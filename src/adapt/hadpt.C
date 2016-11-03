/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: hadpt.C 127 2007-06-07 19:48:25Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#define TARGETPROC -1
//#define FORDEBUG

void refinewrapper(HashTable*HT_Elem_Ptr, HashTable*HT_Node_Ptr, MatProps* matprops_ptr, 
		   ElemPtrList *RefinedList, Element *EmTemp);

extern void  refine(Element*, HashTable*, HashTable*, MatProps* matprops_ptr);

extern void  depchk(Element*, HashTable*, HashTable*, int*, ElemPtrList*);

void update_neighbor_info(HashTable* HT_Elem_Ptr, ElemPtrList* RefinedList, 
			  int myid, int numprocs, 
			  HashTable* HT_Node_Ptr, int h_count );
//extern void  update_neighbor_info(HashTable*, Element*[], int, int,int, HashTable*, int h_count);

extern void  data_com(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		      int myid, int numprocs,
		      int h_count);

extern void  htflush(HashTable*, HashTable*, int);

extern void test_h_refine(HashTable* HT_Elem_Ptr, int myid, int h_count);

extern void all_check(HashTable* eltab, HashTable* ndtab, int myid, int m, double TARGET);


//#define REFINE_THRESHOLD 250000*GEOFLOW_TINY
//#define REFINE_THRESHOLD 10*GEOFLOW_TINY
//#define REFINE_THRESHOLD 35*GEOFLOW_TINY
#define REFINE_THRESHOLD1  5*GEOFLOW_TINY
#define REFINE_THRESHOLD2 15*GEOFLOW_TINY
#define REFINE_THRESHOLD  40*GEOFLOW_TINY


void  H_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int h_count, 
	      double target, MatProps* matprops_ptr, FluxProps *fluxprops, //doesn't need fluxprops
              TimeProps* timeprops_ptr, int num_buffer_layer)
  /*-------------
    scanning element hashtable, if the error of an element is bigger than the target error
    1). check if it has been marked for refinement caused by other element refinement
        if not, store it in the 'refined' array. at the same time, counting the related
	refinement. if it is bigger than a limittation. refuse to do the refinement and
	remove the refinement mark of every related element.
    2). at the same time, if one refinement is required in an element along the interfaces, 
        check its neighbor, if this neighbor also need to be refinement, stop this round of
	checking and refuse to do the refinement. how to make the refinement infomation be
	passed through all the processors is very, very, very difficult.
    3). continue untill every element has been checked


    NOTE:
    this is not the function you should call when a new pile or flux 
    source is placed.  that function is initial_H_adapt.

    *------------------------------------------------------*/

{


  int            myid;
  int            numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  if((myid==TARGETPROC)){//&&(timeprops_ptr->iter==354)){
    printf("entering H_adapt()\n");
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
    printf("After first AssertMeshErrorFree\n");
  }
  
  unsigned elemdebugkey2a[2]={2114123639,2004318068};
  Element* Curr_El=(Element*) HT_Elem_Ptr->lookup(elemdebugkey2a);

  int k, i, j;
  HashEntryPtr   entryp;
  Element*       EmTemp;

  ElemPtrList    RefinedList, TempList;
  int            count=0;
  int            ifg;//--
  int            refine_flag;
  int            htype = 101;//-- 101 is arbitary
  
  htflush(HT_Elem_Ptr, HT_Node_Ptr, 1); 

  int h_begin = 1;
  int h_begin_type = 102;

  delete_unused_elements_nodes(HT_Elem_Ptr, HT_Node_Ptr, myid);

  // must be included to make sure that elements share same side/S_C_CON nodes with neighbors
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  // determine which elements to refine and flag them for refinement
  double geo_target = element_weight(HT_Elem_Ptr, HT_Node_Ptr, myid, numprocs);

  int hash_size=HT_Elem_Ptr->get_no_of_buckets();
  int debug_ref_flag = 0;

  for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
	
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp){
      
      EmTemp = (Element*)(entryp->value);
      assert(EmTemp);
      if(EmTemp->get_adapted_flag()>=NEWSON)
	EmTemp->put_adapted_flag(NOTRECADAPTED);

      entryp=entryp->next;
    }
  }

  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
    
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp){
      
      EmTemp = (Element*)(entryp->value);
      assert(EmTemp);
      //-- this requirement is used to exclude the new elements
      if(((EmTemp->get_adapted_flag()>0)&&
	  (EmTemp->get_adapted_flag()<NEWSON))&&
	 (EmTemp->get_gen()<REFINE_LEVEL)&&
	 ((EmTemp->if_pile_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )>0)||
	  (EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)>0)|| 
	  (EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)>0)|| 
	  (EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD)>0)||
	  (EmTemp->if_source_boundary(HT_Elem_Ptr)>0)||
	  (*(EmTemp->get_el_error())>geo_target)
	  )
	 )
      {
        refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	debug_ref_flag++;
      }
      entryp=entryp->next;
    }
  }
  
  // -h_count for debugging
  //update_neighbor_info(HT_Elem_Ptr, &RefinedList,myid, numprocs, HT_Node_Ptr, h_count);
  refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
		      (void*) &RefinedList,timeprops_ptr);
 
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  if((myid==TARGETPROC)) { //&&(timeprops_ptr->iter==354)){
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
    printf("After third AssertMeshErrorFree\n");
  }
  

  /*************************************************************************/
  /* Add an num_buffer_layer wide buffer layer around the pile/mass-source */
  /*************************************************************************/
  if(num_buffer_layer>=1){
    
    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    //refine where necessary before placing the innermost buffer layer
    for(i=0; i<hash_size; i++){
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){	
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if((EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )==1)||
   	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)==1)||  
   	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)==1)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD)==1)
	   ){
	  refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	  debug_ref_flag++;
	}
	entryp=entryp->next;
      }
    }

    
    refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			(void*) &RefinedList,timeprops_ptr);

    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    //mark the elements in the innermost buffer layer as the BUFFER layer
    for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable	
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if((EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD)>0)
	   )
	  EmTemp->put_adapted_flag(BUFFER);
	entryp=entryp->next;
      }
    }

    for(int ibufferlayer=2;ibufferlayer<=num_buffer_layer;ibufferlayer++)
    {
      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      //refine where necessary before placing the next buffer layer
      for(i=0; i<hash_size; i++)
      {
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp)
        {
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,REFINE_THRESHOLD)==1){
	    refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	    debug_ref_flag++;
	  }
	  entryp=entryp->next;
	}
      }
      
      // -h_count for debugging

      //refine_neigh_update() needs to know new sons are NEWSONs, 
      //can't call them BUFFER until after refine_neigh_update()
      refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			  (void*) &RefinedList,timeprops_ptr);

      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      
      if((myid==TARGETPROC)) {//&&(timeprops_ptr->iter==354)){
	AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
	printf("After fourth AssertMeshErrorFree\n");
      }
      

      //mark the elements just outside the buffer layer as the NEWBUFFER layer
      TempList.trashlist();
      for(i=0; i<hash_size; i++){
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp){
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,
					     REFINE_THRESHOLD)>0) 
	    TempList.add(EmTemp);
	  //EmTemp->put_adapted_flag(NEWBUFFER);
	    

	  entryp=entryp->next;
	}
      }

      //remark the NEWBUFFER elements as BUFFER element (move them from NEWBUFFER to BUFFER)
      for(i=0;i<TempList.get_num_elem();i++) {
	EmTemp=TempList.get(i);
	assert(EmTemp);
	EmTemp->put_adapted_flag(BUFFER);
      }
      TempList.trashlist();

    }
  }
  
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  htflush(HT_Elem_Ptr, HT_Node_Ptr, 2); 

  for(i=0; i<hash_size; i++) {
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;      	    

      switch(EmTemp->get_adapted_flag()) {
      case NEWBUFFER:
	printf("Suspicious element has adapted flag=%d\n aborting",EmTemp->get_adapted_flag());
	assert(0);
	break;
      case BUFFER:
      case NEWSON:
      case NEWFATHER:
      case NOTRECADAPTED:
	//it's an active (non ghost) element
	EmTemp->calc_d_gravity(HT_Elem_Ptr);
	EmTemp->calc_wet_dry_orient(HT_Elem_Ptr);
	break;
      case TOBEDELETED:
	//deleting the refined father elements but not ghost element so don't need to call move_data() again
	EmTemp->void_bcptr();
	HT_Elem_Ptr->remove(EmTemp->pass_key(), 1,stdout,myid,20);
	delete EmTemp;	    	  
	break;
      case -NOTRECADAPTED:
      case -NEWFATHER:
      case -NEWSON:
      case -BUFFER:
	//it's a ghost element, keep these so I don't have to move data again.
	break;
      case OLDFATHER:
      case OLDSON:
	printf("Suspicious element has adapted flag=%d\n aborting",EmTemp->get_adapted_flag());
	assert(0);
	break;
      default:
	//I don't know what kind of Element this is.
	printf("FUBAR element type in H_adapt()!!! key={%u,%u} adapted=%d\naborting.\n",
	       *(EmTemp->pass_key()+0),*(EmTemp->pass_key()+1),
	       EmTemp->get_adapted_flag());
	assert(0);
	break;
      }

    }
  }
  
  return;
}

#ifdef DISABLED
//Keith wrote this because the code block was repeated so many times
void refinewrapper(HashTable*HT_Elem_Ptr, HashTable*HT_Node_Ptr, MatProps* matprops_ptr, 
		   Element **refined, int *count, Element *EmTemp){

  int j, sur = 0; int mii = 0;

  while(refined[mii])
  {
    sur = 1;
    for(j=0;j<KEYLENGTH;j++)
      if(*(refined[mii]->pass_key()+j) != *(EmTemp->pass_key()+j))
	sur = 0;
    if(sur == 1)
      break;
    mii++;
  }
	  
  if(!sur)
  {
    int ifg = 1; 
    j = 0;
    //-- check out the triggered refinement
    depchk(EmTemp, HT_Elem_Ptr, HT_Node_Ptr, &ifg, refined, count);
    if(ifg)
    {
      while(refined[j])
      {
	if(!refined[j]->get_refined_flag())
	  refine(refined[j], HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr); 
	j++; 
      }
    }
  }

  return;
}
#endif
//Keith wrote this because the code block was repeated so many times
void refinewrapper(HashTable*HT_Elem_Ptr, HashTable*HT_Node_Ptr, MatProps* matprops_ptr, 
		   ElemPtrList *RefinedList, Element *EmTemp)
{

  int sur = 0, ifg=1, ielem;

  for(ielem=0;ielem<RefinedList->get_num_elem();ielem++) 
    if(sur=compare_key(RefinedList->get_key(ielem),EmTemp->pass_key()))
      break;

  if(!sur) {	  
    //-- check out the triggered refinement
    depchk(EmTemp,HT_Elem_Ptr,HT_Node_Ptr,&ifg,RefinedList);
    
    if(ifg)
      for(ielem=0;ielem<RefinedList->get_num_elem();ielem++)
	if(!((RefinedList->get(ielem))->get_refined_flag())) {
	  refine(RefinedList->get(ielem),HT_Elem_Ptr,HT_Node_Ptr, 
		 matprops_ptr); 
	  (RefinedList->get(ielem))->put_adapted_flag(OLDFATHER);
	  (RefinedList->get(ielem))->put_refined_flag(1);
	}    
  }
  return;
}

//from init_piles.C
extern void elliptical_pile_height(HashTable* HT_Node_Ptr, Element *EmTemp, 
			    MatProps* matprops_ptr, PileProps* pileprops_ptr);


void  initial_H_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int h_count, 
		      MatProps* matprops_ptr, PileProps *pileprops_ptr,
		      FluxProps *fluxprops_ptr, TimeProps* timeprops_ptr,
		      int num_buffer_layer)
{

  int k, i, j;
  HashEntryPtr   entryp;
  Element*       EmTemp;

  ElemPtrList    RefinedList, TempList;
  int            count=0;
  int            myid;
  int            numprocs;
  int            ifg;//--
  int            refine_flag;
  unsigned       sent_buf[4*KEYLENGTH];//-- 1 source; 2 target; 3 refined flag; 4 generation
  unsigned       recv_buf[4*KEYLENGTH];//-- same
  int            htype = 101;//-- 101 is arbitary
  unsigned       NodeDebugKey[2]={3489660928,         0};

  /*
  if(myid==TARGETPROC){
    printf("Entered initial_H_adapt()\n");
    NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
  }
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  //move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  //printf("myid=%d init_H_adapt 1\n",myid);
  //AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);

  int num_ellipse_centers=0;

  if(timeprops_ptr->iter==0){
    num_ellipse_centers+=pileprops_ptr->numpiles;

    for(int isrc=0; isrc<fluxprops_ptr->no_of_sources; isrc++) {
      if(fluxprops_ptr->start_time[isrc]==0.0)
	num_ellipse_centers++;
    }
  }
  else{
    double begoftimestep=timeprops_ptr->time-timeprops_ptr->dtime;
    double endoftimestep=timeprops_ptr->time;
    for(int isrc=0; isrc<fluxprops_ptr->no_of_sources; isrc++)
      if((begoftimestep<=fluxprops_ptr->start_time[isrc])&&
	 (fluxprops_ptr->start_time[isrc]<endoftimestep))
	num_ellipse_centers++;
  }

  double **xycenter=CAllocD2(num_ellipse_centers,2);
  //printf("num_ellipse_centers=%d\n",num_ellipse_centers);

  int icenter=0;
  if(timeprops_ptr->time==0) {
    for(int ipile=0; ipile<pileprops_ptr->numpiles; ipile++, icenter++){
      xycenter[icenter][0]=pileprops_ptr->xCen[ipile];
      xycenter[icenter][1]=pileprops_ptr->yCen[ipile];
    }

    for(int isrc=0; isrc<fluxprops_ptr->no_of_sources; isrc++)
      if(fluxprops_ptr->start_time[isrc]==0.0) {
	xycenter[icenter][0]=fluxprops_ptr->xCen[isrc];
	//*(matprops_ptr->LENGTH_SCALE);;
	xycenter[icenter][1]=fluxprops_ptr->yCen[isrc];
	//*(matprops_ptr->LENGTH_SCALE);;
	icenter++;
      }
  }
  else{
    double begoftimestep=timeprops_ptr->time-timeprops_ptr->dtime;
    double endoftimestep=timeprops_ptr->time;
    for(int isrc=0; isrc<fluxprops_ptr->no_of_sources; isrc++)
      if((begoftimestep<=fluxprops_ptr->start_time[isrc])&&
	 (fluxprops_ptr->start_time[isrc]<endoftimestep)) {
	xycenter[icenter][0]=fluxprops_ptr->xCen[isrc];
	//*(matprops_ptr->LENGTH_SCALE);
	xycenter[icenter][1]=fluxprops_ptr->yCen[isrc];
	//*(matprops_ptr->LENGTH_SCALE);
	icenter++;
      }	
  }

  htflush(HT_Elem_Ptr, HT_Node_Ptr, 1); 

  //for(k=0;k<297200;k++) refined[k] = 0;
  int h_begin = 1;
  int h_begin_type = 102;
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  /*
  if(myid==0) {
    printf("before AssertMeshErrorFree() 0.8\n");
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
    printf("after AssertMeshErrorFree() 0.8\n");
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  */

  delete_unused_elements_nodes(HT_Elem_Ptr, HT_Node_Ptr, myid);

  /*
  if(myid==TARGETPROC){
    printf("After delete_unused_elements_nodes()\n");
    NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
  }
  */


  //printf("init_H_adapt 2\n");
  //AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);

  // must be included to make sure that elements share same side/S_C_CON nodes with neighbors
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  /*
  if(myid==0) {
    printf("before AssertMeshErrorFree() 0.9\n");
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
    printf("after AssertMeshErrorFree() 0.9\n");
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  */

  //printf("myid=%d init_H_adapt 3\n",myid);
  //AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);

  //MPI_Barrier(MPI_COMM_WORLD); //-- every process arrive here

  int hash_size=HT_Elem_Ptr->get_no_of_buckets();
  int debug_ref_flag = 0;
  
  int intswap;
  int icounter=0;
  int mincentergen;
  do{
    mincentergen=REFINE_LEVEL;
    for(i=0; i<hash_size; i++) {//-- every process begin to scan their own hashtable
	
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){
	
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if(EmTemp->get_adapted_flag()>=NEWSON)
	  EmTemp->put_adapted_flag(NOTRECADAPTED);
	entryp=entryp->next;
      }
    }
    

    TempList.trashlist();
    for(i=0; i<hash_size; i++) {//-- every process begin to scan their own hashtable	
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){
	
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);

	if(EmTemp->get_adapted_flag()>0) {
	  double minx=(*(EmTemp->get_coord()+0)-(*(EmTemp->get_dx()+0))*0.5);
	  //*(matprops_ptr->LENGTH_SCALE);
	  double maxx=(*(EmTemp->get_coord()+0)+(*(EmTemp->get_dx()+0))*0.5);
	  //*(matprops_ptr->LENGTH_SCALE);
	  double miny=(*(EmTemp->get_coord()+1)-(*(EmTemp->get_dx()+1))*0.5);
	  //*(matprops_ptr->LENGTH_SCALE);
	  double maxy=(*(EmTemp->get_coord()+1)+(*(EmTemp->get_dx()+1))*0.5);
	  //*(matprops_ptr->LENGTH_SCALE);
	  
	  //printf("x=[%g,%g] y=[%g,%g]\n",minx,maxx,miny,maxy);

	  for(icenter=0;icenter<num_ellipse_centers;icenter++)
	    if((minx<=xycenter[icenter][0])&&
	       (xycenter[icenter][0]<=maxx)&&
	       (miny<=xycenter[icenter][1])&&
	       (xycenter[icenter][1]<=maxy)){

	      if((EmTemp->get_gen()<REFINE_LEVEL)&&
		 (EmTemp->get_adapted_flag()>0)&&
		 (EmTemp->get_adapted_flag()<NEWSON))
		refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);  
	      Element *EmSon=EmTemp;
	      if(EmTemp->get_adapted_flag()==OLDFATHER)
		for(int ison=0;ison<4;ison++) {
		  EmSon=(Element *) HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*ison);
		  TempList.add(EmSon);
		}
	      else
		TempList.add(EmTemp);
	      
	      if(mincentergen>(EmSon->get_gen())) mincentergen=EmSon->get_gen();
	      break;
	    }
	}

	entryp=entryp->next;
      }
    }

    intswap=mincentergen;
    MPI_Allreduce(&mincentergen,&intswap,1, MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    mincentergen=intswap;


    //printf("icounter=%d\n",icounter);
    //update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);

    //refine_neigh_update() needs to know new sons are NEWSONs, 
    //can't call them BUFFER until after refine_neigh_update()
    refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			(void*) &RefinedList,timeprops_ptr);

    //mark the new buffer elements as BUFFER element, 
    for(i=0;i<TempList.get_num_elem();i++) {
      EmTemp=TempList.get(i);
      assert(EmTemp);
      EmTemp->put_adapted_flag(BUFFER);
    }
    TempList.trashlist();

    /*
    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
    if(myid==TARGETPROC){
      printf("After refine_neigh_update() 1\n");
      NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
    }
    */

    // -h_count for debugging

    if(numprocs>1){
      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      /*
	if(myid==0) {
	  printf("before AssertMeshErrorFree() 1.0\n");
	  AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
	  printf("after AssertMeshErrorFree() 1.0\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
      */
      
      //put in a small buffer layer to facilitate refinement across interprocessor boundaries
      for(int ibufferlayer=0;ibufferlayer<2;ibufferlayer++){    
	
	//refine where necessary and store the next buffer layer in TempList
	TempList.trashlist();
	for(i=0; i<hash_size; i++){
	  entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	  while(entryp){
	    EmTemp = (Element*)(entryp->value);
	    assert(EmTemp);
	    if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,
					       REFINE_THRESHOLD)==1){
	      
	      if(EmTemp->get_gen()<REFINE_LEVEL)
		refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);  
	      if(EmTemp->get_adapted_flag()==OLDFATHER)
		for(int ison=0;ison<4;ison++) {
		  Element *EmSon=(Element *) HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*ison);
		  TempList.add(EmSon);
		}
	      else{
		TempList.add(EmTemp);
	      }	      
	      debug_ref_flag++;
	    }
	    entryp=entryp->next;
	  }
	}
	

	//update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);
	//refine_neigh_update() needs to know new sons are NEWSONs, 
	//can't call them BUFFER until after refine_neigh_update()
        refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			    (void*) &RefinedList,timeprops_ptr);

	//mark the new buffer elements as BUFFER element, 
	for(i=0;i<TempList.get_num_elem();i++) {
	  EmTemp=TempList.get(i);
	  assert(EmTemp);
	  EmTemp->put_adapted_flag(BUFFER);
	}
	TempList.trashlist();



	move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
	/*
	if(myid==0) {
	  printf("before AssertMeshErrorFree() 2.0\n");
	  AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
	  printf("after AssertMeshErrorFree() 2.0\n");

	}

	MPI_Barrier(MPI_COMM_WORLD);
	*/
	/*
	//move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
    if(myid==TARGETPROC){
      printf("After refine_neigh_update() 2\n");
      NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
    }
	*/
      }
    }
  }while(((++icounter)<(REFINE_LEVEL+3))&& //failsafe to prevent infinite loop
	 (mincentergen<REFINE_LEVEL) //the usual criteria
	 );
  //printf("icounter=%d\n",icounter);

  CDeAllocD2(xycenter);

  /*
  if(myid==TARGETPROC){
    printf("After center refined()\n");
    NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
  }
  */

  icounter=0;
  int minboundarygen;

  do{
    /*
    if(myid==0) {
      printf("before AssertMeshErrorFree() 2.2\n");
      AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
      printf("after AssertMeshErrorFree() 2.2\n");
      
    }    
    MPI_Barrier(MPI_COMM_WORLD);
    */

    minboundarygen=REFINE_LEVEL;
    for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
	
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){
	
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if(EmTemp->get_adapted_flag()>=NEWSON)
	  EmTemp->put_adapted_flag(NOTRECADAPTED);
	
	entryp=entryp->next;
      }
    }
    
    /*
    if(myid==0) {
      printf("before AssertMeshErrorFree() 2.3\n");
      AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
      printf("after AssertMeshErrorFree() 2.3\n");
      
    }    
    MPI_Barrier(MPI_COMM_WORLD);
    */

    //if(myid==0) printf("initial_H_adapt %d\n",1); fflush(stdout)
    //if(myid==0) printf("initial_H_adapt %d\n",2); fflush(stdout)

    if(timeprops_ptr->iter==0)
      //mark the piles so we can refine their boundaries
      for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp){	
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  elliptical_pile_height(HT_Node_Ptr, EmTemp, matprops_ptr, pileprops_ptr);
	  entryp=entryp->next;
	}
      }
    
    /*
    if(myid==0) {
      printf("before AssertMeshErrorFree() 2.4\n");
      AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
      printf("after AssertMeshErrorFree() 2.4\n");
      
    }    
    MPI_Barrier(MPI_COMM_WORLD);
    */

    //if(myid==0) printf("initial_H_adapt %d\n",3); fflush(stdout)

    //mark the mass flux sources so we can refine their boundaries
    mark_flux_region(HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr, fluxprops_ptr, 
		     timeprops_ptr);
    
    //if(myid==0) printf("initial_H_adapt %d\n",4); fflush(stdout)


    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    /*
    printf("myid=%d Initial_H_adapt() before first if_pile_boundary()\n",myid);
    if(myid==0) {
      printf("before AssertMeshErrorFree() 3.0\n");
      AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
      printf("after AssertMeshErrorFree() 3.0\n");
      
    }    
    MPI_Barrier(MPI_COMM_WORLD);

    unsigned checkthiselement[2]={2092040192,0};
    unsigned checkthisneighbor[2]={2112662186,2863311530};
    
    ElemBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,checkthiselement,stdout);
    ElemBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,checkthisneighbor,stdout);
    */

    for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){	

	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	//-- this requirement is used to exclude the new elements
	if(((EmTemp->get_adapted_flag()>0)&&
	    (EmTemp->get_adapted_flag()<NEWSON))&&
	   (EmTemp->get_gen()<REFINE_LEVEL)){
	  if(((EmTemp->if_pile_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )>0)|| 
	    //(EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)>0)|| 
	    //(EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)>0)|| 
	      (EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD )>0)|| 
	      (EmTemp->if_source_boundary(HT_Elem_Ptr)>0))
	     ){
	    refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	    debug_ref_flag++;
	  }
	}
	entryp=entryp->next;
      }
    }
    //printf("myid=%d Initial_H_adapt() after first if_pile_boundary()\n",myid);

    //if(myid==0) printf("initial_H_adapt %d\n",5); fflush(stdout)

    // -h_count for debugging
    //update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);

    //refine_neigh_update() needs to know new sons are NEWSONs, 
    //can't call them BUFFER until after refine_neigh_update()
    refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			(void*) &RefinedList,timeprops_ptr);

    /*
    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
    if(myid==TARGETPROC){
      printf("After refine_neigh_update() 3\n");
      NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
    }
    */


    //if(myid==0) printf("initial_H_adapt %d\n",6); fflush(stdout)

    if(timeprops_ptr->iter==0)
      //mark the piles so we can refine their boundaries
      for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp){	
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  elliptical_pile_height(HT_Node_Ptr, EmTemp, matprops_ptr, pileprops_ptr);
	  entryp=entryp->next;
	}
      }

    //if(myid==0) printf("initial_H_adapt %d\n",7); fflush(stdout)

    //mark the mass flux sources so we can refine their boundaries
    mark_flux_region(HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr, fluxprops_ptr, 
		     timeprops_ptr);

    //if(myid==0) printf("initial_H_adapt %d\n",8); fflush(stdout)

    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    /*
    printf("myid=%d Initial_H_adapt() before second if_pile_boundary()\n",myid);
    ElemBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,checkthiselement,stdout);
    ElemBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,checkthisneighbor,stdout);
    */

    for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){	

	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	//-- this requirement is used to exclude the new elements
	if(EmTemp->get_adapted_flag()>0)
	  if(
	     ((EmTemp->if_pile_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )>0)|| 
	      //(EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)>0)|| 
	      //(EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)>0)|| \
	      (EmTemp->if_pile_boundary(HT_Elem_Ptr,REFINE_THRESHOLD )>0)|| 
	      (EmTemp->if_source_boundary(HT_Elem_Ptr)>0))
	     ){
	    EmTemp->put_adapted_flag(BUFFER);
	    if(minboundarygen>EmTemp->get_gen()) minboundarygen=EmTemp->get_gen();
	  }
	
	entryp=entryp->next;
      }
    }
    //printf("myid=%d Initial_H_adapt() after second if_pile_boundary()\n",myid);


    intswap=minboundarygen;
    MPI_Allreduce(&minboundarygen,&intswap,1, MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    minboundarygen=intswap;

    //if(myid==0) printf("initial_H_adapt %d\n",9); fflush(stdout)

    if(numprocs>1){
      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      //if(myid==0) printf("initial_H_adapt %d\n",10); fflush(stdout)

      //put in a small buffer layer to facilitate refinement across interprocessor boundaries
      for(int ibufferlayer=0;ibufferlayer<2;ibufferlayer++){    
	
	//refine where necessary and store the next buffer layer in TempList
	TempList.trashlist();
	for(i=0; i<hash_size; i++){
	  entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	  while(entryp){
	    EmTemp = (Element*)(entryp->value);
	    assert(EmTemp);
	    if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,
					       REFINE_THRESHOLD)==1){
	      
	      if(EmTemp->get_gen()<REFINE_LEVEL)
		refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);  
	      if(EmTemp->get_adapted_flag()==OLDFATHER)
		for(int ison=0;ison<4;ison++) {
		  Element *EmSon=(Element *) HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*ison);
		  TempList.add(EmSon);
		}
	      else{
		TempList.add(EmTemp);
	      }
	      
	      debug_ref_flag++;
	    }
	    entryp=entryp->next;
	  }
	}

	//if(myid==0) printf("initial_H_adapt %d\n",11); fflush(stdout)
	
	//update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);

	//refine_neigh_update() needs to know new sons are NEWSONs, 
	//can't call them BUFFER until after refine_neigh_update()
        refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			    (void*) &RefinedList,timeprops_ptr);
	
	//if(myid==0) printf("initial_H_adapt %d\n",12); fflush(stdout)

	//mark the new buffer elements as BUFFER element 
	for(i=0;i<TempList.get_num_elem();i++) {
	  EmTemp=TempList.get(i);
	  assert(EmTemp);
	  EmTemp->put_adapted_flag(BUFFER);
	}
	TempList.trashlist();
	

	//if(myid==0) printf("initial_H_adapt %d\n",13); fflush(stdout)


	move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

	/*
	//move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
    if(myid==TARGETPROC){
      printf("After refine_neigh_update() 4\n");
      NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
    }
	*/
	//if(myid==0) printf("initial_H_adapt %d\n",14); fflush(stdout)

      }
    }

  }while(((++icounter)<(REFINE_LEVEL+3))&& //failsafe to prevent infinite loop 
	 (minboundarygen<REFINE_LEVEL) //the usual criteria
	 );

  /*
  if(myid==TARGETPROC){
    printf("After boundary refined\n");
    NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
  }
  */

  //if(myid==0) printf("initial_H_adapt %d\n",15); fflush(stdout)

  //printf("icounter=%d\n",icounter);

  //unmark the newsons and buffer layer elements
  for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable
	
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp){
      
      EmTemp = (Element*)(entryp->value);
      assert(EmTemp);
      if(EmTemp->get_adapted_flag()>=NEWSON)
	EmTemp->put_adapted_flag(NOTRECADAPTED);

      entryp=entryp->next;
    }
  }

  //add the real buffer layer
  /*************************************************************************/
  /* Add an num_buffer_layer wide buffer layer around the pile/mass-source */
  /*************************************************************************/
  if(num_buffer_layer>=1){
    
    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    //refine where necessary before placing the innermost buffer layer
    for(i=0; i<hash_size; i++){
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){	
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if((EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )==1)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)==1)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)==1)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD )==1)
	   ){
	  refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	  debug_ref_flag++;
	}
	entryp=entryp->next;
      }
    }
    
    //update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);
    //refine_neigh_update() needs to know new sons are NEWSONs, 
    //can't call them BUFFER until after refine_neigh_update()
    refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			(void*) &RefinedList,timeprops_ptr);

    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

    /*
    //move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
    if(myid==TARGETPROC){
      printf("After refine_neigh_update() 5\n");
      NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
    }
    */

    //mark the elements in the innermost buffer layer as the BUFFER layer
    for(i=0; i<hash_size; i++){//-- every process begin to scan their own hashtable	
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){
	EmTemp = (Element*)(entryp->value);
	assert(EmTemp);
	if((EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,GEOFLOW_TINY     )>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD1)>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD2)>0)||
	   (EmTemp->if_first_buffer_boundary(HT_Elem_Ptr,REFINE_THRESHOLD )>0)	   )
	  EmTemp->put_adapted_flag(BUFFER);
	entryp=entryp->next;
      }
    }

    //increase the width of the buffer layer by one element at a time
    //until it's num_buffer_layer Elements wide
    for(int ibufferlayer=2;ibufferlayer<=num_buffer_layer;ibufferlayer++){    

      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      //refine where necessary before placing the next buffer layer
      for(i=0; i<hash_size; i++){
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp){
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,
					     REFINE_THRESHOLD)==1){
	    refinewrapper(HT_Elem_Ptr,HT_Node_Ptr,matprops_ptr,&RefinedList,EmTemp);
	    debug_ref_flag++;
	  }
	  entryp=entryp->next;
	}
      }
      
      // -h_count for debugging
      //update_neighbor_info(HT_Elem_Ptr, &RefinedList, myid, numprocs, HT_Node_Ptr, h_count);
      refine_neigh_update(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,
			  (void*) &RefinedList,timeprops_ptr);

      move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

      /*
      //move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr); //this move_data() only here for debug
      if(myid==TARGETPROC){
	printf("After refine_neigh_update() 6\n");
	NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
      }
      */

      //the elements just outside the buffer layer will be the next 
      //buffer layer, store them in TempList
      TempList.trashlist();
      for(i=0; i<hash_size; i++){
	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
	while(entryp){
	  EmTemp = (Element*)(entryp->value);
	  assert(EmTemp);
	  if(EmTemp->if_next_buffer_boundary(HT_Elem_Ptr,HT_Node_Ptr,
					     REFINE_THRESHOLD)>0)
	    TempList.add(EmTemp);
	  entryp=entryp->next;
	}
      }

      //mark the new buffer elements as BUFFER element 
      for(i=0;i<TempList.get_num_elem();i++) {
	EmTemp=TempList.get(i);
	assert(EmTemp);
	EmTemp->put_adapted_flag(BUFFER);
      }
      TempList.trashlist();

    }
  }

  /*
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  if(myid==TARGETPROC){
    printf("After buffer layer placed\n");
    NodeBackgroundCheck(HT_Elem_Ptr,HT_Node_Ptr,NodeDebugKey);
  }
  */

  /* transfer ghost elements to proper processors */
 
  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  htflush(HT_Elem_Ptr, HT_Node_Ptr, 2); 

#ifdef FORDEBUG  

  for(i=0; i<hash_size; i++) {
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;      	    

      if((-BUFFER<=EmTemp->get_adapted_flag())&& 
	 (EmTemp->get_adapted_flag()<=-NOTRECADAPTED)) {
	EmTemp->void_bcptr();
	HT_Elem_Ptr->remove(EmTemp->pass_key(), 1,stdout,myid,21);
	delete EmTemp;	    	  
      }
    }
  }

  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);

  Node* NdTemp;
  for(i=0; i<HT_Node_Ptr->get_no_of_buckets(); i++) {
    entryp = *(HT_Node_Ptr->getbucketptr() + i);
    while(entryp) {
      NdTemp=(Node*)(entryp->value);
      entryp=entryp->next;      
      assert(NdTemp);

      NdTemp->put_num_assoc_elem(0);
    }
  }
  int inode;
#endif


  for(i=0; i<hash_size; i++) {
    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;      	    

      switch(EmTemp->get_adapted_flag()) {
      case NEWBUFFER:
	printf("Suspicious element has adapted flag=%d\n aborting",EmTemp->get_adapted_flag());
	assert(0);
	break;
      case BUFFER:
      case NEWSON:
      case NEWFATHER:
      case NOTRECADAPTED:
	//it's an active (non ghost) element
	EmTemp->calc_d_gravity(HT_Elem_Ptr);
	EmTemp->calc_wet_dry_orient(HT_Elem_Ptr);

#ifdef FORDEBUG
	NdTemp=(Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());
	assert(NdTemp);
	NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);
	
	for(inode=0;inode<8;inode++) {
	  NdTemp=(Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+inode*KEYLENGTH);
	  assert(NdTemp);
	  NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);
	}
#endif
	break;
      case TOBEDELETED:
	//deleting the refined father elements but ghost element so don't need to call move_data() again before AssertMeshErrorFree
	EmTemp->void_bcptr();
	HT_Elem_Ptr->remove(EmTemp->pass_key(), 1,stdout,myid,22);
	delete EmTemp;	    	  
	break;
      case -NOTRECADAPTED:
      case -NEWFATHER:
      case -NEWSON:
      case -BUFFER:
	//it's a ghost element, keep these so I don't have to move data again.
	break;
      case OLDFATHER:
      case OLDSON:
	printf("Suspicious element has adapted flag=%d\n aborting",EmTemp->get_adapted_flag());
	assert(0);
	break;
      default:
	//I don't know what kind of Element this is.
	printf("FUBAR element type in H_adapt()!!! key={%u,%u} adapted=%d\naborting.\n",
	       *(EmTemp->pass_key()+0),*(EmTemp->pass_key()+1),
	       EmTemp->get_adapted_flag());
	assert(0);
	break;
      }

    }
  }

#ifdef FORDEBUG
  for(i=0; i<HT_Node_Ptr->get_no_of_buckets(); i++) {
    entryp = *(HT_Node_Ptr->getbucketptr() + i);
    while(entryp) {
      NdTemp=(Node*)(entryp->value);
      entryp=entryp->next;      
      assert(NdTemp);

      if(NdTemp->get_num_assoc_elem()==0){
	HT_Node_Ptr->remove(NdTemp->pass_key(), 0,stdout,myid,23);
	delete NdTemp;	    	  
      }
    }
  }

  if(myid==TARGETPROC){
    printf("in initial_H_adapt() before repartition\n");
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
  }
#endif
  
  if(numprocs>1)
    repartition2(HT_Elem_Ptr, HT_Node_Ptr, timeprops_ptr);

  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);



#ifdef FORDEBUG  
  if(myid==TARGETPROC){
    printf("in initial_H_adapt() after repartition\n");
    AssertMeshErrorFree(HT_Elem_Ptr,HT_Node_Ptr,numprocs,myid,0.0);
  }
#endif  
  return;

  printf("myid=%d reached the end of initial_H_adapt\n",myid);
  exit(1);
}
/***********************************************************************/

void  H_adapt_to_level(HashTable* El_Table, HashTable* NodeTable, 
		       MatProps* matprops_ptr, PileProps* pileprops_ptr,
		       FluxProps *fluxprops_ptr, TimeProps* timeprops_ptr, 
		       int refinelevel) {

  if(refinelevel>REFINE_LEVEL) refinelevel=REFINE_LEVEL;  

  int myid, numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  int minrefinelevel;
  Element* EmTemp;
  ElemPtrList RefinedList(2048);
  int i, generation;
  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;

  htflush(El_Table,NodeTable, 1); 

  do{

    for(i=0; i<num_buck; i++)
      if(*(buck+i)) {
	
	currentPtr = *(buck+i);
	while(currentPtr){
	  
	  EmTemp=(Element*)(currentPtr->value);
	  currentPtr=currentPtr->next;     
	  assert(EmTemp);
	  
	  if(EmTemp->get_adapted_flag()>=NOTRECADAPTED)
	    EmTemp->put_adapted_flag(NOTRECADAPTED);
	  
	}
      }
    
    minrefinelevel=refinelevel;
    RefinedList.trashlist();
    for(i=0; i<num_buck; i++)
      if(*(buck+i)) {
	
	currentPtr = *(buck+i);
	while(currentPtr){
	
	  EmTemp=(Element*)(currentPtr->value);
	  currentPtr=currentPtr->next;     
	  assert(EmTemp);
	  
	  generation=EmTemp->get_gen();
	  
	  if((EmTemp->get_adapted_flag()==NOTRECADAPTED)&&
	     (generation<refinelevel)) {
	    refinewrapper(El_Table,NodeTable,matprops_ptr,&RefinedList,EmTemp);
	    if(generation<minrefinelevel)
	      minrefinelevel=generation;
	  }
	}
      }

    refine_neigh_update(El_Table,NodeTable,numprocs,myid,
			(void*) &RefinedList,timeprops_ptr);

    move_data(numprocs,myid,El_Table,NodeTable,timeprops_ptr);

  }while(minrefinelevel<refinelevel);


  if(timeprops_ptr->iter==0) 
  {
    if(pileprops_ptr->numpiles>0)
    {
      for(i=0; i<num_buck; i++)
	if(*(buck+i))
	{
	  currentPtr = *(buck+i);
	  while(currentPtr){
	    
	    EmTemp=(Element*)(currentPtr->value);
	    currentPtr=currentPtr->next;     
	    assert(EmTemp);
	    
	    if(EmTemp->get_adapted_flag()>=NOTRECADAPTED)
            {
	      EmTemp->put_adapted_flag(NOTRECADAPTED);
	      elliptical_pile_height(NodeTable,EmTemp,matprops_ptr, 
				     pileprops_ptr);
	    }
	    else{
	      EmTemp->void_bcptr();
	      El_Table->remove(EmTemp->pass_key(), 1,stdout,myid,24);
	      delete EmTemp;	   
	    }	    
	  }
	}
    }
    else
    {
      for(i=0; i<num_buck; i++)
	if(*(buck+i)) {
	  
	  currentPtr = *(buck+i);
	  while(currentPtr){
	    
	    EmTemp=(Element*)(currentPtr->value);
	    currentPtr=currentPtr->next;     
	    assert(EmTemp);
	    
	    if(EmTemp->get_adapted_flag()>=NOTRECADAPTED) {
	      EmTemp->put_adapted_flag(NOTRECADAPTED);
	      EmTemp->put_height(0.0);
	    }
	    else{
	      EmTemp->void_bcptr();
	      El_Table->remove(EmTemp->pass_key(), 1,stdout,myid,25);
	      delete EmTemp;	   
	    }	    
	  }
	}
    }
  }

  mark_flux_region(El_Table,NodeTable,matprops_ptr,fluxprops_ptr, 
		   timeprops_ptr);

  if(numprocs>1)
    repartition2(El_Table,NodeTable,timeprops_ptr);

  move_data(numprocs,myid,El_Table,NodeTable,timeprops_ptr);
  for(i=0; i<<num_buck; i++) {
    currentPtr = *(El_Table->getbucketptr() + i);
    while(currentPtr) {
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;      	    
      assert(EmTemp);

      if((EmTemp->get_adapted_flag()>TOBEDELETED)&&
	 (EmTemp->get_adapted_flag()<=BUFFER))
	EmTemp->calc_wet_dry_orient(El_Table);
    }
  }

  return;
}
