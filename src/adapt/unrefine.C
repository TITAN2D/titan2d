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
 * $Id: unrefine.C 150 2007-06-27 20:28:42Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

//#define MIN_GENERATION -1
#define TARGET_PROC -1


int IfMissingElem(HashTable* El_Table, int myid, int iter, int isearch){

  return(0);

  if(iter<4) return(0);

  int i,yada=0;
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)) {
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr) {	      
	Element* Curr_El = (Element*) currentPtr->value;
	currentPtr=currentPtr->next;
	if((*(Curr_El->pass_key()+0)==2110300160)&&
	   (*(Curr_El->pass_key()+1)==0)
	   )
	  {
	    printf("search %d found it, adapted=%d\n",isearch,Curr_El->get_adapted_flag());
	    //printf("search %d found it, adapted=%d, enter an int\n",isearch,Curr_El->get_adapted_flag());
	    //scanf("%d",&yada);
	    return(0);
	  }
      }
    }

  return(1);
}



void unrefine(HashTable* El_Table, HashTable* NodeTable, double target,
	      int myid, int nump, TimeProps* timeprops_ptr, MatProps* matprops_ptr) 
{

  //printf("myid=%d entering unrefine\n",myid);

  int time_step=timeprops_ptr->iter;

  int i, j, k;
  Element* Curr_El;
  ElemPtrList NewFatherList, OtherProcUpdate;

  //-------------------go through all the elements of the subdomain------------------------
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {	      
	    Curr_El = (Element*) currentPtr->value;
	    currentPtr=currentPtr->next;
	    if(Curr_El->get_adapted_flag()==NEWFATHER)
	      Curr_El->put_adapted_flag(NOTRECADAPTED);
	  }
      }


  // start unrefinement
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)) {
      HashEntryPtr currentPtr = *(buck+i);
      if(i==6)
	j = 1;
      while(currentPtr) {	      
	Curr_El = (Element*) currentPtr->value;
	//  need to get currentPtr->next now since currentPtr might get deleted!
	currentPtr=currentPtr->next;
	if(Curr_El->get_adapted_flag()==NOTRECADAPTED) {//if this is a refined element don't involve!!! 
	      
	  // if this if the original element, don't unrefine.  only son 0 checks for unrefinement!
	  if((Curr_El->get_gen()>MIN_GENERATION)&& 
	     (Curr_El->get_which_son()==0)
	     ) {
	    //check to see if currentPtr might get deleted and if it might, find next ptr that won't
	    if(currentPtr != NULL) {
	      int newnext = 0;
	      while(newnext == 0 && currentPtr != NULL) {
		Element* nextelm = (Element*) currentPtr->value;
		if(nextelm->get_which_son() == 0)
		  newnext = 1;
		else 
		  currentPtr = currentPtr->next;
	      }
	    }
	    Curr_El->find_brothers(El_Table, NodeTable,
				   target, myid, matprops_ptr,
				   &NewFatherList, &OtherProcUpdate);

	  }
	}
      }
    }


  int iproc;  
  time_t tic, toc;

  /*
  for(iproc=0;iproc<nump-1;iproc++)
    if(myid>iproc); MPI_Barrier(MPI_COMM_WORLD);

  printf("myid=%d unref before unref_neigh_update\n",myid); fflush(stdout);

  tic=time(NULL);
  toc=tic+2;
  do{
    tic=time(NULL);
  }while(tic<toc);

  for(iproc=1;iproc<nump;iproc++)
    if(myid<iproc); MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  */
  //assert(!IfMissingElem(El_Table, myid, time_step, 0));
  
  unrefine_neigh_update(El_Table,NodeTable,myid,(void*) &NewFatherList);

  //assert(!IfMissingElem(El_Table, myid, time_step, 1));

  unrefine_interp_neigh_update(El_Table,NodeTable,nump,myid,
			       (void*) &OtherProcUpdate);



  /*
  for(iproc=0;iproc<nump-1;iproc++)
    if(myid>iproc); MPI_Barrier(MPI_COMM_WORLD);

  printf("myid=%d unref before delete\n",myid); fflush(stdout);

  tic=time(NULL);
  toc=tic+2;
  do{
    tic=time(NULL);
  }while(tic<toc);

  for(iproc=1;iproc<nump;iproc++)
    if(myid<iproc); MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
  */

  for(k=0;k<NewFatherList.get_num_elem();k++)
    delete_oldsons(El_Table,NodeTable,myid,NewFatherList.get(k));
  /*
  MPI_Barrier(MPI_COMM_WORLD);

  for(iproc=0;iproc<nump-1;iproc++)
    if(myid>iproc); MPI_Barrier(MPI_COMM_WORLD);

  printf("myid=%d unref after delete\n",myid); fflush(stdout);

  tic=time(NULL);
  toc=tic+2;
  do{
    tic=time(NULL);
  }while(tic<toc);

  for(iproc=1;iproc<nump;iproc++)
    if(myid<iproc); MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  */

  //assert(!IfMissingElem(El_Table, myid, time_step, 2));


  
  /*  printf("myid=%d after deleting\n",myid);  
  MPI_Barrier(MPI_COMM_WORLD);
  printf("myid=%d leaving unrefine\n",myid);
  */

  /* 
     char debugfilename[64];
     sprintf(debugfilename,"unref%02d%08d.debug",myid,time_step);
     FILE *fpdebug=fopen(debugfilename,"w");
     fprintf(fpdebug,"%d elements unrefined on process %d ====================================!!!!!!!!!!!\n", unrefined*4, myid);
     fclose(fpdebug);
  */
  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);

  for(i=0; i<El_Table->get_no_of_buckets(); i++) {
    HashEntryPtr currentPtr = *(buck+i);
    while(currentPtr) {
      Curr_El = (Element*) currentPtr->value;
      currentPtr=currentPtr->next;
      if(Curr_El->get_adapted_flag()>TOBEDELETED)
	Curr_El->calc_wet_dry_orient(El_Table);
    }
  }
	
  //printf("myid=%d exiting unrefine\n",myid);
  return;
}

int Element::find_brothers(HashTable* El_Table, HashTable* NodeTable, 
			   double target, int myid, MatProps* matprops_ptr,
			   void *NFL, void *OPU) 
{
  ElemPtrList* NewFatherList  =(ElemPtrList*) NFL;
  ElemPtrList* OtherProcUpdate=(ElemPtrList*) OPU;

  int i = 0, j;
  int unrefine_flag = 1;
  Element* bros[5];
  if(opposite_brother_flag == 0) {
    find_opposite_brother(El_Table);
    if(opposite_brother_flag == 0) 
      return 0;
  }
  while(i<4 && unrefine_flag == 1) {
    Element* EmTemp = (Element*) El_Table->lookup(&brothers[i][0]);
    if(EmTemp == NULL) //|| EmTemp->refined != 0)
      return 0;
    else if(EmTemp->adapted!=NOTRECADAPTED) //this should be sufficient
      return 0;
    bros[i+1] = EmTemp;
    if(bros[i+1]->get_myprocess()!=myid) return 0; //should not be necessary because of "adapted" check
    
    unrefine_flag = EmTemp->check_unrefinement(El_Table,target);
    i++;
  }

  if(unrefine_flag == 1) { // we want to unrefine this element...
    //first we create the father element
    //printf("==============================\n unrefining an element \n===============================\n");
    if(*(bros[1]->getfather()) == (unsigned) 1529353130)
      printf("creating father %u %u from %u %u\n", *(bros[1]->getfather()), *(bros[1]->getfather()+1), *(bros[1]->pass_key()), *(bros[1]->pass_key()+1));
    bros[0] = new Element((bros+1), NodeTable, El_Table, matprops_ptr);
    El_Table->add(bros[0]->pass_key(), bros[0]);
    assert(bros[0]);  // a copy of the parent should always be on the same process as the sons
    NewFatherList->add(bros[0]);

    for(int ineigh=0;ineigh<8;ineigh++)
      if((bros[0]->neigh_proc[ineigh]>=0)&&
	 (bros[0]->neigh_proc[ineigh]!=myid)) {
	OtherProcUpdate->add(bros[0]);
	break;
      }
  }
  
  return unrefine_flag;
}


int Element::check_unrefinement(HashTable* El_Table, double target) {
  int unrefine_flag = 1, i = 0;
  
  if(adapted!=NOTRECADAPTED)
    //This rules out NEWFATHERs, NEWSONs, BUFFERs, GHOSTs, TOBEDELETEDs, and OLDFATERs  
    //This is a redundant check but is is better to be safe than sorry
    return(0); 

  for(int ineigh=0;ineigh<8;ineigh++) {  
    if(((neigh_proc[i]!=myprocess)&& 
	(neigh_proc[i]>=0) &&
	(generation<=0)
	)|| 
       (neigh_gen[i]>generation)
      )  
      return(0); 
    i++;
  }

  return(1);
 }

 

void delete_oldsons(HashTable* El_Table, HashTable* NodeTable,
		    int myid, void *EmFather_in){
  int ison, isonneigh, ineigh, inode, nodeorder[9];
  Element *EmSon, *EmNeigh;
  Node* NdTemp;
  Element *EmFather=(Element *) EmFather_in;
  assert(EmFather);

  /*
  unsigned elemdebugkey[2]={ 695892755,2973438897};
  if(compare_key(EmFather->pass_key(),elemdebugkey)){
    printf("myid=%d unrefine():deleteoldsons():EmFather={%10u,%10u}\n",
	   myid,elemdebugkey[0],elemdebugkey[1]);
    ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey,stdout);
    printf("it's 4 sons are\n");
    for(ison=0;ison<4;ison++) {
      printf("***ison=%d {%10u,%10u}\n",ison,
	     EmFather->son[ison][0],EmFather->son[ison][1]);
      ElemBackgroundCheck(El_Table,NodeTable,EmFather->son[ison],stdout);
    }
  }
  */

  EmFather->refined=0;

  for(inode=4;inode<8;inode++) {
    NdTemp= (Node *) NodeTable->lookup(EmFather->node_key[inode]);
    if(!NdTemp) {
      ElemBackgroundCheck(El_Table,NodeTable,EmFather->key,stdout);
      printf("inode=%d is missing\n",inode); fflush(stdout);
      NodeBackgroundCheck(El_Table,NodeTable,EmFather->node_key[inode],stdout) ; 
    }
    assert(NdTemp);
    nodeorder[inode]=NdTemp->order;
  }
  inode=8;
  NdTemp= (Node *) NodeTable->lookup(EmFather->key);
  assert(NdTemp);
  nodeorder[inode]=NdTemp->order;

  for(ison=0;ison<4;ison++){
    EmSon=(Element *) El_Table->lookup(EmFather->son[ison]);
    if(EmSon==NULL){
      printf("delete_oldsons() null son, ison=%d son={%u,%u}\n",ison,
	     EmFather->son[ison][0],EmFather->son[ison][1]);
      int yada;
      scanf("%d",&yada);
    }
    assert(EmSon);
    assert(EmSon->adapted==OLDSON);
    EmSon->adapted=TOBEDELETED;

    //delete son's bubble nodes
    NdTemp= (Node *) NodeTable->lookup(EmFather->son[ison]);
    assert(NdTemp);
    if(NdTemp->order>nodeorder[8])
      nodeorder[8]=NdTemp->order;
    NodeTable->remove(NdTemp->key, 0,stdout,myid,7);
    delete NdTemp;

    //delete son to son edge nodes
    inode=(ison+1)%4+4;
    NdTemp= (Node *) NodeTable->lookup(EmSon->node_key[inode]);
    assert(NdTemp);
    NodeTable->remove(NdTemp->key, 0,stdout,myid,8);
    delete NdTemp;
      
    //check 2 other edge nodes per son and delete if necessary
    //the first
    isonneigh=ison;
    inode=isonneigh+4;
    ineigh=ison;

    //EmNeigh=(Element *) El_Table->lookup(EmFather->neighbor[ineigh]);
    if((EmFather->neigh_gen[ineigh]==EmFather->generation)||
       (EmFather->neigh_proc[ineigh]==-1)
       ) { 
      NdTemp= (Node *) NodeTable->lookup(EmSon->node_key[inode]);
      if(NdTemp) {
	if(NdTemp->order>nodeorder[inode]) 
	  nodeorder[inode]=NdTemp->order;
	NodeTable->remove(NdTemp->key, 0,stdout,myid,9);
	delete NdTemp;   
      }   
    }
    else if(EmFather->neigh_proc[ineigh]==myid) {
      assert(EmFather->neigh_gen[ineigh]==EmFather->generation+1);
      NdTemp= (Node *) NodeTable->lookup(EmSon->node_key[inode]);
      assert(NdTemp);
      NdTemp->info=S_S_CON;
    }

    //the second
    isonneigh=(ison+3)%4;
    inode=isonneigh+4;
    ineigh=inode;

    //EmNeigh=(Element *) El_Table->lookup(EmFather->neighbor[ineigh]);
    if((EmFather->neigh_gen[ineigh]==EmFather->generation)||
       (EmFather->neigh_proc[ineigh%4]==-1)
       ) { 
      NdTemp= (Node *) NodeTable->lookup(EmSon->node_key[inode]);
      if(NdTemp) {
	if(NdTemp->order>nodeorder[inode]) 
	  nodeorder[inode]=NdTemp->order;

	NodeTable->remove(NdTemp->key, 0,stdout,myid,10);
	delete NdTemp;
      }

      NdTemp= (Node *) NodeTable->lookup(EmFather->node_key[inode]);
      assert(NdTemp);
      NdTemp->info=SIDE;
    }
    else if(EmFather->neigh_proc[ineigh]==myid) {
      NdTemp= (Node *) NodeTable->lookup(EmSon->node_key[inode]);
      assert(NdTemp);
      NdTemp->info=S_S_CON;

      NdTemp= (Node *) NodeTable->lookup(EmFather->node_key[inode]);
      assert(NdTemp);
      NdTemp->info=S_C_CON;      
      nodeorder[inode]=NdTemp->order;
    }

    //Now delete this oldson Element
    El_Table->remove(EmSon->key, 1,stdout,myid,11);
    EmSon->void_bcptr();
    delete EmSon;
  }

  for(inode=4;inode<8;inode++) {
    NdTemp= (Node *) NodeTable->lookup(EmFather->node_key[inode]);
    NdTemp->order=nodeorder[inode];
  }
  inode=8;

  NdTemp= (Node *) NodeTable->lookup(EmFather->key);
  NdTemp->order=nodeorder[inode];
  NdTemp->info=BUBBLE;

  return;
}

void Element::change_neigh_info(unsigned* fth_key, unsigned* ng_key, int neworder, 
				int fth_gen, int fth_proc) {
  int i,j, which_side = -1, same;
  i = 0;
  
  while(i<8 && which_side == -1) {
    if(neigh_proc[i] >= 0) {
      j = 0;
      same = 1;
      while(j<KEYLENGTH && same == 1) {
	if(neighbor[i][j] != ng_key[j]) 
	  same = 0;
	j++;
      }
      if(same == 1) {
	if(i<4)
	  which_side = i;
	else
	  which_side = i -4;
      }
    }
    i++;
  }
  if(!(which_side>=0)){
    int yada;
    printf("which_side=%d \n",which_side);
    scanf("%d",&yada);
  }

  assert(which_side >= 0);
  order[which_side] = neworder;

  for(i=0;i<KEYLENGTH;i++) {
    neighbor[which_side][i] = fth_key[i];
    neighbor[which_side+4][i] = fth_key[i];
  }

  neigh_gen[which_side  ] = fth_gen;
  neigh_gen[which_side+4] = fth_gen;

  neigh_proc[which_side+4] = -2; 

  return;
} 

//make this an element friend function
void unrefine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int myid, void* NFL){

  ElemPtrList* NewFatherList=(ElemPtrList*) NFL;

  int iupdate, ineigh, isonA, isonB, ineighme, ikey;
  Element *EmNeigh, *EmFather;

  //loop through the NEWFATHER elements
  for(iupdate=0;iupdate<NewFatherList->get_num_elem();iupdate++) {

    //I'm a NEWFATHER I'm going to update my neighbors with 
    //my information and if he's a NEWFATHER too I'm going 
    //to and update my information about him
    EmFather=NewFatherList->get(iupdate);
    assert(EmFather); //Help I've been abducted call the FBI!!!
    /*
    unsigned elemdebugkey[2]={ 695892755,2973438897};
    if(compare_key(EmFather->pass_key(),elemdebugkey)){
      printf("myid=%d unrefine():refine_neigh_update: EmFather={%10u,%10u}\n",
	     myid,elemdebugkey[0],elemdebugkey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey,stdout);
      printf("it's 4 sons are\n");
      for(int ison=0;ison<4;ison++) {
	printf("***ison=%d {%10u,%10u} \n",ison,
	       EmFather->getson()+ison*KEYLENGTH+0,
	       EmFather->getson()+ison*KEYLENGTH+1);
	ElemBackgroundCheck(El_Table,NodeTable,EmFather->getson()+ison*KEYLENGTH,stdout);
      }
    }
    */

    for(ineigh=0;ineigh<8;ineigh++)
      if(EmFather->neigh_proc[ineigh]==myid) {
	//only update the information of on processor neighbors in 
	//this function.

	EmNeigh=(Element*) El_Table->lookup(EmFather->neighbor[ineigh]);
	assert(EmNeigh); //Somebody has abducted my neighbor call the FBI!!!

	if(EmNeigh->adapted!=NEWFATHER) {
	  //If I knew my neighbor was a NEWFATHER that means
	  //my neighbor already updated his and my neighbor 
	  //information about each other

	  if(EmNeigh->adapted==OLDSON) {
	    //I am introduced to a NEWFATHER neighbor by his OLDSON
	    EmNeigh=(Element*) El_Table->lookup(EmNeigh->father);
	    assert(EmNeigh); //Somebody has abducted my neighbor call the FBI!!!
	  }

	  //One of my OLDSONs will introduce my neighbor to me
	  isonA=ineigh%4;
	  isonB=(isonA+1)%4;

	  for(ineighme=0;ineighme<4;ineighme++)
	    if(compare_key(EmNeigh->neighbor[ineighme],EmFather->son[isonA])||
	       compare_key(EmNeigh->neighbor[ineighme],EmFather->son[isonB]))
	      break;
	  if(!(ineighme<4))
	    printf("DANGER\n");
	  assert(ineighme<4);
	  
	  //Give my neighbor my contact information
	  EmNeigh->neigh_gen[ineighme]=
	    EmNeigh->neigh_gen[ineighme+4]=
	    EmFather->generation;

	  EmNeigh->neigh_proc[ineighme+4]=-2;

	  if(EmNeigh->adapted==NEWFATHER) {
	    //if my neighbor is a NEWFATHER, I need to update my 
	    //information about him too
	    EmFather->neigh_gen[ineigh]=
	      EmFather->neigh_gen[ineigh+4]=
	      EmNeigh->generation;

	    EmFather->neigh_proc[ineigh+4]=-2;
	  
	    //EmNeigh is inside this if() so only one loop through 
	    //KEYLENGTH is made, it saves on overhead
	    for(ikey=0;ikey<KEYLENGTH;ikey++) {
	      EmFather->neighbor[ineigh][ikey]=
		EmFather->neighbor[ineigh+4][ikey]=
		EmNeigh->key[ikey];
	      EmNeigh->neighbor[ineighme][ikey]=
		EmNeigh->neighbor[ineighme+4][ikey]=
		EmFather->key[ikey];
	    }
	  }else
	    //EmNeigh is inside this if() so only one loop through 
	    //KEYLENGTH is made, it saves on overhead
	    for(ikey=0;ikey<KEYLENGTH;ikey++) 
	      EmNeigh->neighbor[ineighme][ikey]=
		EmNeigh->neighbor[ineighme+4][ikey]=
		EmFather->key[ikey];
	}
      }
  }

  return;
}

//make this a Node and Element friend fucntion
void  unrefine_interp_neigh_update(HashTable* El_Table, HashTable* NodeTable,
				   int nump, int myid, void* OPU)
{
  ElemPtrList* OtherProcUpdate=(ElemPtrList*) OPU;
  
  if(nump<2) return;
  
  int *num_send=CAllocI1(nump), *num_recv=CAllocI1(nump), *isend=CAllocI1(nump); 
  int ierr, iopu, iproc, ineigh, ison, ikey, neigh_proc; //iopu stands for i other processor update
  int send_tag = 061116; //2006 November 16, the day I coded this function.
  MPI_Request* request = new MPI_Request[2*nump];
  Element* EmFather;
  Node *NdTemp;

  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 1.0\n",myid); fflush(stdout);

  
  for(iproc=0; iproc<nump; iproc++) 
    num_send[iproc]=num_recv[iproc]=isend[iproc]=0;
  
  
  for(iopu=0;iopu<OtherProcUpdate->get_num_elem();iopu++) {
    EmFather=OtherProcUpdate->get(iopu);
    assert(EmFather);
    assert(EmFather->myprocess==myid);
    for(ineigh=0;ineigh<8;ineigh++) {
      neigh_proc=EmFather->neigh_proc[ineigh];
      if((neigh_proc>=0)&&
	 (neigh_proc!=myid)
	 )
	num_send[neigh_proc]++;
    }
  }
  num_send[myid]=0; //better paranoid then dead;

  //send to each processor the number of his elements he has to update because of me
  //receive from each processor the number of elements I have to update because of him
  MPI_Alltoall(num_send,1, MPI_INT, num_recv, 1, MPI_INT, MPI_COMM_WORLD);

  if(myid==TARGET_PROC)
    for(iproc=0;iproc<nump;iproc++){
      printf("myid=%2d to/from proc %2d num_send=%6d num_recv=%6d\n",
	     myid,iproc,num_send[iproc],num_recv[iproc]);
      fflush(stdout);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if((myid!=TARGET_PROC)&&(TARGET_PROC<nump)&&(TARGET_PROC>0)) {
    printf("myid=%2d to/from proc %2d num_send=%6d num_recv=%6d\n",
	   myid,TARGET_PROC,num_send[TARGET_PROC],num_recv[TARGET_PROC]);
      fflush(stdout);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  /*
  for(iproc=0;iproc<nump;iproc++)
    printf("myid=%d neighproc=%d numsend=%d numrecv=%d\n",
	   myid,iproc,num_send[iproc],num_recv[iproc]);
  exit(0);
  */
  num_recv[myid]=0; //better paranoid then dead;

  int max_num_send=0,max_num_recv=0;

  for(iproc=0;iproc<nump;iproc++) {
    if(num_send[iproc]>max_num_send) max_num_send=num_send[iproc];
    if(num_recv[iproc]>max_num_recv) max_num_recv=num_recv[iproc];
  }

  unsigned **send,**recv;

  if(max_num_send>0) send=CAllocU2(nump,4*KEYLENGTH*max_num_send);
  if(max_num_recv>0) recv=CAllocU2(nump,4*KEYLENGTH*max_num_recv);

  for(iproc=0;iproc<nump;iproc++)
    if((iproc!=myid)&&(num_recv[iproc]>0))
      ierr=MPI_Irecv((void *) recv[iproc], 4*KEYLENGTH*num_recv[iproc], 
		     MPI_UNSIGNED, iproc, send_tag+iproc, MPI_COMM_WORLD, 
		     (request+nump+iproc));


  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 2.0\n",myid); fflush(stdout);

  for(iopu=0;iopu<OtherProcUpdate->get_num_elem();iopu++) {
    EmFather=OtherProcUpdate->get(iopu);
    assert(EmFather);
    for(ineigh=0;ineigh<8;ineigh++) {
      neigh_proc=EmFather->neigh_proc[ineigh];
      if((neigh_proc>=0)&&
	 (neigh_proc!=myid)
	 ) {

	switch(ineigh){
	case 0:
	case 7:
	  ison=0;
	  break;
	case 1:
	case 4:
	  ison=1;
	  break;
	case 2:
	case 5:
	  ison=2;
	  break;
	case 3:
	case 6:
	  ison=3;
	  break;
	default:
	  assert(0);	  
	}


	NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[ineigh%4+4]);
	assert(NdTemp);

	if(EmFather->neigh_gen[ineigh]-1==EmFather->generation)
	  NdTemp->info=S_C_CON;	

	for(ikey=0;ikey<KEYLENGTH;ikey++) {
	  //The element I want my neighbor to update
	  send[neigh_proc][(4*isend[neigh_proc]+0)*KEYLENGTH+ikey]=
	    EmFather->neighbor[ineigh][ikey];
	  //the OLDSON who will introduce his NEWFATHER to his
	  //neighbor on another processor
	  send[neigh_proc][(4*isend[neigh_proc]+1)*KEYLENGTH+ikey]=
	    EmFather->son[ison][ikey];
	  //the NEWFATHER element
	  send[neigh_proc][(4*isend[neigh_proc]+2)*KEYLENGTH+ikey]=
	    EmFather->key[ikey];	  
	  send[neigh_proc][(4*isend[neigh_proc]+3)*KEYLENGTH+ikey]=
	    NdTemp->key[ikey];
	}

	isend[neigh_proc]++;
      }
    }
  }
  for(iproc=0;iproc<nump;iproc++)
    assert(isend[iproc]==num_send[iproc]);

  CDeAllocI1(isend);

  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 3.0\n",myid); fflush(stdout);

  for(iproc=0;iproc<nump;iproc++)
    if((iproc!=myid)&&(num_send[iproc]>0))
      ierr=MPI_Isend((void *) send[iproc], 4*KEYLENGTH*num_send[iproc], 
		     MPI_UNSIGNED, iproc, send_tag+myid, MPI_COMM_WORLD, 
		     (request+iproc));

  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 4.0\n",myid); fflush(stdout);


  int NumProcsNotRecvd, ifrecvd, nodeorder;
  MPI_Status status;
  Element *EmTemp;
  if(max_num_recv>0) do{

    for(iproc=0;iproc<nump;iproc++)
      if((iproc!=myid)&&(num_recv[iproc]>0)) { 
	if(myid==TARGET_PROC) printf("myid=%d iproc=%2d unref_interp_neigh_update 5.0\n",myid,iproc); fflush(stdout);


	//only check processors I haven't already handled
	ifrecvd=0;
	//printf("myid=%d before Test",myid); fflush(stdout);
	MPI_Test(request+nump+iproc,&ifrecvd,&status);
	//printf("myid=%d after Test",myid); fflush(stdout);
	if(ifrecvd) {
	  if(myid==TARGET_PROC) printf("myid=%d iproc=%2d unref_interp_neigh_update 6.0\n",myid,iproc); fflush(stdout);

	  //I have just received new data from a neighboring processor
	  //I need to update some elements on my interprocessor 
	  //boundary

	  nodeorder=0; //initialize to zero

	  for(iopu=0;iopu<num_recv[iproc];iopu++) {
	    //one by one check my Element's that my neighbor processor
	    //says I need to update

	    //Hi I'm EmTemp
	    EmTemp=(Element*) 
	      El_Table->lookup(&(recv[iproc][(4*iopu+0)*KEYLENGTH]));

	    //my old neighbor will introduce his NEWFATHER to me
	    for(ineigh=0;ineigh<8;ineigh++) 
	      if(compare_key(EmTemp->neighbor[ineigh], 
			     &(recv[iproc][(4*iopu+1)*KEYLENGTH])))
		break;
	    assert(ineigh<8); //I don't know this Element pretending to be
	    //my neighbor, I'm calling the FBI to report a suspicious 
	    //person... err.. suspicious Element
	    int ineighmod4=ineigh%4;

	    if(EmTemp->adapted==OLDSON) {
	      //I'm moving out too so I'll introduce my NEWFATHER to my 
	      //neighbor's NEWFATHER who my neighbor just introduced to me

	      //we can use my father's key directly instead of having to 
	      //use get_father() because we know that my father's key 
	      //was assigned to me in unrefine_elements()
	      EmFather=(Element*) El_Table->lookup(EmTemp->father); 
	      assert(EmFather);
	      
	      //which means my father will only have 1 neighbor on that side
	      EmFather->neigh_proc[ineighmod4+4]=-2;

	      for(ikey=0;ikey<KEYLENGTH;ikey++)
		EmFather->neighbor[ineighmod4][ikey]=
		  EmFather->neighbor[ineighmod4+4][ikey]=
		  recv[iproc][(4*iopu+2)*KEYLENGTH+ikey];	      	      

	      //I know my neighbor on the other processor is the same 
	      //generation as me because we were both unrefined and only 
	      //one of us would have been able to unrefine if we were of
	      //different generations, that means our NEWFATHERs are the
	      //same generation as each other too
	      EmFather->neigh_gen[ineighmod4]=
		EmFather->neigh_gen[ineighmod4+4]=
		EmFather->generation;


	      //all the other neighbor information remains the same but 
	      //must delete 2 nodes I'll do 1 now and the other one was 
	      //just done previous or will be deleted next
	      NdTemp=(Node *) NodeTable->lookup(EmTemp->node_key[ineighmod4+4]);
	      if(NdTemp) {
		if(NdTemp->order>nodeorder) nodeorder=NdTemp->order;
		NodeTable->remove(NdTemp->key, 0,stdout,myid,12);
		delete NdTemp;
	      }
	      NdTemp=(Node *) NodeTable->lookup(EmFather->node_key[ineighmod4+4]);
	      NdTemp->order=nodeorder;
	      NdTemp->info=SIDE;
	      	      
	    }
	    else if(EmTemp->adapted>=NOTRECADAPTED) {
	      nodeorder=0;
	      
	      //my neighbor on the other processor was of either my 
	      //generation or one generation younger/higher than me, 
	      //that makes
	      //his father either one generation lower than me or of
	      //my generation. Either way his FATHER will be the only 
	      //neighbor I have on that side
	      
	      EmTemp->neigh_proc[ineighmod4+4]=-2;
	      
	      for(ikey=0;ikey<KEYLENGTH;ikey++)
		EmTemp->neighbor[ineighmod4+4][ikey]=
		  EmTemp->neighbor[ineighmod4][ikey]=
		  recv[iproc][(4*iopu+2)*KEYLENGTH+ikey];	      	      	      
	      
	      EmTemp->neigh_gen[ineighmod4+4]=
		EmTemp->neigh_gen[ineighmod4]=
		EmTemp->neigh_gen[ineighmod4]-1;
	      //(recv[iproc][(4*iopu+3)*KEYLENGTH+0]?-1:1)*
	      //(recv[iproc][(4*iopu+3)*KEYLENGTH+1]);
	      
	      //all the other neighbor information remains the same but may need 
	      //to delete some nodes and update side node
	      //titan doesn't use node order but for a continuous galerkin code 
	      //you may need to reset the node order here, which would be extra
	      //information to be communicated accross.  If you enforce the 
	      //minimum generation to be zero then you only need 1 unsigned 
	      //number to hold the generation or you could use a union of
	      //int and unsigned an the first bit of the unsigned would be the
	      //sign of the int so you wouldn't need 2 unsigneds to communicate
	      //the generation.
	      NdTemp=(Node *) NodeTable->lookup(EmTemp->node_key[ineighmod4+4]);
	      assert(NdTemp);
	      if(EmTemp->neigh_gen[ineighmod4]==EmTemp->generation)
		NdTemp->info=SIDE;
	      else{//my neighbor was my generation but his father moved in
		NdTemp->info=S_S_CON;

		NdTemp=(Node*) NodeTable->lookup(recv[iproc]+(4*iopu+3)*KEYLENGTH);
		if(!NdTemp){
		  ElemBackgroundCheck(El_Table,NodeTable,recv[iproc]+(4*iopu+3)*KEYLENGTH,stdout);
		  assert(NdTemp);
		}
		NdTemp->info=S_C_CON;
	      }
	    }
	    else
	      //EmTemp is missing... he may have been abducted
	      //Call the FBI to report a missing person! err... make that
	      //Call the FBI to report a missing Element!
	      assert(0);
	    
	  }
	  num_recv[iproc]=0;
	  if(myid==TARGET_PROC) printf("myid=%d iproc=%2d unref_interp_neigh_update 7.0\n",myid,iproc); fflush(stdout);

	}
      }
    if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 8.0\n",myid); fflush(stdout);

    NumProcsNotRecvd=0;
    for(iproc=0;iproc<nump;iproc++)
      if((iproc!=myid)&&(num_recv[iproc]>0)) NumProcsNotRecvd++;
    
  }while(NumProcsNotRecvd>0);
  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 9.0\n",myid); fflush(stdout);

  
  if(max_num_recv>0) CDeAllocU2(recv);
  CDeAllocI1(num_recv);

  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 10.0\n",myid); fflush(stdout);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(max_num_send>0) CDeAllocU2(send);
  CDeAllocI1(num_send);
  delete []request;

  if(myid==TARGET_PROC) printf("myid=%d unref_interp_neigh_update 11.0\n",myid); fflush(stdout);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return;
}
