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
 * $Id: updatenei.C 150 2007-06-27 20:28:42Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "../header/refined_neighbor_info.h"

#define TARGETPROC -1

#define AssertMeshNodeCheck


//! investigate an Element, question his "friends and family" about him.
void ElemBackgroundCheck(HashTable* El_Table, HashTable* NodeTable,
			 unsigned *debugkey, FILE *fp)
{  
  Element* EmDebug=(Element*) El_Table->lookup(debugkey);
  Element* EmTemp;
  ElemPtrList EmDebugNeigh(128);
  ElemPtrList EmDebugFather(16);
  Node* NdTemp;

  int iFather, ison, ineigh, inode;
  int uniqueneigh=8;


  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(int i=0; i<num_buck; i++)
    if(*(buck+i)){
      
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	
	EmTemp=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;     
	assert(EmTemp);
	
	for(ison=0;ison<4;ison++)	  
	  if(compare_key(EmTemp->son[ison],debugkey)){
	    EmDebugFather.add(EmTemp);
	    break;
	  }
	
	for(ineigh=0;ineigh<8;ineigh++)
	  if(compare_key(EmTemp->neighbor[ineigh],debugkey)){
	    EmDebugNeigh.add(EmTemp);
	    break;
	  }
	
      }
    }
  
  if(EmDebug) {
    if(EmDebugFather.get_num_elem()>0) {
      fprintf(fp,"EmDebug={%10u,%10u} has %d Father(s)\n",
	     EmDebug->key[0],
	     EmDebug->key[1],
	     EmDebugFather.get_num_elem());
      for(int iFather=0;iFather<EmDebugFather.get_num_elem();iFather++) {
	fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d adapted=%d which_son=%d iwetnode=%d Awet=%9.6f Swet=%9.6f drypoint={%9.6f,%9.6f}\n",iFather,
		*(EmDebugFather.get_key(iFather)+0),
		*(EmDebugFather.get_key(iFather)+1),
		EmDebugFather.get(iFather)->myprocess,
		EmDebugFather.get(iFather)->generation,
		EmDebugFather.get(iFather)->adapted,
		EmDebugFather.get(iFather)->which_son,
		EmDebugFather.get(iFather)->iwetnode,
		EmDebugFather.get(iFather)->Awet,
		EmDebugFather.get(iFather)->Swet,
		EmDebugFather.get(iFather)->drypoint[0],
		EmDebugFather.get(iFather)->drypoint[1]);
	uniqueneigh=8;
	for(ineigh=0;ineigh<8;ineigh++)
	  if(EmDebugFather.get(iFather)->neigh_proc[ineigh]<0)
	    uniqueneigh--;
	fprintf(fp,"   who has %d unique neighbors\n",uniqueneigh);
	for(ineigh=0;ineigh<8;ineigh++)
	  fprintf(fp,"   %d:   {%10u,%10u}  proc=%d gen=%d\n",ineigh,
		 EmDebugFather.get(iFather)->neighbor[ineigh][0],
		 EmDebugFather.get(iFather)->neighbor[ineigh][1],
		 EmDebugFather.get(iFather)->neigh_proc[ineigh],
		 EmDebugFather.get(iFather)->neigh_gen[ineigh]);	  
      }
    }

    uniqueneigh=8;
    for(ineigh=0;ineigh<8;ineigh++)
      if(EmDebug->neigh_proc[ineigh]<0)
	uniqueneigh--;
    fprintf(fp,"EmDebug={%10u,%10u} ",EmDebug->key[0],EmDebug->key[1]);
    fprintf(fp,"proc=%d ",EmDebug->myprocess);
    fprintf(fp,"gen=%d ",EmDebug->generation);
    fprintf(fp,"adapted=%d ",EmDebug->adapted);
    fprintf(fp,"which_son=%d ",EmDebug->which_son);
    fprintf(fp,"iwetnode=%d ",EmDebug->iwetnode);
    fprintf(fp,"Awet=%9.6f ",EmDebug->Awet);
    fprintf(fp,"Swet=%9.6f ",EmDebug->Swet);
    fprintf(fp,"drypoint={%9.6f,%9.6f} ",
	    EmDebug->drypoint[0],EmDebug->drypoint[1]);
    fprintf(fp,"has neighbors (%d are unique)\n",uniqueneigh);
    for(ineigh=0;ineigh<8;ineigh++)
      fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d\n",ineigh,
	      EmDebug->neighbor[ineigh][0],
	      EmDebug->neighbor[ineigh][1],
	      EmDebug->neigh_proc[ineigh],
	      EmDebug->neigh_gen[ineigh]);
    if(EmDebugNeigh.get_num_elem()>0) {
      fprintf(fp,"The following %d elements have {%10u,%10u} as a neighbor:\n",
	     EmDebugNeigh.get_num_elem(),debugkey[0],debugkey[1]);
      for(ineigh=0;ineigh<EmDebugNeigh.get_num_elem();ineigh++)
	fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d adapted=%d which_son=%d iwetnode=%d Awet=%9.6f Swet=%9.6f drypoint={%9.6f,%9.6f}\n",ineigh,
	       *(EmDebugNeigh.get_key(ineigh)+0),
		*(EmDebugNeigh.get_key(ineigh)+1),
		EmDebugNeigh.get(ineigh)->myprocess,
		EmDebugNeigh.get(ineigh)->generation,
		EmDebugNeigh.get(ineigh)->adapted,
		EmDebugNeigh.get(ineigh)->which_son,
		EmDebugNeigh.get(ineigh)->iwetnode,
		EmDebugNeigh.get(ineigh)->Awet,
		EmDebugNeigh.get(ineigh)->Swet,
		EmDebugNeigh.get(ineigh)->drypoint[0],
		EmDebugNeigh.get(ineigh)->drypoint[1]);
    }
    fprintf(fp,"his 8 non bubble nodes are:\n");
    for(inode=0;inode<8;inode++) {
      NdTemp=(Node*) NodeTable->lookup(EmDebug->node_key[inode]);
      //if(EmDebug->adapted>NOTRECADAPTED)
      //assert(NdTemp);
      if(NdTemp)
	fprintf(fp," %d:   {%10u,%10u}  info=%d\n",inode,
		NdTemp->key[0],NdTemp->key[1],NdTemp->info);
      else
	fprintf(fp," %d:   {%10u,%10u}  GHOST CELL MISSING THIS NODE\n",inode,
		EmDebug->node_key[inode][0],EmDebug->node_key[inode][1]);
    }
  }
  else if(EmDebugNeigh.get_num_elem()>0) {
    fprintf(fp,"Warning Background Check failed!...\nSuspsicious Element {%10u,%10u} is missing!...\nCall the FBI and put out an APB!\n",debugkey[0],debugkey[1]);
    //assert(0);
  }

  return;
}

void ElemBackgroundCheck2(HashTable *El_Table,HashTable *NodeTable,
			  void *EmDebug_in, FILE *fp)
{  
  if(!EmDebug_in) return;

  Element *EmDebug= (Element *) EmDebug_in;
  Element* EmTemp;
  ElemPtrList EmDebugNeigh(128);
  ElemPtrList EmDebugFather(16);
  Node* NdTemp;

  int iFather, ison, ineigh, inode;
  int uniqueneigh=8;


  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(int i=0; i<num_buck; i++)
    if(*(buck+i)){
      
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	
	EmTemp=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;     
	assert(EmTemp);
	
	for(ison=0;ison<4;ison++)	  
	  if(compare_key(EmTemp->son[ison],EmDebug->key)){
	    EmDebugFather.add(EmTemp);
	    break;
	  }
	
	for(ineigh=0;ineigh<8;ineigh++)
	  if(compare_key(EmTemp->neighbor[ineigh],EmDebug->key)){
	    EmDebugNeigh.add(EmTemp);
	    break;
	  }
	
      }
    }
  
  if(EmDebug) {
    if(EmDebugFather.get_num_elem()>0) {
      fprintf(fp,"EmDebug={%10u,%10u} has %d Father(s)\n",
	     EmDebug->key[0],
	     EmDebug->key[1],
	     EmDebugFather.get_num_elem());
      for(int iFather=0;iFather<EmDebugFather.get_num_elem();iFather++) {
	fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d adapted=%d which_son=%d iwetnode=%d Awet=%9.6f Swet=%9.6f drypoint={%9.6f,%9.6f}\n",iFather,
		*(EmDebugFather.get_key(iFather)+0),
		*(EmDebugFather.get_key(iFather)+1),
		EmDebugFather.get(iFather)->myprocess,
		EmDebugFather.get(iFather)->generation,
		EmDebugFather.get(iFather)->adapted,
		EmDebugFather.get(iFather)->which_son,
		EmDebugFather.get(iFather)->iwetnode,
		EmDebugFather.get(iFather)->Awet,
		EmDebugFather.get(iFather)->Swet,
		EmDebugFather.get(iFather)->drypoint[0],
		EmDebugFather.get(iFather)->drypoint[1]);
	uniqueneigh=8;
	for(ineigh=0;ineigh<8;ineigh++)
	  if(EmDebugFather.get(iFather)->neigh_proc[ineigh]<0)
	    uniqueneigh--;
	fprintf(fp,"   who has %d unique neighbors\n",uniqueneigh);
	for(ineigh=0;ineigh<8;ineigh++)
	  fprintf(fp,"   %d:   {%10u,%10u}  proc=%d gen=%d\n",ineigh,
		 EmDebugFather.get(iFather)->neighbor[ineigh][0],
		 EmDebugFather.get(iFather)->neighbor[ineigh][1],
		 EmDebugFather.get(iFather)->neigh_proc[ineigh],
		 EmDebugFather.get(iFather)->neigh_gen[ineigh]);	  
      }
    }

    uniqueneigh=8;
    for(ineigh=0;ineigh<8;ineigh++)
      if(EmDebug->neigh_proc[ineigh]<0)
	uniqueneigh--;
    fprintf(fp,"EmDebug={%10u,%10u} ",EmDebug->key[0],EmDebug->key[1]);
    fprintf(fp,"proc=%d ",EmDebug->myprocess);
    fprintf(fp,"gen=%d ",EmDebug->generation);
    fprintf(fp,"adapted=%d ",EmDebug->adapted);
    fprintf(fp,"which_son=%d ",EmDebug->which_son);
    fprintf(fp,"iwetnode=%d ",EmDebug->iwetnode);
    fprintf(fp,"Awet=%9.6f ",EmDebug->Awet);
    fprintf(fp,"Swet=%9.6f ",EmDebug->Swet);
    fprintf(fp,"drypoint={%9.6f,%9.6f} ",
	    EmDebug->drypoint[0],EmDebug->drypoint[1]);
    fprintf(fp,"has neighbors (%d are unique)\n",uniqueneigh);
    for(ineigh=0;ineigh<8;ineigh++)
      fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d\n",ineigh,
	      EmDebug->neighbor[ineigh][0],
	      EmDebug->neighbor[ineigh][1],
	      EmDebug->neigh_proc[ineigh],
	      EmDebug->neigh_gen[ineigh]);
    if(EmDebugNeigh.get_num_elem()>0) {
      fprintf(fp,"The following %d elements have {%10u,%10u} as a neighbor:\n",
	     EmDebugNeigh.get_num_elem(),EmDebug->key[0],EmDebug->key[1]);
      for(ineigh=0;ineigh<EmDebugNeigh.get_num_elem();ineigh++)
	fprintf(fp," %d:   {%10u,%10u}  proc=%d gen=%d adapted=%d which_son=%d iwetnode=%d Awet=%9.6f Swet=%9.6f drypoint={%9.6f,%9.6f}\n",ineigh,
	       *(EmDebugNeigh.get_key(ineigh)+0),
		*(EmDebugNeigh.get_key(ineigh)+1),
		EmDebugNeigh.get(ineigh)->myprocess,
		EmDebugNeigh.get(ineigh)->generation,
		EmDebugNeigh.get(ineigh)->adapted,
		EmDebugNeigh.get(ineigh)->which_son,
		EmDebugNeigh.get(ineigh)->iwetnode,
		EmDebugNeigh.get(ineigh)->Awet,
		EmDebugNeigh.get(ineigh)->Swet,
		EmDebugNeigh.get(ineigh)->drypoint[0],
		EmDebugNeigh.get(ineigh)->drypoint[1]);
    }
    fprintf(fp,"his 8 non bubble nodes are:\n");
    for(inode=0;inode<8;inode++) {
      NdTemp=(Node*) NodeTable->lookup(EmDebug->node_key[inode]);
      if(EmDebug->adapted>NOTRECADAPTED)
	assert(NdTemp);
      if(NdTemp)
	fprintf(fp," %d:   {%10u,%10u}  info=%d\n",inode,
		NdTemp->key[0],NdTemp->key[1],NdTemp->info);
      else
	fprintf(fp," %d:   {%10u,%10u}  GHOST CELL MISSING THIS NODE\n",inode,
		EmDebug->node_key[inode][0],EmDebug->node_key[inode][1]);
    }
  }

  return;
}



void NodeBackgroundCheck(HashTable *El_Table, HashTable* NodeTable,
			 unsigned *nodedbkey, FILE *fp){
  Node* NdDebug=(Node*) NodeTable->lookup(nodedbkey);
  Element* EmTemp;
  ElemPtrList ElemList(4);

  int inode;

  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(int i=0; i<num_buck; i++)
    if(*(buck+i)){
      
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	
	EmTemp=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;     
	assert(EmTemp);
	
	for(inode=0;inode<8;inode++)	  
	  if(compare_key(EmTemp->node_key[inode],nodedbkey)){
	    ElemList.add(EmTemp);
	    break;
	  }
	
	if(compare_key(EmTemp->key,nodedbkey))
	   ElemList.add(EmTemp);

	
	
      }
    }

  if(NdDebug||ElemList.get_num_elem())
    fprintf(fp,"========================================\n");

  if(NdDebug)
    fprintf(fp,"NdDebug={%10u,%10u} info=%d belongs to %d elements\n",
	    nodedbkey[0],nodedbkey[1],NdDebug->info,ElemList.get_num_elem());
  else
    fprintf(fp,"NdDebug={%10u,%10u} is missing; it belongs to %d elements\n",
	    nodedbkey[0],nodedbkey[1],ElemList.get_num_elem());  

  for(int ielem=0;ielem<ElemList.get_num_elem();ielem++) {
    fprintf(fp,"_%d_:   NdDebug={%10u,%10u} belongs to:\n",
	    ielem,nodedbkey[0],nodedbkey[1]);
    ElemBackgroundCheck(El_Table,NodeTable,ElemList.get_key(ielem),fp);
  }
   
  if(NdDebug||ElemList.get_num_elem())
    fprintf(fp,"========================================\n");
  return;

}


int ifNodeInfoChange(Node* NdDebug, int *NdDebugInfo) {

  if(NdDebug)
    if(NdDebug->getinfo()!=*NdDebugInfo) {
      *NdDebugInfo=NdDebug->getinfo();
      return 1;
    }

  return 0;
}

int ifCheckNode(HashTable* El_Table, HashTable* NodeTable, 
		int NdDebugInfo, int checkstate, unsigned *elemdbkey, 
		Element *EmFather, Element *EmSonA, Element *EmSonB,
		Element *EmNeighNew[4]) 
{
  int ifcheck=0;

  if((checkstate==0)||
     (NdDebugInfo==checkstate)) {
    if(compare_key(EmFather->pass_key(),elemdbkey)||
       compare_key(EmSonA->pass_key(),elemdbkey)||
       compare_key(EmSonB->pass_key(),elemdbkey)||
       compare_key(EmNeighNew[0]->pass_key(),elemdbkey)
     )
      return(1);
  
    for(int ineigh=1;ineigh<4;ineigh++)
      if(EmNeighNew[ineigh])
	if(compare_key(EmNeighNew[ineigh]->pass_key(),elemdbkey))
	  return(1);
  }


  return(0);
}


//this function checks for any and all possible mesh errors, 
//i.e. it checks if the mesh is legal, it says nothing about 
//the quality of a legal mesh, you must have ghost information 
//present before performing this check, WARNING THIS CHECK TAKES
//A LOT OF TIME, ONLY USE IT TO DEBUG.
void AssertMeshErrorFree(HashTable *El_Table, HashTable* NodeTable, 
			 int numprocs, int myid,double loc){
  
  return;

  char filename[256];
  sprintf(filename,"AssertMeshErrorFree%04d.debug",myid);
  FILE *fp=fopen(filename,"w");
  fprintf(fp,"location %g\n",loc);

  Node *NdTemp;
  Element *EmTemp;
  Element *EmNeigh[2];
  HashEntryPtr  *El_Table_bucket0,  El_Table_entry_ptr;
  HashEntryPtr *NodeTable_bucket0, NodeTable_entry_ptr;

  int ielembucket, iside, inode, ineigh, ineighp4, ineighme, i;
  int inodebucket, El_Table_num_buck, NodeTable_num_buck;
  

  assert(myid>=0);
  assert(numprocs>=1);
  assert(myid<numprocs);

  NodeTable_num_buck=NodeTable->get_no_of_buckets();
  NodeTable_bucket0 =NodeTable->getbucketptr();
  for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
    NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
        
    while(NodeTable_entry_ptr) {
	
      NdTemp=(Node*)(NodeTable_entry_ptr->value);
      NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
      assert(NdTemp);
      NdTemp->num_assoc_elem=0;
    }//while(NodeTable_entry_ptr)
  }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++




  //check the element information including missing neighbors
  El_Table_num_buck=El_Table->get_no_of_buckets();
  El_Table_bucket0 =El_Table->getbucketptr();
  for(ielembucket=0; ielembucket<El_Table_num_buck; ielembucket++) {
    El_Table_entry_ptr = *(El_Table_bucket0+ielembucket);
        
    while(El_Table_entry_ptr) {
	
      EmTemp=(Element*)(El_Table_entry_ptr->value);
      El_Table_entry_ptr=El_Table_entry_ptr->next;     
      assert(EmTemp);
      assert(EmTemp->generation>=MIN_GENERATION);
      assert(EmTemp->generation<=REFINE_LEVEL);
      
      assert((EmTemp->get_refined_flag()>=0)||
	     (EmTemp->get_refined_flag()==GHOST));
      
      if(!((EmTemp->get_adapted_flag()>=OLDSON)&&
	   (EmTemp->get_adapted_flag()<=BUFFER)&&
	   (EmTemp->get_adapted_flag()!=-NEWBUFFER)
	   )
	 ){
	ElemBackgroundCheck(El_Table,NodeTable,EmTemp->key,fp);
	fclose(fp);
	assert(0);	
      }
      
      //^ is "exclusive or", not "exclusive or" is "if and only if"
      assert(!((EmTemp->get_refined_flag()>=1)^
	       ((EmTemp->get_adapted_flag()==TOBEDELETED)||
		(EmTemp->get_adapted_flag()==OLDFATHER)||
		(EmTemp->get_adapted_flag()==OLDSON)
		)
	       )
	     );

      //^ is "exclusive or", not "exclusive or" is "if and only if"      
      assert(!((EmTemp->get_refined_flag()==GHOST)^
	       ((EmTemp->get_adapted_flag()>=-BUFFER)&&
		(EmTemp->get_adapted_flag()<=-NOTRECADAPTED)
		)
	       )
	     );
      
      //^ is "exclusive or", not "exclusive or" is "if and only if"
      assert(!((EmTemp->myprocess!=myid)^
	       ((EmTemp->get_adapted_flag()>=-BUFFER)&&
		(EmTemp->get_adapted_flag()<=-NOTRECADAPTED)
		)
	       )
	     );
      
      if(EmTemp->myprocess!=myid)
	assert(numprocs>1);
      
      //^ is "exclusive or", not "exclusive or" is "if and only if"
      assert((EmTemp->get_refined_flag()==0)^
	     !((EmTemp->get_adapted_flag()>=NOTRECADAPTED)&&
	       (EmTemp->get_adapted_flag()<=BUFFER)
	       )
	     );

      if(EmTemp->adapted>=NOTRECADAPTED) {
	NdTemp=(Node*) NodeTable->lookup(EmTemp->key);
	if(EmTemp->adapted>=NOTRECADAPTED) {
	  assert(NdTemp);
	  //if(NdTemp)
	  NdTemp->num_assoc_elem++;}
      }

      double xmax=-HUGE_VAL, xmin=HUGE_VAL, ymax=-HUGE_VAL, ymin=HUGE_VAL;
      double *coord;
      for(int inode=0;inode<8;inode++){
	NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
	if(EmTemp->adapted>=NOTRECADAPTED) {
	  assert(NdTemp);
	  NdTemp->num_assoc_elem++;
	}
	if(NdTemp) {
	  coord=NdTemp->get_coord();
	  if(coord[0]>xmax) xmax=coord[0];
	  if(coord[0]<xmin) xmin=coord[0];
	  if(coord[1]>ymax) ymax=coord[1];
	  if(coord[1]<ymin) ymin=coord[1];	  
	}
      }
      
      if(EmTemp->get_adapted_flag()>=NOTRECADAPTED){


	double tol=(xmax-xmin+ymax-ymin)/2048.0;
	double xmean=0.5*(xmax+xmin);
	double ymean=0.5*(ymax+ymin);
	
	NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[0]);
	assert(NdTemp);
	coord=NdTemp->get_coord();
	assert((coord[0]-tol<=xmin)&&(xmin<=coord[0]+tol)&&
	       (coord[1]-tol<=ymin)&&(ymin<=coord[1]+tol));
	
	NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[2]);
	assert(NdTemp);
	coord=NdTemp->get_coord();
	assert((coord[0]-tol<=xmax)&&(xmax<=coord[0]+tol)&&
	       (coord[1]-tol<=ymax)&&(ymax<=coord[1]+tol));
	
	NdTemp=(Node*) NodeTable->lookup(EmTemp->key);
	assert(NdTemp);
	coord=NdTemp->get_coord();
	assert((coord[0]-tol<=xmean)&&(xmean<=coord[0]+tol)&&
	       (coord[1]-tol<=ymean)&&(ymean<=coord[1]+tol));

      for(iside=0;iside<4;iside++) {
	
	ineigh  =iside;
	ineighp4=iside+4;

	assert(EmTemp->neigh_proc[ineigh]>=-1);
	
	if(EmTemp->neigh_proc[ineigh]>=0) {
	  //if neighbor is not a map boundary

	  EmNeigh[0]=(Element *) 
	    El_Table->lookup(EmTemp->neighbor[ineigh]);
	  if(!EmNeigh[0]) {
	    //if(myid==TARGETPROC) {
	    fprintf(stdout,"neighbor %d is missing\n",ineigh);
	    ElemBackgroundCheck(El_Table,NodeTable,EmTemp->key,stdout);
	    //}
	    fclose(fp);
	    assert(EmNeigh[0]);
	  }
	  
	  EmNeigh[1]=(Element *) 
	    El_Table->lookup(EmTemp->neighbor[ineighp4]);
	  assert(EmNeigh[1]);
	  
	  //^ is "exclusive or", not "exclusive or" is "if and only if"
	  assert(!((EmTemp->neigh_proc[ineighp4]==-2)^ 
		   compare_key(EmTemp->neighbor[ineigh  ],
			       EmTemp->neighbor[ineighp4])
		   )
		 );

	  if(EmTemp->neigh_proc[ineighp4]!=-2)
	    assert(EmTemp->neigh_proc[ineighp4]==
		   EmNeigh[1]->myprocess);
	  
	  
	  //loop over the 2 neighbors _on_this_side_
	  for(i=0;i<2;i++) if(EmNeigh[i]->adapted>=NOTRECADAPTED){
	    

	    //which of my neighbor's neighbors am I
	    for(ineighme=0;ineighme<8;ineighme++)
	      if(compare_key(EmNeigh[i]->neighbor[ineighme],
			     EmTemp->key))
		  break;
	    
	    assert(ineighme<8);
	    assert((ineighme+2)%4==ineigh); //correct sides match up
	    if(!(EmTemp->neigh_gen[ineigh+4*i]==
		 EmNeigh[i]->generation)) {
	      fprintf(fp,"EmTemp={%10u,%10u}\n",
		      *(EmTemp->pass_key()+0),*(EmTemp->pass_key()+1));
	      ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),fp);
	      fprintf(fp,"EmNeigh={%10u,%10u} ineigh=%d\n",
		      *(EmNeigh[i]->pass_key()+0),*(EmNeigh[i]->pass_key()+1),
		      ineigh+4*i);
	      ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[i]->pass_key(),fp);
	      fclose(fp);

	      assert(EmTemp->neigh_gen[ineigh+4*i]==
		     EmNeigh[i]->generation);
	    }
	    if(!(EmNeigh[i]->neigh_gen[ineighme]==
		 EmTemp->generation)) {
	      fprintf(fp,"EmTemp={%10u,%10u}\n",
		      *(EmTemp->pass_key()+0),*(EmTemp->pass_key()+1));
	      ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),fp);
	      fprintf(fp,"EmNeigh={%10u,%10u} ineigh=%d\n",
		      *(EmNeigh[i]->pass_key()+0),*(EmNeigh[i]->pass_key()+1),
		      ineigh+4*i);
	      ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[i]->pass_key(),fp);
	      fclose(fp);
	      assert(EmNeigh[i]->neigh_gen[ineighme]==
		     EmTemp->generation);
	    }
	    
	    if(EmTemp->generation<=EmNeigh[i]->generation){
	      assert(ineighme<4);
	      assert(EmNeigh[i]->neigh_proc[ineighme+4]==-2);
	    }
	 

	    //check to make sure my neighbor shares the correct nodes 
	    //with me and those nodes are of the correct types (infos)
	    //also checks if the 1 irregularity rule has been broken	    
	    if(EmTemp->get_adapted_flag()>=NOTRECADAPTED) {
	      //ghost cells will be missing nodes so don't check them
	      
	      switch(EmTemp->generation-EmNeigh[i]->generation) {
	      case -1:
		//I'm one generation older than my neighbor
		switch(i){ 
		case 0: //neighbor is ineigh=iside
		  assert(compare_key(EmTemp->node_key[ineigh  ],
				     EmNeigh[0]->node_key[(ineighme+1)%4]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineigh ]);
		  assert(NdTemp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		  
		  assert(compare_key(EmTemp->node_key[ineighp4],
				     EmNeigh[0]->node_key[ineighme  ]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineighp4]);
		  assert(NdTemp);

		  if(!(NdTemp->info==S_C_CON)){
		    //printf("myid=%d\n",myid);
		    //if(myid==TARGETPROC)
		    NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
		    fclose(fp);
		      //MPI_Barrier(MPI_COMM_WORLD);
		    assert(NdTemp->info==S_C_CON);		
		  }
		  
		  break;
		case 1: //neighbor is ineighp4=iside+4
		  assert(compare_key(EmTemp->node_key[ineighp4],
				     EmNeigh[1]->node_key[(ineighme+1)%4]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineighp4]);
		  assert(NdTemp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		  
		  assert(compare_key(EmTemp->node_key[(ineigh+1)%4],
				     EmNeigh[1]->node_key[ineighme  ]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[(ineigh+1)%4]);
		  assert(NdTemp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		  
		  break;
		default:
		  assert(0);
		  break;
		}
		break;
	      case 0:
		//I'm the same generation as my neighbor
		if(!compare_key(EmTemp->node_key[ineigh  ],
				EmNeigh[i]->node_key[(ineighme+1)%4])) {
		  fprintf(fp,"iside=%d ineigh=%d ineighme=%d i=%d\n",
			 iside,ineigh,ineighme,i);
		  ElemBackgroundCheck(El_Table,NodeTable,EmTemp->key,fp);
		  fprintf(fp,"[i=0]\n");
		  ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[0]->key,fp);
		  fprintf(fp,"[i=1]\n");
		  ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[1]->key,fp);
		  //fprintf(fp,"stop me\n");
		  fclose(fp);
		  assert(compare_key(EmTemp->node_key[ineigh  ],
				     EmNeigh[i]->node_key[(ineighme+1)%4]));
		}
		NdTemp=(Node*) 
		  NodeTable->lookup(EmTemp->node_key[ineigh]);
		assert(NdTemp);
		if(!((NdTemp->info==CORNER)||
		       (NdTemp->info==S_C_CON)
		     )
		   ){
		  //if(myid==TARGETPROC)
		  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
		  fclose(fp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		}

		assert(compare_key(EmTemp->node_key[ineighp4],
				   EmNeigh[i]->node_key[ineighme+4]));
		NdTemp=(Node*) 
		  NodeTable->lookup(EmTemp->node_key[ineighp4]);
		assert(NdTemp);
		assert(NdTemp->info==SIDE);
		
		assert(compare_key(EmTemp->node_key[(ineigh+1)%4],
				   EmNeigh[i]->node_key[ineighme  ]));
		NdTemp=(Node*) 
		  NodeTable->lookup(EmTemp->node_key[(ineigh+1)%4]);
		assert(NdTemp);

		if(!((NdTemp->info==S_C_CON)||(NdTemp->info==CORNER))){
		  //printf("myid=%d\n",myid);
		  //if(myid==TARGETPROC)
		  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
		    //MPI_Barrier(MPI_COMM_WORLD);  
		    //printf("stop me\n");
		  fclose(fp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));	      
		}
		break;
	      case 1:
		//I'm one generation younger than my neighbor
		switch(ineighme/4){
		case 0:
		  assert(compare_key(EmTemp->node_key[ineigh  ],
				     EmNeigh[i]->node_key[ineighme+4]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineigh]);
		  assert(NdTemp);
		  if(!(NdTemp->info==S_C_CON)){
		    NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
		    fclose(fp);
		    assert(NdTemp->info==S_C_CON);
		  }
		  
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineighp4]);
		  assert(NdTemp);
		  assert(NdTemp->info==S_S_CON);
		  
		  assert(compare_key(EmTemp->node_key[(ineigh+1)%4],
				     EmNeigh[i]->node_key[ineighme  ]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[(ineigh+1)%4]);
		  assert(NdTemp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		  
		  break;
		case 1:
		  if(!compare_key(EmTemp->node_key[ineigh],
				  EmNeigh[i]->node_key[(ineighme+1)%4])){
		    fprintf(fp,"iside=%d ineigh=%d ineighme=%d i=%d\n",
			   iside,ineigh,ineighme,i);
		    ElemBackgroundCheck(El_Table,NodeTable,EmTemp->key,fp);
		    fprintf(fp,"[i=0]\n");
		    ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[0]->key,fp);
		    fprintf(fp,"[i=1]\n");
		    ElemBackgroundCheck(El_Table,NodeTable,EmNeigh[1]->key,fp);
		    fclose(fp);
		    
		    assert(compare_key(EmTemp->node_key[ineigh],
				       EmNeigh[i]->node_key[(ineighme+1)%4]));
		  }
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineigh]);
		  assert(NdTemp);
		  assert((NdTemp->info==CORNER)||
			 (NdTemp->info==S_C_CON));
		  
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[ineighp4]);
		  assert(NdTemp);
		  assert(NdTemp->info==S_S_CON);
		  
		  assert(compare_key(EmTemp->node_key[(ineigh+1)%4],
				     EmNeigh[i]->node_key[ineighme  ]));
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[(ineigh+1)%4]);
		  assert(NdTemp);
		  if(!(NdTemp->info==S_C_CON)){
		    NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
		    fclose(fp);
		    assert(NdTemp->info==S_C_CON);
		  }

		  
		  
		  break;
		default:
		  //for(i=0; i<2; i++) loop has failed
		  assert(0);
		  break;
		}
		break;		
	      default:
		//My neighbor and I are more than 1 generation appart
		//the one irregularity rule has been broken
		assert(0);
		break;
	      }
	    }//if(EmTemp->get_adapted_flag()>=NOTRECADAPTED)
	    else 
	      assert(abs(EmTemp->generation-EmNeigh[i]->generation)<=1);
	  }//for(int i=0;i<2;i++)
	}//if(EmTemp->neigh_proc[ineigh]>=-1);
      }// for(iside=0;iside<4;iside++)
      }//if(EmTemp->get_adapted_flag()>=NOTRECADAPTED)
    }//while(El_Table_entry_ptr)
  }//for(ielembucket=0; ielembucket<El_Table_num_buck; ielembucket++) 

#ifdef AssertMeshNodeCheck
  //check for obvious extra non ghost elements by checking how 
  //many elements each node belongs to.  This won't catch extra 
  //elements with their own extra nodes but the element check for 
  //missing neighbors should catch those.
  for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
    NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
        
    while(NodeTable_entry_ptr) {
      NdTemp=(Node*)(NodeTable_entry_ptr->value);
      NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
      assert(NdTemp);

      //check to see if node is of an allowed type
      assert((NdTemp->info==CORNER)||
	     (NdTemp->info==SIDE)||
	     (NdTemp->info==S_S_CON)||
	     (NdTemp->info==BUBBLE)||
	     (NdTemp->info==S_C_CON)
	     );

      /*
	if(!(NdTemp->num_assoc_elem>=1)){
	NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	fclose(fp);
	assert(0);
      }
      */

      if(NdTemp->info==CORNER)
	if(!(NdTemp->num_assoc_elem<=4)){
	  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	  fclose(fp);
	  assert(0);
	}  

      if(NdTemp->info==SIDE)
	if(!(NdTemp->num_assoc_elem<=2)){
	  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	  fclose(fp);
	  assert(0);
	}  

      if(NdTemp->info==S_S_CON)
	if(!(NdTemp->num_assoc_elem<=2)){
	  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	  fclose(fp);
	  assert(0);
	}  


      if(NdTemp->info==BUBBLE)
	if(!(NdTemp->num_assoc_elem<=1)){
	  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	  fclose(fp);
	  assert(0);
	}  

      if(NdTemp->info==S_C_CON)
	if(!(NdTemp->num_assoc_elem<=3)){
	  NodeBackgroundCheck(El_Table,NodeTable,NdTemp->key,fp);
	  fclose(fp);
	  assert(0);
	}  

    }//while(NodeTable_entry_ptr)
  }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++)
#endif

  fclose(fp);

  return;
}

int IfNeighProcChange(HashTable* El_Table, HashTable* NodeTable, int myid, Element* EmDebug, Element* EmTemp) {
  if(myid==2)
    if(*(EmDebug->get_neigh_proc()+1)==2) {
      printf("changed neigh_proc element\n");
      ElemBackgroundCheck(El_Table,NodeTable,EmDebug->pass_key(),stdout);
      printf("element who mistakenly changed the neigh_proc\n");
      ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),stdout);
      return 1;
    }

  return 0;
}

void update_neighbor_info(HashTable* HT_Elem_Ptr, ElemPtrList* RefinedList, 
			  int myid, int numprocs, 
			  HashTable* HT_Node_Ptr, int h_count );

//#define NEWCODEDOESNOTWORKYET
#ifdef NEWCODEDOESNOTWORKYET
void refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, 
			 int numprocs, int myid, void* RL,
			 TimeProps* timeprops_ptr){

  update_neighbor_info(El_Table, (ElemPtrList*) RL, myid, numprocs, 
		       NodeTable, 0);
  return;
}

#else
void refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, 
			 int nump, int myid, void* RL,
			 TimeProps* timeprops_ptr)
{

  Element* EmFather;
  Element* EmSon[4];
  Element* EmNeighNew[4];
  Element* EmNeighOld[2];
  ElemPtrList* RefinedList=(ElemPtrList*) RL;
  Node*  NdTemp;
  int ifather, iside, ineigh, ineighp4, isonA, isonB;
  int ineighme, ineighmep4, ineighson, ikey, inewcase, inode;

  //interproc update only variables
  int *num_send, *num_recv, *isend;
  int iproc, ineighm4, ierr, send_tag=061201*2, neigh_proc;
  unsigned **send, **recv;
  MPI_Request* request = new MPI_Request[2*nump];

  //printf("myid=%d, nump=%d\n",myid,nump);

  //unsigned ElemDebugKey[2]={1897922560,         0};

  //unsigned nodedebugkey[2]={2115132074,2863311530};
  //unsigned NodeDebugKey[2]={1895825408,         0};
  //int NodeDebugInfo=-99999;
  //Node* NodeDebug=(Node*) NodeTable->lookup(NodeDebugKey);
  FILE *fpbg=stdout;
  char fname2[256];
  sprintf(fname2,"refine_neigh_update%04d.debug",myid);
  FILE *fpdb2;

  unsigned ElemDebugKey[2]={ 695804849, 991146299};
  Element* ElemDebug=(Element*) El_Table->lookup(ElemDebugKey);
  unsigned ElemDebugNeighKey[2]={ 695876266,2863311530};
  Element* ElemDebugNeigh=(Element*) El_Table->lookup(ElemDebugNeighKey);
  unsigned ElemDebugFatherKey[2]={0,0};
  Element* ElemDebugFather=NULL;
  unsigned ElemDebugFatherNeighKey[2]={0,0};
  Element* ElemDebugFatherNeigh=NULL;
  unsigned ElemDebugKidnapperKey[2]={ 351592920,3079262749};
  Element* ElemDebugKidnapper=(Element*) El_Table->lookup(ElemDebugKidnapperKey);
  
  if(0&&myid==2) {
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d iter=%d\n",myid,timeprops_ptr->iter);
    if(timeprops_ptr->iter==19) {
      fprintf(fpdb2,"**********************************************************\n");
      if(ElemDebug) {
	fprintf(fpdb2,"ElemDebug={%10u,%10u}============================\n",
		ElemDebugKey[0],ElemDebugKey[1]);
	ElemBackgroundCheck(El_Table,NodeTable,ElemDebugKey,fpdb2);
	ElemDebugFatherKey[0]=ElemDebug->father[0];
	ElemDebugFatherKey[1]=ElemDebug->father[1];
	ElemDebugFather=(Element*) El_Table->lookup(ElemDebugFatherKey);
	
      }
      if(ElemDebugNeigh) {
	fprintf(fpdb2,"ElemDebugNeigh={%10u,%10u}=======================\n",
		ElemDebugNeighKey[0],ElemDebugNeighKey[1]);
	ElemBackgroundCheck(El_Table,NodeTable,ElemDebugNeighKey,fpdb2);      
      }
      if(ElemDebugFather) {
	fprintf(fpdb2,"ElemDebugFather={%10u,%10u}======================\n",
		ElemDebugFatherKey[0],ElemDebugFatherKey[1]);
	ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherKey,fpdb2);
	ineigh=-1;
	if(ElemDebug->which_son==3) ineigh=3;
	else if(ElemDebug->which_son==0) ineigh=7;
	if(ineigh!=-1) {
	  ElemDebugFatherNeighKey[0]=ElemDebug->neighbor[ineigh][0];
	  ElemDebugFatherNeighKey[1]=ElemDebug->neighbor[ineigh][1];
	  ElemDebugFatherNeigh=(Element*) 
	    El_Table->lookup(ElemDebugFatherNeighKey);
	}
	fprintf(fpdb2,"Father ineigh=%d\n",ineigh);
      }
      if(ElemDebugFatherNeigh) {
	fprintf(fpdb2,"ElemDebugFatherNeigh={%10u,%10u} ineigh=%d========\n",
		ElemDebugFatherNeighKey[0],ElemDebugFatherNeighKey[1],ineigh);
	ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherNeighKey,fpdb2);
      }
      fprintf(fpdb2,"**********************************************************\n");
      fclose(fpdb2);
    }
  }

  //printf("myid=%d nump=%d\n",myid,nump);


  /*************************************************************/
  /* so I won't have to waste time waiting for neighbor update */
  /* information from other processors later, we all send the  */
  /* information first, i.e. before we do on processor updates */
  /* so the other processor's information will already be here */
  /* when I'm ready for it                                     */
  /*************************************************************/
  if(nump>1) {
    int ineighother;  //only need this variable inside this if statement 
    //so declare here rather than up to and it will disapear when the loop
    //exits

    num_send=CAllocI1(nump);
    num_recv=CAllocI1(nump);
    isend   =CAllocI1(nump); 

    for(iproc=0;iproc<nump;iproc++)
      num_send[iproc]=num_recv[iproc]=isend[iproc]=0;

    /*
    for(iproc=0;iproc<nump;iproc++)
      printf("A myid=%d iproc=%d num_send=%d, num_recv=%d\n",
	     myid,iproc,num_send[iproc],num_recv[iproc]);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    */

    for(ifather=RefinedList->get_inewstart();
	ifather<RefinedList->get_num_elem();
	ifather++){
    
      EmFather=RefinedList->get(ifather); //Hello I'm the OLDFATHER
      assert(EmFather); //Help I've been abducted call the FBI!!!
      assert(EmFather->adapted==OLDFATHER); //sanity check
      EmSon[0]=EmSon[1]=EmSon[2]=EmSon[3]=NULL;
      
      for(ineigh=0;ineigh<8;ineigh++){
	neigh_proc=EmFather->neigh_proc[ineigh];
	assert(neigh_proc<nump);
	ineighm4=ineigh%4;
	ineighp4=ineighm4+4;
	if(ineigh==ineighm4)
	  ineighother=ineighp4;
	else
	  ineighother=ineighm4;
	if((neigh_proc>=0)&&(neigh_proc!=myid)) {
	  num_send[neigh_proc]++;
	  switch(ineigh){
	  case 0:
	    isonA=0;
	    isonB=1;
	    break;
	  case 7:
	    isonA=0;
	    isonB=3;
	    break;
	  case 1:
	    isonA=1;
	    isonB=2;
	    break;
	  case 4:
	    isonA=1;
	    isonB=0;
	    break;
	  case 2:
	    isonA=2;
	    isonB=3;
	    break;
	  case 5:
	    isonA=2;
	    isonB=1;
	    break;
	  case 3:
	    isonA=3;
	    isonB=0;
	    break;
	  case 6:
	    isonA=3;
	    isonB=2;
	    break;
	  default:
	    assert(0);
	    break;
	  }//switch(ineigh)
	  
	  if(!EmSon[isonA]) {//to save a little work, only lookup this
	    //son if I don't already know him
	    EmSon[isonA]=(Element*) El_Table->lookup(EmFather->son[isonA]);
	    assert(EmSon[isonA]);
	  }

	  if(!EmSon[isonB]) {//to save a little work, only lookup "other"
	    //son if I don't already know him
	    EmSon[isonB]=(Element*) El_Table->lookup(EmFather->son[isonB]);
	    assert(EmSon[isonB]);
	  }

	  switch(EmFather->generation-EmFather->neigh_gen[ineigh]){
	  case -1:
	    for(ikey=0;ikey<KEYLENGTH;ikey++)
	      EmSon[  isonA]->neighbor[ineighm4][ikey]=
		EmSon[isonA]->neighbor[ineighp4][ikey]=
		EmFather->neighbor[ineigh][ikey];

	    EmSon[  isonA]->neigh_gen[ ineighm4]=
	      EmSon[isonA]->neigh_gen[ ineighp4]=
	      EmSon[isonA]->generation;

	    EmSon[  isonA]->neigh_proc[ineighm4]=neigh_proc;
	    EmSon[  isonA]->neigh_proc[ineighp4]=-2;

	    NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[ineighp4]);
	    assert(NdTemp);
	    NdTemp->info=SIDE;
	    //printf("node update yada 1\n");
	    
	    if(ineigh<4) {
	      //printf("node update yada 2\n");
	      NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[ineighp4]);
	      assert(NdTemp);
	      NdTemp->info=CORNER;
	    }

	    if(EmFather->neigh_proc[ineighother]==myid){
	      //one neighbor on this side belongs to another processor and
	      //has just been taken care of
	      //however the other neighbor on this side belongs to my
	      //processor and I need to update him now... the logic is simpler
	      //and cleaner to do it here than with the on processor update
	      //
	      //Since the neighbors on this side are already younger than me
	      //I know they could not have been adapted already.  That means
	      //all I need to do are update this neighbor and myself, in other
	      //words I know the other neighbor on this side can not possibly 
	      //be an OLDFATHER.

	      //EmNeighOld[0] is only being assigned here for clarity, i.e.
	      //to show the old neighbor is the same as the newneighbor
	      EmNeighNew[0]=EmNeighOld[0]=(Element*) 
		El_Table->lookup(EmFather->neighbor[ineighother]);

	      assert(EmNeighNew[0]);
	      
	      for(ineighme=0;ineighme<4;ineighme++)
		if(compare_key(EmFather->key,
			       EmNeighNew[0]->neighbor[ineighme]))
		  break;
	  
	      assert(ineighme<4);
	      ineighmep4=ineighme+4;

	      for(ikey=0;ikey<KEYLENGTH;ikey++) {
		EmSon[  isonB]->neighbor[ineighm4][ikey]=
		  EmSon[isonB]->neighbor[ineighp4][ikey]=
		  EmFather->neighbor[ ineighother][ikey];
		EmNeighNew[  0]->neighbor[ineighme  ][ikey]=
		  EmNeighNew[0]->neighbor[ineighmep4][ikey]=
		  EmSon[ isonB]->key[ikey];
	      }
	      
	      
	      EmSon[   isonB]->neigh_gen[ ineighm4]=
		EmSon[ isonB]->neigh_gen[ ineighp4]=
		EmNeighNew[0]->generation;
	      

	      EmNeighNew[  0]->neigh_gen[ineighme  ]=
		EmNeighNew[0]->neigh_gen[ineighmep4]=
		EmSon[ isonB]->generation;

	      EmSon[  isonB]->neigh_proc[ineighm4]=myid;
	      EmSon[  isonB]->neigh_proc[ineighp4]=-2;

	      NdTemp=(Node*) 
		NodeTable->lookup(EmSon[isonB]->node_key[ineighp4]);
	      assert(NdTemp);
	      NdTemp->info=SIDE;
	      //printf("node update yada 1\n");
	    
	      if(ineigh>=4) {
		//printf("node update yada 2\n");
		NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[ineighp4]);
		assert(NdTemp);
		NdTemp->info=CORNER;
	      }
	    }

	    break;
	  case  0: //old neighbor was same generation as OLDFATHER and is 
	    //one generation older than NEWSON's so, because of the "new 
	    //difference" in generation we know that OLDFATHER's 2 corner 
	    //nodes on this side are actually CORNER's and not S_C_CON's
	    if(ineigh<4) {//only do once, this if should not be necessary 
	      //but better safe than sorry

	      isonB=(isonA+1)%4;
	      if(!EmSon[isonB]) {//to save a little work, only lookup this
		//son if I don't already know him
		EmSon[isonB]=(Element*) El_Table->lookup(EmFather->son[isonB]);
		assert(EmSon[isonB]);
	      }	      

	      for(ikey=0;ikey<KEYLENGTH;ikey++)
		EmSon[  isonA]->neighbor[ineighm4][ikey]=
		  EmSon[isonA]->neighbor[ineighp4][ikey]=
		  EmSon[isonB]->neighbor[ineighm4][ikey]=
		  EmSon[isonB]->neighbor[ineighp4][ikey]=
		  EmFather->neighbor[ineigh][ikey];

	      EmSon[  isonA]->neigh_gen[ ineighm4]=
		EmSon[isonA]->neigh_gen[ ineighp4]=
		EmSon[isonB]->neigh_gen[ ineighm4]=
		EmSon[isonB]->neigh_gen[ ineighp4]=
		EmFather->generation;
	      
	      EmSon[  isonA]->neigh_proc[ineighm4]=
		EmSon[isonB]->neigh_proc[ineighm4]=
		neigh_proc;

	      EmSon[  isonA]->neigh_proc[ineighp4]=
		EmSon[isonB]->neigh_proc[ineighp4]=
		-2;
	      
	      
	      NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[ineigh]);
	      assert(NdTemp);
	      NdTemp->info=CORNER;

	      NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[ineighp4]);
	      assert(NdTemp);
	      NdTemp->info=S_S_CON;

	      NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[ineighp4]);
	      assert(NdTemp);
	      NdTemp->info=S_C_CON;

	      NdTemp=(Node*) NodeTable->lookup(EmSon[isonB]->node_key[ineighp4]);
	      assert(NdTemp);
	      NdTemp->info=S_S_CON;

	      NdTemp=(Node*) NodeTable->lookup(EmSon[isonB]->node_key[(ineigh+1)%4]);
	      assert(NdTemp);
	      NdTemp->info=CORNER;

	    }
	    break;
	  default:
	    printf("FUBAR in refine_neigh_update() interproc part 1\n");
	    assert(0);
	    break;
	  }//switch(EmFather->generation-EmFather->neigh_gen[ineigh])
	}//if((neigh_proc>=0)&&(neigh_proc!=myid))

      }//for(ineigh=0;ineigh<8;ineigh++)
    }//for(ifather=RefinedList->get_inewstart();...

    MPI_Alltoall(num_send,1,MPI_INT,num_recv,1,MPI_INT,MPI_COMM_WORLD);
    num_send[myid]=num_recv[myid]=0;

    int max_num_send=0,max_num_recv=0;

    for(iproc=0;iproc<nump;iproc++) {
      if(num_send[iproc]>max_num_send) max_num_send=num_send[iproc];
      if(num_recv[iproc]>max_num_recv) max_num_recv=num_recv[iproc];
      //printf("B myid=%d iproc=%d num_send=%d, num_recv=%d\n",
      //myid,iproc,num_send[iproc],num_recv[iproc]);
    }
    fflush(stdout);

    send=CAllocU2(nump,4*KEYLENGTH*max_num_send);
    recv=CAllocU2(nump,4*KEYLENGTH*max_num_recv);

    for(iproc=0;iproc<nump;iproc++)
      if((iproc!=myid)&&(num_recv[iproc]>0))       
	ierr=MPI_Irecv((void *) recv[iproc], 4*KEYLENGTH*num_recv[iproc], 
		       MPI_UNSIGNED, iproc, send_tag+iproc, MPI_COMM_WORLD, 
		       (request+nump+iproc));

    for(ifather=RefinedList->get_inewstart();
	ifather<RefinedList->get_num_elem();
	ifather++){
    
      EmFather=RefinedList->get(ifather); //Hello I'm the OLDFATHER
      assert(EmFather); //Help I've been abducted call the FBI!!!
      assert(EmFather->adapted==OLDFATHER); //sanity check
      
      for(iside=0;iside<4;iside++){
	isonA =iside;
	isonB =(iside+1)%4;

	ineigh=iside;
	neigh_proc  =EmFather->neigh_proc[ineigh  ];
	if((neigh_proc>=0)&&
	   (neigh_proc!=myid)
	   ) {
	  
	  for(ikey=0;ikey<KEYLENGTH;ikey++) {
	    //The element I want my neighbor proc to update
	    send[neigh_proc][(4*isend[neigh_proc]+0)*KEYLENGTH+ikey]=
	      EmFather->neighbor[ineigh][ikey];
	    //ME 
	    send[neigh_proc][(4*isend[neigh_proc]+1)*KEYLENGTH+ikey]=
	      EmFather->key[ikey];
	  }

	  switch(EmFather->generation-EmFather->neigh_gen[ineigh]) {
	  case -1: //I'm one generation older than my old neighbor so 
	    //I know he couldn't have just refined and I know I have to
	    //introduce him to sonA

	    for(ikey=0;ikey<KEYLENGTH;ikey++)
	      send[  neigh_proc][(4*isend[neigh_proc]+2)*KEYLENGTH+ikey]=
		send[neigh_proc][(4*isend[neigh_proc]+3)*KEYLENGTH+ikey]=
		EmFather->son[isonA][ikey];

	    break;
	  case  0: //I'm the same generation as my old neighbor so I 
	    //know I have to introduce him to both my sons
	    for(ikey=0;ikey<KEYLENGTH;ikey++) {
	      send[neigh_proc][(4*isend[neigh_proc]+2)*KEYLENGTH+ikey]=
		EmFather->son[isonB][ikey];
	      send[neigh_proc][(4*isend[neigh_proc]+3)*KEYLENGTH+ikey]=
		EmFather->son[isonA][ikey];
	    }
	    break;
	  default:
	    printf("FUBAR\n");
	    assert(0);
	    break; 
	  }

	  isend[neigh_proc]++;
	}
	
	ineigh=iside+4;
	neigh_proc  =EmFather->neigh_proc[ineigh  ];	
	if((neigh_proc>=0)&&
	   (neigh_proc!=myid)
	   ) {
	  //I know I'm one generation older than my old neighbor so
	  //I know he couldn't have just refined and I know I have to
	  //introduce them to sonB, because he's neighbor 4<=ineigh<8
	  
	  for(ikey=0;ikey<KEYLENGTH;ikey++) {
	    //The element I want my neighbor proc to update
	    send[neigh_proc][(4*isend[neigh_proc]+0)*KEYLENGTH+ikey]=
	      EmFather->neighbor[ineigh][ikey];
	    //ME 
	    send[neigh_proc][(4*isend[neigh_proc]+1)*KEYLENGTH+ikey]=
	      EmFather->key[ikey];

	    send[  neigh_proc][(4*isend[neigh_proc]+2)*KEYLENGTH+ikey]=
	      send[neigh_proc][(4*isend[neigh_proc]+3)*KEYLENGTH+ikey]=
	      EmFather->son[isonB][ikey];
	  }  

	  isend[neigh_proc]++;	  
	}
      }//for(iside=0;iside<4;iside++)

    }//for(ifather=RefinedList->get_inewstart();...

    //send the update neighbor information to the other processors
    for(iproc=0;iproc<nump;iproc++)
      if((iproc!=myid)&&(num_send[iproc]>0)) 
	ierr=MPI_Isend((void *) send[iproc], 4*KEYLENGTH*num_send[iproc], 
		       MPI_UNSIGNED, iproc, send_tag+myid, MPI_COMM_WORLD, 
		       (request+iproc));
    CDeAllocI1(isend);
  }//if(nump>1)


  /*************************************************************/
  /* now do the on processor updates while I'm waiting to      */
  /* receive neighbor update information from other processors */
  /*************************************************************/

  for(ifather=RefinedList->get_inewstart();
      ifather<RefinedList->get_num_elem();
      ifather++){

    EmFather=RefinedList->get(ifather); //Hello I'm the OLDFATHER
    assert(EmFather); //Help I've been abducted call the FBI!!!
    assert(EmFather->adapted==OLDFATHER); //sanity check
    
    NdTemp=(Node*) NodeTable->lookup(EmFather->key);
    assert(NdTemp);
    NdTemp->info=CORNER;

    //These are my sons, I'm going to introduce them to my neighbors
    for(isonA=0;isonA<4;isonA++) {
      EmSon[isonA]=(Element*) El_Table->lookup(EmFather->son[isonA]);
      assert(EmSon[isonA]); //MY son has been abducted call the FBI!!!

      NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->key);
      assert(NdTemp);
      NdTemp->info=BUBBLE;

      NdTemp=(Node*) 
	NodeTable->lookup(EmSon[isonA]->node_key[(isonA+1)%4+4]);
      assert(NdTemp);
      NdTemp->info=SIDE;
    }

    //visit my neighbors on each side
    for(iside=0;iside<4;iside++) {

      ineigh  =iside;
      ineighp4=ineigh+4;
      isonA   =ineigh;
      isonB   =(ineighp4+1)%4;

      if(EmFather->neigh_proc[ineigh]==-1) {
	//handle map boundary special
	for(ikey=0;ikey<KEYLENGTH;ikey++){
	  EmSon[  isonA]->neighbor[ineigh  ][ikey]=
	    EmSon[isonA]->neighbor[ineighp4][ikey]=
	    EmSon[isonB]->neighbor[ineigh  ][ikey]=
	    EmSon[isonB]->neighbor[ineighp4][ikey]=
	    0;
	}	
	EmSon[  isonA]->neigh_gen[ineigh  ]=
	  EmSon[isonA]->neigh_gen[ineighp4]=
	  EmSon[isonB]->neigh_gen[ineigh  ]=
	  EmSon[isonB]->neigh_gen[ineighp4]=
	  0;

	EmSon[  isonA]->neigh_proc[ineigh]=
	  EmSon[isonB]->neigh_proc[ineigh]=
	  -1;

	EmSon[  isonA]->neigh_proc[ineighp4]=
	  EmSon[isonB]->neigh_proc[ineighp4]=
	  -2;
	//if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
      }
      else if((EmFather->neigh_proc[ineigh  ]==myid)&&
	      ((EmFather->neigh_proc[ineighp4]==myid)||
	       (EmFather->neigh_proc[ineighp4]==-2))
	      ) {
	//case where one neighbor on this side is on my proc while the other 
	//is on another proc has already been handled up above, when packing 
	//the information to send to the other proc.
						
	//knock knock, Hello Neighbors
	EmNeighOld[0]=(Element*) El_Table->lookup(EmFather->neighbor[ineigh  ]);
	/*
	if(!EmNeighOld[0]){
	  printf("father element {%u,%u} is missing neighbor %d/{%u,%u}\n",
		 *(EmFather->pass_key()+0),*(EmFather->pass_key()+1),
		 ineigh,
		 EmFather->neighbor[ineigh][0],
		 EmFather->neighbor[ineigh][1]);
	  ElemBackgroundCheck(El_Table,NodeTable,EmFather->pass_key(),stdout);
	  fflush(stdout);
	  printf("***********\n");
	  ElemBackgroundCheck(El_Table,NodeTable,EmFather->neighbor[ineigh],
			      stdout);
	  fflush(stdout);
	}
	*/
	assert(EmNeighOld[0]);



	EmNeighOld[1]=(Element*) El_Table->lookup(EmFather->neighbor[ineighp4]);
	assert(EmNeighOld[1]);
	EmNeighNew[0]=EmNeighNew[1]=EmNeighNew[2]=EmNeighNew[3]=NULL;
		
	for(ineighme=0;ineighme<8;ineighme++){
	  if(compare_key(EmFather->key,EmNeighOld[1]->neighbor[ineighme]))
	    break;

	}
	if(!(ineighme<8)) {
	  printf("FUBAR 0 detected in refine_neigh_update\nEmFather={%10u,%10u}\nEmNeighOld[0]={%10u,%10u} ineigh=%d isonA=%d isonB=%d\n",EmFather->key[0],EmFather->key[1],EmNeighOld[0]->key[0],EmNeighOld[0]->key[1],ineigh,isonA,isonB);
 	  printf("aborting!!");

	  assert(ineighme<8);
	}	    

	
	//There are 5 cases I need to worry about about
	//A: my old neighbor is one generation older than me, only 
	//   possible if we've both been refined
	//B: my old neighbor is of my generation and hasn't been refined
	//C: my old neighbor is of my generation and has been refined
	//D: my old neighbor is one generation younger than me and hasn't 
	//   been refined
	//E: my old neighbor is one generation younger than me and has
	//   been refined
	//
	//I'm going to compress this into one of 3 cases I need to handle
	//0: case A: same as case E but from the other side, only need 
	//   to do once so don't do anything this time
	//1: my new neighbor is my generation
	//2: my new neighbor is one generation younger than me
	//3: my new neighbor is two generations younger than me
	
	  
	//this switch loop compresses cases
	switch(EmNeighOld[0]->generation-EmFather->generation) {
	case -1:
	  //this is a case A
	  inewcase=0;
	  assert(EmNeighOld[0]->adapted==OLDFATHER);	    
	  ineighme=ineighmep4=-1; //for sanity check
	  break;
	case 0:
	  assert(ineighme<4);
	  ineighmep4=ineighme+4;
	  
	  if(EmNeighOld[0]->adapted==OLDFATHER) {
	    //this is a case C 
	    inewcase=2;
	    
	    EmNeighNew[0]=(Element*) 
	      El_Table->lookup(EmNeighOld[0]->son[(ineighme+1)%4]);
	    assert(EmNeighNew[0]);
	    
	    EmNeighNew[1]=(Element*) 
	      El_Table->lookup(EmNeighOld[0]->son[ineighme]);
	    assert(EmNeighNew[1]);
	    EmNeighNew[2]=EmNeighNew[3]=NULL;
	  }
	  else{
	    //this is a case B
	    if(EmNeighOld[0]->adapted>TOBEDELETED){
	      inewcase=1;
	      EmNeighNew[0]=EmNeighOld[0];
	      EmNeighNew[1]=EmNeighNew[2]=EmNeighNew[3]=NULL;
	    }
	    else inewcase=0;

	  }
	  break;
	case 1:
	  assert(ineighme<4);
	  ineighmep4=ineighme+4;

	  if((EmNeighOld[0]->adapted==OLDFATHER)||
	     (EmNeighOld[1]->adapted==OLDFATHER)
	     ){
	    //this is a case E
	    inewcase=3;

	    if(EmNeighOld[0]->adapted==OLDFATHER) {
	      EmNeighNew[0]=(Element*) 
		El_Table->lookup(EmNeighOld[0]->son[(ineighme+1)%4]);
	      assert(EmNeighNew[0]);

	      EmNeighNew[1]=(Element*) 
		El_Table->lookup(EmNeighOld[0]->son[ineighme]);
	      assert(EmNeighNew[1]);
	    }
	    else
	      EmNeighNew[1]=EmNeighNew[0]=EmNeighOld[0];
	    

	    if(EmNeighOld[1]->adapted==OLDFATHER) {
	      EmNeighNew[2]=(Element*) 
		El_Table->lookup(EmNeighOld[1]->son[(ineighme+1)%4]);
	      assert(EmNeighNew[2]);
	      
	      EmNeighNew[3]=(Element*) 
		El_Table->lookup(EmNeighOld[1]->son[ineighme]);
	      assert(EmNeighNew[3]);
	    }
	    else
	      EmNeighNew[3]=EmNeighNew[2]=EmNeighOld[1];
	  }
	  else{
	    //this is a case D
	    inewcase=2;
	    
	    EmNeighNew[0]=EmNeighOld[0];
	    EmNeighNew[1]=EmNeighOld[1];
	    EmNeighNew[2]=EmNeighNew[3]=NULL;
	  }
	  break;
	default:
	  inewcase=-1;
	  
	  printf("FUBAR 1 detected in refine_neigh_update! aborting.\n");
	  assert(0);
	  break;
	} //switch based on difference in generation between me and my old neighbor, this is used to reduce the number of cases from 5 to 3 (based on new neighbor generation)
	
	//sanity check
	assert((ineigh>=0)&&(ineigh<4));
	assert(ineighp4==ineigh+4);
	if(inewcase){
	  assert((ineighme>=0)&&(ineighme<4)); 
	  assert(ineighmep4==ineighme+4);
	}
		  
	//now only deal with the new cases, and yes I know that I
	//am resetting neighbor information in ghost cells but 
	//not neighbor information of the original cells on other 
	//processors, I'm going to fix that in a minute
	switch(inewcase) {
	case 0: 
	  //case A
	  break;
	case 1: 
	  //case B
	  //new neighbor generation is my (the OLDFATHER) generation
	  for(ikey=0;ikey<KEYLENGTH;ikey++) {
	    EmNeighNew[ 0]->neighbor[ineighme  ][ikey]=
	      EmSon[isonB]->key[ikey];
	    EmNeighNew[ 0]->neighbor[ineighmep4][ikey]=
	      EmSon[isonA]->key[ikey];
	    
	    EmSon[   isonA]->neighbor[ineigh  ][ikey]=
	      EmSon[ isonA]->neighbor[ineighp4][ikey]=
	      EmSon[ isonB]->neighbor[ineigh  ][ikey]=
	      EmSon[ isonB]->neighbor[ineighp4][ikey]=
	      EmNeighNew[0]->key[ikey];
	  }

	  EmNeighNew[  0]->neigh_gen[ineighme  ]=
	    EmNeighNew[0]->neigh_gen[ineighmep4]=
	    EmSon[ isonA]->generation;
	  
	  EmSon[   isonA]->neigh_gen[ineigh  ]=
	    EmSon[ isonA]->neigh_gen[ineighp4]=
	    EmSon[ isonB]->neigh_gen[ineigh  ]=
	    EmSon[ isonB]->neigh_gen[ineighp4]=
	    EmNeighNew[0]->generation;
	  
	  EmSon[   isonA]->neigh_proc[ineighp4]=
	    EmSon[ isonB]->neigh_proc[ineighp4]=
	    -2;

	  EmNeighNew[  0]->neigh_proc[ineighme  ]=
	    EmNeighNew[0]->neigh_proc[ineighmep4]=
	    EmFather->myprocess;

	  //update the nodes on this side
	  //The new difference in generation tells me the OLDFATHER's 
	  //2 corner nodes on this side are actually CORNER's and not
	  //S_C_CON's
	  
	  inode=ineigh;
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);

	  inode=ineighp4;
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(S_S_CON);
	  
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(S_C_CON);
	  
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(S_S_CON);

	  inode=(ineigh+1)%4;
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);
	  
	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);

	  break;
	case 2:
	  //cases C & D
	  //new neighbor generation is my son's generation

	  for(ikey=0;ikey<KEYLENGTH;ikey++) {
	    
	    EmNeighNew[  0]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[0]->neighbor[ineighmep4][ikey]=
	      EmSon[ isonA]->key[ikey];
	    EmNeighNew[  1]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[1]->neighbor[ineighmep4][ikey]=
	      EmSon[ isonB]->key[ikey];
	    
	    EmSon[   isonA]->neighbor[ineigh  ][ikey]=
	      EmSon[ isonA]->neighbor[ineighp4][ikey]=
	      EmNeighNew[0]->key[ikey];
	    EmSon[   isonB]->neighbor[ineigh  ][ikey]=
	      EmSon[ isonB]->neighbor[ineighp4][ikey]=
	      EmNeighNew[1]->key[ikey];
	  }
	  
	  EmNeighNew[  0]->neigh_gen[ineighme  ]=
	    EmNeighNew[0]->neigh_gen[ineighmep4]=
	    EmNeighNew[1]->neigh_gen[ineighme  ]=
	    EmNeighNew[1]->neigh_gen[ineighmep4]=
	    EmSon[ isonA]->generation;
	  
	  EmNeighNew[  0]->neigh_proc[ineighmep4]=
	    EmNeighNew[1]->neigh_proc[ineighmep4]=
	    EmSon[ isonA]->neigh_proc[ineighp4  ]=
	    EmSon[ isonB]->neigh_proc[ineighp4  ]=
	    -2;
	  
	  EmSon[   isonA]->neigh_gen[ineigh  ]=
	    EmSon[ isonA]->neigh_gen[ineighp4]=
	    EmSon[ isonB]->neigh_gen[ineigh  ]=
	    EmSon[ isonB]->neigh_gen[ineighp4]=
	    EmNeighNew[0]->generation;
	  
	  EmSon[   isonA]->neigh_proc[ineigh]=
	    EmNeighNew[0]->myprocess;
	  
	  EmSon[   isonB]->neigh_proc[ineigh]=
	    EmNeighNew[1]->myprocess;
	  
	  EmNeighNew[  0]->neigh_proc[ineighme]=
	    EmNeighNew[1]->neigh_proc[ineighme]=
	    EmFather->myprocess;

	  //update the nodes on this side
	  //don't update my corner nodes because they could be S_C_CON's
	  //if they should be S_C_CON's and I reset them to CORNERs I
	  //will no longer conserve mass/volume in a dramatically 
	  //observable fashion

	  if(EmSon[isonA]->neigh_gen[(ineigh+3)%4]==
	     EmSon[isonA]->generation) {
	    //neighbor before (tested here) and after this (the ineigh)
	    //corner (i.e. the ineigh neighbor) are the same generation 
	    //as me, therefor this (the ineigh) node is a CORNER and not
	    //an S_C_CON node
	    inode=ineigh;
	    NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(CORNER);
	  }

	  inode=ineigh+4;
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(SIDE);
	  
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);
	  
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(SIDE);
	  
	  if(EmSon[isonB]->neigh_gen[(ineigh+1)%4]==
	     EmSon[isonB]->generation) {
	    //neighbor before (i.e. the ineigh neighbor) and after 
	    //(tested here) this (the (ineigh+1)%4) corner are the
	    //the same generation as me, therefore this (the 
	    //(ineigh+1)%4) node is a CORNER and not an S_C_CON node
	    inode=(ineigh+1)%4;
	    NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(CORNER);
	  }
	  
	  break;
	case 3:
	  //case E

	  //update the nodes on this side	  

	  inode=ineigh; //father corner node
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);

	  inode=ineighp4; //father edge node
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);

	  inode=(ineigh+1)%4; //father corner node
	  NdTemp=(Node*) NodeTable->lookup(EmFather->node_key[inode]);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);
	  
	  for(ikey=0;ikey<KEYLENGTH;ikey++) {
	    
	    EmNeighNew[  0]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[0]->neighbor[ineighmep4][ikey]=
	      EmNeighNew[1]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[1]->neighbor[ineighmep4][ikey]=
	      EmSon[ isonA]->key[ikey];
	    
	    EmNeighNew[  2]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[2]->neighbor[ineighmep4][ikey]=
	      EmNeighNew[3]->neighbor[ineighme  ][ikey]=
	      EmNeighNew[3]->neighbor[ineighmep4][ikey]=
	      EmSon[ isonB]->key[ikey];
	    
	    EmSon[   isonA]->neighbor[ineigh  ][ikey]=
	      EmNeighNew[0]->key[ikey];
	    EmSon[   isonA]->neighbor[ineighp4][ikey]=
	      EmNeighNew[1]->key[ikey];
	    
	    EmSon[   isonB]->neighbor[ineigh  ][ikey]=
	      EmNeighNew[2]->key[ikey];
	    EmSon[   isonB]->neighbor[ineighp4][ikey]=
	      EmNeighNew[3]->key[ikey];
	  }

	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);

	    
	  EmNeighNew[  0]->neigh_gen[ineighme  ]=
	    EmNeighNew[0]->neigh_gen[ineighmep4]=
	    EmNeighNew[1]->neigh_gen[ineighme  ]=
	    EmNeighNew[1]->neigh_gen[ineighmep4]=
	    EmNeighNew[2]->neigh_gen[ineighme  ]=
	    EmNeighNew[2]->neigh_gen[ineighmep4]=
	    EmNeighNew[3]->neigh_gen[ineighme  ]=
	    EmNeighNew[3]->neigh_gen[ineighmep4]=
	    EmSon[ isonA]->generation;
	  
	  EmSon[   isonA]->neigh_gen[ineigh  ]=
	    EmSon[ isonA]->neigh_gen[ineighp4]=
	    EmNeighNew[0]->generation;
	  
	  EmSon[   isonB]->neigh_gen[ineigh  ]=
	    EmSon[ isonB]->neigh_gen[ineighp4]=
	    EmNeighNew[2]->generation;
	  
	  EmNeighNew[  0]->neigh_proc[ineighmep4]=
	    EmNeighNew[1]->neigh_proc[ineighmep4]=
	    EmNeighNew[2]->neigh_proc[ineighmep4]=
	    EmNeighNew[3]->neigh_proc[ineighmep4]=
	    -2;

	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);

	  EmSon[  isonA]->neigh_proc[ineigh  ]=
	    EmFather->neigh_proc[ineigh];
	  
	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);


	  EmSon[  isonB]->neigh_proc[ineigh  ]=
	    EmFather->neigh_proc[ineighp4];
	  
	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);

	  inode=ineighp4; //sonA edge node 
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonA]->node_key[inode]);
	  assert(NdTemp);
	  if(compare_key(EmSon[isonA]->neighbor[ineigh  ],
			 EmSon[isonA]->neighbor[ineighp4])) {
	    EmSon[isonA]->neigh_proc[ineighp4]=-2;
	    NdTemp->putinfo(SIDE);
	  }
	  else{
	    EmSon[  isonA]->neigh_proc[ineighp4]=
	      EmSon[isonA]->neigh_proc[ineigh  ];
	    
	    NdTemp->putinfo(S_C_CON);
	    
	    inode=ineighmep4;

	    NdTemp=(Node*) NodeTable->lookup(EmNeighNew[0]->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(S_S_CON);

	    NdTemp=(Node*) NodeTable->lookup(EmNeighNew[1]->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(S_S_CON);
	  }
	  
	  //if(Curr_El) if(IfNeighProcChange(El_Table,NodeTable,myid,Curr_El,EmFather)) assert(0);
	  inode=ineighp4; //sonB edge node
	  NdTemp=(Node*) NodeTable->lookup(EmSon[isonB]->node_key[inode]);
	  assert(NdTemp);
	  if(compare_key(EmSon[isonB]->neighbor[ineigh  ],
			 EmSon[isonB]->neighbor[ineighp4])) {
	    EmSon[isonB]->neigh_proc[ineighp4]=-2;
	    NdTemp->putinfo(SIDE);
	  }
	  else{
	    EmSon[  isonB]->neigh_proc[ineighp4]=
	      EmSon[isonB]->neigh_proc[ineigh  ];
	    
	    NdTemp->putinfo(S_C_CON);

	    inode=ineighmep4;
	    
	    NdTemp=(Node*) NodeTable->lookup(EmNeighNew[2]->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(S_S_CON);
	    
	    NdTemp=(Node*) NodeTable->lookup(EmNeighNew[3]->node_key[inode]);
	    assert(NdTemp);
	    NdTemp->putinfo(S_S_CON);
	  }
	    
	  break;
	default:
	  printf("FUBAR 2 detected in refine_neigh_update! aborting.\n");
	  assert(0);

	  break;
	} //switch(inewcase), case based on generation of my new neighbor

      } //else: not a map boundary

    } //iside loop
      
  } //ifather loop
  
  /*************************************************************/
  /* The interprocessor update information should be here by   */
  /* now or at the very least I won't have to weight very long */
  /* receive neighbor update information from other processors */
  /*************************************************************/

  //now complete the interprocessor update
  if(nump>1){
    int NumProcsNotRecvd, ifrecvd, nodeorder, irecvd;
    MPI_Status status;
    Element *EmTemp, *EmNeigh;

    do{

      for(iproc=0;iproc<nump;iproc++){
	//printf("stop me\n");
	if(num_recv[iproc]>0) { 
	  //only check processors I haven't already handled
	  ifrecvd=0;
	  //printf("myid=%d before Test",myid); fflush(stdout);
	  MPI_Test(request+nump+iproc,&ifrecvd,&status);
	  //printf("myid=%d after Test",myid); fflush(stdout);
	  if(ifrecvd) {

	    for(irecvd=0;irecvd<num_recv[iproc];irecvd++) {
	      //one by one check my Element's that my neighbor processor
	      //says I need to update

	      //Hi I'm EmTemp
	      EmTemp=(Element*) 
		El_Table->lookup(&(recv[iproc][(4*irecvd+0)*KEYLENGTH]));
	      /*
	      if(compare_key(EmTemp->key,ElemDebugKey)){
		ElemBackgroundCheck(El_Table,NodeTable,ElemDebugKey,fpbg);
		printf("stop me\n");		
	      }
	      */
	      if(!EmTemp){
		printf("refine_neigh_update(): myid=%d receiving from iproc=%d, %dth element being received is {%10u,%10u} but this off processor neighbor is missing\n",myid,iproc,irecvd,
		       recv[iproc][(4*irecvd+0)*KEYLENGTH+0],
		       recv[iproc][(4*irecvd+0)*KEYLENGTH+1]);
		ElemBackgroundCheck(El_Table,NodeTable,&(recv[iproc][(4*irecvd+0)*KEYLENGTH]),stdout);
		assert(EmTemp);
	      }

	      //my old neighbor will introduce his NEWSONs to me
	      for(ineigh=0;ineigh<4;ineigh++) 
		if(compare_key(EmTemp->neighbor[ineigh], 
			       &(recv[iproc][(4*irecvd+1)*KEYLENGTH])))
		  break;
	      assert(ineigh<4); //I don't know this Element pretending
	      //to be my neighbor, I'm calling the FBI to report a 
	      //suspicious person... err.. suspicious Element
	      //note it is impossible for my old neighbor to be younger
	      //than me because there is no direct triggered refinement
	      //across the interprocessor boundary, the "adapted" flag
	      //value -BUFFER is an _instruction_ by my neighbor on 
	      //another processor telling me to refine, but he can't 
	      //refine unless he is the same age as me or older than me.
	      ineighp4=ineigh+4;

	      EmNeigh=(Element*) 
		El_Table->lookup(&(recv[iproc][(4*irecvd+1)*KEYLENGTH]));
	      if(EmNeigh){
		if(!compare_key(&(recv[iproc][(4*irecvd+1)*KEYLENGTH]),EmNeigh->pass_key())) {
		  printf("refine_neigh_update(): myid=%d receiving from iproc=%d, %dth element being received is {%10u,%10u}\n That element is supposed to have an onprocessor element {%10u,%10u} as a neighbor but that \"on processor element\" thinks it key is {%10u,%10u}.\n",myid,iproc,irecvd,
			 recv[iproc][(4*irecvd+0)*KEYLENGTH+0],
			 recv[iproc][(4*irecvd+0)*KEYLENGTH+1],
			 recv[iproc][(4*irecvd+1)*KEYLENGTH+0],
			 recv[iproc][(4*irecvd+1)*KEYLENGTH+1],
			 *(EmNeigh->pass_key()+0),*(EmNeigh->pass_key()+1));
		  ElemBackgroundCheck(El_Table,NodeTable,&(recv[iproc][(4*irecvd+0)*KEYLENGTH]),stdout);
		  ElemBackgroundCheck(El_Table,NodeTable,&(recv[iproc][(4*irecvd+1)*KEYLENGTH]),stdout);
		  ElemBackgroundCheck(El_Table,NodeTable,EmNeigh->pass_key(),stdout);
		  assert(EmNeigh);
		}

		El_Table->remove(EmNeigh->key, 1,stdout,myid,13);
		EmNeigh->void_bcptr();
		delete EmNeigh;	    
	      }


	      //EmFather=(Element*) El_Table->lookup(EmTemp->father);

	      if(EmTemp->adapted==OLDFATHER) {
		//I know my neighbor was my generation, and we both refined
		//so the corners I share with this old neighbor are CORNERs
		//and not S_C_CONs
		isonA=ineigh;
		isonB=(ineigh+1)%4;

		inode=isonA;
		NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		assert(NdTemp);
		NdTemp->info=CORNER;

		inode=isonB;
		NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		assert(NdTemp);
		NdTemp->info=CORNER;

		EmSon[0]=EmSon[1]=EmSon[2]=EmSon[3]=NULL;
		EmSon[isonA]=(Element*) 
		  El_Table->lookup(EmTemp->son[isonA]); 
		assert(EmSon[isonA]);

		EmSon[isonB]=(Element*) 
		  El_Table->lookup(EmTemp->son[isonB]); 
		assert(EmSon[isonA]);

		for(ikey=0;ikey<KEYLENGTH;ikey++) {
		  EmSon[  isonA]->neighbor[ineigh  ][ikey]=
		    EmSon[isonA]->neighbor[ineighp4][ikey]=
		    recv[iproc][(4*irecvd+2)*KEYLENGTH+ikey];


		  //bob
		  EmSon[  isonB]->neighbor[ineigh  ][ikey]=
		    EmSon[isonB]->neighbor[ineighp4][ikey]=
		    recv[iproc][(4*irecvd+3)*KEYLENGTH+ikey];
		}
		
		EmSon[  isonA]->neigh_gen[ineigh  ]=
		  EmSon[isonA]->neigh_gen[ineighp4]=
		  EmSon[isonB]->neigh_gen[ineigh  ]=
		  EmSon[isonB]->neigh_gen[ineighp4]=
		  EmSon[isonA]->generation;

		EmSon[  isonA]->neigh_proc[ineigh  ]=
		  EmSon[isonB]->neigh_proc[ineigh  ]=
		  iproc;

		EmSon[  isonA]->neigh_proc[ineighp4]=
		  EmSon[isonB]->neigh_proc[ineighp4]=
		  -2;

		NdTemp=(Node*) 
		  NodeTable->lookup(EmTemp->node_key[ineighp4]);
		assert(NdTemp);
		NdTemp->info=CORNER;

		NdTemp=(Node*) 
		  NodeTable->lookup(EmSon[isonA]->node_key[ineighp4]);
		assert(NdTemp);
		NdTemp->info=SIDE;

		NdTemp=(Node*) 
		  NodeTable->lookup(EmSon[isonB]->node_key[ineighp4]);
		assert(NdTemp);
		NdTemp->info=SIDE;
	      }
	      else{ //my oldneighbor was either a generation older than me
		//or he was my generation, 

		if(EmTemp->neigh_gen[ineigh]==EmTemp->generation){
		  //he was my generation so the two corner nodes we 
		  //shared are CORNER and not S_C_CONs
		  inode=ineigh;
		  NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		  assert(NdTemp);
		  NdTemp->info=CORNER;

		  inode=(ineigh+1)%4;
		  NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		  assert(NdTemp);
		  NdTemp->info=CORNER;
		}
		else{
		  //he was a generation older than me (his son is now my 
		  //generation so if my neighbor element's who share the 
		  //corner's I share with the new son are my generation, 
		  //then these shared corners are CORNERs and not S_C_CONs 
		  if(EmTemp->neigh_gen[(ineigh+3)%4]==EmTemp->generation) {
		    inode=ineigh;
		    NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		    assert(NdTemp);
		    NdTemp->info=CORNER;
		  }

		  if(EmTemp->neigh_gen[(ineigh+1)%4]==EmTemp->generation) {
		    inode=(ineigh+1)%4;
		    NdTemp=(Node*) NodeTable->lookup(EmTemp->node_key[inode]);
		    assert(NdTemp);
		    NdTemp->info=CORNER;
		  }
		}


		for(ikey=0;ikey<KEYLENGTH;ikey++) {
		  EmTemp->neighbor[ineigh  ][ikey]=
		    recv[iproc][(4*irecvd+2)*KEYLENGTH+ikey];

		  EmTemp->neighbor[ineighp4][ikey]=
		    recv[iproc][(4*irecvd+3)*KEYLENGTH+ikey];
		}
		EmTemp->neigh_gen[ineighp4]=
		  EmTemp->neigh_gen[ineigh]=
		  EmTemp->neigh_gen[ineigh]+1;

		EmTemp->neigh_proc[ineigh]=iproc;
		
		NdTemp=(Node*) 
		  NodeTable->lookup(EmTemp->node_key[ineighp4]);
		assert(NdTemp);

		if(compare_key(EmTemp->neighbor[ineigh  ],
			       EmTemp->neighbor[ineighp4])) {
		  EmTemp->neigh_proc[ineighp4]=-2;
		  NdTemp->info=SIDE;
		}
		else{
		  EmTemp->neigh_proc[ineighp4]=iproc;
		  NdTemp->info=S_C_CON;
		  /*
		  if(compare_key(EmTemp->key,ElemDebugKey)){
		    ElemBackgroundCheck(El_Table,NodeTable,ElemDebugKey,fpbg);
		    printf("stop me 2\n");		
		  }
		  */
		}

		if(EmTemp->generation>EmTemp->neigh_gen[ineigh]) {
		  NdTemp=(Node*) 
		    NodeTable->lookup(EmTemp->node_key[EmTemp->get_which_son()]);
		  assert(NdTemp);
		  NdTemp->info=CORNER;
		}
		//ElemBackgroundCheck(El_Table,NodeTable,EmTemp->key,fpbg);
		
	      }//else

	    }//for(irecvd=0;irecvd<num_recv[iproc];irecvd++)
	    num_recv[iproc]=0;	      
	  }//if(ifrecvd)
	}//if(num_recv[iproc]>0)
	
      }//for(iproc=0;iproc<nump;iproc++)
      NumProcsNotRecvd=0;
      for(iproc=0;iproc<nump;iproc++)
	if(num_recv[iproc]>0) NumProcsNotRecvd++;

    }while(NumProcsNotRecvd>0);
    
    CDeAllocU2(recv);
    CDeAllocI1(num_recv);
  }//if(nump>1)

  Element *EmTemp;
  HashEntryPtr* bucketzero = El_Table->getbucketptr();
  int ibucketdebugneigh=El_Table->hash(ElemDebugNeighKey);
  HashEntryPtr entryp;

  /*
  if(timeprops_ptr->iter==2389) {


    ElemDebugFatherNeigh=ElemDebugFather=NULL;
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    if(ElemDebug) {
      printf("ElemDebug={%10u,%10u}============================\n",
	     ElemDebugKey[0],ElemDebugKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugKey,stdout);
      ElemDebugFatherKey[0]=ElemDebug->father[0];
      ElemDebugFatherKey[1]=ElemDebug->father[1];
      ElemDebugFather=(Element*) El_Table->lookup(ElemDebugFatherKey);

    }
    if(ElemDebugNeigh) {
      printf("ElemDebugNeigh={%10u,%10u}=======================\n",
	     ElemDebugNeighKey[0],ElemDebugNeighKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugNeighKey,stdout);      
    }
    if(ElemDebugFather) {
      printf("ElemDebugFather={%10u,%10u}======================\n",
	     ElemDebugFatherKey[0],ElemDebugFatherKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherKey,stdout);
      ineigh=-1;
      if(ElemDebug->which_son==3) ineigh=3;
      else if(ElemDebug->which_son==0) ineigh=7;
      if(ineigh!=-1) {
	ElemDebugFatherNeighKey[0]=ElemDebug->neighbor[ineigh][0];
	ElemDebugFatherNeighKey[1]=ElemDebug->neighbor[ineigh][1];
	ElemDebugFatherNeigh=(Element*) 
	  El_Table->lookup(ElemDebugFatherNeighKey);
      }
      printf("Father ineigh=%d\n",ineigh);
    }
    if(ElemDebugFatherNeigh) {
      printf("ElemDebugFatherNeigh={%10u,%10u} ineigh=%d========\n",
	     ElemDebugFatherNeighKey[0],ElemDebugFatherNeighKey[1],ineigh);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherNeighKey,stdout);
    }
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");


    int ibucket, ilink;
    int icounter=0;
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
    for(ifather=0;ifather<RefinedList->get_num_elem();ifather++) {
      if(compare_key(RefinedList->get_key(ifather),ElemDebugNeighKey)) {
	printf("Lookup by key\n");
	ElemBackgroundCheck(El_Table,NodeTable,ElemDebugNeighKey,stdout);
	printf("Stored in RefinedList\n");
	ElemBackgroundCheck2(El_Table,NodeTable,RefinedList->get(ifather),stdout);
      }
      else if(RefinedList->get(ifather)==ElemDebugNeigh) {
	printf("RefinedList ifather=%d has same address as ElemDebugNeigh\n");

      }
      else{
	icounter++;

	ibucket=El_Table->hash(RefinedList->get_key(ifather));
	if(ibucket==ibucketdebugneigh)
	  printf("missing elem={%10u,%10u} is in same bucket as oldfather={%10u,%10u}\n",
		 ElemDebugNeighKey[0],ElemDebugNeighKey[1],
		 *(RefinedList->get_key(ifather)+0),
		 *(RefinedList->get_key(ifather)+1));
	fflush(stdout);

	entryp=*(bucketzero+ibucket);
	ilink=0;
	while(entryp) {
	  EmTemp=(Element*) (entryp->value);	  
	  if(compare_key(EmTemp->key,RefinedList->get_key(ifather)))
	    break;
	  ilink++;
	  entryp=entryp->next;	  
	}

	if(entryp) {
	  if(ilink>0) {
	    EmTemp=(Element*) (entryp->value);	  
	    if(compare_key(EmTemp->key,ElemDebugNeighKey)) 
	      printf("missing elem={%10u,%10u} is hashentryp previous to oldfather={%10u,%10u}\n",
		     ElemDebugNeighKey[0],ElemDebugNeighKey[1],
		     *(RefinedList->get_key(ifather)+0),
		     *(RefinedList->get_key(ifather)+1));
	  }

	  if(entryp->next) {
	    EmTemp=(Element*) ((entryp->next)->value);
	    if(compare_key(EmTemp->key,ElemDebugNeighKey)) 
	      printf("missing elem={%10u,%10u} is hashentryp after oldfather={%10u,%10u}\n",
		     ElemDebugNeighKey[0],ElemDebugNeighKey[1],
		     *(RefinedList->get_key(ifather)+0),
		     *(RefinedList->get_key(ifather)+1));
	    
	  }
	}
      }
    }
    printf("icounter=%d num refined elements=%d\n",icounter,
	   RefinedList->get_num_elem());	   
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
  }
  fflush(stdout);
  
  entryp=*(bucketzero+ibucketdebugneigh);
  int ilink=0;
  while(entryp) {
    EmTemp=(Element*) (entryp->value);	  
    if(compare_key(EmTemp->key,ElemDebugNeighKey))
      break;
    ilink++;
    entryp=entryp->next;	  
  }
  if(!entryp) ilink=-1;
  int ilinkdebugneigh=ilink;
  HashEntryPtr ElemDebugNeighHEP=entryp;


  int ibucketkidnapper=El_Table->hash(ElemDebugKidnapperKey);
  entryp=*(bucketzero+ibucketkidnapper);
  ilink=0;
  while(entryp) {
    EmTemp=(Element*) (entryp->value);	  
    if(compare_key(EmTemp->key,ElemDebugKidnapperKey))
      break;
    ilink++;
    entryp=entryp->next;	  
  }
  HashEntryPtr ElemKidnapperHEP=entryp;

  int ilinkyada=ilink;
  if(!entryp) ilink=-1;
  int ilinkkidnapper=ilink;
  ilink=ilinkyada;
  while(entryp) {
    EmTemp=(Element*) (entryp->value);	  
    if(compare_key(EmTemp->key,ElemDebugNeighKey))
      break;
    ilink++;
    entryp=entryp->next;
  }
  int ilinkyada2=ilink;
  if(!entryp) ilink=-1;
  ilinkyada=ilink;

  //must do this in a seperate loop
  for(ifather=0; //RefinedList->get_inewstart();
      ifather<RefinedList->get_num_elem();
      ifather++){

    EmFather=RefinedList->get(ifather); //Hello I'm the OLDFATHER
    assert(EmFather); //Help I've been abducted call the FBI!!!
    if(!(EmFather->adapted==OLDFATHER)){
      printf("ifather=%d**************************************\n",ifather);
      ElemBackgroundCheck(El_Table,NodeTable,EmFather->pass_key(),stdout);
    }
    assert(EmFather->adapted==OLDFATHER); //sanity check
    assert(EmFather->refined==1);
    EmFather->adapted=TOBEDELETED; //I've lived a good life, it's my time to die
    El_Table->remove(EmFather->key,1,stdout,myid,14);
    if(ElemDebugNeigh){
      EmTemp=(Element*) El_Table->lookup(ElemDebugNeighKey);
      if(!EmTemp) {
	printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	printf("Element={%10u,%10u} removed from hashtable at same time as refined father={%10u,%10u} ifather=%d\n",
	       ElemDebugNeighKey[0],ElemDebugNeighKey[1],
	       EmFather->key[0],EmFather->key[1],ifather);
	fflush(stdout);
	printf("MissingElem ibucket=%d ilink=%d address=%x ilinkyada=%d\n",
	       ibucketdebugneigh,ilinkdebugneigh,ElemDebugNeigh,ilinkyada);
	ElemBackgroundCheck2(El_Table,NodeTable,ElemDebugNeigh,stdout);
	printf("KidnapperElem ibucket=%d ilink=%d address=%x\n",
	       ibucketkidnapper,ilinkkidnapper,ElemDebugKidnapper);
	ElemBackgroundCheck2(El_Table,NodeTable,EmFather,stdout);
	printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	ElemDebugNeigh=NULL;
      }
    }
    EmFather->void_bcptr();
    delete EmFather;

  }
  */

  for(ifather=0; //RefinedList->get_inewstart();
      ifather<RefinedList->get_num_elem();
      ifather++) {

    EmFather=RefinedList->get(ifather); //Hello I'm the OLDFATHER
    assert(EmFather); //Help I've been abducted call the FBI!!!
    assert(EmFather->adapted==OLDFATHER); //sanity check
    assert(EmFather->refined==1);
    EmFather->adapted=TOBEDELETED; //I've lived a good life, it's my time to die
    EmFather->void_bcptr();
    El_Table->remove(EmFather->key,1,stdout,myid,15);
    delete EmFather;
  }
  //RefinedList->set_inewstart(RefinedList->get_num_elem());

  /*
  if(myid==TARGETPROC)
    AssertMeshErrorFree(El_Table,NodeTable,nump,myid);
  */
  
  //clear the refined list
  RefinedList->trashlist();

  if(nump>1) {
    MPI_Barrier(MPI_COMM_WORLD);
    CDeAllocU2(send);
    CDeAllocI1(num_send);
  }

  /*
  if(timeprops_ptr->iter==2389) {
    ElemDebugFatherNeigh=ElemDebugFather=NULL;
    printf("//////////////////////////////////////////////////////////\n");
    if(ElemDebug) {
      printf("ElemDebug={%10u,%10u}============================\n",
	     ElemDebugKey[0],ElemDebugKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugKey,stdout);
      ElemDebugFatherKey[0]=ElemDebug->father[0];
      ElemDebugFatherKey[1]=ElemDebug->father[1];
      ElemDebugFather=(Element*) El_Table->lookup(ElemDebugFatherKey);

    }
    if(ElemDebugNeigh) {
      printf("ElemDebugNeigh={%10u,%10u}=======================\n",
	     ElemDebugNeighKey[0],ElemDebugNeighKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugNeighKey,stdout);      
    }
    if(ElemDebugFather) {
      printf("ElemDebugFather={%10u,%10u}======================\n",
	     ElemDebugFatherKey[0],ElemDebugFatherKey[1]);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherKey,stdout);
      ineigh=-1;
      if(ElemDebug->which_son==3) ineigh=3;
      else if(ElemDebug->which_son==0) ineigh=7;
      if(ineigh!=-1) {
	ElemDebugFatherNeighKey[0]=ElemDebug->neighbor[ineigh][0];
	ElemDebugFatherNeighKey[1]=ElemDebug->neighbor[ineigh][1];
	ElemDebugFatherNeigh=(Element*) 
	  El_Table->lookup(ElemDebugFatherNeighKey);
      }
      printf("Father ineigh=%d\n",ineigh);
    }
    if(ElemDebugFatherNeigh) {
      printf("ElemDebugFatherNeigh={%10u,%10u} ineigh=%d========\n",
	     ElemDebugFatherNeighKey[0],ElemDebugFatherNeighKey[1],ineigh);
      ElemBackgroundCheck(El_Table,NodeTable,ElemDebugFatherNeighKey,stdout);
    }
    printf("//////////////////////////////////////////////////////////\n");
  }
  */

  return;
}


#endif


extern void update_neighbor_interprocessor(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
					   refined_neighbor* refined_start, int myid, int nump);

void update_neighbor_info(HashTable* HT_Elem_Ptr, ElemPtrList* RefinedList, 
			  int myid, int nump, 
			  HashTable* HT_Node_Ptr, int h_count )/*-h_count for  debugging-*/
{
  int son_order[2];
  Element* EmTemp;
  Element* Neighbor;
  Element* SonTemp;
  Node*    NdTemp;
  unsigned* orig_neighbors;
  int* orig_neigh_proc;
  int which_side;
  unsigned NewNeighbor[2][KEYLENGTH];
  
  int SIDE_SONS[4][2]={{0, 1}, {1, 2}, {2, 3}, {3, 0}};
 
  refined_neighbor* refined_start=new refined_neighbor();
  refined_neighbor* refined_current=refined_start;
  refined_neighbor* refined_new;

  //printf("H_ADAPT  IN  SUBDOMAIN %d  ELEMENTS REFINED: %d\n",myid, count);  


  unsigned* Sons;
  unsigned* NeighSons;
  int MyGeneration;
  int NeighGeneration;
  int NeighRefined;
  unsigned* Mykey;
  int reg;
  int a, b, k;
    

  for(int i=0; i<RefinedList->get_num_elem(); i++) {
      
    EmTemp=RefinedList->get(i);//-- element ready for refinement


    orig_neighbors=EmTemp->get_neighbors();
    orig_neigh_proc=EmTemp->get_neigh_proc();
    Sons=EmTemp->getson();
    MyGeneration=EmTemp->get_gen();
    Mykey=EmTemp->pass_key();

    assert((EmTemp->get_adapted_flag()==OLDFATHER)||
	   (EmTemp->get_adapted_flag()==TOBEDELETED));

    if(EmTemp->get_adapted_flag()==OLDFATHER) {
      EmTemp->put_adapted_flag(TOBEDELETED); //mark old father for deletion

      //update the neighbors of the new sons
      for(int j = 0; j<8; j++)
	{
	  if((*(orig_neigh_proc+j)>=0)&&(*(orig_neigh_proc+j)==myid)) //neighbor is at the same proc.
	    {
	      Neighbor = (Element*) HT_Elem_Ptr->lookup(orig_neighbors+j*KEYLENGTH);	
	      assert(Neighbor);
	      NeighGeneration=Neighbor->get_gen();
	      NeighRefined=Neighbor->get_refined_flag();

	      /*-------------------------------------------------
		  case 1: neighbor has the same gen. and was not refined in this step
	          case 2: neighbor has the same gen. and was refined in this step
		  case 3: neighbor is older and was refined
		  case 4: neighbor is younger and was refined
		  case 5: neighbor is younger and was not refined
		  other : wrong request
	       *-----------------------------------------------*/
	      if(!NeighRefined)
		{
		  assert(NeighGeneration>=MyGeneration);
		  which_side=Neighbor->which_neighbor(Mykey);
		  assert(which_side<4);

		  if(NeighGeneration == MyGeneration)//--case 1
		    {	
		      a=j+1;
		      if (a==4) a=0; 
		      b = j;
		      reg = 1;
		    }
		  else//Neighbor is 1 smaller
		    {
		      a = EmTemp->which_neighbor(orig_neighbors+j*KEYLENGTH);//--case 5
		      if(a == 7) a = 0;
		      else if(a>=4) a = a-3;
		      b = a;
		      reg = 6;
		    }

		  for(k=0; k<KEYLENGTH; k++)
		    {		      		     
		      NewNeighbor[0][k]=*(Sons+a*KEYLENGTH+k);
		      NewNeighbor[1][k]=*(Sons+b*KEYLENGTH+k);
		    }
		  Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);  
		}

	      else//-- neighbor is refined
		{	      
		  if(NeighGeneration == (MyGeneration +1))//--case 4
		    {
		      
		      a=j;/*-- a: my son's number--*/
		      if(j>=4) a=j-3;
		      if(j==7) a=0;
		      
		      which_side=Neighbor->which_neighbor(Mykey); 
		      NdTemp = (Node*)(HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side+4)*KEYLENGTH));
		      NdTemp->putinfo(S_C_CON);
		      //NdTemp->putorder
		      
		      int b1 = which_side;
		      int b2=which_side+1;/*-- b1, b2: neighbor son's number--*/
		      if(b2==4) b2=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(Sons+a*KEYLENGTH+k);
		      reg = 6;//changed from 4
		      
		      unsigned* NeiSon = Neighbor->getson();

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b1*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b2*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);
		      
		    }

		  else if(NeighGeneration == MyGeneration )//-- case 2
		    {
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=j+1;
		      if(a==4) a=0;
		      b=which_side+1;
		      if(b==4) b=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+which_side*KEYLENGTH+k);
		      
		      Element* EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+a*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		  
		      for(k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+b*KEYLENGTH+k);

		      EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+j*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		    }

		  else if(NeighGeneration<MyGeneration)//-- case 3  
		    { 
		      if(j>=4) {
			printf("neigh_proc[%d] = %d\n",j,orig_neigh_proc[j]);
			j = j-4;
		      }
			
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=which_side;
		      if(which_side>=4) a=which_side-3;
		      if(which_side==7) a=0;

		      b = j+1;
		      if(b==4) b = 0;
		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeighSons+a*KEYLENGTH));

		      int edge_of_son;
		      if(a == which_side) edge_of_son = which_side;
		      else edge_of_son = which_side - 4;
		      
		      for(int k=0; k<KEYLENGTH; k++)
			{		      		     
			  NewNeighbor[0][k]=*(Sons+b*KEYLENGTH+k);
			  NewNeighbor[1][k]=*(Sons+j*KEYLENGTH+k);
			}
		      Neighbor->change_neighbor(&NewNeighbor[0][0], edge_of_son, myid, 3);  
		      
		    }		  		  
		}
	      
	    }//end of if same proc
	  
	  else if(*(orig_neigh_proc+j)!=myid && *(orig_neigh_proc+j) >= 0)
	    {
	      refined_new=new refined_neighbor();
	      refined_current->next=refined_new;
	      
	      NeighGeneration=*(EmTemp->get_neigh_gen()+j);

	      if(NeighGeneration==MyGeneration)
		{
		  assert(j<4);
		  refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH), 
					      Mykey, Sons, SIDE_SONS[j], MyGeneration, 1);
		}
	      
	      else
		refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH),
					    Mykey, Sons, &SIDE_SONS[j%4][j/4], MyGeneration, 2);//the neighbor is younger...

	      refined_current=refined_new;	    

	      /*for the neighbor information...assume that the neighbor was not 
		refined...if yes then it will be changed in update_inter*/
	      if(j<4)
		{
		  if(NeighGeneration==MyGeneration)
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(S_C_CON);
		      NdTemp->put_order(*(EmTemp->get_order()+j));
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(S_S_CON);
			}
		    }

		  else//NeighGen < MyGen
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(CORNER);
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(SIDE);
			  NdTemp->put_order(*(SonTemp->get_order()+j));
			}
		      
		    }
		}

	    }//end of is different proc
	  
	} 

    } 
  }

  MPI_Barrier(MPI_COMM_WORLD);
  update_neighbor_interprocessor(HT_Elem_Ptr, HT_Node_Ptr, refined_start, myid, nump);
  delete refined_start;
  //cout<<"ready with update: "<<myid<<"\n\n"<<flush;
}

#ifdef DISABLED
void update_neighbor_info(HashTable* HT_Elem_Ptr, Element* refined[], 
			  int count, int myid, int nump, 
			  HashTable* HT_Node_Ptr, int h_count )/*-h_count for  debugging-*/
{
  int son_order[2];
  Element* EmTemp;
  Element* Neighbor;
  Element* SonTemp;
  Node*    NdTemp;
  unsigned* orig_neighbors;
  int* orig_neigh_proc;
  int which_side;
  unsigned NewNeighbor[2][KEYLENGTH];
  
  int SIDE_SONS[4][2]={{0, 1}, {1, 2}, {2, 3}, {3, 0}};
 
  refined_neighbor* refined_start=new refined_neighbor();
  refined_neighbor* refined_current=refined_start;
  refined_neighbor* refined_new;

  //printf("H_ADAPT  IN  SUBDOMAIN %d  ELEMENTS REFINED: %d\n",myid, count);  


    unsigned* Sons;
    unsigned* NeighSons;
    int MyGeneration;
    int NeighGeneration;
    int NeighRefined;
    unsigned* Mykey;
    int reg;
    int a, b, k;
    

  for(int i=0; i<count; i++) {
      
    EmTemp=refined[i];//-- element ready for refinement


    orig_neighbors=EmTemp->get_neighbors();
    orig_neigh_proc=EmTemp->get_neigh_proc();
    Sons=EmTemp->getson();
    MyGeneration=EmTemp->get_gen();
    Mykey=EmTemp->pass_key();

    if(EmTemp->get_adapted_flag()==OLDFATHER) {
      EmTemp->put_adapted_flag(TOBEDELETED); //mark old father for deletion

      //update the neighbors of the new sons
      for(int j = 0; j<8; j++)
	{
	  if((*(orig_neigh_proc+j)>=0)&&(*(orig_neigh_proc+j)==myid)) //neighbor is at the same proc.
	    {
	      Neighbor = (Element*) HT_Elem_Ptr->lookup(orig_neighbors+j*KEYLENGTH);	
	      assert(Neighbor);
	      NeighGeneration=Neighbor->get_gen();
	      NeighRefined=Neighbor->get_refined_flag();

	      /*-------------------------------------------------
		  case 1: neighbor has the same gen. and was not refined in this step
	          case 2: neighbor has the same gen. and was refined in this step
		  case 3: neighbor is older and was refined
		  case 4: neighbor is younger and was refined
		  case 5: neighbor is younger and was not refined
		  other : wrong request
	       *-----------------------------------------------*/
	      if(!NeighRefined)
		{
		  assert(NeighGeneration>=MyGeneration);
		  which_side=Neighbor->which_neighbor(Mykey);
		  assert(which_side<4);

		  if(NeighGeneration == MyGeneration)//--case 1
		    {	
		      a=j+1;
		      if (a==4) a=0; 
		      b = j;
		      reg = 1;
		    }
		  else//Neighbor is 1 smaller
		    {
		      a = EmTemp->which_neighbor(orig_neighbors+j*KEYLENGTH);//--case 5
		      if(a == 7) a = 0;
		      else if(a>=4) a = a-3;
		      b = a;
		      reg = 6;
		    }

		  for(k=0; k<KEYLENGTH; k++)
		    {		      		     
		      NewNeighbor[0][k]=*(Sons+a*KEYLENGTH+k);
		      NewNeighbor[1][k]=*(Sons+b*KEYLENGTH+k);
		    }
		  Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);  
		}

	      else//-- neighbor is refined
		{	      
		  if(NeighGeneration == (MyGeneration +1))//--case 4
		    {
		      
		      a=j;/*-- a: my son's number--*/
		      if(j>=4) a=j-3;
		      if(j==7) a=0;
		      
		      which_side=Neighbor->which_neighbor(Mykey); 
		      NdTemp = (Node*)(HT_Node_Ptr->lookup(Neighbor->getNode()+(which_side+4)*KEYLENGTH));
		      NdTemp->putinfo(S_C_CON);
		      //NdTemp->putorder
		      
		      int b1 = which_side;
		      int b2=which_side+1;/*-- b1, b2: neighbor son's number--*/
		      if(b2==4) b2=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(Sons+a*KEYLENGTH+k);
		      reg = 6;//changed from 4
		      
		      unsigned* NeiSon = Neighbor->getson();

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b1*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);

		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeiSon+b2*KEYLENGTH));
		      Neighbor->change_neighbor(&NewNeighbor[0][0], which_side, myid, reg);
		      
		    }

		  else if(NeighGeneration == MyGeneration )//-- case 2
		    {
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=j+1;
		      if(a==4) a=0;
		      b=which_side+1;
		      if(b==4) b=0;

		      for(int k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+which_side*KEYLENGTH+k);
		      
		      Element* EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+a*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		  
		      for(k=0; k<KEYLENGTH; k++)
			NewNeighbor[0][k]=NewNeighbor[1][k]=*(NeighSons+b*KEYLENGTH+k);

		      EmSon= (Element*) HT_Elem_Ptr->lookup(Sons+j*KEYLENGTH);
		      
		      EmSon->change_neighbor(&NewNeighbor[0][0], j, myid, 6);//changed to 6 from 2 lend
		    }

		  else if(NeighGeneration<MyGeneration)//-- case 3  
		    { 
		      if(j>=4) {
			printf("neigh_proc[%d] = %d\n",j,orig_neigh_proc[j]);
			j = j-4;
		      }
			
		      which_side=Neighbor->which_neighbor(Mykey);
		      NeighSons=Neighbor->getson();

		      a=which_side;
		      if(which_side>=4) a=which_side-3;
		      if(which_side==7) a=0;

		      b = j+1;
		      if(b==4) b = 0;
		      Neighbor = (Element*)(HT_Elem_Ptr->lookup(NeighSons+a*KEYLENGTH));

		      int edge_of_son;
		      if(a == which_side) edge_of_son = which_side;
		      else edge_of_son = which_side - 4;
		      
		      for(int k=0; k<KEYLENGTH; k++)
			{		      		     
			  NewNeighbor[0][k]=*(Sons+b*KEYLENGTH+k);
			  NewNeighbor[1][k]=*(Sons+j*KEYLENGTH+k);
			}
		      Neighbor->change_neighbor(&NewNeighbor[0][0], edge_of_son, myid, 3);  
		      
		    }		  		  
		}
	      
	    }//end of if same proc
	  
	  else if(*(orig_neigh_proc+j)!=myid && *(orig_neigh_proc+j) >= 0)
	    {
	      refined_new=new refined_neighbor();
	      refined_current->next=refined_new;
	      
	      NeighGeneration=*(EmTemp->get_neigh_gen()+j);

	      if(NeighGeneration==MyGeneration)
		{
		  assert(j<4);
		  refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH), 
					      Mykey, Sons, SIDE_SONS[j], MyGeneration, 1);
		}
	      
	      else
		refined_new->set_parameters(*(orig_neigh_proc+j), (orig_neighbors+j*KEYLENGTH),
					    Mykey, Sons, &SIDE_SONS[j%4][j/4], MyGeneration, 2);//the neighbor is younger...

	      refined_current=refined_new;	    

	      /*for the neighbor information...assume that the neighbor was not 
		refined...if yes then it will be changed in update_inter*/
	      if(j<4)
		{
		  if(NeighGeneration==MyGeneration)
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(S_C_CON);
		      NdTemp->put_order(*(EmTemp->get_order()+j));
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(S_S_CON);
			}
		    }

		  else//NeighGen < MyGen
		    {
		      NdTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(j+4)*KEYLENGTH);
		      NdTemp->putinfo(CORNER);
		      for(int k=0; k<2; k++)
			{
			  SonTemp=(Element*)HT_Elem_Ptr->lookup(EmTemp->getson()+KEYLENGTH*SIDE_SONS[j][k]);
			  NdTemp=(Node*)HT_Node_Ptr->lookup(SonTemp->getNode()+(j+4)*KEYLENGTH);
			  assert(NdTemp);
			  NdTemp->putinfo(SIDE);
			  NdTemp->put_order(*(SonTemp->get_order()+j));
			}
		      
		    }
		}

	    }//end of is different proc
	  
	} 

    } 
  }

  MPI_Barrier(MPI_COMM_WORLD);
  update_neighbor_interprocessor(HT_Elem_Ptr, HT_Node_Ptr, refined_start, myid, nump);
  delete refined_start;
  //cout<<"ready with update: "<<myid<<"\n\n"<<flush;
}
#endif
