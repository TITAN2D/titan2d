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
 * $Id: move_data.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "../header/geoflow.h"

//#define PRINT_MOVE

void move_data(int numprocs, int myid, HashTable* El_Table, HashTable* NodeTable,
	       TimeProps* timeprops_ptr)
{

  if(numprocs<2) return;

#ifdef PRINT_MOVE
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid == 0)    printf("========================================================================================\n");
#endif
  int ibuck, iproc, inode, ierr, ineigh, ivar, ielem;
  /* assume that no elements share a neighboring element on another processor */
  int* num_send_recv = new int[numprocs];
  int* IfSendDone = new int[numprocs]; //initially holds current count of elements 
  //packed to send then "morphs" into whether or not the send has completed.
  int* IfRecvDone = new int[numprocs]; //says whether or not (1 or 0) the receive 
  //has completed, note if I'm not going to receive from a particular element 
  //then I'm already done receiving from them.
  MPI_Request* RequestSend= new MPI_Request[numprocs];
  MPI_Request* RequestRecv= new MPI_Request[numprocs];

  for(iproc=0;iproc<numprocs;iproc++)
    num_send_recv[iproc]=IfSendDone[iproc]=0;

  char filename[256];
  sprintf(filename,"move_data_elem.%04d",myid);
  FILE *fpelem1; 
  //if(timeprops_ptr->iter==279) fpelem1=fopen(filename,"w");
  sprintf(filename,"yada%d.txt",myid);  
  FILE *fpelem2; 
  //if(timeprops_ptr->iter==279) fpelem2=fopen(filename,"w");

  unsigned elemdebugkey2a[2]={2114123639,2004318068};
  Element* EmTemp=(Element*) El_Table->lookup(elemdebugkey2a);
  Node* NdTemp;

  /*
  if(timeprops_ptr->iter==49)
    if(EmTemp){
      if(EmTemp->get_adapted_flag()>=NOTRECADAPTED)
	ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,stdout);
    }
    else if(myid==2)
      ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,stdout);
  */

  int num_elem_on_proc=0, num_neighbors=0;
  /* count how many elements we should send and receive from other procs */
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr entryp;
  int *neigh_proc;
  int ifprint, numprocsrint=0; //for debug

  for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++) {

    entryp = *(buck+ibuck);
    while(entryp) {	
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;
      
      if((EmTemp->get_refined_flag()==0)&&
	 (EmTemp->get_adapted_flag()>0)
	 ) {
	//if this is a refined element don't involve!!!
	
	num_elem_on_proc++;
	
	neigh_proc = EmTemp->get_neigh_proc();
	ifprint=0;
	for(ineigh=0;ineigh<8;ineigh++) { 
	  if(neigh_proc[ineigh]>=0)
	    num_neighbors++;
	  
	  if((neigh_proc[ineigh]>=0)&&(neigh_proc[ineigh]!=myid)) {
	    num_send_recv[neigh_proc[ineigh]] += 1;
	    if((neigh_proc[ineigh]==0)||(neigh_proc[ineigh]==1))
	      ifprint=1;
	  }
	}
	if(0&&ifprint&&(timeprops_ptr->iter==279)){
	  fprintf(fpelem1,"___%d___==============================\n",numprocsrint);
	  ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),fpelem1);
	  fprintf(fpelem2,"%03d:  : {%10u,%10u}:",numprocsrint,
		  *(EmTemp->pass_key()+0),
		  *(EmTemp->pass_key()+1));
	  for(ineigh=0;ineigh<8;ineigh++)
	    if((neigh_proc[ineigh]>=0)&&(neigh_proc[ineigh]!=myid)&&
	       ((neigh_proc[ineigh]==0)||(neigh_proc[ineigh]==1))
	       ) 
	      fprintf(fpelem2," {%10u,%10u}",
		      *(EmTemp->get_neighbors()+ineigh*KEYLENGTH+0),
		      *(EmTemp->get_neighbors()+ineigh*KEYLENGTH+1));
	  fprintf(fpelem2,"\n");
	  
	  numprocsrint++;
	}//if(0&&ifprint&&(timeprops_ptr->iter==279))
	
      }//if((EmTemp->get_refined_flag()==0)&& ...
    }//while(entryp)
  }//for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++)

  /*
 if(timeprops_ptr->iter==279) {
    fclose(fpelem1);
    fclose(fpelem2);
  }
  */
  num_send_recv[myid] = 0;  // don't need to send info to myself
  





  //printf("myid=%d move_data() %d\n",myid,1); fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);  
  //printf("myid=%d move_data() %d\n",myid,2); fflush(stdout);
  

#ifdef NONONONO  
  //if(timeprops_ptr->iter==279) 
    {
      char FILENAME[256];
      sprintf(FILENAME,"num_send_receive%06d.%04d",timeprops_ptr->iter,myid);
      FILE *fpdebug;
      fpdebug=fopen(FILENAME,"a");
      fprintf(fpdebug,"myid=%d iter=%d time=%g [secs] num_elem_on_proc=%d num_neighbors=%d\n",
	      myid,timeprops_ptr->iter,timeprops_ptr->timesec(),num_elem_on_proc,num_neighbors);
      for(iproc=0;iproc<numprocs;iproc++)
	fprintf(fpdebug,"to/from proc %d: #elems=%d\n",
		iproc, num_send_recv[iproc]);
      fclose(fpdebug);
    }
#endif

  int send_tag = 22674; //original value Keith didn't change
  ElemPack** send_array = new ElemPack*[numprocs];
  ElemPack** recv_array = new ElemPack*[numprocs];

  for(iproc=0;iproc<numprocs;iproc++) {
    IfRecvDone[iproc]=!num_send_recv[iproc];
    if(num_send_recv[iproc]>0) {
      send_array[iproc] = new ElemPack[num_send_recv[iproc]];
      recv_array[iproc] = new ElemPack[num_send_recv[iproc]];      

      ierr= MPI_Irecv((void*) recv_array[iproc], num_send_recv[iproc],  
		      ELEMTYPE, iproc, send_tag+iproc, MPI_COMM_WORLD, 
		      &(RequestRecv[iproc]));
    }//if(num_send_recv[iproc] != 0)
  }//for(iproc=0;iproc<numprocs;iproc++)

  /* put (GHOST) elements to be moved in the proper arrays */
  for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++) {

    entryp = *(buck+ibuck);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;
      
      if((EmTemp->get_refined_flag()==0)&&
	 (EmTemp->get_adapted_flag()>0)
	 ) {
	//if this element should not be on this processor don't involve!!!
	
	neigh_proc = EmTemp->get_neigh_proc();
	for(ineigh=0;ineigh<8;ineigh++) {
	  iproc=neigh_proc[ineigh];
	  if((iproc!=myid) && (iproc>=0)) {
	    assert(IfSendDone[iproc]<num_send_recv[iproc]);
	    
	    Pack_element(EmTemp,(send_array[iproc]+IfSendDone[iproc]),NodeTable,myid);
	    
	    (send_array[iproc]+IfSendDone[iproc])->refined=GHOST;
	    (send_array[iproc]+IfSendDone[iproc])->adapted=-(EmTemp->get_adapted_flag());
	    
	    IfSendDone[iproc]++;
	  }//if((iproc!=myid) && (iproc>=0))
	}//for(ineigh=0;ineigh<8;ineigh++)
      }//if((EmTemp->get_refined_flag()==0)&& ...
    }//while(entryp) 
  }//for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++)
 
  for(iproc=0;iproc<numprocs;iproc++) {

    //better paranoid than dead
    assert(IfSendDone[iproc]==num_send_recv[iproc]); 

    IfSendDone[iproc]=!IfSendDone[iproc]; //"It's morphing time!" 
    //IfSendDone now holds whether or not (1 or 0) I'm done sending
    //to each processor, if I'm not going to send to a processor then
    //I'm already "done" sending to them.

    if(num_send_recv[iproc]>0) 
      //send ghost elements to processor iproc
      ierr=MPI_Isend((void*) send_array[iproc], num_send_recv[iproc], 
		     ELEMTYPE, iproc, send_tag+myid, MPI_COMM_WORLD, 
		     &(RequestSend[iproc]));    
  }
    

  //printf("myid=%d move_data() %d\n",myid,3); fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);  
  //printf("myid=%d move_data() %d\n",myid,4); fflush(stdout);
  

#ifdef PRINT_MOVE
  for(iproc=0;iproc<numprocs;iproc++)
    printf("proc %d is sending/receiving %d from %d\n", myid, num_send_recv[iproc], iproc);
#endif

  int NumNotRecvd, IfSentRecvd;
  MPI_Status status;  
  Element *elm, *new_elm;
  double not_used, *dPtr, *d2Ptr;
  int add_counter = 0, update_counter = 0;
    
  //wait for incomming data from each processor and incorporate each
  //processor's data as soon as I get it.
  do{
    NumNotRecvd=0;
    for(iproc=0;iproc<numprocs;iproc++) 
      if(!IfRecvDone[iproc]) {
	MPI_Test(&(RequestRecv[iproc]),&IfSentRecvd,&status);   

	if(IfSentRecvd) {
	  
	  for(ielem=0;ielem<num_send_recv[iproc];ielem++) {
	    elm = (Element*) (El_Table->lookup((recv_array[iproc]+ielem)->key));
	    if(elm == NULL) { // this elm doesn't exist on this proc
	      new_elm = new Element();

	      construct_el(new_elm, (recv_array[iproc]+ielem), NodeTable, myid, &not_used);
	      if((new_elm->get_adapted_flag()<0)&&
		 (new_elm->get_adapted_flag()>=-BUFFER))
		new_elm->put_myprocess(iproc);
	      El_Table->add(new_elm->pass_key(), new_elm);
	      add_counter++;
	    } //if(elm == NULL)
	    else{ 
	      //this elm is already on this proc, rather than delete old copy 
	      //and allocate space for a new one, save time by only copying the 
	      //new element data to the old element.
	      construct_el(    elm, (recv_array[iproc]+ielem), NodeTable, myid, &not_used);
	      if((elm->get_adapted_flag()<0)&&
		 (elm->get_adapted_flag()>=-BUFFER))
		elm->put_myprocess(iproc);
	      update_counter++;		
	    }//else
	  }//for(ielem=0;ielem<num_send_recv[iproc];ielem++)
	
	  IfRecvDone[iproc]=1;
	  delete [](recv_array[iproc]);
	}//if(IfSentRecvd)
	
	else
	  NumNotRecvd++;
      }//if(!IfRecvDone[iproc])
    
  }while(NumNotRecvd>0);
  delete []recv_array;
  delete []IfRecvDone;
  delete []RequestRecv;
  delete []num_send_recv;  
  
    
  //wait for sends to complete, delete sent arrays as soon as possible
  int NumNotSent;
  do{

    NumNotSent=0;
    for(iproc=0;iproc<numprocs;iproc++) 
      if(!IfSendDone[iproc]) {
	MPI_Test(&(RequestSend[iproc]),&IfSentRecvd,&status);   

	if(IfSentRecvd) {
	  IfSendDone[iproc]=1;
	  delete [](send_array[iproc]);
	}
	else
	  NumNotSent++;
      }
  }while(NumNotSent>0);
  delete []send_array;
  delete []IfSendDone;
  delete []RequestSend;


#ifdef PRINT_MOVE
  printf("proc %d has added %d  and updated %d ghost elements \n", myid, add_counter, update_counter);
#endif

  //shouldn't need this Barrier but "better paranoid than dead"
  MPI_Barrier(MPI_COMM_WORLD);  
  
  return;
    }


/* delete the ghost elements that were put in the element hashtable */
void delete_ghost_elms(HashTable* El_Table, int myid) 
{
  int ibuck;
  int delete_counter = 0;
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr entryp;
  Element *EmTemp;
  for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++) {
    
    entryp = *(buck+ibuck);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      assert(EmTemp);
      entryp=entryp->next;
	    
      if((EmTemp->get_refined_flag()==GHOST)//||
	 //((EmTemp->get_adapted_flag()<0)&&
	 //(EmTemp->get_adapted_flag()>=-BUFFER))
	 ) 
	{ //this is a GHOST element
	  EmTemp->void_bcptr();
	  El_Table->remove(EmTemp->pass_key(),1,stdout,myid,26);
	  delete EmTemp;
	  delete_counter++;
	}
    }//while(currentPtr)
  }//for(ibuck=0;ibuck<El_Table->get_no_of_buckets();ibuck++)
#ifdef PRINT_MOVE
  printf("proc %d has deleted %d ghost elms \n",myid, delete_counter);
#endif
  return;
}

void create_delete_memory() {

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  double* junk = new double[100000];
  junk[0] = 0;
  delete []junk;

  return;
}
