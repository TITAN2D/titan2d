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
 * $Id: repartition_BSFC.C 135 2007-06-07 20:15:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "./repartition_BSFC.h"

/* decent values for this partitioning scheme, the user can
   set them in the as parameters to tune for better performance */
#define BINS_PER_PROC 50 /* minimum amount of coarse bins on each processor */
#define SUBBINS_PER_BIN 50 /* amount of subbins a bin is divided into */
#define MAX_REFINEMENT_LEVEL 10 /* amount of refinement of the bins */

#define MIN_NUM_2_SEND 10 //if you're going to send fewer than this number of 
//elements to the processor before or after you on the space filling curve
//don't bother because the communication isn't worth it, this number is 
//machine dependent, should probably be much larger than 10

// stuff for the load-balancing weights used
#define NON_EMPTY_CELL    1.40 //the original and still best value
#define EMPTY_BUFFER_CELL 1.35 //Keith added
#define EMPTY_CELL        0.95 //cells not involved in computation
//#define EMPTY_CELL      1.00 //orignal, from before there were buffer cells

#define DEBUG_REPART2
#define DEBUG_REPART2C
//(mdj)#define DEBUG_ITER 179
#define DEBUG_ITER 1000000

int SequentialSend(int numprocs, int myid, 
		   HashTable* El_Table, HashTable* NodeTable,
		   TimeProps* timeprops_ptr, 
		   double *NewProcDoubleKeyBoundaries,int iseqsend);

void NonSequentialSendAndUpdateNeigh(int numprocs, int myid, 
				     HashTable* El_Table, HashTable* NodeTable,
				     TimeProps* timeprops_ptr, 
				     double *NewProcDoubleKeyBoundaries);



void BSFC_create_refinement_info(int* number_of_cuts, 
				 float* global_actual_work_allocated,
				 float total_weight, 
				 float* work_percent_array,
				 unstructured_communication verts_in_cuts_info,
				 float**, int myid, int numprocs);

void BSFC_create_bins(int num_local_objects,
		      BSFC_VERTEX_PTR sfc_vert_ptr, 
		      int* amount_of_bits_used, int size_of_unsigned,
		      float* global_actual_work_allocated, 
		      float *work_percent_array, float* total_weight_ptr,
		      int* balanced_flag, unstructured_communication* verts_in_cuts_info,
		      int* number_of_cuts,  
		      int bins_per_proc,
		      int myid, int numprocs);

void BSFC_update_element_proc(int myid, int numprocs, 
			      HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			      BSFC_VERTEX_PTR sfc_vert_ptr);

/* Space filling curve (BSFC) partioning routine */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void repartition(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int time_step)
{
  int ierr, i, j, k;                     /* local variables */
  int num_local_objects;              /* the number of objects this processor owns */
  BSFC_VERTEX_PTR sfc_vert_ptr;        /* array that stores the sfc objects */
  float* objs_wgt = NULL;             /* array of objects weights */
  float *global_actual_work_allocated = NULL; /* cumulative actual work allocated */
  float total_weight = 0;   /* sum of the work i.e. the amount of work for all objects */
  float *work_percent_array = NULL;   /* the cumulative percent of work each
					 processor should ideally get */
  int balanced_flag;                  /* flag to indicate if all partitions
					 are balanced */
  float* wgts_in_cut_ptr = NULL;      /* array of weights for sfc objects in 
					 a cut bin */
  int number_of_cuts = 0; /* maximum amount of cuts in a coarse bin on this processor */
  int amount_of_bits_used = 0;        /* amount of bits used in calculating the
					 bin an sfc object belongs to */
  int refinement_level_counter = 0;   /* counter to keep track of how many 
					 levels of bin refinement have been performed */
  int myid, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int size_of_unsigned = sizeof(unsigned);
  //printf("proc %d has entered the repartitioning scheme...\n",myid);
  
  /* get application data (number of objects, ids, weights, and coords */
  int elm_counter = 0;  //used to keep count how many active elements are on this processor
  int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
  HashEntryPtr entryp;
  Element* EmTemp;
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag()) {
	    if(*(EmTemp->get_state_vars()) > GEOFLOW_TINY)
	      EmTemp->put_lb_weight(NON_EMPTY_CELL);
	    else if(EmTemp->get_adapted_flag()==BUFFER)
	      EmTemp->put_lb_weight(EMPTY_BUFFER_CELL);
	    else
	      EmTemp->put_lb_weight(EMPTY_CELL);

	    EmTemp->put_myprocess(-1);
	    
	    total_weight += EmTemp->get_lb_weight();
	    EmTemp->put_new_old(BSFC_NEW);
	    elm_counter++;
	  }
	  entryp = entryp->next;
	}
    }
  //printf("there are %d active elements on proc %d\n",elm_counter,myid);
  //elements that share constrained nodes are grouped into 1 objects (since they cannot be separated onto different processors)  
  num_local_objects = 0; 
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag() && EmTemp->get_new_old() == BSFC_NEW) {
	    //check for constrained nodes on the vertex nodes
	    j = 4;
	    EmTemp->copy_key_to_lb_key();
	    num_local_objects++;
	  }
	  entryp = entryp->next;
	}
    }

  
  sfc_vert_ptr = (BSFC_VERTEX_PTR) malloc(num_local_objects * sizeof(BSFC_VERTEX));
  // fill up the sfc_vert_ptr array which stores all the necessary info about the load-balancing objects
  j = 0; 
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag() && EmTemp->get_new_old() > 0) {
	    sfc_vert_ptr[j].lb_weight = EmTemp->get_lb_weight();
	    unsigned* elem_key = EmTemp->get_lb_key();
	    for(k=0;k<KEYLENGTH;k++)
	      sfc_vert_ptr[j].sfc_key[k] = elem_key[k];
	 
	    j++;
	  }
	  entryp = entryp->next;
	}
    }    
  assert(j == num_local_objects);
  global_actual_work_allocated=(float*) malloc(sizeof(float)*(numprocs+1));

  work_percent_array = (float*) malloc(sizeof(float) * numprocs);
  unstructured_communication verts_in_cut_info;
  verts_in_cut_info.used_flag = 0;
  /*create bins, fill global weight vector and perform initial partition of the bins*/
  BSFC_create_bins(num_local_objects, sfc_vert_ptr, 
		   &amount_of_bits_used, size_of_unsigned,
		   global_actual_work_allocated, work_percent_array, 
		   &total_weight, &balanced_flag, &verts_in_cut_info,
		   &number_of_cuts, BINS_PER_PROC, myid, numprocs);
 
  if(balanced_flag != BSFC_BALANCED) { 
    int* local_balanced_flag_array; /* used to indicate which partitions on this 
				       processor are already balanced - useful
				       for when more than one cut in a bin */
    int* ll_bins_head; /* used to indicate the beginning of the linklist */
    float* work_prev_allocated = NULL; /* stores the weights of all 
					  objects before a cut */
    
    if(verts_in_cut_info.recv_count == 0 || myid == 0) 
      balanced_flag = BSFC_BALANCED;
    BSFC_create_refinement_info(&number_of_cuts, 
				global_actual_work_allocated, 
				total_weight, work_percent_array,
				verts_in_cut_info, &work_prev_allocated,
				myid, numprocs);
    
    ll_bins_head = (int*) malloc(sizeof(int) * (1+number_of_cuts));
    
    if(number_of_cuts == 0)
      balanced_flag = BSFC_BALANCED;
    
    if(ll_bins_head != NULL)
      ll_bins_head[number_of_cuts] = 0;  /* the first link list starts off
					    at array location number_of_cuts-1 ! */
    for(i=0;i<number_of_cuts;i++)
      ll_bins_head[i] = -1;
    
    local_balanced_flag_array = (int*) malloc(sizeof(int) * (1+number_of_cuts));
    
    for(i=0;i<number_of_cuts;i++)
      local_balanced_flag_array[i] = BSFC_BALANCED;
    local_balanced_flag_array[number_of_cuts] = balanced_flag;
    
    /* refine bins until a satisfactory partition tolerance is obtained */
    while(balanced_flag != BSFC_BALANCED &&
	  refinement_level_counter < MAX_REFINEMENT_LEVEL) { 
      BSFC_refine_partition(&balanced_flag, 
			    &amount_of_bits_used, verts_in_cut_info.recv_count,
			    verts_in_cut_info.recv_sfc_vert,  
			    work_percent_array, total_weight,
			    global_actual_work_allocated, 
			    number_of_cuts,
			    ll_bins_head, work_prev_allocated, 
			    SUBBINS_PER_BIN, local_balanced_flag_array, myid, numprocs);
      refinement_level_counter++;
    }
    // printf("proc %d has refined %d levels===============\n",myid,refinement_level_counter);
    free(local_balanced_flag_array);
    free(work_prev_allocated);
    free(ll_bins_head);
  }
  
  /* if the load-balancing objects were moved to different processors,
     we need to move them back now */
 
  if(verts_in_cut_info.used_flag != 0) {  // move back and free up the space...
    int recv_count = 0;
    sfc_vertex* send_sfc_vert = new sfc_vertex[verts_in_cut_info.send_count];  //this array is actually for receiving the info...
    // fill up the send array...
    int* proc_counter = new int[numprocs];
    proc_counter[0] = 0;
    for(i=1;i<numprocs;i++)
      proc_counter[i] = proc_counter[i-1] + verts_in_cut_info.send_procs_ptr[i-1];
    //done filling up the send array
    int tag = 21504;
    
    recv_count = 0;
    MPI_Request* send_request = new MPI_Request[numprocs];
    MPI_Request* recv_request = new MPI_Request[numprocs];
    for(i=0;i<numprocs;i++) {
      if(i!= myid) {
	//  receive necessary info here...
	if(verts_in_cut_info.send_procs_ptr[i] != 0) {
	  j = MPI_Irecv((send_sfc_vert+proc_counter[i]), verts_in_cut_info.send_procs_ptr[i],
			LB_VERT_TYPE, i, tag, MPI_COMM_WORLD, (send_request+i));
	  
	  
	}
	// send necessary info here...
	if(verts_in_cut_info.recv_procs_ptr[i] != 0) {
	  j = MPI_Isend(&(verts_in_cut_info.recv_sfc_vert[recv_count]), verts_in_cut_info.recv_procs_ptr[i],
			LB_VERT_TYPE, i, tag, MPI_COMM_WORLD, (recv_request+i));
	}
      }
      recv_count += verts_in_cut_info.recv_procs_ptr[i];
    }
    // wait until the info is sent and received...
    for(i=0;i<numprocs;i++)
      if(i!= myid)
	{
	  if(verts_in_cut_info.send_procs_ptr[i] != 0) {
	    MPI_Status status;
	    j = MPI_Wait(&(send_request[i]), &status);
	  }
	  if(verts_in_cut_info.recv_procs_ptr[i] != 0) {
	    MPI_Status status;
	    j = MPI_Wait(&(recv_request[i]), &status);
	  }
	}
	  
    recv_count = 0;
    for(i=0;i<myid;i++)
      recv_count += verts_in_cut_info.recv_procs_ptr[i];
    for(i=0;i<num_local_objects;i++)
      if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT) {
	if(sfc_vert_ptr[i].sfc_key[0] == (unsigned) 3157082553)
	  j = 0;
	if(sfc_vert_ptr[i].destination_proc != myid) {
	  j = sfc_vert_ptr[i].destination_proc;
	  sfc_vert_ptr[i].destination_proc = 
	    send_sfc_vert[proc_counter[sfc_vert_ptr[i].destination_proc]].destination_proc;
	  proc_counter[j] += 1;
	}
	else { // if i need to send to myself... 
	  sfc_vert_ptr[i].destination_proc = verts_in_cut_info.recv_sfc_vert[recv_count].destination_proc;
	  recv_count++;
	}
      }
    delete []proc_counter;
    delete []send_request;
    delete []recv_request;
    delete []send_sfc_vert;
  }
  free(global_actual_work_allocated);
  free(work_percent_array);  

  BSFC_update_element_proc(myid, numprocs, HT_Elem_Ptr, HT_Node_Ptr, sfc_vert_ptr);
    
  //debug stuff...
/*  unsigned* max = new unsigned[numprocs];
  for(i=0;i<numprocs;i++)
    max[i] = 0;
  unsigned* min = new unsigned[numprocs];
  for(i=0;i<numprocs;i++)
    min[i] = ~0;
  unsigned* ustore = new unsigned[numprocs];
  for(i=0;i<num_local_objects;i++) {
    if(max[sfc_vert_ptr[i].destination_proc] < sfc_vert_ptr[i].sfc_key[0])
      max[sfc_vert_ptr[i].destination_proc] = sfc_vert_ptr[i].sfc_key[0];
    if(min[sfc_vert_ptr[i].destination_proc] > sfc_vert_ptr[i].sfc_key[0])
      min[sfc_vert_ptr[i].destination_proc] = sfc_vert_ptr[i].sfc_key[0];
  }
  i = MPI_Allreduce(min, ustore, numprocs, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
  if(myid == 0) {
   // printf(" time step %d ****************** the minimums are: ", time_step);
    for(i=0;i<numprocs;i++)
      // printf("%u ",ustore[i]);
    //printf("\n");
  }
  delete []min;
  i = MPI_Allreduce(max, ustore, numprocs, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  if(myid == 0) {
    //printf(" time step %d ****************** the maximums are: ", time_step);
    for(i=0;i<numprocs;i++)
     // printf("%u ",ustore[i]);
    //printf("\n");
  }
  delete []max;
  delete []ustore;  */
  

  // more debugging stuff...
/*  double* proc_load = new double[numprocs];
  for(i=0;i<numprocs;i++)
    proc_load[i] = 0;
  //meshplotter(HT_Elem_Ptr, HT_Node_Ptr, time_step+10000);
  for(i=0;i<no_of_buckets;i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{
	  EmTemp = (Element*)(entryp->value);
	  if(!EmTemp->get_refined_flag()) {
	    if(EmTemp->get_new_old() > 0)
	      proc_load[EmTemp->get_myprocess()] += EmTemp->get_lb_weight();
	  }
	  entryp = entryp->next;
	}
    }
  double* procload2 = new double[numprocs];
  i = MPI_Allreduce(proc_load, procload2, numprocs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(myid == 0)
    for(i=0;i<numprocs;i++)
     // printf("proc %d load is %e total weight is %e\n",i, procload2[i], total_weight);
  delete []procload2;
  delete []proc_load;*/ 
  //done debug stuff
  free(sfc_vert_ptr);

  BSFC_update_and_send_elements(myid, numprocs, HT_Elem_Ptr, HT_Node_Ptr, time_step);

  
  return;
}

/* create info before starting the multi-level refinement of the bins 
 */

void BSFC_create_refinement_info(int* number_of_cuts, 
				 float* global_actual_work_allocated,
				 float total_weight, 
				 float* work_percent_array,
				 unstructured_communication verts_in_cut_info,
				 float** work_prev_allocated_ptr,
				 int myid, int numprocs)
{
  float my_work, work2;
  int i = 0, j;
  //printf("the total weight is %e\n",total_weight);
  /* find out how many cuts are in this bin. */
  my_work = global_actual_work_allocated[myid];
  while(my_work > total_weight * work_percent_array[myid-i]) {
    i++;
  }
  *number_of_cuts = i;
  if(verts_in_cut_info.recv_count == 0 || myid ==0) {
    *number_of_cuts = 0;
    return;
  }

  /* create link list for objects in the array.  link list
     is set up so that the objects in the array are
     traversed consecutively */
  for(i=0;i<(verts_in_cut_info.recv_count-1);i++)
    verts_in_cut_info.recv_sfc_vert[i].next_sfc_vert_index = i+1;
  verts_in_cut_info.recv_sfc_vert[verts_in_cut_info.recv_count-1].next_sfc_vert_index = -1;
     
  /* update work previously allocated to include work in all bins
     with higher keys than this bin */
  work2 = 0;
  for(i=0;i<verts_in_cut_info.recv_count;i++) 
    work2 += verts_in_cut_info.recv_sfc_vert[i].lb_weight;
  
  work2 = global_actual_work_allocated[myid] - work2;
  
  // note that work_prev_allocated data will only be correct for the
  // actual processors that are getting examined
  *work_prev_allocated_ptr = 
    (float*) malloc(sizeof(float)  * numprocs);
  for(i=myid-*number_of_cuts+1;i<=myid;i++)
    *(*work_prev_allocated_ptr+i) = work2;
  
  return;
}


//#define NEWCODEDOESNOTWORKYET
#ifdef NEWCODEDOESNOTWORKYET
void repartition2(HashTable* El_Table, HashTable* NodeTable,
		  TimeProps* timeprops_ptr) {
  repartition(El_Table,NodeTable,timeprops_ptr->iter);
  return;
}
#else
/*! Keith wrote this repartitioning function to make it work with any
 *  hash function, to not bother with constraining nodes (which is only
 *  useful in Continuous Galerkin method, titan uses a finite-difference/
 *  finite-volume predictor/corrector scheme), and to not bother with 
 *  keeping one brother of every element on the processor, now the 
 *  unrefinement Element constructor computes the key of the opposite 
 *  brother and that's all you need.  
 *
 *  The 2 purposes of repartition2() are to 
 *  1) balance the work between processors
 *  2) remove any overlap of key ranges caused by refinement. 
 *  with this implementation the second purpose is the more critical, but
 *  to understand why you first have to understand the problem.
 *
 *  In Titan elements and nodes are "stored" in hash tables, which is just 
 *  a one dimensional array of buckets and in each bucket is a linked list 
 *  of elements or nodes.  the "hash function" turns a "key" into the 
 *  indice of the bucket the element or node resides in and then you
 *  search the bucket's linked-list for the element/node with the key you
 *  are looking for.  Array access is fast, linked list access is slow.
 *  and there is an entire art form to minimizing the length of each linked 
 *  list, i.e. making sure the elements/nodes are as equally distributed 
 *  among the buckets as possible, in order to decrease the average time 
 *  it takes to retrieve and element/node, that is the true goal.
 *
 *  another complimentary way to achieve this goal is to ensure that 
 *  elements/nodes are "preloaded" in to cache, which can be achieved if
 *  elements/nodes close to each other in physical space are stored close 
 *  to each other in memory.  The way titan accomplishes this is to 
 *  organize data according to its position on (distance from the 
 *  beginning of) a space filling curve.  Essentially the keys are nothing
 *  more than the position of the node (or center node of the element) on
 *  the space filling curve.  A space filling curve is simply a curve that 
 *  travels to EVERY __POINT__ (not just every element or node) in 
 *  physical space, and visits all the points close to each other before 
 *  moving on and then never comes back to the same region.
 *
 *  basically it's a "power of 2" thing.  The normalized map is a unit 
 *  square ranging from (0,0) to (1,1).  if you divide this square into 
 *  4 sub squares (by dividing each dimension in half), the space filling 
 *  curve will visit all the points in one sub square before moving on to 
 *  the next sub square.  if you divide a sub square into 4 sub sub 
 *  squares the same holds true, and this relationship is infinitely 
 *  recursive all the way down to a single point in theory and down to 
 *  the last bit in the (currently) 8 byte key in practice.
 *
 *  however the physical dimensions of the map make it a rectangle not
 *  a square and we want elements/cells to be squares in PHYSICAL space 
 *  which means that each dimension of the map will be divided into a 
 *  different, and usually not a power of 2, integer number of elements.
 *  the problem arises when a "father" element is divided into it's 4
 *  "son" elements and some of the son elements are on a different "sub 
 *  square" or different "sub sub square" or different "sub sub sub (you 
 *  get the idea) square."  Since each processor owns one continuous 
 *  segment of the space filling curve, this means that refinement can 
 *  result in some of the "sons" having keys that should be on another 
 *  processor. 
 *
 *  why does it matter which processor an element belongs to?  Each 
 *  element needs its neighbor's information to update itself, which 
 *  during multiprocessor simulations means processors have to 
 *  communicate with each other, which means they have to know which 
 *  processors they need to send information to and receive information 
 *  from.  So if an element is on the wrong processor it's game over.  
 *  Luckily, elements "remember" which processors it's neighbors belong 
 *  to so this grants a _temporary_ reprieve but during repartitioning, 
 *  when elements are being moved from one processor to another, this 
 *  information needs to be correctly reset.  And it is a whole lot 
 *  easier and cheaper in terms of communication (which is slow and 
 *  hence you want to minimize it) for all the elements to belong to 
 *  the processors that owns the section of the space filling curve 
 *  they're on.
 *
 *  that is why it is absolutely essential to fix key range overlap 
 *  during repartitioning (at least for the way that I have implemented
 *  repartitioning).  A slight load imbalance can be tolerated but key 
 *  range overlap can not be.
 *
 *  repartition2 has 2 steps a sequential send (sending/receiving 
 *  elements to/from the processor(s) immediately before and/or after
 *  you on the space filling curve) and a non sequential send that 
 *  fixes any remaining key range overlap and then updates the 
 *  neighbor information of every element it owns.
 *
 *  I (Keith) implemented the sequential send in an "intelligent" way
 *  or at least intelligent enough that it can actually be "confused"
 *  by a pathological case.  The "sequential send" determines how many
 *  elements it would need to send and receive from its 2 neighbors on 
 *  the space filling curve to
 *  1) achieve load (computational work) balance
 *  2) to eliminate it's key range overlap with its neighbor
 *  since communication is expensive it is preferable to do a one way
 *  only send/receive, that is to send OR receive enough elements to 
 *  fix BOTH load balance and key range overlap in just one send OR 
 *  receive.  If you have to sacrifice a little load balance to ensure 
 *  the key range overlap is fixed that's okay because the slight load 
 *  imbalance will be corrected by the next repartitioning so it's not 
 *  a big deal.
 *
 *  the pathological case is when a processor has to send away _all_ of 
 *  its elements to fix the key range overlap, and possibly give away
 *  the same element(s) to BOTH of its neighbor on the space filling
 *  curve. Yes that actually happened and it caused titan to crash, 
 *  which is why I had to rewrite repartition2 to build in a failsafe to 
 *  protect against that.  The failsafe is to make each processor refuse 
 *  to send away more than half (actually total number of elements minus 
 *  one divided by 2) of its elements to either neighbor, and if then if 
 *  it is necessary, repeat the sequential send until each processor's 
 *  "maximum key" is greater than its "minimum key".  This is the reason 
 *  the sequential send is inside a while loop.  Once the maximum key is 
 *  greater than than the minimum key, the non sequential send can fix 
 *  the remaining key range overlap without any difficulty.
 *
 *  note because of some tricks I've played with the initial grid 
 *  generator (to exploit the space filling curve's "power of 2" effect)
 *  the sequential send will usually occur only once per repartitioning
 *  and the non sequential send will usually not occur at all.  In fact
 *  these will only occur when there are far too few elements per 
 *  processor OR when the exceedingly vast majority of the elements are
 *  very close together on the space filling curve and very few exist 
 *  elsewhere.  This means communication will usually be minimal and 
 *  thus the code will be fast.
 *
 *  I (Keith) spent a minor amount of work to make this fast (or at 
 *  least not "as dumb as a post" slow), by doing the intelligent 
 *  (usually one way) sequential send and overlapping computation and 
 *  communication to keep the CPU busy while it's waiting to exchange 
 *  elements with other processors so no time will be wasted then.
 *
 *  The sequential send and non sequential send do all of their own 
 *  communication (do not rely on other non MPI functions to do it)
 *  but a move_data() is required immediately after repartitioning to 
 *  create the layer of "ghost" cells/elements around the processor's 
 *  collection of elements.
 *
 */
void repartition2(HashTable* El_Table, HashTable* NodeTable,
		  TimeProps* timeprops_ptr) {

  int myid, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(!(numprocs>1)) return;

  double *NewProcDoubleKeyBoundaries=CAllocD1(numprocs+1);

  int counter=0;
  int ifrepeat=0;
  do{
    ifrepeat=
      SequentialSend(numprocs, myid, El_Table, NodeTable, timeprops_ptr, 
		     NewProcDoubleKeyBoundaries, counter);
    //printf("myid=%d counter=%d MyFirstAndLastDoubleKey=%g %g\n",myid,counter,MyFirstAndLastDoubleKey[0],MyFirstAndLastDoubleKey[1]);
    if(counter++>4*numprocs){
      printf("myid=%d repartition2() did %d SequentialSend()'s when there were only %d processors.  This means we are probably in an infinite loop, so we are aborting\n",myid,counter,numprocs);
      assert(0);
    }
    
  }while(ifrepeat);

  //printf("myid=%d counter=%d\n",myid,counter);

  NonSequentialSendAndUpdateNeigh(numprocs, myid, El_Table, NodeTable, 
				  timeprops_ptr, NewProcDoubleKeyBoundaries);

  CDeAllocD1(NewProcDoubleKeyBoundaries);

  return;
}
#endif


void checkelemnode(HashTable *El_Table, HashTable *NodeTable, 
		   int myid, FILE *fpdebug, double loc){
  unsigned elemdebugkey2a[2]={ 695804849, 991146299};
  unsigned nodedebugkey2a[2]={ 695852110,3303820997};

  return;

  if(fpdebug==NULL) fpdebug=stdout;

  fprintf(fpdebug,"**************************\nmyid=%d location=%g\n",
	  myid,loc);
  ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,fpdebug);
  NodeBackgroundCheck(El_Table,NodeTable,nodedebugkey2a,fpdebug);

  fflush(fpdebug);

  return;
}
int SequentialSend(int numprocs, int myid, 
		   HashTable* El_Table, HashTable* NodeTable,
		   TimeProps* timeprops_ptr, 
		   double *NewProcDoubleKeyBoundaries, int iseqsend) {

  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  int ibuck, ielem, ikey, iproc, ineigh, ierr;
  Element *EmTemp;
  unsigned nullkey[2]={0,0};
  double MyFirstAndLastDoubleKey[2]={-1.0,-2.0};



  char fname2[256];
  sprintf(fname2,"seqsend%04d.debug",myid);
  FILE *fpdb2; 
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"w");
    fprintf(fpdb2,"myid=%d iter=%d iseqsend=%d\n",myid,timeprops_ptr->iter,iseqsend);
    fclose(fpdb2);
    fpdb2=fopen(fname2,"a");
  }

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 0.0);
  /*
  if((timeprops_ptr->iter==279)&&((myid==0)||(myid==1))) {
    printf("checking proc myid=%d for element ={%10u,%10u}\n",
	   myid,elemdebugkey2a[0],elemdebugkey2a[1]);
    ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,stdout);
  }
  */
  //unsigned debugkey1[2]={2111546026,2863311530};
  //unsigned debugkey2[2]={2111570602,2863311530};
  //unsigned debugkey3[2]={2111566506,2863311530};
  //unsigned debugkey4[2]={2111563776,         0};

  /*
  if(myid==1) MPI_Barrier(MPI_COMM_WORLD);
  if(timeprops_ptr->iter==9){
    printf("***********************************************\n");
    printf("myid=%d ***************************************\n",myid);
    printf("***********************************************\n");
    ElemBackgroundCheck(El_Table,NodeTable,debugkey1,stdout);
    ElemBackgroundCheck(El_Table,NodeTable,debugkey2,stdout);
    ElemBackgroundCheck(El_Table,NodeTable,debugkey3,stdout);
    ElemBackgroundCheck(El_Table,NodeTable,debugkey4,stdout);      
    printf("***********************************************\n");
    NodeBackgroundCheck(El_Table,NodeTable,debugkey1,stdout);
    NodeBackgroundCheck(El_Table,NodeTable,debugkey2,stdout);
    NodeBackgroundCheck(El_Table,NodeTable,debugkey3,stdout);
    NodeBackgroundCheck(El_Table,NodeTable,debugkey4,stdout);          
  }
  if(myid==0) MPI_Barrier(MPI_COMM_WORLD);
  if(myid==1){
    printf("***********************************************\n");
    printf("***********************************************\n");
    printf("***********************************************\n");
  }
  */

  /*
  Element* Curr_El=(Element*) El_Table->lookup(elemdebugkey2a);
  
  printf("myid=%d iter=%d\n",myid,timeprops_ptr->iter);
  if(timeprops_ptr->iter==49) {
    if(Curr_El){
      if((Curr_El->get_adapted_flag()>=NOTRECADAPTED)||(myid==2)) {
	printf("myid=%d repartition2\n",myid);
	ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,stdout);
      }
    }
    else if(myid==2) {
      printf("myid=%d repartition2\n",myid);
      ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,stdout);
    }
  }
  */

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 1.0\n",myid);
    fclose(fpdb2);
  }
#endif

  int num_elem = 0;
  for(ibuck=0; ibuck<num_buck; ibuck++) {
    currentPtr = *(buck+ibuck);

    while(currentPtr){
	
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;     
      assert(EmTemp);
      
      if(EmTemp->get_adapted_flag()>=NOTRECADAPTED){
	//this is an active element

	//weight the active elements for repartitioning
	if(*(EmTemp->get_state_vars())>GEOFLOW_TINY)
	  //this has pile
	  EmTemp->put_lb_weight(NON_EMPTY_CELL);
	else if(EmTemp->get_adapted_flag()==BUFFER)
	  //this might have pile before the next adaptation
	  EmTemp->put_lb_weight(EMPTY_BUFFER_CELL);
	else
	  //this will not have pile before the next adaptation
	  //but might have pile before the next repartitioning
	  EmTemp->put_lb_weight(EMPTY_CELL);
	
	EmTemp->put_new_old(BSFC_NEW); //legacy from "original" 
	//repartition(), I don't know the consequences of not 
	//having it so included it just to be safe

	num_elem++; //count the active elements
      }
      else{
	//delete the non active elements
	EmTemp->void_bcptr();
	El_Table->remove(EmTemp->pass_key(),1,stdout,myid,1);
	delete EmTemp;
      }
    }//while(currentPtr)
  }//for(ibuck=0; ibuck<num_buck; ibuck++)

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 2.0\n",myid); 
    fclose(fpdb2);
  }
#endif

  double *LoadBalArray=(double*)   malloc(num_elem*sizeof(double));
  double *DoubleKeyArray=(double*)   malloc(num_elem*sizeof(double));
  Element **ElemArray=(Element**) malloc(num_elem*sizeof(void* ));
  double doublekeyrange1=*(El_Table->get_doublekeyrange()+1);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 at 3.0\n",myid);
    fclose(fpdb2);
  }

#endif

  //store pointers to all the active elements in an array and store
  //double precision versions of the elements' keys in an array
  ielem=0;
  for(ibuck=0; ibuck<num_buck; ibuck++) {
    currentPtr = *(buck+ibuck);
    
    while(currentPtr) {	
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;     
      assert(EmTemp);
      
      if(EmTemp->get_adapted_flag()>=NOTRECADAPTED) {
	DoubleKeyArray[ielem]=*(EmTemp->pass_key()+0);
	for(ikey=1;ikey<KEYLENGTH;ikey++)
	  DoubleKeyArray[ielem]=DoubleKeyArray[ielem]*doublekeyrange1+
	    *(EmTemp->pass_key()+ikey);
	
	ElemArray[ielem]=EmTemp;
	
	ielem++;
      }
    }//while(currentPtr)
  }//for(ibuck=0; ibuck<num_buck; ibuck++)
  
  assert(ielem==num_elem); //sanity check

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 at 4.0\n",myid); 
    fclose(fpdb2);
  }
#endif

  //quicksort into ascending order of keys, the array of pointers to
  //elements
  q_sort_data(DoubleKeyArray,(void **) ElemArray, 0, num_elem-1);
  
  /*
  if(timeprops_ptr->iter==99&&myid==2) {
    printf("last element on proc myid=%d is {%10u,%10u}\n",myid,
	   *(ElemArray[num_elem-1]->pass_key()+0),
	   *(ElemArray[num_elem-1]->pass_key()+1));
    ElemBackgroundCheck(El_Table,NodeTable,ElemArray[num_elem-1]->pass_key(),stdout);
  }
  */

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 1.0);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 5.0\n",myid); 
    fclose(fpdb2);
  }
#endif

  unsigned MyFirstAndLastKey[2*KEYLENGTH];

  //store my first key
  for(ikey=0;ikey<KEYLENGTH;ikey++)
    MyFirstAndLastKey[0*KEYLENGTH+ikey]=
      *(ElemArray[0]->pass_key()+ikey);

  //DoubleArray now becomes the cumulative sum of the load balancing
  //weights (lb_weight);
  LoadBalArray[0]=ElemArray[0]->get_lb_weight();
  for(ielem=1;ielem<num_elem;ielem++)
    LoadBalArray[ielem]=LoadBalArray[ielem-1]+
      ElemArray[ielem]->get_lb_weight();

  //store my last key, do this after finding the cumulative lb_weight
  //sum instead of before it, because I know ElemArray[num_elem-1] 
  //will be in cache now because I just used him, so save a SMALL 
  //amount of time  
  for(ikey=0;ikey<KEYLENGTH;ikey++)
    MyFirstAndLastKey[1*KEYLENGTH+ikey]=
      *(ElemArray[num_elem-1]->pass_key()+ikey);

  /*
  if(myid==1) MPI_Barrier(MPI_COMM_WORLD);
  printf("myid=%d A) MyFirstAndLastKey={{%10u,%10u},{%10u,%10u}}\n",myid,
	 MyFirstAndLastKey[0],MyFirstAndLastKey[1],
	 MyFirstAndLastKey[2],MyFirstAndLastKey[3]);
  if(myid==0) MPI_Barrier(MPI_COMM_WORLD);
  */

  double *cum_sum_proc_lb_weight=CAllocD1(numprocs);
  double *temprecvdouble=CAllocD1(3*numprocs);
  double **AllFirstAndLastDoubleKeys=CAllocD2(numprocs,2);


  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 2.0);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 6.0\n",myid);
    fclose(fpdb2);
  }
#endif

  double tempsenddouble[3]=
    {(double) MyFirstAndLastKey[0],
     (double) MyFirstAndLastKey[KEYLENGTH],
     LoadBalArray[num_elem-1]};

  for(ikey=1;ikey<KEYLENGTH;ikey++) {
    tempsenddouble[0]=tempsenddouble[0]*doublekeyrange1+
      MyFirstAndLastKey[          ikey];

    tempsenddouble[1]=tempsenddouble[1]*doublekeyrange1+
      MyFirstAndLastKey[KEYLENGTH+ikey];
  }

  

  // **********************************
  // ** section for debugging output **
  // **********************************
#ifdef DEBUG_REPART2
  //MPI_Barrier(MPI_COMM_WORLD);
  //for(iproc=0;iproc<numprocs-1;iproc++)
  //if(myid>iproc) MPI_Barrier(MPI_COMM_WORLD);

  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 MyF&LKey={{%10u,%10u},{%10u,%10u}} doublekeyrange1=%20.14g tempsenddouble={%20.14g,%20.14g,%20.14g}\n",myid,
	    MyFirstAndLastKey[0],MyFirstAndLastKey[1],
	    MyFirstAndLastKey[2],MyFirstAndLastKey[3],
	    doublekeyrange1,
	    tempsenddouble[0],tempsenddouble[1],tempsenddouble[2]); fflush(stdout);
    
    fclose(fpdb2);
  }


  //for(iproc=1;iproc<numprocs;iproc++)
  //if(myid<iproc) MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
#endif

  MPI_Allgather(tempsenddouble,3,MPI_DOUBLE,
		temprecvdouble,3,MPI_DOUBLE,MPI_COMM_WORLD); 
  
#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 at 7.0\n",myid);
    fclose(fpdb2);
  }

#endif
  
  for(iproc=0;iproc<numprocs;iproc++) {
    AllFirstAndLastDoubleKeys[iproc][0]=temprecvdouble[iproc*3+0];
    AllFirstAndLastDoubleKeys[iproc][1]=temprecvdouble[iproc*3+1];
    cum_sum_proc_lb_weight[   iproc]   =temprecvdouble[iproc*3+2];
  }
  CDeAllocD1(temprecvdouble);

  //used in first-send criteria 1
  for(iproc=1;iproc<numprocs;iproc++)
    cum_sum_proc_lb_weight[iproc]+=cum_sum_proc_lb_weight[iproc-1];

  double lb_weight_per_proc=cum_sum_proc_lb_weight[numprocs-1]/numprocs;
  double mynewstart=lb_weight_per_proc*myid;
  double mynewstop=mynewstart+lb_weight_per_proc;


  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 3.0);

  // **********************************
  // ** section for debugging output **
  // **********************************
#ifdef DEBUG_REPART2
  //MPI_Barrier(MPI_COMM_WORLD);
  //for(iproc=0;iproc<numprocs-1;iproc++)
  //if(myid>iproc) MPI_Barrier(MPI_COMM_WORLD);

  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d F&LDKeys=\n",myid);
    for(iproc=0;iproc<numprocs;iproc++)
      fprintf(fpdb2,"    iproc=%d {%20.14g,%20.14g}\n",iproc,
	      AllFirstAndLastDoubleKeys[iproc][0], 
	      AllFirstAndLastDoubleKeys[iproc][1]);
    fclose(fpdb2);
  }
  //fflush(stdout);

  //for(iproc=1;iproc<numprocs;iproc++)
  //if(myid<iproc) MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  

  //printf("myid=%d repart2 at 8.0\n",myid); fflush(stdout);
#endif

  // *********************************************************************
  // ** Because the keys are generated from normalized (0 to 1) x and y **
  // ** coordinates, there is only 1 space filling curve ever used for  **
  // ** titan, regardless of the physical size of the domain. However,  **
  // ** the number and position on the space filing curve of points that**
  // ** represent nodes will be different from simulation to simulation **
  // ** Let DX be the entire length of the physical domain in the X     **
  // ** direction and DY be the entire length of the physical domain in **
  // ** the Y direction.  Further let dX=DX/2^N and dY=DY/2^N.  The     **
  // ** nature of the space filling curve is such that every subdomain  **
  // ** that is dX by dY long, and starts and ends at integer multiples **
  // ** of dX and dY (from the lower left corner of the domain) will    **
  // ** contain exactly one (continuous) segment of the space filling   **
  // ** curve.  Any subdomain that does not start and end at integer    **
  // ** multiples of dX and dY (from the lower left corner of the       **
  // ** domain) will contain multiple discontinous segments of the      **
  // ** space filling curve.  Consequently sufficient refinement of     **
  // ** elements at boundaries of such (discontinuous segment)          **
  // ** subdomains will result in son elements whose keys should belong **
  // ** to the next subdomain over. In other words the key ranges of    **
  // ** discontinuous segment subdomains will overlap, which can cause  **
  // ** all kinds of havoc.  On two processor runs this isn't so bad,   **
  // ** you just need 2 criteria for determining which elements to send **
  // ** to your only neighbor.  The first criteria based on load        **
  // ** balancing weight, the second criteria based on how many of your **
  // ** keys overlap with your neighbor key range.  This "first-send"   **
  // ** will fix the key boundaries of a single segment of the space    **
  // ** filling curve that belongs on each processor.  However, on more **
  // ** than 2 processors, you will end up with the situation where     **
  // ** refinement creates elements that belong, according to their     **
  // ** keys, to processors that can be located anywhere, i.e. not      **
  // ** immediately before or after you, on the space filling curve.    **
  // ** This necessesitates a "second-send" of elements to arbirtrary   **
  // ** processors, and takes up a lot of communication time.  The good **
  // ** news is there is a simple initial gridding trick you can play   **
  // ** that will minimize the ammount of non-sequential repartitioning **
  // ** if you have numprocs processors then divide your 2D normalized  **
  // ** grid into 4^ceil(log4(numprocs)) square subdomains, that is     **
  // ** subdomains that are dX=DX/2^N by dY=DY/2^N in dimensional       **
  // ** length. As we said above these subdomains will have exactly one **
  // ** (continuous) segment of the space filling curve, and each       **
  // ** processor will own a whole subdomain's segment of space filling **
  // ** curve (which will not have key ranges that over lap with other  **
  // ** processors) or it will own parts of the segment of space        **
  // ** filling curve contained by one subdomain shared with other      **
  // ** processors.  In this second case there can be overlapping key   **
  // ** ranges, but the key ranges will only overlap with other         **
  // ** processors who share the same subdomain with you, and hence are **
  // ** located very close to you on the space filling curve, which as  **
  // ** a result minimizes non sequential interprocessor repartitioning **
  // ** Just to be clear, the repartitioning code must handle non       **
  // ** sequential processor repartitioning (hence the "second-send")   **
  // ** to be failsafe, but the initial gridding trick will restore the **
  // ** lost performance to you.  The moral of the story is to not      **
  // ** screw with the preprocessor/grid generator I (Keith) wrote, if  **
  // ** you don't want to severely degrade your performance.  You have  **
  // ** been warned!!!                                                  **
  // *********************************************************************




  // *********************************************************************
  // ** determine the number of elements you need to repartition during **
  // ** the "first-send" according to both criteria                     **
  // *********************************************************************

  //{-1,num_elem} initialization provides failsafe, this failsafe 
  //is ALWAYS used for myid==0 and myid==numprocs-1
  int isend1[2]={-1,num_elem}; //criteria 1 of first-send
  int isend2[2]={-1,num_elem}; //criteria 2 of first-send

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");


    fprintf(fpdb2,"myid=%d repart2 A) isend1={%d,%d}\n",myid,isend1[0],isend1[1]); 
    fclose(fpdb2);
  }
#endif

  if(myid>0)
    for(ielem=0;ielem<num_elem;ielem++)
      LoadBalArray[ielem]+=cum_sum_proc_lb_weight[myid-1];

  CDeAllocD1(cum_sum_proc_lb_weight);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 LoadBalArray[%d]=%20.14g\n",
	 myid,num_elem-1,LoadBalArray[num_elem-1]);
    fclose(fpdb2);
  }
  //printf("myid=%d repart2 at 9.0\n",myid); fflush(stdout);
#endif

  // **********************************************************************
  // find, according to criteria 1, the first element to send the processor 
  // after me. do this before finding last element to send to process 
  // before me so in case of conflict more work goes to higher rank 
  // processes because rank 0 process (myid==0) has a little extra work 
  // because it outputs more data
  // **********************************************************************
  if(myid<numprocs-1)
    for(ielem=8;ielem<num_elem;ielem++) 
      //starts at ielem=8 so I keep a minumum of 8 elements for myself
      if(LoadBalArray[ielem]<mynewstop)
	isend1[1]=ielem+1;
      else
	break;

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 isend1[1]=%d ielem=%d num_elem=%d\n",
	 myid,isend1[1],ielem,num_elem);
    fclose(fpdb2);
  }
#endif

  //find the last element to send the processor before me
  if(myid>0)
    for(ielem=0;ielem<isend1[1]-8;ielem++)
      //stops at isend1[1]-8 so I keep a minumum of 8 elements for myself
      if(LoadBalArray[ielem]>=mynewstart) {
	isend1[0]=ielem-1;
	break; 
      }

  // **********************************
  // ** section for debugging output **
  // **********************************
#ifdef DEBUG_REPART2
  //MPI_Barrier(MPI_COMM_WORLD);
  //for(iproc=0;iproc<numprocs-1;iproc++)
  //if(myid>iproc) MPI_Barrier(MPI_COMM_WORLD);

  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 B) isend1={%d,%d}\n",myid,isend1[0],isend1[1]); 
    fclose(fpdb2);
  }


  //for(iproc=1;iproc<numprocs;iproc++)
  //if(myid<iproc) MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
#endif


  // **********************************************************************
  // find, according to criteria 2, the number of elements that overlap 
  // with the processors before and after you, you will either have to send
  // or receive the overlapping elements.
  // **********************************************************************
  double doublekey;
  if(myid<numprocs-1)
    for(ielem=0;ielem<num_elem;ielem++) {

      //compute the double precision equivalent of this element's key
      doublekey=*(ElemArray[ielem]->pass_key()+0);
      for(ikey=1;ikey<KEYLENGTH;ikey++)
	doublekey=doublekey*doublekeyrange1+
	  *(ElemArray[ielem]->pass_key()+ikey);
      
      if(doublekey>=AllFirstAndLastDoubleKeys[myid+1][0]) {
	isend2[1]=ielem;
	break;
      }
    }

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 4.0);


#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 at 10.0\n",myid);
    fclose(fpdb2);
  }

#endif

  if(myid>0) {
    for(ielem=0;ielem<num_elem;ielem++) {

      //compute the double precision equivalent of this element's key
      doublekey=*(ElemArray[ielem]->pass_key()+0);
      for(ikey=1;ikey<KEYLENGTH;ikey++)
	doublekey=doublekey*doublekeyrange1+
	  *(ElemArray[ielem]->pass_key()+ikey);
      
      if(doublekey>AllFirstAndLastDoubleKeys[myid-1][1]) {
	isend2[0]=ielem-1;
	break;
      }
    }
    if(ielem==num_elem) isend2[0]=num_elem-1;
  }

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 C) isend2={%d,%d}\n",myid,isend2[0],isend2[1]); 
    fclose(fpdb2);
  }

#endif

  // **********************************************************************
  // now combine criteria 1 and 2 to decide how many elements will I send 
  // or receive from the processers before/after me on the space filling 
  // curve, the default is to not send anything to keep communication cost 
  // down.
  // num_send[0] is the number of elements you send to the processor 
  //             before you on the space filling curve
  // num_send[1] is the number of elements you send to the processor
  //             after you on the space filling curve
  // num_recv[0] is the number of elements you receive from the processor 
  //             before you on the space filling curve
  // num_recv[1] is the number of elements you receive from the processor
  //             after you on the space filling curve
  // **********************************************************************
  int NumMy2WayElem=(num_elem-1)/2, NumYour2WayElem[2]={-1,-1};
  int num_send[ 2]={          0,                 0}; //final
  int num_send1[2]=//criteria 1 capped as failsafe
    {(isend1[0]+1       <NumMy2WayElem)?isend1[0]+1       :NumMy2WayElem,
     (num_elem-isend1[1]<NumMy2WayElem)?num_elem-isend1[1]:NumMy2WayElem};
  int num_send2[2]={isend2[0]+1,num_elem-isend2[1]}; //criteria 2
  int num_sendmax[2]= //max of criteria 1 and criteria 2
    {(num_send1[0]>num_send2[0])?num_send1[0]:num_send2[0],
     (num_send1[1]>num_send2[1])?num_send1[1]:num_send2[1]};
  
  int num_recv[   2]={0,0}; //final
  int num_recv1[  2]={0,0}; //criteria 1
  int num_recv2[  2]={0,0}; //criteria 2
  int num_recvmax[2]={0,0}; //max of criteria 1 and criteria 2
  
  //keep cheap communication simple, this short blocking send won't
  //cost much time because the processors should be synced since I 
  //just did a (necessary) MPI_Allgather() a little before now
  int tempsendint[5]={num_send1[0],num_send1[1],num_send2[0],num_send2[1],num_elem};
  int *temprecvint=CAllocI1(5*numprocs);
  int num_neigh_elem[2]={-1,-1};
  int If2WaySend[2]={0,0};


  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 5.0);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 num_send1={%d,%d} num_send2={%d,%d}\n",myid,
	    num_send1[0],num_send1[1],num_send2[0],num_send2[1]); 
    fclose(fpdb2);
  }

#endif

  MPI_Allgather(tempsendint,5,MPI_INT,
		temprecvint,5,MPI_INT,MPI_COMM_WORLD); 

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 at 11.0\n",myid);
    fclose(fpdb2);
  }

#endif

  //figure out how many elements I send-to/receive-from the processor 
  //before me on the space filling curve, the default is zero
  if(myid>0){
    num_recv1[0]=temprecvint[(myid-1)*5+1];
    num_recv2[0]=temprecvint[(myid-1)*5+3];
    num_recvmax[0]=(num_recv1[0]>num_recv2[0])?num_recv1[0]:num_recv2[0];
    num_neigh_elem[0]=temprecvint[(myid-1)*5+4];
    //NumYour2WayElem[0]=(num_neigh_elem[0]-1)/2;


    if(((num_send1[0]<MIN_NUM_2_SEND)&&
	(num_send2[0]==0            )
	)||
       (((num_recv1[0]>=MIN_NUM_2_SEND)||
	 (num_recv2[0]>0))&&
	(num_recv1[0]>=num_recv2[0]))
       ){
      num_send[0]=0;
      if((num_recv1[0]>=MIN_NUM_2_SEND)||
	 (num_recv2[0]>0))
	num_recv[0]=num_recv1[0];
      else
	num_recv[0]=0;
    }
    else if(num_send1[0]>=num_send2[0]){
      num_send[0]=num_send1[0];
      num_recv[0]=0;
    }
    else{
      /*
      num_send[0]=(num_send2[0]<NumMy2WayElem[0])?num_send2[0]:
	NumMy2WayElem[0];
      num_recv[0]=(num_recv2[0]<NumYour2WayElem[0])?num_recv2[0]:
	NumYour2WayElem[0];
      */
      If2WaySend[0]=1;
    }


    /*
    //"<" not "<=" in (num_sendmax[0]<num_recvmax[0]) because want to 
    //favor sending elements to processors with higher rank
    if(((num_send1[0]>MIN_NUM_2_SEND)||
	(num_send2[0]>0             )
	)&&
       ((num_sendmax[0]<num_recvmax[0])||
	(num_recvmax[0]==0))
       ){
      num_send[0]=num_sendmax[0];
    }
    else if((num_recv1[0]>MIN_NUM_2_SEND)||
	    (num_recv2[0]>0             )){
      num_recv[0]=num_recvmax[0];
    }

    if((num_send[0]>(num_elem-1)/2)||
       ((num_neigh_elem[0]!=-1)&&
	(num_recv[0]>(num_neigh_elem[0]-1)/2))
       )
      If2WaySend[0]=1;
    */
  }
   

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 6.0);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 12.0\n",myid);
    fclose(fpdb2);
  }

#endif

  //figure out how many elements I send-to/receive-from the processor 
  //after me on the space filling curve, the default is zero
  if(myid<numprocs-1) {
    num_recv1[1]=temprecvint[(myid+1)*5+0];
    num_recv2[1]=temprecvint[(myid+1)*5+2];
    num_recvmax[1]=(num_recv1[1]>num_recv2[1])?num_recv1[1]:num_recv2[1];
    num_neigh_elem[1]=temprecvint[(myid+1)*5+4];
    //NumYour2WayElem[1]=(num_neigh_elem[1]-1)/2;


    if(((num_send1[1]<MIN_NUM_2_SEND)&&
	(num_send2[1]==0            )
	)||
       (((num_recv1[1]>=MIN_NUM_2_SEND)||
	 (num_recv2[1]>0))&&
	(num_recv1[1]>=num_recv2[1]))
       ){
      num_send[1]=0;
      if((num_recv1[1]>=MIN_NUM_2_SEND)||
	 (num_recv2[1]>0))
	num_recv[1]=num_recv1[1];
      else
	num_recv[1]=0;
    }
    else if(num_send1[1]>=num_send2[1]){
      num_send[1]=num_send1[1];
      num_recv[1]=0;
    }
    else{
      /*
      num_send[1]=(num_send2[1]<NumMy2WayElem[1])?num_send2[1]:
	NumMy2WayElem[1];
      num_recv[1]=(num_recv2[1]<NumYour2WayElem[1])?num_recv2[1]:
	NumYour2WayElem[1];
      */
      If2WaySend[1]=1;
    }
      
    /*
    //"<=" not "<" in (num_sendmax[1]<=num_recvmax[1]) because want to 
    //favor sending elements to processors with higher rank
    if(((num_send1[1]>MIN_NUM_2_SEND)||
	(num_send2[1]>0             )
	)&&
       ((num_sendmax[1]<=num_recvmax[1])||
	(num_recvmax[1]==0))
       ) {
      num_send[1]=num_sendmax[1]; 
    }
    else if((num_recv1[1]>MIN_NUM_2_SEND)||
	    (num_recv2[1]>0             )){
      num_recv[1]=num_recvmax[1];
    }

    if((num_send[1]>(num_elem-1)/2)||
       ((num_neigh_elem[1]!=-1)&&
	(num_recv[1]>(num_neigh_elem[1]-1)/2))
       )
      If2WaySend[1]=1;
    */

  }

  CDeAllocI1(temprecvint);  //I don't need you anymore

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 7.0);

  /*
  char filename[256];
  sprintf(filename,"sequentialsend%04d.debug",myid);
  FILE *fpdebug=fopen(filename,"w");
  */


  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    FILE *fpdebug=fpdb2;

    fprintf(fpdebug,"myid=%d SequentialSend() A If2WaySend[0]=%d\n  (num_elem      -1)/2=%6d num_send1[0]=%6d num_send2[0]=%6d num_send[0]=%6d\n  (num_neigh_elem-1)/2=%6d num_recv1[0]=%6d num_recv2[0]=%6d num_recv[0]=%6d\n",myid,If2WaySend[0],(num_elem-1)/2,num_send1[0],num_send2[0],num_send[0],(num_neigh_elem[0]-1)/2,num_recv1[0],num_recv2[0],num_recv[0]);
    fprintf(fpdebug,"myid=%d SequentialSend() B If2WaySend[1]=%d\n  (num_elem      -1)/2=%6d num_send1[1]=%6d num_send2[1]=%6d num_send[1]=%6d\n  (num_neigh_elem-1)/2=%6d num_recv1[1]=%6d num_recv2[1]=%6d num_recv[1]=%6d\n",myid,If2WaySend[1],(num_elem-1)/2,num_send1[1],num_send2[1],num_send[1],(num_neigh_elem[1]-1)/2,num_recv1[1],num_recv2[1],num_recv[1]);
    
    fclose(fpdb2);
  }
  
  /*
  fclose(fpdebug);
  fpdebug=fopen(filename,"a");
  */

  MPI_Request Request2WaySend[2];
  MPI_Request Request2WayRecv[2];

  int i2way, send_tag1=070202*6;
  double *DoubleKeyMy2WaySend0, *DoubleKeyYour2WaySend0;
  double *DoubleKeyMy2WaySend1, *DoubleKeyYour2WaySend1;

  if(If2WaySend[0]){
    NumYour2WayElem[0]=(num_neigh_elem[0]-1)/2;
    DoubleKeyMy2WaySend0  =CAllocD1(NumMy2WayElem);
    DoubleKeyYour2WaySend0=CAllocD1(NumYour2WayElem[0]);

    for(ielem=0;ielem<NumMy2WayElem;ielem++)
      DoubleKeyMy2WaySend0[ielem]=DoubleKeyArray[ielem];
    
    ierr=MPI_Irecv((void *) DoubleKeyYour2WaySend0,
		   NumYour2WayElem[0],MPI_DOUBLE,myid-1,send_tag1+myid-1,
		   MPI_COMM_WORLD,&(Request2WayRecv[0]));


    ierr=MPI_Isend((void *) DoubleKeyMy2WaySend0,
		   NumMy2WayElem,MPI_DOUBLE,myid-1,send_tag1+myid,
		   MPI_COMM_WORLD,&(Request2WaySend[0]));
  }

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 8.0);

  if(If2WaySend[1]){
    NumYour2WayElem[1]=(num_neigh_elem[1]-1)/2;
    DoubleKeyMy2WaySend1  =CAllocD1(NumMy2WayElem);
    DoubleKeyYour2WaySend1=CAllocD1(NumYour2WayElem[1]);

    for(ielem=num_elem-NumMy2WayElem,i2way=0;ielem<num_elem;ielem++,i2way++)
      DoubleKeyMy2WaySend1[i2way]=DoubleKeyArray[ielem];
    
    ierr=MPI_Irecv((void *) DoubleKeyYour2WaySend1,
		   NumYour2WayElem[1],MPI_DOUBLE,myid+1,send_tag1+myid+1,
		   MPI_COMM_WORLD,&(Request2WayRecv[1]));


    ierr=MPI_Isend((void *) DoubleKeyMy2WaySend1,
		   NumMy2WayElem,MPI_DOUBLE,myid+1,send_tag1+myid,
		   MPI_COMM_WORLD,&(Request2WaySend[1]));
  }

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 9.0);

  /*
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    FILE* fpdebug=fpdb2;

    fprintf(fpdebug,"myid=%d SequentialSend() C If2WaySend[0]=%d\n  (num_elem      -1)/2=%6d num_send1[0]=%6d num_send2[0]=%6d num_send[0]=%6d\n  (num_neigh_elem-1)/2=%6d num_recv1[0]=%6d num_recv2[0]=%6d num_recv[0]=%6d\n",myid,If2WaySend[0],(num_elem-1)/2,num_send1[0],num_send2[0],num_send[0],(num_neigh_elem[0]-1)/2,num_recv1[0],num_recv2[0],num_recv[0]);
    fprintf(fpdebug,"myid=%d SequentialSend() D If2WaySend[1]=%d\n  (num_elem      -1)/2=%6d num_send1[1]=%6d num_send2[1]=%6d num_send[1]=%6d\n  (num_neigh_elem-1)/2=%6d num_recv1[1]=%6d num_recv2[1]=%6d num_recv[1]=%6d\n",myid,If2WaySend[1],(num_elem-1)/2,num_send1[1],num_send2[1],num_send[1],(num_neigh_elem[1]-1)/2,num_recv1[1],num_recv2[1],num_recv[1]);

    fclose(fpdb2);
  }
  */
  /*
  fclose(fpdebug);
  */

  //If I need to, do a two way send we are only going to correct the key 
  //range overlap and leave load balancing for the next repartitioning
  //The idea is to keep it simple by keeping the same number of elements 
  //on the 2 processors.  Note this is not guaranteed to fix the key range
  //overlap, and the new "maximum key" can end up being less than the new
  //"minimum key."  If that happens for any processor then the sequential 
  //send for all processors repeats until it is corrected.  If at the end 
  //of that my minimum element is less than the maximum element of the 
  //processor before me on the space filling curve or my maximum key is 
  //more than the minimum key on the processor after me on the space 
  //filling curve I will need to do a non-sequential send to fix this.
  //this bit of code counts how many I need to send and receive during the 
  //two way sequential send(s)

  int IfSentRecvd;
  MPI_Status  status;
  int ime, iyou;

  //fpdebug=fopen(filename,"w");
  //fclose(fpdebug);
  //fpdebug=fopen(filename,"a");
  if(If2WaySend[0]||If2WaySend[1])
    //fprintf(fpdebug,"******************\nmyid=%d If2WaySend[0/1]=%d/%d NumYour2WayElem[0/1]=%d/%d NumMy2WayElem=%d\n",myid,If2WaySend[0],If2WaySend[1],NumYour2WayElem[0],NumYour2WayElem[1],NumMy2WayElem);

  while(If2WaySend[0]||If2WaySend[1]) {
    
    if(If2WaySend[0]){
      MPI_Test(Request2WayRecv+0,&IfSentRecvd,&status);      
      if(IfSentRecvd) {
	
	ime=iyou=num_send[0]=num_recv[0]=0;
	for(i2way=0;i2way<NumMy2WayElem+NumYour2WayElem[0];i2way++) {
	  /*
	  if((ime>NumMy2WayElem)||
	     (iyou>NumYour2WayElem[0])||
	     ((ime==NumMy2WayElem)&&(iyou==NumYour2WayElem[0]))
	     ) {
	    fprintf(fpdebug,"myid=%d SequentialSend() E If2WaySend[0]=%d: i2way=%d/%d ime=%d/%d iyou=%d/%d\n",myid,i2way,NumMy2WayElem+NumYour2WayElem[0],ime,NumMy2WayElem,iyou,NumYour2WayElem[0]);
	    fclose(fpdebug);
	    fpdebug=fopen(filename,"a");
	  }
	  */
	  if(ime==NumMy2WayElem) {
	    iyou++;
	    if(i2way>=NumYour2WayElem[0])
	      num_recv[0]++;
	  }
	  else if(iyou==NumYour2WayElem[0]) {
	    ime++;
	    if(i2way<NumYour2WayElem[0])
	      num_send[0]++;
	  }
	  else if(DoubleKeyMy2WaySend0[ime]<DoubleKeyYour2WaySend0[iyou]) {
	    ime++;
	    if(i2way<NumYour2WayElem[0])
	      num_send[0]++;
	  }
	  else{
	    iyou++;
	    if(i2way>=NumYour2WayElem[0])
	      num_recv[0]++;
	  }	    
	}
	//checkelemnode(El_Table, NodeTable, myid, fpdb2, 10.0);
	
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  FILE *fpdebug=fpdb2;

	  fprintf(fpdebug,"myid=%d SequentialSend() EE If2WaySend[0]=%d ime=%d iyou=%d\n  (num_elem      -1)/2=%6d num_send1[0]=%6d num_send2[0]=%6d num_send[0]=%6d firstdk=%20.14g\n  (num_neigh_elem-1)/2=%6d num_recv1[0]=%6d num_recv2[0]=%6d num_recv[0]=%6d lastdk =%20.14g\n",myid,If2WaySend[0],ime,iyou,(num_elem-1)/2,num_send1[0],num_send2[0],num_send[0],DoubleKeyMy2WaySend0[0],(num_neigh_elem[0]-1)/2,num_recv1[0],num_recv2[0],num_recv[0],DoubleKeyYour2WaySend0[(num_neigh_elem[0]-1)/2-1]);

	  fclose(fpdb2);
	}
	  /*	
	    fclose(fpdebug);
	    fpdebug=fopen(filename,"a");
	  */
	//bob

	CDeAllocD1(DoubleKeyYour2WaySend0);
	If2WaySend[0]=0;
      }
    }

    if(If2WaySend[1]){
      MPI_Test(Request2WayRecv+1,&IfSentRecvd,&status);      
      if(IfSentRecvd) {
	
	ime=iyou=num_send[1]=num_recv[1]=0;
	for(i2way=0;i2way<NumMy2WayElem+NumYour2WayElem[1];i2way++) {	  
	  /*
	  if((ime>NumMy2WayElem)||
	     (iyou>NumYour2WayElem[1])||
	     ((ime==NumMy2WayElem)&&(iyou==NumYour2WayElem[1]))
	     ) {
	    fprintf(fpdebug,"myid=%d SequentialSend() F If2WaySend[1]=%d: i2way=%d/%d ime=%d/%d iyou=%d/%d\n",myid,i2way,NumMy2WayElem+NumYour2WayElem[1],ime,NumMy2WayElem,iyou,NumYour2WayElem[1]);
	    fclose(fpdebug);
	    fpdebug=fopen(filename,"a");
	  }
	  */
	  if(ime==NumMy2WayElem) {
	    iyou++;
	    if(i2way<NumMy2WayElem)
	      num_recv[1]++;
	  }
	  else if(iyou==NumYour2WayElem[1]) {
	    ime++;
	    if(i2way>=NumMy2WayElem)
	      num_send[1]++;
	  }
	  else if(DoubleKeyMy2WaySend1[ime]<DoubleKeyYour2WaySend1[iyou]) {
	    ime++;
	    if(i2way>=NumMy2WayElem)
	      num_send[1]++;
	  }
	  else{
	    iyou++;
	    if(i2way<NumMy2WayElem)
	      num_recv[1]++;
	  }
	}

	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  FILE *fpdebug=fpdb2;

	  fprintf(fpdebug,"myid=%d SequentialSend() FF If2WaySend[1]=%d ime=%d iyou=%d\n  (num_elem      -1)/2=%6d num_send1[1]=%6d num_send2[1]=%6d num_send[1]=%6d lastdk =%20.14g\n  (num_neigh_elem-1)/2=%6d num_recv1[1]=%6d num_recv2[1]=%6d num_recv[1]=%6d firstdk=%20.14g\n",myid,If2WaySend[1],ime,iyou,(num_elem-1)/2,num_send1[1],num_send2[1],num_send[1],DoubleKeyMy2WaySend1[(num_elem-1)/2-1],(num_neigh_elem[1]-1)/2,num_recv1[1],num_recv2[1],num_recv[1],DoubleKeyYour2WaySend1[0]);
	  
	  fclose(fpdb2);
	}

	/*
	fclose(fpdebug);
	fpdebug=fopen(filename,"a");

	checkelemnode(El_Table, NodeTable, myid, fpdb2, 11.0);
	*/

	CDeAllocD1(DoubleKeyYour2WaySend1);
	If2WaySend[1]=0;
      }
    }
  }
  //fprintf(fpdebug,"GG\n");
  //fclose(fpdebug);


  if(num_send[0]&&num_recv[0]) If2WaySend[0]=1;
  if(num_send[1]&&num_recv[1]) If2WaySend[1]=1;


  //find the indices of the last/first element I send to the processor
  //before/after me on the space filling curve
  int isend[2]={num_send[0]-1,num_elem-num_send[1]}; //final
  
#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d repart2 at 13.0\n",myid);
    fclose(fpdb2);
  }

#endif

  // **********************************
  // ** section for debugging output **
  // **********************************

#ifdef DEBUG_REPART2  
    //MPI_Barrier(MPI_COMM_WORLD);
    //for(iproc=0;iproc<numprocs-1;iproc++)
    //if(myid>iproc) MPI_Barrier(MPI_COMM_WORLD);

  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 num_elem=%d\n   num_send={%d,%d} num_send1={%d,%d} num_send2={%d,%d}\n   num_recv={%d,%d} num_recv1={%d,%d} num_recv2={%d,%d}\n",myid,num_elem,
	 num_send[ 0],num_send[ 1],
	 num_send1[0],num_send1[1],
	 num_send2[0],num_send2[1],
	 num_recv[ 0],num_recv[ 1],
	 num_recv1[0],num_recv1[1],
	 num_recv2[0],num_recv2[1]);
    fclose(fpdb2);

  }
#endif

  // **************************************************
  // ** Actually send the "first-send" elements away **
  // **************************************************
  MPI_Request requestsent[2]={0,0};
  MPI_Request requestrecv[2]={0,0};
  int send_tag2=061205*4; //YYMMDD date of coding (times 4 for spacing)

  ElemPack *send_array0, *send_array1;
  ElemPack *recv_array0, *recv_array1;

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 12.0);
  if(num_send[0]>0) {
    send_array0=(ElemPack*) malloc(num_send[0]*sizeof(ElemPack));
    //my first key is changing, update my record
    for(ikey=0;ikey<KEYLENGTH;ikey++)
      MyFirstAndLastKey[0*KEYLENGTH+ikey]=
	*(ElemArray[isend[0]+1]->pass_key()+ikey);
    
    for(ielem=0;ielem<=isend[0];ielem++) {
      Pack_element(ElemArray[ielem],send_array0+ielem,NodeTable,myid-1);
      assert(!compare_key(send_array0[ielem].key,nullkey));
    }
    assert(ielem=num_send[0]);

    
#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

      if(num_send[0]==0)
	fprintf(fpdb2,"myid=%d yadaA isend[0]=%d num_send[0]=%d num_elem=%d ielem=%d\n",
	       myid,isend[0],num_send[0],num_elem,ielem);
      else
	fprintf(fpdb2,"myid=%d yadaA isend[0]=%d num_send[0]=%d num_elem=%d ielem=%d {%10u,%10u}\n",
	       myid,isend[0],num_send[0],num_elem,ielem,
	       send_array0[isend[0]].key[0],send_array0[isend[0]].key[1]);

      fclose(fpdb2);
    
    }
#endif

    //send the "first-send" elements to the processor before me
    ierr=MPI_Isend((void*) send_array0,num_send[0],ELEMTYPE,myid-1,send_tag2+myid  ,
		   MPI_COMM_WORLD,requestsent+0);
  }
  
  if(num_recv[0]>0) {
    //receive the "first-send" elements from the processor before me
#ifdef DEBUG_REPART2
    if(timeprops_ptr->iter==DEBUG_ITER){
      fpdb2=fopen(fname2,"a");
      
      fprintf(fpdb2,"myid=%d repart2 at 14.9 num_recv[0]=%d\n",myid,num_recv[0]);
      fclose(fpdb2);
    }
#endif
    recv_array0=(ElemPack*) malloc(num_recv[0]*sizeof(ElemPack));
    ierr=MPI_Irecv((void *) recv_array0,num_recv[0],ELEMTYPE,myid-1,
		   send_tag2+myid-1,MPI_COMM_WORLD,requestrecv+0);
    
    //checkelemnode(El_Table, NodeTable, myid, fpdb2, 14.0);
  }
  
#ifdef DEBUG_REPART2
  //printf("myid=%d repart2 at 15.0\n",myid); fflush(stdout);
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2  isend[1]=%d num_send[1]=%d num_elem=%d\n",
	    myid,isend[1],num_send[1],num_elem);
    fclose(fpdb2);
  }
  
#endif
  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 15.0);

  if(num_send[1]>0) {
    send_array1=(ElemPack*) malloc(num_send[1]*sizeof(ElemPack));
    //my last key is changing, update my record
    for(ikey=0;ikey<KEYLENGTH;ikey++)
      MyFirstAndLastKey[1*KEYLENGTH+ikey]=
	*(ElemArray[isend[1]-1]->pass_key()+ikey);
    
    for(ielem=isend[1];ielem<num_elem;ielem++) {
      Pack_element(ElemArray[ielem],send_array1+ielem-isend[1],NodeTable,myid+1);
      assert(!compare_key(send_array1[ielem-isend[1]].key,nullkey));
    }
    assert(ielem-isend[1]==num_send[1]);
    
#ifdef DEBUG_REPART2    
    if(timeprops_ptr->iter==DEBUG_ITER){
      fpdb2=fopen(fname2,"a");
      
      if(num_send[1]==0)
	fprintf(fpdb2,"myid=%d yadaB isend[1]=%d num_send[1]=%d num_elem=%d ielem=%d\n",
		myid,isend[1],num_send[1],num_elem,ielem);
      else
	fprintf(fpdb2,"myid=%d yadaB isend[1]=%d num_send[1]=%d num_elem=%d ielem=%d {%10u,%10u}\n",
		myid,isend[1],num_send[1],num_elem,ielem,
		send_array1[0].key[0],send_array1[0].key[1]);
      
      fclose(fpdb2);
      
    }
#endif    

    //checkelemnode(El_Table, NodeTable, myid, fpdb2, 16.0);
    
    //send the "first-send" elements to the processor after me
    ierr=MPI_Isend((void*) send_array1,num_send[1],ELEMTYPE,myid+1,send_tag2+myid  ,
		   MPI_COMM_WORLD,requestsent+1);  
  }

  if(num_recv[1]) {
    //receive the "first-send" elements from the processor after me
#ifdef DEBUG_REPART2
    if(timeprops_ptr->iter==DEBUG_ITER){
      fpdb2=fopen(fname2,"a");
      
      fprintf(fpdb2,"myid=%d repart2 at 15.9 num_recv[1]=%d\n",myid,num_recv[1]);
      fclose(fpdb2);
    }
#endif
    recv_array1=(ElemPack*) malloc(num_recv[1]*sizeof(ElemPack));
    ierr=MPI_Irecv((void *) recv_array1,num_recv[1],ELEMTYPE,myid+1,
		   send_tag2+myid+1,MPI_COMM_WORLD,requestrecv+1);
    
    //checkelemnode(El_Table, NodeTable, myid, fpdb2, 17.0);
  }
  
#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 at 16.0\n",myid);
    fclose(fpdb2);
  }
#endif
  
  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 18.0);
  
  // **********************************
  // ** section for debugging output **
  // **********************************

#ifdef DEBUG_REPART2
  //MPI_Barrier(MPI_COMM_WORLD);
  //for(iproc=0;iproc<numprocs-1;iproc++)
  //if(myid>iproc) MPI_Barrier(MPI_COMM_WORLD);

  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d B) MyFirstAndLastKey={{%10u,%10u},{%10u,%10u}}\n",myid,
	    MyFirstAndLastKey[0],MyFirstAndLastKey[1],
	    MyFirstAndLastKey[2],MyFirstAndLastKey[3]); fflush(stdout);

    //for(iproc=1;iproc<numprocs;iproc++)
    //if(myid<iproc) MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);

    fprintf(fpdb2,"myid=%d repart2 at 16.1\n",myid);
    fclose(fpdb2);
  }
#endif
  

  // **************************************************
  // ** do stuff to pass the time until "first-send" **
  // ** elements arrive                              **
  // **************************************************

  //Delete the elements that I just sent away, recall all 
  //other non-active elements were deleted at the beginning of
  //repartition2()
  
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"repartition2(): myid=%d\nabout to delete elements %d through %d out of %d because I sent them to processor %d\n",myid,0,isend[0],num_elem,myid-1);
    for(ielem=0;ielem<=isend[0];ielem++)
      fprintf(fpdb2,"    deleting %dth element: key={%u,%u}\n",ielem,
	      *(ElemArray[ielem]->pass_key()+0),
	      *(ElemArray[ielem]->pass_key()+1));
    
    fprintf(fpdb2,"\nabout to delete elements %d through %d out of %d because I sent the to processor %d\n",isend[1],num_elem-1,num_elem,myid+1);
    for(ielem=isend[1];ielem<num_elem;ielem++)
      fprintf(fpdb2,"    deleting %dth element: key={%u,%u}\n",ielem,
	     *(ElemArray[ielem]->pass_key()+0),
	     *(ElemArray[ielem]->pass_key()+1));

    fclose(fpdb2);
  }

    /*

    fflush(stdout);
  }
  */

  if(num_send[0]>0)
    for(ielem=0;ielem<=isend[0];ielem++) {
      ElemArray[ielem]->void_bcptr();
      El_Table->remove(ElemArray[ielem]->pass_key(),1,stdout,myid,2);
      delete ElemArray[ielem];
    }
  
  if(num_send[1]>0){
    for(ielem=isend[1];ielem<num_elem;ielem++) {
      El_Table->remove(ElemArray[ielem]->pass_key(),1,stdout,myid,3);
      delete ElemArray[ielem];
    }
  }
  free(ElemArray);
  free(DoubleKeyArray);
  free(LoadBalArray);

  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 19.0);

  //All the non active elements have already been deleted. Now
  //delete unused nodes now to give MPI more time to send, and 
  //we won't have to do this in H_adapt() now which also saves 
  //computation
  int NodeTable_num_buck=NodeTable->get_no_of_buckets();
  HashEntryPtr *NodeTable_bucket0=NodeTable->getbucketptr();
  HashEntryPtr NodeTable_entry_ptr;
  int inodebucket, inode;
  Node *NdTemp;

  //zero the number of elems each node is associated with
  for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
    NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
        
    while(NodeTable_entry_ptr) {
	
      NdTemp=(Node*)(NodeTable_entry_ptr->value);
      NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
      assert(NdTemp);
      NdTemp->put_num_assoc_elem(0);
    }//while(NodeTable_entry_ptr)
  }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++

  
  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 20.0);

  //loop through elements to see which nodes are associated with them
  for(ibuck=0; ibuck<num_buck; ibuck++) {
    currentPtr = *(buck+ibuck);
    
    while(currentPtr){	
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;     
      assert(EmTemp);
      NdTemp=(Node*) NodeTable->lookup(EmTemp->pass_key());
      assert(NdTemp);
      NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);

      for(inode=0;inode<8;inode++){
	NdTemp=(Node*) NodeTable->lookup(EmTemp->getNode()+inode*KEYLENGTH);
	if(NdTemp==NULL){
	  //char fname2[256];
	  //sprintf(fname2,"seqsendmissnode%04d.debug",myid);
	  fpdb2=fopen(fname2,"a");
	  //checkelemnode(El_Table, NodeTable, myid, fpdb2, 21.0);
	  //fpdb2=stdout;
	  fprintf(fpdb2,"myid=%d iter=%d iseqsend=%d node inode=%d is missing\n",myid,timeprops_ptr->iter,iseqsend,inode);
	  ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),fpdb2);
	  fclose(fpdb2);
	  fpdb2=fopen(fname2,"a");
	  NodeBackgroundCheck(El_Table,NodeTable,EmTemp->getNode()+inode*KEYLENGTH,fpdb2);
	  //unsigned tempkey[2]={695892755,2973438897};
	  //printf(

	  fclose(fpdb2);

	}

	assert(NdTemp);
	NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);
      }
    }//while(currentPtr)
  }//for(ibuck=0; ibuck<num_buck; ibuck++)


  //if a node is not associated with any elements delete it
  for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
    NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
        
    while(NodeTable_entry_ptr) {
	
      NdTemp=(Node*)(NodeTable_entry_ptr->value);
      NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
      assert(NdTemp);

      if(NdTemp->get_num_assoc_elem()==0){
	NodeTable->remove(NdTemp->pass_key(), 0,stdout,myid,4);
	delete NdTemp;
      }
    }//while(NodeTable_entry_ptr)
  }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"myid=%d repart2 at 17.0\n",myid);
    fclose(fpdb2);
  }
#endif



  // *******************************************
  // ** NOW ACCEPT THE FIRST-RECEIVE ELEMENTS **
  // *******************************************

  do{
    if(num_recv[0]>0) {
      MPI_Test(requestrecv+0,&IfSentRecvd,&status);      
      if(IfSentRecvd) {
#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.1... received \"sequential send\" elements from processor %d\n",myid,myid-1);
	  fclose(fpdb2);
	}
#endif
	IncorporateNewElements(El_Table,NodeTable,myid,num_recv[0],
			       recv_array0,timeprops_ptr);

#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");


	  fprintf(fpdb2,"myid=%d sender id=%d num_recv[0]=%d IncomingFirstLastKey={{%10u,%10u} {%10u,%10u}}\n",myid,myid-1,num_recv[0],
		  recv_array0[            0].key[0],
		  recv_array0[            0].key[1],
		  recv_array0[num_recv[0]-1].key[0],
		  recv_array0[num_recv[0]-1].key[1]);
	  fclose(fpdb2);
	}
#endif

	doublekey=recv_array0[0].key[0];
	for(ikey=1;ikey<KEYLENGTH;ikey++)
	  doublekey=doublekey*doublekeyrange1+
	    recv_array0[0].key[ikey];

	if(doublekey<AllFirstAndLastDoubleKeys[myid][0])
	  //my first key has changed, update my record
	  for(ikey=0;ikey<KEYLENGTH;ikey++)
	    MyFirstAndLastKey[0*KEYLENGTH+ikey]=
	      recv_array0[0].key[ikey];

	free(recv_array0);
	num_recv[0]=0;
      }
    }

    if(num_recv[1]>0) {
      MPI_Test(requestrecv+1,&IfSentRecvd,&status);
      if(IfSentRecvd) {
#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.3... received \"sequential send\" elements from processor %d\n",myid,myid+1);
	  fclose(fpdb2);
	}
#endif
	IncorporateNewElements(El_Table,NodeTable,myid,num_recv[1],
			       recv_array1,timeprops_ptr);
#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");


	  fprintf(fpdb2,"myid=%d sender id=%d num_recv[1]=%d IncomingFirstLastKey={{%10u,%10u} {%10u,%10u}}\n",myid,myid+1,num_recv[1],
		  recv_array1[            0].key[0],
		  recv_array1[            0].key[1],
		  recv_array1[num_recv[1]-1].key[0],
		  recv_array1[num_recv[1]-1].key[1]);
	  fclose(fpdb2);
	}
#endif
	doublekey=recv_array1[num_recv[1]-1].key[0];
	for(ikey=1;ikey<KEYLENGTH;ikey++)
	  doublekey=doublekey*doublekeyrange1+
	    recv_array1[num_recv[1]-1].key[ikey];

	if(AllFirstAndLastDoubleKeys[myid][1]<doublekey)
	  //my last key has changed, update my record
	  for(ikey=0;ikey<KEYLENGTH;ikey++)
	    MyFirstAndLastKey[1*KEYLENGTH+ikey]=
	      recv_array1[num_recv[1]-1].key[ikey];

	free(recv_array1);
	num_recv[1]=0;
      }
    }

    if(num_send[0]>0) {
      MPI_Test(requestsent+0,&IfSentRecvd,&status);
      if(IfSentRecvd){
#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.5 about to free element sent to proc %d\n",myid,myid-1);
	  fclose(fpdb2);
	}
#endif

	free(send_array0);
	num_send[0]=0;

#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.6 I have freed element sent to proc %d\n",myid,myid-1);
	  fclose(fpdb2);
	}
#endif


      }
    }

    if(num_send[1]>0) {
      MPI_Test(requestsent+1,&IfSentRecvd,&status);
      if(IfSentRecvd){
#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.7 about to free element sent to proc %d\n",myid,myid+1);
	  fclose(fpdb2);
	}
#endif

	free(send_array1);
	num_send[1]=0;

#ifdef DEBUG_REPART2
	if(timeprops_ptr->iter==DEBUG_ITER){
	  fpdb2=fopen(fname2,"a");
	  
	  fprintf(fpdb2,"myid=%d repart2 at 17.8 I have freed element sent to proc %d\n",myid,myid+1);
	  fclose(fpdb2);
	}
#endif
      }
    }

  }while(num_send[0]+num_send[1]+num_recv[0]+num_recv[1]);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    
    fprintf(fpdb2,"myid=%d repart2 at 18.0\n",myid);
    fclose(fpdb2);
  }
#endif

  //CDeAllocD2(AllFirstAndLastDoubleKeys);

  //these are guaranteed to have arrived at their destination by this 
  //point so we don't need to do an MPI wait before DeAllocating them
  if(If2WaySend[0])
    CDeAllocD1(DoubleKeyMy2WaySend0);

  if(If2WaySend[1])  
    CDeAllocD1(DoubleKeyMy2WaySend1);

  MyFirstAndLastDoubleKey[0]=MyFirstAndLastKey[0];
  MyFirstAndLastDoubleKey[1]=MyFirstAndLastKey[KEYLENGTH];
  for(ikey=1;ikey<KEYLENGTH;ikey++) {
    MyFirstAndLastDoubleKey[0]=
      MyFirstAndLastDoubleKey[0]*doublekeyrange1+
      MyFirstAndLastKey[ikey];
    MyFirstAndLastDoubleKey[1]=
      MyFirstAndLastDoubleKey[1]*doublekeyrange1+
      MyFirstAndLastKey[KEYLENGTH+ikey];        
  }


  //determine the interpartition boundaries on the spacefilling curve
  if(myid==0)
    MyFirstAndLastDoubleKey[0]=0.0;
  
  if(myid==numprocs-1) {
    MyFirstAndLastDoubleKey[1]=doublekeyrange1;
    for(ikey=1;ikey<KEYLENGTH;ikey++) 
      MyFirstAndLastDoubleKey[1]=
	MyFirstAndLastDoubleKey[1]*doublekeyrange1+
	doublekeyrange1;
  }
  
  //double **AllFirstAndLastDoubleKeys=CAllocD2(numprocs,2);
  
  MPI_Allgather(  MyFirstAndLastDoubleKey ,2,MPI_DOUBLE,
		  *AllFirstAndLastDoubleKeys,2,MPI_DOUBLE,MPI_COMM_WORLD); 
  
  
  NewProcDoubleKeyBoundaries[0]=-1.0;

  int ifrepeat=0;

  int jproc;
  for(iproc=0;iproc<numprocs;iproc++) {
    if(AllFirstAndLastDoubleKeys[iproc][0]>
       AllFirstAndLastDoubleKeys[iproc][1]) {
      ifrepeat=1;
      break;
    }
    for(jproc=iproc+1;jproc<numprocs;jproc++)
      if((AllFirstAndLastDoubleKeys[iproc][0]>
	  AllFirstAndLastDoubleKeys[jproc][0])||
	 (AllFirstAndLastDoubleKeys[iproc][1]>
	  AllFirstAndLastDoubleKeys[jproc][1])) {
	ifrepeat=1;
	break;
      }
    if(ifrepeat) break;
  }

 
  if(!ifrepeat) 
    for(iproc=0;iproc<numprocs;iproc++)     
      NewProcDoubleKeyBoundaries[iproc+1]=
	AllFirstAndLastDoubleKeys[iproc][1];
  

  CDeAllocD2(AllFirstAndLastDoubleKeys);

#ifdef DEBUG_REPART2
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");

    fprintf(fpdb2,"************************\n************************\n************************\ndone\n************************\n************************\n************************\n");
    fclose(fpdb2);
  }
#endif

  return ifrepeat;
}









//bob

void NonSequentialSendAndUpdateNeigh(int numprocs, int myid, 
				     HashTable* El_Table, HashTable* NodeTable,
				     TimeProps* timeprops_ptr, 
				     double *NewProcDoubleKeyBoundaries) {


  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  int ibuck, ielem, ikey, iproc, ineigh, ierr, inodebucket, inode;
  Element *EmTemp;
  unsigned nullkey[2]={0,0};
  double doublekeyrange1=*(El_Table->get_doublekeyrange()+1);
  double doublekey;
  int NodeTable_num_buck=NodeTable->get_no_of_buckets();
  HashEntryPtr *NodeTable_bucket0=NodeTable->getbucketptr();
  HashEntryPtr NodeTable_entry_ptr;
  Node *NdTemp;
  int IfSentRecvd;
  MPI_Status status;


  char fname2[256];
  sprintf(fname2,"nonseqsend%04d.debug",myid);
  FILE *fpdb2; 

#ifdef DEBUG_REPART2C
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"w");
    fprintf(fpdb2,"myid=%d iter=%d\n",myid,timeprops_ptr->iter);
    for(iproc=0;iproc<=numprocs;iproc++)
      fprintf(fpdb2,"   NewProcDoubleKeyBoundaries[%3d]=%20.14g\n",
	      iproc,NewProcDoubleKeyBoundaries[iproc]);
    

    fclose(fpdb2);
  }
#endif
  
  // ***********************************************************************
  // ***********************************************************************
  // ** In case there was non sequential processor overlap we now have to **
  // ** so a second-send (only a second, and not a third or more, send    **
  // ** because the partion boundaries were finalized by the first send)  **
  // ***********************************************************************
  // ***********************************************************************




  double MyFirstAndLastDoubleKey[2];

  MyFirstAndLastDoubleKey[0]=
    NewProcDoubleKeyBoundaries[myid];
  MyFirstAndLastDoubleKey[1]=
    NewProcDoubleKeyBoundaries[myid+1];

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 19.0\n",myid); fflush(stdout);
#endif 

  // **************************************************************
  // now that I know all the boundaries of all processors I am able
  // to determine 
  // 1) which elements I "own" but should not
  // 2) which processors these elements should actually belong to
  // **************************************************************
  ElemPtrList NotMyElem(256);  //will deallocate self with destructor
  int *NumToSecondSend=CAllocI1(numprocs);
  int *NumToSecondRecv=CAllocI1(numprocs);

  for(iproc=0;iproc<numprocs;iproc++)
    NumToSecondSend[iproc]=NumToSecondRecv[iproc]=0;

  //figure out which elements I shouldn't own and how many elements I
  //should send to each of the other processors.
  for(ibuck=0; ibuck<num_buck; ibuck++) {
    currentPtr = *(buck+ibuck);

    while(currentPtr){
	
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;     
      assert(EmTemp);

#ifdef DEBUG_REPART2C
      if((timeprops_ptr->iter==DEBUG_ITER)&&
	 (EmTemp->get_adapted_flag()<=-NOTRECADAPTED)&& 
	 (EmTemp->get_adapted_flag()>=-BUFFER)) {	

	doublekey=*(EmTemp->pass_key()+0);
	for(ikey=1;ikey<KEYLENGTH;ikey++) 
	  doublekey=doublekey*doublekeyrange1+
	    *(EmTemp->pass_key()+ikey);

	fpdb2=fopen(fname2,"a");
	fprintf(fpdb2,"myid=%d iter=%d MyFLDK={%20.14g,%20.14g} NPDKB={%20.14g,%20.14g} elemdoublekey=%20.14g\n",
		myid,timeprops_ptr->iter,
		MyFirstAndLastDoubleKey[0],
		MyFirstAndLastDoubleKey[1],
		NewProcDoubleKeyBoundaries[myid],
		NewProcDoubleKeyBoundaries[myid+1],
		doublekey);
	fclose(fpdb2);
      }
#endif     

      if(EmTemp->get_adapted_flag()>=NOTRECADAPTED) {
	doublekey=*(EmTemp->pass_key()+0);
	for(ikey=1;ikey<KEYLENGTH;ikey++) 
	  doublekey=doublekey*doublekeyrange1+
	    *(EmTemp->pass_key()+ikey);
	
	if((doublekey<MyFirstAndLastDoubleKey[0])||
	   (doublekey>MyFirstAndLastDoubleKey[1])
	   ) {
	  NotMyElem.add(EmTemp);
	  for(iproc=0;iproc<numprocs;iproc++)
	    if((NewProcDoubleKeyBoundaries[iproc  ]< doublekey)&&
	       (NewProcDoubleKeyBoundaries[iproc+1]>=doublekey)
	       ) {
	      NumToSecondSend[iproc]++;
	      break;
	    }
	  /*
	  printf("myid=%d iter=%d: element {%10u,%10u} does not belong to me but says it does\n",
		 myid,timeprops_ptr->iter,
		 *(EmTemp->pass_key()+0),*(EmTemp->pass_key()+1));
	  ElemBackgroundCheck(El_Table,NodeTable,EmTemp->pass_key(),stdout);
	  assert(0);
	  */
	}
      }
    }
  }

  assert(NumToSecondSend[myid]==0);

  //tell all the other processsors how many elements I will second send to them
  //also find out how many elements I will second receive from each of them
  MPI_Alltoall(NumToSecondSend,1, MPI_INT, NumToSecondRecv, 1, MPI_INT, MPI_COMM_WORLD);
  
#ifdef DEBUG_REPART2C
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d iter=%d\n",myid,timeprops_ptr->iter);
    for(iproc=0;iproc<numprocs;iproc++)
      fprintf(fpdb2,"     to/from iproc=%4d  #send=%5d #recv=%5d\n",
	      iproc,NumToSecondSend[iproc],NumToSecondRecv[iproc]);
    fclose(fpdb2);
  }
#endif


  assert(NumToSecondRecv[myid]==0);

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 20.0\n",myid); fflush(stdout);
#endif

  //now that I know how many elements I will second-send/second-receive from each
  //of the other processors I can allocate space for some send/receive info/data
  ElemPack **SecondSendArray=(ElemPack**) malloc(numprocs*sizeof(ElemPack*));
  ElemPack **SecondRecvArray=(ElemPack**) malloc(numprocs*sizeof(ElemPack*));
  int *IfSecondSendDone=CAllocI1(numprocs); //this will initially hold the 
  //current element place number in each of the second-send arrays and then
  //will "morph" into whether or not I'm done sending to them.
  int *IfSecondRecvDone=CAllocI1(numprocs); //this will hold whether or not
  //(1 or 0) I'm done second-receiving from a given processor
  MPI_Request *RequestSecondSend=(MPI_Request*) malloc(numprocs*sizeof(MPI_Request));
  MPI_Request *RequestSecondRecv=(MPI_Request*) malloc(numprocs*sizeof(MPI_Request));
  
  int send_tag3=070102*5; //YYMMDD date of coding (times 5 for spacing)

  //allocate space to second-send/second-receive and initiate the non blocking second 
  //receive
  for(iproc=0;iproc<numprocs;iproc++) {
    if(NumToSecondRecv[iproc]>0) {
      IfSecondRecvDone[iproc]=0;
      SecondRecvArray[iproc]=(ElemPack*) malloc(NumToSecondRecv[iproc]*sizeof(ElemPack));
      ierr=MPI_Irecv((void *) SecondRecvArray[iproc],NumToSecondRecv[iproc],ELEMTYPE,iproc,
		     send_tag3+iproc,MPI_COMM_WORLD,&(RequestSecondRecv[iproc]));
    }
    else //if I'm not going to second-receive from them I'm already done 
      //second-receiving from them
      IfSecondRecvDone[iproc]=1;

    IfSecondSendDone[iproc]=0; //initialize current element place in 
    //second-send arrays to zero
    if(NumToSecondSend[iproc]>0)
      SecondSendArray[iproc]=(ElemPack*) malloc(NumToSecondSend[iproc]*sizeof(ElemPack));    
  }


  //go through the short list of Elements that shouldn't belong to me and pack them up
  //to send them to the elements they should belong to
  for(ielem=0;ielem<NotMyElem.get_num_elem();ielem++) {
    doublekey=*(NotMyElem.get_key(ielem));
    for(ikey=1;ikey<KEYLENGTH;ikey++) 
      doublekey=doublekey*doublekeyrange1+
	*(NotMyElem.get_key(ielem)+ikey);

    for(iproc=0;iproc<numprocs;iproc++)
      if((NewProcDoubleKeyBoundaries[iproc  ]< doublekey)&&
	 (NewProcDoubleKeyBoundaries[iproc+1]>=doublekey)
	 ) {
	assert(NumToSecondSend[iproc]>IfSecondSendDone[iproc]);
	Pack_element(NotMyElem.get(ielem),
		     SecondSendArray[iproc]+IfSecondSendDone[iproc],
		     NodeTable,iproc);

	IfSecondSendDone[iproc]++;
	break;
      }
  }

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 21.0\n",myid); fflush(stdout);
#endif
  
  int iSecondSend;
  for(iproc=0;iproc<numprocs;iproc++) {
#ifdef DEBUG_REPART2C
    if(timeprops_ptr->iter==DEBUG_ITER){
      fpdb2=fopen(fname2,"a");
      fprintf(fpdb2,"myid=%d iter=%d iproc=%d IfSecondSendDone=%d NumToSecondSend=%d NumNotMyElem=%d\n",myid,timeprops_ptr->iter,iproc,IfSecondSendDone[iproc],NumToSecondSend[iproc],NotMyElem.get_num_elem());
      fclose(fpdb2);
    }
#endif

    assert(IfSecondSendDone[iproc]==NumToSecondSend[iproc]);

    IfSecondSendDone[iproc]=!IfSecondSendDone[iproc]; //It's morphing time!
    //IfSecondSendDone now says whether or not I'm done second-sending to 
    //each processor, If I'm not going to second-send to a given processor
    //I'm already done second-sending to them
    
    if(NumToSecondSend[iproc]>0) {
      //I need to second-send to processor iproc

      //make sure I'm not sending any empty elements
      for(iSecondSend=0;iSecondSend<NumToSecondSend[iproc];iSecondSend++)
	assert(!compare_key(SecondSendArray[iproc][iSecondSend].key,nullkey));
      
      //non blocking second-send to processor iproc
      ierr=MPI_Isend((void *) SecondSendArray[iproc],NumToSecondSend[iproc],ELEMTYPE,iproc,
		     send_tag3+myid,MPI_COMM_WORLD,&(RequestSecondSend[iproc]));
#ifdef DEBUG_REPART2C
      if(timeprops_ptr->iter==DEBUG_ITER){
	fpdb2=fopen(fname2,"a");
	fprintf(fpdb2,"myid=%d iter=%d did a non blocking nonsequential send to proc %d  #elem sent=%d\n",myid,timeprops_ptr->iter,iproc,NumToSecondSend[iproc]);
	fclose(fpdb2);
      }
#endif
      
    }
  }

#ifdef DEBUG_REPART2C
  if(timeprops_ptr->iter==DEBUG_ITER){
    fpdb2=fopen(fname2,"a");
    fprintf(fpdb2,"myid=%d iter=%d all non blocking nonsequential sends have been sent\n",myid,timeprops_ptr->iter);
    fclose(fpdb2);
      }
#endif



  
  // ******************************************************************
  // pass a little bit of time while waiting to receive the second-send 
  // elements
  // ******************************************************************

  // remove second-send elements from my hashtable
  for(ielem=0;ielem<NotMyElem.get_num_elem();ielem++) {
    NotMyElem.get(ielem)->void_bcptr();
    El_Table->remove(NotMyElem.get_key(ielem),1,stdout,myid,5);
    delete NotMyElem.get(ielem);
  }

  
  if(NotMyElem.get_num_elem()>0) {
    //delete the extra nodes associated with the second-send elements
    
    //zero the number of elems each node is associated with
    for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
      NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
      
      while(NodeTable_entry_ptr) {
	
	NdTemp=(Node*)(NodeTable_entry_ptr->value);
	NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
	assert(NdTemp);
	NdTemp->put_num_assoc_elem(0);
      }//while(NodeTable_entry_ptr)
    }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++
    
  
    //loop through elements to see which nodes are associated with them
    for(ibuck=0; ibuck<num_buck; ibuck++) {
      currentPtr = *(buck+ibuck);
      
      while(currentPtr){	
	EmTemp=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;     
	assert(EmTemp);
	NdTemp=(Node*) NodeTable->lookup(EmTemp->pass_key());
	assert(NdTemp);
	NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);
      
	for(inode=0;inode<8;inode++){
	  NdTemp=(Node*) NodeTable->lookup(EmTemp->getNode()+inode*KEYLENGTH);
	  assert(NdTemp);
	  NdTemp->put_num_assoc_elem(NdTemp->get_num_assoc_elem()+1);
	}
      }//while(currentPtr)
    }//for(ibuck=0; ibuck<num_buck; ibuck++)
    

    //if a node is not associated with any elements delete it
    for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++) {
      NodeTable_entry_ptr = *(NodeTable_bucket0+inodebucket);
      
      while(NodeTable_entry_ptr) {
	
	NdTemp=(Node*)(NodeTable_entry_ptr->value);
	NodeTable_entry_ptr=NodeTable_entry_ptr->next;     
	assert(NdTemp);
	
	if(NdTemp->get_num_assoc_elem()==0){
	  NodeTable->remove(NdTemp->pass_key(), 0,stdout,myid,6);
	  delete NdTemp;
	}
      }//while(NodeTable_entry_ptr)
    }//for(inodebucket=0; inodebucket<NodeTable_num_buck; inodebucket++
  }//if(NotMyElem.get_num_elem()>0)



  // ***************************************************************
  // I've passed as much time as I can, now receive the second-send 
  // elements and incorporate them
  // ***************************************************************

  int NumNotRecvd;
  do{

    NumNotRecvd=0;
    for(iproc=0;iproc<numprocs;iproc++) 
      if(!IfSecondRecvDone[iproc]) {
	MPI_Test(&(RequestSecondRecv[iproc]),&IfSentRecvd,&status);   

	if(IfSentRecvd) {
	  IncorporateNewElements(El_Table,NodeTable,myid,NumToSecondRecv[iproc],
				 SecondRecvArray[iproc],timeprops_ptr);     
	  IfSecondRecvDone[iproc]=1;
	  free(SecondRecvArray[iproc]);
	}
	else
	  NumNotRecvd++;
      }
  }while(NumNotRecvd>0);
  free(SecondRecvArray);
  free(RequestSecondRecv);
  CDeAllocI1(NumToSecondRecv);
  CDeAllocI1(IfSecondRecvDone);



#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 22.0\n",myid); fflush(stdout);
#endif

  // ****************************************************************
  // When debugging it is useful to compare the actual unsigned keys
  // to their double precision versions this section prints that info
  // out 
  // ****************************************************************

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 before AllFirstAndLastKeys\n",myid); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  unsigned *AllFirstAndLastKeys=CAllocU1(2*KEYLENGTH*numprocs);

  MPI_Allgather( MyFirstAndLastKey ,2*KEYLENGTH,MPI_UNSIGNED,
                AllFirstAndLastKeys,2*KEYLENGTH,MPI_UNSIGNED,
		MPI_COMM_WORLD);

  FILE *fp;
  char filename[256];
  if(0&&timeprops_ptr->iter==99) 
  {
    //printf("myid=%d repart2 numprocs=%d\n",myid,numprocs);

    sprintf(filename,"firstandlastkeys%04d.%06d",myid,timeprops_ptr->iter);
    fp=fopen(filename,"w");
  
  
    for(iproc=0;iproc<numprocs;iproc++)
      fprintf(fp,"iproc=%d {%10u,%10u} {%10u,%10u}\n",iproc,
	      AllFirstAndLastKeys[iproc*2*KEYLENGTH+0],
	      AllFirstAndLastKeys[iproc*2*KEYLENGTH+1],
	      AllFirstAndLastKeys[iproc*2*KEYLENGTH+2],
	      AllFirstAndLastKeys[iproc*2*KEYLENGTH+3]);

    fprintf(fp,"doublekeyrange1=%20.14g\nNPDKB={%20.14g,%20.14g,%20.14g",
	    doublekeyrange1,	 
	    NewProcDoubleKeyBoundaries[0],
	    NewProcDoubleKeyBoundaries[1],
	    NewProcDoubleKeyBoundaries[2]);
    
    for(iproc=2;iproc<numprocs;iproc++)
      fprintf(fp,",%20.14g",NewProcDoubleKeyBoundaries[iproc+1]);
    
    fprintf(fp,"}\n");
    
    fclose(fp);

    CDeAllocU1(AllFirstAndLastKeys);
  }

  printf("myid=%d repart2 at 23.0\n",myid); fflush(stdout);
#endif

  // ****************************************************************
  // While I am waiting for my second-send to complete, I will update
  // the neigh_proc of EVERY element in my hashtable. EVERY because 
  // the processes of elements neighboring the elements I now own 
  // could have change.  I can tell which processor every element 
  // will belong to after the second-send/second-receive completes 
  // just by comparing their keys to NewProcDoubleKeyBoundaries
  // ****************************************************************

  //Now scan the hashtable
  int num_elem=0;
  for(ibuck=0; ibuck<num_buck; ibuck++) {
    currentPtr = *(buck+ibuck);

    while(currentPtr){
	
      EmTemp=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;     
      assert(EmTemp);
      assert(EmTemp->get_adapted_flag()>=NOTRECADAPTED);
      num_elem++;
      EmTemp->put_myprocess(myid);
      
      for(ineigh=0;ineigh<8;ineigh++)
	if(*(EmTemp->get_neigh_proc()+ineigh)>=0) {
	  
	  //make a double precision version of the neighbor's key
	  doublekey=*(EmTemp->get_neighbors()+ineigh*KEYLENGTH);
	  for(ikey=1;ikey<KEYLENGTH;ikey++)
	    doublekey=doublekey*doublekeyrange1+
	      *(EmTemp->get_neighbors()+ineigh*KEYLENGTH+ikey);

	  //check which processor the neighbor's key belongs to
	  for(iproc=0;iproc<numprocs;iproc++)
	    if((NewProcDoubleKeyBoundaries[iproc  ]< doublekey)&&
	       (NewProcDoubleKeyBoundaries[iproc+1]>=doublekey)
	       ) {
	      EmTemp->put_neigh_proc(ineigh,iproc);
	      break;
	    }

	}//if(EmTemp->get_neigh_proc(ineigh)>=0) 
    }//while(currentPtr)	
  }//for(ibuck=0; ibuck<num_buck; ibuck++)


  //congratulations we've updated every element's information 
  //about which processors its neighbors belong too.

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 24.0\n",myid); fflush(stdout);
#endif

  /***************************************************************/
  /** I've done everything possible to kill time while I waited **/
  /** for second-send to complete.  Now I will wait some more   **/
  /** and deallocate space as soon as I'm allowed to.           **/
  /***************************************************************/
  int NumNotSent;
  do{

    NumNotSent=0;
    for(iproc=0;iproc<numprocs;iproc++) 
      if(!IfSecondSendDone[iproc]) {
	MPI_Test(&(RequestSecondSend[iproc]),&IfSentRecvd,&status);   

	if(IfSentRecvd) {
	  IfSecondSendDone[iproc]=1;
	  free(SecondSendArray[iproc]);
	}
	else
	  NumNotSent++;
      }
  }while(NumNotSent>0);
  free(SecondSendArray);
  free(RequestSecondSend);
  CDeAllocI1(NumToSecondSend);
  CDeAllocI1(IfSecondSendDone);

#ifdef DEBUG_REPART2B
  printf("myid=%d repart2 at 25.0\n",myid); fflush(stdout);
#endif

  return;
}



void IncorporateNewElements(HashTable* El_Table, HashTable* NodeTable,
			    int myid, int num_recv, ElemPack *recv_array,
			    TimeProps* timeprops_ptr) {
  int ielem;
  Element *EmTemp;

  unsigned nullkey[2]={0,0};
  
  for(ielem=0;ielem<num_recv;ielem++) {
    if(compare_key(recv_array[ielem].key,nullkey)){
      printf("myid=%d num_recv=%d recv_array[%d].key=={0,0}\n",
	     myid,num_recv,ielem);
      assert(0); 
    }
    EmTemp = (Element*) (El_Table->lookup(recv_array[ielem].key));

    assert(EmTemp==NULL);  //this forces a deletion of ghost elements 
    //before repartitioning, which is done within repartition2(), if
    //deletion of ghost cells is skipped it causes drastic problems.
    //for _starters_ ghost elements might be missing nodes and new 
    //nodes get created when new elements are created, so if we 
    //create new elements we can be sure that they will have all 
    //their node, secondly the presence of ghost cells screws up the
    //neighbor updating.
    
    Element* EmNew = new Element();
    double not_used;
    //if((timeprops_ptr->iter==119)&&(myid==0))
    //printf("myid=%d num_recv=%d ielem=%d\n",myid,num_recv,ielem);

    construct_el(EmNew,recv_array+ielem,NodeTable,myid,&not_used);
    El_Table->add(EmNew->pass_key(),EmNew);
  }

  return;
}

void q_sort_data(double *numbers, void **data, int left, int right)
{
  int int_pivot, l_hold, r_hold;
  double pivot;
  //unsigned pivot_key[2];
  void *pivot_data;

  l_hold = left;
  r_hold = right;
  pivot        = numbers[left]; 
  pivot_data   = data[   left];
  //pivot_key[0] = keys[left][0]; 
  //pivot_key[1] = keys[left][1];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      data[   left] = data[   right];
      //keys[left][0] = keys[right][0];
      //keys[left][1] = keys[right][1];

      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      data[   right] = data[   left];
      //keys[right][0] = keys[left][0];
      //keys[right][1] = keys[left][1];

      right--;
    }
  }
  numbers[left] = pivot;
  data[   left] = pivot_data;
  //keys[left][0] = pivot_key[0];
  //keys[left][1] = pivot_key[1];

  int_pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < int_pivot)
    q_sort_data(numbers, data, left, int_pivot-1);
  if (right > int_pivot)
    q_sort_data(numbers, data, int_pivot+1, right);
}

void q_sort(double *numbers, int left, int right)
{
  int int_pivot, l_hold, r_hold;
  double pivot;

  l_hold = left;
  r_hold = right;
  pivot        = numbers[left]; 
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      right--;
    }
  }
  numbers[left] = pivot;

  int_pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < int_pivot)
    q_sort(numbers, left, int_pivot-1);
  if (right > int_pivot)
    q_sort(numbers, int_pivot+1, right);
}
