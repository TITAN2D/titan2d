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
 * $Id: BSFC_refine_partition.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "./repartition_BSFC.h"

void BSFC_refine_partition(int* local_balanced_flag, 
			   int *amount_of_used_bits, int num_vert_in_cut,
			   BSFC_VERTEX_PTR vert_in_cut_ptr, 
			   float* work_percent_array, float total_weight,
			   float* global_actual_work_allocated, 
			   int number_of_cuts, 
			   int* ll_bins_head, float* work_prev_allocated,
			   int subbins_per_bin, int* local_balanced_flag_array,
			   int myid, int numprocs)
		      
{
  //printf("proc %d is refining the partition level\n",myid);
  int i=0, j=0, k, current_proc;
  int amount_of_bits;
  float* binned_weight_array;
  int* bin_proc_array;
  int* ll_prev_bins;
  int ll_counter, ll_location, *ll_bins_head_copy;

  /* amount of sub-bins in a bin, probably want this as a passed in parameter */
  int number_of_bins = subbins_per_bin;

  /* check to see that all of the bits of the sfc key 
     have not already been used */
  if(*amount_of_used_bits >= sizeof(unsigned) * KEYLENGTH * 8) {
    //printf("No more refinement is possible in the repartitioning on proc %d.\n", myid);
    *local_balanced_flag = BSFC_BALANCED;
    return;
  }
  
  /*  assume initially that all the partitions on this processor are balanced.
      we will check later on whether any are not balanced */
  *local_balanced_flag = BSFC_BALANCED;

  /* if there are a lot of cuts in a bin, we want the amount of bins
     to be greater than the amount of cuts */
  if(number_of_cuts >= number_of_bins)
    number_of_bins = number_of_cuts + 1;

  /*increase sub-bins so that there is a power of 2 */
  i=0;
  while(number_of_bins > BSFC_pow(2,i))
    i++;
  amount_of_bits = i;
  if(amount_of_bits + *amount_of_used_bits > 8*sizeof(unsigned) * KEYLENGTH)
    amount_of_bits = 8*sizeof(unsigned) * KEYLENGTH - *amount_of_used_bits;
  number_of_bins = BSFC_pow(2,i);
    
  ll_prev_bins = (int*) malloc(sizeof(int) * (number_of_cuts+1));
  
  ll_bins_head_copy = (int*) malloc(sizeof(int) * (number_of_cuts+1));

  for(i=0;i<=number_of_cuts;i++)
    ll_bins_head_copy[i] = -1;
  
  /* loop over all bins that have a cut in them using linklist 
     to find objects in the cut bins */
  for(ll_counter=0;ll_counter<=number_of_cuts;ll_counter++) 
    if((local_balanced_flag_array[ll_counter]==BSFC_NOT_BALANCED)
       && ll_bins_head[ll_counter] != -1) {

      /* calculate new bin numbers for objects that are in a cut bin */
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	vert_in_cut_ptr[ll_location].my_bin = 
	  BSFC_get_array_location(number_of_bins, amount_of_bits, 
				  *amount_of_used_bits, (vert_in_cut_ptr+ll_location));
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }  
      
      binned_weight_array = (float*) malloc(number_of_bins*sizeof(float));
      for(i=0;i<number_of_bins;i++)
	binned_weight_array[i] = 0;
      
      /* fill up the weight array */
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {  
	binned_weight_array[vert_in_cut_ptr[ll_location].my_bin] += 
	  vert_in_cut_ptr[ll_location].lb_weight;
	  ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }
      
      bin_proc_array = (int*) malloc(sizeof(int) * number_of_bins);
      i = number_of_bins-1;

      current_proc = vert_in_cut_ptr[ll_bins_head[ll_counter]].destination_proc;
      float scanned_work_prev_allocated = work_prev_allocated[current_proc];
      while(i >= 0 && current_proc > -1) {
	scanned_work_prev_allocated += binned_weight_array[i];
	bin_proc_array[i] = current_proc;
	if(scanned_work_prev_allocated > work_percent_array[current_proc]* total_weight) {
	  global_actual_work_allocated[current_proc] = scanned_work_prev_allocated;      
	  bin_proc_array[i] = -current_proc;
	  current_proc--;
	  // the while statement is if there is more than 1 cut in a bin
	  while(current_proc >= 0 &&
		scanned_work_prev_allocated > work_percent_array[current_proc]*total_weight) {
	    global_actual_work_allocated[current_proc] = scanned_work_prev_allocated;      
	    current_proc--;
	  }
	}
	i--;
      }
                  
      /* specify which processor an object belongs to,
	 we will know this because we know what bin an object 
	 belongs to and we know what processor a bin belongs to */
      
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	if(bin_proc_array[vert_in_cut_ptr[ll_location].my_bin] >= 0) {
	  vert_in_cut_ptr[ll_location].cut_bin_flag = BSFC_NO_CUT;
	  vert_in_cut_ptr[ll_location].destination_proc = 
	    bin_proc_array[vert_in_cut_ptr[ll_location].my_bin];
	  ll_location =  vert_in_cut_ptr[ll_location].next_sfc_vert_index;
	}
	else {	/* if this object is in a bin with a cut... */
	  vert_in_cut_ptr[ll_location].cut_bin_flag = BSFC_CUT;
	  vert_in_cut_ptr[ll_location].destination_proc = 
	    -bin_proc_array[vert_in_cut_ptr[ll_location].my_bin];
	  if(ll_bins_head_copy[number_of_cuts-myid+ //array continued on the next line
			      vert_in_cut_ptr[ll_location].destination_proc] != -1) {
	    int ll_next = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
	    vert_in_cut_ptr[ll_location].next_sfc_vert_index = 
	      ll_bins_head_copy[number_of_cuts-myid+vert_in_cut_ptr[ll_location].destination_proc];
	    ll_bins_head_copy[number_of_cuts-myid+vert_in_cut_ptr[ll_location].destination_proc] = ll_location;
	    ll_location = ll_next;
	  }
	  else {
	    int ll_next = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
	    vert_in_cut_ptr[ll_location].next_sfc_vert_index = -1;
	    ll_bins_head_copy[number_of_cuts-myid+vert_in_cut_ptr[ll_location].destination_proc] = ll_location;
	    /* calculate work_prev_allocated for this new partition */
	    work_prev_allocated[vert_in_cut_ptr[ll_location].destination_proc] = 
	      work_prev_allocated[ll_counter-number_of_cuts+myid];
	    for(i=vert_in_cut_ptr[ll_location].my_bin+1;i<number_of_bins;i++)
	      work_prev_allocated[vert_in_cut_ptr[ll_location].destination_proc] 
		+= binned_weight_array[i];
	    ll_location = ll_next;
	  }
	}
      }
      free(binned_weight_array);
      free(bin_proc_array);
    }
  
  for(i=0;i<=number_of_cuts;i++) 
    ll_bins_head[i] = ll_bins_head_copy[i];

  free(ll_prev_bins);
  free(ll_bins_head_copy);

  *amount_of_used_bits += amount_of_bits;
  
  /* check which partitions that are not balanced */
  for(i=0;i<=number_of_cuts;i++) 
    if(ll_bins_head[i] != -1 || local_balanced_flag_array[i] != BSFC_BALANCED)
      local_balanced_flag_array[i] =
	BSFC_find_imbalance(work_percent_array, 
			    global_actual_work_allocated[(myid+i-number_of_cuts)],
			    total_weight, myid+i-number_of_cuts, numprocs);
  
  /* check if any of the partitions are not balanced */
  *local_balanced_flag = BSFC_BALANCED;
  i=0;
  while(*local_balanced_flag == BSFC_BALANCED && i<=number_of_cuts) {
    *local_balanced_flag = local_balanced_flag_array[i];
    i++;
  }
  
  /* check the partitions to see if any more improvement can be made on them */
  if(*local_balanced_flag == BSFC_NOT_BALANCED) {
    for(i=0;i<=number_of_cuts;i++)
      if(ll_bins_head[i] != -1) {
	/* check if there is only 1 object in this bin. if there is, 
	   no further bin refinement will improve load-balance */
	if(vert_in_cut_ptr[ll_bins_head[i]].next_sfc_vert_index == -1) {
	  ll_bins_head[i] = -1;
	  local_balanced_flag_array[i] = BSFC_BALANCED;
	  //printf("Bin refinement cannot improve load balance on proc %d\n", myid);
	}
      }
    /* check again if any of the partitions are not balanced */
    *local_balanced_flag = BSFC_BALANCED;
    j=0;
    while(*local_balanced_flag == BSFC_BALANCED && j<=number_of_cuts) {
      *local_balanced_flag = local_balanced_flag_array[j];
      j++;
    }
  }
  
  return;
}

