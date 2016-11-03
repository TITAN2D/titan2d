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
 * $Id: element_weight.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

/*! element_weight() cycles through the element Hashtable (listing of all 
 *  elements) and for each element (that has not been refined this iteration 
 *  and is not a ghost_element) calls Element member function 
 *  Element::calc_flux_balance() (which returns a double precision value 
 *  representing the weight that an element is assigned based on the 
 *  magnitude of its net mass/momentum fluxes). Note that this value is 
 *  adjusted to give non-zero weight even to elements with zero pile-heights
 *  The cumulative weights (along with a count of the evaluated elements) 
 *  are stored in sub_weight[]; based on this, the return value for this 
 *  function is calculated and stored in global_weight[] (i.e. the sum of 
 *  sub_weight[] from all processors).
 */
double element_weight(HashTable* El_Table, HashTable* NodeTable, int myid, int nump) {
  int i,j,k, counter;
  double tiny = GEOFLOW_TINY;
  int el_counter = 0;
  double evalue = 1;
  double sub_weight[2] = {0,0}; // second number is to keep track of the number of objects
  
  //-------------------go through all the elements of the subdomain and  
  //-------------------find the edge states

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    if((Curr_El->get_adapted_flag()>0)) //Keith added this
            {
	      //if this element doesn't belong on this processor don't involve!!! 
	      Curr_El->calc_flux_balance(NodeTable);
	      if(*(Curr_El->get_state_vars())>GEOFLOW_TINY)
              {
		sub_weight[1] += 1;
		sub_weight[0] += *(Curr_El->get_el_error());
	      }
	      else if(Curr_El->get_adapted_flag()==BUFFER){
		sub_weight[1] += 0.1;
		sub_weight[0] += *(Curr_El->get_el_error())*0.1;
	      }	      
	    }
	    
	    currentPtr=currentPtr->next;      	    
	  }
      }

  double global_weight[2];
  i = MPI_Allreduce(sub_weight, global_weight, 2,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  //global_weight[0] = (global_weight[0]-global_weight[1])/global_weight[1];
  if(global_weight[1]>0.0)
    global_weight[0]=(global_weight[0])/global_weight[1]; //to protect from division by zero
  else
    global_weight[0]=1.0; //just to make it nonzero

  return global_weight[0];
}
