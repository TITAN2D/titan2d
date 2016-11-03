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
 * $Id: setup_geoflow.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"


void setup_geoflow(HashTable* El_Table, HashTable* NodeTable, int myid, 
		   int nump,MatProps* matprops_ptr,TimeProps *timeprops_ptr) 
{

  int i;
  int num_buckets = El_Table->get_no_of_buckets();
  int num_node_buckets = NodeTable->get_no_of_buckets();
  /* zero out the fluxes for all of the nodes */
  HashEntryPtr* buck = NodeTable->getbucketptr();
  for(i=0; i<num_node_buckets; i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Node* Curr_Node=(Node*)(currentPtr->value);
	    Curr_Node->zero_flux();
	      
	    currentPtr=currentPtr->next;      	    
	  }
      }

  
  /* put the coord for the center node in the element */
  buck = El_Table->getbucketptr();
  for(i=0; i<num_buckets; i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    int refined = Curr_El->get_refined_flag();
 	    if(Curr_El->get_adapted_flag()>0)//if this is a refined element don't involve!!!
	      {
		Curr_El->find_positive_x_side(NodeTable);
		Curr_El->calculate_dx(NodeTable); 
		Curr_El->calc_topo_data(matprops_ptr);
		Curr_El->calc_gravity_vector(matprops_ptr);
	      }
	      
	    currentPtr=currentPtr->next;      	    
	  }
      }


  /* transfer ghost elements to proper processors */
  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);

  /* calculate d_gravity array for each element */
  buck = El_Table->getbucketptr();
  for(i=0; i<num_buckets; i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
	    int refined = Curr_El->get_refined_flag();
 	    if(Curr_El->get_adapted_flag()>0)//if this is a refined element don't involve!!!
	      {
		Curr_El->calc_d_gravity(El_Table);
	      }
	    
	    currentPtr=currentPtr->next;      	    
	  }
      }

  return;
}
