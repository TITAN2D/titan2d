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
 * $Id: BSFC_combine_elements.C 135 2007-06-07 20:15:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "./repartition_BSFC.h"

// routine figures out how to "bunch" together elements that cannot be put on different processors because of
// a constrained node 
// if destination_proc > -1 , update the element proc, otherwise combine the elements together
void BSFC_combine_elements(int side, Element *EmTemp, HashTable *HT_Elem_Ptr, 
			   HashTable *HT_Node_Ptr, int destination_proc)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // printf("%u %u constrained for combine on %d\n", *(EmTemp->pass_key()), *(EmTemp->pass_key()+1), myid);
  
  Element* elm_bunch[4];
  int orig_side_number = side;
  Node* nd;
  int i, j;
  for(i=0;i<4;i++)
    elm_bunch[i] = NULL;

  elm_bunch[EmTemp->which_son] = EmTemp;
  
  if(side == EmTemp->which_son)
    switch(EmTemp->which_son) {
    case 0:
      elm_bunch[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
      break;
    case 1:
      elm_bunch[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
      break;
    case 2:
      elm_bunch[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
      break;
    case 3:
      elm_bunch[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
      break;
    }
  else 
    switch(EmTemp->which_son) {
    case 0:
      elm_bunch[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
      break;
    case 1:
      elm_bunch[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
      break;
    case 2:
      elm_bunch[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
      break;
    case 3:
      elm_bunch[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
      break;
    }
  // check the next side
  int next_side = side + 1;
  if(EmTemp->which_son == 0 && next_side == 1)
    next_side = 3;
  if(next_side <4) {
    nd = (Node*) HT_Node_Ptr->lookup(&(EmTemp->node_key[next_side+4][0]));
    if(nd->getinfo() == S_S_CON) {
      if(next_side == EmTemp->which_son)
	switch(EmTemp->which_son) {
	case 0:
	  elm_bunch[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	  break;
	case 1:
	  elm_bunch[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	  break;
	case 2:
	  elm_bunch[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	  break;
	case 3:
	  elm_bunch[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	  break;
	}
      else 
	switch(EmTemp->which_son) {
	case 0:
	  elm_bunch[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	  break;
	case 1:
	  elm_bunch[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	  break;
	case 2:
	  elm_bunch[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	  break;
	case 3:
	  elm_bunch[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	  break;
	}
    }
  }
  // look at "opposite" element and check for constrained nodes there
  Element* orig_elem = EmTemp; // save for later
  int next_elm = EmTemp->get_which_son()+2;
  if(next_elm > 3)
    next_elm -= 4;
  Element* elm_bunch2[4];
  for(i=0;i<4;i++)
    elm_bunch2[i] = NULL;
  if(orig_elem->get_opposite_brother_flag() == 0)
    orig_elem->find_opposite_brother(HT_Elem_Ptr);
  if(orig_elem->get_opposite_brother_flag() == 1)
    EmTemp = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[next_elm][0]));
  else 
    EmTemp = NULL;

  if(EmTemp != NULL && EmTemp->get_refined_flag() == 0) {
    side = -1;
    i = 4;
    while(i<8 && side == -1) {
      nd = (Node*) HT_Node_Ptr->lookup(&(EmTemp->node_key[i][0]));
      if(nd->getinfo() == S_S_CON)
	side = i-4;
      i++;
    }

    if(side != -1) {
      elm_bunch2[EmTemp->which_son] = EmTemp;
      
      if(side == EmTemp->which_son)
	switch(EmTemp->which_son) {
	case 0:
	  elm_bunch2[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	  break;
	case 1:
	  elm_bunch2[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	  break;
	case 2:
	  elm_bunch2[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	  break;
	case 3:
	  elm_bunch2[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	  break;
	}
      else 
	switch(EmTemp->which_son) {
	case 0:
	  elm_bunch2[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	  break;
	case 1:
	  elm_bunch2[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	  break;
	case 2:
	  elm_bunch2[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	  break;
	case 3:
	  elm_bunch2[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	  break;
	}
      // check the next side
      next_side = side + 1;
      if(EmTemp->which_son == 0 && next_side == 1)
	next_side = 3;
      if(next_side <4) {
	nd = (Node*) HT_Node_Ptr->lookup(&(EmTemp->node_key[next_side+4][0]));
      	if(nd->getinfo() == S_S_CON) {
	  if(next_side == EmTemp->which_son)
	    switch(EmTemp->which_son) {
	    case 0:
	      elm_bunch2[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	      break;
	    case 1:
	      elm_bunch2[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	      break;
	    case 2:
	      elm_bunch2[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	      break;
	    case 3:
	      elm_bunch2[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	      break;
	    }
	  else 
	    switch(EmTemp->which_son) {
	    case 0:
	      elm_bunch2[3] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[3][0]));
	      break;
	    case 1:
	      elm_bunch2[0] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[0][0]));
	      break;
	    case 2:
	      elm_bunch2[1] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[1][0]));
	      break;
	    case 3:
	      elm_bunch2[2] = (Element*) HT_Elem_Ptr->lookup(&(EmTemp->brothers[2][0]));
	      break;
	    }
	}
      }
    }
  }

  //combine the elements into one "big" element
  j = 0;
  i = 0;
  while(i<4 && j == 0) {
    if(elm_bunch[i] != NULL && elm_bunch2[i] != NULL)
      j = 1;
    i++;
  }
  if(j == 1)
    for(i=0;i<4;i++)
      if(elm_bunch[i] == NULL)
	elm_bunch[i] = elm_bunch2[i];
  
  if(destination_proc >= 0) {  //done load-balancing -- now just update the element process (myprocess)
    for(i=0;i<4;i++)
      if(elm_bunch[i] != NULL) {
	elm_bunch[i]->myprocess = destination_proc;;
      }
  }
  else {
    // combine the info
    double sum_lb_weights = 0;
    j = 0;
    for(i=0;i<4;i++)
      if(elm_bunch[i] != NULL) {
	//printf("new_old value for %u %u\n", *(elm_bunch[i]->pass_key()), *(elm_bunch[i]->pass_key()+1));
	j++;
	sum_lb_weights += elm_bunch[i]->lb_weight;
	assert(elm_bunch[i]->new_old== BSFC_NEW);
	elm_bunch[i]->new_old=BSFC_OLD;
      }
    orig_elem->lb_weight = sum_lb_weights;
    orig_elem->new_old = j;
    // put in lb_key representing this "bunch" of elements 
    if(j<2)
      assert(j>1);
    if(j==2) {  // 2 elements get shared side node key as the representative key...
      for(i=0;i<KEYLENGTH;i++)
	orig_elem->lb_key[i] = orig_elem->node_key[orig_side_number+4][i];
    }
    else { // 3 or more elements get parent's bubble key as the representative key...
      if(elm_bunch[0] != NULL) {
	for(i=0;i<KEYLENGTH;i++)
	  orig_elem->lb_key[i] = elm_bunch[0]->node_key[2][i];
      }
      else {
	for(i=0;i<KEYLENGTH;i++)
	  orig_elem->lb_key[i] = elm_bunch[1]->node_key[3][i];
      }
    }
  }
 
  return;
}

