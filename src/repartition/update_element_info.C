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
 * $Id: update_element_info.C 153 2007-06-27 20:30:26Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
/*
  destroy_element
  create_element
  same_proc
  diff_proc
  construct_el
  check_neighbor_info

 */


#include "../header/hpfem.h"


void same_proc(Element *r_element, HashTable* HT_Elem_Ptr, 
	       int target_proc, int side);
void diff_proc(Element* r_element, HashTable* HT_Elem_Ptr, 
	       int new_proc, int side,  ELinkPtr* EL_head);

//! construct_el is a friend function of the Element class that fills an element with information it receives in a variable of the ElemPack class from an MPI call
void construct_el(Element* newelement, ElemPack* elem2, 
		  HashTable* HT_Node_Ptr, int myid, 
		  double* e_error);

void check_neighbor_info(Element* newelement, HashTable* HT_Elem_Ptr, 
			 int myid);


/*! destroy_element() is a friend function of the Element and Node classes that does the following
 *
 *  1.  Update the neighbor_proc of the neighbors
 *  1.1 if neighbor is at the same proc---->ok
 *  1.2 if neighbor is at the target proc-->done when the element 
 *     is created in its new subdomain
 *  1.3 if neighbor is at a 3rd proc------->these elements are linked 
 *     for later communication
 *
 *  2.  Remove element from the hashtable
 *
 *  3.  Remove some nodes..........later not now
 */
void destroy_element(void *r_element_in, 
		     HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		     int target_proc, ELinkPtr* EL_head){

  int myid, numprocs, i;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  Element *r_element = (Element *) r_element_in;

  if(!r_element->get_refined_flag())//if parent don't care about neighbor info
    {
      for(i=0; i<8; i++)
	{
	  if(r_element->neigh_proc[i] == myid)
	    same_proc(r_element, HT_Elem_Ptr,target_proc, i);
	  
	  else if(r_element->neigh_proc[i] != myid && 
		  r_element->neigh_proc[i] != target_proc && r_element->neigh_proc[i] >= 0)
	    diff_proc(r_element,HT_Elem_Ptr, 
		      target_proc, i, EL_head);
	}
    }

  HT_Elem_Ptr->remove(r_element->key,1,stdout,myid,26);


}


//! create_element() is a friend function of the Element and Node classes. After receiving an ElemPack, create_element() instances a new element, calls construct_el() to transfer data from ElemPack to the new element, and inserts the new element into the Hashtable. Don't call this if s_flag is 0 (original repartitioning scheme)
void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, 
		    HashTable* HT_Node_Ptr, int myid, double* e_error)
{

  Element* newelement = new Element();
  
  construct_el(newelement, elem2, HT_Node_Ptr, myid, e_error);

  HT_Elem_Ptr->add(newelement->pass_key(), newelement);

  if(!newelement->get_refined_flag())//if parent .... don't care about neighbor info
    check_neighbor_info(newelement, HT_Elem_Ptr, myid);

}

// bsfc repartitioning scheme
void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, 
		    HashTable* HT_Node_Ptr, int myid)
{
  Element* newelement = new Element();
  double e_error = 0;
  construct_el(newelement, elem2, HT_Node_Ptr, myid, &e_error);
  
  Element* EmTemp = (Element*) HT_Elem_Ptr->lookup(newelement->pass_key());
  if(EmTemp != NULL) {// update this element
    // first check that it is the same element
    // if the generation, son number, and parent keys are the same then will
    // be the same element unless the initial mesh was too fine and 
    // multiple elements had the same initial key
    assert(elem2->generation == EmTemp->get_gen());
    if(elem2->generation > 0) {
      assert(elem2->which_son == EmTemp->get_which_son());
    }
    // the same element...
    HT_Elem_Ptr->remove(elem2->key,1,stdout,myid,27);
    delete EmTemp;
  }
  if(!((*(newelement->pass_key()+0)==0)&&
       (*(newelement->pass_key()+1)==0)))
    HT_Elem_Ptr->add(newelement->pass_key(), newelement);
  else delete newelement;
  //printf("processor %d just added element %u %u\n",myid, elem2->key[0], elem2->key[1]);
  if(myid == 3 && elem2->key[0] == 0)
    e_error = 0;

  return;
}


void same_proc(Element* r_element, HashTable* HT_Elem_Ptr, 
	       int target_proc, int side)
{
  Element* Neighbor;
  Neighbor=(Element*)HT_Elem_Ptr->lookup(r_element->get_neighbors()+side*KEYLENGTH);
  //r_element->put_neigh_gen(side, Neighbor->get_gen());// added by jp 9.30

  int which = Neighbor->which_neighbor(r_element->pass_key());
  int gen = r_element->get_gen();// added by jp 9.30

  Neighbor->change_neighbor_process(which, target_proc);
  Neighbor->put_neigh_gen(which, gen);// added by jp 9.30

}

void diff_proc(Element* r_element, HashTable* HT_Elem_Ptr, 
	       int new_proc, int side,  ELinkPtr* EL_head)
{

  //create a linked list
  ELinkPtr EL_new;
  ELinkPtr EL_temp;
  
  EL_new = new ElementLink(r_element->pass_key(), 
			   (r_element->get_neighbors()+side*KEYLENGTH),
			   *(r_element->getassoc()+side), new_proc);

  if(!(*EL_head)) *EL_head = EL_new;

  else
    {
      EL_new->next = *EL_head;
      (*EL_head)->pre = EL_new;      
      *EL_head = EL_new;
    }
}

void construct_el(Element* newelement, ElemPack* elem2, 
		  HashTable* HT_Node_Ptr, int myid, double* e_error)
{ 
  Node* node;
  int i, j;
  newelement->myprocess = myid;
  newelement->generation = elem2->generation;
  newelement->opposite_brother_flag = elem2->opposite_brother_flag;
  newelement->material=elem2->material;

  for(i=0; i<8; i++)
    {
      newelement->neigh_proc[i]=elem2->neigh_proc[i];
      newelement->neigh_gen[i]=elem2->neigh_gen[i];
    }
  for(i=0; i<5; i++)
    newelement->order[i] = elem2->order[i];

  newelement->ndof = elem2->ndof;
  newelement->no_of_eqns = elem2->no_of_eqns;
  newelement->refined = elem2->refined;
  newelement->adapted = elem2->adapted;
  newelement->which_son = elem2->which_son;
  newelement->new_old = elem2->new_old;
  for(i=0;i<4;i++)
    for(j=0;j<KEYLENGTH;j++)
      newelement->brothers[i][j] = elem2->brothers[i][j];

  for(i=0; i<KEYLENGTH; i++)
    newelement->key[i] = elem2->key[i];
    
  for(i=0; i<8; i++)
    for(int j=0; j<KEYLENGTH; j++)
      {
	newelement->node_key[i][j] = elem2->node_key[i][j];
	newelement->neighbor[i][j] = elem2->neighbor[i][j];
	if(i<4)newelement->son[i][j] = elem2->son[i][j];

      }
  for(i=0; i<EQUATIONS; i++)
    {
      newelement->el_error[i] = elem2->el_error[i];
      newelement->el_solution[i] = elem2->el_solution[i];
    }

  //and the node info -- ignore some info if this is just getting a parent from another processor...
  for(i=0; i<8; i++) {
    if(elem2->n_coord[i][0]*elem2->n_coord[i][1]==0) {
      printf("myid=%d elem2->key={%u,%u} elem2->coord=(%20g,%20g) inode=%d node->key={%u,%u} node->coord=(%20g,%20g)\n",myid,elem2->key[0],elem2->key[1],elem2->n_coord[8][0],elem2->n_coord[8][1],i,elem2->node_key[i][0],elem2->node_key[i][1],elem2->n_coord[i][0],elem2->n_coord[i][1]);
    }

      node = (Node*) HT_Node_Ptr->lookup(elem2->node_key[i]);
      if(!node) {
	  node = new Node(elem2->node_key[i], elem2->n_coord[i],
			  elem2->n_info[i], elem2->n_order[i], elem2->node_elevation[i],i);
	  
	  HT_Node_Ptr->add(elem2->node_key[i], node);
	}
      else {
	//because of storing all the node but not updating the 
	//info and order if the node was not previously in the subdomain
	//check if the sfc is screwed
	if(*(node->get_coord())!=elem2->n_coord[i][0] ||
	   *(node->get_coord()+1)!=elem2->n_coord[i][1])  {
	  printf("myid=%d\n  pack  elem(x,y)=(%20g,%20g)\n exist elem(x,y)=(%20g,%20g)\n\n",myid,elem2->n_coord[i][0],elem2->n_coord[i][1],*(node->get_coord()),*(node->get_coord()+1)); fflush(stdout);
	    int screwd=0;
	    assert(screwd);
	  }
	if(newelement->refined == 0)  // only update if this is from an active element
	  node->set_parameters(elem2->n_info[i], elem2->n_order[i]);
      }
  }
  
  node = (Node*) HT_Node_Ptr->lookup(elem2->key);
  if(!node) {
      node = new Node(elem2->key, elem2->n_coord[8],
		      elem2->n_info[8], elem2->n_order[8], elem2->node_elevation[8],8);
      
      HT_Node_Ptr->add(elem2->key, node);
    }
  else if(newelement->refined != 0) // only update if this is from an active element
    node->set_parameters(elem2->n_info[8], elem2->n_order[8]);
    
  //The boundary conditions

  if((elem2->bc)!=0)
    {
      BC* newbc = new BC;
      for(i=0; i<4; i++)
	{
	  newbc->type[i]=elem2->bc_type[i];
	  for(int j=0; j<2; j++)
	   for(int k=0; k<2; k++)
	    newbc->value[i][j][k] = elem2->bc_value[i][j][k];
	}
       newelement->bcptr = newbc;
    }

  else newelement->bcptr = 0;

  *e_error = newelement->el_error[0];

  //geoflow stuff
  newelement->positive_x_side = elem2->positive_x_side;
  newelement->elevation = elem2->elevation;
  for(i=0;i<DIMENSION;i++) {
    newelement->coord[i]       = elem2->n_coord[8][i];
    newelement->dx[i]          = elem2->dx[i];
    newelement->eigenvxymax[i] = elem2->eigenvxymax[i];
    newelement->kactxy[i]      = elem2->kactxy[i];
    newelement->zeta[i]        = elem2->zeta[i];
    newelement->curvature[i]   = elem2->curvature[i];
    newelement->d_gravity[i]   = elem2->d_gravity[i];
  }
  for(i=0;i<NUM_STATE_VARS;i++) {
    newelement->state_vars[i] = elem2->state_vars[i];
    newelement->Influx[    i] = elem2->Influx[    i];
  }
  for(i=0;i<3;i++)
    newelement->gravity[i] = elem2->gravity[i];
  for(i=0;i<DIMENSION*NUM_STATE_VARS;i++)
    newelement->d_state_vars[i] = elem2->d_state_vars[i];
  newelement->shortspeed = elem2->shortspeed;
  newelement->lam        = elem2->lam;
  newelement->lb_weight  = elem2->lb_weight;
  newelement->elm_loc[0] = elem2->elm_loc[0];
  newelement->elm_loc[1] = elem2->elm_loc[1];
  
  newelement->iwetnode   = elem2->iwetnode;
  newelement->Awet       = elem2->Awet;
  newelement->Swet       = elem2->Swet;
  newelement->drypoint[0]= elem2->drypoint[0];
  newelement->drypoint[1]= elem2->drypoint[1];

  return;
}

void check_neighbor_info(Element* newelement, HashTable* HT_Elem_Ptr, int myid)
{

  int* neigh_proc = newelement->getassoc();
  unsigned* neigh_key = newelement->get_neighbors();
  Element* neighbor;
  int which;

  for(int i=0; i<8; i++)
    {
      if(*(neigh_proc+i) == myid)
	{
	  neighbor = (Element*)HT_Elem_Ptr->lookup(neigh_key+i*KEYLENGTH);
	  which = neighbor->which_neighbor(newelement->pass_key());
	  neighbor->change_neighbor_process(which, myid);
	}

    }

}






//for the 3rd party involved in updating
void diff_proc1_2(int counter, NeighborPack packed_neighbor_info[], HashTable* HT_Elem_Ptr)
{

  Element* element;
  int which;
  for(int i=0; i<counter; i++)
    {
      element = (Element*)HT_Elem_Ptr->lookup(packed_neighbor_info[i].targetkey);
      assert(element);
      which = element->which_neighbor(packed_neighbor_info[i].elkey);
      element->change_neighbor_process(which, packed_neighbor_info[i].new_proc);

    }

}
