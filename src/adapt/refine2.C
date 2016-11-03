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
 * $Id: refine2.C 127 2007-06-07 19:48:25Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include"../header/hpfem.h"

extern void fhsfc2d_(double, unsigned, unsigned);
extern void hsfc2d(unsigned* , unsigned* , unsigned* );
extern void create_new_node(int, int, int, HashTable*, Node*[], 
			    unsigned[][2], int, int*, int, int, MatProps*);


//()---new node numbering


//  3---(14)--6---(15)--2
//  |         |         |
//  |         |         |
// (11) (3)  (12) (2)  (13)
//  |         |         |
//  |         |         |
//  7---(9)---E---(10)--5
//  |         |         |
//  |         |         |
// (6)  (0)  (7)  (1)  (8)
//  |         |         |
//  |         |         |
//  0---(4)---4---(5)---1


//max & min coordinates need to be passed later now only for L-shape!!!
//only 4 one step because of the info FLAG!!!
//if the new node is on INTERFACE flag will be -1

void refine(Element* EmTemp, HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, MatProps* matprops_ptr)
{ 
  //printf("refining element %u %u \n",*(EmTemp->pass_key()), *(EmTemp->pass_key()+1));
  int which;
  Node *n1, *n2, *n3, *n4;
  
  unsigned* KeyTemp;
  Node* NodeTemp[9];
  unsigned NewNodeKey[16][KEYLENGTH];
  Element* Quad9P;
  int numprocs, myid, i;
  Element* neigh_elm;
  unsigned* neigh_node_key;
  int RefinedNeigh=0;
  int info;
  int other_proc=0;
  int boundary;
  int order[5];
  int NewOrder[4][5];
  
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int elm_loc[2], my_elm_loc[2];
  int* fth_elm_loc = EmTemp->get_elm_loc();
  elm_loc[0] = 2*fth_elm_loc[0];
  elm_loc[1] = 2*fth_elm_loc[1];  
  
  KeyTemp = EmTemp->getNode();
  
  for(i=0; i<8; i++)//-- corners and sides
    {
      NodeTemp[i]=(Node*) HT_Node_Ptr->lookup(KeyTemp+i*KEYLENGTH);
      assert(NodeTemp[i] );
    }
  
  for(i=0; i<5; i++)
    order[i]=*(EmTemp->get_order()+i);
  
  /*filling up the new order array
    str: side orders remain;
    newsides get the higher order of the already existing sides
    bubbles get the order of the old bubble
  */
  for(i=0; i<4; i++)//-- orders for the 4 sons
    {
      int a=i-1;
      if(a==-1) 
	a=3;
      NewOrder[i][a]=order[a];
      NewOrder[i][i]=order[i];
      NewOrder[i][4]=order[4];
    }
  
  //for the new internal sides they get the order of the max side order
  for(i=0; i<4; i++)
    {
      int a=i+1;
      int b=i+2;
      int c=i-1;      
      if(i==2) b=0;
      if(i==3) {a=0; b=1;}//in the case of the 3rd element
      if(c==-1) c=3;//in the case of the 0th element
      
      NewOrder[i][a]=NewOrder[i][b]=order[0];

      for(int j=1; j<4; j++)
	if(order[j]>NewOrder[i][a]) NewOrder[i][a]=NewOrder[i][b]=order[j];
    }

  NodeTemp[8]=(Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());//-- bubble

  //SIDE 0
  if(*(EmTemp->get_neigh_proc())==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()) <= EmTemp->get_gen()) {
    RefinedNeigh = 0;
    info=S_S_CON;
    if(*(EmTemp->get_neigh_proc())!=myid) 
      {
	other_proc=1;
	info=-1;
      }

    which=4;
    //---Fourth new node---
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
    
    create_new_node(which, 0, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[0], matprops_ptr);
    
    //---Fourth old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[4]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[4]->putinfo(-1);
    else 
      NodeTemp[4]->putinfo(S_C_CON);
    
    //---Fifth new node---
    which=5;
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
    
    create_new_node(which, 1, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh,
		    boundary, order[0], matprops_ptr);
  }
  else {

    //Keith Added this if
    if((*(EmTemp->get_neigh_proc())!=myid)||
       ((*(EmTemp->get_neigh_proc()+4)!=myid)&&
	(*(EmTemp->get_neigh_proc()+4)!=-2)) 
       )
      other_proc=1;
    else
      other_proc=0;

    // fourth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors());
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[4][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[4]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
    //fourth old node
    NodeTemp[4]->putinfo(CORNER);
    if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
      NodeTemp[4]->putinfo(-1);
    // fifth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+4*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[5][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[5]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);


  }

//+++++++++++++++++++++++++++SIDE1

  if(*(EmTemp->get_neigh_proc()+1)==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()+1) <= EmTemp->get_gen()) {
    RefinedNeigh = 0;
    info=S_S_CON;  
    if(*(EmTemp->get_neigh_proc()+1)!=myid)// && *(EmTemp->get_neigh_proc()+1)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;

    //---Eight new node---
    which=8;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
    
    create_new_node(which, 1, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[1], matprops_ptr);
    
    //---Fifth old node---
    if(RefinedNeigh || boundary)
      NodeTemp[5]->putinfo(CORNER);
    else { 
      if(other_proc) 
	NodeTemp[5]->putinfo(info);
      else 
	NodeTemp[5]->putinfo(S_C_CON);
    }
    
    //---Thirteenth new node---
    which=13;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
    
    create_new_node(which, 2, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		    &RefinedNeigh, boundary, order[1], matprops_ptr);
  }
  else {
    //Keith Added this if
    if((*(EmTemp->get_neigh_proc()+1)!=myid)||
       ((*(EmTemp->get_neigh_proc()+5)!=myid)&&
	(*(EmTemp->get_neigh_proc()+5)!=-2)) 
       )
      other_proc=1;
    else
      other_proc=0;

    // eighth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[8][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
    // fifth old node
    NodeTemp[5]->putinfo(CORNER);
    if(other_proc) //ERROR: other_proc is set based on side 0 neigbor not being more refined or never set, we never checked to see if the more refined neighbor was on another processor 
      NodeTemp[5]->putinfo(-1);
    // thirteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+5*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[13][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
  }

  //+++++++++++++++++++++++++++SIDE2

  if(*(EmTemp->get_neigh_proc()+2)==-1) 
    boundary=1;
  else 
    boundary=0;
  
  if(boundary == 1 || *(EmTemp->get_neigh_gen()+2) <= EmTemp->get_gen()) {
    info=S_S_CON;
    
    if(*(EmTemp->get_neigh_proc()+2)!=myid)// && *(EmTemp->get_neigh_proc()+2)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;

    RefinedNeigh=0;
    
    //---Fourteenth new node---
    which=14;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
    
    create_new_node(which, 3, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[2], matprops_ptr);
    
    //---Sixth old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[6]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[6]->putinfo(-1);
    else
      NodeTemp[6]->putinfo(S_C_CON);
    
    //---Fifteenth new node---
    which=15;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
    
    create_new_node(which, 2, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[2], matprops_ptr);
  }
  else {
    //Keith Added this if
    if((*(EmTemp->get_neigh_proc()+2)!=myid)||
       ((*(EmTemp->get_neigh_proc()+6)!=myid)&&
	(*(EmTemp->get_neigh_proc()+6)!=-2)) 
       )
      other_proc=1;
    else
      other_proc=0;


    // fourteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+6*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[14][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[14]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
    // sixth old node
    NodeTemp[6]->putinfo(CORNER);
    if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
      NodeTemp[6]->putinfo(-1);
    // fifteenth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+2*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[15][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[15]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
  }

  //+++++++++++++++++++++++++++SIDE 3

  if(*(EmTemp->get_neigh_proc()+3)==-1) 
    boundary=1;
  else 
    boundary=0;

  if(boundary == 1 || *(EmTemp->get_neigh_gen()+3) <= EmTemp->get_gen()) {
    info=S_S_CON;
    
    if(*(EmTemp->get_neigh_proc()+3)!=myid)  //&& *(EmTemp->get_neigh_proc()+3)>0) 
      {
	other_proc=1;
	info=-1;
      }
    else 
      other_proc=0;  
    
    RefinedNeigh=0;
    
    //---Sixth new node----
    which=6;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
    
    create_new_node(which, 0, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[3], matprops_ptr);
    
    //---Seventh old node---
    if(RefinedNeigh || boundary) 
      NodeTemp[7]->putinfo(CORNER);
    else if(other_proc) 
      NodeTemp[7]->putinfo(-1);
    else 
      NodeTemp[7]->putinfo(S_C_CON);
    
    //---Eleventh new node---
    which=11;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
    
    create_new_node(which, 3, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		    &RefinedNeigh, boundary, order[3], matprops_ptr);
  }
  else {
    //Keith Added this if
    if((*(EmTemp->get_neigh_proc()+3)!=myid)||
       ((*(EmTemp->get_neigh_proc()+7)!=myid)&&
	(*(EmTemp->get_neigh_proc()+7)!=-2)) 
       )
      other_proc=1;
    else
      other_proc=0;

    // sixth new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+7*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[6][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
    // seventh old node
    NodeTemp[7]->putinfo(CORNER);
    if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
      NodeTemp[7]->putinfo(-1);
    // eleventh new node
    neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+3*KEYLENGTH);
    i = 0;
    which = -1;
    while(i<4 && which == -1) {
      if(compare_key(neigh_elm->get_neighbors()+i*KEYLENGTH, EmTemp->pass_key()))
	which = i;
      i++;
    }
    assert(which != -1);
    neigh_node_key = neigh_elm->getNode();
    for(i=0;i<KEYLENGTH;i++)
      NewNodeKey[11][i] = neigh_node_key[i+(which+4)*KEYLENGTH];
    n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
    if(neigh_elm->get_refined_flag() == 0 || neigh_elm->get_refined_flag() == GHOST)
      n1->putinfo(SIDE);
    //else if(neigh_elm->get_refined_flag()==GHOST)
    //n1->putinfo(-1);
    else
      n1->putinfo(S_C_CON);
  }
  //++++++++++++++++INTERNAL SIDE NODES 7, 8OLD, 12, 9, 10

  RefinedNeigh=0;
  boundary=0;
  info=SIDE;

  //---Seventh new node---

  which=7;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
  n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());
 
  create_new_node(which, 4, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		  &RefinedNeigh, boundary, NewOrder[0][1], matprops_ptr);

  NodeTemp[8]->putinfo(CORNER);//changing the old bubble


  //---Twelwth new node---

  which=12;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);

  create_new_node(which, 6, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[3][1], matprops_ptr);


  //---Ninth new node---

  which=9;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 7, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[0][2], matprops_ptr);

  //---Tenth new node---

  which=10;
  // geoflow info
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);

  create_new_node(which, 5, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[1][2], matprops_ptr);


  //+++++++++++++++++++THE NEW BUBBLES 0, 1, 2, 3

  info=BUBBLE;

  //---0th new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
  which=0;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode());
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*4);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[0][4], matprops_ptr);
  


  //---1st new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
  which=1;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[1][4], matprops_ptr);



  //---2nd new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
  which=2;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*2);
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*5);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info,
		  &RefinedNeigh, boundary, NewOrder[2][4], matprops_ptr);


  //---3rd new node---

  NodeTemp[0]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
  NodeTemp[1]= (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
  which=3;
  n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*6);
  n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*3);
  n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->getNode()+KEYLENGTH*7);

  create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, 
		  &RefinedNeigh, boundary, NewOrder[3][4], matprops_ptr);


  //---NEW ELEMENTS---

  unsigned nodes[9][KEYLENGTH];
  unsigned neigh[8][KEYLENGTH];
  unsigned* orig_neighbors=EmTemp->get_neighbors();
  int* orig_neigh_proc=EmTemp->get_neigh_proc();
  int neigh_proc[8];
  BC* bcptr=NULL;
  BC* orig_bcptr=EmTemp->get_bcptr();
  int generation=EmTemp->get_gen()+1;
  int* orig_neigh_gen=EmTemp->get_neigh_gen();
  int neigh_gen[4];
  int material=EmTemp->get_material();

  double coord[DIMENSION];
  //---0th new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+i);
      nodes[1][i]=*(KeyTemp+4*KEYLENGTH+i);
      nodes[2][i]=*((EmTemp->pass_key())+i);
      nodes[3][i]=*(KeyTemp+7*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[4][i];
      nodes[5][i]=NewNodeKey[7][i];
      nodes[6][i]=NewNodeKey[9][i];
      nodes[7][i]=NewNodeKey[6][i];
      nodes[8][i]=NewNodeKey[0][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);
  //neighbors
 
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=*(orig_neighbors+i); //Why is this ok if not ok for 3 down
      neigh[1][i]=neigh[5][i]=NewNodeKey[1][i];
      neigh[2][i]=neigh[6][i]=NewNodeKey[3][i];
      if(*(EmTemp->get_neigh_proc()+7)!=-2) 
	neigh[3][i]=neigh[7][i]=*(orig_neighbors+7*KEYLENGTH+i); //This should be okay no matter what
      else
	neigh[3][i]=neigh[7][i]=*(orig_neighbors+3*KEYLENGTH+i); //This is only ok if neigh_proc==-2
    }

  //process of the neighbors

  neigh_proc[0]=*(orig_neigh_proc);
  neigh_proc[1]=myid;
  neigh_proc[2]=myid;
  if(*(orig_neigh_proc+7)!=-2) 
    neigh_proc[3]=*(orig_neigh_proc+7);//depending if the neighboring element is already refined
  else 
    neigh_proc[3]=*(orig_neigh_proc+3);

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=*orig_neigh_gen;
  neigh_gen[1]=generation;
  neigh_gen[2]=generation;
  neigh_gen[3]=*(orig_neigh_gen+3);


  //boundary conditions
  if(orig_bcptr && (orig_bcptr->type[0] || orig_bcptr->type[3])) //else bcptr is a NULL pointer by default, ERROR this should crash if orig_bcptr==NULL
    {
      bcptr=new BC;
      bcptr->type[0]=orig_bcptr->type[0];
      bcptr->type[3]=orig_bcptr->type[3];
      for(i=0; i<2; i++) 
	for(int j=0; j<2; j++)
	  {
	    bcptr->value[0][i][j]=orig_bcptr->value[0][i][j];
	    bcptr->value[3][i][j]=orig_bcptr->value[3][i][j];
	  }	
    }
 
  double err = (*(EmTemp->get_el_error()))*.5; //added by jp oct11
  double sol = (*(EmTemp->get_el_solution()))*.5;//added by jp oct11
  // son 0 can use elm_loc
  int iwetnodefather=EmTemp->get_iwetnode();
  double Awetfather=EmTemp->get_Awet();
  double dpson[2];
  dpson[0]=*(EmTemp->get_drypoint()+0)*2+0.5;
  dpson[1]=*(EmTemp->get_drypoint()+1)*2+0.5;  
  Quad9P=new Element(nodes, neigh, neigh_proc, bcptr, generation, elm_loc, 
		     &NewOrder[0][0], neigh_gen, material, EmTemp, coord,
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr,
		     iwetnodefather, Awetfather, dpson);
  double* state_vars=Quad9P->get_state_vars();

  Quad9P->put_which_son(0);//--by jp, 0 means son 0

  Quad9P->putel_sq(sol, err);//added by jp oct11
  Element* old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    old_elm->put_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
    old_elm->void_bcptr();
    HT_Elem_Ptr->remove(old_elm->pass_key(),1,stdout,myid,16);
    delete old_elm;
  }
    
  HT_Elem_Ptr->add(nodes[8], Quad9P);



  //---1st new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+4*KEYLENGTH+i);
      nodes[1][i]=*(KeyTemp+1*KEYLENGTH+i);
      nodes[2][i]=*(KeyTemp+5*KEYLENGTH+i);
      nodes[3][i]=*((EmTemp->pass_key())+i);
      nodes[4][i]=NewNodeKey[5][i];
      nodes[5][i]=NewNodeKey[8][i];
      nodes[6][i]=NewNodeKey[10][i];
      nodes[7][i]=NewNodeKey[7][i];
      nodes[8][i]=NewNodeKey[1][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
 
  for(i=0; i<KEYLENGTH; i++)
    {
      if(*(EmTemp->get_neigh_proc()+4)!=-2)
	neigh[0][i]=neigh[4][i]=*(orig_neighbors+4*KEYLENGTH+i); //this should be ok now matter what
      else
	neigh[0][i]=neigh[4][i]=*(orig_neighbors+0*KEYLENGTH+i); //this is only ok if neigh_proc==-2
      neigh[1][i]=neigh[5][i]=*(orig_neighbors+1*KEYLENGTH+i); 
      neigh[2][i]=neigh[6][i]=NewNodeKey[2][i];
      neigh[3][i]=neigh[7][i]=NewNodeKey[0][i];
     
    }


  //process of the neighbors

  neigh_proc[0]=(*(orig_neigh_proc+4)!=-2) ? *(orig_neigh_proc+4) : *(orig_neigh_proc);
  neigh_proc[1]=*(orig_neigh_proc+1);
  neigh_proc[2]=myid;
  neigh_proc[3]=myid;

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=*orig_neigh_gen;
  neigh_gen[1]=*(orig_neigh_gen+1);
  neigh_gen[2]=generation;
  neigh_gen[3]=generation;

  bcptr=NULL;
  //boundary conditions
  if(orig_bcptr && (orig_bcptr->type[0] || orig_bcptr->type[1])) //else bcptr is a NULL pointer by default
    {
      bcptr=new BC;
      bcptr->type[0]=orig_bcptr->type[0];
      bcptr->type[1]=orig_bcptr->type[1];
      for(i=0; i<2; i++) 
	for(int j=0; j<2; j++)
	  {
	    bcptr->value[0][i][j]=orig_bcptr->value[0][i][j];
	    bcptr->value[1][i][j]=orig_bcptr->value[1][i][j];
	  }
	
    }
  my_elm_loc[0] = elm_loc[0]+1;
  my_elm_loc[1] = elm_loc[1];
  dpson[0]=*(EmTemp->get_drypoint()+0)*2-0.5;
  dpson[1]=*(EmTemp->get_drypoint()+1)*2+0.5;  
  Quad9P=new Element(nodes, neigh, neigh_proc, bcptr, generation, my_elm_loc, 
		     &NewOrder[1][0], neigh_gen, material, EmTemp, coord,
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr,
		     iwetnodefather, Awetfather, dpson);
  state_vars=Quad9P->get_state_vars();

  Quad9P->put_which_son(1);//--by jp

  Quad9P->putel_sq(sol, err);//added by jp oct11
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    old_elm->put_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
    old_elm->void_bcptr();
    HT_Elem_Ptr->remove(old_elm->pass_key(),1,stdout,myid,17);
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);


  //---2nd new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*((EmTemp->pass_key())+i);
      nodes[1][i]=*(KeyTemp+5*KEYLENGTH+i);
      nodes[2][i]=*(KeyTemp+2*KEYLENGTH+i);
      nodes[3][i]=*(KeyTemp+6*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[10][i];
      nodes[5][i]=NewNodeKey[13][i];
      nodes[6][i]=NewNodeKey[15][i];
      nodes[7][i]=NewNodeKey[12][i];
      nodes[8][i]=NewNodeKey[2][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
  

 
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=NewNodeKey[1][i];
      if(*(EmTemp->get_neigh_proc()+5)!=-2)
	neigh[1][i]=neigh[5][i]=*(orig_neighbors+5*KEYLENGTH+i); //This should be ok no matter what
      else
	neigh[1][i]=neigh[5][i]=*(orig_neighbors+1*KEYLENGTH+i); //this is only ok is neigh_proc==-2
      neigh[2][i]=neigh[6][i]=*(orig_neighbors+2*KEYLENGTH+i);
      neigh[3][i]=neigh[6][i]=NewNodeKey[3][i];

    }


  //process of the neighbors

  neigh_proc[0]=myid;
  neigh_proc[1]=(*(orig_neigh_proc+5)!=-2) ? *(orig_neigh_proc+5) : *(orig_neigh_proc+1);
  neigh_proc[2]=*(orig_neigh_proc+2);
  neigh_proc[3]=myid;

  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=generation;
  neigh_gen[1]=*(orig_neigh_gen+1);
  neigh_gen[2]=*(orig_neigh_gen+2);
  neigh_gen[3]=generation;

  bcptr=NULL;
  //boundary conditions
  if(orig_bcptr && (orig_bcptr->type[1] || orig_bcptr->type[2])) //else bcptr is a NULL pointer by default
    {
      bcptr=new BC;
      bcptr->type[1]=orig_bcptr->type[1];
      bcptr->type[2]=orig_bcptr->type[2];
      for(i=0; i<2; i++) 
	for(int j=0; j<2; j++)
	  {
	    bcptr->value[1][i][j]=orig_bcptr->value[1][i][j];
	    bcptr->value[2][i][j]=orig_bcptr->value[2][i][j];
	  }
      
    }
  my_elm_loc[0] = elm_loc[0] +1;
  my_elm_loc[1] = elm_loc[1] +1;
  dpson[0]=*(EmTemp->get_drypoint()+0)*2-0.5;
  dpson[1]=*(EmTemp->get_drypoint()+1)*2-0.5;  
  Quad9P=new Element(nodes, neigh, neigh_proc, bcptr, generation, my_elm_loc,
		     &NewOrder[2][0], neigh_gen, material, EmTemp, coord,
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr,
		     iwetnodefather, Awetfather, dpson);
  state_vars=Quad9P->get_state_vars();

  Quad9P->put_which_son(2);//--by jp

  Quad9P->putel_sq(sol, err);//added by jp oct11
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    old_elm->put_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
    old_elm->void_bcptr();
    HT_Elem_Ptr->remove(old_elm->pass_key(),1,stdout,myid,18);
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);


  //---3rd new element---


  //the nodes
  
  for(i=0; i<KEYLENGTH; i++)
    {
      nodes[0][i]=*(KeyTemp+7*KEYLENGTH+i);
      nodes[1][i]=*((EmTemp->pass_key())+i);
      nodes[2][i]=*(KeyTemp+6*KEYLENGTH+i);
      nodes[3][i]=*(KeyTemp+3*KEYLENGTH+i);
      nodes[4][i]=NewNodeKey[9][i];
      nodes[5][i]=NewNodeKey[12][i];
      nodes[6][i]=NewNodeKey[14][i];
      nodes[7][i]=NewNodeKey[11][i];
      nodes[8][i]=NewNodeKey[3][i];
    }
  n1 = (Node*) HT_Node_Ptr->lookup(&(nodes[8][0]));
  for(i=0;i<DIMENSION;i++)
    coord[i] = *(n1->get_coord()+i);

  //neighbors
  for(i=0; i<KEYLENGTH; i++)
    {
      neigh[0][i]=neigh[4][i]=NewNodeKey[0][i];
      neigh[1][i]=neigh[5][i]=NewNodeKey[2][i];
      if(*(EmTemp->get_neigh_proc()+6)!=-2)
	neigh[2][i]=neigh[6][i]=*(orig_neighbors+6*KEYLENGTH+i);
      else
	neigh[2][i]=neigh[6][i]=*(orig_neighbors+2*KEYLENGTH+i);
      
      neigh[3][i]=neigh[6][i]=*(orig_neighbors+3*KEYLENGTH+i);    
      
    }


  //process of the neighbors

  neigh_proc[0]=myid;
  neigh_proc[1]=myid;
  neigh_proc[2]=(*(orig_neigh_proc+6)!=-2) ? *(orig_neigh_proc+6) : *(orig_neigh_proc+2);
  neigh_proc[3]=*(orig_neigh_proc+3);
 
  neigh_proc[4]=neigh_proc[5]=neigh_proc[6]=neigh_proc[7]=-2;

  neigh_gen[0]=generation;
  neigh_gen[1]=generation;
  neigh_gen[2]=*(orig_neigh_gen+2);
  neigh_gen[3]=*(orig_neigh_gen+3);
 
  bcptr=NULL;
 //boundary conditions
  if(orig_bcptr && (orig_bcptr->type[2] || orig_bcptr->type[3])) //else bcptr is a NULL pointer by default
    {
      bcptr=new BC;
      bcptr->type[2]=orig_bcptr->type[2];
      bcptr->type[3]=orig_bcptr->type[3];
      for(i=0; i<2; i++) 
	for(int j=0; j<2; j++)
	  {
	    bcptr->value[2][i][j]=orig_bcptr->value[2][i][j];
	    bcptr->value[3][i][j]=orig_bcptr->value[3][i][j];
	  }
	
    }

  my_elm_loc[0] = elm_loc[0];
  my_elm_loc[1] = elm_loc[1] +1;
  dpson[0]=*(EmTemp->get_drypoint()+0)*2+0.5;
  dpson[1]=*(EmTemp->get_drypoint()+1)*2-0.5;  
  Quad9P=new Element(nodes, neigh, neigh_proc, bcptr, generation, my_elm_loc, 
		     &NewOrder[3][0], neigh_gen, material, EmTemp, coord,
		     HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr,
		     iwetnodefather, Awetfather, dpson);
  state_vars=Quad9P->get_state_vars();

  Quad9P->put_which_son(3);//--by jp

  Quad9P->putel_sq(sol, err);//added by jp oct11
  old_elm = (Element*) HT_Elem_Ptr->lookup(Quad9P->pass_key());
  if(old_elm != NULL) {
    old_elm->put_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
    old_elm->void_bcptr();
    HT_Elem_Ptr->remove(old_elm->pass_key(),1,stdout,myid,19);
    delete old_elm;
  }

  HT_Elem_Ptr->add(nodes[8], Quad9P);


  //---CHANGING THE FATHER---
  EmTemp->putson(&NewNodeKey[0][0]);
  // putting in brother info
  for(i=0;i<4;i++) {
    EmTemp = (Element*) HT_Elem_Ptr->lookup(&NewNodeKey[i][0]);
    EmTemp->putbrothers(&NewNodeKey[0][0]);  //was  EmTemp->putbrothers(&NewNodeKey[i][0]);
  }

  return;
}


void create_new_node(int which, int Node1, int Node2, HashTable* HT_Node_Ptr, 
		     Node* NodeTemp[], unsigned NewNodeKey[][KEYLENGTH], int info,
		     int* RefNe, int boundary, int order, MatProps* matprops_ptr)
{
  double NewNodeCoord[2];
  double norm_coord[2];
  unsigned u_norm_coord[2];
  unsigned nkey=2;
  unsigned key[KEYLENGTH];
  Node* NewNode;
  Node* p;
  static double XRange[2]; 
  static double YRange[2];
  int i;
 
  for(i=0; i<2; i++)
    {
      XRange[i]=*(HT_Node_Ptr->get_Xrange()+i);
      YRange[i]=*(HT_Node_Ptr->get_Yrange()+i);
    }
 
  for(i=0; i<2; i++)    
    NewNodeCoord[i]=(*(NodeTemp[Node1]->get_coord()+i) + *(NodeTemp[Node2]->get_coord()+i))*.5; 
  
  norm_coord[0]=(NewNodeCoord[0]-XRange[0])/(XRange[1]-XRange[0]);
  norm_coord[1]=(NewNodeCoord[1]-YRange[0])/(YRange[1]-YRange[0]);

  fhsfc2d_(norm_coord, &nkey, key);

  for(i=0; i<KEYLENGTH; i++)

    NewNodeKey[which][i]=key[i];


  p=(Node*) HT_Node_Ptr->lookup(key);
  
  if(!p)
    
    {
      NewNode = new Node(key, NewNodeCoord, info, order, matprops_ptr);

      HT_Node_Ptr->add(key, NewNode);

      p=NewNode;

    }

  else if(*(p->get_coord())!=NewNodeCoord[0] || *(p->get_coord()+1)!=NewNodeCoord[1])
   {
     short same_key=0;
     assert(same_key);
   }

  else
    *RefNe=1;      
  
  if(*RefNe || boundary)
    p->putinfo(SIDE);
      
  return;
}








