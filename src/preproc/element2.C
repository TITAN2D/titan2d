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
 * $Id: element2.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include "element.h"
#include "../header/FileFormat.h"
#include "useful_lib.h"


Element::Element(){
  int i;
  for(i=0;i<4;i++)
    neighbor[i] = NULL;
  opposite_brother = NULL;
}


void Element::setparameters(int di, Node* nodes[], int mat, int* elm_loc_in)
{
  int i;
  elementid=di;
  for(i=0; i<8; i++)
    element_nodes[i]=nodes[i];

  for(i=0; i<4; i++)
    {
      neighbor[i]=0;
      for(int j=0; j<2; j++)
	boundary[i][j]=0;
    }

  material=mat;
  elm_loc[0] = elm_loc_in[0];
  elm_loc[1] = elm_loc_in[1];
  if(elm_loc[0] %2 == 0) {
    if(elm_loc[1] %2 == 0)
      which_son = 0;
    else
      which_son = 3;
  }
  else  {
    if(elm_loc[1] %2 == 0)
      which_son = 1;
    else
      which_son = 2;
  }
}


void Element::order_nodes()
{
  double v[3][2];
  int count0=0, count1=0, count2=0, count3=0, count4=0, count5=0;

      for(int i=0; i<3; i++)
	for(int j=0; j<2; j++)

	  v[i][j]=element_nodes[i+1]->node_coord[j]-element_nodes[0]->node_coord[j];
    
      double v01=v[0][0]*v[1][1]-v[0][1]*v[1][0];
      double v12=v[1][0]*v[2][1]-v[1][1]*v[2][0];//The cross products
      double v20=v[2][0]*v[0][1]-v[2][1]*v[0][0];
 
      if(v01<0 && v12<0 && v20>0){ case5(); count5++;}
      
      else if(v01>0 && v12>0 && v20<0) count0++;

      else if(v01>0 && v12<0 && v20<0){case1(); count1++;}
      
      else if(v01<0 && v12>0 && v20<0){case2(); count2++;}
      
      else if(v01>0 && v12<0 && v20>0) {case3();count3++;}
      
      else if(v01<0 && v12>0 && v20>0){case4();count4++;}
       
}


void Element::case5()
{

  Node* temp;

  temp=element_nodes[1];
  element_nodes[1]=element_nodes[3];
  element_nodes[3]=temp;
  temp=element_nodes[5];
  element_nodes[5]=element_nodes[7];
  element_nodes[7]=temp;

}   


void Element::case1()
{
  Node* temp;

  temp=element_nodes[2];
  element_nodes[2]=element_nodes[3];
  element_nodes[3]=temp;
	  
}

void Element::case2()
{

  Node* temp;

  temp=element_nodes[1];
  element_nodes[1]=element_nodes[2];
  element_nodes[2]=temp;
}

void Element::case3()
{
  Node* temp;

  temp=element_nodes[1];
  element_nodes[1]=element_nodes[2];
  element_nodes[2]=element_nodes[3];
  element_nodes[3]=temp;

}

void Element::case4()
{
  Node* temp;

  temp=element_nodes[1];
  element_nodes[1]=element_nodes[3];
  element_nodes[3]=element_nodes[2];
  element_nodes[2]=temp;
}


void Element::create_m_node(double* max, double *min)
{

  Node* newPtr=new Node();
  
  for(int j=0; j<2; j++)
   {
    newPtr->node_coord[j]=0;
    for(int i=0; i<4; i++)
      newPtr->node_coord[j]+=element_nodes[i]->node_coord[j]/4;
   }
  newPtr->written=0;
  
  element_nodes[8]=newPtr;
  element_nodes[8]->nodeid=-elementid;

  if(element_nodes[8]->node_coord[0]>*max) *max=element_nodes[8]->node_coord[0];
  if(element_nodes[8]->node_coord[1]>*(max+1)) *(max+1)=element_nodes[8]->node_coord[1];
  if(element_nodes[8]->node_coord[0]<*min) *min=element_nodes[8]->node_coord[0];
  if(element_nodes[8]->node_coord[1]<*(min+1)) *(min+1)=element_nodes[8]->node_coord[1];
  
  return;
}


/*void Element::determine_the_key(unsigned nkey, double* max, double *min)
{

  extern void fhsfc2d_(double* , unsigned* , unsigned* );
  double norm_coord[2];
  norm_coord[0]=(element_nodes[8]->node_coord[0] - *min)/(*max-*min);
  norm_coord[1]=(element_nodes[8]->node_coord[1] - *(min+1))/(*(max+1)-*(min+1));

  fhsfc2d_(norm_coord, &nkey, key);

}*/

void Element::determine_neighbors(int array_loc, Element* el)
{
  int i, neigh_loc;

  for(i=4;i<8;i++)
    {
      neigh_loc = element_nodes[i]->get_element_array_loc(array_loc);
      if(neigh_loc >= 0)
	neighbor[i-4] = (el+neigh_loc);
    }
}

//*************************finding the opposite brother***********************
void Element::determine_opposite_brother()
{
  if(which_son != 3) {
    if(neighbor[which_son+1] != NULL ) {
      if(which_son != 2)
	opposite_brother = neighbor[which_son+1]->get_neighbors(which_son+2);
      else
	opposite_brother = neighbor[which_son+1]->get_neighbors(0);
    }
  }
  else {
    if(neighbor[0] != NULL)
      opposite_brother = neighbor[0]->get_neighbors(1);
  }


}

Element* Element::get_neighbors(int side) 
{
  return neighbor[side];
}

void Element::myproc(int np, int i, int ec)
{
  myprocess=i/(ec/np);

  if(i>=(ec/np)*np) 
    myprocess=np-1;
}

void Element::write_element_data_bin(FILE *fp){
  int i;
  for(i=0; i<9; i++){
    fwriteU(fp,element_nodes[i]->key[0]);
    fwriteU(fp,element_nodes[i]->key[1]);}

  for(i=0; i<4; i++)
    if(neighbor[i]!=0){
      fwriteI(fp,neighbor[i]->myprocess);
      fwriteU(fp,neighbor[i]->element_nodes[8]->key[0]);
      fwriteU(fp,neighbor[i]->element_nodes[8]->key[1]);}
    else fwriteI(fp,-1);

  for(i=0; i<4; i++)
    if(boundary[i][0]!=0){
      fwriteI(fp,boundary[i][0]->type);
      /* note the boundaries SHOULD be written as floats even when
	 you __DON'T__ want to WRITEDOUBLEASFLOAT because in boundary.h
	 for hpfem the values to be read into are FLOATs not doubles */
      fwriteF(fp,boundary[i][0]->x_value);
      fwriteF(fp,boundary[i][0]->y_value);}
    else fwriteI(fp,-1);

  for(i=0; i<4; i++)
    if(boundary[i][1]!=0){
      fwriteI(fp,boundary[i][1]->type);
      /* note the boundaries SHOULD be written as floats even when
	 you __DON'T__ want to WRITEDOUBLEASFLOAT because in boundary.h
	 for hpfem the values to be read into are FLOATs not doubles */
      fwriteF(fp,boundary[i][1]->x_value);
      fwriteF(fp,boundary[i][1]->y_value);}
    else fwriteI(fp,-1);
   
  fwriteI(fp,elm_loc[0]);
  fwriteI(fp,elm_loc[1]);

  if(opposite_brother == NULL) {
    fwriteI(fp,0);
    fwriteI(fp,0);}
  else{
    fwriteU(fp,*(opposite_brother->pass_key()));
    fwriteU(fp,*(opposite_brother->pass_key()+1));}

  fwriteI(fp,material);

  return;
}


void Element::write_element_data(ofstream* out){
  int i;
  for(i=0; i<9; i++)
    //    *out<<elementid<<'\n';
    *out<<element_nodes[i]->key[0]<<' '<<element_nodes[i]->key[1]<<' ';
  *out<<'\n';
  for(i=0; i<4; i++)
    {
    if(neighbor[i]!=0) *out<<neighbor[i]->myprocess<<' '<<neighbor[i]->element_nodes[8]->key[0]<<' '<<neighbor[i]->element_nodes[8]->key[1]
		      <<' ';
    else *out<<-1<<' ';
    }
  /*  *out<<elementid<<'\n';
  for(int i=0; i<4; i++)
    {
    if(neighbor[i]!=0) *out<<neighbor[i]->elementid
		      <<' '<<neighbor[i]->myprocess<<' ';
    else *out<<-1<<' ';
    }*/ //Just 4 test

  *out<<'\n';
  for(i=0; i<4; i++)
    {
      if(boundary[i][0]!=0) *out<<boundary[i][0]->type<<' '<<boundary[i][0]->x_value<<' '<<boundary[i][0]->y_value<<' ';
      else *out<<-1<<' ';
    }
  *out<<'\n';

  for(i=0; i<4; i++)
    {
      if(boundary[i][1]!=0) *out<<boundary[i][1]->type<<' '<<boundary[i][1]->x_value<<' '<<boundary[i][1]->y_value<<' ';
      else *out<<-1<<' ';
    }

  *out<<'\n';
  *out<<elm_loc[0]<<' '<<elm_loc[1]<<'\n';
  if(opposite_brother == NULL) {
    *out<<0<<' '<<0<<'\n';
  }
  else
    *out<<*(opposite_brother->pass_key())<<' '<<*(opposite_brother->pass_key()+1)<<'\n';
  *out<<material<<endl;
      
}


void Element::reset_written_flag()
{

  for(int i=0; i<8; i++) if(element_nodes[i]->written!=0) element_nodes[i]->written=0;

}

void Element::set_boundary(Boundary *b)
{
  for(int j=0; j<4; j++)
     if(b->node->nodeid==element_nodes[j+4]->nodeid)
       {
	 if(b->type==-2)
	   boundary[j][0]=b;
	 
	 if(b->type==-3)
	   boundary[j][1]=b;
       }
  return;
}

void Element::find_boundary(int cc, int fc, Boundary b[])
{
  for(int i=0; i<cc+fc; i++)
    for(int j=4; j<8; j++)
      if(b[i].node->nodeid==element_nodes[j]->nodeid)
	{
	  if(b[i].type==-2)
	    boundary[j-4][0]=&b[i];

	  if(b[i].type==-3)
	    boundary[j-4][1]=&b[i];

	  break;
	}    
}

unsigned* Element::pass_key() 
{
  unsigned* mykey = element_nodes[8]->get_key();

  return mykey;
}
