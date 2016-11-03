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
 * $Id: node.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <stdio.h>
#include "node.h"
#include "../header/FileFormat.h"
#include "useful_lib.h"

Node::Node(){

 for(int i=0; i<2; i++)
   node_coord[i]=0;


}

void Node::setparameters(int di, double* co)
{

  for(int i=0; i<2; i++)
      node_coord[i]=*(co+i);
  
  nodeid=di;
  written=0;
}

void Node::determine_max_min(double max[], double min[])

{
  if(node_coord[0]>max[0]) max[0]=node_coord[0];
  if(node_coord[1]>max[1]) max[1]=node_coord[1];
  if(node_coord[0]<min[0]) min[0]=node_coord[0];
  if(node_coord[1]<min[1]) min[1]=node_coord[1];

}


void Node::determine_the_key(unsigned nkey, double* max, double *min, unsigned ma[], unsigned mi[])
{

  extern void fhsfc2d_(double* , unsigned* , unsigned* );

  double norm_coord[2];
  norm_coord[0]=(node_coord[0] - *min)/(*max-*min);
  norm_coord[1]=(node_coord[1] - *(min+1))/(*(max+1)-*(min+1));

  fhsfc2d_(norm_coord, &nkey, key);


  if(key[0]>ma[0] || (key[0]==ma[0] && key[1]>ma[1])) {ma[0]=key[0]; ma[1]=key[1];}
  if(key[0]<mi[0] || (key[0]==mi[0] && key[1]<mi[1])) {mi[0]=key[0]; mi[1]=key[1];}
}

void Node::write_node_data_bin(FILE *fp)
{
  if (!written){
    fwriteU(fp,key[0]);
    fwriteU(fp,key[1]);
#ifdef WRITEDOUBLEASFLOAT
    fwriteF(fp,node_coord[0]);
    fwriteF(fp,node_coord[1]);
#else
    fwriteD(fp,node_coord[0]);
    fwriteD(fp,node_coord[1]);
#endif
  }
  written=1;
}

void Node::write_node_data(ofstream* out)
{
  if (!written)
    *out<<key[0]<<' '<<key[1]<<' '<< setprecision(20) <<node_coord[0]<<' '<< setprecision(20) << node_coord[1]<<'\n';
  written=1;
}

void Node::clear_written_flag()
{
  written=0;
}

void Node::put_element_array_loc(int in)
{
  if(in < 0) {
    element_array_loc[0] = element_array_loc[1] = -1;
  }
  else {
    if(element_array_loc[0] == -1)
      element_array_loc[0] = in;
    else if(element_array_loc[1] == -1)
      element_array_loc[1] = in;
    else 
      cout<<"Something wrong with node.element_array_loc!"<<endl;
    
  }
  return;
}

int Node::get_element_array_loc(int in)
{
  if(in == element_array_loc[0] && element_array_loc[1] >= 0)
    return(element_array_loc[1]);
  if(in == element_array_loc[1] && element_array_loc[0] >= 0)
    return(element_array_loc[0]);
  
  return(-2);
}
