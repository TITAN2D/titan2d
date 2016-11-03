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
 * $Id: node.C 156 2007-06-27 20:33:08Z dkumar $ 
 */
//#define DEBUG_SAVE_NODE

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#ifndef NODE_C
#define NODE_C
#include "../header/node.h"
#include "../gisapi/GisApi.h"
#include "../header/properties.h"

#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include <assert.h>

Node:: Node() {
  printf("creating a node without setting its values\n");
}

Node:: Node(unsigned* keyi, double* coordi, MatProps* matprops_ptr)
{
  int      i;
  id     =0;  /* keith added this so save_node wouldn't write an uninitialized
		 variable and valgrind wouldn't flag an error.  id is used in
		 ../repartition/BSFC_update_and_send_elements.C */
  order  =0;  /* keith added this for no good reason, so if you don't know 
		 that this is right find out, keith isn't responsible for 
		 any errors it may cause because right now the node order 
		 is unused, it shouldn't be used, but I wanted to see if 
		 assigning it zero gets rid of what looks like a "self 
		 contained" memory error
	      */
  nextptr =0;
  preptr = 0;
  sol    = 0;
  sol_deleted = 0;
  
  info   = INIT;
	 
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  zero_flux();
  // find the max resolution of the GIS info and then get the elevation at this node
  double resolution = 0;
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation(%d) r=%g (x,y)=(%g,%g) e=%g\n",i,resolution,xcoord,ycoord,elevation);
    exit(1);
  }    
  elevation = elevation/matprops_ptr->LENGTH_SCALE;
/*  if((unsigned) 1548032885 == key[0])
    printf("created the missing node...\n"); */

/*
  if((coord[0]==0)||(coord[1]==0)){
    printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
    assert(coord[0]*coord[1]);
  }
*/
}

Node::Node(unsigned* keyi, double* coordi, int inf, int ord, MatProps* matprops_ptr)  //for refined
{
  int      i;
  id     =0;  /* keith added this so save_node wouldn't write an uninitialized
		 variable and valgrind wouldn't flag an error.  id is used in
		 ../repartition/BSFC_update_and_send_elements.C */

  nextptr =0;
  preptr = 0;
  sol    = 0;
  sol_deleted = 0;
  
  info   = INIT;
  
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  
  info=inf;
  order=ord;
  //geoflow stuff
  zero_flux();
  // find the max resolution of the GIS info and then get the elevation at this node
  double resolution = 0;
  i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation\n");
    exit(1);
  }  
  elevation = elevation/matprops_ptr->LENGTH_SCALE;
/*  if((unsigned) 1548032885 == key[0])
    printf("created the missing node111111...\n");*/

/*
  if((coord[0]==0)||(coord[1]==0)){
    printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
    assert(coord[0]*coord[1]);
  }
*/
  return;
}

Node:: Node(unsigned* keyi, double* coordi, int inf, 
	    int ord, double elev,int yada)
{
  int      i;
  id     =0;  /* keith added this so save_node wouldn't write an uninitialized
		 variable and valgrind wouldn't flag an error.  id is used in
		 ../repartition/BSFC_update_and_send_elements.C */  

  nextptr =0;
  preptr = 0;
  sol    = 0;
  sol_deleted = 0;
  
  info   = INIT;
  
  for(i=0; i<DIMENSION; i++)
    coord[i] = *(coordi+i);
  
  for(i=0; i<KEYLENGTH; i++)
    key[i] = *(keyi+i);
  dof[0] = INIT;
  dof[1] = INIT;
  
  info=inf;
  order=ord;
  //geoflow stuff
  zero_flux();
  elevation = elev;
  /*
  if((coord[0]==0)||(coord[1]==0)){
    printf("inode=%d node={%u,%u} coord=(%20g,%20g)\n",yada,key[0],key[1],coord[0],coord[1]);
    assert(coord[0]*coord[1]);
  }
  */
  return;
}


void Node::set_parameters(int inf, int ord)
{
  info = inf;
  order = ord;
/*  if(key[0] == (unsigned) 3197207111) {
    int mmmyid;
    MPI_Comm_rank(MPI_COMM_WORLD, &mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    printf("changing info and ord %u %u on %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",key[0], key[1], mmmyid);
    printf("======================================================================\n");
    printf("======================================================================\n");

  }*/
}
void Node::putinfo(int in)
{
  info = in;
/*  if(key[0] == (unsigned) 2962355296) {
    int mmmyid;
    MPI_Comm_rank(MPI_COMM_WORLD, &mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    printf("changing info %u %u on %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",key[0], key[1], mmmyid);
    printf("?????????????????????????????????????????????????????? \n");
    printf("?????????????????????????????????????????????????????? \n");
    
  }*/
}


Node:: ~Node() { 
  if(sol != NULL) 
    delete []sol;
}

void Node::zero_flux() {
  int i;
  for(i=0;i<NUM_STATE_VARS;i++)
    flux[i] = refinementflux[i]=0.0;
}


void Node::set_elevation(MatProps* matprops_ptr){  
  double resolution = 0;
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  int i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation\n");
    exit(1);
  }    
  elevation = elevation/matprops_ptr->LENGTH_SCALE;
  

}

void Node::save_node(FILE* fp) {

  FourBytes  temp4;
  EightBytes temp8;
  unsigned writespace[13];

  int Itemp=0, itemp;
  for(itemp=0;itemp<KEYLENGTH;itemp++) {
    writespace[Itemp++]=key[itemp];    
  }
  assert(Itemp==2);
#ifdef DEBUG_SAVE_NODE
  FILE *fpdb=fopen("save_node.debug","w");
  fprintf(fpdb,"key=%u %u\n",key[0],key[1]); 
#endif

  for(itemp=0;itemp<DIMENSION;itemp++) {
    temp8.d=coord[itemp];
    writespace[Itemp++]=temp8.u[0];
    writespace[Itemp++]=temp8.u[1];
  }
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"coord=%g %g\n",coord[0],coord[1]); 
#endif
  assert(Itemp==6);

  temp4.i=id;
  writespace[Itemp++]=temp4.u;
  assert(Itemp==7);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"id=%d\n",id); 
#endif

  temp4.i=info;
  writespace[Itemp++]=temp4.u;
  assert(Itemp==8);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"info=%d\n",info); 
#endif

  /* these are Legacy and are not used
  temp4.i=order;
  writespace[Itemp++]=temp4.u;
  assert(Itemp==9);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"order=%d\n",order); 
#endif

  temp4.i=dof[0];
  writespace[Itemp++]=temp4.u;
  temp4.i=dof[1];
  writespace[Itemp++]=temp4.u;
  assert(Itemp==11);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"dof=%d %d\n",dof[0],dof[1]); 
#endif

  temp4.i=glnum;
  writespace[Itemp++]=temp4.u;
  assert(Itemp==12);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"glnum=%d\n",glnum); 
#endif

  temp4.i=reconstructed;
  writespace[Itemp++]=temp4.u;
  assert(Itemp==13);
#ifdef DEBUG_SAVE_NODE
  fprintf(fpdb,"reconstructed=%d\n",reconstructed); 
#endif
  */

#ifdef DEBUG_SAVE_NODE
  fclose(fpdb);
#endif
  fwrite(writespace,sizeof(unsigned),Itemp,fp);

  return;
}

Node::Node(FILE* fp, MatProps* matprops_ptr) {

  sol=0;           //never USED anywhere in TITAN except when initialized
  sol_deleted=0;   //never USED anywhere in TITAN except when initialized
  nextptr=preptr=0;//never USED anywhere in TITAN except when initialized

  FourBytes  temp4;
  EightBytes temp8;
  //unsigned readspace[13];
  unsigned readspace[8];
  int Itemp=0, itemp;

  //fread(readspace,sizeof(unsigned),13,fp);
  fread(readspace,sizeof(unsigned),8,fp);

  //KEYLENGTH should be 2 but put it in a loop to make it generic.
  for(itemp=0;itemp<KEYLENGTH;itemp++) {
    key[itemp]=readspace[Itemp++];
  }
  assert(Itemp==2);

  //DIMENSION should be 2 but put it in a loop to make it generic.
  for(itemp=0;itemp<DIMENSION;itemp++) {
    temp8.u[0]=readspace[Itemp++];
    temp8.u[1]=readspace[Itemp++];
    coord[itemp]=temp8.d;
  }
  assert(Itemp==6);

  temp4.u=readspace[Itemp++];
  id=temp4.i;
  assert(Itemp==7);

  temp4.u=readspace[Itemp++];
  info=temp4.i;
  assert(Itemp==8);

  /* these are legacy and are not used
  temp4.u=readspace[Itemp++];
  order=temp4.i;
  assert(Itemp==9);

  temp4.u=readspace[Itemp++];
  dof[0]=temp4.i;
  temp4.u=readspace[Itemp++];
  dof[1]=temp4.i;
  assert(Itemp==11);

  temp4.u=readspace[Itemp++];
  glnum=temp4.i;
  assert(Itemp==12);

  temp4.u=readspace[Itemp++];
  reconstructed=temp4.i;
  assert(Itemp==13);
  */
  // find the max resolution of the GIS info and then get the elevation at this node
  double resolution = 0;
  double xcoord = coord[0]*(matprops_ptr->LENGTH_SCALE);
  double ycoord = coord[1]*(matprops_ptr->LENGTH_SCALE);
  int i = Get_max_resolution(&resolution);
  if(i != 0) {
    printf("error in Get_max_resolution\n");
    exit(1);
  }    
  i = Get_elevation(resolution, xcoord, ycoord, &elevation);
  if(i != 0) {
    printf("error in Get_elevation\n");
    exit(1);
  }    
  elevation = elevation/matprops_ptr->LENGTH_SCALE;

  return;
}


#endif





