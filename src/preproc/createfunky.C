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
 * new file created 2003/08/28 17:20:48 by kdalbey 
 *
 */


#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "boundary.h"
#include "element.h"
#include "../header/FileFormat.h"
#include "GisApi.h"
#include "node.h"
#include "useful_lib.h"

#define MIN_NONSEQ_REPART

/* found in preprocess.C */
void Read_material_data(int *material_count, char ***materialnames,
			double **lambda, double **mu);


void createfunky(int NumProc, char *GISDbase, char *location, 
		 char *mapset, char *topomap, 
		 int havelimits, double limits[4],
		 int *node_count, Node **node, 
		 int *element_count, Element **element,
		 int *force_count, int *constraint_count, 
		 Boundary **boundary,
		 int *material_count, char ***materialnames, 
		 double **lambda, double **mu){

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

  if(!Initialize_GIS_data(GISDbase, location, mapset, topomap)){
    double xmin, xmax, ymin, ymax;
    double res;
    int havematmap=0; //default is we don't have one
    FILE *fp;

    /* read in the material names and properties from file "frict.data"
       this data is otherwise read at the same time as funky.bin or funky.dat 
       see the function Read_material_data() in preprocess.C */
    Read_material_data(material_count, materialnames, lambda, mu);

    Get_max_resolution(&res);
    Get_xmax(res, &xmax);
    Get_xmin(res, &xmin);
    Get_ymax(res, &ymax);
    Get_ymin(res, &ymin);  
    
    printf("The GIS region is:  x = %e to %e   y = %e to %e.\n",xmin, xmax, ymin, ymax);


    double matres; /* need to put this outside this if statement so if 
		      downstream will see it */
    if(*material_count>1){
      char *matmap=(char *) malloc((strlen(topomap)+5)*sizeof(char));
      strcpy(matmap,topomap);
      strcat(matmap+strlen(topomap),"_Mat");

      if(Initialize_Raster_data(GISDbase, location, mapset, matmap)){
	printf("createfunky.C cannot open material property map \"%s\".\n  Check if the file exists.\n",matmap);
	exit(1);}

      free(matmap); //just the name not the map
      
      double xminmat, xmaxmat, yminmat, ymaxmat;
      Get_raster_resolution(&matres);
      Get_raster_xmax(matres, &xmaxmat);
      Get_raster_xmin(matres, &xminmat);
      Get_raster_ymax(matres, &ymaxmat);
      Get_raster_ymin(matres, &yminmat);  
      
      printf("The material property GIS region is:  x = %e to %e   y = %e to %e.\n",xminmat, xmaxmat, yminmat, ymaxmat);
     
      /* compare the limits of the two GIS maps ... 
	 the sanctioned method is that the 2 maps are the same size
	 but just in case they're not */
      if(xminmat>xmin) xmin=xminmat;
      if(yminmat>ymin) ymin=yminmat;
      if(xmaxmat>xmax) xmax=xmaxmat;
      if(ymaxmat>ymax) ymax=ymaxmat;
    }

    if(havelimits){
      /* check the values for user input min and max coordinates ... */
      if((limits[0]>xmin)&&(limits[0]<xmax)) xmin = limits[0];
      if((limits[1]>ymin)&&(limits[1]<ymax)) ymin = limits[1];
      if((limits[2]>xmin)&&(limits[2]<xmax)) xmax = limits[2];
      if((limits[3]>ymin)&&(limits[3]<ymax)) ymax = limits[3];}

    printf("The simulation region is:  x = %e to %e   y = %e to %e.\n",xmin, xmax, ymin, ymax);   
    
    double xlength = 0.98*(xmax-xmin);
    double ylength = 0.98*(ymax-ymin);
    xmin    = xmin+0.1*xlength;
    xmax    = xmax-0.1*xlength;
    xlength = xmax-xmin;
    ymin    = ymin+0.1*ylength;
    ymax    = ymax-0.1*ylength;
    ylength = ymax-ymin;

    int NumDim=2; //x,y

#ifdef MIN_NONSEQ_REPART
    //the number of subdomain divisions per dimension is such that each
    //processor will originaly own exactly 2 or 4 subdomains in two
    //dimensions (or in three dimensions 2, 4, or 8 subdomains per proc). 
    //(mdj)2007-04-11 int NumDimDiv=(int) pow(2,floor(log(pow(NumProc,1.0/NumDim))/log(2.0))+1);
    int NumDimDiv=(int) pow(2,floor(log(pow(NumProc,1.0/NumDim))/log(2.0))+3);
    int NumCellPerSub=pow(2,NumDim);
#else
    //disable the initial gridding trick that minimizes non-sequential 
    //processor repartitioning

    int NumDimDiv=1;
    int NumCellPerSub=8*NumProc;
#endif

    double dX=(xmax-xmin)/NumDimDiv; //subdomain X length
    double dY=(ymax-ymin)/NumDimDiv; //subdomain Y length

    int nx=2;
    double dx=dX/nx;
    int ny=2;
    double dy=dY/ny;

    if(dy>=dx) {
      ny=ceil(dY/dx);
      dy=dY/ny; }
    else{
      nx=ceil(dX/dy);
      dx=dX/nx; 
    }

    double aspect_ratio_toler=0.05;

    while(fabs(1.0-min(dx,dy)/max(dx,dy))>aspect_ratio_toler) {
      
      ny++;
      dy=dY/ny;
      nx=(int) (dX/dy+0.5); //round to nearest whole number
      dx=dX/nx;
    }
    nx*=NumDimDiv;
    ny*=NumDimDiv;
  
  
    //create funky here
    //int nx;
    int i,j, icounter;
    int **elem_node, *elem_mat, **elem_loc;
    int *bound_id;
    double **bound_xy;
    double *x, *y, **xy;

    /*
    nx = (int) (ny*(0.00001+xlength/ylength));
    if(nx <= 0){
      nx = 10;
      ny = (int) (nx*(ylength/xlength));}
    */

    //set up corner, side and bubble nodes

    //first x nodes
    x=CAllocD1(2*nx+1);
    for(i=0; i<2*nx+1; i++)
      x[i] = xlength*i/(2.0*nx)+xmin;
    
    //now y nodes
    y=CAllocD1(2*ny+1);
    for(i=0; i<2*ny+1; i++)
      y[i] = ylength*i/(2.0*ny)+ymin;

    //determine how many nodes, elements, bc's and materials there are
    *node_count   =(2*nx+1)*(2*ny+1)-nx*ny;
    *element_count=nx*ny;
    *constraint_count=2*(nx+ny);
    *force_count=0;


    //generate the node list
    xy=CAllocD2(*node_count,2);

    *node=(Node *) calloc(*node_count,sizeof(Node));
    *boundary=(Boundary *) calloc(2*(nx+ny),sizeof(Boundary));

    
    int ibc=0, ielem, inode=0;
    int *ibc2inode=CAllocI1(2*(nx+ny));
    int *ibc2ielem=CAllocI1(2*(nx+ny));
    for(j=0; j<2*ny+1; j++)
      for(i=0; i<2*nx+1; i++)
	if(!(i%2)||!(j%2)){ //only output side and corner nodes, not bubble
	  xy[inode][0]=x[i];
	  xy[inode][1]=y[j];

	  (*node)[inode].setparameters(inode+1,xy[inode]);

	  //if on a boundary
	  if((((j==0)||(j==2*ny))&&(i%2))||
	     (((i==0)||(i==2*nx))&&(j%2))) {

	    ielem=((j<2*ny)?j/2:ny-1)*nx+((i<2*nx)?i/2:nx-1);
	    ibc2ielem[ibc]=ielem;
	    ibc2inode[ibc]=inode;	    

	    (*boundary)[ibc++].setparameters(&((*node)[inode]),0.0,0.0,-3);
	  }
	  inode++;
	}


    CDeAllocD1(x);
    CDeAllocD1(y);


    //generate the element lists
    elem_node=CAllocI2(*element_count,8);
    elem_mat=CAllocI1(*element_count);
    elem_loc=CAllocI2(*element_count,2);

    ibc=0;
    ielem=0;
    if(*material_count==1) //only 1 material 
      for(j=1; j<=ny; j++)
	for(i=1; i<=nx; i++){
	  elem_node[ielem][0]=2*i-1+(j-1)*(3*nx+2);
	  elem_node[ielem][1]=2*i+1+(j-1)*(3*nx+2);
	  elem_node[ielem][2]=2*i+1+j*(3*nx+2);
	  elem_node[ielem][3]=2*i-1+j*(3*nx+2);
	  elem_node[ielem][4]=2*i+(j-1)*(3*nx+2);
	  elem_node[ielem][5]=i+1+j*(2*nx+1)+(j-1)*(nx+1);
	  elem_node[ielem][6]=2*i+j*(3*nx+2);
	  elem_node[ielem][7]=i+j*(2*nx+1)+(j-1)*(nx+1);
	  elem_mat[ielem]=1; //only 1 material so don't need to get from map
	  elem_loc[ielem  ][0]=i-1;
	  elem_loc[ielem++][1]=j-1;}
    else //more than 1 material
      for(j=1; j<=ny; j++)
	for(i=1; i<=nx; i++){
	  elem_node[ielem][0]=2*i-1+(j-1)*(3*nx+2);
	  elem_node[ielem][1]=2*i+1+(j-1)*(3*nx+2);
	  elem_node[ielem][2]=2*i+1+j*(3*nx+2);
	  elem_node[ielem][3]=2*i-1+j*(3*nx+2);
	  elem_node[ielem][4]=2*i+(j-1)*(3*nx+2);
	  elem_node[ielem][5]=i+1+j*(2*nx+1)+(j-1)*(nx+1);
	  elem_node[ielem][6]=2*i+j*(3*nx+2);
	  elem_node[ielem][7]=i+j*(2*nx+1)+(j-1)*(nx+1);
	  Get_raster_id(matres,xy[ielem][0],xy[ielem][1],
			&(elem_mat[ielem]));
	  elem_loc[ielem  ][0]=i-1;
	  elem_loc[ielem++][1]=j-1;}

    /******************************************/
    /* pass the grid directly to preprocess.C */
    /******************************************/

    int w; //name is legacy, it might stand for "which" or "while"
    Node* address[8];


    *element=(Element *) calloc(*element_count,sizeof(Element));

    //store the elements
    for(ielem=0;ielem<*element_count;ielem++){
      for(j=0; j<8; j++){
	w = elem_node[ielem][j] - 1;
	if((*node)[w].get_nodeid()!=elem_node[ielem][j]) {
	  w = 0;
	  while((*node)[w].get_nodeid()!=elem_node[ielem][j]) w++;
	  printf("found the node the long way\n");fflush(stdout);}
	address[j]=&((*node)[w]);}
      (*element)[ielem].setparameters(ielem+1,address,
					 elem_mat[ielem]-1,
					 elem_loc[ielem]);}

    //store the boundary conditions (bc's)
    for(ibc=0;ibc<2*(nx+ny);ibc++)
      (*element)[ibc2ielem[ibc]].set_boundary(&((*boundary)[ibc]));
    

    
    //deallocate the data list
    CDeAllocD2(xy);
    CDeAllocI2(elem_node);
    CDeAllocI1(elem_mat);
    CDeAllocI2(elem_loc);
    CDeAllocI1(ibc2inode);
    CDeAllocI1(ibc2ielem);

    Delete_GIS_data();
    if(havematmap) Delete_Raster_data();}
  else{
    printf("Couldn't initialize the GIS information.\n");
    exit(1);}

  return;
}
