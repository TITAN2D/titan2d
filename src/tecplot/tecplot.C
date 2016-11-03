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
 * $Id: tecplot.C 223 2011-05-25 11:08:47Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#define TECPLOTASCII
#define TARGETPROCA -1
#define TARGETPROCB -1

#ifndef TECPLOTASCII
//generates a corrupt binary tecplot file... it doesn't work
#define num_var 8
#endif

void get_elem_orient(Element *EmTemp, int *xm, int *xp, int *ym, int *yp);
int get_elem_elev(HashTable *NodeTable, MatProps *matprops,
		  Element *EmTemp, double *elevation);


int get_ll_polygon(HashTable* El_Table, HashTable* NodeTable,
		   int myid, Element* EmArray[4]);
int get_ur_tri(HashTable* El_Table, HashTable* NodeTable, 
	       int myid, Element* EmArray[4]);


int print_bubble_node(FILE *fp, HashTable* NodeTable, MatProps* matprops,
		      Element *EmTemp);

void DumpString(FILE *fp, char *str);

void testkey(HashTable* El_Table, Element* EmTemp) {
  unsigned key[2];
  key[0]=*(EmTemp->pass_key()+0);
  key[1]=*(EmTemp->pass_key()+1);
  Element *EmTemp2=(Element*) El_Table->lookup(key);
  if(EmTemp!=EmTemp2) {
    printf("Element can't look up self with key");
    printf("\n");
  }
  if((key[0]==486561007)&&(key[1]==1010586963)) {
    printf("danger key present");
    printf("\n");
  }

  return;
}
    
void testkey2(HashTable* El_Table) {
  unsigned key[2];
  key[0]=486561007;
  key[1]=1010586963;
  /* if(((Element*) El_Table->lookup(key))==0) {
    printf("danger Element went bad");
    printf("\n");
    } */

  return;
}

void tecplotter(HashTable* El_Table, HashTable* NodeTable, 
		MatProps* matprops, TimeProps* timeprops, MapNames* mapnames,
		double v_star)
{
  int i_buck, i_neigh;   //indices
  int xp, xm, yp, ym;    //x plus, x minus, y plus, y minus directions
  int gen, *neigh_gen;   //level of grid refinement
  int numprocs;          //number of processes
  int myid, *neigh_proc; //process id's

  int done = 1; 
  unsigned key[2];
  int TECTAG = 123;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  Element *EmTemp, *EmArray[4];
  HashEntry *entryp;
  
  int e_buckets=El_Table->get_no_of_buckets();
  
  /*There are 2 ways to draw triangles
    1) as a separate ZONE of triangles 
    2) as quads with duplicate points
    I chose the 2nd option.  Why?  Because I want all the elements
    from a given processor to be in a single zone.  Also duplicate 
    points for every element that uses them have been eliminated.
    That makes the "boundary" option in tecplot work right, but 
    more importantly the nodes the triangles use aren't repeated 
    (you get the triangles for free apart from the element 
    conectivity data) and nodes for quads are reduced by about a
    factor of 4.  This reduces file the size by about a factor of 4 
    when grid adaptation is used, and by about a factor of 3 when 
    it's not.  

    MODIFICATION: 20070115: Now Each And Every BUBBLE Node of Active 
    and Ghost Elements are printed once in the order they appear in 
    the Element HashTable.  This eliminates the need to sort (quick
    or otherwise) the elements speeding up tecplot, at the cost of 
    slightly increasing the file size on multiprocessor runs 
    (adding BUBBLE nodes of GHOST elements along the right and upper
    edges of the collection of elements this processor "owns")
  */
   
  //count the number of BUBBLE nodes, and tecplot quads and tecplot 
  //triangles I need to draw
  int num_tec_node=0, num_tec_tri=0, num_tec_quad=0;
  for(i_buck=0;i_buck<e_buckets;i_buck++) {
    entryp = *(El_Table->getbucketptr()+i_buck);
    while(entryp) {	
      EmArray[0]=EmTemp=(Element*)entryp->value;
      assert(EmTemp);
      entryp=entryp->next;
      if((EmTemp->get_adapted_flag()!=TOBEDELETED)&&
	 (abs(EmTemp->get_adapted_flag())<=BUFFER)
	 ) {
	num_tec_node++;
	EmTemp->put_ithelem(num_tec_node);           

	if(EmTemp->get_adapted_flag()>TOBEDELETED) {
	  switch(get_ll_polygon(El_Table,NodeTable,myid,EmArray)) {
	  case 4: //we should draw a lower left quad
	    num_tec_quad++;
	    break;
	  case 3: //we should draw a lower left triangle
	    num_tec_tri++;
	    break;
	  case 0: //we should NOT draw a lower left polygon
	    break;
	  default: //something wierd happened report this to the "FBI"
	    assert(0);
	  }//switch(get_ll_polygon(El_Table,NodeTable,myid,EmArray))

	  //match the lower left triangle at the lower left corner of 
	  //a different subdomain with an upper right triangle in this 
	  //subdomain
	  if(get_ur_tri(El_Table,NodeTable,myid,EmArray)) num_tec_tri++;  
	  
	  neigh_proc = EmTemp->get_neigh_proc();
	  gen=EmTemp->get_gen();
	  neigh_gen=EmTemp->get_neigh_gen();
	  //one triangle is added to an element for each side whose
	  //neighboring element is more refined than current element 
	  //provided the triangle wouldn't fall on a domain boundary
	  for(i_neigh=0; i_neigh<4; i_neigh++)
	    if((neigh_gen[i_neigh]>gen)&&(neigh_proc[i_neigh]!=INIT)) 
	      num_tec_tri++;
	}//if(EmTemp->get_adapted_flag()>TOBEDELETED)
      }//if((EmTemp->get_adapted_flag()!=TOBEDELETED)&&
    }//while(entryp)
  }//for(i_buck=0;i_buck<e_buckets;i_buck++)

  int num_tec_elem=num_tec_quad+num_tec_tri;
  int num_tec_elem2=0;

  char filename[20];
  sprintf(filename,"tecpl%02d%08d.tec",myid,timeprops->iter);
  FILE *fp=fopen(filename,"w");


  //print the tecplot header
  int hrs, mins; double secs;
  timeprops->chunktime(&hrs,&mins,&secs);

  fprintf ( fp, "TITLE= \" %s: time %d:%02d:%g (hrs:min:sec), V*=%g\"\n",
	    mapnames->gis_map,hrs,mins,secs,v_star);
  fprintf ( fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"PILE_HEIGHT\","
                 "\"VOL_FRACT\", \"SOLID_X_MOMENTUM\", \"SOLID_Y_MOMENTUM\","
                 "\"FLUID_X_MOMENTUM\", \"FLUID_Y_MOMENTUM\","
                 "\"ELEVATION\", \"SOLID_SPEED\", \"FLUID_SPEED\" \n");
  
  fprintf(fp,"\n");
  fprintf(fp,"ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
	  num_tec_node, num_tec_elem);


  //print the element/BUBBLE-node information

  int num_missing_bubble_node=0; //for legacy debugging, I (Keith) 
  //believe the missing BUBBLE node problem has been solved, it 
  //doesn't seem to be a problem anymore

  for(i_buck=0;i_buck<e_buckets;i_buck++) {
    entryp = *(El_Table->getbucketptr()+i_buck);
    while(entryp) {	
      EmTemp=(Element*)entryp->value;
      assert(EmTemp);
      entryp=entryp->next;
      if((EmTemp->get_adapted_flag()!=TOBEDELETED)&&
	 (abs(EmTemp->get_adapted_flag())<=BUFFER)
	 )
	num_missing_bubble_node+=
	  print_bubble_node(fp, NodeTable, matprops, EmTemp);
    }//while(entryp)
  }//for(i_buck=0;i_buck<e_buckets;i_buck++)
  

  //print the tecplot element connectivity data
  int tec_nodes_in_tec_quad[4];
  for(i_buck=0;i_buck<e_buckets;i_buck++) {
    entryp = *(El_Table->getbucketptr()+i_buck);
    while(entryp) {	
      EmArray[0]=EmTemp=(Element*)entryp->value;
      assert(EmTemp);
      entryp=entryp->next;
      if((EmTemp->get_adapted_flag()>TOBEDELETED)&&
	 (EmTemp->get_adapted_flag()<=BUFFER)
	 ) {

	neigh_proc = EmTemp->get_neigh_proc();
	gen=EmTemp->get_gen();  
	neigh_gen=EmTemp->get_neigh_gen();
		
	switch(get_ll_polygon(El_Table,NodeTable,myid,EmArray)) {
	case 4: //we should draw a lower left quad
	  tec_nodes_in_tec_quad[0]=EmArray[2]->get_ithelem();
	  tec_nodes_in_tec_quad[1]=EmArray[0]->get_ithelem();
	  tec_nodes_in_tec_quad[2]=EmArray[1]->get_ithelem();
	  tec_nodes_in_tec_quad[3]=EmArray[3]->get_ithelem();
	  
	  fprintf(fp,"%d %d %d %d\n",
		  tec_nodes_in_tec_quad[0],
		  tec_nodes_in_tec_quad[1],
		  tec_nodes_in_tec_quad[2],
		  tec_nodes_in_tec_quad[3]);

	  num_tec_elem2++; //counter of number of drawing elements
	  break;
	case 3: //we should draw a lower left triangle
	  tec_nodes_in_tec_quad[  0]=EmArray[2]->get_ithelem();
	  tec_nodes_in_tec_quad[  1]=EmArray[0]->get_ithelem();
	  tec_nodes_in_tec_quad[  2]=
	    tec_nodes_in_tec_quad[3]=EmArray[1]->get_ithelem();

	  fprintf(fp,"%d %d %d %d\n",
		  tec_nodes_in_tec_quad[0],
		  tec_nodes_in_tec_quad[1],
		  tec_nodes_in_tec_quad[2],
		  tec_nodes_in_tec_quad[3]);

	  num_tec_elem2++; //counter of number of drawing elements
	  break;
	case 0: //we should NOT draw a lower left polygon
	  break;
	default: //something weird happened report this to the "FBI"
	  assert(0);
	}//switch(get_ll_polygon(El_Table, myid, EmArray))
	
	/**********************************************************/
	/* match the lower left triangle at the lower left corner */
	/* of a different subdomain with an upper right one       */
	/**********************************************************/

	if(get_ur_tri(El_Table,NodeTable,myid,EmArray)) {
	  tec_nodes_in_tec_quad[  0]=EmArray[1]->get_ithelem();
	  tec_nodes_in_tec_quad[  1]=EmArray[2]->get_ithelem();
	  tec_nodes_in_tec_quad[  2]=
	    tec_nodes_in_tec_quad[3]=EmArray[0]->get_ithelem();

	  fprintf(fp,"%d %d %d %d\n",
		  tec_nodes_in_tec_quad[0],
		  tec_nodes_in_tec_quad[1],
		  tec_nodes_in_tec_quad[2],
		  tec_nodes_in_tec_quad[3]);

	  num_tec_elem2++; //counter of number of drawing elements
	}//if(get_ur_tri(El_Table,myid,EmArray))
	
	/**********************************************************/
	/* add triangles for elements with more refined neighbors */
	/* (aka split neighbors) this fills in the holes          */
	/**********************************************************/
	tec_nodes_in_tec_quad[0]=EmTemp->get_ithelem();
	for(i_neigh=0; i_neigh<4; i_neigh++)  
	  if((neigh_gen[i_neigh]>gen)&&(neigh_proc[i_neigh]!=INIT)){
	    /* draw triangles as quads with the last corner of
	       the triangle repeated.  The three corners are
	       1) this element
	       2) the primary neighbor (one of sides 0->3) 
	       3) the secondary neighbor (the primary +4 side) */
	
	    EmArray[1]=EmArray[2]=EmArray[3]=NULL;

	    //The Primary Neighbor
	    EmArray[1]=(Element*) El_Table->lookup(EmTemp->get_neighbors()+i_neigh*KEYLENGTH);
	    assert(EmArray[1]);
	    tec_nodes_in_tec_quad[  1]=EmArray[1]->get_ithelem();

	    //The Secondary Neighbor
	    EmArray[2]=(Element*) El_Table->lookup(EmTemp->get_neighbors()+(i_neigh+4)*KEYLENGTH);
	    assert(EmArray[2]);
	    tec_nodes_in_tec_quad[  2]=
	      tec_nodes_in_tec_quad[3]=EmArray[2]->get_ithelem();

	    fprintf(fp,"%d %d %d %d\n",
		    tec_nodes_in_tec_quad[0],
		    tec_nodes_in_tec_quad[1],
		    tec_nodes_in_tec_quad[2],
		    tec_nodes_in_tec_quad[3]);
	    	    
	    num_tec_elem2++; //counter of number of drawing elements
	  }//if((neigh_gen[i_neigh]>gen)&&(neigh_proc[i_neigh]!=INIT))
      }//if((EmTemp->get_adapted_flag()>TOBEDELETED)&& ...
    }//while(entryp)
  }//for(i_buck=0;i_buck<e_buckets;i_buck++)
  fflush(fp);
  fclose(fp);
  //printf("num_tec_elem {1/2}={%d/%d}\n",num_tec_elem,num_tec_elem2);
  assert(num_tec_elem2==num_tec_elem);

  MPI_Barrier(MPI_COMM_WORLD);

  return;
}

  
/*
 ***********************************************************
 ***********************************************************
 *********** get orientation of an element *****************
 ***********************************************************
 ***********************************************************
 */
void get_elem_orient(Element *EmTemp, int *xm, int *xp, int *ym, int *yp){
  //this does the same thing as Andy's switch statement and is shorter
  *xp = EmTemp->get_positive_x_side();
  *xm = (2+*xp)%4;
  *yp = (1+*xp)%4;
  *ym = (3+*xp)%4;
  return;
}

/*
 ***********************************************************
 ***********************************************************
 ************* get elevation for an element's **************
 ********************** bubble node ************************
 ***********************************************************
 ***********************************************************
 */

int get_elem_elev(HashTable *NodeTable, MatProps *matprops, 
		   Element *EmTemp, double *elevation){
  int j;
  double resolution, xcoord, ycoord;
  Node *NdTemp;


  NdTemp = (Node*) NodeTable->lookup(EmTemp->pass_key());

  if(NdTemp == NULL){
    //printf("\"start\" of get_elem_elev return(1)\n"); fflush(stdout);
    j = Get_max_resolution(&resolution);
    if(j != 0) {
      printf("error in Get_max_resolution\n");
      exit(1);
    }    

    xcoord = *(EmTemp->get_coord())*(matprops->LENGTH_SCALE);
    ycoord = *(EmTemp->get_coord()+1)*(matprops->LENGTH_SCALE);
    //printf("elev_ptr=%d, ",(int) elevation);
    j = Get_elevation(resolution, xcoord, ycoord, elevation);
    //printf("elev_ptr=%d\n",(int) elevation);
    if(j != 0) {
      printf("error in Get_elevation\n");
      exit(1);
    }   
    //printf("End of get_elem_elev return(1)\n"); fflush(stdout);
    return(1);
  }
  else{
    *elevation = NdTemp->get_elevation()*(matprops->LENGTH_SCALE);
    return(0);
  }
}


/* Keith Wrote this 20061120
 ************************************************************************
 * get lower, left, and lower-left neighbors: store them in EmArray[1..3]
 * returns 4 if should draw a lower left quad
 * returns 3 if should draw a lower left triangle
 * returns 0 if shouldn't draw a lower left polygon
 ************************************************************************
*/
int get_ll_polygon(HashTable* El_Table, HashTable* NodeTable,
		   int myid, Element* EmArray[4]) {
  
  int   xm,   xp,   ym,   yp; //this element's  left right down up 
  int zmxm, zmxp, zmym, zmyp; //left/down neighbor's left right down up
  int gen=EmArray[0]->get_gen();
  int* neigh_gen=EmArray[0]->get_neigh_gen();
  int* neigh_proc=EmArray[0]->get_neigh_proc(); 

  int yada=0;


  for(int ielem=1;ielem<4; ielem++)
    EmArray[ielem]=NULL;

  //yada++;
  /*
  if((*(EmArray[0]->pass_key()+0)==905969664)&&
     (*(EmArray[0]->pass_key()+1)==0))
    printf("stop me\n");
  */
  //yada++;

  //determine which ways are left (xm) and down (ym)
  get_elem_orient(EmArray[0],&xm,&xp,&ym,&yp);

  if(((neigh_proc[xm]==INIT)||(neigh_proc[ym]==INIT)
      )||
     (((neigh_gen[xm]<gen)||(neigh_gen[ym]<gen))&&
      (EmArray[0]->get_which_son()!=0)
      )
     )
    return 0;


  //yada++;

  EmArray[1]=
    (Element*) El_Table->lookup(EmArray[0]->get_neighbors()+(xm+4)*KEYLENGTH);
  EmArray[2]=
    (Element*) El_Table->lookup(EmArray[0]->get_neighbors()+ym*KEYLENGTH);

  if(!(((neigh_proc[xm+4]==myid)||
	((neigh_proc[xm+4]==-2)&&(neigh_proc[xm]==myid)))||
       (neigh_proc[ym]==myid))
     )
    //can't get to lower left neighbor by taking either route
    //because lower and left neighbors are both ghost cells
    //so draw a lower left triangle instead of a lower left quad
    return 3;

  if((neigh_proc[xm+4]==myid)||
     ((neigh_proc[xm+4]==-2)&&(neigh_proc[xm]==myid))) {
    //can get to lower left neighbor by going left then down
    assert(EmArray[1]);
    get_elem_orient(EmArray[1],&zmxm,&zmxp,&zmym,&zmyp);
    if(*(EmArray[1]->get_neigh_proc()+zmym)==-1){
      EmArray[1]=EmArray[2]=EmArray[3]=NULL;
      return 0;
    }
    EmArray[3]=(Element*) El_Table->lookup(EmArray[1]->get_neighbors()+(zmym+4)*KEYLENGTH);
    //printf("left then down\n");
    //yada++;

  }
  else{ //neigh_proc[ym]==myid
    //printf("down then left\n");
    //can get to lower left neighbor by going down then left
    assert(EmArray[2]);
    get_elem_orient(EmArray[2],&zmxm,&zmxp,&zmym,&zmyp);
    if(*(EmArray[2]->get_neigh_proc()+zmxm)==-1){
      EmArray[1]=EmArray[2]=EmArray[3]=NULL;
      return 0;
    }
    EmArray[3]=(Element*) El_Table->lookup(EmArray[2]->get_neighbors()+zmxm*KEYLENGTH);
    //printf("down then left\n");
    //yada++;
  }
  //yada++;
  /*
  if((*(EmArray[0]->pass_key()+0)==905969664)&&
     (*(EmArray[0]->pass_key()+1)==0))
    printf("stop me\n");
  */
  if(!(EmArray[3])){
    printf("myid=%d xm=%d xp=%d ym=%d yp=%d zmxm=%d zmxp=%d zmym=%d zmyp=%d EmArray[0]->key={%u,%u}\n",myid,xm,xp,ym,yp,zmxm,zmxp,zmym,zmyp,*(EmArray[0]->pass_key()+0),*(EmArray[0]->pass_key()+1));
    scanf("%d",&yada);
    //yada++;
  }
  /*
  myid+=0;
  xm+=0;
  xp+=0;
  ym+=0;
  yp+=0;
  zmxm+=0;
  zmxp+=0;
  zmym+=0;
  zmyp+=0;  
  yada++;
  */
  assert(EmArray[3]); //make this full proof

  return 4;
}


/* Keith Wrote this 20061120
 ************************************************************************
 * get right and upper neighbors: store them in EmArray[1..2]
 * returns 3 if should draw an upper right triangle
 * returns 0 otherwise
 ************************************************************************
*/


int get_ur_tri(HashTable* El_Table, HashTable* NodeTable, 
	       int myid, Element* EmArray[4]){
  int   xm,   xp,   ym,   yp; //this  element's  left right down up 
  int xpxm, xpxp, xpym, xpyp; //right neighbor's left right down up
  int ypxm, ypxp, ypym, ypyp; //up    neighbor's left right down up
  int gen, *neigh_gen, *neigh_neigh_gen; 
  int *neigh_proc, *xneigh_neigh_proc, *yneigh_neigh_proc;

  assert(EmArray[0]);
  EmArray[1]=EmArray[2]=EmArray[3]=NULL;

  //determine which ways are left (xm) and down (ym)
  get_elem_orient(EmArray[0],&xm,&xp,&ym,&yp);
  gen=EmArray[0]->get_gen();
  neigh_gen=EmArray[0]->get_neigh_gen();
  int xpfirst=xp;
  if(neigh_gen[xp]>gen) xp+=4; //account for more refined neighbors

  neigh_proc = EmArray[0]->get_neigh_proc(); 

  //make sure it's not a domain boundary
  if((neigh_proc[xpfirst]!=INIT)&&(neigh_proc[yp]!=INIT)){

    EmArray[2]=(Element*) El_Table->lookup(EmArray[0]->get_neighbors()+yp*KEYLENGTH);
    if(!EmArray[2]) {
      ElemBackgroundCheck(El_Table,NodeTable,EmArray[0]->pass_key(),stdout);
      ElemBackgroundCheck(El_Table,NodeTable,
			  EmArray[0]->get_neighbors()+yp*KEYLENGTH,stdout);
      
    }
    assert(EmArray[2]);

    get_elem_orient(EmArray[2],&ypxm,&ypxp,&ypym,&ypyp);
    yneigh_neigh_proc=EmArray[2]->get_neigh_proc(); 

    EmArray[1]=(Element*) El_Table->lookup(EmArray[0]->get_neighbors()+xp*KEYLENGTH);
    if(!EmArray[1]) {
      ElemBackgroundCheck(El_Table,NodeTable,EmArray[0]->pass_key(),stdout);
      ElemBackgroundCheck(El_Table,NodeTable,
			  EmArray[0]->get_neighbors()+xp*KEYLENGTH,stdout);
    }
    assert(EmArray[1]);
    get_elem_orient(EmArray[1],&xpxm,&xpxp,&xpym,&xpyp);
    neigh_neigh_gen=EmArray[1]->get_neigh_gen();
    //account for more refined neighbor's neighbors
    if(neigh_neigh_gen[xpyp]>neigh_gen[xp]) xpyp+=4; 

    xneigh_neigh_proc=EmArray[1]->get_neigh_proc(); 

    //check if it's a lower left corner of another subdomain
    if((yneigh_neigh_proc[ypxp]!=INIT)&&
       (yneigh_neigh_proc[ypxp]!=neigh_proc[xp])&&
       (yneigh_neigh_proc[ypxp]!=neigh_proc[yp])&&
       (yneigh_neigh_proc[ypxp]==xneigh_neigh_proc[xpyp]))
      return(3);}

  EmArray[1]=EmArray[2]=EmArray[3]=NULL;
  
  return(0);
}


/*
 ***********************************************************
 ***********************************************************
 ************ print the element's bubble node **************
 ***********************************************************
 ***********************************************************
 */
//for use when writing ascii tecplot files
int print_bubble_node(FILE *fp, HashTable* NodeTable, MatProps* matprops,
		      Element *EmTemp){
  int num_missing_bubble_node;
  double velocity_scale, momentum_scale, elevation;

  velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities
  momentum_scale = matprops->HEIGHT_SCALE * velocity_scale; // scaling factor for the momentums

  num_missing_bubble_node=
    get_elem_elev(NodeTable,matprops,EmTemp,&elevation);

  double Vel[4];
  double volf ;
  if (*(EmTemp->get_state_vars()) > GEOFLOW_TINY)
    volf=(*(EmTemp->get_state_vars()+1)/(*(EmTemp->get_state_vars())));
  else
    volf=0;
  EmTemp->eval_velocity(0.0,0.0,Vel);
  fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e\n",
	  (*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE,
	  (*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE,
	  elevation+(*(EmTemp->get_state_vars()))*(matprops)->HEIGHT_SCALE, 
	  *(EmTemp->get_state_vars())*(matprops)->HEIGHT_SCALE,
          volf,
	  *(EmTemp->get_state_vars()+2)*momentum_scale, 
	  *(EmTemp->get_state_vars()+3)*momentum_scale, 
	  *(EmTemp->get_state_vars()+4)*momentum_scale, 
	  *(EmTemp->get_state_vars()+5)*momentum_scale, 
	  elevation, 
	  sqrt(Vel[0]*Vel[0]+Vel[1]*Vel[1])*velocity_scale,
	  sqrt(Vel[2]*Vel[2]+Vel[3]*Vel[3])*velocity_scale);

  return(num_missing_bubble_node);
}

/*
 ***********************************************************
 ***********************************************************
 ******** write a binary string the way tecplot does *******
 ***********************************************************
 ***********************************************************
 */

/* the name DumpString is what tecplot uses in 
   tecplot9/util/preplot/preplot.c
   see useful_lib.C for fwriteI() */
void DumpString(FILE *fp, char *str){ 
  int i = 0;
  while (str[i] != '\0')
    fwriteI(fp,(int)str[i++]);
  fwriteI(fp,0);
  return;
}

/*
 ***********************************************************
 ***********************************************************
 *********** geoflow visualization output ******************
 ***********************************************************
 ***********************************************************
 */

void viz_output(HashTable* El_Table, HashTable* NodeTable,
		int myid, int numprocs, MatProps* matprops, 
		TimeProps* timeprops, MapNames* mapnames)
{
  int i, k;
  int element_counter = 0;
  Element* EmTemp;
  HashEntry* entryp;
  char filename[18] = "viz_outputxxx.out";
  int e_buckets=El_Table->get_no_of_buckets();

  FILE*    fp;

  // global descriptor -- list of file names
  if(myid == 0 && timeprops->ifstart()) {

    char gl_filename[19] = "viz_filenames.out";
    fp = fopen ( gl_filename, "w" );
    fprintf(fp, "%d    !number of files\n",numprocs);
    
    for(i=0;i<numprocs;i++) {
      filename[10] = (i % 1000)/100 + 48;
      filename[11] = (i % 100)/10 + 48;
      filename[12] = (i % 10) + 48;

      fprintf( fp, "%s\n", filename);
    }
    fprintf ( fp, "\n" );
    fclose ( fp );
  }

  filename[10] = (myid % 1000)/100 + 48;
  filename[11] = (myid % 100)/10 + 48;
  filename[12] = (myid % 10) + 48;
  double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities

  // actual data output
  for(i=0; i<e_buckets; i++)
    {
      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0) 
	      element_counter++;
	  
	  entryp = entryp->next;
	}
    }
  //information below will probably be changed later
  int points = element_counter;
  int variables = 7;
  int inputs = 12;
  int parameters = 2;
 
  if(timeprops->ifstart()){
    fp = fopen ( filename, "w" );
    /*    fprintf(fp,"%s\n%s\n%s\n%s\n",mapnames->gis_main,mapnames->gis_sub,
	    mapnames->gis_mapset,mapnames->gis_map);
    fprintf(fp,"%d %g\n",matprops->material_count,matprops->intfrict);
    //material i=0 is NOT YET USED, and is reserved for something like water
    for(i=1;i<=matprops->material_count;i++)
      fprintf(fp,"%s\n%g\n",matprops->matnames[i],matprops->bedfrict[i]); */}
  else 
    fp = fopen ( filename, "a+" );
  fprintf ( fp, "%f %d %d %d %d\n", timeprops->time, parameters, points, variables, inputs );
  fprintf ( fp, "\"key unsigned[%d]\", \"coordinate double[%d]\", \"terrain_slope double[%d]\", \"pile_height double[1]\", \"pile_velocity double[%d]\", \"processor int[1]\", \"FOM double[1]\" \n",
	    2, 3, 2, 2);
  // output parameters
  fprintf(fp, "time_step int[1] %d \n", timeprops->iter);
  fprintf(fp, "number_of_processors int[1] %d \n", numprocs);

  //output field data
  for(i=0; i<e_buckets; i++)
    {
      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    {
	      double* state_vars = EmTemp->get_state_vars();
	      Node* NdTemp = (Node*) NodeTable->lookup(EmTemp->pass_key());
	      if(state_vars[0] < GEOFLOW_TINY) {
		double zero = 0;
		fprintf(fp, "%u %u %f %f %f %f %f %f %f %f %d %f \n",
			*(EmTemp->pass_key()), 
			*(EmTemp->pass_key()+1),
			(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE, 
			(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			(NdTemp->get_elevation())*(matprops)->LENGTH_SCALE,
			*(EmTemp->get_zeta()), 
			*(EmTemp->get_zeta()+1), 
			zero, zero, zero, myid, zero);
	      }
	      else {
		double Vel[4];
		EmTemp->eval_velocity(0.0,0.0,Vel);
		fprintf(fp, "%u %u %f %f %f %f %f %f %f %f %d %f \n",
			*(EmTemp->pass_key()), 
			*(EmTemp->pass_key()+1),
			(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE, 
			(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			(NdTemp->get_elevation())*(matprops)->LENGTH_SCALE,
			*(EmTemp->get_zeta()), 
			*(EmTemp->get_zeta()+1), 
			state_vars[0]*(matprops)->HEIGHT_SCALE, 
			velocity_scale*Vel[0], //*state_vars[1]/state_vars[0],
			velocity_scale*Vel[1], //state_vars[2]/state_vars[0], 
			myid, (double) (state_vars[0]/20.));  //state_vars[0]/20 is just a filler item
	      }
	    }
	  entryp = entryp->next; 
	  
	} 
    } 

  fprintf ( fp, "\n" );
  fclose ( fp );

  return;
} 
/*************************************
**************************************
**************************************
*  output mesh only
**************************************
**************************************
*************************************/
void meshplotter(HashTable* El_Table, HashTable* NodeTable, 
		 MatProps* matprops, TimeProps* timeprops, MapNames* mapnames,
		 double v_star)
{
  int myid, i;
  int numprocs;
  int material;
  int done = 1; 
  int TECTAG = 123;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int element_counter = 0;
  Element* EmTemp;
  Node* NodeTemp;
  HashEntry* entryp;
  unsigned* nodes;
  char filename[256];

  if(myid==TARGETPROCB) {printf("at meshplotter 1.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);

  sprintf(filename,"mshpl%02d%08d.tec",myid,timeprops->iter);
  int order;
  int e_buckets=El_Table->get_no_of_buckets();

  double momentum_scale = matprops->HEIGHT_SCALE * sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums

  FILE*    fp;
 
  //  if ( myid == 0 )
  // {
  fp = fopen ( filename, "w" );


  if(myid==TARGETPROCB) {printf("at meshplotter 2.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);

  int hours, minutes; double seconds;
  timeprops->chunktime(&hours,&minutes,&seconds);
  
  fprintf ( fp, "TITLE= \" %s (MESH OUTPUT) time %d:%02d:%g (hrs:min:sec), V*=%g\"\n",
	    mapnames->gis_map,hours,minutes,seconds,v_star);

  //fprintf ( fp, "TITLE= \"MESH OUTPUT\"\n" );

  fprintf ( fp, "VARIABLES = \"X\" \"Y\" \"NODAL_ELEVATION\" \"PROC\" \"PILE_HEIGHT\" \"XMOMENTUM\" \"YMOMENTUM\" \"KEY0\" \"KEY1\" \"GENERATION\" \"SON\" \"ELM_ELEVATION\" \"XSLOPE\" \"YSLOPE\" \"XCURVATURE\" \"YCURVATURE\" \"ELMLOC1\" \"ELMLOC2\"" );

  if(myid==TARGETPROCB) {printf("at meshplotter 3.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);


  for(i=0; i<e_buckets; i++)
    {
      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    element_counter++;
	  entryp = entryp->next;
	}
    }

  if(myid==TARGETPROCB) {printf("at meshplotter 4.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);

  fprintf ( fp, "\n" );
  fprintf ( fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", 
            element_counter*4, element_counter );

  int elements = El_Table->get_no_of_buckets();
  if(myid==TARGETPROCB) {printf("at meshplotter 4.1\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);

  for(i=0; i<elements; i++)
    {
      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    {
	      order = 1;
	      for ( int k = 0; k < 5; k++ )
		{
		  int help_order=*(EmTemp->get_order()+k);
		  if ( help_order > order ) order = help_order;
		}
	      nodes = EmTemp->getNode();
	      double* state_vars=EmTemp->get_state_vars();
	      double err = sqrt(*(EmTemp->get_el_error()));
	      for(int j=0; j<4; j++)
		{
		  NodeTemp = (Node*) NodeTable->lookup(nodes+j*KEYLENGTH);
		  assert(NodeTemp);
		  //int* dof = NodeTemp->getdof();
		  int jj = j;
		  if(1) {//NodeTemp->getinfo() != S_C_CON) {
		    fprintf ( fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n", 
			      (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE,
			      (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			      NodeTemp->get_elevation() * (matprops->LENGTH_SCALE), 
			      myid, state_vars[0]*(matprops)->HEIGHT_SCALE,
			      state_vars[2]*momentum_scale, state_vars[3]*momentum_scale,
			      *(EmTemp->pass_key()), 
			      *(EmTemp->pass_key()+1), 
			      EmTemp->get_gen(), 
			      EmTemp->get_which_son(), EmTemp->get_elevation() * (matprops->LENGTH_SCALE),
			      *(EmTemp->get_zeta()), *(EmTemp->get_zeta()+1),
			      *(EmTemp->get_curvature())/(matprops->LENGTH_SCALE), 
			      *(EmTemp->get_curvature()+1)/(matprops->LENGTH_SCALE),
			      *(EmTemp->get_elm_loc()), *(EmTemp->get_elm_loc()+1));
		  }
		  else { // S_C_CON will have a discontinuity in the elevation so fix that by interpolation
		    double elev;
		    int neighside, mynode;
		    if(j == EmTemp->get_which_son() + 1) {
		      mynode = j - 1;
		      neighside = j;
		    }
		    else if(j == EmTemp->get_which_son() - 1) {
		      mynode = j+1;
		      neighside = j-1;
		      if(neighside < 0)
			neighside = 3;
		    }
		    else if(EmTemp->get_which_son() == 0) {
		      mynode = 0;
		      neighside = 2;
		    }
		    else if(EmTemp->get_which_son() == 3) {
		      mynode = 3;
		      neighside = 0;
		    }
		    else {
			mynode = 1;
			neighside = 1;
		      //assert(0);
		    }
		    Node* NodeTemp2 = (Node*) NodeTable->lookup(nodes+mynode*KEYLENGTH);
		    assert(NodeTemp2);

		    elev = .5 * NodeTemp2->get_elevation();
		    Element* EmTemp2 = (Element*) El_Table->lookup((EmTemp->get_neighbors()+KEYLENGTH*neighside));
		    assert(EmTemp2);

		    NodeTemp2 = (Node*) NodeTable->lookup(EmTemp2->getNode()+j*KEYLENGTH);
		    elev += .5 * NodeTemp2->get_elevation();
		    fprintf ( fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n", 
			      (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE,
			      (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			      elev * (matprops->LENGTH_SCALE), 
			      myid, state_vars[0]*(matprops)->HEIGHT_SCALE,
			      state_vars[2]*momentum_scale, state_vars[3]*momentum_scale,
			      *(EmTemp->pass_key()), 
			      *(EmTemp->pass_key()+1), 
			      EmTemp->get_gen(), 
			      EmTemp->get_which_son(), EmTemp->get_elevation() * (matprops->LENGTH_SCALE),
			      *(EmTemp->get_zeta()), *(EmTemp->get_zeta()+1),
			      *(EmTemp->get_curvature())/(matprops->LENGTH_SCALE), 
			      *(EmTemp->get_curvature()+1)/(matprops->LENGTH_SCALE),
			      *(EmTemp->get_elm_loc()), *(EmTemp->get_elm_loc()+1));
		  }
		    
		}  
	    } 
	  
	  entryp = entryp->next; 
	  
	} 
    } 

  if(myid==TARGETPROCB) {printf("at meshplotter 5.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);


  //outData<<'\n'; 
  fprintf ( fp, "\n" );

  for(i=0; i<element_counter;i++) 
    { 
      for(int j=0; j<4; j++) 
	//outData<<i*4+j+1<<' '; 
	fprintf (fp, "%d ", i*4+j+1);

      //outData<<'\n'; 
      fprintf ( fp, "\n" );
    }
 
  if(myid==TARGETPROCB) {printf("at meshplotter 6.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);


  fclose ( fp );

  if(myid==TARGETPROCB) {printf("at meshplotter 7.0\n"); fflush(stdout);}
  MPI_Barrier(MPI_COMM_WORLD);


  //  if(myid != numprocs-1) MPI_Send(&done, 1, MPI_INT, myid+1, TECTAG, MPI_COMM_WORLD); 


  return;
} 

/*************************************
**************************************
**************************************
*  output mesh for geoflow viz people
**************************************
**************************************
*************************************/
void vizplotter(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops, TimeProps* timeprops)
{
  int myid, i;
  int numprocs;
  int material;
  int done = 1; 
  int TECTAG = 123;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int element_counter = 0;
  Element* EmTemp;
  Node* NodeTemp;
  HashEntry* entryp;
  unsigned* nodes;
  char filename[20]; //="vizplotxxxxxxxx.plt";
  sprintf(filename,"vizplot%08d.plt",timeprops->iter);

  /*filename[14] = (which % 10) + 48;
  filename[13] = (which % 100)/10 + 48;
  filename[12] = (which % 1000)/100 + 48;
  filename[11] = (which % 10000)/1000 + 48;
  filename[10] = (which % 100000)/10000 + 48;
  filename[9] = (which % 1000000)/100000 + 48;
  filename[8] = (which % 10000000)/1000000 + 48;
  filename[7] = (which % 100000000)/10000000 + 48; */
  int order;
  int e_buckets=El_Table->get_no_of_buckets();
  double momentum_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums


  FILE*    fp;
 
  if ( myid == 0 )
    {
      fp = fopen ( filename, "w" );
      fprintf ( fp, "TITLE= \"MESH OUTPUT\"\n" );
      fprintf ( fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"PROC\" \"PILE_HEIGHT\" \"XMOMENTUM\" \"YMOMENTUM\" \"KEY0\" \"KEY1\" \"GENERATION\" \"SON\"" );
    }
  else 
    {
      MPI_Recv(&done, 1, MPI_INT, myid-1, TECTAG, MPI_COMM_WORLD, &status);
      //outData.open(filename, ios::app);
      fp = fopen ( filename, "a+" );
    }


  for(i=0; i<e_buckets; i++)
    {

      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    element_counter++;
	  entryp = entryp->next;
	  
	}
    }


  //outData<<'\n';
  fprintf ( fp, "\n" );

  //outData<<"ZONE N="<<element_counter*4<<", E="<<element_counter<<", F=FEPOINT, ET=QUADRILATERAL"<<'\n';
  fprintf ( fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", element_counter*4, element_counter );

  int elements = El_Table->get_no_of_buckets();

  for(i=0; i<elements; i++)
    {

      entryp = *(El_Table->getbucketptr() + i);
      while(entryp)
	{	
	  
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    {
	      order = 1;
	      for ( int k = 0; k < 5; k++ )
		{
		  int help_order=*(EmTemp->get_order()+k);
		  if ( help_order > order ) order = help_order;
		}

	      nodes = EmTemp->getNode();
	      double* state_vars=EmTemp->get_state_vars();
	      double err = sqrt(*(EmTemp->get_el_error()));
	      for(int j=0; j<4; j++)
		{
		  NodeTemp = (Node*) NodeTable->lookup(nodes+j*KEYLENGTH);
		  //int* dof = NodeTemp->getdof();
		  fprintf ( fp, "%e %e %e %d %e %e %e %u %u %d %d\n", 
			    (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE,
			    (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE, 
			    (NodeTemp->get_elevation())*(matprops)->LENGTH_SCALE, 
			    myid, state_vars[0]*(matprops)->HEIGHT_SCALE,
			    state_vars[2]*momentum_scale, state_vars[3]*momentum_scale,
			    *(EmTemp->pass_key()), 
			    *(EmTemp->pass_key()+1), 
			    EmTemp->get_gen(), 
			    EmTemp->get_which_son());

		}  
	    } 
	  
	  entryp = entryp->next; 
	  
	} 
    } 

  //outData<<'\n'; 
  fprintf ( fp, "\n" );

  for(i=0; i<element_counter;i++) 
    { 
      for(int j=0; j<4; j++) 
	//outData<<i*4+j+1<<' '; 
	fprintf (fp, "%d ", i*4+j+1);

      //outData<<'\n'; 
      fprintf ( fp, "\n" );
    }
 
  fclose ( fp );

  if(myid != numprocs-1) MPI_Send(&done, 1, MPI_INT, myid+1, TECTAG, MPI_COMM_WORLD); 

} 

