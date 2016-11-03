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
 * $Id: triplot.C 136 2007-06-07 20:18:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
void incr_tri_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		     int myid, int numprocs, MatProps* matprops, 
		     TimeProps* timeprops, double v_star){
  int i, k;
  int element_counter = 0, element_counter2=0;
  int istri1edge, istri2edge;
  Element *EmTemp, *EmTemp2;
  HashEntry* entryp;
  //  char filename[18] = "tri_outputxxx.out";
  char filename[18] = "tri_outputxxx.bin";
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();
  Node *NodeTemp, *NodeTemp2;
  unsigned *nodes;
  FILE*    fp;

  sprintf(filename,"tri_output%03d.bin",myid);
  //filename[10] = (myid % 1000)/100 + 48;
  //filename[11] = (myid % 100)/10 + 48;
  //filename[12] = (myid % 10) + 48;
  double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities

  double momentum_scale = matprops->HEIGHT_SCALE * sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums
  double min_height = matprops->MAX_NEGLIGIBLE_HEIGHT;
  double *state_vars, *state_vars2;
  double pile_max=0.;
  double pile_min=1.E30;
  double xmom_max=0.;
  double xmom_min=1.E30;
  double ymom_max=0.;
  double ymom_min=1.E30;
  double v_max=0.0;
  double velocity2;
  double snd_array[4];
  double rcv_array[4];
  // actual data output

  //printf("alive on %d\n",myid);
  for(i=0; i<e_buckets; i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp){

	EmTemp = (Element*)entryp->value;
	assert(EmTemp);
	state_vars=EmTemp->get_state_vars();
	
	if(EmTemp->get_adapted_flag()>0)	  
	  if (timeprops->ifstart() || 
	      (state_vars[0]*(matprops)->HEIGHT_SCALE>GEOFLOW_TINY)){
	    
	    element_counter++;
	    
	    
	    if (pile_max < (state_vars[0]*(matprops)->HEIGHT_SCALE))
	      pile_max = (state_vars[0]*(matprops)->HEIGHT_SCALE);
	    
	    if (pile_min > (state_vars[0]*(matprops)->HEIGHT_SCALE))
	      pile_min = (state_vars[0]*(matprops)->HEIGHT_SCALE);
	    
	    if(state_vars[0]>min_height){
	      if (xmom_max < (state_vars[1]*momentum_scale))
		xmom_max = (state_vars[1]*momentum_scale);
	      
	      if (xmom_min > (state_vars[1]*momentum_scale))
		xmom_min = (state_vars[1]*momentum_scale);
	      
	      if (ymom_max < (state_vars[2]*momentum_scale))
		ymom_max = (state_vars[2]*momentum_scale);
	      
	      if (ymom_min > (state_vars[2]*momentum_scale))
		ymom_min = (state_vars[2]*momentum_scale);
	      
	      double VxVy[2];
	      EmTemp->eval_velocity(0.0,0.0,VxVy);
	      velocity2 = VxVy[0]*VxVy[0]+VxVy[1]*VxVy[1];
	      if(v_max < velocity2) v_max=velocity2;}
	    
	  }
	
	entryp = entryp->next;
	
      }
    }
  snd_array[0]=pile_max; snd_array[1]=xmom_max; snd_array[2]=ymom_max;
  snd_array[3]=v_max;
  int ier=MPI_Allreduce (snd_array,rcv_array,4, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  pile_max=rcv_array[0]; xmom_max=rcv_array[1]; ymom_max=rcv_array[2];

  v_max=sqrt(rcv_array[3])*velocity_scale;

  snd_array[0]=pile_min; snd_array[1]=xmom_min; snd_array[2]=ymom_min;
  ier=MPI_Allreduce (snd_array,rcv_array,3, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  pile_min=rcv_array[0]; xmom_min=rcv_array[1]; ymom_min=rcv_array[2];
  //printf("done comm %d\n",myid);

  if(timeprops->ifstart())
    fp = fopen_bin(filename,"w");
  else 
    fp = fopen_bin(filename,"a");
  

  /* file version: it's the date the file format was established, this will
     allow old and new versions of the tri_outputxxx.bin file to be read
     by the same code (with if statements) */
  if(timeprops->ifstart()) fwriteI(fp,20030721); 
  
  //printf("ready to write data  for %d elements \n",element_counter);
  double xmin,xmax,ymin,ymax,resolution,elevmax,elevmin;
  i=Get_window(&xmin,&xmax,&ymin,&ymax);
  i=Get_max_resolution(&resolution);
  i = Get_elev_max(resolution, &elevmax);
  i = Get_elev_min(resolution, &elevmin);

  /* fprintf ( fp, "No. of Triangles=%d, generated from %d, elements, time step %d, time %e, Pile Max = %e, Pile Min = %e, XMom mx= %e, XMom Min= %e YMom Max= %e YMom Min= %e Xmax= %e, Xmin =%e, Ymax= %e, Ymin= %e,Zmax =%e, Zmin= %e\n", element_counter*2, element_counter,time_step,time, pile_max,pile_min,xmom_max, xmom_min,ymom_max,ymom_min,xmax,xmin,ymax,ymin,elevmax,elevmin );
     fprintf ( fp, "\n" ); */
  
     
  //write this time step's header
  fwriteI(fp,element_counter*2); //number of triangles
  fwriteI(fp,timeprops->iter);   //the time step indice

  fwriteF(fp,timeprops->time*timeprops->TIME_SCALE); //the time at this time step
  fwriteF(fp,pile_max); 
  fwriteF(fp,pile_min);
  fwriteF(fp,xmom_max);
  fwriteF(fp,xmom_min);
  fwriteF(fp,ymom_max);
  fwriteF(fp,ymom_min);
  fwriteF(fp,xmax);
  fwriteF(fp,xmin);
  fwriteF(fp,ymax);
  fwriteF(fp,ymin);
  fwriteF(fp,elevmax);
  fwriteF(fp,elevmin);
  fwriteF(fp,v_max);


  int elements = HT_Elem_Ptr->get_no_of_buckets();

  for(i=0; i<elements; i++){

    entryp = *(HT_Elem_Ptr->getbucketptr() + i);
    while(entryp){

      EmTemp = (Element*)entryp->value;
      assert(EmTemp);
      if(EmTemp->get_adapted_flag()>0){

	nodes = EmTemp->getNode();
	double* state_vars=EmTemp->get_state_vars();
	double err = sqrt(*(EmTemp->get_el_error()));
	
	// declare some data Abani for generating the triangles

	double xco[5],yco[5],z_elev[5],pile_ht,xmom,ymom;
	//
	for(int j=0; j<4; j++){
	  NodeTemp = (Node*) HT_Node_Ptr->lookup(nodes+j*KEYLENGTH);
	  //int* dof = NodeTemp->getdof();
	  if(NodeTemp->getinfo() != S_C_CON) {
	    
	    xco[j]=   (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE;
	    yco[j]=  (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
	    z_elev[j]=NodeTemp->get_elevation() * (matprops->LENGTH_SCALE);}
	  else { // S_C_CON will have a discontinuity in the elevation so fix that by interpolation
	    double elev;
	    int neighside, mynode;
	    if(j == EmTemp->get_which_son() + 1) {
	      mynode = j - 1;
	      neighside = j;}
	    else if(j == EmTemp->get_which_son() - 1) {
	      mynode = j+1;
	      neighside = j-1;
	      if(neighside < 0)
		neighside = 3;}
	    else if(EmTemp->get_which_son() == 0) {
	      mynode = 0;
	      neighside = 2;}
	    else if(EmTemp->get_which_son() == 3) {
	      mynode = 3;
	      neighside = 0;}
	    else {
	      mynode = 1;
	      neighside = 1;
	      //assert(0);
	    }
	    NodeTemp2 = (Node*) HT_Node_Ptr->lookup(nodes+mynode*KEYLENGTH);
	    elev = .5 * NodeTemp2->get_elevation();
	    EmTemp2 = (Element*) HT_Elem_Ptr->lookup((EmTemp->get_neighbors()+KEYLENGTH*neighside));
	    NodeTemp2 = (Node*) HT_Node_Ptr->lookup(EmTemp2->getNode()+j*KEYLENGTH);
	    elev += .5 * NodeTemp2->get_elevation();
	    
	    xco[j]=   (*(NodeTemp->get_coord()))*(matprops)->LENGTH_SCALE;
	    yco[j]=  (*(NodeTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
	    z_elev[j]=       elev * (matprops->LENGTH_SCALE) ;


	  }
	}//end of nodes loop
	if (timeprops->ifstart() || state_vars[0]*(matprops)->HEIGHT_SCALE>GEOFLOW_TINY ){
	  element_counter2++;

	  //check both triangles for a neighbor with "no" height
	  
	  //first the triangle with sides 0 and 1 (points 0,1,2)
	  istri1edge=0;  //default is it's not an edge
	  EmTemp2=(Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+0*KEYLENGTH); //side 0
	  if(EmTemp2==NULL) istri1edge=1;
	  else if(*(EmTemp2->get_state_vars())*(matprops)->HEIGHT_SCALE<GEOFLOW_TINY)
	    istri1edge=1;
	  else{
	    EmTemp2=(Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+1*KEYLENGTH); //side 1
	    if(EmTemp2==NULL) istri1edge=1;
	    else if(*(EmTemp2->get_state_vars())*(matprops)->HEIGHT_SCALE<GEOFLOW_TINY)
	      istri1edge=1;}
	  
	  //next the triangle with sides 2 and 3 (points 0,2,3)
	  istri2edge=0;  //default is it's not an edge
	  EmTemp2=(Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+2*KEYLENGTH); //side 2
	  if(EmTemp2==NULL) istri2edge=1;
	  else if(*(EmTemp2->get_state_vars())*(matprops)->HEIGHT_SCALE<GEOFLOW_TINY)
	    istri2edge=1;
	  else{
	    EmTemp2=(Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+3*KEYLENGTH); //side 3
	    if(EmTemp2==NULL) istri2edge=1;
	    else if(*(EmTemp2->get_state_vars())*(matprops)->HEIGHT_SCALE<GEOFLOW_TINY)
	      istri2edge=1;}

		      
	
	  /* fprintf(fp, "%u %u %d %e %e %e %e %e %e %e %e %e %e %e %e %d \n",
	   *(EmTemp->pass_key()),*(EmTemp->pass_key()+1),EmTemp->get_gen(),xco[0],yco[0],z_elev[0],xco[1],yco[1],z_elev[1],xco[2],yco[2],z_elev[2],state_vars[0]*(matprops)->HEIGHT_SCALE,state_vars[1]*momentum_scale,state_vars[2]*momentum_scale,istri1edge); */

	  
	  //write the first triangle's data
	  fwriteI(fp,*(EmTemp->pass_key()));
	  fwriteI(fp,*(EmTemp->pass_key()+1));
	  fwriteI(fp,EmTemp->get_gen());
	  fwriteF(fp,xco[0]);
	  fwriteF(fp,yco[0]);
	  fwriteF(fp,z_elev[0]);
	  fwriteF(fp,xco[1]);
	  fwriteF(fp,yco[1]);
	  fwriteF(fp,z_elev[1]);
	  fwriteF(fp,xco[2]);
	  fwriteF(fp,yco[2]);
	  fwriteF(fp,z_elev[2]);
	  fwriteF(fp,state_vars[0]*(matprops)->HEIGHT_SCALE);
	  fwriteF(fp,state_vars[1]*momentum_scale);
	  fwriteF(fp,state_vars[2]*momentum_scale);
	  fwriteI(fp,istri1edge);
	  
	  
	  /* fprintf(fp, "%u %u %d %e %e %e %e %e %e %e %e %e %e %e %e %d \n",
	   *(EmTemp->pass_key()),*(EmTemp->pass_key()+1),EmTemp->get_gen(),xco[2],yco[2],z_elev[2],xco[3],yco[3],z_elev[3],xco[0],yco[0],z_elev[0],state_vars[0]*(matprops)->HEIGHT_SCALE,state_vars[1]*momentum_scale,state_vars[2]*momentum_scale,istri2edge); 
	   */
	   
	  //write the second triangle's data
	  fwriteI(fp,*(EmTemp->pass_key()));
	  fwriteI(fp,*(EmTemp->pass_key()+1));
	  fwriteI(fp,EmTemp->get_gen());
	  fwriteF(fp,xco[2]);
	  fwriteF(fp,yco[2]);
	  fwriteF(fp,z_elev[2]);
	  fwriteF(fp,xco[3]);
	  fwriteF(fp,yco[3]);
	  fwriteF(fp,z_elev[3]);
	  fwriteF(fp,xco[0]);
	  fwriteF(fp,yco[0]);
	  fwriteF(fp,z_elev[0]);
	  fwriteF(fp,state_vars[0]*(matprops)->HEIGHT_SCALE);
	  fwriteF(fp,state_vars[1]*momentum_scale);
	  fwriteF(fp,state_vars[2]*momentum_scale);
	  fwriteI(fp,istri2edge);
	  
	  
	}
      }
      entryp = entryp->next;
      
    }
  }
  
  //printf("element_counter(1,2)=%d,%d\n",element_counter,element_counter2);
  //outData<<'\n';
  //fprintf ( fp, "\n" );



  fclose ( fp );
  return;
}


