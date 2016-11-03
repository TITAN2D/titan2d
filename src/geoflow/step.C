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
 * $Id: step.C 164 2007-07-03 15:27:22Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#define APPLY_BC

void step(HashTable* El_Table, HashTable* NodeTable, int myid, int nump,
	  MatProps* matprops_ptr, TimeProps* timeprops_ptr, 
	  PileProps *pileprops_ptr, FluxProps *fluxprops, 
	  StatProps* statprops_ptr, int* order_flag, 
	  OutLine* outline_ptr, DISCHARGE* discharge, int adaptflag)
{
  /* 
   * PREDICTOR-CORRECTED based on Davis' Simplified Godunov Method 
   */

  /* pass off proc data here (really only need state_vars for off-proc neighbors) */
  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);

  slopes(El_Table, NodeTable, matprops_ptr);

  // get coefficients, eigenvalues, hmax and calculate the time step 
  double dt = get_coef_and_eigen(El_Table, NodeTable, matprops_ptr, 
                                 fluxprops, timeprops_ptr,0);

  timeprops_ptr->incrtime(&dt); //also reduces dt if necessary
  
  // assign influxes and then if any new sources are activating in 
  // current time step refine and re-mark cells 
  adapt_fluxsrc_region(El_Table,NodeTable,matprops_ptr,pileprops_ptr,fluxprops,
                       timeprops_ptr,dt,myid,adaptflag);

  int i;
  HashEntryPtr* buck = El_Table->getbucketptr();
  HashEntryPtr currentPtr;
  Element* Curr_El;

  double dt2 = .5*dt; // dt2 is set as dt/2 !

  /*
   *  predictor step
   */
  int j,k, counter;
  double tiny = GEOFLOW_TINY;
  double flux_src_coef=0;

#ifdef SECOND_ORDER
  //-------------------go through all the elements of the subdomain and  
  //-------------------calculate the state variables at time .5*delta_t
  /* mdj 2007-04 */
  int IF_STOPPED;
  double curr_time,influx[3],*d_uvec; //VxVy[2];
  double VxVy[2];
  Node* nd;
#pragma omp parallel for                                                \
private(currentPtr,Curr_El,IF_STOPPED,influx,j,k,curr_time,flux_src_coef,VxVy)
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
    {
      currentPtr = *(buck+i);
      while(currentPtr){

	Curr_El=(Element*)(currentPtr->value);

	  influx[3];
	  influx[0]=*(Curr_El->get_influx()+0);
	  influx[1]=*(Curr_El->get_influx()+1);
	  influx[2]=*(Curr_El->get_influx()+2);

	  if(!(influx[0]>=0.0)){
	    printf("negative influx=%g\n",influx[0]);
	    assert(0);
	  }


	if(Curr_El->get_adapted_flag()>0){

	  d_uvec = Curr_El->get_d_state_vars();
	  nd = (Node*) NodeTable->lookup(Curr_El->pass_key());

  // -- calc contribution of flux source
          flux_src_coef=0;
          curr_time=(timeprops_ptr->time)*(timeprops_ptr->TIME_SCALE);

	  //VxVy[2]; 
	  if(*(Curr_El->get_state_vars()+0)>GEOFLOW_TINY) {
	    VxVy[0]=*(Curr_El->get_state_vars()+1)/ *(Curr_El->get_state_vars()+0);
	    VxVy[1]=*(Curr_El->get_state_vars()+2)/ *(Curr_El->get_state_vars()+0);
	  }
	  else
	    VxVy[0]=VxVy[1]=0.0;

#ifdef STOPCRIT_CHANGE_SOURCE
	  IF_STOPPED=Curr_El->get_stoppedflags();
#else
	  IF_STOPPED=!(!(Curr_El->get_stoppedflags()));
#endif

	  predict_(Curr_El->get_state_vars(), d_uvec, (d_uvec+NUM_STATE_VARS),
		   Curr_El->get_prev_state_vars(), &tiny, 
		   Curr_El->get_kactxy(), &dt2, Curr_El->get_gravity(), 
		   Curr_El->get_curvature(),
		   &(matprops_ptr->bedfrict[Curr_El->get_material()]), 
		   &(matprops_ptr->intfrict),
		   Curr_El->get_d_gravity(), &(matprops_ptr->frict_tiny), 
		   order_flag, VxVy, 
		   &IF_STOPPED,influx);

	  /* apply bc's */
#ifdef APPLY_BC
	  for(j=0;j<4;j++)
	    if(*(Curr_El->get_neigh_proc()+j) == INIT)   // this is a boundary!
	      for(k=0;k<NUM_STATE_VARS;k++)
		*(Curr_El->get_state_vars()+k) = 0;
#endif
	}
	currentPtr=currentPtr->next;      	    
      }
    }
  /* finished predictor step */
  /* really only need to share dudx, state_vars, and kactxy */
  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);

  /* calculate the slopes for the new (half-time step) state variables */
  slopes(El_Table, NodeTable, matprops_ptr);
#endif  //SECOND_ORDER
  
  /* really only need to share dudx, state_vars, and kactxy */
  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);

  /* calculate kact/pass */
  double dt_not_used = get_coef_and_eigen(El_Table, NodeTable, matprops_ptr,
                       fluxprops, timeprops_ptr, 1);

  /*
   * calculate edge states
   */
  double outflow=0.0;  //shouldn't need the =0.0 assignment but just being cautious.
  //printf("step: before calc_edge_states\n"); fflush(stdout);
  calc_edge_states(El_Table,NodeTable,matprops_ptr,timeprops_ptr,myid,
		   order_flag,&outflow);
  outflow*=dt;

  /*
    * corrector step and b.c.s
    */

  //for comparison of magnitudes of forces in slumping piles
  double forceint=0.0, elemforceint;
  double forcebed=0.0, elemforcebed;
  double eroded=0.0, elemeroded;
  double deposited=0.0, elemdeposited;
  double realvolume=0.0;

  buck = El_Table->getbucketptr();
  // mdj 2007-04 this loop has pretty much defeated me - there is
  //             a dependency in the Element class that causes incorrect
  //             results
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	 {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_adapted_flag()>0) { //if this is a refined element don't involve!!!
	      
	      double *dxy=Curr_El->get_dx();
              // if calculations are first-order, predict is never called
              // ... so we need to update prev_states
              if ( *order_flag == 1 )
                Curr_El->update_prev_state_vars();

	      void *Curr_El_out= (void *) Curr_El;
	      correct(NodeTable, El_Table, dt, matprops_ptr,
		      fluxprops, timeprops_ptr,
		      Curr_El_out,
		      &elemforceint,&elemforcebed,
		      &elemeroded,&elemdeposited);

	      forceint+=fabs(elemforceint);
	      forcebed+=fabs(elemforcebed);
	      realvolume+=dxy[0]*dxy[1]**(Curr_El->get_state_vars());
	      eroded+=elemeroded;
	      deposited+=elemdeposited;

	      double *coord=Curr_El->get_coord();	      
	      //update the record of maximum pileheight in the area covered by this element
	      double hheight=*(Curr_El->get_state_vars());
	      if(hheight>0 && hheight<0);
	      double pfheight[6];


#ifdef MAX_DEPTH_MAP
	      outline_ptr->update(coord[0]-0.5*dxy[0],coord[0]+0.5*dxy[0],
				  coord[1]-0.5*dxy[1],coord[1]+0.5*dxy[1],
				  hheight,pfheight);      
#endif

#ifdef APPLY_BC
	      for(j=0;j<4;j++)
		if(*(Curr_El->get_neigh_proc()+j) == INIT)   // this is a boundary!
		  for(k=0;k<NUM_STATE_VARS;k++)
		    *(Curr_El->get_state_vars()+k) = 0;
#endif
	    }
	    currentPtr=currentPtr->next;      	    
	  }
      }


  //update the orientation of the "dryline" (divides partially wetted cells
  //into wet and dry parts solely based on which neighbors currently have 
  //pileheight greater than GEOFLOW_TINY
  for(i=0; i<El_Table->get_no_of_buckets(); i++) {
    HashEntryPtr currentPtr = *(buck+i);

    while(currentPtr) {      
      Element* Curr_El=(Element*)(currentPtr->value);
      currentPtr=currentPtr->next;      	    
      if(Curr_El->get_adapted_flag()>0) //if this is a refined element don't involve!!!
	Curr_El->calc_wet_dry_orient(El_Table);
    }
  }

  /* finished corrector step */

  //printf("EXIT at %s: %d",__FILE__, __LINE__);
  //exit(1);

  calc_stats(El_Table, NodeTable, myid, matprops_ptr, timeprops_ptr, 
	     statprops_ptr, discharge, dt);

  //  if(statprops_ptr->timereached!=-1.0)
  //    printf("timereached=%g after calc_stats()\n",statprops_ptr->timereached);


  double tempin[6], tempout[6];
  tempin[0]=outflow;    //volume that flew out the boundaries this iteration
  tempin[1]=eroded;     //volume that was eroded this iteration
  tempin[2]=deposited;  //volume that is currently deposited
  tempin[3]=realvolume; //"actual" volume within boundaries
  tempin[4]=forceint;   //internal friction force
  tempin[5]=forcebed;   //bed friction force

  MPI_Reduce(tempin,tempout,6,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  statprops_ptr->outflowvol+=tempout[0]*(matprops_ptr->HEIGHT_SCALE)*
    (matprops_ptr->LENGTH_SCALE)*(matprops_ptr->LENGTH_SCALE);
  statprops_ptr->erodedvol+=tempout[1]*(matprops_ptr->HEIGHT_SCALE)*
    (matprops_ptr->LENGTH_SCALE)*(matprops_ptr->LENGTH_SCALE);
  statprops_ptr->depositedvol=tempout[2]*(matprops_ptr->HEIGHT_SCALE)*
    (matprops_ptr->LENGTH_SCALE)*(matprops_ptr->LENGTH_SCALE);
  statprops_ptr->realvolume=tempout[3]*(matprops_ptr->HEIGHT_SCALE)*
    (matprops_ptr->LENGTH_SCALE)*(matprops_ptr->LENGTH_SCALE);

  statprops_ptr->forceint=tempout[4]/tempout[3]*matprops_ptr->GRAVITY_SCALE;
  statprops_ptr->forcebed=tempout[5]/tempout[3]*matprops_ptr->GRAVITY_SCALE;

  //calc_volume(El_Table, myid, matprops_ptr, timeprops_ptr, dt, v_star, nz_star);

  return;
}


/***********************************************************************/
/* calc_volume():                                                      */
/* calculates volume to verify mass conservation                       */
/* determines the maximum height                                       */
/* determines the maximum velocity                                     */
/* determines and returns v_star the non-dimensional stopping velocity */
/***********************************************************************/

void calc_volume(HashTable* El_Table, int myid, MatProps* matprops_ptr, 
		 TimeProps* timeprops_ptr, double d_time, double *v_star, 
		 double *nz_star)
{
  int i,j,k, counter, imax = 0;
  double tiny = GEOFLOW_TINY;
  //-------------------go through all the elements of the subdomain and  
  //-------------------calculate the state variables at time .5*delta_t
  double volume = 0, volume2 =0, max_height = 0;
  double v_ave=0, gl_v_ave;
  double g_ave=0, gl_g_ave;
  double v_max=0, gl_v_max;
  double min_height=matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
  register double temp;

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
      {
	HashEntryPtr currentPtr = *(buck+i);
	while(currentPtr)
	  {
	    Element* Curr_El=(Element*)(currentPtr->value);
 	    if(Curr_El->get_adapted_flag()>0){
	      
	      double* state_vars = Curr_El->get_state_vars();
	      double dx = *(Curr_El->get_dx());
	      double dy = *(Curr_El->get_dx()+1);

	      if(state_vars[0] > max_height) 
              {
		max_height = state_vars[0];
		imax = i;
              }

	      // rule out non physical fast moving thin layers
	      if(state_vars[0] > min_height)
              {
		temp=sqrt(state_vars[2]*state_vars[2]
			  +state_vars[3]*state_vars[3]);
		v_ave+=temp*dx*dy;
		temp/=state_vars[1];  

		double dvol=state_vars[0]*dx*dy;
		g_ave+=*(Curr_El->get_gravity()+2)*dvol;
		volume2+=dvol;

		if(temp>v_max) v_max=temp;
              }
	      volume += state_vars[0]*dx*dy;
	    }
	    currentPtr=currentPtr->next;      	    
	  }
      } 
  
  double gl_volume = 0, gl_volume2=0, gl_max_height;
  double send[5], receive[5]={0.0,0.0,0.0,0.0,0.0};
  send[0]=volume;
  send[1]=volume2;
  send[2]=v_ave;
  send[3]=g_ave;
  i = MPI_Reduce(send, receive, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gl_volume =receive[0];
  gl_volume2=receive[1];
  gl_v_ave  =receive[2];
  *nz_star  =receive[3];
  
  send[0]=max_height;
  send[1]=v_max;

  i = MPI_Reduce(send, receive, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  gl_max_height=receive[0];
  gl_v_max     =receive[1];


  if(myid == 0) {
    *nz_star=*nz_star/gl_volume2*matprops_ptr->GRAVITY_SCALE/9.8;
    //dimensionalize
    gl_v_ave=gl_v_ave/gl_volume2*sqrt(matprops_ptr->LENGTH_SCALE * 
				      matprops_ptr->GRAVITY_SCALE);
    gl_v_max=gl_v_max*sqrt(matprops_ptr->LENGTH_SCALE * 
			   (matprops_ptr->GRAVITY_SCALE));

    gl_volume = gl_volume*(matprops_ptr->LENGTH_SCALE)*
                (matprops_ptr->LENGTH_SCALE)*(matprops_ptr->HEIGHT_SCALE);
    gl_max_height = gl_max_height * (matprops_ptr->HEIGHT_SCALE);


    d_time*=timeprops_ptr->TIME_SCALE;

    /* v_star is the nondimensional global average velocity by v_slump
       once v_slump HAS BEEN CALIBRATED (not yet done see ../main/datread.C) 
       the calculation will terminate when v_star reaches 1 */
    *v_star=gl_v_ave/matprops_ptr->Vslump;

    //chunk time
    int hours, minutes; double seconds;
    timeprops_ptr->chunktime(&hours, &minutes, &seconds);

    printf("At the end of time step %d the time is %d:%02d:%g (hrs:min:sec),\n"
           "time step length is %g [sec], volume is %g [m^3],\n"
           "max height is %g [m], max velocity is %g [m/s],\n"
           "ave velocity is %g [m/s], v* = %g\n\n",
	   timeprops_ptr->iter, hours, minutes, seconds, d_time,
	   gl_volume, gl_max_height, gl_v_max,gl_v_ave,*v_star);
  }

  return;
}

/***********************************************************************/
/* the get_max_momentum function was put here because it is similar to */
/* calc volume... which is admittedly not the best reason so if you    */
/* can think of a better place to put it go ahead                      */
/***********************************************************************/

double get_max_momentum(HashTable* El_Table, MatProps* matprops_ptr){
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  double mom2, max_mom=0, gl_max_mom;
  double min_height=matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
  int i;

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i))
    {
      HashEntryPtr entryp = *(buck+i);
      while(entryp)
      {
	Element* EmTemp=(Element*)(entryp->value);
	if(EmTemp->get_adapted_flag()>0)
        {
	  double* state_vars = EmTemp->get_state_vars();
	  //eliminate fast moving very thin pile from consideration
	  if(state_vars[0]>=min_height)
          {
	    mom2=(state_vars[2]*state_vars[2]+state_vars[3]*state_vars[3]);
	    /* mom2 is not a mistake... only need to take the root of 
	       the maximum value */
	    if(mom2>max_mom) max_mom=mom2;
          }
	}
	entryp=entryp->next;      	    
      }
    }

  max_mom=sqrt(max_mom);

  if(numprocs>1){
    if(myid==0) printf("get_max_momentum()1: i=%d\n",i); fflush(stdout);
    i=MPI_Reduce(&max_mom, &gl_max_mom, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(myid==0) printf("get_max_momentum()2: i=%d\n",i); fflush(stdout);}
  else gl_max_mom=max_mom;

  return(gl_max_mom * matprops_ptr->HEIGHT_SCALE * 
	 sqrt(matprops_ptr->LENGTH_SCALE * (matprops_ptr->GRAVITY_SCALE))); 

}


/**********************************************************************/
/* the sim_end_warning function was put here because it is similar to */
/* calc volume... which is admittedly not the best reason so if you   */
/* can think of a better place to put it go ahead                     */
/**********************************************************************/

void sim_end_warning(HashTable* El_Table, MatProps* matprops_ptr,
		     TimeProps* timeprops_ptr, double v_star){
  FILE *fp;
  int myid, numprocs;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  /****************************************/
  /* print out the final dimensional time */
  /****************************************/
  //chunk time
  int hours, minutes; double seconds;
  timeprops_ptr->chunktime(&hours, &minutes, &seconds);
  
  if(myid==0){
    //print to screen
    printf("\nTitan2D performed %d time steps before the calculation ended.\n",
	   timeprops_ptr->iter);
    printf("%d:%02d:%g (hrs:min:sec) of time was simulated.\n",
	   hours,minutes,seconds);

    //print to file
    fp=fopen("sim_end_warning.readme","w");
    fprintf(fp,
	    "Titan2D performed %d time steps before the calculation ended.\n",
	    timeprops_ptr->iter);
    fprintf(fp,"%d:%02d:%g (hrs:min:sec) of time was simulated.\n",
	    hours,minutes,seconds);}

  /*****************************************/
  /* find and print maximum final velocity */
  /*****************************************/
  //double send[2], receive[2];
  double velocity2;
  double v_max=0;
  double xy_v_max[2];
  double min_height=matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
  int i;
  struct{ //for use with MPI_MAXLOC
    double val;
    int    rank;
  } send, receive;


  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){

      HashEntryPtr entryp = *(buck+i);
      while(entryp)
      {
	Element* EmTemp=(Element*)(entryp->value);
	if(EmTemp->get_adapted_flag()>0)
        {  
	  double* state_vars = EmTemp->get_state_vars();
	  //eliminate fast moving very thin pile from consideration
	  if(state_vars[0]>=min_height)
          {
	    velocity2=(state_vars[2]*state_vars[2]+state_vars[3]*state_vars[3])
	      /(state_vars[1]*state_vars[1]);

	    if(velocity2>v_max)
            {
	      /* velocity2 is not a mistake... only need to take the root of 
		 the maximum value */
	      v_max=velocity2;
	      xy_v_max[0]=*(EmTemp->get_coord());
	      xy_v_max[1]=*(EmTemp->get_coord()+1);}}
	}
	entryp=entryp->next;      	    
      }
    }
  v_max=sqrt(v_max);

  /* get the max value accross all processors */
  send.val=v_max;
  send.rank=myid;

  if(numprocs>1){
    MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    v_max=receive.val;

    if(receive.rank!=0){ /* don't send location if it's already on the 
			    root processor */
      if(receive.rank==myid)
	MPI_Send(xy_v_max,2,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      else if(myid==0)
	MPI_Recv(xy_v_max,2,MPI_DOUBLE,receive.rank,0,MPI_COMM_WORLD,&status);}
  }

  // print the rest of the warning
  if(myid == 0){
    //print to screen
    printf("the final v* = v/v_slump = %g\n",v_star);
    printf("The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
	   v_max*sqrt(matprops_ptr->LENGTH_SCALE * 
		      (matprops_ptr->GRAVITY_SCALE)),
	   xy_v_max[0]*matprops_ptr->LENGTH_SCALE,
	   xy_v_max[1]*matprops_ptr->LENGTH_SCALE);
    
    //print to file
    fprintf(fp,"the final v* = v/v_slump = %g\n",v_star);
    fprintf(fp,"The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
	    v_max*sqrt(matprops_ptr->LENGTH_SCALE * 
		       (matprops_ptr->GRAVITY_SCALE)),
	    xy_v_max[0]*matprops_ptr->LENGTH_SCALE,
	    xy_v_max[1]*matprops_ptr->LENGTH_SCALE);
    fclose(fp);}
  
  return;
}

