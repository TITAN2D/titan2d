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
 * $Id: outsum.C 136 2007-06-07 20:18:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#define ADAM_HEIGHT_FRAC 0.02

void output_summary(TimeProps* timeprops, StatProps* statprops, int savefileflag) {
  //FILE *fp=fopen("output_summary.readme","a");
  char filename[256];
  sprintf(filename,"output_summary.%06d",statprops->runid);
  FILE *fp=fopen(filename,"a");

  if(timeprops->ifstart())
    fprintf(fp,"Summary of Output from Titan2D\n\n");

  int hours, minutes; double seconds;
  timeprops->chunktime(&hours,&minutes,&seconds);

  fprintf(fp,"%d:%02d:%g (hrs:min:sec)   iter=%8d   Vave=%g [m/s]   hmax=%g [m]   Xcen=%g [m]   Ycen=%g [m]   <(X-<X>)^2>=%g [m^2]   <(Y-<Y>)^2>=%g [m^2]   area=%g [m^2]   stat volume=%g [m^3]   real volume=%g [m^3]   eroded volume=%g [m^3]   deposited volume=%g [m^3]   outflow volume=%g [m^3]   forceint=%g [m/s^2]   forcebed=%g [m/s^2]   slope=%g   savefile=%d\n",
	  hours, minutes, seconds, timeprops->iter,
	  statprops->vmean, statprops->hmax, statprops->xcen, statprops->ycen,
	  statprops->xvar, statprops->yvar,
	  statprops->area, statprops->statvolume, statprops->realvolume,
	  statprops->erodedvol, statprops->depositedvol, statprops->outflowvol,
	  statprops->forceint, statprops->forcebed,
	  statprops->slopemean,savefileflag);

  /*  fprintf(fp,"%d:%02d:%g (hrs:min:sec)   iter=%8d   V*=%g   ave slope=%g   hmax =%g [m]   area=%g [m^2]   stat volume=%g [m^3]   forceint=%g [m/s^2]   forcebed=%g [m/s^2]\n",
	  hours,minutes,seconds,timeprops->iter,statprops->vstar,
	  statprops->slopemean, statprops->hmax,  
	  statprops->area, statprops->statvolume,
	  statprops->forceint,statprops->forcebed); */

  fclose(fp);

  return;
}


void output_stoch_stats(MatProps* matprops, StatProps* statprops) {
  char filename[128];
  sprintf(filename,"statout_lhs.%02d",statprops->lhs.refnum);
  FILE *fp=fopen(filename,"w");
  fprintf(fp,"%2d %6d %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n",
	  statprops->lhs.refnum,statprops->lhs.runid,statprops->timereached,
	  statprops->vstar*matprops->Vslump,statprops->hmax,statprops->xcen,
	  statprops->ycen,statprops->xyminmax[1]-statprops->xyminmax[0],
	  statprops->xyminmax[3]-statprops->xyminmax[2]);
  fflush(fp);
  fclose(fp);
  return;
}

void output_discharge(MatProps* matprops, TimeProps* timeprops, 
		      DISCHARGE* discharge, int myid) {
  double *send, *receive, doubleswap;
  int iplane, num_planes=discharge->num_planes;

  if(num_planes>0) {

    send=CAllocD1(num_planes);
    receive=CAllocD1(num_planes);
    for(iplane=0;iplane<num_planes;iplane++)
      send[iplane]=discharge->planes[iplane][9];
    
    MPI_Reduce(send,receive,num_planes,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(myid==0) { 
      doubleswap=(matprops->LENGTH_SCALE)*(matprops->LENGTH_SCALE)*
	(matprops->HEIGHT_SCALE);
      //doubleswap=(matprops->LENGTH_SCALE)*(matprops->HEIGHT_SCALE)*
      //	sqrt((matprops->LENGTH_SCALE)*(matprops->GRAVITY_SCALE));
      
      for(iplane=0;iplane<num_planes;iplane++)
	receive[iplane]*=doubleswap;
       
      
      FILE* fp=fopen("discharge.out","a");
      if(timeprops->iter==0) {
	fprintf(fp,"number of discharge planes = %g\n",num_planes);
	for(iplane=0;iplane<num_planes;iplane++)
	  fprintf(fp,"plane %d endpoints are (%16.10g,%16.10g) (%16.10g,%16.10g)\n",iplane+1,discharge->planes[iplane][0]*(matprops->LENGTH_SCALE),discharge->planes[iplane][2]*(matprops->LENGTH_SCALE),discharge->planes[iplane][1]*(matprops->LENGTH_SCALE),discharge->planes[iplane][3]*(matprops->LENGTH_SCALE));
      }
      
      fprintf(fp,"\n%16.10g  %16.10g",timeprops->timesec(),receive[0]);
      for(iplane=1;iplane<num_planes;iplane++) 
	fprintf(fp,"  %16.10g",receive[iplane]);

      fclose(fp);
    }

    CDeAllocD1(send);
    CDeAllocD1(receive);    
  }
  
  return;
}

/*************************************/
/* adam's function: written by keith */
/*************************************/

void OUTPUT_ADAM_STATS(HashTable* El_Table, MatProps* matprops_ptr,
		       TimeProps* timeprops_ptr, StatProps* statprops_ptr) {
  int myid, numprocs;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  double velocity2=0.0;
  double vmax=0, hmax=0;
  double xy_cen[2]={0.0,0.0}, vh_cen[2]={0.0,0.0};
  double xyh_vmax[3]={0.0,0.0,0.0};
  double xyv_hmax[3]={0.0,0.0,0.0};
  double masscenterdist2=0.0, masscentermindist2=HUGE_VAL, xycen[2]={0.0,0.0};
  double vmax_min_height=
    matprops_ptr->MAX_NEGLIGIBLE_HEIGHT*512.0*ADAM_HEIGHT_FRAC;
  int i;
  struct{ //for use with MPI_MAXLOC
    double val;
    int    rank;
  } send, receive;


  xy_cen[0]=statprops_ptr->xcen/(matprops_ptr->LENGTH_SCALE);
  xy_cen[1]=statprops_ptr->ycen/(matprops_ptr->LENGTH_SCALE);

  double VxVy[2];

  HashEntryPtr* buck = El_Table->getbucketptr();
  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)) {

      HashEntryPtr entryp = *(buck+i);
      while(entryp){

	Element* EmTemp=(Element*)(entryp->value);
	if(EmTemp->get_adapted_flag()>0){
	  
	  double* state_vars = EmTemp->get_state_vars();
	  EmTemp->eval_velocity(0.0,0.0,VxVy);
	  velocity2=VxVy[0]*VxVy[0]+VxVy[1]*VxVy[1];

	  //get v and h at center of mass
	  double* xy = EmTemp->get_coord();
	  masscenterdist2=
	    (xy[0]-xy_cen[0])*(xy[0]-xy_cen[0])+
	    (xy[1]-xy_cen[1])*(xy[1]-xy_cen[1]);
	  if(masscenterdist2<masscentermindist2){
	    masscentermindist2=masscenterdist2;
	    vh_cen[0]=velocity2;
	    vh_cen[1]=state_vars[0];
	    xycen[0]=xy[0];
	    xycen[1]=xy[1];
	  }

	  //eliminate fast moving very thin pile from consideration
	  if(state_vars[0]>=vmax_min_height){

	    if(velocity2>vmax){
	      /* velocity2 is not a mistake... only need to take the root of 
		 the maximum value */
	      vmax=velocity2;

	      xyh_vmax[0]=*(EmTemp->get_coord());
	      xyh_vmax[1]=*(EmTemp->get_coord()+1);
	      xyh_vmax[2]=state_vars[0];}}

	  if(state_vars[0]>hmax) {
	    
	    hmax=state_vars[0];
	    
	    xyv_hmax[0]=*(EmTemp->get_coord());
	    xyv_hmax[1]=*(EmTemp->get_coord()+1);
	    xyv_hmax[2]=velocity2;

	  }
	}
	entryp=entryp->next;      	    
      }
    }
  /*
  printf("myid=%d xycen=%16.10g %16.10g v=%16.10g h=%16.10g\n",myid,
	 xycen[0]*matprops_ptr->LENGTH_SCALE,
	 xycen[1]*matprops_ptr->LENGTH_SCALE,
	 vh_cen[0]*sqrt(matprops_ptr->LENGTH_SCALE*(matprops_ptr->GRAVITY_SCALE)),
	 vh_cen[1]*matprops_ptr->HEIGHT_SCALE);
  */
  vh_cen[0]=sqrt(vh_cen[0]);
  vmax=sqrt(vmax);
  xyv_hmax[2]=sqrt(xyv_hmax[2]);

  /* get the max value accross all processors */


  if(numprocs>1){
    send.rank=myid;


    //at center of mass
    send.val=masscentermindist2;
    MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if(receive.rank!=0){ /* don't send location if it's already on the 
			    root processor */
      if(receive.rank==myid)
	MPI_Send(vh_cen,2,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      else if(myid==0)
	MPI_Recv(vh_cen,2,MPI_DOUBLE,receive.rank,0,MPI_COMM_WORLD,&status);}



    //at location of vmax
    send.val=vmax;
    MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    vmax=receive.val;

    if(receive.rank!=0){ /* don't send location if it's already on the 
			    root processor */
      if(receive.rank==myid)
	MPI_Send(xyh_vmax,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      else if(myid==0)
	MPI_Recv(xyh_vmax,3,MPI_DOUBLE,receive.rank,0,MPI_COMM_WORLD,&status);}


    //at location of hmax
    send.val=hmax;
    MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    hmax=receive.val;

    if(receive.rank!=0){ /* don't send location if it's already on the 
			    root processor */
      if(receive.rank==myid)
	MPI_Send(xyv_hmax,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      else if(myid==0)
	MPI_Recv(xyv_hmax,3,MPI_DOUBLE,receive.rank,0,MPI_COMM_WORLD,&status);}
  }

  if(myid == 0) {
    FILE *fp=fopen("flow_dynamics.stats","a");

    fprintf(fp,"%16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g\n",
	    timeprops_ptr->timesec(), //time in seconds

	    statprops_ptr->vmean, //average velocity

	    //x,y,v,h at center of mass
	    statprops_ptr->xcen, 
	    statprops_ptr->ycen, 
	    vh_cen[0]*sqrt(matprops_ptr->LENGTH_SCALE * 
		       (matprops_ptr->GRAVITY_SCALE)),
	    vh_cen[1]*(matprops_ptr->HEIGHT_SCALE),

	    //x,y,v,h at location of vmax
	    xyh_vmax[0]*matprops_ptr->LENGTH_SCALE,
	    xyh_vmax[1]*matprops_ptr->LENGTH_SCALE,
	    vmax*sqrt(matprops_ptr->LENGTH_SCALE * 
		       (matprops_ptr->GRAVITY_SCALE)),
	    xyh_vmax[2]*(matprops_ptr->HEIGHT_SCALE),

	    //x,y,v,h at location of hmax
	    xyv_hmax[0]*matprops_ptr->LENGTH_SCALE,
	    xyv_hmax[1]*matprops_ptr->LENGTH_SCALE,
	    xyv_hmax[2]*sqrt(matprops_ptr->LENGTH_SCALE * 
		       (matprops_ptr->GRAVITY_SCALE)),
	    hmax*(matprops_ptr->HEIGHT_SCALE));
	   

    fclose(fp);
  }
  return;
}
