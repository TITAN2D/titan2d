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
 * $Id: get_coef_and_eigen.C 143 2007-06-25 17:58:08Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
//#define DEBUGINHERE
 
#include "../header/hpfem.h"
#include "../header/geoflow.h"

double get_coef_and_eigen(HashTable* El_Table, HashTable* NodeTable,
			  MatProps* matprops_ptr, FluxProps* fluxprops_ptr,
			  TimeProps* timeprops_ptr, int ghost_flag)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  double min_distance = 1000000, max_evalue = GEOFLOW_TINY, doubleswap;

  int ibuck, ierr;
  double tiny = GEOFLOW_TINY, min_dx_dy_evalue = 10000000, hmax = 0;
  double evalue = 1.0;  //might need to change this
  //-------------------go through all the elements of the subdomain and get 
  //-------------------the coefficients and eigenvalues and calculate the time step
  double global_dt[3], dt[3]={0.0,0.0,HUGE_VAL};



#ifdef DEBUGINHERE
  FILE *fp=fopen("junk.debug","w");
#endif


  HashEntryPtr* elem_bucket_zero = El_Table->getbucketptr();
  HashEntryPtr entryp;
  int num_elem_buckets=El_Table->get_no_of_buckets();
  Element *EmTemp;


  //beginning of section that SHOULD ____NOT___ be openmp'd
  double maxinflux=0.0;  
  if((maxinflux=fluxprops_ptr->MaxInfluxNow(matprops_ptr,timeprops_ptr)*
      (matprops_ptr->epsilon))>0.0) {
    double mindx=-1.0;;

    for(ibuck=0; ibuck<num_elem_buckets; ibuck++) {
      entryp = *(elem_bucket_zero+ibuck);
      while(entryp) {
	EmTemp=(Element*)(entryp->value);
	entryp=entryp->next;
	
	if((EmTemp->get_adapted_flag()>0)||
	   (EmTemp->get_adapted_flag()<0)) {
	  mindx=((*(EmTemp->get_dx()+0)<*(EmTemp->get_dx()+1))?
		 *(EmTemp->get_dx()+0):*(EmTemp->get_dx()+1))*
	    pow(0.5,REFINE_LEVEL-EmTemp->get_gen());
	  break;
	}
      }
      if(mindx>0.0) break;
    }
  
    double dttemp=mindx/maxinflux; //equivalent to dx/eigen_speed
    
    //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
    double dttemp2=0.81*maxinflux*(matprops_ptr->GRAVITY_SCALE)/9.8;
    
    dt[2]=c_dmin1(dttemp,dttemp2); 
  } //end of section that SHOULD ____NOT___ be openmp'd



  double u_vec_alt[3];
  double* d_uvec, *dx_ptr;
  int intswap;
  double *curve, maxcurve;
  int ifanynonzeroheight=0;
  double VxVy[2];
  for(ibuck=0; ibuck<num_elem_buckets; ibuck++) {
    entryp = *(elem_bucket_zero+ibuck);
    while(entryp) {
      EmTemp=(Element*)(entryp->value);
      entryp=entryp->next;
      
      if((EmTemp->get_adapted_flag()>0)||
	 ((EmTemp->get_adapted_flag()<0)&&
	  (ghost_flag == 1))){
	//if this element does not belong on this processor don't involve!!!
	
	if(*(EmTemp->get_state_vars())>GEOFLOW_TINY) {
	  ifanynonzeroheight=1;
	  
	  /* calculate hmax */
	  if(hmax < *(EmTemp->get_state_vars()))
	    hmax = *(EmTemp->get_state_vars());
	  
	  d_uvec = EmTemp->get_d_state_vars();
	  dx_ptr = EmTemp->get_dx();
#ifdef SUNOS
	  gmfggetcoef_(EmTemp->get_state_vars(), d_uvec, 
		       (d_uvec+NUM_STATE_VARS), dx_ptr, 
		       &(matprops_ptr->bedfrict[EmTemp->get_material()]), 
		       &(matprops_ptr->intfrict), EmTemp->get_kactxy(),
		       (EmTemp->get_kactxy()+1), &tiny, 
		       &(matprops_ptr->epsilon));
	  EmTemp->calc_stop_crit(matprops_ptr);
	  intswap=EmTemp->get_stoppedflags();
#ifdef DEBUGINHERE
	  fprintf(fp,"%d\n",intswap);
#endif
	  if((intswap<0)||(intswap>2)) printf("get_coef_and_eigen stopped flag=%d\n",intswap);
	  
	  //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's 
	  //rule speed if it is smaller) because underestimating speed (which 
	  //results in over estimating the timestep) is fatal to stability...
	  VxVy[0]=(*(EmTemp->get_state_vars()+1))/
	    (*(EmTemp->get_state_vars()));
	  VxVy[1]=(*(EmTemp->get_state_vars()+2))/
	    (*(EmTemp->get_state_vars()));

	  //eigen_(EmTemp->eval_state_vars(u_vec_alt),
	  eigen_(EmTemp->get_state_vars(),
		 (EmTemp->get_eigenvxymax()),
		 (EmTemp->get_eigenvxymax()+1), &evalue, &tiny, 
		 EmTemp->get_kactxy(), EmTemp->get_gravity(),
		 VxVy);
	  
#endif
	  
	  //printf("evalue=%g\n",evalue);
	  
	  // ***********************************************************
	  // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
	  // ***********************************************************
	  doubleswap = c_dmin1(dx_ptr[0],dx_ptr[1]);
	  if(doubleswap/evalue < min_dx_dy_evalue) {
	    min_distance = doubleswap;
	    max_evalue = evalue;}
	  
	  if(evalue > 1000000000.) { 
	    curve=EmTemp->get_curvature();
	    maxcurve=(dabs(curve[0])>dabs(curve[1]))?curve[0]:curve[1];
	    
	    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
		   evalue, myid, *(EmTemp->get_state_vars()+1), 
		   *(EmTemp->get_state_vars()+2),
		   *(EmTemp->get_state_vars()),maxcurve,
		   *(EmTemp->get_coord()),*(EmTemp->get_coord()+1));
	  }
	  
	  min_dx_dy_evalue = c_dmin1(c_dmin1(dx_ptr[0],dx_ptr[1])/evalue,
				     min_dx_dy_evalue);
	}
	else{
	  EmTemp->put_stoppedflags(2);
#ifdef DEBUGINHERE
	  fprintf(fp,"%d\n",EmTemp->get_stoppedflags()); 
#endif
	}
	
      } //(EmTemp->get_adapted_flag()>0)||...
      
    }//while(entryp)
  }
  
#ifdef DEBUGINHERE
  fclose(fp);
#endif

  //if(!ifanynonzeroheight) min_dx_dy_evalue=0.0;

  dt[0] = 0.5*min_dx_dy_evalue;
  //printf("hmax=%g epsilon=%g\n",hmax,matprops_ptr->epsilon);
  dt[1] = -0.9*sqrt(hmax*(matprops_ptr->epsilon)*(matprops_ptr->GRAVITY_SCALE)/9.8); //find the negative of the max not the positive min

  //printf("myid=%d, dt={%g, %g, %g}\n",myid,dt[0],-dt[1],dt[2]);

  ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  /* if(myid==0)
    printf("myid=0, global_dt={%g, %g, %g}\n",global_dt[0],-global_dt[1],global_dt[2]);
  */

  dt[0] = 0.5*c_dmin1( global_dt[0], -global_dt[1]);
  if(dt[0]==0.0) dt[0]=0.5*global_dt[2];
  
  //printf("====== for proc %d time step should be %e from distance %e and evalue %e  =========\n",
  //	 myid, min_distance/max_evalue, min_distance, max_evalue);

  //exit(0);
  
  return dt[0];
}

