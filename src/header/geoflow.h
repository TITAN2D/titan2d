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
 * $Id: geoflow.h 152 2007-06-27 20:29:54Z dkumar $ 
 */

#ifndef __GEOFLOW
#define __GEOFLOW

/* geoflow header file */
#define WEIGHT_ADJUSTER 1


#define NUM_FREEFALLS_2_STOP 2 //stopping criteria parameter
//#define STOPCRIT_CHANGE_FLUX
//#define STOPCRIT_CHANGE_BED
//#define STOPCRIT_CHANGE_SOURCE
//#define DO_EROSION

//#define REFINE_LEVEL 3
extern int REFINE_LEVEL; //make REFINE_LEVEL a global variable that can be changed set it in  Read_grid() (datread.C) or loadrun() (restart.C)
//(mdj)2007-04-11 #define MIN_GENERATION -1 //minimum refinement level
#define MIN_GENERATION -3 //minimum refinement level

//! non member C++ function that wraps the fortran correct_() function
void correct(HashTable* NodeTable, HashTable* El_Table,
	     double dt, MatProps* matprops_ptr, 
	     FluxProps *fluxprops, TimeProps *timeprops,
	     void *EmTemp,
	     double *forceint, double *forcebed, 
	     double *eroded, double *deposited);

//! this function is legacy, the prototype exists but the function is not defined
void checknodesol(HashTable*);

//! This function assigns a global_weight to the collection of elements based on the sum of their element_weight
double element_weight(HashTable* El_Table, HashTable*, int myid, int nump);

//! This function calculates the vast majority of statistics used for output, including most of what appears in output_summary.######, the friction body forces however are not calculated in here, Keith wrote this to replace calc_volume()
void calc_stats(HashTable* El_Table, HashTable* NodeTable, int myid, 
		MatProps* matprops, TimeProps* timeprops, 
		StatProps* statprops, DISCHARGE* discharge, double d_time);

//! calc_volume() has been replaced by calc_stats(), calc_volume() is out of date legacy code, the function is still defined in step.C but it is not called.
void calc_volume(HashTable* El_Table, int myid, MatProps* matprops_ptr, 
		 TimeProps* timeprops_ptr, double d_time, double* v_star, double* nz_star);

//! get_max_momentum() is legacy, it has been replaced by calc_stats()
double get_max_momentum(HashTable* El_Table, MatProps* matprops_ptr);

//! this function prints a warning message at end of the simulation to say if the flow is still moving and thus should be run longer before using the data to make decisions
void sim_end_warning(HashTable* El_Table, MatProps* matprops_ptr, 
		     TimeProps* timeprops_ptr,double v_star);

//! this function outputs final stats for one run in a collection of stochastic/probabilistic runs
void out_final_stats(TimeProps* timeprops_ptr,StatProps* statprops_ptr);

//! this function loops through the nodes zeroing the fluxes, then loops through the elements and finds the positive x direction of the element, calculates element size, calculates local terrain elevation, slopes, and curvatures, and calculates the gravity vector in local coordinates.
void setup_geoflow(HashTable* El_Table, HashTable* NodeTable, int myid, 
		   int nump, MatProps* matprops_ptr, TimeProps *timeprops_ptr);

//! this function calculates the spatial derivatives of the state variables
void slopes(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr);

//! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model) calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum allowable timestep for this iteration.
double get_coef_and_eigen(HashTable* El_Table, HashTable* NodeTable, 
			  MatProps* matprops_ptr, FluxProps* fluxprops_ptrs,
			  TimeProps* timeprops_ptr, int ghost_flag);

//! this function transfers information during events such as ghost element data exchange and repartitioning
void move_data(int nump, int myid, HashTable* El_Table, HashTable* NodeTable, TimeProps* timeprops_ptr);

//! this function deletes the current ghost elements
void delete_ghost_elms(HashTable* El_Table, int myid) ;


//! This function loops through all the non-ghost current elements and calls the Element member function Element::calc_edge_states() which calculates the Riemann fluxes between elements and stores the Riemann fluxes in the edge nodes. 
void calc_edge_states(HashTable* El_Table, HashTable* NodeTable,
		      MatProps* matprops_ptr, TimeProps* timeprops_ptr,
		      int myid, int* order_flag, double *outflow);

//! c++ sgn function 
inline double c_sgn(double zz)
{
  double sgn;
  if(zz > GEOFLOW_TINY)
    sgn = 1;
  else if(zz < -GEOFLOW_TINY)
    sgn = -1;
  else 
    sgn = 0;

  return sgn;
}

//! c++ dmin1 function
inline double c_dmin1(double d1, double d2) {
  
  if(d1> d2)
    d1 = d2;

  return d1;
}

//! another c++ dmin1 function
inline double c_dmin1(double d1, double d2, double d3) {
  
  if(d1> d2)
    d1 = d2;
  if(d1 > d3)
    d1 = d3;

  return d1;
}

//! a c++ dmax1 function
inline double c_dmax1(double d1, double d2) {
  
  if(d1 < d2)
    d1 = d2;

  return d1;
}

//! another c++ dmax1 function
inline double c_dmax1(double d1, double d2, double d3) {
  
  if(d1 < d2)
    d1 = d2;
  if(d1 < d3)
    d1 = d3;

  return d1;
}

// a c++ dabs function 
inline double dabs(double dd) {
  if(dd < 0)
    dd = -dd;

  return dd;
}
  
/* fortran calls */
#ifdef SUNOS 
//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
extern "C" void gmfggetcoef_(double*, double*, double*, double*, double*,
			     double*, double*, double*, double*, double*);

//! the actual calculation of wave speeds (eigen vectors of the flux jacoboians) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void eigen_(double *Uvec, double *eigenvxmax, double *eigenvymax, 
		       double *evalue, double *tiny, double *kactxy, 
		       double *gravity, double *VxVy);

//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void predict_(double *Uvec, double *dUdx, double *dUdy, 
			 double *Uprev, double *tiny, double *kactxy,
			 double *dt2, double *g, double *curv,
			 double *bedfrictang, double *intfrictang, 
			 double *dgdx, double *frict_tiny, int *order_flag, 
			 double *VxVy, int *if_stopped, double *fluxcoef);

//! the actual corrector timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void correct_(double *Uvec, double *Uprev, double *fluxxp, 
			 double *fluxyp, double *fluxxm, double *fluxym, 
			 double *tiny, double *dtdx, double *dtdy, double *dt,
			 double *dUdx, double *dUdy, double *xslope, 
			 double *yslope, double *curv, double *intfrictang, 
			 double *bedfrictang, double *g, double *kactxy, 
			 double *dgdx, double *frict_tiny, double *forceint, 
			 double *forcebed, int *do_erosion, double *eroded, 
			 double *VxVy, int *if_stopped, double *fluxcoef);
#endif
#ifdef IBMSP
extern "C" void gmfggetcoef(double*, double*, double*, double*, double*, 
			    double*, double*, double*, double*, double*);
extern "C" void eigen(double *Uvec, double *eigenvxmax, double *eigenvymax, 
		      double *evalue, double *tiny, double *kactxy, 
		      double *gravity, double *VxVy);

#endif


#endif
