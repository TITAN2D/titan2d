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

#include "constant.h"

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
void correct(ElementType elementType,NodeHashTable* NodeTable, ElementsHashTable* El_Table, double dt, MatProps* matprops_ptr, FluxProps *fluxprops,
             TimeProps *timeprops, Integrator *integrator, void *EmTemp, double *forceint, double *forcebed, double *eroded, double *deposited);

//! this function is legacy, the prototype exists but the function is not defined
void checknodesol(NodeHashTable*);

//! calc_volume() has been replaced by calc_stats(), calc_volume() is out of date legacy code, the function is still defined in step.C but it is not called.
void calc_volume(ElementType elementType,ElementsHashTable* El_Table, int myid, MatProps* matprops_ptr, TimeProps* timeprops_ptr, double d_time,
                 double* v_star, double* nz_star);

//! get_max_momentum() is legacy, it has been replaced by calc_stats()
double get_max_momentum(ElementType elementType,ElementsHashTable* El_Table, MatProps* matprops_ptr);

//! this function prints a warning message at end of the simulation to say if the flow is still moving and thus should be run longer before using the data to make decisions
void sim_end_warning(ElementType elementType, ElementsHashTable* El_Table, MatProps* matprops_ptr, TimeProps* timeprops_ptr, double v_star);

//! this function outputs final stats for one run in a collection of stochastic/probabilistic runs
void out_final_stats(TimeProps* timeprops_ptr, StatProps* statprops_ptr);

//! this function loops through the nodes zeroing the fluxes, then loops through the elements and finds the positive x direction of the element, calculates element size, calculates local terrain elevation, slopes, and curvatures, and calculates the gravity vector in local coordinates.
void setup_geoflow(ElementsHashTable* El_Table, NodeHashTable* NodeTable, int myid, int nump, MatProps* matprops_ptr,
                   TimeProps *timeprops_ptr);

//! this function calculates the spatial derivatives of the state variables
void slopes(ElementsHashTable* El_Table, NodeHashTable* NodeTable, MatProps* matprops_ptr);


//! this function transfers information during events such as ghost element data exchange and repartitioning
void move_data(int nump, int myid, ElementsHashTable* El_Table, NodeHashTable* NodeTable, TimeProps* timeprops_ptr);

//! this function deletes the current ghost elements
void delete_ghost_elms(ElementsHashTable* El_Table, int myid);

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
//!pretty cool branchless sign
inline double signum(double zz){return (0.0<zz)-(zz<0.0);}
inline double sgn_tiny(double zz,double tiny){return (tiny<zz)-(zz<-tiny);}


//! c++ dmin1 function
inline double c_dmin1(const double d1, const double d2)
{
    if(d1 > d2)
        return d2;
    else
        return d1;
}

//! another c++ dmin1 function
inline double c_dmin1(double d1, double d2, double d3)
{
    
    if(d1 > d2)
        d1 = d2;
    if(d1 > d3)
        d1 = d3;
    
    return d1;
}

//! a c++ dmax1 function
inline double c_dmax1(const double d1, const double d2)
{
    if(d1 < d2)
        return d2;
    else
        return d1;
}

//! another c++ dmax1 function
inline double c_dmax1(double d1, double d2, double d3)
{
    
    if(d1 < d2)
        d1 = d2;
    if(d1 < d3)
        d1 = d3;
    
    return d1;
}

// a c++ dabs function 
inline double dabs(double dd)
{
    if(dd < 0)
        dd = -dd;
    
    return dd;
}

//! This function updates phi variable
//void update_phi(HashTable *El_Table, double *update);

void find_min_dx(ElementsHashTable* El_Table, double* mindx);

void reinitialization(ElementType elementType, NodeHashTable* NodeTable, ElementsHashTable* El_Table, MatProps* matprops_ptr,
    TimeProps *timeprops, PileProps *pileprops_ptr, int nump, int rank);

void initialization(ElementType elementType, NodeHashTable* NodeTable, ElementsHashTable* El_Table, MatProps* matprops_ptr,
    TimeProps *timeprops, PileProps *pileprops_ptr, int nump, int rank);

/* fortran calls */
#ifdef SUNOS 
//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
/*extern "C" void gmfggetcoef_(double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*);

//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
extern "C" void gmfggetcoef2ph_(double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*);*/

//! the actual calculation of wave speeds (eigen vectors of the flux jacoboians) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
/*extern "C" void eigen_(double *Uvec, double *eigenvxmax, double *eigenvymax,
        double *evalue, double *tiny, double *kactxy,
        double *gravity, double *VxVy);

//! the actual calculation of wave speeds (eigen vectors of the flux jacoboians) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void eigen2ph_(double *Uvec, double *eigenvxmax, double *eigenvymax,
        double *evalue, double *tiny, double *kactxy,
        double *gravity, double *Vs, double *Vf, double *eps, int *);*/

//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
/*extern "C" void predict_(double *Uvec, double *dUdx, double *dUdy,
        double *Uprev, double *tiny, double *kactxy,
        double *dt2, double *g, double *curv,
        double *bedfrictang, double *intfrictang,
        double *dgdx, double *frict_tiny, int *order_flag,
        double *VxVy, int *if_stopped, double *fluxcoef);*/

//! the actual corrector timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
/*extern "C" void correct_(double *Uvec, double *Uprev, double *fluxxp,
        double *fluxyp, double *fluxxm, double *fluxym,
        double *tiny, double *dtdx, double *dtdy, double *dt,
        double *dUdx, double *dUdy, double *xslope,
        double *yslope, double *curv, double *intfrictang,
        double *bedfrictang, double *g, double *kactxy,
        double *dgdx, double *frict_tiny, double *forceint,
        double *forcebed, int *do_erosion, double *eroded,
        double *VxVy, int *if_stopped, double *fluxcoef);*/

//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
/*extern "C" void predict2ph_(double *Uvec, double *dUdx, double *dUdy,
        double *Uprev, double *tiny, double *kactxy,
        double *dt2, double *g, double *curv,
        double *bedfrictang, double *intfrictang,
        double *dgdx, double *frict_tiny, int *order_flag,
        double *VxVy, int *if_stopped, double *fluxcoef);*/

//! the actual corrector timestep update
/*extern "C" void correct2ph_(double *Uvec, double *Uprev, double *fluxxp,
        double *fluxyp, double *fluxxm, double *fluxym,
        double *tiny, double *dtdx, double *dtdy, double *dt,
        double *dUdx, double *dUdy, double *xslope,
        double *yslope, double *curv, double *intfrictang,
        double *bedfrictang, double *g, double *kactxy,
        double *frict_tiny, double *forceint,
        double *forcebed, int *do_erosion, double *eroded,
        double *Vsolid, double *Vfluid, double *den_solid,
        double *den_fluid, double *terminal_vel, double *eps,
        int *if_stopped, double *fluxcoef);*/

#endif
#ifdef IBMSP
extern "C" void gmfggetcoef(double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*);
extern "C" void eigen(double *Uvec, double *eigenvxmax, double *eigenvymax,
        double *evalue, double *tiny, double *kactxy,
        double *gravity, double *VxVy);

#endif

#endif
