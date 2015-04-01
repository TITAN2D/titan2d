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
 * $Id: get_coef_and_eigen.C,v 1.4 2004/08/11 15:58:46 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
//#define DEBUGINHERE
//#define DO_EROSION
 
#include "../header/hpfem.h"
#include "../header/geoflow.h"
void correct(HashTable* NodeTable, HashTable* El_Table,
	     double dt, MatProps* matprops_ptr, 
	     FluxProps *fluxprops, TimeProps *timeprops,
	     void *EmTemp_in,double *forceint, double *forcebed, 
	     double *eroded, double *deposited)
{
  Element *EmTemp=(Element *) EmTemp_in;
  double *dx=EmTemp->get_dx();
  double dtdx = dt/dx[0];
  double dtdy = dt/dx[1];
  

  double tiny = GEOFLOW_TINY;
  int xp=EmTemp->get_positive_x_side();
  int yp=(xp+1)%4, xm=(xp+2)%4, ym=(xp+3)%4; 

  int ivar, j, k;

  double fluxxp[3], fluxyp[3], fluxxm[3], fluxym[3];

  Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode()+(xp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxp[ivar] = nxp->flux[ivar];
    
  Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode()+(yp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxyp[ivar] = nyp->flux[ivar];
    
  Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode()+(xm+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxm[ivar] = nxm->flux[ivar];
    
  Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode()+(ym+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxym[ivar] = nym->flux[ivar];
    

#ifdef SUNOS
  /* the values being passed to correct are for a SINGLE element, NOT a
     region, as such the only change that having variable bedfriction 
     requires is to pass the bedfriction angle for the current element 
     rather than the only bedfriction 

     I wonder if this is legacy code, it seems odd that it is only called
     for the SUN Operating System zee ../geoflow/correct.f */

#ifdef DO_EROSION
  int do_erosion=1;
#else
  int do_erosion=0;
#endif

#ifdef STOPCRIT_CHANGE_SOURCE
  int IF_STOPPED=EmTemp->get_stoppedflags();
#else
  int IF_STOPPED=!(!EmTemp->get_stoppedflags());
#endif

  double *state_vars=EmTemp->get_state_vars();
  double *prev_state_vars=EmTemp->get_prev_state_vars();
  double *d_state_vars=EmTemp->get_d_state_vars();
  double *gravity=EmTemp->get_gravity();
  double *d_gravity=EmTemp->get_d_gravity();
  double *zeta=EmTemp->get_zeta();
  double *curvature=EmTemp->get_curvature();
  double effect_bedfrict=EmTemp->get_effect_bedfrict();
  double *effect_kactxy=EmTemp->get_effect_kactxy(); 
  double *Influx=EmTemp->get_influx();


  double VxVy[2];
  if(state_vars[0]>GEOFLOW_TINY) {
    VxVy[0]=state_vars[1]/state_vars[0];
    VxVy[1]=state_vars[2]/state_vars[0];
  }
  else{VxVy[0]=VxVy[1]=0.0;}

  
  EmTemp->convect_dryline(VxVy,dt); //this is necessary

  correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym,
	   &tiny,&dtdx,&dtdy,&dt,d_state_vars,(d_state_vars+NUM_STATE_VARS), 
	   &(zeta[0]),&(zeta[1]),curvature,
	   &(matprops_ptr->intfrict), &effect_bedfrict,
	   gravity, effect_kactxy, d_gravity, &(matprops_ptr->frict_tiny),
	   forceint,forcebed,&do_erosion,eroded,VxVy,
	   &IF_STOPPED,Influx);

  *forceint*=dx[0]*dx[1];
  *forcebed*=dx[0]*dx[1];
  *eroded*=dx[0]*dx[1];

#endif 

  if(EmTemp->get_stoppedflags()==2) 
    *deposited=state_vars[0]*dx[0]*dx[1];
  else
    *deposited=0.0;

  if(EmTemp->get_stoppedflags()) 
    *eroded=0.0;
  
  EmTemp->calc_shortspeed(1.0/dt);




  return;
}
