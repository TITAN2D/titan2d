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

#define DO_EROSION
 
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
  double kactxy[DIMENSION]; 

  double tiny = GEOFLOW_TINY;
  int xp=EmTemp->get_positive_x_side();
  int yp=(xp+1)%4, xm=(xp+2)%4, ym=(xp+3)%4; 

  int ivar,i, j, k;
  double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS]; 
  double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];

  Node* nxp = (Node*) NodeTable->lookup(EmTemp->getNode()+(xp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxp[ivar] = nxp->flux[ivar];
    
  Node* nyp = (Node*) NodeTable->lookup(EmTemp->getNode()+(yp+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxyp[ivar] = nyp->flux[ivar];
    
  Node* nxm = (Node*) NodeTable->lookup(EmTemp->getNode()+(xm+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxxm[ivar] = nxm->flux[ivar];
    
  Node* nym = (Node*) NodeTable->lookup(EmTemp->getNode()+(ym+4)*2);
  for(ivar=0;ivar<NUM_STATE_VARS;ivar++) fluxym[ivar] = nym->flux[ivar];
    

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
  double bedfrict=EmTemp->get_effect_bedfrict();
  double *Influx=EmTemp->get_influx();
  double solid_den=matprops_ptr->den_solid;
  double fluid_den=matprops_ptr->den_fluid;
  double terminal_vel=matprops_ptr->v_terminal;


  double Vfluid[DIMENSION];
  double volf;
  if ( state_vars[0] > GEOFLOW_TINY)
  {
    for (i=0; i<DIMENSION; i++)
     kactxy[i]=*(EmTemp->get_effect_kactxy()+i);

    // fluid velocities
    Vfluid[0]=state_vars[4]/state_vars[0];
    Vfluid[1]=state_vars[5]/state_vars[0];

    // volume fractions
    volf = state_vars[1]/state_vars[0];
  }
  else
  {
    for (i=0; i<DIMENSION; i++)
    {
      kactxy[i]=matprops_ptr->epsilon;
      Vfluid[i]=0.;
    }
    volf=1.;
    bedfrict=matprops_ptr->bedfrict[EmTemp->get_material()];
  }

  double Vsolid[DIMENSION];
  if(state_vars[1]>GEOFLOW_TINY) 
  {
    Vsolid[0]=state_vars[2]/state_vars[1];
    Vsolid[1]=state_vars[3]/state_vars[1];
  }
  else
  {
    Vsolid[0]=Vsolid[1]=0.0;
  }

  double V_avg[DIMENSION];
  V_avg[0] = Vsolid[0]*volf + Vfluid[0]*(1.-volf);
  V_avg[1] = Vsolid[1]*volf + Vfluid[1]*(1.-volf);
  EmTemp->convect_dryline(V_avg,dt); //this is necessary

  correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym,
	   &tiny, &dtdx, &dtdy, &dt, d_state_vars, (d_state_vars+NUM_STATE_VARS), 
	   &(zeta[0]), &(zeta[1]), curvature,
	   &(matprops_ptr->intfrict), &bedfrict,
	   gravity, kactxy, &(matprops_ptr->frict_tiny),
	   forceint, forcebed, &do_erosion, eroded, Vsolid, Vfluid,
	   &solid_den, &fluid_den, &terminal_vel,
           &(matprops_ptr->epsilon), &IF_STOPPED, Influx);

  *forceint*=dx[0]*dx[1];
  *forcebed*=dx[0]*dx[1];
  *eroded*=dx[0]*dx[1];

   bool print_vars=false;
   for (i=0; i<NUM_STATE_VARS; i++)
     if ( isnan(state_vars[i]) )
        print_vars=true;

  if ( print_vars )
  {
    printf("ElemKey: %u\n", *EmTemp->pass_key());
    printf("Kactxy = %10.5f%10.5f\n", kactxy[0], kactxy[1]);
    printf("BedFrict: %10.5f: IntFrict: %10.5f\n", bedfrict, matprops_ptr->intfrict);
    printf("state_vars: \n");
    for (i=0; i<NUM_STATE_VARS; i++)
      printf("%10.5f", state_vars[i]);
    printf("\n");
    printf("prev_state_vars: \n");
    for (i=0; i<NUM_STATE_VARS; i++)
      printf("%10.5f", prev_state_vars[i]);
    printf("\n");
    printf("fluxes: \n");
    for (i=0; i<NUM_STATE_VARS; i++)
      printf("%10.5f%10.5f%10.5f%10.5f\n",
             fluxxp[i],fluxxm[i],fluxyp[i],fluxym[i]);
  }

  if(EmTemp->get_stoppedflags()==2) 
    *deposited=state_vars[0]*dx[0]*dx[1];
  else
    *deposited=0.0;

  if(EmTemp->get_stoppedflags()) 
    *eroded=0.0;
  
  EmTemp->calc_shortspeed(1.0/dt);
  return;
}
