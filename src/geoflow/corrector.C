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

//
//dUdx[0] dh_dx
//dUdx[1] dhVx_dx
//dUdx[2] dhVy_dx
//dUdy[0] dh_dy
//dUdy[1] dhVx_dy
//dUdy[2] dhVy_dy
//
void correct1ph(Element *Elm, const double *fluxxp,
        const double *fluxyp, const double *fluxxm, const double *fluxym,
        const double tiny, const double dtdx, const double dtdy, const double dt,
        const double dh_dx, const double dhVy_dx, 
        const double dh_dy, const double dhVx_dy,
        const double xslope, const double yslope, 
        const double curv_x, const double curv_y, const double intfrictang,
        const double bedfrictang, const double *g, const double kactxy,
        const double *dgdx, const double frict_tiny, double &forceint,
        double &forcebed, const int DO_EROSION, double &eroded,
        const double *VxVy, const int IF_STOPPED)
{
    double speed;
    double forceintx, forceinty;
    double forcebedx, forcebedy;
    double forcebedmax, forcebedequil, forcegrav;
    double unitvx, unitvy;
    double tanbed;
    double Ustore[3];

    double h_inv;
    double sgn_dudy, sgn_dvdx, tmp;
    double slope;
    double es, totalShear;
    const double threshold = 1.0E-02, erosion_rate = 0.1;

    // parameter(erosion_rate=0.015)(data for this is in ...DATA dir)
    // parameter(erosion_rate=0.125).....2.983063e+05.........5.718950e+05.....dec_2_erosion_4 
    // parameter(erosion_rate=0.250) ....2.983063e+05.........2.582245e+06.....dec_2_erosion_5
    // parameter(erosion_rate=0.19)......2.983063e+05........ 1.183863e+06.....dec_10_erosion_2.
    // parameter(erosion_rate=0.17)......2.983063e+05........ 9.570268e+05.....dec_10_erosion_3.
    // parameter(erosion_rate=0.0) ......2.983063e+05.........2.983062e+05.....dec_11_erosion1.
    // parameter(erosion_rate=0.16)......2.983063e+05.........9.306538e+05.....dec_11_erosion2.
    // parameter(erosion_rate=0.14)......2.983063e+05.........6.282101e+05.....dec_11_erosion3.
    // parameter(erosion_rate=0.15)
    // ......2.983063e+05.........7.039735e+05.....dec_11_erosion4.
    // parameter(erosion_rate=0.155).....2.983063e+05.........7.390704e+05.....dec_11_erosion5.
    // parameter(erosion_rate=0.156)
    // external sgn_tiny            

    slope = sqrt(xslope * xslope + yslope * yslope);
    //
    // threshold=1.278820338*dcos(slope)*
    //     &     c_dmax1((1-dtan(slope)/dtan(intfrictang)),0.0)**2

    Ustore[0] = Elm->prev_state_vars(0)
            - dtdx * (fluxxp[0] - fluxxm[0])
            - dtdy * (fluxyp[0] - fluxym[0])
            + dt * Elm->Influx(0);
    Ustore[0] = c_dmax1(Ustore[0], 0.0);

    Ustore[1] = Elm->prev_state_vars(1)
            - dtdx * (fluxxp[1] - fluxxm[1])
            - dtdy * (fluxyp[1] - fluxym[1])
            + dt * Elm->Influx(1);

    Ustore[2] = Elm->prev_state_vars(2)
            - dtdx * (fluxxp[2] - fluxxm[2])
            - dtdy * (fluxyp[2] - fluxym[2])
            + dt * Elm->Influx(2);

    // initialize to zero
    forceintx = 0.0;
    forcebedx = 0.0;
    forceinty = 0.0;
    forcebedy = 0.0;
    unitvx = 0.0;
    unitvy = 0.0;
    eroded = 0.0;

    if (Elm->state_vars(0) > tiny) {
        // S terms
        // here speed is speed squared
        speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];
        if (speed > 0.0) {
            // here speed is speed
            speed = sqrt(speed);
            unitvx = VxVy[0] / speed;
            unitvy = VxVy[1] / speed;
        }
        else {
            unitvx = 0.0;
            unitvy = 0.0;
        }
        tanbed = tan(bedfrictang);
        h_inv = 1.0 / Elm->state_vars(0);

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        // x direction source terms
        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        // the gravity force in the x direction
        forcegrav = g[0] * Elm->state_vars(0);

        // the internal friction force
        tmp = h_inv * (dhVx_dy - VxVy[0] * dh_dy);
        sgn_dudy = sgn_tiny(tmp, frict_tiny);
        forceintx = sgn_dudy * Elm->state_vars(0) * kactxy * (g[2] * dh_dy + dgdx[1] * Elm->state_vars(0)) * sin(intfrictang);

        // the bed friction force for fast moving flow 
        forcebedx = unitvx * c_dmax1(g[2] * Elm->state_vars(0) + VxVy[0] * Elm->state_vars(1) * curv_x, 0.0) * tanbed;

        if (IF_STOPPED == 2 && 1 == 0) {
            // the bed friction force for stopped or nearly stopped flow

            // the static friction force is LESS THAN or equal to the friction
            // coefficient times the normal force but it can NEVER exceed the 
            // NET force it is opposing

            // maximum friction force the bed friction can support
            forcebedmax = g[2] * Elm->state_vars(0) * tanbed;

            // the NET force the bed friction force is opposing 
            forcebedequil = forcegrav - forceintx;
            // $           -kactxy*g[2]*EmTemp->state_vars(0)*dh_dx

            // the "correct" stopped or nearly stopped flow bed friction force 
            // (this force is not entirely "correct" it will leave a "negligible"
            // (determined by stopping criteria) amount of momentum in the cell
            forcebedx = sgn_tiny(forcebedequil, c_dmin1(forcebedmax, fabs(forcebedx) + fabs(forcebedequil)));
            // forcebedx=sgn_tiny(forcebed2,dmin1(forcebed1,fabs(forcebed2)))

            // not really 1 but this makes friction statistics accurate
            unitvx = 1.0;
            // else

        }

        Ustore[1] = Ustore[1] + dt * (forcegrav - forcebedx - forceintx);

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        // y direction source terms
        //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        // the gravity force in the y direction
        forcegrav = g[1] * Elm->state_vars(0);

        // the internal friction force
        tmp = h_inv * (dhVy_dx - VxVy[1] * dh_dx);
        sgn_dvdx = sgn_tiny(tmp, frict_tiny);
        forceinty = sgn_dvdx * Elm->state_vars(0) * kactxy * (g[2] * dh_dx + dgdx[0] * Elm->state_vars(0)) * sin(intfrictang);

        // the bed friction force for fast moving flow 
        forcebedy = unitvy * c_dmax1(g[2] * Elm->state_vars(0) + VxVy[1] * Elm->state_vars(2) * curv_y, 0.0) * tanbed;
        if (IF_STOPPED == 2 && 1 == 0) {
            // the bed friction force for stopped or nearly stopped flow

            // the NET force the bed friction force is opposing 
            forcebedequil = forcegrav - forceinty;
            // $           -kactxy*g[2]*EmTemp->state_vars(0)*dh_dy

            // the "correct" stopped or nearly stopped flow bed friction force 
            // (this force is not entirely "correct" it will leave a "negligible"
            // (determined by stopping criteria) amount of momentum in the cell
            forcebedy = sgn_tiny(forcebedequil, c_dmin1(forcebedmax, fabs(forcebedy) + fabs(forcebedequil)));

            // not really 1 but this makes friction statistics accurate
            unitvy = 1.0;
            //    else
        }

        Ustore[2] = Ustore[2] + dt * (forcegrav - forcebedy - forceinty);

        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        // (erosion terms) this is Camil's logic, Keith changed some variable 
        //names for clarity
        //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if ((false) && (DO_EROSION != 0) && (IF_STOPPED == 0)) {
            totalShear = sqrt(forcebedx * forcebedx + forcebedy * forcebedy);
            if ((totalShear > threshold) && (Elm->state_vars(0) > 0.004)) {

                es = erosion_rate * sqrt(fabs(totalShear - threshold));
                eroded = dt*es;
                Ustore[0] = Ustore[0] + eroded;
                Ustore[1] = Ustore[1] + eroded * VxVy[0];
                Ustore[2] = Ustore[2] + eroded * VxVy[1];
                //write (*,*) 'Doing Keith Erosion Model'
            }
        }
        if ((DO_EROSION != 0) && (Elm->state_vars(0) > threshold)) {
            es = erosion_rate * sqrt(Elm->state_vars(1) * Elm->state_vars(1) + Elm->state_vars(2) * Elm->state_vars(2)) / Elm->state_vars(0);
            Ustore[0] = Ustore[0] + dt*es;
            Ustore[1] = Ustore[1] + dt * es * Ustore[1];
            Ustore[2] = Ustore[2] + dt * es * Ustore[2];
            //write (*,*) 'Doing Camil Erosion Model'
        }
    }

    // computation of magnitude of friction forces for statistics
    forceint = unitvx * forceintx + unitvy*forceinty;
    forcebed = unitvx * forcebedx + unitvy*forcebedy;

    // update the state variables
    Elm->state_vars(0, Ustore[0]);
    Elm->state_vars(1, Ustore[1]);
    Elm->state_vars(2, Ustore[2]);
}
void calc_drag_force(Element *Elm, const double *vsolid, const double *vfluid,
        const double den_solid, const double den_fluid, const double vterminal,
        double *drag, const double tiny)
{
      double temp, volf, denfrac;
      double delv[2];
      const double exponant = 3.0;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!     model in pitman-le paper
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      volf = Elm->state_vars(1)/Elm->state_vars(0);
      temp = Elm->state_vars(1)*pow(1.-volf,1.0-exponant)/vterminal;
      denfrac = den_fluid/den_solid;

      for(int i=0;i<4;++i)drag[i]=0.0;
 
      //fluid vel - solid vel
      if(Elm->state_vars(1) > tiny)
      {
        for(int i=0;i<2;++i)
            delv[i] = vfluid[i]-vsolid[i];

//compute individual drag-forces
        drag[0] = (1.0-denfrac)*temp*delv[0];
        drag[1] = (1.0-denfrac)*temp*delv[1];
        drag[2] = (1.0-denfrac)*temp*delv[0]/denfrac;
        drag[3] = (1.0-denfrac)*temp*delv[1]/denfrac;
      }
}


//dUdx[0] dh_dx
//dUdx[1] dh_dx_liq
//dUdx[3] dhVy_dx_sol
//dUdy[0] dh_dy
//dUdy[1] dh_dy_liq
//dUdy[2] dh_dy_liq
void correct2ph(Element *Elm, const double *fluxxp,
        const double *fluxyp, const double *fluxxm, const double *fluxym,
        const double tiny, const double dtdx, const double dtdy, const double dt,
        const double dh_dx, const double dh_dx_liq, const double dhVy_dx_sol, 
        const double dh_dy, const double dh_dy_liq, const double dhVx_dy_sol,
        const double xslope,const double yslope, const double curv_x, const double curv_y, const double intfrictang,
        const double bedfrictang, const double *g, const double kactxy,
        const double frict_tiny, double &forceint,
        double &forcebed, const int DO_EROSION, double &eroded,
        const double *v_solid, const double *v_fluid,
        const double den_solid, const double den_fluid, const double terminal_vel, 
        const double eps, const int IF_STOPPED)
{
      double speed;
      double forceintx, forceinty;
      double forcebedx, forcebedy;
      double forcebedmax, forcebedequil, forcegrav;
      double unitvx, unitvy;
      double den_frac;
      double alphaxx, alphayy, alphaxy, alphaxz, alphayz;
      double tanbed;

      double Ustore[6];
      double h_inv, hphi_inv;
      double sgn_dudy, sgn_dvdx, tmp;
      double slope;
      double t1, t2, t3, t4, t5;
      double es,totalShear;
      double drag[4];

      int i;
      const double threshold=1.0E-02,erosion_rate=0.1;

// initialize to zero
      forceintx=0.0;
      forcebedx=0.0;
      forceinty=0.0;
      forcebedy=0.0;
      unitvx=0.0;
      unitvy=0.0;
      eroded=0.0;

      slope=sqrt(xslope*xslope+yslope*yslope);
      den_frac = den_fluid/den_solid;
      for(i=0;i<6;++i)
        Ustore[i]=Elm->prev_state_vars(i)+dt*Elm->Influx(i)-dtdx*(fluxxp[i]-fluxxm[i])-dtdy*(fluxyp[i]-fluxym[i]);

      if(Ustore[0] > tiny) {
// Source terms ...
// here speed is speed squared
         speed=v_solid[0]*v_solid[0]+v_solid[1]*v_solid[1];
         if(speed>0.0) {
// here speed is speed
            speed=sqrt(speed);
            unitvx=v_solid[0]/speed;
            unitvy=v_solid[1]/speed;
         }
         else{
            unitvx=0.0;
            unitvy=0.0;
         }
         tanbed=tan(bedfrictang);
         h_inv = 1.0/Elm->state_vars(0);
         hphi_inv = 1.0/Elm->state_vars(1);
         alphaxx = kactxy;
         alphayy = kactxy;
         den_frac = den_fluid/den_solid;
         calc_drag_force(Elm, v_solid, v_fluid, den_solid, den_fluid, terminal_vel, drag, tiny);
 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    solid fraction x-direction source terms
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    alphaxy -- see pitman-le (2005)
         tmp = hphi_inv*(dhVx_dy_sol-v_solid[0]*dh_dy_liq);
         sgn_dudy = sgn_tiny(tmp, frict_tiny);
         alphaxy = sgn_dudy*sin(intfrictang)*kactxy;

//    alphaxz (includes centrifugal effects)
         alphaxz = -unitvx*tanbed*(1.0+(v_solid[0]*v_solid[0])*curv_x/g[2]);

//    evaluate t1 
         t1 = (1.0-den_frac)*(-alphaxx*xslope-alphaxy*yslope + alphaxz)*Elm->state_vars(1)*g[2];
//    evaluate t2
         t2 = eps*den_frac*Elm->state_vars(1)*g[2]*dh_dx;
//    evaluate t3
         t3 = eps*den_frac*Elm->state_vars(1)*g[2]*xslope;
//    evaluate t4
         t4 = Elm->state_vars(1)*g[0];
//    evaluate drag
         t5 = drag[0];
//    update Ustore
         Ustore[2] = Ustore[2] + dt*(t1-t2-t3+t4+t5);

//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// solid fraction y-direction source terms
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    alphaxy -- see pitman-le (2005) for definitions
         tmp = hphi_inv*(dhVy_dx_sol-v_solid[1]*dh_dx_liq);
         sgn_dvdx = sgn_tiny(tmp, frict_tiny);
         alphaxy = sgn_dvdx*sin(intfrictang)*kactxy;

//    alphayz
         alphayz = -unitvy*tanbed*(1.0+(v_solid[1]*v_solid[1])*curv_y/g[2]);

//    evaluate t1
         t1 = (1.0-den_frac)*(-alphaxy*xslope-alphayy*yslope + alphayz)*Elm->state_vars(1)*g[2];
//    evaluate t2
         t2 = eps*den_frac*Elm->state_vars(1)*dh_dy;
//    evaluate t3
         t3 = eps*den_frac*Elm->state_vars(1)*g[2]*yslope;
//    evaluate t4 ( gravity along y-dir )
         t4 = Elm->state_vars(1)*g[1];
//    drag term
         t5 = drag[1];
         Ustore[3] = Ustore[3] + dt*(t1-t2-t3+t4+t5);

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    fluid fraction x-direction source terms
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    gravity on fluid
         t4 = Elm->state_vars(0)*g[0];
//    drag force on fluid
         t5 = drag[2];
         Ustore[4] = Ustore[4] + dt*(t4 - t5);

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    fluid fraction y-direction source terms
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//    gravity on fluid
         t4 = Elm->state_vars(0)*g[1];
//    drag force on fluid
         t5 = drag[3];
         Ustore[5] = Ustore[5] + dt*(t4 - t5);
      }

// computation of magnitude of friction forces for statistics
      forceint=unitvx*forceintx+unitvy*forceinty;
      forcebed=unitvx*forcebedx+unitvy*forcebedy;

// update the state variables
      for(i=0;i<6;++i)
         Elm->state_vars(i, Ustore[i]);

}
void correct(ElementType elementType,HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr, FluxProps *fluxprops,
             TimeProps *timeprops, void *EmTemp_in, double *forceint, double *forcebed, double *eroded,
             double *deposited)
{
    MatPropsTwoPhases* matprops2_ptr{nullptr};
    if(elementType == ElementType::TwoPhases)
    {
        matprops2_ptr=static_cast<MatPropsTwoPhases*>(matprops_ptr);
    }
    Element *EmTemp = (Element *) EmTemp_in;
    double dxdy = EmTemp->dx(0)*EmTemp->dx(1);
    double dtdx = dt / EmTemp->dx(0);
    double dtdy = dt / EmTemp->dx(1);
    
    double tiny = GEOFLOW_TINY;
    int xp = EmTemp->positive_x_side();
    int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
    
    int ivar, j, k;
    
    double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
    double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];
    
    Node** nodes = EmTemp->getNodesPtrs();
    
    Node* nxp = nodes[xp + 4]; //(Node*) NodeTable->lookup(EmTemp->getNode() + (xp + 4) * 2);
    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxxp[ivar] = nxp->flux(ivar);
    
    Node* nyp = nodes[yp + 4]; //(Node*) NodeTable->lookup(EmTemp->getNode() + (yp + 4) * 2);
    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxyp[ivar] = nyp->flux(ivar);
    
    Node* nxm = nodes[xm + 4]; //(Node*) NodeTable->lookup(EmTemp->getNode() + (xm + 4) * 2);
    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxxm[ivar] = nxm->flux(ivar);
    
    Node* nym = nodes[ym + 4]; //(Node*) NodeTable->lookup(EmTemp->getNode() + (ym + 4) * 2);
    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxym[ivar] = nym->flux(ivar);
    

    /* the values being passed to correct are for a SINGLE element, NOT a
     region, as such the only change that having variable bedfriction
     requires is to pass the bedfriction angle for the current element
     rather than the only bedfriction

     I wonder if this is legacy code, it seems odd that it is only called
     for the SUN Operating System zee ../geoflow/correct.f */

#ifdef DO_EROSION
    int do_erosion=1;
#else
    int do_erosion = 0;
#endif
    
#ifdef STOPCRIT_CHANGE_SOURCE
    int IF_STOPPED=EmTemp->stoppedflags();
#else
    int IF_STOPPED = !(!EmTemp->stoppedflags());
#endif
    
    double gravity[3]{EmTemp->gravity(0),EmTemp->gravity(1),EmTemp->gravity(2)};
    double d_gravity[3]{EmTemp->d_gravity(0),EmTemp->d_gravity(1),EmTemp->d_gravity(2)};
    
    if(elementType == ElementType::TwoPhases)
    {
        int i;
        double kactxy[DIMENSION];
        double bedfrict = EmTemp->effect_bedfrict();
        double solid_den = matprops2_ptr->den_solid;
        double fluid_den = matprops2_ptr->den_fluid;
        double terminal_vel = matprops2_ptr->v_terminal;

        double Vfluid[DIMENSION];
        double volf;
        if(EmTemp->state_vars(0) > GEOFLOW_TINY)
        {
            for(i = 0; i < DIMENSION; i++)
                kactxy[i] = EmTemp->effect_kactxy(i);

            // fluid velocities
            Vfluid[0] = EmTemp->state_vars(4) / EmTemp->state_vars(0);
            Vfluid[1] = EmTemp->state_vars(5) / EmTemp->state_vars(0);

            // volume fractions
            volf = EmTemp->state_vars(1) / EmTemp->state_vars(0);
        }
        else
        {
            for(i = 0; i < DIMENSION; i++)
            {
                kactxy[i] = matprops2_ptr->epsilon;
                Vfluid[i] = 0.;
            }
            volf = 1.;
            bedfrict = matprops2_ptr->bedfrict[EmTemp->material()];
        }

        double Vsolid[DIMENSION];
        if(EmTemp->state_vars(1) > GEOFLOW_TINY)
        {
            Vsolid[0] = EmTemp->state_vars(2) / EmTemp->state_vars(1);
            Vsolid[1] = EmTemp->state_vars(3) / EmTemp->state_vars(1);
        }
        else
        {
            Vsolid[0] = Vsolid[1] = 0.0;
        }

        double V_avg[DIMENSION];
        V_avg[0] = Vsolid[0] * volf + Vfluid[0] * (1. - volf);
        V_avg[1] = Vsolid[1] * volf + Vfluid[1] * (1. - volf);
        EmTemp->convect_dryline(V_avg[0],V_avg[1], dt); //this is necessary

        /*correct2ph_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt, d_state_vars,
                 (d_state_vars + NUM_STATE_VARS), &(zeta[0]), &(zeta[1]), curvature, &(matprops2_ptr->intfrict), &bedfrict,
                 gravity, kactxy, &(matprops2_ptr->frict_tiny), forceint, forcebed, &do_erosion, eroded, Vsolid, Vfluid,
                 &solid_den, &fluid_den, &terminal_vel, &(matprops2_ptr->epsilon), &IF_STOPPED, Influx);*/
        
        correct2ph(EmTemp, fluxxp, fluxyp, fluxxm, fluxym, tiny, dtdx, dtdy, dt,
                EmTemp->dh_dx(), EmTemp->dh_dx_liq(), EmTemp->dhVy_dx_sol(), 
                EmTemp->dh_dy(), EmTemp->dh_dy_liq(), EmTemp->dhVx_dy_sol(),
                EmTemp->zeta(0), EmTemp->zeta(1), EmTemp->curvature(0),EmTemp->curvature(1), matprops2_ptr->intfrict, bedfrict,
                gravity, kactxy[0], matprops2_ptr->frict_tiny, *forceint, *forcebed, do_erosion, *eroded, Vsolid, Vfluid,
                solid_den, fluid_den, terminal_vel, matprops2_ptr->epsilon, IF_STOPPED);
        
        
        bool print_vars = false;
        for(i = 0; i < NUM_STATE_VARS; i++)
            if(isnan(EmTemp->state_vars(i)))
                print_vars = true;
        
        if(print_vars)
        {
            cout<<"ElemKey: "<<EmTemp->key()<<endl;
            printf("Kactxy = %10.5f%10.5f\n", kactxy[0], kactxy[1]);
            printf("BedFrict: %10.5f: IntFrict: %10.5f\n", bedfrict, matprops2_ptr->intfrict);
            printf("state_vars: \n");
            for(i = 0; i < NUM_STATE_VARS; i++)
                printf("%10.5f", EmTemp->state_vars(i));
            printf("\n");
            printf("prev_state_vars: \n");
            for(i = 0; i < NUM_STATE_VARS; i++)
                printf("%10.5f", EmTemp->prev_state_vars(i));
            printf("\n");
            printf("fluxes: \n");
            for(i = 0; i < NUM_STATE_VARS; i++)
                printf("%10.5f%10.5f%10.5f%10.5f\n", fluxxp[i], fluxxm[i], fluxyp[i], fluxym[i]);
        }
    }
    if(elementType == ElementType::SinglePhase)
    {
        double VxVy[2];
        if(EmTemp->state_vars(0) > GEOFLOW_TINY)
        {
            VxVy[0] = EmTemp->state_vars(1) / EmTemp->state_vars(0);
            VxVy[1] = EmTemp->state_vars(2) / EmTemp->state_vars(0);
        }
        else
        {
            VxVy[0] = VxVy[1] = 0.0;
        }

        EmTemp->convect_dryline(VxVy[0], VxVy[1], dt); //this is necessary
        
        /*void correct1ph(double *Uvec, const double *Uprev, const double *fluxxp,
        const double *fluxyp, const double *fluxxm, const double *fluxym,
        const double tiny, const double dtdx, const double dtdy, const double dt,
        const double *dUdx, const double *dUdy, const double xslope,
        const double yslope, const double *curv, const double intfrictang,
        const double bedfrictang, const double *g, const double kactxy,
        const double *dgdx, const double frict_tiny, double &forceint,
        double &forcebed, const int DO_EROSION, double &eroded,
        const double *VxVy, const int IF_STOPPED, const double *fluxsrc);*/
        
        correct1ph(EmTemp, fluxxp, fluxyp, fluxxm, fluxym, tiny, dtdx, dtdy, dt,
                EmTemp->dh_dx(), EmTemp->dhVy_dx(), 
                EmTemp->dh_dy(), EmTemp->dhVx_dy(),
                EmTemp->zeta(0), EmTemp->zeta(1), EmTemp->curvature(0),EmTemp->curvature(1), matprops_ptr->intfrict,
                EmTemp->effect_bedfrict(), gravity, EmTemp->effect_kactxy(0), d_gravity, matprops_ptr->frict_tiny, *forceint, *forcebed,
                do_erosion, *eroded, VxVy, IF_STOPPED);

        /*correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt, d_state_vars,
                 (d_state_vars + NUM_STATE_VARS), &(zeta[0]), &(zeta[1]), curvature, &(matprops_ptr->intfrict),
                 &effect_bedfrict, gravity, effect_kactxy, d_gravity, &(matprops_ptr->frict_tiny), forceint, forcebed,
                 &do_erosion, eroded, VxVy, &IF_STOPPED, Influx);*/
    }
    
    *forceint *= dxdy;
    *forcebed *= dxdy;
    *eroded *= dxdy;

    
    if(EmTemp->stoppedflags() == 2)
        *deposited = EmTemp->h() * dxdy;
    else
        *deposited = 0.0;
    
    if(EmTemp->stoppedflags())
        *eroded = 0.0;
    
    EmTemp->calc_shortspeed(1.0 / dt);
    
    return;
}
