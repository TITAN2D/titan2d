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
 * $Id: step.C 225 2012-02-06 21:05:29Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#define APPLY_BC

#include "../header/titan2d_utils.h"

#include "../header/outline.h"

//dUdx[0] dh_dx
//dUdx[1] dhVx_dx
//dUdx[2] dhVy_dx
//dUdy[0] dh_dy
//dUdy[1] dhVx_dy
//dUdy[2] dhVy_dy

//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
void predict(Element *Elm,
        const double dh_dx, const double dhVx_dx, const double dhVy_dx, 
        const double dh_dy, const double dhVx_dy, const double dhVy_dy,
        const double tiny, const double kactx,
        const double dt2, const double *g, const double curv_x, const double curv_y,
        const double bedfrictang, const double intfrictang,
        const double *dgdx, const double frict_tiny, const int order_flag,
        double *VxVy, const int IF_STOPPED, double *fluxsrc)
{
    //NOTE:  d(g[2]*Elm->state_vars(0))/dx is approximated by g[2]*dUvec[0]/dx !!!
    double c_sq;
    double h_inv;
    double tanbed;
    double VxVyS[2];
    double VxVyB[2];
    double unitvx, unitvy;
    double tmp, sgn_dudy,sgn_dvdx;
    double forcegrav;
    double forceintx,forceinty, forcebedx,forcebedy;
    double forcebedequil, forcebedmax;
    double speed;
    //    curv := inverse of radius of curvature = second derivative of 
    //    position normal to tangent with respect to distance along tangent,
    //    if dz/dx=0 curve=d2z/dx2, otherwise rotate coordinate system so 
    //    dz/dx=0, that is mathematical definition of curvature I believe
    //    laercio returns d2z/dx2 whether or not dz/dx=0 in his GIS functions

    //TEST********** order_flag = 1 
    //     write(*,*) "order_flag in predict=", order_flag 

    Elm->prev_state_vars(0, Elm->state_vars(0));
    Elm->prev_state_vars(1, Elm->state_vars(1));
    Elm->prev_state_vars(2, Elm->state_vars(2));

    if (IF_STOPPED == 2) {
        VxVy[0] = 0.0;
        VxVy[1] = 0.0;
        VxVyS[0] = 0.0;
        VxVyS[1] = 0.0;
    }
    else {
        VxVy[0] = VxVyB[0];
        //Elm->state_vars(1)/Elm->state_vars(0);
        VxVy[1] = VxVyB[1];
        //Elm->state_vars(2)/Elm->state_vars(0);
        VxVyS[0] = VxVyB[0];
        VxVyS[1] = VxVyB[1];
    }

    if (order_flag == 2) {

        c_sq = kactx * g[2] * Elm->state_vars(0);
        //h_inv := 1/Elm->state_vars(0)                         

        Elm->state_vars(0, Elm->state_vars(0) - dt2 * (dhVx_dx + dhVy_dy + fluxsrc[0]));
        Elm->state_vars(0, c_dmax1(Elm->state_vars(0), 0.0));

        //dF/dU, dG/dU and S terms if Elm->state_vars(0) > TINY !
        if (Elm->prev_state_vars(0) > tiny) {
            h_inv = 1.0 / Elm->prev_state_vars(0);
            tanbed = tan(bedfrictang);

            //here speed is speed squared
            speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];
            if (speed > 0.0) {
                //here speed is speed
                speed = sqrt(speed);
                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else {
                unitvx = 0.0;
                unitvy = 0.0;
            }

            //dnorm=dsqrt(Uprev[1]**2+Uprev[2]**2+tiny**2)

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //****** X-dir ******
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // dF/dU and dG/dU terms
            Elm->state_vars(1, Elm->state_vars(1) -
                    dt2 * ((c_sq - VxVy[0] * VxVy[0]) * dh_dx +
                    2.0 * VxVy[0] * dhVx_dx -
                    VxVy[0] * VxVy[1] * dh_dy +
                    VxVy[1] * dhVx_dy +
                    VxVy[0] * dhVy_dy +
                    fluxsrc[1]));

            // x direction source terms

            // the gravity force in the x direction
            forcegrav = g[0] * Elm->state_vars(0);

            // the internal friction force
            tmp = h_inv * (dhVx_dy - VxVyS[0] * dh_dy);
            sgn_dudy = sgn_tiny(tmp, frict_tiny);
            forceintx = sgn_dudy * Elm->state_vars(0) * kactx * (g[2] * dh_dy + dgdx[1] * Elm->state_vars(0)) * sin(intfrictang);

            // the bed friction force for fast moving flow 
            forcebedx = unitvx * c_dmax1(g[2] * Elm->state_vars(0) + VxVyS[0] * Elm->state_vars(1) * curv_x, 0.0) * tanbed;

            if (IF_STOPPED == 2 && 1 == 0) {
                //the bed friction force for stopped or nearly stopped flow

                //the static friction force is LESS THAN or equal to the friction
                //coefficient times the normal force but it can NEVER exceed the 
                //NET force it is opposing

                //maximum friction force the bed friction can support
                forcebedmax = c_dmax1(g[2] * Elm->state_vars(0) + VxVyS[0] * Elm->state_vars(1) * curv_x, 0.0) * tanbed;

                //     the NET force the bed friction force is opposing 
                forcebedequil = forcegrav - forceintx;
                //                   -kactxy*g[2]*Elm->state_vars(0)*dh_dx


                // the "correct" stopped or nearly stopped flow bed friction force 
                // (this force is not entirely "correct" it will leave a "negligible"
                // (determined by stopping criteria) amount of momentum in the cell
                forcebedx = sgn_tiny(forcebedequil, c_dmin1(forcebedmax,
                        dabs(forcebedx) + dabs(forcebedequil)));
                //  forcebedx=
                //                   sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))
                //            else
            }

            // all the x source terms
            Elm->state_vars(1, Elm->state_vars(1) + dt2 * (forcegrav - forcebedx - forceintx));
            //            write(*,*) 'int', forceintx, 'bed', forcebedx

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //****** Y-dir ******
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //dF/dU and dG/dU terms
            Elm->state_vars(2, Elm->state_vars(2) -
                    dt2 * ((c_sq - VxVy[1] * VxVy[1]) * dh_dy +
                    2.0 * VxVy[1] * dhVy_dy -
                    VxVy[0] * VxVy[1] * dh_dx +
                    VxVy[1] * dhVx_dx +
                    VxVy[0] * dhVy_dx +
                    fluxsrc[2]));

            //the gravity force in the y direction
            forcegrav = g[1] * Elm->state_vars(0);

            //the internal friction force
            tmp = h_inv * (dhVy_dx - VxVyS[1] * dh_dx);
            sgn_dvdx = sgn_tiny(tmp, frict_tiny);
            forceinty = sgn_dvdx * Elm->state_vars(0) * kactx * (g[2] *
                    dh_dx + dgdx[0] * Elm->state_vars(0)) * sin(intfrictang);

            //the bed friction force for fast moving flow 
            forcebedy = unitvy *
                    c_dmax1(g[2] * Elm->state_vars(0) + VxVyS[1] * Elm->state_vars(2) * curv_y, 0.0)
                    * tan(bedfrictang);

            if (IF_STOPPED == 2 && 1 == 0) {
                //the bed friction force for stopped or nearly stopped flow

                forcebedmax =
                        c_dmax1(g[2] * Elm->state_vars(0) + VxVyS[1] * Elm->state_vars(2) * curv_y, 0.0)
                        * tanbed;

                //the NET force the bed friction force is opposing 
                forcebedequil = forcegrav
                        //     $              -kactxy*g[2]*Elm->state_vars(0)*dh_dy
                        - forceinty;

                //the "correct" stopped or nearly stopped flow bed friction force 
                //(this force is not entirely "correct" it will leave a "negligible"
                //(determined by stopping criteria) amount of momentum in the cell
                forcebedy = sgn_tiny(forcebedequil, c_dmin1(forcebedmax,
                        fabs(forcebedy) + fabs(forcebedequil)));

                //          forcebedy=sgn(forcebed2,dmin1(forcebed1,dabs(forcebed2)))
                //           else
            }

            // all the y source terms
            Elm->state_vars(2, Elm->state_vars(2) + dt2 * (forcegrav - forcebedy - forceinty));

        }
    }
}
//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
void predict2ph(Element *Elm, const double *dUdx, const double *dUdy,
        double *Uprev, const double tiny, const double kactx,
        const double dt2, const double *g, const double curv_x, const double curv_y,
        const double bedfrictang, const double intfrictang,
        const double *dgdx, const double frict_tiny, const int order_flag,
        double *VxVy, const int IF_STOPPED, double *fluxsrc)
/*void predict(double *Uvec,dUdx,dUdy,Uprev,tiny,kactxy,dt2, g, 
          curv, bedfrictang, intfrictang,
          dgdx, frict_tiny, order_flag, VxVyB, 
          IF_STOPPED,fluxsrc)*/
{
    
}


void step(ElementType elementType,ElementsHashTable* El_Table, NodeHashTable* NodeTable, int myid, int nump, MatProps* matprops_ptr,
          TimeProps* timeprops_ptr, PileProps *pileprops_ptr, FluxProps *fluxprops, StatProps* statprops_ptr,
          int* order_flag, OutLine* outline_ptr, DischargePlanes* discharge, int adaptflag)
{
    double t_start,t_start2;
    
    /*El_Table->updateAllEntries();
     El_Table->updateAllLocalEntries();
     update_elements_pointers(El_Table, NodeTable);*/

    ASSERT3(El_Table->checkPointersToNeighbours("check_elements_pointers_StepStart")==0);
    ASSERT3(El_Table->ckeckLocalElementsPointers("ckeckAllLocalEntriesPointers_StepStart")==0);
    
    /* 
     * PREDICTOR-CORRECTED based on Davis' Simplified Godunov Method 
     */

    /* pass off proc data here (really only need state_vars for off-proc neighbors) */
    move_data(nump, myid, El_Table, NodeTable, timeprops_ptr);
    

    
    t_start = MPI_Wtime();
    slopes(El_Table, NodeTable, matprops_ptr);
    titanTimings.slopesCalcTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.slopesCalcTime += MPI_Wtime() - t_start;

    // get coefficients, eigenvalues, hmax and calculate the time step 
    double dt = get_coef_and_eigen(elementType, El_Table, NodeTable, matprops_ptr, fluxprops, timeprops_ptr, 0);
    //printf("step(): iter=%d %g+",timeprops_ptr->iter,timeprops_ptr->time);
    
    timeprops_ptr->incrtime(&dt); //also reduces dt if necessary
    //printf("step: oldtime+%g=%g\n",dt,timeprops_ptr->time);
    
    // assign influxes and then if any new sources are activating in current time step refine and re-mark cells 
    adapt_fluxsrc_region(El_Table, NodeTable, matprops_ptr, pileprops_ptr, fluxprops, timeprops_ptr, dt, myid,
                         adaptflag);
    
    int i;
    Element* Curr_El;
    double VxVy[2];
    
    //  move_data(nump, myid, El_Table, NodeTable,timeprops_ptr);
    
    //slopes(El_Table, NodeTable, matprops_ptr);
    
    /* get coefficients, eigenvalues, hmax and calculate the time step */
    //dt = get_coef_and_eigen(El_Table, NodeTable, matprops_ptr, 
    //			  fluxprops, timeprops_ptr,0);
    //timeprops_ptr->incrtime(&dt); //also reduces dt if necessary
    double dt2 = .5 * dt; // dt2 is set as dt/2 !
            
    /*
     *  predictor step
     */
    t_start = MPI_Wtime();
    int j, k, counter;
    double tiny = GEOFLOW_TINY;
    double flux_src_coef = 0;

// in TWO PHASES #ifdef SECOND_ORDER

    //-------------------go through all the elements of the subdomain and  
    //-------------------calculate the state variables at time .5*delta_t
    /* mdj 2007-04 */
    int IF_STOPPED;
    double curr_time, influx[3]; //VxVy[2];

    int Nelms = El_Table->getNumberOfLocalElements();
    Element** Elms = (Element**) El_Table->getLocalElementsValues();
    //Node* nd;
//#pragma omp parallel for                                                \
//private(currentPtr,Curr_El,IF_STOPPED,influx,j,k,curr_time,flux_src_coef,VxVy)
    for(i = 0; i < Nelms; i++)
    {
        
        Curr_El = Elms[i];
        
        influx[0] = Curr_El->Influx(0);
        influx[1] = Curr_El->Influx(1);
        influx[2] = Curr_El->Influx(2);
        
        //note, now there is no check for fluxes from non-local elements
        if(!(influx[0] >= 0.0))
        {
            printf("negative influx=%g\n", influx[0]);
            assert(0);
        }
        
        /*
         if((timeprops_ptr->timesec()>30.0)&&
         (influx[0]>0.0)){
         printf("flux not zeroed 30 seconds\n");
         assert(0);
         }
         */

        //nd = (Node*) NodeTable->lookup(Curr_El->pass_key());
        
        // -- calc contribution of flux source
        flux_src_coef = 0;
        curr_time = (timeprops_ptr->cur_time) * (timeprops_ptr->TIME_SCALE);
        
        /*
         double influx[3];

         influx[0]=*(Curr_El->get_influx()+0);
         influx[1]=*(Curr_El->get_influx()+1);
         influx[2]=*(Curr_El->get_influx()+2);

         if((timeprops_ptr->timesec()>30.0)&&
         (influx[0]>0.0)){
         printf("flux not zeroed 30 seconds\n");
         assert(0);
         }
         */

        //VxVy[2];
        if(Curr_El->state_vars(0) > GEOFLOW_TINY)
        {
            VxVy[0] = Curr_El->state_vars(1) / Curr_El->state_vars(0);
            VxVy[1] = Curr_El->state_vars(2) / Curr_El->state_vars(0);
        }
        else
            VxVy[0] = VxVy[1] = 0.0;
        
#ifdef STOPCRIT_CHANGE_SOURCE
        IF_STOPPED=Curr_El->stoppedflags();
#else
        IF_STOPPED = !(!(Curr_El->stoppedflags()));
#endif
        double gravity[3]{Curr_El->gravity(0),Curr_El->gravity(1),Curr_El->gravity(2)};
        double d_gravity[3]{Curr_El->d_gravity(0),Curr_El->d_gravity(1),Curr_El->d_gravity(2)};
        
        if(elementType == ElementType::TwoPhases)
        {
            //nothing there
            /*predict2ph(Curr_El, d_uvec, (d_uvec + NUM_STATE_VARS), Curr_El->get_prev_state_vars(),
                     tiny, Curr_El->kactxy(0), dt2, gravity, Curr_El->curvature(0), Curr_El->curvature(1),
                     matprops_ptr->bedfrict[Curr_El->material()], matprops_ptr->intfrict,
                     d_gravity, matprops_ptr->frict_tiny, *order_flag, VxVy, IF_STOPPED, influx);*/
        }
        if(elementType == ElementType::SinglePhase)
        {
            /*predict_(Curr_El->get_state_vars(), d_uvec, (d_uvec + NUM_STATE_VARS), Curr_El->get_prev_state_vars(),
                     &tiny, Curr_El->get_kactxy(), &dt2, Curr_El->get_gravity(), Curr_El->get_curvature(),
                     &(matprops_ptr->bedfrict[Curr_El->material()]), &(matprops_ptr->intfrict),
                     Curr_El->get_d_gravity(), &(matprops_ptr->frict_tiny), order_flag, VxVy, &IF_STOPPED, influx);*/
            
            predict(Curr_El, 
                    Curr_El->dh_dx(), Curr_El->dhVx_dx(), Curr_El->dhVy_dx(), 
                    Curr_El->dh_dy(), Curr_El->dhVx_dy(), Curr_El->dhVy_dy(),
                    tiny, Curr_El->kactxy(0), dt2, gravity, Curr_El->curvature(0), Curr_El->curvature(1),
                    matprops_ptr->bedfrict[Curr_El->material()], matprops_ptr->intfrict,
                    d_gravity, matprops_ptr->frict_tiny, *order_flag, VxVy, IF_STOPPED, influx);
        }

        

        /* apply bc's */
#ifdef APPLY_BC
        for(j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(k = 0; k < NUM_STATE_VARS; k++)
                    Curr_El->state_vars(k,0.0);
#endif
        
    }
    titanTimings.predictorStepTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.predictorStepTime += MPI_Wtime() - t_start;
    /* finished predictor step */

    /* really only need to share dudx, state_vars, and kactxy */
    move_data(nump, myid, El_Table, NodeTable, timeprops_ptr);
    
    /* calculate the slopes for the new (half-time step) state variables */
    t_start = MPI_Wtime();
    slopes(El_Table, NodeTable, matprops_ptr);

    // in TWO PHASES #endif  //SECOND_ORDER


    
    /* really only need to share dudx, state_vars, and kactxy */
    move_data(nump, myid, El_Table, NodeTable, timeprops_ptr);
    
    /* calculate kact/pass */
    double dt_not_used = get_coef_and_eigen(elementType, El_Table, NodeTable, matprops_ptr, fluxprops, timeprops_ptr, 1);
    
    /*
     * calculate edge states
     */
    double outflow = 0.0;  //shouldn't need the =0.0 assignment but just being cautious.
    //printf("step: before calc_edge_states\n"); fflush(stdout);
    calc_edge_states(El_Table, NodeTable, matprops_ptr, timeprops_ptr, myid, order_flag, &outflow);
    outflow *= dt;
    //printf("step: after calc_edge_states\n"); fflush(stdout);
    
    /*
     * corrector step and b.c.s
     */
    t_start = MPI_Wtime();
    //for comparison of magnitudes of forces in slumping piles
    double forceint = 0.0, elemforceint;
    double forcebed = 0.0, elemforcebed;
    double eroded = 0.0, elemeroded;
    double deposited = 0.0, elemdeposited;
    double realvolume = 0.0;
    
    // mdj 2007-04 this loop has pretty much defeated me - there is
    //             a dependency in the Element class that causes incorrect
    //             results
    for(i = 0; i < Nelms; i++)
    {
        Curr_El = Elms[i];
        double dx = Curr_El->dx(0);
        double dy = Curr_El->dx(1);
        if(elementType == ElementType::TwoPhases)
        {
            // if calculations are first-order, predict is never called
            // ... so we need to update prev_states
            if(*order_flag == 1)
                Curr_El->update_prev_state_vars();

        }
        void *Curr_El_out = (void *) Curr_El;
        correct(elementType, NodeTable, El_Table, dt, matprops_ptr, fluxprops, timeprops_ptr, Curr_El_out, &elemforceint,
                &elemforcebed, &elemeroded, &elemdeposited);
        
        forceint += fabs(elemforceint);
        forcebed += fabs(elemforcebed);
        realvolume += dx * dy * Curr_El->state_vars(0);
        eroded += elemeroded;
        deposited += elemdeposited;
        
#ifdef APPLY_BC
        for(j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(k = 0; k < NUM_STATE_VARS; k++)
                    Curr_El->state_vars(k, 0.0);
#endif
    }
    t_start2=MPI_Wtime();
    outline_ptr->update(El_Table, NodeTable);
   /* for(i = 0; i < Nelms; i++)
    {
        Curr_El = Elms[i];
        double dx = Curr_El->dx(0);
        double dy = Curr_El->dx(1);
        //update the record of maximum pileheight in the area covered by this element
        double hheight = Curr_El->state_vars(0);
        //if (hheight > 0 && hheight < 0)
        //	;
        
#ifdef MAX_DEPTH_MAP
        double momenta[4];
        if(elementType == ElementType::TwoPhases)
        {
            momenta[0]=Curr_El->state_vars(2);
            momenta[1]=Curr_El->state_vars(3);
            momenta[2]=Curr_El->state_vars(4);
            momenta[3]=Curr_El->state_vars(5);
        }
        if(elementType == ElementType::SinglePhase)
        {
            momenta[0]=Curr_El->state_vars(1);
            momenta[1]=Curr_El->state_vars(2);
        }

        outline_ptr->update(Curr_El->coord(0) - 0.5 * dx, Curr_El->coord(0) + 0.5 * dx, Curr_El->coord(1) - 0.5 * dy,
                            Curr_El->coord(1) + 0.5 * dy, hheight, momenta);

#endif
    }*/
    titanTimings.outlineStepTime += MPI_Wtime() - t_start2;
    titanTimingsAlongSimulation.outlineStepTime += MPI_Wtime() - t_start2;
   /* for(i = 0; i < Nelms; i++)
    {
        Curr_El = Elms[i];
#ifdef APPLY_BC
        for(j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(k = 0; k < NUM_STATE_VARS; k++)
                    Curr_El->state_vars(k, 0.0);
#endif
        
    }*/
    titanTimings.correctorStepTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.correctorStepTime += MPI_Wtime() - t_start;
    //update the orientation of the "dryline" (divides partially wetted cells
    //into wet and dry parts solely based on which neighbors currently have 
    //pileheight greater than GEOFLOW_TINY
    for(i = 0; i < Nelms; i++)
    {
        Elms[i]->calc_wet_dry_orient(El_Table);
    }
    
    /* finished corrector step */

    calc_stats(elementType, El_Table, NodeTable, myid, matprops_ptr, timeprops_ptr, statprops_ptr, discharge, dt);
    
    double tempin[6], tempout[6];
    tempin[0] = outflow;    //volume that flew out the boundaries this iteration
    tempin[1] = eroded;     //volume that was eroded this iteration
    tempin[2] = deposited;  //volume that is currently deposited
    tempin[3] = realvolume; //"actual" volume within boundaries
    tempin[4] = forceint;   //internal friction force
    tempin[5] = forcebed;   //bed friction force
            
    MPI_Reduce(tempin, tempout, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    statprops_ptr->outflowvol += tempout[0] * (matprops_ptr->HEIGHT_SCALE) * (matprops_ptr->LENGTH_SCALE)
                                 * (matprops_ptr->LENGTH_SCALE);
    statprops_ptr->erodedvol += tempout[1] * (matprops_ptr->HEIGHT_SCALE) * (matprops_ptr->LENGTH_SCALE)
                                * (matprops_ptr->LENGTH_SCALE);
    statprops_ptr->depositedvol = tempout[2] * (matprops_ptr->HEIGHT_SCALE) * (matprops_ptr->LENGTH_SCALE)
                                  * (matprops_ptr->LENGTH_SCALE);
    statprops_ptr->realvolume = tempout[3] * (matprops_ptr->HEIGHT_SCALE) * (matprops_ptr->LENGTH_SCALE)
                                * (matprops_ptr->LENGTH_SCALE);
    
    statprops_ptr->forceint = tempout[4] / tempout[3] * matprops_ptr->GRAVITY_SCALE;
    statprops_ptr->forcebed = tempout[5] / tempout[3] * matprops_ptr->GRAVITY_SCALE;
    
#ifdef DEBUG_EXTRA_CHECKING
    El_Table->checkPointersToNeighbours("check_elements_pointers_StepEnd");
    El_Table->ckeckLocalElementsPointers("ckeckAllLocalEntriesPointers_StepEnd");
#endif
    return;
}

/***********************************************************************/
/* calc_volume():                                                      */
/* calculates volume to verify mass conservation                       */
/* determines the maximum height                                       */
/* determines the maximum velocity                                     */
/* determines and returns v_star the non-dimensional stopping velocity */
/***********************************************************************/

void calc_volume(ElementType elementType,ElementsHashTable* El_Table, int myid, MatProps* matprops_ptr, TimeProps* timeprops_ptr, double d_time,
                 double *v_star, double *nz_star)
{
    int i, j, k, counter, imax = 0;
    double tiny = GEOFLOW_TINY;
    //-------------------go through all the elements of the subdomain and  
    //-------------------calculate the state variables at time .5*delta_t
    double volume = 0, volume2 = 0, max_height = 0;
    double v_ave = 0, gl_v_ave;
    double g_ave = 0, gl_g_ave;
    double v_max = 0, gl_v_max;
    double min_height = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
    register double temp;
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(Curr_El->adapted_flag() > 0)
            {
                double dx = Curr_El->dx(0);
                double dy = Curr_El->dx(1);

                if(Curr_El->state_vars(0) > max_height)
                {
                    max_height = Curr_El->state_vars(0);
                    imax = i;
                }

                // rule out non physical fast moving thin layers
                if(Curr_El->state_vars(0) > min_height)
                {
                    if(elementType == ElementType::TwoPhases)
                    {
                        temp = sqrt(Curr_El->state_vars(2) * Curr_El->state_vars(2) + Curr_El->state_vars(3) * Curr_El->state_vars(3));
                        v_ave += temp * dx * dy;
                        temp /= Curr_El->state_vars(1);
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        temp = sqrt(Curr_El->state_vars(1) * Curr_El->state_vars(1) + Curr_El->state_vars(2) * Curr_El->state_vars(2));
                        v_ave += temp * dx * dy;
                        temp /= Curr_El->state_vars(0);
                    }


                    double dvol = Curr_El->state_vars(0) * dx * dy;
                    g_ave += Curr_El->gravity(2) * dvol;
                    volume2 += dvol;

                    if(temp > v_max)
                        v_max = temp;
                }

                /* dxg = c_dmax1(dx, dxg);
                 dyg = c_dmax1(dy, dyg);*/

                volume += Curr_El->state_vars(0) * dx * dy;
            }
        }
    }
    
    double gl_volume = 0, gl_volume2 = 0, gl_max_height;
    double send[5], receive[5] =
    { 0.0, 0.0, 0.0, 0.0, 0.0 };
    send[0] = volume;
    send[1] = volume2;
    send[2] = v_ave;
    send[3] = g_ave;
    i = MPI_Reduce(send, receive, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    gl_volume = receive[0];
    gl_volume2 = receive[1];
    gl_v_ave = receive[2];
    *nz_star = receive[3];
    /*
     i = MPI_Reduce(&volume, &gl_volume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     i = MPI_Reduce(&volume2, &gl_volume2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     i = MPI_Reduce(&v_ave, &gl_v_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     i = MPI_Reduce(&g_ave, &gl_g_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
    send[0] = max_height;
    send[1] = v_max;
    
    i = MPI_Reduce(send, receive, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    gl_max_height = receive[0];
    gl_v_max = receive[1];
    /*  i = MPI_Reduce(&max_height, &gl_max_height, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     i = MPI_Reduce(&v_max, &gl_v_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     */
    if(myid == 0)
    {
        *nz_star = *nz_star / gl_volume2 * matprops_ptr->GRAVITY_SCALE / 9.8;
        //dimensionalize
        gl_v_ave = gl_v_ave / gl_volume2 * sqrt(matprops_ptr->LENGTH_SCALE * matprops_ptr->GRAVITY_SCALE);
        gl_v_max = gl_v_max * sqrt(matprops_ptr->LENGTH_SCALE * (matprops_ptr->GRAVITY_SCALE));
        
        gl_volume = gl_volume * (matprops_ptr->LENGTH_SCALE) * (matprops_ptr->LENGTH_SCALE)
                    * (matprops_ptr->HEIGHT_SCALE);
        gl_max_height = gl_max_height * (matprops_ptr->HEIGHT_SCALE);
        
        d_time *= timeprops_ptr->TIME_SCALE;
        
        /* v_star is the nondimensional global average velocity by v_slump
         once v_slump HAS BEEN CALIBRATED (not yet done see ../main/datread.C) 
         the calculation will terminate when v_star reaches 1 */
        *v_star = gl_v_ave / matprops_ptr->Vslump;
        
        //chunk time
        int hours, minutes;
        double seconds;
        timeprops_ptr->chunktime(&hours, &minutes, &seconds);
        
        printf("At the end of time step %d the time is %d:%02d:%g (hrs:min:sec),\ntime step length is %g [sec], volume is %g [m^3],\nmax height is %g [m], max velocity is %g [m/s],\n ave velocity is %g [m/s], v* = %g\n\n",
               timeprops_ptr->iter, hours, minutes, seconds, d_time, gl_volume, gl_max_height, gl_v_max, gl_v_ave,
               *v_star);
    }
    
    return;
}

/***********************************************************************/
/* the get_max_momentum function was put here because it is similar to */
/* calc volume... which is admittedly not the best reason so if you    */
/* can think of a better place to put it go ahead                      */
/***********************************************************************/

double get_max_momentum(ElementType elementType,ElementsHashTable* El_Table, MatProps* matprops_ptr)
{
    int numprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    double mom2, max_mom = 0, gl_max_mom;
    double min_height = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
    int i;
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
                
            if(EmTemp->adapted_flag() > 0)
            {

                //eliminate fast moving very thin pile from consideration
                if(EmTemp->state_vars(0) >= min_height)
                {
                    if(elementType == ElementType::TwoPhases)
                    {
                        mom2 = (EmTemp->state_vars(2) * EmTemp->state_vars(2) + EmTemp->state_vars(3) * EmTemp->state_vars(3));
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        mom2 = (EmTemp->state_vars(1) * EmTemp->state_vars(1) + EmTemp->state_vars(2) * EmTemp->state_vars(2));
                    }
                    /* mom2 is not a mistake... only need to take the root of 
                     the maximum value */
                    if(mom2 > max_mom)
                        max_mom = mom2;
                }

            }
        }
    }
    
    max_mom = sqrt(max_mom);
    
    if(numprocs > 1)
    {
        if(myid == 0)
            printf("get_max_momentum()1: i=%d\n", i);
        fflush(stdout);
        i = MPI_Reduce(&max_mom, &gl_max_mom, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(myid == 0)
            printf("get_max_momentum()2: i=%d\n", i);
        fflush(stdout);
    }
    else
        gl_max_mom = max_mom;
    
    return (gl_max_mom * matprops_ptr->HEIGHT_SCALE * sqrt(matprops_ptr->LENGTH_SCALE * (matprops_ptr->GRAVITY_SCALE)));
    
}

/**********************************************************************/
/* the sim_end_warning function was put here because it is similar to */
/* calc volume... which is admittedly not the best reason so if you   */
/* can think of a better place to put it go ahead                     */
/**********************************************************************/

void sim_end_warning(ElementType elementType,ElementsHashTable* El_Table, MatProps* matprops_ptr, TimeProps* timeprops_ptr, double v_star)
{
    FILE *fp;
    int myid, numprocs;
    MPI_Status status;
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    /****************************************/
    /* print out the final dimensional time */
    /****************************************/
    //chunk time
    int hours, minutes;
    double seconds;
    timeprops_ptr->chunktime(&hours, &minutes, &seconds);
    
    if(myid == 0)
    {
        //print to screen
        printf("\nTitan2D performed %d time steps before the calculation ended.\n", timeprops_ptr->iter);
        printf("%d:%02d:%g (hrs:min:sec) of time was simulated.\n", hours, minutes, seconds);
        
        //print to file
        fp = fopen("sim_end_warning.readme", "w");
        fprintf(fp, "Titan2D performed %d time steps before the calculation ended.\n", timeprops_ptr->iter);
        fprintf(fp, "%d:%02d:%g (hrs:min:sec) of time was simulated.\n", hours, minutes, seconds);
    }
    
    /*****************************************/
    /* find and print maximum final velocity */
    /*****************************************/
    //double send[2], receive[2];
    double velocity2;
    double v_max = 0;
    double xy_v_max[2];
    double min_height = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT;
    int i;
    struct
    { //for use with MPI_MAXLOC
        double val;
        int rank;
    } send, receive;
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;

    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
            {

                //eliminate fast moving very thin pile from consideration
                if(EmTemp->state_vars(0) >= min_height)
                {
                    if(elementType == ElementType::TwoPhases)
                    {
                        velocity2 = (EmTemp->state_vars(2) * EmTemp->state_vars(2) + EmTemp->state_vars(3) * EmTemp->state_vars(3))
                            / (EmTemp->state_vars(1) * EmTemp->state_vars(1));
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        velocity2 = (EmTemp->state_vars(1) * EmTemp->state_vars(1) + EmTemp->state_vars(2) * EmTemp->state_vars(2))
                            / (EmTemp->state_vars(0) * EmTemp->state_vars(0));
                    }


                    if(velocity2 > v_max)
                    {
                        /* velocity2 is not a mistake... only need to take the root of 
                         the maximum value */
                        v_max = velocity2;

                        xy_v_max[0] = EmTemp->coord(0);
                        xy_v_max[1] = EmTemp->coord(1);
                    }
                }
            }
        }
    }
    
    v_max = sqrt(v_max);
    
    /* get the max value accross all processors */
    send.val = v_max;
    send.rank = myid;
    
    if(numprocs > 1)
    {
        MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        v_max = receive.val;
        
        if(receive.rank != 0)
        { /* don't send location if it's already on the 
         root processor */
            if(receive.rank == myid)
                MPI_Send(xy_v_max, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            else if(myid == 0)
                MPI_Recv(xy_v_max, 2, MPI_DOUBLE, receive.rank, 0, MPI_COMM_WORLD, &status);
        }
    }
    
    // print the rest of the warning
    if(myid == 0)
    {
        //print to screen
        printf("the final v* = v/v_slump = %g\n", v_star);
        printf("The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
               v_max * sqrt(matprops_ptr->LENGTH_SCALE * (matprops_ptr->GRAVITY_SCALE)),
               xy_v_max[0] * matprops_ptr->LENGTH_SCALE, xy_v_max[1] * matprops_ptr->LENGTH_SCALE);
        
        //print to file
        fprintf(fp, "the final v* = v/v_slump = %g\n", v_star);
        fprintf(fp, "The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
                v_max * sqrt(matprops_ptr->LENGTH_SCALE * (matprops_ptr->GRAVITY_SCALE)),
                xy_v_max[0] * matprops_ptr->LENGTH_SCALE, xy_v_max[1] * matprops_ptr->LENGTH_SCALE);
        fclose(fp);
    }
    
    return;
}

