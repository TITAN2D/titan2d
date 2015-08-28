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

#include "../header/hpfem.h"
#include "../header/geoflow.h"
#include <math.h>

//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
void gmfggetcoef(const double h,const double hVx,const double hVy, const double *dUdx, const double *dUdy,
        const double bedfrictang, const double intfrictang, 
        double &Kactx, double &Kacty, const double tiny, 
	const double epsilon)
{
    //vel is used to determine if velocity gradients are converging or diverging
    double vel;
    
    //COEFFICIENTS
    double hSQ=h*h;
    double cosphiSQ = cos(intfrictang);
    double tandelSQ = tan(bedfrictang);
    cosphiSQ*=cosphiSQ;
    tandelSQ*=tandelSQ;
    
    
    
    if(h > tiny)
    {
         vel=dUdx[1]/h - hVx*dUdx[0]/hSQ+
             dUdy[2]/h - hVy*dUdy[0]/hSQ;
         Kactx=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
             sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;
         Kacty=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
             sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;

         //if there is no yielding...
         if(fabs(hVx/h) < tiny && fabs(hVy/h) < tiny)
         {
            Kactx = 1.0;
            Kacty = 1.0;
         }
    }
    else
    {
	vel = 0.0;
	Kactx = 1.0;
	Kacty = 1.0;
    }
    Kactx = epsilon * Kactx;
    Kacty = epsilon * Kacty;
}
//! the actual calculation of k active passive is done by a fortran call this should be ripped out and rewritten as a C++ Element member function
void gmfggetcoef2ph(const double h_liq,const double hVx_sol,const double hVy_sol, const double *dUdx, const double *dUdy,
        const double bedfrictang, const double intfrictang, 
        double &Kactx, double &Kacty, const double tiny, 
	const double epsilon)
{
    //vel is used to determine if velocity gradients are converging or diverging
    double vel;
    
    //COEFFICIENTS
    double h_liqSQ=h_liq*h_liq;
    double cosphiSQ = cos(intfrictang);
    double tandelSQ = tan(bedfrictang);
    cosphiSQ*=cosphiSQ;
    tandelSQ*=tandelSQ;
    
    if(h_liq > tiny)
    {
         vel=dUdx[2]/h_liq - hVx_sol*dUdx[1]/h_liqSQ+
             dUdy[3]/h_liq - hVy_sol*dUdy[1]/h_liqSQ;
         Kactx=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
             sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;
         Kacty=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
             sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;

         //if there is no yielding...
         if(fabs(hVx_sol/h_liq) < tiny && fabs(hVy_sol/h_liq) < tiny)
         {
            Kactx = 1.0;
            Kacty = 1.0;
         }
    }
    else
    {
	vel = 0.0;
	Kactx = 1.0;
	Kacty = 1.0;
    }
    Kactx = epsilon * Kactx;
    Kacty = epsilon * Kacty;
}
//! the actual calculation of wave speeds (eigen vectors of the flux jacoboians) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
void eigen( const double h, double &eigenvxmax, double &eigenvymax, double &evalue, 
        const double tiny, double &kactx, const double gravity_z, const double *VxVy) 
{
    if (h > tiny)
    {
        //     iverson and denlinger
        if (kactx < 0.0)
        {
            //negative kactxy
            kactx = -kactx;
        }
        eigenvxmax = fabs(VxVy[0]) + sqrt(kactx * gravity_z * h);
        eigenvymax = fabs(VxVy[1]) + sqrt(kactx * gravity_z * h);

    }
    else
    {
        eigenvxmax = tiny;
        eigenvymax = tiny;
    }

    evalue = c_dmax1(eigenvxmax, eigenvymax);
}
//@TODO is it really h2
void eigen2ph( const double h_sol, const double h_liq, double &eigenvxmax, double &eigenvymax, double &evalue, 
        const double tiny, double &kactx, const double gravity_z, 
        const double *v_solid, const double *v_fluid, const int flowtype) 

{
    double sound_speed;
    if (h_sol > tiny) {
        //iverson and denlinger
        if (kactx < 0.0) {
            kactx = -kactx;
        }

        if (flowtype == 1)
            sound_speed = sqrt(h_sol * kactx * gravity_z);
        else if (flowtype == 2)
            sound_speed = sqrt(h_sol * gravity_z);
        else
            sound_speed = sqrt(h_liq * gravity_z * kactx+(h_sol - h_liq) * gravity_z);

        //x-direction
        eigenvxmax = c_dmax1(fabs(v_solid[0] + sound_speed), fabs(v_fluid[0] + sound_speed));

        //y-direction
        eigenvymax = c_dmax1(fabs(v_solid[1] + sound_speed), fabs(v_fluid[1] + sound_speed));
    }
    else {
        eigenvxmax = tiny;
        eigenvymax = tiny;
    }
    evalue = c_dmax1(eigenvxmax, eigenvymax);
}


double get_coef_and_eigen(ElementType elementType, HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, FluxProps* fluxprops_ptr,
                          TimeProps* timeprops_ptr, int ghost_flag)
{
    MatPropsTwoPhases* matprops2_ptr{nullptr};
    if(elementType == ElementType::TwoPhases)
    {
        matprops2_ptr=static_cast<MatPropsTwoPhases*>(matprops_ptr);
    }
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double min_distance = 1000000, max_evalue = GEOFLOW_TINY, doubleswap;
    
    int ibuck, ierr;
    double tiny = GEOFLOW_TINY, min_dx_dy_evalue = 10000000, hmax = 0;
    double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] =
    { 0.0, 0.0, HUGE_VAL };
    
    HashEntryPtr* elem_bucket_zero = El_Table->getbucketptr();
    HashEntryPtr entryp;
    int num_elem_buckets = El_Table->get_no_of_buckets();
    Element *EmTemp;
    
    //beginning of section that SHOULD ____NOT___ be openmp'd
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * (matprops_ptr->epsilon)) > 0.0)
    {
        double mindx = -1.0;
        ;
        
        for(ibuck = 0; ibuck < num_elem_buckets; ibuck++)
        {
            entryp = *(elem_bucket_zero + ibuck);
            while (entryp)
            {
                EmTemp = (Element*) (entryp->value);
                entryp = entryp->next;
                
                if((EmTemp->adapted_flag() > 0) || (EmTemp->adapted_flag() < 0))
                {
                    mindx = (EmTemp->dx(0) < EmTemp->dx(1)) ?
                                    EmTemp->dx(0) : EmTemp->dx(1)
                            * pow(0.5, REFINE_LEVEL - EmTemp->generation());
                    break;
                }
            }
            if(mindx > 0.0)
                break;
        }
        
        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed
                
        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (matprops_ptr->GRAVITY_SCALE) / 9.8;
        
        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd
    
    double u_vec_alt[3];
    double* d_uvec;
    int intswap;
    double maxcurve;
    int ifanynonzeroheight = 0;
    //for TWO PHASES
    double Vsolid[2], Vfluid[2];
    //for single PHASE
    double VxVy[2];

    for(ibuck = 0; ibuck < num_elem_buckets; ibuck++)
    {
        entryp = *(elem_bucket_zero + ibuck);
        while (entryp)
        {
            EmTemp = (Element*) (entryp->value);
            entryp = entryp->next;
            
            if((EmTemp->adapted_flag() > 0) || ((EmTemp->adapted_flag() < 0) && (ghost_flag == 1)))
            {
                //if this element does not belong on this processor don't involve!!!
                
                if(EmTemp->state_vars(0) > GEOFLOW_TINY)
                {
                    ifanynonzeroheight = 1;
                    
                    /* calculate hmax */
                    if(hmax < EmTemp->state_vars(0))
                        hmax = EmTemp->state_vars(0);
                    
                    d_uvec = EmTemp->get_d_state_vars();

                    if(elementType == ElementType::TwoPhases)
                    {
                        gmfggetcoef2ph(EmTemp->state_vars(1),EmTemp->state_vars(2),EmTemp->state_vars(3), d_uvec, (d_uvec + NUM_STATE_VARS),
                                     matprops_ptr->bedfrict[EmTemp->material()], matprops_ptr->intfrict,
                                     EmTemp->kactxy_ref(0), EmTemp->kactxy_ref(1), tiny, matprops_ptr->epsilon);
                        
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        gmfggetcoef(EmTemp->h(), EmTemp->hVx(), EmTemp->hVy(), d_uvec, (d_uvec + NUM_STATE_VARS),
                                     matprops_ptr->bedfrict[EmTemp->material()], matprops_ptr->intfrict,
                                     EmTemp->kactxy_ref(0), EmTemp->kactxy_ref(1), tiny, matprops_ptr->epsilon);
                        
                        /*gmfggetcoef_(EmTemp->get_state_vars(), d_uvec, (d_uvec + NUM_STATE_VARS), dx_ptr,
                                     &(matprops_ptr->bedfrict[EmTemp->material()]), &(matprops_ptr->intfrict),
                                     EmTemp->get_kactxy(), (EmTemp->get_kactxy() + 1), &tiny, &(matprops_ptr->epsilon));*/
                    }

                    EmTemp->calc_stop_crit(matprops_ptr);
                    intswap = EmTemp->stoppedflags();

                    if((intswap < 0) || (intswap > 2))
                        printf("get_coef_and_eigen stopped flag=%d\n", intswap);
                    
                    //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                    //rule speed if it is smaller) because underestimating speed (which
                    //results in over estimating the timestep) is fatal to stability...
                    if(elementType == ElementType::TwoPhases)
                    {
                        Vsolid[0] = (EmTemp->state_vars(2)) / (EmTemp->state_vars(1));
                        Vsolid[1] = (EmTemp->state_vars(3)) / (EmTemp->state_vars(1));

                        Vfluid[0] = (EmTemp->state_vars(4)) / (EmTemp->state_vars(0));
                        Vfluid[1] = (EmTemp->state_vars(5)) / (EmTemp->state_vars(0));

                        //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                        eigen2ph(EmTemp->h(), EmTemp->h2(), EmTemp->eigenvxymax_ref(0),
                                  EmTemp->eigenvxymax_ref(1), evalue, tiny, EmTemp->kactxy_ref(0),
                                  EmTemp->gravity(2), Vsolid, Vfluid,
                                  matprops2_ptr->flow_type);
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        VxVy[0] = (EmTemp->state_vars(1)) / (EmTemp->state_vars(0));
                        VxVy[1] = (EmTemp->state_vars(2)) / (EmTemp->state_vars(0));

                        //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                        eigen(EmTemp->h(), EmTemp->eigenvxymax_ref(0),
                                  EmTemp->eigenvxymax_ref(1),
                               evalue, tiny, EmTemp->kactxy_ref(0), EmTemp->gravity(2), VxVy);
                    }
                    
                    //printf("evalue=%g\n",evalue);
                    
                    // ***********************************************************
                    // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                    // ***********************************************************
                    doubleswap = c_dmin1(EmTemp->dx(0), EmTemp->dx(1));
                    if(doubleswap / evalue < min_dx_dy_evalue)
                    {
                        min_distance = doubleswap;
                        max_evalue = evalue;
                    }
                    
                    if(evalue > 1000000000.)
                    {
                        maxcurve = (dabs(EmTemp->curvature(0)) > dabs(EmTemp->curvature(1)) ? EmTemp->curvature(0) : EmTemp->curvature(1));
                        if(elementType == ElementType::TwoPhases)
                        {
                            fprintf(stderr,
                                    "eigenvalue is %e for procd %d momentums are:\n \
                     solid :(%e, %e) \n \
                     fluid :(%e, %e) \n \
                     for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                                    evalue, myid, EmTemp->state_vars(2), EmTemp->state_vars(3),
                                    EmTemp->state_vars(4), EmTemp->state_vars(5),
                                    EmTemp->state_vars(0), maxcurve, EmTemp->coord(0),
                                    EmTemp->coord(1));
                            exit(1);
                        }
                        if(elementType == ElementType::SinglePhase)
                        {
                            printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                                   evalue, myid, EmTemp->state_vars(1), EmTemp->state_vars(2),
                                   EmTemp->state_vars(0), maxcurve, EmTemp->coord(0),
                                   EmTemp->coord(1));
                            exit(1);
                        }
                    }
                    
                    min_dx_dy_evalue = c_dmin1( c_dmin1(EmTemp->dx(0), EmTemp->dx(1)) / evalue, min_dx_dy_evalue);
                }
                else
                {
                    if(elementType == ElementType::TwoPhases)
                    {
                        EmTemp->calc_stop_crit(matprops2_ptr); // ensure decent values of kactxy
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        EmTemp->set_stoppedflags(2);
                    }
                }
                
            } //(EmTemp->get_adapted_flag()>0)||...
            
        } //while(entryp)
    }
    
    //if(!ifanynonzeroheight) min_dx_dy_evalue=0.0;
    
    dt[0] = 0.5 * min_dx_dy_evalue;
    //printf("hmax=%g epsilon=%g\n",hmax,matprops_ptr->epsilon);
    dt[1] = -0.9 * sqrt(hmax * (matprops_ptr->epsilon) * (matprops_ptr->GRAVITY_SCALE) / 9.8); //find the negative of the max not the positive min
            
    //printf("myid=%d, dt={%g, %g, %g}\n",myid,dt[0],-dt[1],dt[2]);
    
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    /* if(myid==0)
     printf("myid=0, global_dt={%g, %g, %g}\n",global_dt[0],-global_dt[1],global_dt[2]);
     */

    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];
    
    //printf("====== for proc %d time step should be %e from distance %e and evalue %e  =========\n",
    //	 myid, min_distance/max_evalue, min_distance, max_evalue);
    
    //exit(0);
    
    return dt[0];
}

