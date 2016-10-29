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
# include <titan_config.h>
#endif

#include "../header/hpfem.h"
#include "../header/geoflow.h"
#include "../header/integrators.h"
#include <math.h>


double Integrator_SinglePhase_Heuristic_Coulomb::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_C(h[ndx], hVx[ndx], hVy[ndx],
                        dh_dx[ndx], dhVx_dx[ndx],
                        dh_dy[ndx], dhVy_dy[ndx],
                        matprops_ptr->bedfrict[material_[ndx]], int_frict,
                        kactxy_[0][ndx], kactxy_[1][ndx], tiny, scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_C(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_SinglePhase_LevelSet_Coulomb::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if((adapted_[ndx]> 0) || (adapted_[ndx]< 0))
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
            if (mindx > 0.0)
            		break;
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_C(h[ndx], hVx[ndx], hVy[ndx],
                        dh_dx[ndx], dhVx_dx[ndx],
                        dh_dy[ndx], dhVy_dy[ndx],
                        matprops_ptr->bedfrict[material_[ndx]], int_frict,
                        kactxy_[0][ndx], kactxy_[1][ndx], tiny, scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_C(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_SinglePhase_Heuristic_Voellmy_Salm::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_VS(kactxy_[0][ndx], kactxy_[1][ndx], scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_VS(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_SinglePhase_LevelSet_Voellmy_Salm::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_VS(kactxy_[0][ndx], kactxy_[1][ndx], scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_VS(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_SinglePhase_Heuristic_Pouliquen_Forterre::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_PF(kactxy_[0][ndx], kactxy_[1][ndx], scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_PF(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_SinglePhase_LevelSet_Pouliquen_Forterre::get_coef_and_eigen(int ghost_flag)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int ierr;
    double min_dx_dy_evalue = 10000000.0, hmax = 0.0;
    //double evalue = 1.0;  //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };

    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * scale_.epsilon) > 0.0)
    {
        double mindx = -1.0;

        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }

        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed

        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;

        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!

            if(h[ndx] > GEOFLOW_TINY)
            {
                double VxVy[2];
                double evalue;

                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];

                gmfggetcoef_PF(kactxy_[0][ndx], kactxy_[1][ndx], scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                VxVy[0] = hVx[ndx] / h[ndx];
                VxVy[1] = hVy[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen_PF(h[ndx], eigenvxymax_[0][ndx],eigenvxymax_[1][ndx],
                        evalue, tiny, kactxy_[0][ndx], gravity_[2][ndx], VxVy);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    printf(" eigenvalue is %e for procd %d momentums are %e %e for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx[ndx], hVy[ndx],
                           h[ndx], maxcurve, coord_[0][ndx],
                           coord_[1][ndx]);
                    assert(0);
                }

                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                stoppedflags_[ndx]=2;
            }
        }
    }

    dt[0] = 0.5 * min_dx_dy_evalue;

    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min

#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];

    return dt[0];
}
double Integrator_TwoPhases::get_coef_and_eigen(int ghost_flag)
{
    MatPropsTwoPhases* matprops2_ptr=static_cast<MatPropsTwoPhases*>(matprops_ptr);
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    int ierr;
    double min_dx_dy_evalue = 10000000, hmax = 0;
    //might need to change this
    //-------------------go through all the elements of the subdomain and get
    //-------------------the coefficients and eigenvalues and calculate the time step
    double global_dt[3], dt[3] = { 0.0, 0.0, HUGE_VAL };
    
    
    //beginning of section that SHOULD ____NOT___ be openmp'd
    //why the first element is good enough?
    double maxinflux = 0.0;
    if((maxinflux = fluxprops_ptr->MaxInfluxNow(matprops_ptr, timeprops_ptr) * (matprops_ptr->scale.epsilon)) > 0.0)
    {
        double mindx = -1.0;
        
        for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        {
            if(adapted_[ndx]!=0)
            {
                mindx = (dx_[0][ndx] < dx_[1][ndx]) ? dx_[0][ndx] : dx_[1][ndx]
                        * pow(0.5, REFINE_LEVEL - generation_[ndx]);
                break;
            }
        }
        
        double dttemp = mindx / maxinflux; //equivalent to dx/eigen_speed
                
        //equivalent to 0.9*sqrt(hmax/g) where hmax=maxinflux*dt when xVel and yVel =0
        double dttemp2 = 0.81 * maxinflux * (scale_.gravity) / 9.8;
        
        dt[2] = c_dmin1(dttemp, dttemp2);
    } //end of section that SHOULD ____NOT___ be openmp'd
    


    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(min:min_dx_dy_evalue) reduction(max:hmax)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if((adapted_[ndx] > 0) || ((adapted_[ndx] < 0) && (ghost_flag == 1)))
        {
            //if this element does not belong on this processor don't involve!!!
            
            if(h[ndx] > GEOFLOW_TINY)
            {
                double Vsolid[2], Vfluid[2];
                double evalue;
                /* calculate hmax */
                if(hmax < h[ndx])
                    hmax = h[ndx];
                

                gmfggetcoef2ph(h_liq[ndx],hVx_sol[ndx],hVy_sol[ndx],
                        dh_liq_dx[ndx],dhVx_sol_dx[ndx],
                        dh_liq_dy[ndx],dhVy_sol_dy[ndx],
                        matprops_ptr->bedfrict[material_[ndx]], int_frict,
                        kactxy_[0][ndx], kactxy_[1][ndx], tiny, scale_.epsilon);

                elements_[ndx].calc_stop_crit(matprops_ptr, this);

                if((stoppedflags_[ndx] < 0) || (stoppedflags_[ndx] > 2))
                    printf("get_coef_and_eigen stopped flag=%d\n", stoppedflags_[ndx]);

                //must use hVx/h and hVy/h rather than eval_velocity (L'Hopital's
                //rule speed if it is smaller) because underestimating speed (which
                //results in over estimating the timestep) is fatal to stability...
                Vsolid[0] = hVx_sol[ndx] / h_liq[ndx];
                Vsolid[1] = hVy_sol[ndx] / h_liq[ndx];

                Vfluid[0] = hVx_liq[ndx] / h[ndx];
                Vfluid[1] = hVy_liq[ndx] / h[ndx];

                //eigen_(EmTemp->eval_state_vars(u_vec_alt),
                eigen2ph(h[ndx], h_liq[ndx], eigenvxymax_[0][ndx],
                          eigenvxymax_[1][ndx], evalue, tiny, kactxy_[0][ndx],
                          gravity_[2][ndx], Vsolid, Vfluid,
                          matprops2_ptr->flow_type);

                // ***********************************************************
                // !!!!!!!!!!!!!!!!!!!!!check dx & dy!!!!!!!!!!!!!!!!!!!!!!!!
                // ***********************************************************
                if(evalue > 1000000000.)
                {
                    double maxcurve = (dabs(curvature_[0][ndx]) > dabs(curvature_[1][ndx]) ? curvature_[0][ndx] : curvature_[1][ndx]);
                    fprintf(stderr,
                            "eigenvalue is %e for procd %d momentums are:\n \
             solid :(%e, %e) \n \
             fluid :(%e, %e) \n \
             for pile height %e curvature=%e (x,y)=(%e,%e)\n",
                            evalue, myid, hVx_sol[ndx],hVy_sol[ndx],
                            hVx_liq[ndx],hVy_liq[ndx],
                            h[ndx], maxcurve, coord_[0][ndx], coord_[1][ndx]);
                    assert(0);
                }
                
                min_dx_dy_evalue = min( min(dx_[0][ndx], dx_[1][ndx]) / evalue, min_dx_dy_evalue);
            }
            else
            {
                elements_[ndx].calc_stop_crit(matprops2_ptr, this); // ensure decent values of kactxy
            }
            
        }
    }
    
    dt[0] = 0.5 * min_dx_dy_evalue;
    
    dt[1] = -0.9 * sqrt(hmax * scale_.epsilon * scale_.gravity / 9.8); //find the negative of the max not the positive min
            
#ifdef USE_MPI
    ierr = MPI_Allreduce(dt, global_dt, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else //USE_MPI
    global_dt[0]=dt[0];
    global_dt[1]=dt[1];
    global_dt[2]=dt[2];
#endif //USE_MPI
    
    dt[0] = 0.5 * c_dmin1(global_dt[0], -global_dt[1]);
    if(dt[0] == 0.0)
        dt[0] = 0.5 * global_dt[2];
    
    
    return dt[0];
}


void find_min_dx(ElementsHashTable* El_Table, double* mindx) {

	int no_of_buckets = El_Table->get_no_of_buckets();
	vector<HashEntryLine> &bucket = El_Table->bucket;
	tivector<Element> &elenode_ = El_Table->elenode_;

	//@ElementsBucketDoubleLoop
	for (int ibuck = 0; ibuck < no_of_buckets; ibuck++) {
		for (int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++) {
			Element* EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);

			if ((EmTemp->adapted_flag() > 0) || (EmTemp->adapted_flag() < 0)) {
				*mindx = (
						((EmTemp->dx(0)) < (EmTemp->dx(1))) ?
								(EmTemp->dx(0)) : (EmTemp->dx(1)))
						* pow(0.5, REFINE_LEVEL - EmTemp->generation());
				break;
			}
		}
		if (*mindx > 0.0)
			break;
	}
}
