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
# include <titan_config.h>
#endif

#include "../header/integrators.h"

#include "../header/hpfem.h"

#include "../header/titan2d_utils.h"
#include "../header/titan_simulation.h"
#include "../header/outline.h"

//#include <advisor-annotate.h>double

Integrator::Integrator(cxxTitanSimulation *_titanSimulation):
    EleNodeRef(_titanSimulation->ElemTable,_titanSimulation->NodeTable),
    timeprops_ptr(_titanSimulation->get_timeprops()),
    matprops_ptr(_titanSimulation->get_matprops()),
    fluxprops_ptr(_titanSimulation->get_fluxprops()),
    pileprops_ptr(_titanSimulation->get_pileprops()),
    statprops_ptr(_titanSimulation->get_statprops()),
    outline_ptr(_titanSimulation->get_outline()),
    discharge_ptr(_titanSimulation->get_discharge_planes()),
	localquants_ptr(_titanSimulation->get_local_quants()),
    elementType(_titanSimulation->get_element_type()),
    adapt(_titanSimulation->adapt),
    TiScalableObject(_titanSimulation->scale_),
    ElemProp(ElemTable, NodeTable)
{

    tiny = GEOFLOW_TINY;

    dt=0.0;

    forceint = 0.0;
    forcebed = 0.0;
    eroded = 0.0;
    deposited = 0.0;
    realvolume = 0.0;

	force_gx = 0.0;
	force_gy = 0.0;
	force_bx = 0.0;
	force_by = 0.0;
	force_bcx = 0.0;
	force_bcy = 0.0;
	force_rx = 0.0;
	force_ry = 0.0;
	power_g = 0.0;
	power_b = 0.0;
	power_bc = 0.0;
	power_r = 0.0;

	Tforce_gx = 0.0;
	Tforce_bx = 0.0;
	Tforce_bcx = 0.0;
	Tforce_rx = 0.0;
	Tforce_gy = 0.0;
	Tforce_by = 0.0;
	Tforce_bcy = 0.0;
	Tforce_ry = 0.0;
	Tpower_g = 0.0;
	Tpower_b = 0.0;
	Tpower_bc = 0.0;
	Tpower_r = 0.0;


    int_frict = 37.0;
    frict_tiny = 0.1;

    order=1;
}
Integrator::~Integrator()
{

}
bool Integrator::scale()
{
    if(TiScalableObject::scale())
    {
        frict_tiny=scale_.frict_tiny;
        int_frict = int_frict * PI / 180.0;
        return true;
    }
    return false;
}
bool Integrator::unscale()
{
    if(TiScalableObject::unscale())
    {
        frict_tiny=scale_.frict_tiny;
        int_frict = int_frict * 180.0/PI;
        return true;
    }
    return false;
}
void Integrator::print0(int spaces)
{
    printf("%*cIntegrator base\n", spaces,' ');
    printf("%*corder: %d\n", spaces+4,' ',order);
    printf("%*cfrict_tiny: %.3f\n", spaces+4,' ',frict_tiny);
}

void Integrator::step()
{
    assert(ElemTable->all_elenodes_are_permanent);
    assert(NodeTable->all_elenodes_are_permanent);

    ASSERT2(ElemTable->check_that_all_elenodes_are_permanent());
    ASSERT2(NodeTable->check_that_all_elenodes_are_permanent());
    ASSERT3(ElemTable->checkPointersToNeighbours("check_elements_pointers_StepStart")==0);

    TIMING1_DEFINE(t_start);
    TIMING1_DEFINE(t_start2);

    PROFILING3_DEFINE(pt_start);
    PROFILING3_START(pt_start);

    //reset class members which can be reset outside
    frict_tiny=scale_.frict_tiny;

    /*
     * PREDICTOR-CORRECTED based on Davis' Simplified Godunov Method
     */

    /* pass off proc data here (really only need state_vars for off-proc neighbors) */
    move_data(numprocs, myid, ElemTable, NodeTable, timeprops_ptr);
    PROFILING3_STOPADD_RESTART(step_other,pt_start);


    TIMING1_START(t_start);
    ElemProp.slopes(matprops_ptr);
    TIMING1_STOPADD(slopesCalcTime, t_start);
    PROFILING3_STOPADD_RESTART(step_slopesCalc,pt_start);

    // get coefficients, eigenvalues, hmax and calculate the time step
    dt = get_coef_and_eigen(0);
    PROFILING3_STOPADD_RESTART(step_get_coef_and_eigen,pt_start);

    timeprops_ptr->incrtime(&dt); //also reduces dt if necessary
    PROFILING3_STOPADD_RESTART(step_other,pt_start);

    // assign influxes and then if any new sources are activating in current time step refine and re-mark cells
    adapt_fluxsrc_region(ElemTable, NodeTable, matprops_ptr, pileprops_ptr, fluxprops_ptr, timeprops_ptr, dt, myid,
                         adapt);
    PROFILING3_STOPADD_RESTART(step_adapt_fluxsrc_region,pt_start);

    int i;
    Element* Curr_El;



    /*
     *  predictor step
     */
    TIMING1_START(t_start);
    PROFILING3_STOPADD_RESTART(step_other,pt_start);

    //-------------------go through all the elements of the subdomain and
    //-------------------calculate the state variables at time .5*delta_t
    if(order!=1)predictor();
    PROFILING3_STOPADD_RESTART(step_predict,pt_start);
    TIMING1_STOPADD(predictorStepTime, t_start);
    /* finished predictor step */

    /* really only need to share dudx, state_vars, and kactxy */
    move_data(numprocs, myid, ElemTable, NodeTable, timeprops_ptr);
    PROFILING3_STOPADD_RESTART(step_other,pt_start);

    /* calculate the slopes for the new (half-time step) state variables */
    TIMING1_START(t_start);
    if(order!=1)ElemProp.slopes(matprops_ptr);
    TIMING1_STOPADD(slopesCalcTime, t_start);
    PROFILING3_STOPADD_RESTART(step_slopesCalc,pt_start);
    // in TWO PHASES #endif  //SECOND_ORDER



    /* really only need to share dudx, state_vars, and kactxy */
    move_data(numprocs, myid, ElemTable, NodeTable, timeprops_ptr);
    PROFILING3_STOPADD_RESTART(step_other,pt_start);

    /* calculate kact/pass */
    if(order!=1||numprocs>1)
        double dt_not_used = get_coef_and_eigen(1);
    PROFILING3_STOPADD_RESTART(step_get_coef_and_eigen,pt_start);
    /*
     * calculate edge states
     */
    double outflow = 0.0;  //shouldn't need the =0.0 assignment but just being cautious.
    //printf("step: before calc_edge_states\n"); fflush(stdout);
    ElemProp.calc_edge_states(matprops_ptr, timeprops_ptr, this, myid, order, outflow);
    PROFILING3_STOPADD_RESTART(step_calc_edge_states,pt_start);

    outflow *= dt;
    //printf("step: after calc_edge_states\n"); fflush(stdout);

    /*
     * corrector step and b.c.s
     */
    TIMING1_START(t_start);
    PROFILING3_STOPADD_RESTART(step_other,pt_start);
    corrector();
    TIMING1_STOPADD(correctorStepTime,t_start);
    PROFILING3_STOPADD_RESTART(step_corrector,pt_start);


    // Recording Local and Global flow characteristics
    flowrecords();

    //statistics, etc.
    TIMING1_START(t_start2);
    if(outline_ptr->enabled)outline_ptr->update();
    TIMING1_STOPADD(outlineStepTime,t_start2);
    PROFILING3_STOPADD_RESTART(step_outline,pt_start);


    //update the orientation of the "dryline" (divides partially wetted cells
    //into wet and dry parts solely based on which neighbors currently have
    //pileheight greater than GEOFLOW_TINY
    ElemProp.calc_wet_dry_orient();
    PROFILING3_STOPADD_RESTART(step_calc_wet_dry_orient,pt_start);

    /* finished corrector step */

    statprops_ptr->calc_stats(myid, matprops_ptr, timeprops_ptr, discharge_ptr, dt);

    double FORCE_SCALE = scale_.gravity * scale_.height * scale_.length * scale_.length;
    double POWER_SCALE = sqrt(scale_.gravity * scale_.length) * FORCE_SCALE;

    double tempin[18], tempout[18];
    tempin[0] = outflow;    //volume that flew out the boundaries this iteration
    tempin[1] = eroded;     //volume that was eroded this iteration
    tempin[2] = deposited;  //volume that is currently deposited
    tempin[3] = realvolume; //"actual" volume within boundaries
    tempin[4] = forceint;   //internal friction force
    tempin[5] = forcebed;   //bed friction force
    tempin[6] = force_gx;
    tempin[7] = force_gy;
    tempin[8] = force_bx;
    tempin[9] = force_by;
    tempin[10] = force_bcx;
    tempin[11] = force_bcy;
    tempin[12] = force_rx;
    tempin[13] = force_ry;
    tempin[14] = power_g;
    tempin[15] = power_b;
    tempin[16] = power_bc;
    tempin[17] = power_r;

#ifdef USE_MPI
    MPI_Reduce(tempin, tempout, 18, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else //USE_MPI
    for(int i=0;i<18;++i)tempout[i]=tempin[i];
#endif //USE_MPI

    statprops_ptr->outflowvol += tempout[0] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length)
                                 * (matprops_ptr->scale.length);
    statprops_ptr->erodedvol += tempout[1] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length)
                                * (matprops_ptr->scale.length);
    statprops_ptr->depositedvol = tempout[2] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length)
                                  * (matprops_ptr->scale.length);
    statprops_ptr->realvolume = tempout[3] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length)
                                * (matprops_ptr->scale.length);

    statprops_ptr->forceint = tempout[4] / tempout[3] * matprops_ptr->scale.gravity;
    statprops_ptr->forcebed = tempout[5] / tempout[3] * matprops_ptr->scale.gravity;

    statprops_ptr->force_gx = tempout[6] * FORCE_SCALE;
    statprops_ptr->force_gy = tempout[7] * FORCE_SCALE;
    statprops_ptr->force_bx = tempout[8] * FORCE_SCALE;
    statprops_ptr->force_by = tempout[9] * FORCE_SCALE;
    statprops_ptr->force_bcx = tempout[10] * FORCE_SCALE;
    statprops_ptr->force_bcy = tempout[11] * FORCE_SCALE;
    statprops_ptr->force_rx = tempout[12] * FORCE_SCALE;
    statprops_ptr->force_ry = tempout[13] * FORCE_SCALE;
    statprops_ptr->power_g = tempout[14] * POWER_SCALE;
    statprops_ptr->power_b = tempout[15] * POWER_SCALE;
    statprops_ptr->power_bc = tempout[16] * POWER_SCALE;
    statprops_ptr->power_r = tempout[17] * POWER_SCALE;

    PROFILING3_STOPADD_RESTART(step_calc_stats,pt_start);

    return;
}
void Integrator::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator","Type");
    TiH5_writeDoubleAttribute(group,  int_frict);
    TiH5_writeDoubleAttribute(group,  frict_tiny);
    TiH5_writeIntAttribute(group,  order);
    TiH5_writeDoubleAttribute(group,  tiny);
}
void Integrator::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));
    TiH5_readDoubleAttribute(group,  int_frict);
    TiH5_readDoubleAttribute(group,  frict_tiny);
    TiH5_readIntAttribute(group,  order);
    TiH5_readDoubleAttribute(group,  tiny);
}
Integrator* Integrator::createIntegrator(const H5::CommonFG *parent, cxxTitanSimulation *_titanSimulation, const  string group_name)
{
    Integrator* integrator=nullptr;
    string integratorType;
    H5::Group group(parent->openGroup(group_name));
    TiH5_readStringAttribute__(group,integratorType,"Type");

    if (integratorType == "Integrator_SinglePhase_Coulomb")
    {
        integrator = new Integrator_SinglePhase_Coulomb(_titanSimulation);
    }
    else if (integratorType == "Integrator_SinglePhase_Voellmy_Salm")
    {
        integrator = new Integrator_SinglePhase_Voellmy_Salm(_titanSimulation);
    }
    else if (integratorType == "Integrator_SinglePhase_Pouliquen_Forterre")
    {
        integrator = new Integrator_SinglePhase_Pouliquen_Forterre(_titanSimulation);
    }
    else if (integratorType == "Integrator_TwoPhases_Coulomb")
    {
        integrator = new Integrator_TwoPhases_Coulomb(_titanSimulation);
    }
    else
    {
        cout << "ERROR: Unknown type of integrator:" << integratorType << "\n";
        assert(integrator != nullptr);
    }
    integrator->h5read(parent,group_name);
    return integrator;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_SinglePhase::Integrator_SinglePhase(cxxTitanSimulation *_titanSimulation):
        Integrator(_titanSimulation),
        h(state_vars_[0]),
        hVx(state_vars_[1]),
        hVy(state_vars_[2]),
        dh_dx(d_state_vars_[0]),
        dh_dy(d_state_vars_[NUM_STATE_VARS]),
        dhVx_dx(d_state_vars_[1]),
        dhVx_dy(d_state_vars_[NUM_STATE_VARS+1]),
        dhVy_dx(d_state_vars_[2]),
        dhVy_dy(d_state_vars_[NUM_STATE_VARS+2])
{
    assert(elementType==ElementType::SinglePhase);

    threshold = 1.0E-02;
    erosion_rate = 0.1;
    do_erosion = 0;
}

void Integrator_SinglePhase::predictor()
{
}

void Integrator_SinglePhase::corrector()
{

}

void Integrator_SinglePhase::flowrecords()
{

}

void Integrator_SinglePhase::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_SinglePhase","Type");
    TiH5_writeDoubleAttribute(group,  threshold);
    TiH5_writeDoubleAttribute(group,  erosion_rate);
    TiH5_writeIntAttribute(group, do_erosion);
}
void Integrator_SinglePhase::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_readDoubleAttribute(group,  threshold);
    TiH5_readDoubleAttribute(group,  erosion_rate);
    TiH5_readIntAttribute(group, do_erosion);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_SinglePhase_Coulomb::Integrator_SinglePhase_Coulomb(cxxTitanSimulation *_titanSimulation):
        Integrator_SinglePhase(_titanSimulation)
{
    assert(elementType==ElementType::SinglePhase);
    assert(order==1);

    stopping_criteria=0;
    thr=0.05;
    //intfrictang=matprops_ptr->intfrict;
    //frict_tiny=matprops_ptr->frict_tiny;
}
bool Integrator_SinglePhase_Coulomb::scale()
{
    if(Integrator_SinglePhase::scale())
    {
        thr = thr/scale_.height;
        return true;
    }
    return false;
}
bool Integrator_SinglePhase_Coulomb::unscale()
{
    if(Integrator_SinglePhase::unscale())
    {
        thr = thr*scale_.height;
        return true;
    }
    return false;
}
void Integrator_SinglePhase_Coulomb::print0(int spaces)
{
    printf("%*cIntegrator: single phase, Coulomb model\n", spaces,' ');
    printf("%*cint_frict:%.3f\n", spaces+4,' ',scaled?int_frict*180.0/PI:int_frict);
    Integrator_SinglePhase::print0(spaces+4);
}

void Integrator_SinglePhase_Coulomb::predictor()
{
    //@TODO OMP me
    if(order==1)return;
    /* mdj 2007-04 */
    int IF_STOPPED;
    double curr_time, influx[3];
    double VxVy[2];
    double dt2 = .5 * dt; // dt2 is set as dt/2 !

    Element* Curr_El;
    //#pragma omp parallel for                                                \
    //private(currentPtr,Curr_El,IF_STOPPED,influx,j,k,curr_time,flux_src_coef,VxVy)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

        Curr_El = &(elements_[ndx]);
        elements_[ndx].update_prev_state_vars();

        influx[0] = Influx_[0][ndx];
        influx[1] = Influx_[1][ndx];
        influx[2] = Influx_[2][ndx];

        //note, now there is no check for fluxes from non-local elements
        if(!(influx[0] >= 0.0))
        {
            printf("negative influx=%g\n", influx[0]);
            assert(0);
        }

        // -- calc contribution of flux source
        curr_time = (timeprops_ptr->cur_time) * (timeprops_ptr->TIME_SCALE);


        //VxVy[2];
        if(h[ndx] > GEOFLOW_TINY)
        {
            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];
        }
        else
        {
            VxVy[0] = VxVy[1] = 0.0;
        }

#ifdef STOPCRIT_CHANGE_SOURCE
        IF_STOPPED=Curr_El->stoppedflags();
#else
        IF_STOPPED = !(!(Curr_El->stoppedflags()));
#endif
        double g[3]{gravity_[0][ndx],gravity_[1][ndx],gravity_[2][ndx]};
        double d_g[3]{d_gravity_[0][ndx],d_gravity_[1][ndx],d_gravity_[2][ndx]};

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //predictor itself

        //NOTE:  d(g[2]*Elm->state_vars(0))/dx is approximated by g[2]*dUvec[0]/dx !!!
        double c_sq;
        double h_inv;
        double tanbed;
        double VxVyS[2];
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

        if (IF_STOPPED == 2) {
            VxVy[0] = 0.0;
            VxVy[1] = 0.0;
            VxVyS[0] = 0.0;
            VxVyS[1] = 0.0;
        }
        else {
            //VxVy[0] = VxVy[0];
            //Elm->state_vars(1)/Elm->state_vars(0);
            //VxVy[1] = VxVy[1];
            //Elm->state_vars(2)/Elm->state_vars(0);
            VxVyS[0] = VxVy[0];
            VxVyS[1] = VxVy[1];
        }


        c_sq = kactxy_[0][ndx] * g[2] * h[ndx];
        //h_inv := 1/h[ndx];

        h[ndx]=h[ndx] - dt2 * (dhVx_dx[ndx] + dhVy_dy[ndx] + influx[0]);
        h[ndx]=c_dmax1(h[ndx], 0.0);

        //dF/dU, dG/dU and S terms if h[ndx] > TINY !
        if (h[ndx] > tiny) {
            h_inv = 1.0 / h[ndx];
            tanbed = tan(matprops_ptr->bedfrict[material_[ndx]]);

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
            hVx[ndx]=hVx[ndx] -
                    dt2 * ((c_sq - VxVy[0] * VxVy[0]) * dh_dx[ndx] +
                    2.0 * VxVy[0] * dhVx_dx[ndx] -
                    VxVy[0] * VxVy[1] * dh_dy[ndx] +
                    VxVy[1] * dhVx_dy[ndx] +
                    VxVy[0] * dhVy_dy[ndx] +
                    influx[1]);

            // x direction source terms

            // the gravity force in the x direction
            forcegrav = g[0] * h[ndx];

            // the internal friction force
            tmp = h_inv * (dhVx_dy[ndx] - VxVyS[0] * dh_dy[ndx]);
            sgn_dudy = sgn_tiny(tmp, frict_tiny);
            forceintx = sgn_dudy * h[ndx] * kactxy_[0][ndx] * (g[2] * dh_dy[ndx] + d_g[1] * h[ndx]) * sin(int_frict);

            // the bed friction force for fast moving flow
            forcebedx = unitvx * c_dmax1(g[2] * h[ndx] + VxVyS[0] * hVx[ndx] * curvature_[0][ndx], 0.0) * tanbed;

            if (IF_STOPPED == 2 && 1 == 0) {
                //the bed friction force for stopped or nearly stopped flow

                //the static friction force is LESS THAN or equal to the friction
                //coefficient times the normal force but it can NEVER exceed the
                //NET force it is opposing

                //maximum friction force the bed friction can support
                forcebedmax = c_dmax1(g[2] * h[ndx] + VxVyS[0] * hVx[ndx] * curvature_[0][ndx], 0.0) * tanbed;

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
            hVx[ndx]=hVx[ndx] + dt2 * (forcegrav - forcebedx - forceintx);
            //            write(*,*) 'int', forceintx, 'bed', forcebedx

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //****** Y-dir ******
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //dF/dU and dG/dU terms
            hVy[ndx]=hVy[ndx] -
                    dt2 * ((c_sq - VxVy[1] * VxVy[1]) * dh_dy[ndx] +
                    2.0 * VxVy[1] * dhVy_dy[ndx] -
                    VxVy[0] * VxVy[1] * dh_dx[ndx] +
                    VxVy[1] * dhVx_dx[ndx] +
                    VxVy[0] * dhVy_dx[ndx] +
                    influx[2]);

            //the gravity force in the y direction
            forcegrav = g[1] * h[ndx];

            //the internal friction force
            tmp = h_inv * (dhVy_dx[ndx] - VxVyS[1] * dh_dx[ndx]);
            sgn_dvdx = sgn_tiny(tmp, frict_tiny);
            forceinty = sgn_dvdx * h[ndx] * kactxy_[0][ndx] * (g[2] * dh_dx[ndx] + d_g[0] * h[ndx]) * sin(int_frict);

            //the bed friction force for fast moving flow
            forcebedy = unitvy *
                    c_dmax1(g[2] * h[ndx] + VxVyS[1] * hVy[ndx] * curvature_[1][ndx], 0.0)
                    * tanbed;

            if (IF_STOPPED == 2 && 1 == 0) {
                //the bed friction force for stopped or nearly stopped flow

                forcebedmax =
                        c_dmax1(g[2] * h[ndx] + VxVyS[1] * hVy[ndx] * curvature_[1][ndx], 0.0)
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
            hVy[ndx]=hVy[ndx] + dt2 * (forcegrav - forcebedy - forceinty);

        }


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // apply bc's
        for(int j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(int k = 0; k < NUM_STATE_VARS; k++)
                    state_vars_[k][ndx]=0.0;
    }
}

void Integrator_SinglePhase_Coulomb::corrector()
{
    //for comparison of magnitudes of forces in slumping piles
    double m_forceint = 0.0;
    double m_forcebed = 0.0;
    double m_eroded = 0.0;
    double m_deposited = 0.0;
    double m_realvolume = 0.0;

    const double sin_intfrictang=sin(int_frict);

    //convinience ref
    tivector<double> *g=gravity_;
    tivector<double> *dgdx=d_gravity_;
    tivector<double> &kactxy=effect_kactxy_[0];
    tivector<double> &bedfrictang=effect_bedfrict_;

    // mdj 2007-04 this loop has pretty much defeated me - there is
    //             a dependency in the Element class that causes incorrect
    //             results
    //ANNOTATE_SITE_BEGIN(ISPC_cor);
    //ANNOTATE_TASK_BEGIN(Integrator_SinglePhase_Coulomb_FirstOrder_corrector_loop);
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_MIDIUM_CHUNK) \
        reduction(+: m_forceint, m_forcebed, m_eroded, m_deposited, m_realvolume)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        //ANNOTATE_ITERATION_TASK(ISPC_cor_iter);
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!
        //if first order states was not updated as there is no predictor
        if(order==1)
        {
            for (int i = 0; i < NUM_STATE_VARS; i++)
                prev_state_vars_[i][ndx]=state_vars_[i][ndx];
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double elem_forceint;
        double elem_forcebed;
        double elem_eroded;
        double elem_deposited;

        double dxdy = dx_[0][ndx] * dx_[1][ndx];
        double dtdx = dt / dx_[0][ndx];
        double dtdy = dt / dx_[1][ndx];

        int xp = positive_x_side_[ndx];
        int yp = (xp + 1) % 4;
        int xm = (xp + 2) % 4;
        int ym = (xp + 3) % 4;

        int ivar, j, k;

        double fluxxp[MAX_NUM_STATE_VARS], fluxyp[MAX_NUM_STATE_VARS];
        double fluxxm[MAX_NUM_STATE_VARS], fluxym[MAX_NUM_STATE_VARS];


        ti_ndx_t nxp = node_key_ndx_[xp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxp[ivar] = node_flux_[ivar][nxp];

        ti_ndx_t nyp = node_key_ndx_[yp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxyp[ivar] = node_flux_[ivar][nyp];

        ti_ndx_t nxm = node_key_ndx_[xm + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxm[ivar] = node_flux_[ivar][nxm];

        ti_ndx_t nym = node_key_ndx_[ym + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxym[ivar] = node_flux_[ivar][nym];


        /* the values being passed to correct are for a SINGLE element, NOT a
         region, as such the only change that having variable bedfriction
         requires is to pass the bedfriction angle for the current element
         rather than the only bedfriction
         I wonder if this is legacy code, it seems odd that it is only called
         for the SUN Operating System zee ../geoflow/correct.f */

#ifdef STOPPED_FLOWS
    #ifdef STOPCRIT_CHANGE_SOURCE
        int IF_STOPPED=stoppedflags_[ndx];
    #else
        int IF_STOPPED = !(!stoppedflags_[ndx]);
    #endif
#endif


        double VxVy[2];
        if(h[ndx] > tiny)
        {
            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];
        }
        else
        {
            VxVy[0] = VxVy[1] = 0.0;
        }

        elements_[ndx].convect_dryline(VxVy[0], VxVy[1], dt); //this is necessary

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //corrector itself

        double speed;
        double forceintx, forceinty;
        double forcebedx, forcebedy;
        double forcebedx_curv, forcebedy_curv;
        double forcebedmax, forcebedequil, forcegravx , forcegravy;
        double unitvx, unitvy;
        double tanbed;
        double Ustore[3];

        double h_inv;
        double sgn_dudy, sgn_dvdx, tmp;
        double es, totalShear;

        Ustore[0] = prev_state_vars_[0][ndx]
                - dtdx * (fluxxp[0] - fluxxm[0])
                - dtdy * (fluxyp[0] - fluxym[0])
                + dt * Influx_[0][ndx];
        Ustore[0] = c_dmax1(Ustore[0], 0.0);

        Ustore[1] = prev_state_vars_[1][ndx]
                - dtdx * (fluxxp[1] - fluxxm[1])
                - dtdy * (fluxyp[1] - fluxym[1])
                + dt * Influx_[1][ndx];

        Ustore[2] = prev_state_vars_[2][ndx]
                - dtdx * (fluxxp[2] - fluxxm[2])
                - dtdy * (fluxyp[2] - fluxym[2])
                + dt * Influx_[2][ndx];

        // initialize to zero
        forceintx = 0.0;
        forcebedx = 0.0;
        forcebedx_curv = 0.0;
        forcebedy_curv = 0.0;
        forceinty = 0.0;
        forcebedy = 0.0;
        unitvx = 0.0;
        unitvy = 0.0;
        elem_eroded = 0.0;

        if(h[ndx] > tiny)
        {
            double inertial_x,inertial_y,drag_x, drag_y;
            // S terms
            // here speed is speed squared
            speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];
            if (speed > 0.0)
            {
                // here speed is speed
                speed = sqrt(speed);
                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }
            tanbed = tan(bedfrictang[ndx]);
            h_inv = 1.0 / h[ndx];

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // x direction source terms
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // the gravity force in the x direction
            forcegravx = g[0][ndx] * h[ndx];

            // the internal friction force
            tmp = h_inv * (dhVx_dy[ndx] - VxVy[0] * dh_dy[ndx]);
            sgn_dudy = sgn_tiny(tmp, frict_tiny);
            forceintx = sgn_dudy * h[ndx]* kactxy[ndx] * (g[2][ndx] * dh_dy[ndx] + dgdx[1][ndx] * h[ndx]) * sin_intfrictang;

            // the bed friction force for fast moving flow
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedx = unitvx * tanbed * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * tanbed * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }
#ifdef STOPPED_FLOWS
            if (IF_STOPPED == 2 && 1 == 0) {
                // the bed friction force for stopped or nearly stopped flow

                // the static friction force is LESS THAN or equal to the friction
                // coefficient times the normal force but it can NEVER exceed the
                // NET force it is opposing

                // maximum friction force the bed friction can support
                forcebedmax = g[2][ndx] * h[ndx] * tanbed;

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
#endif

            tmp = Ustore[1] + dt * (forcegravx - forcebedx - forcebedx_curv - forceintx);
            //STOPPING CRITERIA
            if(stopping_criteria==1)
            {
                inertial_x = fabs(Ustore[1] + dt * forcegravx);
                drag_x = fabs(dt * (forcebedx + forcebedx_curv + forceintx) );

                if (inertial_x <= drag_x)
                    tmp = 0.0;
            }
            Ustore[1] = tmp;

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // y direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // the gravity force in the y direction
            forcegravy = g[1][ndx] * h[ndx];

            // the internal friction force
            tmp = h_inv * (dhVy_dx[ndx] - VxVy[1] * dh_dx[ndx]);
            sgn_dvdx = sgn_tiny(tmp, frict_tiny);
            forceinty = sgn_dvdx * h[ndx] * kactxy[ndx] * (g[2][ndx] * dh_dx[ndx] + dgdx[0][ndx] * h[ndx]) * sin_intfrictang;

            // the bed friction force for fast moving flow
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedy = unitvy * tanbed * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * tanbed * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }
#ifdef STOPPED_FLOWS
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
#endif
            tmp = Ustore[2] + dt * (forcegravy - forcebedy - forcebedy_curv - forceinty);
            //STOPPING CRITERIA
            if(stopping_criteria==1)
            {
                inertial_y = fabs(Ustore[2] + dt * forcegravy);
                drag_y = fabs(dt * (forcebedy + forcebedy_curv + forceinty) );

                if (inertial_y <= drag_y)
                    tmp = 0.0;
            }
            Ustore[2] = tmp;


#ifdef STOPPED_FLOWS
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // (erosion terms) this is Camil's logic, Keith changed some variable
            //names for clarity
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if ((false) && (do_erosion != 0) && (IF_STOPPED == 0)) {
                totalShear = sqrt(forcebedx * forcebedx + forcebedy * forcebedy);
                if ((totalShear > threshold) && (h[ndx] > 0.004)) {

                    es = erosion_rate * sqrt(fabs(totalShear - threshold));
                    elem_eroded = dt*es;
                    Ustore[0] = Ustore[0] + elem_eroded;
                    Ustore[1] = Ustore[1] + elem_eroded * VxVy[0];
                    Ustore[2] = Ustore[2] + elem_eroded * VxVy[1];
                    //write (*,*) 'Doing Keith Erosion Model'
                }
            }
#endif
            if ((do_erosion != 0) && (h[ndx] > threshold)) {
                es = erosion_rate * sqrt(hVx[ndx] * hVx[ndx] + hVy[ndx] * hVy[ndx]) / h[ndx];
                Ustore[0] = Ustore[0] + dt * es;
                Ustore[1] = Ustore[1] + dt * es * Ustore[1];
                Ustore[2] = Ustore[2] + dt * es * Ustore[2];
                //write (*,*) 'Doing Camil Erosion Model'
            }

        }


        // computation of magnitude of friction forces for statistics
        elem_forceint = unitvx * forceintx + unitvy*forceinty;
        elem_forcebed = unitvx * forcebedx + unitvy*forcebedy;

        // update the state variables
        h[ndx]=Ustore[0];
        hVx[ndx]=Ustore[1];
        hVy[ndx]=Ustore[2];

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        elem_forceint *= dxdy;
        elem_forcebed *= dxdy;
        elem_eroded *= dxdy;


        if(stoppedflags_[ndx] == 2)
            elem_deposited = h[ndx] * dxdy;
        else
            elem_deposited = 0.0;

        if(stoppedflags_[ndx])
            elem_eroded = 0.0;

        elements_[ndx].calc_shortspeed(1.0 / dt);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		m_forceint += fabs(elem_forceint);
		m_forcebed += fabs(elem_forcebed);
		m_realvolume += dxdy * h[ndx];
		m_eroded += elem_eroded;
		m_deposited += elem_deposited;

		// apply bc's
		for (int j = 0; j < 4; j++)
			if (neigh_proc_[j][ndx] == INIT)   // this is a boundary!
				for (int k = 0; k < NUM_STATE_VARS; k++)
					state_vars_[k][ndx] = 0.0;
	}

	forceint = m_forceint;
	forcebed = m_forcebed;
	eroded = m_eroded;
	deposited = m_deposited;
	realvolume = m_realvolume;

    //ANNOTATE_TASK_END(Integrator_SinglePhase_Coulomb_FirstOrder_corrector_loop);
    //ANNOTATE_SITE_END(Integrator_SinglePhase_Coulomb_FirstOrder_corrector);
}
void Integrator_SinglePhase_Coulomb::flowrecords()
{
	// Updating spatial derivatives of State Variables
	ElemProp.slopes(matprops_ptr);

	// Temporary reduction variables for recording QoIs globally
    double m_force_gx = 0.0;
    double m_force_gy = 0.0;
    double m_force_bx = 0.0;
    double m_force_by = 0.0;
    double m_force_bcx = 0.0;
    double m_force_bcy = 0.0;
    double m_force_rx = 0.0;
    double m_force_ry = 0.0;
    double m_power_g = 0.0;
    double m_power_b = 0.0;
    double m_power_bc = 0.0;
    double m_power_r = 0.0;

    //convinience ref
    tivector<double> *g=gravity_;
    tivector<double> *dgdx=d_gravity_;
    tivector<double> &bedfrictang=effect_bedfrict_;

    double Kactxy;
    const double sin_intfrictang=sin(int_frict);
    double cosphiSQ = cos(int_frict);
    cosphiSQ*=cosphiSQ;

    // Initializing temporary arrays before searching for locations
    if (localquants_ptr->no_locations > 0)
    	for (int iloc = 0; iloc < localquants_ptr->no_locations; iloc++)
    	{
    		localquants_ptr->temps[iloc].resize(0);
    		localquants_ptr->TimeInts[iloc].resize(0);
    	}

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_MIDIUM_CHUNK) \
        reduction(+: m_force_gx, m_force_gy, m_force_bx, m_force_by, m_force_bcx, m_force_bcy, m_force_rx, m_force_ry) \
		reduction(+: m_power_g, m_power_b, m_power_bc, m_power_r)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

        if(h[ndx] > thr)
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double dxdy = dx_[0][ndx] * dx_[1][ndx];

            double VxVy[2];
            double tanbed = tan(bedfrictang[ndx]);
            double tandelSQ = tanbed * tanbed;

            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];

            //vel is used to determine if velocity gradients are converging or diverging
            double vel;

            //COEFFICIENTS
            double hSQ = h[ndx] * h[ndx];

            vel = dhVx_dx[ndx]/h[ndx] - hVx[ndx]*dh_dx[ndx]/hSQ +
            		dhVy_dy[ndx]/h[ndx] - hVy[ndx]*dh_dy[ndx]/hSQ;

            Kactxy = (2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny) *
            		sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;

            //if there is no yielding...
            if(fabs(VxVy[0]) < tiny && fabs(VxVy[1]) < tiny)
               	Kactxy = 1.0;

            Kactxy *= scale_.epsilon;
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            double speed, h_inv;
            double forceintx, forceinty;
            double forcebedx, forcebedy;
            double forcebedx_curv, forcebedy_curv;
            double forcegravx , forcegravy;
            double unitvx, unitvy;
            double sgn_dudy, sgn_dvdx, tmp;

            // initialize to zero
            forceintx = forcebedx = 0.0;
            forcebedx_curv = forcebedy_curv = 0.0;
            forceinty = forcebedy = 0.0;
            unitvx = unitvy = 0.0;

            // S terms
            // here speed is speed squared
            speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];

            if (speed > 0.0)
            {
                // here speed is speed
                speed = sqrt(speed);
                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = unitvy = 0.0;
            }
            h_inv = 1.0 / h[ndx];

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // x direction source terms
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // the gravity force in the x direction
            forcegravx = g[0][ndx] * h[ndx];

            // the internal friction force
            tmp = h_inv * (dhVx_dy[ndx] - VxVy[0] * dh_dy[ndx]);
            sgn_dudy = sgn_tiny(tmp, frict_tiny);
            forceintx = sgn_dudy * h[ndx] * Kactxy * (g[2][ndx] * dh_dy[ndx] + dgdx[1][ndx] * h[ndx]) * sin_intfrictang;

            // the bed friction force for fast moving flow
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedx = unitvx * tanbed * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * tanbed * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // y direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // the gravity force in the y direction
            forcegravy = g[1][ndx] * h[ndx];

            // the internal friction force
            tmp = h_inv * (dhVy_dx[ndx] - VxVy[1] * dh_dx[ndx]);
            sgn_dvdx = sgn_tiny(tmp, frict_tiny);
            forceinty = sgn_dvdx * h[ndx] * Kactxy * (g[2][ndx] * dh_dx[ndx] + dgdx[0][ndx] * h[ndx]) * sin_intfrictang;

            // the bed friction force for fast moving flow
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedy = unitvy * tanbed * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * tanbed * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }

            /////////////////////////////////////////////////////////////////////////////
            // Recording QoIs globally from updated cells
			m_force_gx += (forcegravx * dxdy);
			m_force_gy += (forcegravy * dxdy);
			m_force_bx -= (forcebedx * dxdy);
			m_force_by -= (forcebedy * dxdy);
			m_force_bcx -= (forcebedx_curv * dxdy);
			m_force_bcy -= (forcebedy_curv * dxdy);
			m_force_rx -= (forceintx * dxdy);
			m_force_ry -= (forceinty * dxdy);

			m_power_g += (forcegravx * VxVy[0] + forcegravy * VxVy[1]) * dxdy;
			m_power_b -= (forcebedx * VxVy[0] + forcebedy * VxVy[1]) * dxdy;
			m_power_bc -= (forcebedx_curv * VxVy[0] + forcebedy_curv * VxVy[1]) * dxdy;
			m_power_r -= (forceintx * VxVy[0] + forceinty * VxVy[1]) * dxdy;

			// Searching user-defined locations to record QoIs
			if (localquants_ptr->no_locations > 0)
				localquants_ptr->FindElement(dt, dx_[0][ndx], dx_[1][ndx],
						coord_[0][ndx], coord_[1][ndx], h[ndx], hVx[ndx], hVy[ndx],
						forcegravx, forcegravy, -forcebedx, -forcebedy,
						-forcebedx_curv, -forcebedy_curv, -forceintx, -forceinty,
						zeta_[0][ndx], zeta_[1][ndx]);
		}
	}

    // Storing QoIs globally from updated cells at current time
	force_gx = m_force_gx;
	force_gy = m_force_gy;
	force_bx = m_force_bx;
	force_by = m_force_by;
	force_bcx = m_force_bcx;
	force_bcy = m_force_bcy;
	force_rx = m_force_rx;
	force_ry = m_force_ry;
	power_g = m_power_g;
	power_b = m_power_b;
	power_bc = m_power_bc;
	power_r = m_power_r;

	// Recording time-integration of QoIs globally from updated cells
	Tforce_gx += m_force_gx * dt;
	Tforce_bx += m_force_bx * dt;
	Tforce_bcx += m_force_bcx * dt;
	Tforce_rx += m_force_rx * dt;
	Tforce_gy += m_force_gy * dt;
	Tforce_by += m_force_by * dt;
	Tforce_bcy += m_force_bcy * dt;
	Tforce_ry += m_force_ry * dt;
	Tpower_g += m_power_g * dt;
	Tpower_b += m_power_b * dt;
	Tpower_bc += m_power_bc * dt;
	Tpower_r += m_power_r * dt;

	// Storing recorded specified QoIs at user-defined locations
    if (localquants_ptr->no_locations > 0)
     	localquants_ptr->StoreQuant(matprops_ptr, timeprops_ptr);
}
void Integrator_SinglePhase_Coulomb::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator_SinglePhase::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_SinglePhase_Coulomb","Type");
}
void Integrator_SinglePhase_Coulomb::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator_SinglePhase::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_SinglePhase_Voellmy_Salm::Integrator_SinglePhase_Voellmy_Salm(cxxTitanSimulation *_titanSimulation):
        Integrator_SinglePhase(_titanSimulation)
{
    assert(elementType==ElementType::SinglePhase);
    assert(order==1);

    mu = 0.5;
    xi = 120.0;
    thr = 0.05;
}
bool Integrator_SinglePhase_Voellmy_Salm::scale()
{
    if(Integrator_SinglePhase::scale())
    {
        xi = xi/scale_.gravity;
        thr = thr/scale_.height;
        return true;
    }
    return false;
}
bool Integrator_SinglePhase_Voellmy_Salm::unscale()
{
    if(Integrator_SinglePhase::unscale())
    {
        xi = xi*scale_.gravity;
        thr = thr*scale_.height;
        return true;
    }
    return false;
}
void Integrator_SinglePhase_Voellmy_Salm::print0(int spaces)
{
    printf("%*cIntegrator: single phase, Voellmy_Salm model, first order\n", spaces,' ');
    printf("%*cmu:%.3f\n", spaces+4,' ',mu);
    printf("%*cxi:%.3f\n", spaces+4,' ',scaled?xi*scale_.gravity:xi);
    Integrator_SinglePhase::print0(spaces+4);
}
void Integrator_SinglePhase_Voellmy_Salm::predictor()
{
}

void Integrator_SinglePhase_Voellmy_Salm::corrector()
{
    //for comparison of magnitudes of forces in slumping piles
    double m_forceint = 0.0;
    double m_forcebed = 0.0;
    double m_eroded = 0.0;
    double m_deposited = 0.0;
    double m_realvolume = 0.0;

    double inv_xi= 1.0/xi;

    //convinience ref
    tivector<double> *g=gravity_;

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) \
        reduction(+: m_forceint, m_forcebed, m_eroded, m_deposited, m_realvolume)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!
        //if first order states was not updated as there is no predictor
        if(order==1)
        {
            for (int i = 0; i < NUM_STATE_VARS; i++)
                prev_state_vars_[i][ndx]=state_vars_[i][ndx];
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double elem_forceint;
        double elem_forcebed;
        double elem_eroded;
        double elem_deposited;

        double dxdy = dx_[0][ndx] * dx_[1][ndx];
        double dtdx = dt / dx_[0][ndx];
        double dtdy = dt / dx_[1][ndx];

        int xp = positive_x_side_[ndx];
        int yp = (xp + 1) % 4;
        int xm = (xp + 2) % 4;
        int ym = (xp + 3) % 4;

        int ivar, j, k;

        double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
        double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];


        ti_ndx_t nxp = node_key_ndx_[xp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxp[ivar] = node_flux_[ivar][nxp];

        ti_ndx_t nyp = node_key_ndx_[yp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxyp[ivar] = node_flux_[ivar][nyp];

        ti_ndx_t nxm = node_key_ndx_[xm + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxm[ivar] = node_flux_[ivar][nxm];

        ti_ndx_t nym = node_key_ndx_[ym + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxym[ivar] = node_flux_[ivar][nym];


        double VxVy[2];
        if(h[ndx] > tiny)
        {
            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];
        }
        else
        {
            VxVy[0] = VxVy[1] = 0.0;
        }

        elements_[ndx].convect_dryline(VxVy[0], VxVy[1], dt); //this is necessary


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //corrector itself
        double speed, speed_squared;
        double forceintx, forceinty;
        double forcebedx, forcebedy;
        double forcebedx_curv, forcebedy_curv;
        double forcegravx,forcegravy;
        double unitvx, unitvy, tmp;
        double Ustore[3];
        double inertial_x, inertial_y;
        double drag_x, drag_y;

        Ustore[0] = prev_state_vars_[0][ndx]
                - dtdx * (fluxxp[0] - fluxxm[0])
                - dtdy * (fluxyp[0] - fluxym[0])
                + dt * Influx_[0][ndx];
        Ustore[0] = c_dmax1(Ustore[0], 0.0);

        Ustore[1] = prev_state_vars_[1][ndx]
                - dtdx * (fluxxp[1] - fluxxm[1])
                - dtdy * (fluxyp[1] - fluxym[1])
                + dt * Influx_[1][ndx];

        Ustore[2] = prev_state_vars_[2][ndx]
                - dtdx * (fluxxp[2] - fluxxm[2])
                - dtdy * (fluxyp[2] - fluxym[2])
                + dt * Influx_[2][ndx];

        // initialize to zero
        forceintx = 0.0;
        forcebedx = 0.0;
        forcebedx_curv = 0.0;
        forcebedy_curv = 0.0;
        forceinty = 0.0;
        forcebedy = 0.0;
        unitvx = 0.0;
        unitvy = 0.0;
        elem_eroded = 0.0;

        if(h[ndx] > tiny)
        {
            // S terms
        	speed_squared = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];

            if (speed_squared > 0.0)
            {

                speed = sqrt(speed_squared);

                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             // x direction source terms
             //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            // The gravity force in the x direction
            forcegravx = g[0][ndx] * h[ndx];

            // The basal type friction force in x direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedx = unitvx * mu * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * mu * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }

            // The velocity-dependent resisting force in x direction
            forceintx = unitvx * speed_squared * inv_xi / scale_.epsilon;

            //STOPPING CRITERIA
            inertial_x = fabs( Ustore[1] + dt * forcegravx );

            drag_x = fabs( dt * ( forceintx + forcebedx + forcebedx_curv) );

            if ( inertial_x > drag_x )
            	Ustore[1] = Ustore[1] + dt * (forcegravx - forcebedx - forcebedx_curv - forceintx);
            else
            	Ustore[1] = 0.0;
             //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             // y direction source terms
             //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            // The gravity force in the y direction
            forcegravy = g[1][ndx] * h[ndx];

            // The basal friction force  in y direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedy = unitvy * mu * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * mu * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }

            // The velocity-dependent resisting force in y direction
            forceinty = unitvy * speed_squared * inv_xi / scale_.epsilon;

            //STOPPING CRITERIA
            inertial_y = fabs( Ustore[2] + dt * forcegravy );

            drag_y = fabs( dt * ( forceinty + forcebedy + forcebedy_curv) );

            if ( inertial_y > drag_y )
            	Ustore[2] = Ustore[2] + dt * (forcegravy - forcebedy - forcebedy_curv - forceinty);
    	    else
    	    	Ustore[2] = 0.0;

        }

        // computation of magnitude of friction forces for statistics
        elem_forceint = speed_squared * inv_xi / scale_.epsilon;
        elem_forcebed = unitvx * forcebedx + unitvy*forcebedy;

        // update the state variables
        h[ndx]=Ustore[0];
        hVx[ndx]=Ustore[1];
        hVy[ndx]=Ustore[2];

        //end of correct
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        elem_forceint *= dxdy;
        elem_forcebed *= dxdy;
        elem_eroded *= dxdy;


        if(stoppedflags_[ndx] == 2)
            elem_deposited = h[ndx] * dxdy;
        else
            elem_deposited = 0.0;

        if(stoppedflags_[ndx])
            elem_eroded = 0.0;

        elements_[ndx].calc_shortspeed(1.0 / dt);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        m_forceint += fabs(elem_forceint);
        m_forcebed += fabs(elem_forcebed);
        m_realvolume += dxdy * h[ndx];
        m_eroded += elem_eroded;
        m_deposited += elem_deposited;

		// apply bc's
		for (int j = 0; j < 4; j++)
			if (neigh_proc_[j][ndx] == INIT)   // this is a boundary!
				for (int k = 0; k < NUM_STATE_VARS; k++)
					state_vars_[k][ndx] = 0.0;
	}

	forceint = m_forceint;
	forcebed = m_forcebed;
	eroded = m_eroded;
	deposited = m_deposited;
	realvolume = m_realvolume;
}
void Integrator_SinglePhase_Voellmy_Salm::flowrecords()
{
	// Updating spatial derivatives of State Variables
	ElemProp.slopes(matprops_ptr);

	// Temporary reduction variables for recording QoIs globally
    double m_force_gx = 0.0;
    double m_force_gy = 0.0;
    double m_force_bx = 0.0;
    double m_force_by = 0.0;
    double m_force_bcx = 0.0;
    double m_force_bcy = 0.0;
    double m_force_rx = 0.0;
    double m_force_ry = 0.0;
    double m_power_g = 0.0;
    double m_power_b = 0.0;
    double m_power_bc = 0.0;
    double m_power_r = 0.0;

    double inv_xi= 1.0/xi;

    //convinience ref
    tivector<double> *g=gravity_;

    // Initializing temporary arrays before searching for locations
    if (localquants_ptr->no_locations > 0)
    	for (int iloc = 0; iloc < localquants_ptr->no_locations; iloc++)
    	{
    		localquants_ptr->temps[iloc].resize(0);
    		localquants_ptr->TimeInts[iloc].resize(0);
    	}

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_MIDIUM_CHUNK) \
        reduction(+: m_force_gx, m_force_gy, m_force_bx, m_force_by, m_force_bcx, m_force_bcy, m_force_rx, m_force_ry) \
		reduction(+: m_power_g, m_power_b, m_power_bc, m_power_r)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

        if(h[ndx] > thr)
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double dxdy = dx_[0][ndx] * dx_[1][ndx];

            double VxVy[2];

            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            double speed, speed_squared, tmp;
            double forceintx, forceinty;
            double forcebedx, forcebedy;
            double forcebedx_curv, forcebedy_curv;
            double forcegravx , forcegravy;
            double unitvx, unitvy;

            // initialize to zero
            forceintx = forcebedx = 0.0;
            forcebedx_curv = forcebedy_curv = 0.0;
            forceinty = forcebedy = 0.0;
            unitvx = unitvy = 0.0;

            // S terms
            // here speed is speed squared
        	speed_squared = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];

            if (speed_squared > 0.0)
            {

                speed = sqrt(speed_squared);

                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // x direction source terms
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            // The gravity force in the x direction
            forcegravx = g[0][ndx] * h[ndx];

            // The basal friction force in x direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0) {
            	forcebedx = unitvx * mu * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * mu * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }

            // The velocity-dependent resisting force in x direction
            forceintx = unitvx * speed_squared * inv_xi / scale_.epsilon;

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // y direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            // The gravity force in the y direction
            forcegravy = g[1][ndx] * h[ndx];

            // The basal friction force in y direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0) {
            	forcebedy = unitvy * mu * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * mu * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }

            // The velocity-dependent resisting force in y direction
            forceinty = unitvy * speed_squared * inv_xi / scale_.epsilon;

            /////////////////////////////////////////////////////////////////////////////
            // Recording QoIs globally from updated cells at current time
			m_force_gx += (forcegravx * dxdy);
			m_force_gy += (forcegravy * dxdy);
			m_force_bx -= (forcebedx * dxdy);
			m_force_by -= (forcebedy * dxdy);
			m_force_bcx -= (forcebedx_curv * dxdy);
			m_force_bcy -= (forcebedy_curv * dxdy);
			m_force_rx -= (forceintx * dxdy);
			m_force_ry -= (forceinty * dxdy);

			m_power_g += (forcegravx * VxVy[0] + forcegravy * VxVy[1]) * dxdy;
			m_power_b -= (forcebedx * VxVy[0] + forcebedy * VxVy[1]) * dxdy;
			m_power_bc -= (forcebedx_curv * VxVy[0] + forcebedy_curv * VxVy[1]) * dxdy;
			m_power_r -= (forceintx * VxVy[0] + forceinty * VxVy[1]) * dxdy;

			// Searching user-defined locations to record QoIs
			if (localquants_ptr->no_locations > 0)
				localquants_ptr->FindElement(dt, dx_[0][ndx], dx_[1][ndx],
						coord_[0][ndx], coord_[1][ndx], h[ndx], hVx[ndx], hVy[ndx],
						forcegravx, forcegravy, -forcebedx, -forcebedy,
						-forcebedx_curv, -forcebedy_curv, -forceintx, -forceinty,
						zeta_[0][ndx], zeta_[1][ndx]);
		}
	}

    // Storing QoIs globally from updated cells at current time
	force_gx = m_force_gx;
	force_gy = m_force_gy;
	force_bx = m_force_bx;
	force_by = m_force_by;
	force_bcx = m_force_bcx;
	force_bcy = m_force_bcy;
	force_rx = m_force_rx;
	force_ry = m_force_ry;
	power_g = m_power_g;
	power_b = m_power_b;
	power_bc = m_power_bc;
	power_r = m_power_r;

	// Recording time-integration of QoIs globally from updated cells
	Tforce_gx += m_force_gx * dt;
	Tforce_bx += m_force_bx * dt;
	Tforce_bcx += m_force_bcx * dt;
	Tforce_rx += m_force_rx * dt;
	Tforce_gy += m_force_gy * dt;
	Tforce_by += m_force_by * dt;
	Tforce_bcy += m_force_bcy * dt;
	Tforce_ry += m_force_ry * dt;
	Tpower_g += m_power_g * dt;
	Tpower_b += m_power_b * dt;
	Tpower_bc += m_power_bc * dt;
	Tpower_r += m_power_r * dt;

	// Storing recorded specified QoIs at user-defined locations
    if (localquants_ptr->no_locations > 0)
     	localquants_ptr->StoreQuant(matprops_ptr, timeprops_ptr);
}
void Integrator_SinglePhase_Voellmy_Salm::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator_SinglePhase::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_SinglePhase_Voellmy_Salm","Type");
    TiH5_writeDoubleAttribute(group, mu);
    TiH5_writeDoubleAttribute(group, xi);
    TiH5_writeDoubleAttribute(group, thr);
}
void Integrator_SinglePhase_Voellmy_Salm::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator_SinglePhase::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_readDoubleAttribute(group, mu);
    TiH5_readDoubleAttribute(group, xi);
    TiH5_readDoubleAttribute(group, thr);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_SinglePhase_Pouliquen_Forterre::Integrator_SinglePhase_Pouliquen_Forterre(cxxTitanSimulation *_titanSimulation):
        Integrator_SinglePhase(_titanSimulation)
{
    assert(elementType==ElementType::SinglePhase);
    assert(order==1);

    phi1=32.9;//in degrees, will convert to rad on scale (0.41887902 rad);
    phi2=42.0;//in degrees, will convert to rad on scale (0.523598776 rad);
    phi3=33.9;
    Beta=0.65;
    L_material=1.0E-3;
    thr=0.05;
}

bool Integrator_SinglePhase_Pouliquen_Forterre::scale()
{
    if(Integrator_SinglePhase::scale())
    {
        phi1*=PI/180.0;
        phi2*=PI/180.0;
        phi3*=PI/180.0;
        L_material = L_material/scale_.height;
        thr = thr/scale_.height;
        return true;
    }
    return false;
}
bool Integrator_SinglePhase_Pouliquen_Forterre::unscale()
{
    if(Integrator_SinglePhase::unscale())
    {
        phi1*=180.0/PI;
        phi2*=180.0/PI;
        phi3*=180.0/PI;
        L_material = L_material*scale_.height;
        thr = thr*scale_.height;
        return true;
    }
    return false;
}
void Integrator_SinglePhase_Pouliquen_Forterre::print0(int spaces)
{
    printf("%*cIntegrator: single phase, Pouliquen_Forterre model, first order\n", spaces,' ');
    printf("%*cphi1:%.3f\n", spaces+4,' ',scaled?phi1*180.0/PI:phi1);
    printf("%*cphi2:%.3f\n", spaces+4,' ',scaled?phi2*180.0/PI:phi2);
    printf("%*cphi3:%.3f\n", spaces+4,' ',scaled?phi3*180.0/PI:phi3);
    printf("%*cBeta:%.3f\n", spaces+4,' ',Beta);
    printf("%*cL_material:%.3e\n", spaces+4,' ',scaled?L_material*scale_.height:L_material);
    Integrator_SinglePhase::print0(spaces+4);
}

void Integrator_SinglePhase_Pouliquen_Forterre::predictor()
{
}

void Integrator_SinglePhase_Pouliquen_Forterre::corrector()
{
    Element* Curr_El;

    //for comparison of magnitudes of forces in slumping piles
    double m_forceint = 0.0;
    double m_forcebed = 0.0;
    double m_eroded = 0.0;
    double m_deposited = 0.0;
    double m_realvolume = 0.0;

    double mu_1 = tan(phi1);
    double mu_2 = tan(phi2);
    double mu_3 = tan(phi3);

    //convinience ref
    tivector<double> *g=gravity_;

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) \
        reduction(+: m_forceint, m_forcebed, m_eroded, m_deposited, m_realvolume)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!
        //if first order states was not updated as there is no predictor
        if(order==1)
        {
            for (int i = 0; i < NUM_STATE_VARS; i++)
                prev_state_vars_[i][ndx]=state_vars_[i][ndx];
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double elem_forceint;
        double elem_forcebed;
        double elem_eroded;
        double elem_deposited;

        double dxdy = dx_[0][ndx] * dx_[1][ndx];
        double dtdx = dt / dx_[0][ndx];
        double dtdy = dt / dx_[1][ndx];

        int xp = positive_x_side_[ndx];
        int yp = (xp + 1) % 4;
        int xm = (xp + 2) % 4;
        int ym = (xp + 3) % 4;

        int ivar, j, k;

        double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
        double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];


        ti_ndx_t nxp = node_key_ndx_[xp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxp[ivar] = node_flux_[ivar][nxp];

        ti_ndx_t nyp = node_key_ndx_[yp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxyp[ivar] = node_flux_[ivar][nyp];

        ti_ndx_t nxm = node_key_ndx_[xm + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxm[ivar] = node_flux_[ivar][nxm];

        ti_ndx_t nym = node_key_ndx_[ym + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxym[ivar] = node_flux_[ivar][nym];


        double VxVy[2];
        if(h[ndx] > tiny)
        {
            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];
        }
        else
        {
            VxVy[0] = VxVy[1] = 0.0;
        }

        elements_[ndx].convect_dryline(VxVy[0], VxVy[1], dt); //this is necessary

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //corrector itself

        double speed;
        double forceintx, forceinty;
        double forcebedx, forcebedy;
        double forcebedx_curv, forcebedy_curv;
        double forcegravx, forcegravy;
        double unitvx, unitvy, tmp;
        double Ustore[3];
        double mu_bed, Local_Fr;

        Ustore[0] = prev_state_vars_[0][ndx]
                - dtdx * (fluxxp[0] - fluxxm[0])
                - dtdy * (fluxyp[0] - fluxym[0])
                + dt * Influx_[0][ndx];
        Ustore[0] = c_dmax1(Ustore[0], 0.0);

        Ustore[1] = prev_state_vars_[1][ndx]
                - dtdx * (fluxxp[1] - fluxxm[1])
                - dtdy * (fluxyp[1] - fluxym[1])
                + dt * Influx_[1][ndx];

        Ustore[2] = prev_state_vars_[2][ndx]
                - dtdx * (fluxxp[2] - fluxxm[2])
                - dtdy * (fluxyp[2] - fluxym[2])
                + dt * Influx_[2][ndx];

        // initialize to zero
        forceintx = 0.0;
        forcebedx = 0.0;
        forcebedx_curv = 0.0;
        forcebedy_curv = 0.0;
        forceinty = 0.0;
        forcebedy = 0.0;
        unitvx = 0.0;
        unitvy = 0.0;
        elem_eroded = 0.0;
        mu_bed = 0.0;
        
        if(h[ndx] > tiny)
        {
        	double inertial_x,inertial_y,drag_x, drag_y;
            // S terms
            // here speed is speed squared
            speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];
            if (speed > 0.0)
            {
                // here speed is speed
                speed = sqrt(speed);
                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }
            // Calculation of Froude number
            Local_Fr = speed / sqrt( g[2][ndx] * h[ndx] * scale_.epsilon);

            //ccccccccccccccc Calculation of mu_bed(Local_Fr,h) ccccccccccccccccc

            //Dynamic flow regime
			if ( Local_Fr >= Beta )
				mu_bed = mu_1 + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] * Beta / ( L_material * Local_Fr ) );

            //Intermediate flow regime
			else if ( ( Local_Fr < Beta ) && ( Local_Fr > 0.0 ) )
				mu_bed = mu_3 + pow( ( Local_Fr / Beta ), 0.001 ) * ( mu_1 - mu_3 ) + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] / L_material );

            //Static regime
			else if ( Local_Fr == 0.0 )
				mu_bed = mu_3 + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] / L_material);

			//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			// x direction source terms
			//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

			// The gravity force in the x direction
			forcegravx = g[0][ndx] * h[ndx];

			// The basal friction forces in x direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedx = unitvx * mu_bed * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * mu_bed * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }

            // The resisting forces due to flow thickness gradient in x direction
			forceintx = h[ndx] * g[2][ndx] * dh_dx[ndx] * scale_.epsilon;

			//STOPPING CRITERIA
			inertial_x = fabs( Ustore[1] + dt * forcegravx );

			drag_x = fabs( dt * ( forcebedx + forcebedx_curv + forceintx ) );

			if ( inertial_x > drag_x )
				Ustore[1] = Ustore[1] + dt * ( forcegravx - forcebedx - forcebedx_curv - forceintx );
			else
				Ustore[1] = 0.0;

			//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			// y direction source terms
			//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

			// The gravity force in the y direction
			forcegravy = g[1][ndx] * h[ndx];

			// The basal friction forces in y direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedy = unitvy * mu_bed * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * mu_bed * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }

            // The resisting forces due to flow thickness gradient in y direction
			forceinty = h[ndx] * g[2][ndx] * dh_dy[ndx] * scale_.epsilon;

			//STOPPING CRITERIA
			inertial_y = fabs( Ustore[2] + dt * forcegravy );

			drag_y = fabs( dt * ( forcebedy + forcebedy_curv + forceinty ) );

			if ( inertial_y > drag_y )
				Ustore[2] = Ustore[2] + dt * ( forcegravy - forcebedy - forcebedy_curv - forceinty );
			else
				Ustore[2] = 0.0;
        }


        // computation of magnitude of friction forces for statistics
        elem_forceint = unitvx * forceintx + unitvy*forceinty;
        elem_forcebed = unitvx * forcebedx + unitvy*forcebedy;

        // update the state variables
        h[ndx]=Ustore[0];
        hVx[ndx]=Ustore[1];
        hVy[ndx]=Ustore[2];

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        elem_forceint *= dxdy;
        elem_forcebed *= dxdy;
        elem_eroded *= dxdy;


        if(stoppedflags_[ndx] == 2)
            elem_deposited = h[ndx] * dxdy;
        else
            elem_deposited = 0.0;

        if(stoppedflags_[ndx])
            elem_eroded = 0.0;

        elements_[ndx].calc_shortspeed(1.0 / dt);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		m_forceint += fabs(elem_forceint);
		m_forcebed += fabs(elem_forcebed);
		m_realvolume += dxdy * h[ndx];
		m_eroded += elem_eroded;
		m_deposited += elem_deposited;

		// apply bc's
		for (int j = 0; j < 4; j++)
			if (neigh_proc_[j][ndx] == INIT)   // this is a boundary!
				for (int k = 0; k < NUM_STATE_VARS; k++)
					state_vars_[k][ndx] = 0.0;
	}

	forceint = m_forceint;
	forcebed = m_forcebed;
	eroded = m_eroded;
	deposited = m_deposited;
	realvolume = m_realvolume;
}
void Integrator_SinglePhase_Pouliquen_Forterre::flowrecords()
{
	// Updating spatial derivatives of State Variables
	ElemProp.slopes(matprops_ptr);

    // Temporary reduction variables for recording QoIs globally
    double m_force_gx = 0.0;
    double m_force_gy = 0.0;
    double m_force_bx = 0.0;
    double m_force_by = 0.0;
    double m_force_bcx = 0.0;
    double m_force_bcy = 0.0;
    double m_force_rx = 0.0;
    double m_force_ry = 0.0;
    double m_power_g = 0.0;
    double m_power_b = 0.0;
    double m_power_bc = 0.0;
    double m_power_r = 0.0;

    double mu_1 = tan(phi1);
    double mu_2 = tan(phi2);
    double mu_3 = tan(phi3);

    //convinience ref
    tivector<double> *g=gravity_;

    // Initializing temporary arrays before searching for locations
    if (localquants_ptr->no_locations > 0)
    	for (int iloc = 0; iloc < localquants_ptr->no_locations; iloc++)
    	{
    		localquants_ptr->temps[iloc].resize(0);
    		localquants_ptr->TimeInts[iloc].resize(0);
    	}

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_MIDIUM_CHUNK) \
        reduction(+: m_force_gx, m_force_gy, m_force_bx, m_force_by, m_force_bcx, m_force_bcy, m_force_rx, m_force_ry) \
		reduction(+: m_power_g, m_power_b, m_power_bc, m_power_r)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

        if(h[ndx] > thr)
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double dxdy = dx_[0][ndx] * dx_[1][ndx];

            double VxVy[2];

            VxVy[0] = hVx[ndx] / h[ndx];
            VxVy[1] = hVy[ndx] / h[ndx];

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double speed;
            double forceintx, forceinty;
            double forcebedx, forcebedy;
            double forcebedx_curv, forcebedy_curv;
            double forcegravx , forcegravy;
            double unitvx, unitvy, tmp;
            double mu_bed, Local_Fr;

            // initialize to zero
            forceintx = forcebedx = 0.0;
            forcebedx_curv = forcebedy_curv = 0.0;
            forceinty = forcebedy = 0.0;
            unitvx = unitvy = 0.0;

            // S terms
            // here speed is speed squared
        	speed = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];

            if (speed > 0.0)
            {
                speed = sqrt(speed);

                unitvx = VxVy[0] / speed;
                unitvy = VxVy[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }
            Local_Fr = speed / sqrt( g[2][ndx] * h[ndx] * scale_.epsilon);

            //ccccccccccccccc Calculation of mu_bed(Local_Fr,h) ccccccccccccccccc

            //Dynamic flow regime
			if ( Local_Fr >= Beta )
				mu_bed = mu_1 + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] * Beta / ( L_material * Local_Fr ) );

            //Intermediate flow regime
			else if ( ( Local_Fr < Beta ) && ( Local_Fr > 0.0 ) )
				mu_bed = mu_3 + pow( ( Local_Fr / Beta ), 0.001 ) * ( mu_1 - mu_3 ) + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] / L_material );

            //Static regime
			else if ( Local_Fr == 0.0 )
				mu_bed = mu_3 + ( mu_2 - mu_1 ) / ( 1.0 + h[ndx] / L_material);
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // x direction source terms
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

			// The gravity force in the x direction
			forcegravx = g[0][ndx] * h[ndx];

			// The basal friction forces in x direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[0] * hVx[ndx] * curvature_[0][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedx = unitvx * mu_bed * g[2][ndx] * h[ndx];
            	forcebedx_curv = unitvx * mu_bed * VxVy[0] * hVx[ndx] * curvature_[0][ndx];
            }

            // The resisting forces due to flow thickness gradient in x direction
			forceintx = h[ndx] * g[2][ndx] * dh_dx[ndx] * scale_.epsilon;

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // y direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

			// The gravity force in the y direction
			forcegravy = g[1][ndx] * h[ndx];

			// The basal friction forces in y direction
            tmp = c_dmax1(g[2][ndx] * h[ndx] + VxVy[1] * hVy[ndx] * curvature_[1][ndx], 0.0);
            if (tmp > 0.0){
            	forcebedy = unitvy * mu_bed * g[2][ndx] * h[ndx];
            	forcebedy_curv = unitvy * mu_bed * VxVy[1] * hVy[ndx] * curvature_[1][ndx];
            }

            // The resisting forces due to flow thickness gradient in y direction
			forceinty = h[ndx] * g[2][ndx] * dh_dy[ndx] * scale_.epsilon;

            /////////////////////////////////////////////////////////////////////////////
			// Recording QoIs globally from updated cells
			m_force_gx += (forcegravx * dxdy);
			m_force_gy += (forcegravy * dxdy);
			m_force_bx -= (forcebedx * dxdy);
			m_force_by -= (forcebedy * dxdy);
			m_force_bcx -= (forcebedx_curv * dxdy);
			m_force_bcy -= (forcebedy_curv * dxdy);
			m_force_rx -= (forceintx * dxdy);
			m_force_ry -= (forceinty * dxdy);

			m_power_g += (forcegravx * VxVy[0] + forcegravy * VxVy[1]) * dxdy;
			m_power_b -= (forcebedx * VxVy[0] + forcebedy * VxVy[1]) * dxdy;
			m_power_bc -= (forcebedx_curv * VxVy[0] + forcebedy_curv * VxVy[1]) * dxdy;
			m_power_r -= (forceintx * VxVy[0] + forceinty * VxVy[1]) * dxdy;

			// Searching user-defined locations to record QoIs
			if (localquants_ptr->no_locations > 0)
				localquants_ptr->FindElement(dt, dx_[0][ndx], dx_[1][ndx],
						coord_[0][ndx], coord_[1][ndx], h[ndx], hVx[ndx], hVy[ndx],
						forcegravx, forcegravy, -forcebedx, -forcebedy,
						-forcebedx_curv, -forcebedy_curv, -forceintx, -forceinty,
						zeta_[0][ndx], zeta_[1][ndx]);
		}
	}

    // Storing QoIs globally from updated cells at current time
	force_gx = m_force_gx;
	force_gy = m_force_gy;
	force_bx = m_force_bx;
	force_by = m_force_by;
	force_bcx = m_force_bcx;
	force_bcy = m_force_bcy;
	force_rx = m_force_rx;
	force_ry = m_force_ry;
	power_g = m_power_g;
	power_b = m_power_b;
	power_bc = m_power_bc;
	power_r = m_power_r;

	// Recording time-integration of QoIs globally from updated cells
	Tforce_gx += m_force_gx * dt;
	Tforce_bx += m_force_bx * dt;
	Tforce_bcx += m_force_bcx * dt;
	Tforce_rx += m_force_rx * dt;
	Tforce_gy += m_force_gy * dt;
	Tforce_by += m_force_by * dt;
	Tforce_bcy += m_force_bcy * dt;
	Tforce_ry += m_force_ry * dt;
	Tpower_g += m_power_g * dt;
	Tpower_b += m_power_b * dt;
	Tpower_bc += m_power_bc * dt;
	Tpower_r += m_power_r * dt;

	// Storing recorded specified QoIs at user-defined locations
    if (localquants_ptr->no_locations > 0)
     	localquants_ptr->StoreQuant(matprops_ptr, timeprops_ptr);
}
void Integrator_SinglePhase_Pouliquen_Forterre::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator_SinglePhase::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_SinglePhase_Pouliquen_Forterre","Type");
    TiH5_writeDoubleAttribute(group, phi1);
    TiH5_writeDoubleAttribute(group, phi2);
    TiH5_writeDoubleAttribute(group, phi3);
    TiH5_writeDoubleAttribute(group, Beta);
    TiH5_writeDoubleAttribute(group, L_material);
    TiH5_writeDoubleAttribute(group, thr);
}
void Integrator_SinglePhase_Pouliquen_Forterre::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator_SinglePhase::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_readDoubleAttribute(group, phi1);
    TiH5_readDoubleAttribute(group, phi2);
    TiH5_readDoubleAttribute(group, phi3);
    TiH5_readDoubleAttribute(group, Beta);
    TiH5_readDoubleAttribute(group, L_material);
    TiH5_readDoubleAttribute(group, thr);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_TwoPhases::Integrator_TwoPhases(cxxTitanSimulation *_titanSimulation):
        Integrator(_titanSimulation),
        h(state_vars_[0]),
        h_liq(state_vars_[1]),
        hVx_sol(state_vars_[2]),
        hVy_sol(state_vars_[3]),
        hVx_liq(state_vars_[4]),
        hVy_liq(state_vars_[5]),
        dh_dx(d_state_vars_[0]),
        dh_dy(d_state_vars_[NUM_STATE_VARS]),
        dh_liq_dx(d_state_vars_[1]),
        dh_liq_dy(d_state_vars_[NUM_STATE_VARS+1]),
        dhVx_sol_dx(d_state_vars_[2]),
        dhVx_sol_dy(d_state_vars_[NUM_STATE_VARS+2]),
        dhVy_sol_dx(d_state_vars_[3]),
        dhVy_sol_dy(d_state_vars_[NUM_STATE_VARS+3]),
        dhVx_liq_dx(d_state_vars_[4]),
        dhVx_liq_dy(d_state_vars_[NUM_STATE_VARS+4]),
        dhVy_liq_dx(d_state_vars_[5]),
        dhVy_liq_dy(d_state_vars_[NUM_STATE_VARS+5])
{
    assert(elementType==ElementType::TwoPhases);
}

void Integrator_TwoPhases::predictor()
{
}

void Integrator_TwoPhases::corrector()
{

}
void Integrator_TwoPhases::flowrecords()
{

}
void Integrator_TwoPhases::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_TwoPhases","Type");
}
void Integrator_TwoPhases::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Integrator_TwoPhases_Coulomb::Integrator_TwoPhases_Coulomb(cxxTitanSimulation *_titanSimulation):
        Integrator_TwoPhases(_titanSimulation)
{
    assert(elementType==ElementType::TwoPhases);
    assert(order==1);


}
bool Integrator_TwoPhases_Coulomb::scale()
{
    if(Integrator_TwoPhases::scale())
    {
        //anything to scale?
        return true;
    }
    return false;
}
bool Integrator_TwoPhases_Coulomb::unscale()
{
    if(Integrator_TwoPhases::unscale())
    {
        //anything to unscale?
        return true;
    }
    return false;
}
void Integrator_TwoPhases_Coulomb::print0(int spaces)
{
    printf("%*cIntegrator: Two Phases Pitman-Le model, first order\n", spaces,' ');
    Integrator_TwoPhases::print0(spaces+4);
}
void Integrator_TwoPhases_Coulomb::predictor()
{

}
void Integrator_TwoPhases_Coulomb::corrector()
{
    MatPropsTwoPhases* matprops2_ptr=static_cast<MatPropsTwoPhases*>(matprops_ptr);

    double den_solid = matprops2_ptr->den_solid;
    double den_fluid = matprops2_ptr->den_fluid;
    double terminal_vel = matprops2_ptr->v_terminal;

    //for comparison of magnitudes of forces in slumping piles
    double m_forceint = 0.0;
    double m_forcebed = 0.0;
    double m_eroded = 0.0;
    double m_deposited = 0.0;
    double m_realvolume = 0.0;


    //intfrictang=matprops_ptr->intfrict;
    //frict_tiny=matprops_ptr->frict_tiny;
    double sin_intfrictang=sin(int_frict);
    double epsilon=scale_.epsilon;

    //convinience ref
    //tivector<double> *g=gravity_;
    //tivector<double> *dgdx=d_gravity_;
    //tivector<double> &kactxy=effect_kactxy_[0];
    //tivector<double> &bedfrictang=effect_bedfrict_;

    // mdj 2007-04 this loop has pretty much defeated me - there is
    //             a dependency in the Element class that causes incorrect
    //             results
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) \
        reduction(+: m_forceint, m_forcebed, m_eroded, m_deposited, m_realvolume)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!
        //if first order states was not updated as there is no predictor
        if(order==1)
        {
            for (int i = 0; i < NUM_STATE_VARS; i++)
                prev_state_vars_[i][ndx]=state_vars_[i][ndx];
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double elem_forceint;
        double elem_forcebed;
        double elem_eroded;
        double elem_deposited;

        double dxdy = dx_[0][ndx] * dx_[1][ndx];
        double dtdx = dt / dx_[0][ndx];
        double dtdy = dt / dx_[1][ndx];

        int xp = positive_x_side_[ndx];
        int yp = (xp + 1) % 4;
        int xm = (xp + 2) % 4;
        int ym = (xp + 3) % 4;

        int ivar, j, k;

        double fluxxp[NUM_STATE_VARS], fluxyp[NUM_STATE_VARS];
        double fluxxm[NUM_STATE_VARS], fluxym[NUM_STATE_VARS];


        ti_ndx_t nxp = node_key_ndx_[xp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxp[ivar] = node_flux_[ivar][nxp];

        ti_ndx_t nyp = node_key_ndx_[yp + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxyp[ivar] = node_flux_[ivar][nyp];

        ti_ndx_t nxm = node_key_ndx_[xm + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxxm[ivar] = node_flux_[ivar][nxm];

        ti_ndx_t nym = node_key_ndx_[ym + 4][ndx];
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            fluxym[ivar] = node_flux_[ivar][nym];


        /* the values being passed to correct are for a SINGLE element, NOT a
         region, as such the only change that having variable bedfriction
         requires is to pass the bedfriction angle for the current element
         rather than the only bedfriction
         I wonder if this is legacy code, it seems odd that it is only called
         for the SUN Operating System zee ../geoflow/correct.f */

//#ifdef STOPPED_FLOWS
    #ifdef STOPCRIT_CHANGE_SOURCE
        int IF_STOPPED=stoppedflags_[ndx];
    #else
        int IF_STOPPED = !(!stoppedflags_[ndx]);
    #endif
//#endif
        double g[3]{gravity_[0][ndx],gravity_[1][ndx],gravity_[2][ndx]};
        double d_g[3]{d_gravity_[0][ndx],d_gravity_[1][ndx],d_gravity_[2][ndx]};


        int i;
        double kactxy[DIMENSION];
        double bedfrict = effect_bedfrict_[ndx];
        double Vfluid[DIMENSION];
        double volf;


        if(h[ndx] > GEOFLOW_TINY)
        {
            for(i = 0; i < DIMENSION; i++)
                kactxy[i] = effect_kactxy_[i][ndx];

            // fluid velocities
            Vfluid[0] = hVx_liq[ndx] / h[ndx];
            Vfluid[1] = hVy_liq[ndx] / h[ndx];

            // volume fractions
            volf = h_liq[ndx] / h[ndx];
        }
        else
        {
            for(i = 0; i < DIMENSION; i++)
            {
                kactxy[i] = matprops2_ptr->scale.epsilon;
                Vfluid[i] = 0.;
            }
            volf = 1.;
            bedfrict = matprops2_ptr->bedfrict[material_[ndx]];
        }

        double Vsolid[DIMENSION];
        if(h_liq[ndx] > GEOFLOW_TINY)
        {
            Vsolid[0] = hVx_sol[ndx] / h_liq[ndx];
            Vsolid[1] = hVy_sol[ndx] / h_liq[ndx];
        }
        else
        {
            Vsolid[0] = Vsolid[1] = 0.0;
        }

        double V_avg[DIMENSION];
        V_avg[0] = Vsolid[0] * volf + Vfluid[0] * (1. - volf);
        V_avg[1] = Vsolid[1] * volf + Vfluid[1] * (1. - volf);
        elements_[ndx].convect_dryline(V_avg[0],V_avg[1], dt); //this is necessary

        double curv_x=curvature_[0][ndx];
        double curv_y=curvature_[1][ndx];
        double xslope=zeta_[0][ndx];
        double yslope=zeta_[1][ndx];

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //corrector
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
        double es, totalShear;
        double drag[4];

        // initialize to zero
        forceintx = 0.0;
        forcebedx = 0.0;
        forceinty = 0.0;
        forcebedy = 0.0;
        unitvx = 0.0;
        unitvy = 0.0;
        elem_eroded = 0.0;

        slope = sqrt(xslope * xslope + yslope * yslope);
        den_frac = den_fluid / den_solid;
        for (i = 0; i < 6; ++i)
            Ustore[i] = prev_state_vars_[i][ndx] + dt * Influx_[i][ndx]
                    - dtdx * (fluxxp[i] - fluxxm[i])
                    - dtdy * (fluxyp[i] - fluxym[i]);

        if (Ustore[0] > tiny)
        {
            // Source terms ...
            // here speed is speed squared
            speed = Vsolid[0] * Vsolid[0] + Vsolid[1] * Vsolid[1];
            if (speed > 0.0)
            {
                // here speed is speed
                speed = sqrt(speed);
                unitvx = Vsolid[0] / speed;
                unitvy = Vsolid[1] / speed;
            }
            else
            {
                unitvx = 0.0;
                unitvy = 0.0;
            }
            tanbed = tan(bedfrict);
            h_inv = 1.0 / h[ndx];
            hphi_inv = 1.0 / h_liq[ndx];
            alphaxx = kactxy[0];
            alphayy = kactxy[0];
            den_frac = den_fluid / den_solid;
            calc_drag_force(ndx, Vsolid, Vfluid, den_solid, den_fluid,terminal_vel, drag);

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    solid fraction x-direction source terms
            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    alphaxy -- see pitman-le (2005)
            tmp = hphi_inv * (dhVx_sol_dy[ndx] - Vsolid[0] * dh_liq_dy[ndx]);
            sgn_dudy = sgn_tiny(tmp, frict_tiny);
            alphaxy = sgn_dudy * sin_intfrictang * kactxy[0];

            //    alphaxz (includes centrifugal effects)
            alphaxz = -unitvx * tanbed
                    * (1.0 + (Vsolid[0] * Vsolid[0]) * curv_x / g[2]);

            //    evaluate t1
            t1 = (1.0 - den_frac)
                    * (-alphaxx * xslope - alphaxy * yslope + alphaxz)
                    * h_liq[ndx] * g[2];
            //    evaluate t2
            t2 = epsilon * den_frac * h_liq[ndx] * g[2] * dh_dx[ndx];
            //    evaluate t3
            t3 = epsilon * den_frac * h_liq[ndx] * g[2] * xslope;
            //    evaluate t4
            t4 = h_liq[ndx] * g[0];
            //    evaluate drag
            t5 = drag[0];
            //    update Ustore
            Ustore[2] = Ustore[2] + dt * (t1 - t2 - t3 + t4 + t5);

            //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            // solid fraction y-direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    alphaxy -- see pitman-le (2005) for definitions
            tmp = hphi_inv * (dhVy_sol_dx[ndx] - Vsolid[1] * dh_liq_dx[ndx]);
            sgn_dvdx = sgn_tiny(tmp, frict_tiny);
            alphaxy = sgn_dvdx * sin_intfrictang * kactxy[0];

            //    alphayz
            alphayz = -unitvy * tanbed
                    * (1.0 + (Vsolid[1] * Vsolid[1]) * curv_y / g[2]);

            //    evaluate t1
            t1 = (1.0 - den_frac)
                    * (-alphaxy * xslope - alphayy * yslope + alphayz)
                    * h_liq[ndx] * g[2];
            //    evaluate t2
            t2 = epsilon * den_frac * h_liq[ndx] * dh_dy[ndx];
            //    evaluate t3
            t3 = epsilon * den_frac * h_liq[ndx] * g[2] * yslope;
            //    evaluate t4 ( gravity along y-dir )
            t4 = h_liq[ndx] * g[1];
            //    drag term
            t5 = drag[1];
            Ustore[3] = Ustore[3] + dt * (t1 - t2 - t3 + t4 + t5);

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    fluid fraction x-direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    gravity on fluid
            t4 = h[ndx] * g[0];
            //    drag force on fluid
            t5 = drag[2];
            Ustore[4] = Ustore[4] + dt * (t4 - t5);

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    fluid fraction y-direction source terms
            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            //    gravity on fluid
            t4 = h[ndx] * g[1];
            //    drag force on fluid
            t5 = drag[3];
            Ustore[5] = Ustore[5] + dt * (t4 - t5);
        }

        // computation of magnitude of friction forces for statistics
        elem_forceint = unitvx * forceintx + unitvy * forceinty;
        elem_forcebed = unitvx * forcebedx + unitvy * forcebedy;

        // update the state variables
        for (i = 0; i < 6; ++i)
            state_vars_[i][ndx]=Ustore[i];

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        elem_forceint *= dxdy;
        elem_forcebed *= dxdy;
        elem_eroded *= dxdy;


        if(stoppedflags_[ndx] == 2)
            elem_deposited = h[ndx] * dxdy;
        else
            elem_deposited = 0.0;

        if(stoppedflags_[ndx])
            elem_eroded = 0.0;

        elements_[ndx].calc_shortspeed(1.0 / dt);


        //correct(elementType, NodeTable, ElemTable, dt, matprops_ptr, fluxprops_ptr, timeprops_ptr, this, Curr_El_out, &elemforceint,
        //        &elemforcebed, &elemeroded, &elemdeposited);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


        m_forceint += fabs(elem_forceint);
        m_forcebed += fabs(elem_forcebed);
        m_realvolume += dxdy * h[ndx];
        m_eroded += elem_eroded;
        m_deposited += elem_deposited;

        // apply bc's
        for(int j = 0; j < 4; j++)
            if(neigh_proc_[j][ndx] == INIT)   // this is a boundary!
                for(int k = 0; k < NUM_STATE_VARS; k++)
                    state_vars_[k][ndx]=0.0;
    }
    forceint = m_forceint;
    forcebed = m_forcebed;
    eroded = m_eroded;
    deposited = m_deposited;
    realvolume = m_realvolume;
}
void Integrator_TwoPhases_Coulomb::flowrecords()
{

}
void Integrator_TwoPhases_Coulomb::h5write(H5::CommonFG *parent, string group_name) const
{
    Integrator_TwoPhases::h5write(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"Integrator_TwoPhases_Coulomb","Type");
}
void Integrator_TwoPhases_Coulomb::h5read(const H5::CommonFG *parent, const  string group_name)
{
    Integrator_TwoPhases::h5read(parent,group_name);
    H5::Group group(parent->openGroup(group_name));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
    double min_height = matprops_ptr->scale.max_negligible_height;
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

#ifdef USE_MPI
    i = MPI_Reduce(send, receive, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else //USE_MPI
    for(int i=0;i<4;++i)receive[i]=send[i];
#endif //USE_MPI

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
    
#ifdef USE_MPI
    i = MPI_Reduce(send, receive, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else //USE_MPI
    for(int i=0;i<2;++i)receive[i]=send[i];
#endif //USE_MPI
    gl_max_height = receive[0];
    gl_v_max = receive[1];
    /*  i = MPI_Reduce(&max_height, &gl_max_height, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     i = MPI_Reduce(&v_max, &gl_v_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
     */
    if(myid == 0)
    {
        *nz_star = *nz_star / gl_volume2 * matprops_ptr->scale.gravity / 9.8;
        //dimensionalize
        gl_v_ave = gl_v_ave / gl_volume2 * sqrt(matprops_ptr->scale.length * matprops_ptr->scale.gravity);
        gl_v_max = gl_v_max * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity));
        
        gl_volume = gl_volume * (matprops_ptr->scale.length) * (matprops_ptr->scale.length)
                    * (matprops_ptr->scale.height);
        gl_max_height = gl_max_height * (matprops_ptr->scale.height);
        
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
    double min_height = matprops_ptr->scale.max_negligible_height;
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
    
#ifdef USE_MPI
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
#else //USE_MPI
    gl_max_mom = max_mom;
#endif //USE_MPI
    return (gl_max_mom * matprops_ptr->scale.height * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)));
    
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
    IF_MPI(MPI_Status status);
    
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
    double min_height = matprops_ptr->scale.max_negligible_height;
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
#ifdef USE_MPI
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
#endif //USE_MPI
    // print the rest of the warning
    if(myid == 0)
    {
        //print to screen
        printf("the final v* = v/v_slump = %g\n", v_star);
        printf("The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
               v_max * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)),
               xy_v_max[0] * matprops_ptr->scale.length, xy_v_max[1] * matprops_ptr->scale.length);
        
        //print to file
        fprintf(fp, "the final v* = v/v_slump = %g\n", v_star);
        fprintf(fp, "The maximum final velocity of %g [m/s] \noccured at the UTM coordinates (%g,%g)\n",
                v_max * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)),
                xy_v_max[0] * matprops_ptr->scale.length, xy_v_max[1] * matprops_ptr->scale.length);
        fclose(fp);
    }
    
    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 0
//! the actual predictor half timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void predict_1ph_coul_(double *Uvec, double *dUdx, double *dUdy,
        double *Uprev, double *tiny, double *kactxy,
        double *dt2, double *g, double *curv,
        double *bedfrictang, double *int_frict,
        double *dgdx, double *frict_tiny, int *order_flag,
        double *VxVy, int *if_stopped, double *fluxcoef);
//! the actual corrector timestep update (finite difference predictor finite volume corrector) is done by a fortran call, this should be ripped out and rewritten as a C++ Element member function
extern "C" void correct_1ph_coul_(double *Uvec, double *Uprev, double *fluxxp,
        double *fluxyp, double *fluxxm, double *fluxym,
        double *tiny, double *dtdx, double *dtdy, double *dt,
        double *dUdx, double *dUdy, double *xslope,
        double *yslope, double *curv, double *int_frict,
        double *bedfrictang, double *g, double *kactxy,
        double *dgdx, double *frict_tiny, double *forceint,
        double *forcebed, int *do_erosion, double *eroded,
        double *VxVy, int *if_stopped, double *fluxcoef);
void correct_legacy_c_call(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
    FluxProps *fluxprops, TimeProps *timeprops, Element *EmTemp, double *forceint, double *forcebed,
    double *eroded, double *deposited)
{
    double dx[2] = {EmTemp->dx(0),EmTemp->dx(1)};
    double dtdx = dt / dx[0];
    double dtdy = dt / dx[1];
    double tiny = GEOFLOW_TINY;
    int xp = EmTemp->positive_x_side();
    int yp = (xp + 1) % 4, xm = (xp + 2) % 4, ym = (xp + 3) % 4;
    int ivar, j, k;
    double fluxxp[3], fluxyp[3], fluxxm[3], fluxym[3];
    Node** nodes=EmTemp->getNodesPtrs();
    Node* nxp = nodes[xp+4];//(Node*) NodeTable->lookup(EmTemp->getNode() + (xp + 4) * 2);
    for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxxp[ivar] = nxp->flux[ivar];
    Node* nyp = nodes[yp+4];//(Node*) NodeTable->lookup(EmTemp->getNode() + (yp + 4) * 2);
    for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxyp[ivar] = nyp->flux[ivar];
    Node* nxm = nodes[xm+4];//(Node*) NodeTable->lookup(EmTemp->getNode() + (xm + 4) * 2);
    for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxxm[ivar] = nxm->flux[ivar];
    Node* nym = nodes[ym+4];//(Node*) NodeTable->lookup(EmTemp->getNode() + (ym + 4) * 2);
    for (ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        fluxym[ivar] = nym->flux[ivar];
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
    int IF_STOPPED=EmTemp->get_stoppedflags();
#else
    int IF_STOPPED = !(!EmTemp->get_stoppedflags());
    double *state_vars = EmTemp->get_state_vars();
    double *prev_state_vars = EmTemp->get_prev_state_vars();
    double *d_state_vars = EmTemp->get_d_state_vars();
    double *gravity = EmTemp->get_gravity();
    double *d_gravity = EmTemp->get_d_gravity();
    double *zeta = EmTemp->get_zeta();
    double *curvature = EmTemp->get_curvature();
    double effect_bedfrict = EmTemp->get_effect_bedfrict();
    double *effect_kactxy = EmTemp->get_effect_kactxy();
    double *Influx = EmTemp->get_influx();
    double VxVy[2];
    if (state_vars[0] > GEOFLOW_TINY) {
        VxVy[0] = state_vars[1] / state_vars[0];
        VxVy[1] = state_vars[2] / state_vars[0];
    } else {
        VxVy[0] = VxVy[1] = 0.0;
    }
    EmTemp->convect_dryline(VxVy, dt); //this is necessary
    correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt,
        d_state_vars, (d_state_vars + NUM_STATE_VARS), &(zeta[0]), &(zeta[1]), curvature,
        &(matprops_ptr->intfrict), &effect_bedfrict, gravity, effect_kactxy, d_gravity,
        &(matprops_ptr->frict_tiny), forceint, forcebed, &do_erosion, eroded, VxVy, &IF_STOPPED,
        Influx);
    *forceint *= dx[0] * dx[1];
    *forcebed *= dx[0] * dx[1];
    *eroded *= dx[0] * dx[1];
#endif
    if (EmTemp->get_stoppedflags() == 2)
        *deposited = state_vars[0] * dx[0] * dx[1];
    else
        *deposited = 0.0;
    if (EmTemp->get_stoppedflags())
        *eroded = 0.0;
    EmTemp->calc_shortspeed(1.0 / dt);
    return;
}
void Integrator_Legacy::predictor()
{
    /* mdj 2007-04 */
    int IF_STOPPED;
    double curr_time, influx[3];
    double VxVy[2];
    double dt2 = .5 * dt; // dt2 is set as dt/2 !
    Element* Curr_El;
    //#pragma omp parallel for                                                \
    //private(currentPtr,Curr_El,IF_STOPPED,influx,j,k,curr_time,flux_src_coef,VxVy)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        Curr_El = &(elements_[ndx]);
        elements_[ndx].update_prev_state_vars();
        influx[0] = Curr_El->Influx(0);
        influx[1] = Curr_El->Influx(1);
        influx[2] = Curr_El->Influx(2);
        //note, now there is no check for fluxes from non-local elements
        if(!(influx[0] >= 0.0))
        {
            printf("negative influx=%g\n", influx[0]);
            assert(0);
        }
        // -- calc contribution of flux source
        curr_time = (timeprops_ptr->cur_time) * (timeprops_ptr->TIME_SCALE);
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
        }
        if(elementType == ElementType::SinglePhase)
        {
            predict(Curr_El,
                    Curr_El->dh_dx(), Curr_El->dhVx_dx(), Curr_El->dhVy_dx(),
                    Curr_El->dh_dy(), Curr_El->dhVx_dy(), Curr_El->dhVy_dy(),
                    tiny, Curr_El->kactxy(0), dt2, gravity, Curr_El->curvature(0), Curr_El->curvature(1),
                    matprops_ptr->bedfrict[Curr_El->material()], matprops_ptr->intfrict,
                    d_gravity, matprops_ptr->frict_tiny, order, VxVy, IF_STOPPED, influx);
        }
        // apply bc's
        for(int j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(int k = 0; k < NUM_STATE_VARS; k++)
                    Curr_El->state_vars(k,0.0);
    }
}
void Integrator_Legacy::corrector()
{
    Element* Curr_El;
    //for comparison of magnitudes of forces in slumping piles
    forceint = 0.0;
    forcebed = 0.0;
    eroded = 0.0;
    deposited = 0.0;
    realvolume = 0.0;
    double elemforceint;
    double elemforcebed;
    double elemeroded;
    double elemdeposited;
    // mdj 2007-04 this loop has pretty much defeated me - there is
    //             a dependency in the Element class that causes incorrect
    //             results
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        Curr_El = &(elements_[ndx]);
        double dx = Curr_El->dx(0);
        double dy = Curr_El->dx(1);
        //if first order states was not updated as there is no predictor
        if(order==1)
            elements_[ndx].update_prev_state_vars();
        void *Curr_El_out = (void *) Curr_El;
        correct(elementType, NodeTable, ElemTable, dt, matprops_ptr, fluxprops_ptr, timeprops_ptr, Curr_El_out, &elemforceint,
                &elemforcebed, &elemeroded, &elemdeposited);
        forceint += fabs(elemforceint);
        forcebed += fabs(elemforcebed);
        realvolume += dx * dy * Curr_El->state_vars(0);
        eroded += elemeroded;
        deposited += elemdeposited;
        // apply bc's
        for(int j = 0; j < 4; j++)
            if(Curr_El->neigh_proc(j) == INIT)   // this is a boundary!
                for(int k = 0; k < NUM_STATE_VARS; k++)
                    Curr_El->state_vars(k, 0.0);
    }
}
#endif
