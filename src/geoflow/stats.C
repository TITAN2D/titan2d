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
 * $Id: stats.C 232 2012-03-27 00:33:41Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/stats.hpp"
#include "../header/properties.h"

#include <cmath>

//#include <advisor-annotate.h>

/* STAT_VOL_FRAC is only here temporarily until a good value for it
 is found (because changing geoflow.h requires recompiling all of
 titan), it will then be moved to ../header/geoflow.h */
#define STAT_VOL_FRAC 0.95

StatProps::StatProps(ElementsHashTable *_ElemTable, NodeHashTable* _NodeTable):
    EleNodeRef(_ElemTable, _NodeTable)
{
    timereached = -1.0;
    xcen = ycen = xvar = yvar = rmean = area = vmean = vxmean = vymean = slopemean = vstar = 0.0;
    realvolume = statvolume = outflowvol = erodedvol = depositedvol = cutoffheight = 0.0;
    piler = hmax = vmax = forceint = forcebed = 0.0;
    heightifreach = xyifreach[0] = xyifreach[1] = timereached = 0.0;
    xyminmax[0] = xyminmax[1] = xyminmax[2] = xyminmax[3] = hxyminmax = 0.0;
    lhs.refnum = lhs.runid = -1;
    runid = -1;
}
StatProps::~StatProps()
{
}
void StatProps::set(const double edge_height, const double test_height, const double test_location_x, const double test_location_y)
{
    hxyminmax = edge_height;

    heightifreach=test_height;
    if(heightifreach < 0.0)
    {
        xyifreach[0]=test_location_x;
        xyifreach[1]=test_location_y;
    }
    else
    {
        heightifreach = xyifreach[0] = xyifreach[1] = HUGE_VAL;
    }

    //to get rid on uninitiallized memory error in saverun() (restart.C)
    forceint = forcebed = 0.0;
}
void StatProps::scale(const MatProps* matprops_ptr)
{
    if(hxyminmax < 0.0)
    {
        hxyminmax = matprops_ptr->scale.max_negligible_height * 10.0;

    }
    hxyminmax /= matprops_ptr->scale.height;

    if(heightifreach > -1.9)
    {

        //default test height is 10 time the maximum negligible height
        if(heightifreach > -1.1 && heightifreach < -0.9)
            heightifreach = matprops_ptr->scale.max_negligible_height * 10.0;

        heightifreach /= matprops_ptr->scale.height;

        xyifreach[0] /= matprops_ptr->scale.length;
        xyifreach[1] /= matprops_ptr->scale.length;
    }
    else
    {
        heightifreach = xyifreach[0] = xyifreach[1] =
        HUGE_VAL;
    }

    //to get rid on uninitiallized memory error in saverun() (restart.C)
    forceint = forcebed = 0.0;
}

void StatProps::calc_stats(int myid, MatProps* matprops, TimeProps* timeprops,
                DischargePlanes* discharge, double d_time)
{
    assert(ElemTable->all_elenodes_are_permanent);

    const ElementType elementType=ElemTable->elementType_;

    int i, iproc;
    double m_area = 0.0, m_max_height = 0.0;
    //double m_cutoffvolume;
    /* the desired volume of material to take
     statistics from, this is portion is the one
     with the largest heights, for example
     sampling the 95% (by volume) of the pile with
     the largest heights, the fraction is set by
     STAT_VOL_FRAC define statement  */
    double m_cutoffheight = 0.0; /* a pile height criteria equivalent to the
     cutoffvolume, this criteria is found at each
     iteration */
    double m_statvolume=1.0; /* volume where pile thickness >= cutoffheight
     statvolume >= m_cutoffvolume */
    double m_realvolume = 0.0; /* the total volume on this processor or across
     processors */
    //double m_slopevolume = 0.0;
    /* volume used to determine the average slope
     in the direction of velocity,
     slopevolume<=statvolume because can't count
     cells with zero velocity */
    const double testpointheight = heightifreach;
    const double testpointx = xyifreach[0];
    const double testpointy = xyifreach[1];

    double testpointmindist2;
    testpointmindist2 = pow(2.0, 30.0); //HUGE_VAL;
    int testpointreach = 0;

    double testvolume = 0.0;
    //double m_slope_ave = 0.0;
    double m_v_max = 0.0, m_v_ave = 0.0, m_vx_ave = 0.0, m_vy_ave = 0.0;

    double xC = 0.0, yC = 0.0, rC = 0.0, m_piler2 = 0.0;
    double m_xVar = 0.0, m_yVar = 0.0;
    //assume that mean starting location is at (x,y) = (1,1)
    const double m_xCen = 1.2; //0; //1.0/(matprops->scale.length);
    const double m_yCen = 0.3; //0; //1.0/(matprops->scale.length);
    const double min_height = matprops->scale.max_negligible_height;
    int numproc;
    double m_x_min = HUGE;
    double m_x_max = -HUGE;
    double m_y_min = HUGE;
    double m_y_max = -HUGE;

    double resolution = 0;
    Get_max_resolution(&resolution);




    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    //define reference with hints for compiler
    double * RESTRICT h=&(state_vars_[0][0]);
    double * RESTRICT hVx=&(state_vars_[1][0]);
    double * RESTRICT hVy=&(state_vars_[2][0]);

    if(elementType == ElementType::TwoPhases){
        hVx=&(state_vars_[2][0]);
        hVy=&(state_vars_[3][0]);
    }

    TI_ASSUME_ALIGNED(h);
    TI_ASSUME_ALIGNED(hVx);
    TI_ASSUME_ALIGNED(hVy);

    const int N = ElemTable->elenode_.size();

    //initiate m_xyminmax
    for(ti_ndx_t ndx = 0; ndx < N; ndx++)
    {
        if(h[ndx] > min_height && h[ndx] >= hxyminmax && adapted_[ndx] > 0 && myid == myprocess_[ndx])
        {
            m_x_min = coord_[0][ndx];
            m_x_max = coord_[0][ndx];
            m_y_min = coord_[1][ndx];
            m_y_max = coord_[1][ndx];
            break;
        }
    }
    /**************************************************/
    /****** TADA!!!!!!!!! we have finally found  ******/
    /****** this iteration's cut off height (and ******/
    /****** the global volumes too) now we can   ******/
    /****** calculate the rest of the stats in a ******/
    /****** straight forward manner              ******/
    /**************************************************/
    int state_vars_bad_values=0;
    //ANNOTATE_SITE_BEGIN(StatProps_calc_stats);
    //ANNOTATE_TASK_BEGIN(StatProps_calc_stats_loop);
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) \
        reduction(min:m_x_min,m_y_min,testpointmindist2) \
        reduction(max:m_x_max,m_y_max,m_v_max,testpointreach) \
        reduction(max:m_max_height) \
        reduction(+:xC,yC,rC,m_area,m_v_ave,m_vx_ave,m_vy_ave)\
        reduction(+:m_piler2,testvolume,m_xVar,m_yVar) /*m_slope_ave,m_slopevolume*/
    for(ti_ndx_t ndx = 0; ndx < N; ndx++)
    {
        if(adapted_[ndx] > 0 && myid == myprocess_[ndx])
        {
            Element* Curr_El = &(ElemTable->elenode_[ndx]);
            if(discharge->num_planes>0)
            {
                //calculate volume passing through "discharge planes"
                double nodescoord[9][2];
                Node* node;

                for(int inode = 0; inode < 8; inode++)
                {
                    nodescoord[inode][0] = coord_[0][node_key_ndx_[inode][ndx]];
                    nodescoord[inode][1] = coord_[1][node_key_ndx_[inode][ndx]];
                }
                nodescoord[8][0] = coord_[0][ndx];
                nodescoord[8][1] = coord_[1][ndx];

                discharge->update(nodescoord, hVx[ndx], hVy[ndx], d_time);
            }

            // rule out non physical fast moving thin layers
            //if(state_vars_[0][ndx] >= cutoffheight){

            if(h[ndx] > min_height)
            {

                if(h[ndx] > m_max_height)
                    m_max_height = h[ndx];

                if(h[ndx] >= hxyminmax)
                {
                    m_x_min = min(m_x_min,coord_[0][ndx]);
                    m_x_max = max(m_x_max,coord_[0][ndx]);
                    m_y_min = min(m_y_min,coord_[1][ndx]);
                    m_y_max = max(m_y_max,coord_[1][ndx]);
                }

                //to test if pileheight of depth testpointheight
                //has reached the testpoint
                double testpointdist2 = (coord_[0][ndx] - testpointx) * (coord_[0][ndx] - testpointx)
                        + (coord_[1][ndx] - testpointy) * (coord_[1][ndx] - testpointy);

                if(testpointdist2 < testpointmindist2)
                {
                    testpointmindist2 = testpointdist2;
                    testpointreach = h[ndx] >= testpointheight;
                }

                double dA = dx_[0][ndx] * dx_[1][ndx];
                m_area += dA;
                double dVol = state_vars_[0][ndx] * dA;
                testvolume += dVol;
                xC += coord_[0][ndx] * dVol;
                yC += coord_[1][ndx] * dVol;

                m_xVar += coord_[0][ndx] * coord_[0][ndx] * dVol;
                m_yVar += coord_[1][ndx] * coord_[1][ndx] * dVol;
                m_piler2 += (coord_[0][ndx] * coord_[0][ndx] + coord_[1][ndx] * coord_[1][ndx]) * dVol;
                rC += sqrt((coord_[0][ndx] - m_xCen) * (coord_[0][ndx] - m_xCen) + (coord_[1][ndx] - m_yCen) * (coord_[1][ndx] - m_yCen)) * dVol;


                m_v_ave += sqrt(hVx[ndx] * hVx[ndx] + hVy[ndx] * hVy[ndx]) * dA;

                double VxVy[2];
                Curr_El->eval_velocity(0.0, 0.0, VxVy);

                if(elementType == ElementType::TwoPhases)
                {
                    if(std::isnan(m_v_ave)||std::isnan(state_vars_[0][ndx])||std::isnan(state_vars_[1][ndx])||std::isnan(state_vars_[2][ndx])
                        ||std::isnan(state_vars_[3][ndx])||std::isnan(state_vars_[4][ndx])||std::isnan(state_vars_[5][ndx]))
                    {
                        //v_ave is NaN
                        cout<<"calc_stats(): NaN detected in element={"<<ElemTable->key_[ndx]<<"} at iter="<<timeprops->iter<<"\n";
                        printf("prevu={%12.6g,%12.6g,%12.6g,%12.6g,%12.6g,%12.6g}\n",
                               prev_state_vars_[0][ndx], prev_state_vars_[1][ndx],
                               prev_state_vars_[2][ndx], prev_state_vars_[3][ndx],
                               prev_state_vars_[4][ndx], prev_state_vars_[5][ndx]);
                        printf("  u={%12.6g,%12.6g,%12.6g,%12.6g,%12.6g,%12.6g}\n", state_vars_[0][ndx], state_vars_[1][ndx],
                               state_vars_[2][ndx], state_vars_[3][ndx], state_vars_[4][ndx], state_vars_[5][ndx]);
                        printf("prev {Vx_s, Vy_s, Vx_f, Vy_f}={%12.6g,%12.6g,%12.6g,%12.6g}\n",
                               prev_state_vars_[2][ndx] / (prev_state_vars_[1][ndx]),
                               prev_state_vars_[3][ndx] / (prev_state_vars_[1][ndx]),
                               prev_state_vars_[4][ndx] / (prev_state_vars_[0][ndx]),
                               prev_state_vars_[5][ndx] / (prev_state_vars_[0][ndx]));
                        printf("this {Vx_s, Vy_s, Vx_f, Vy_f}={%12.6g,%12.6g,%12.6g,%12.6g}\n",
                               state_vars_[2][ndx] / state_vars_[1][ndx], state_vars_[3][ndx] / state_vars_[1][ndx],
                               state_vars_[4][ndx] / state_vars_[0][ndx], state_vars_[5][ndx] / state_vars_[0][ndx]);
                        ElemBackgroundCheck2(ElemTable, NodeTable, &(ElemTable->elenode_[ndx]), stdout);
                        assert(0);
                    }
                }
                if(elementType == ElementType::SinglePhase)
                {
                    if(std::isnan(m_v_ave)||std::isnan(state_vars_[0][ndx])||std::isnan(state_vars_[1][ndx])||std::isnan(state_vars_[2][ndx]))
                    {
                        //v_ave is NaN

                        cout<<"calc_stats(): NaN detected in element={"<<ElemTable->key_[ndx]<<"} at iter="<<timeprops->iter<<"\n";
                        printf("prevu={%12.6g,%12.6g,%12.6g}\n", prev_state_vars_[0][ndx],
                               prev_state_vars_[1][ndx], prev_state_vars_[2][ndx]);
                        printf("    u={%12.6g,%12.6g,%12.6g}\n", state_vars_[0][ndx], state_vars_[1][ndx], state_vars_[2][ndx]);
                        printf("prev {hVx/h,hVy/h}={%12.6g,%12.6g}\n",
                               prev_state_vars_[1][ndx] / prev_state_vars_[0][ndx],
                               prev_state_vars_[2][ndx] / prev_state_vars_[0][ndx]);
                        printf("this {hVx/h,hVy/h}={%12.6g,%12.6g}\n", state_vars_[1][ndx] / state_vars_[0][ndx],
                               state_vars_[2][ndx] / state_vars_[0][ndx]);
                        ElemBackgroundCheck2(ElemTable, NodeTable, &(ElemTable->elenode_[ndx]), stdout);
                        assert(0);
                    }
                }


                double temp = sqrt(VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1]);
                m_v_max = max(m_v_max,temp);
                m_vx_ave += hVx[ndx] * dA;
                m_vy_ave += hVy[ndx] * dA;



                //these are garbage, Bin Yu wanted them when he was trying to come up
                //with a global stopping criteria (to stop the calculation, not the pile)
                //volume averaged slope in the direction of velocity
                //a negative number means the flow is headed uphill
                /*double xslope = 0, yslope = 0;

                Get_slope(resolution, coord_[0][ndx] * matprops->scale.length,
                          coord_[1][ndx] * matprops->scale.length, xslope, yslope);
                if(temp > GEOFLOW_TINY)
                {
                    m_slope_ave += -(hVx[ndx] * xslope + hVy[ndx] * yslope) * dA / temp;
                    m_slopevolume += dVol;
                }*/
            }
        }
    }
    //ANNOTATE_TASK_END(StatProps_calc_stats_loop);
    //ANNOTATE_SITE_END(StatProps_calc_stats);

    MPI_Reduce(&m_x_min, xyminmax, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_x_max, xyminmax+1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_y_min, xyminmax+2, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_y_max, xyminmax+3, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(myid == 0)
    {
        xyminmax[0] *= matprops->scale.length;
        xyminmax[1] *= matprops->scale.length;
        xyminmax[2] *= matprops->scale.length;
        xyminmax[3] *= matprops->scale.length;
    }

    int inttempout;
    double tempin[14], tempout[14], temp2in[2], temp2out[2];

    //find the minimum distance (squared) to the test point
    MPI_Allreduce(&testpointmindist2, tempout, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    //if this processor isn't the closest to the test point it doesn't count as it's flow reaching the point
    if(tempout[0] < testpointmindist2)
        testpointreach = 0;

    //did the closest point to the test point get reached by the flow?
    MPI_Reduce(&testpointreach, &inttempout, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    testpointreach = inttempout;

    tempin[0] = xC;
    tempin[1] = yC;
    tempin[2] = rC;
    tempin[3] = m_area;
    tempin[4] = m_v_ave;
    tempin[5] = m_vx_ave;
    tempin[6] = m_vy_ave;
    tempin[7] = 0.0;//m_slope_ave;
    tempin[8] = m_piler2;
    tempin[9] = 0.0;//m_slopevolume;
    tempin[10] = testvolume;
    tempin[11] = m_xVar;
    tempin[12] = m_yVar;
    tempin[13] = ElemTable->get_no_of_entries();

    i = MPI_Reduce(tempin, tempout, 14, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    temp2in[0] = m_max_height;
    temp2in[1] = m_v_max;
    i = MPI_Reduce(temp2in, temp2out, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myid == 0)
    {
        if(testpointreach && (timereached < 0.0))
            timereached = timeprops->timesec();

        double VELOCITY_SCALE = sqrt(matprops->scale.length * matprops->scale.gravity);
        //dimensionalize
        xcen = tempout[0] * (matprops->scale.length) / tempout[10];
        ycen = tempout[1] * (matprops->scale.length) / tempout[10];
        xvar = tempout[11] * (matprops->scale.length) * (matprops->scale.length) / tempout[10]
                - (xcen) * (xcen);
        yvar = tempout[12] * (matprops->scale.length) * (matprops->scale.length) / tempout[10]
                - (ycen) * (ycen);
        rmean = tempout[2] * (matprops->scale.length) / tempout[10];
        area = tempout[3] * (matprops->scale.length) * (matprops->scale.length);
        vmean = tempout[4] * VELOCITY_SCALE / tempout[10];
        vxmean = tempout[5] * VELOCITY_SCALE / tempout[10];
        vymean = tempout[6] * VELOCITY_SCALE / tempout[10];

        //slopemean = (tempout[9] > 0) ? tempout[7] / tempout[9] : 0.0;

        realvolume = m_realvolume * (matprops->scale.length) * (matprops->scale.length)
                                * (matprops->scale.height);

        //statvolume is really testvolume which is statvolume if it's not disabled
        statvolume = tempout[10] * (matprops->scale.length) * (matprops->scale.length)
                                * (matprops->scale.height);

        cutoffheight = m_cutoffheight * (matprops->scale.height);//@TODO m_cutoffheight is not init
        testvolume = tempout[10] / m_statvolume;//@TODO m_statvolume is not init

        /* the factor of 3^0.5 is a safety factor, this value was chosen because
         * it makes the "radius" of a uniformly distributed line equal to half
         * the line length
         */
        //3 standard deviations out ~ 99.5% of the material
        piler = 3.0 * sqrt(xvar + yvar);
        hmax = temp2out[0] * (matprops->scale.height);
        vmax = temp2out[1] * VELOCITY_SCALE;

        /* v_star is the nondimensional global average velocity by v_slump
         once v_slump HAS BEEN CALIBRATED (not yet done see ../main/datread.C)
         the calculation will terminate when v_star reaches 1 */
        vstar = vmean / matprops->Vslump;

        /******************/
        /* output section */
        /******************/

        /* output Center Of Mass and x and y components of mean velocity to
         assist the dynamic gis update daemon */
        if(elementType == ElementType::SinglePhase)
        {
            FILE* fp2 = fopen("com.up", "w");
            fprintf(fp2, "%d %g %g %g %g %g %g\n", timeprops->iter, timeprops->timesec(), xcen, ycen,
                    vxmean, vymean, piler);
            fclose(fp2);
        }
        /* standard to screen output */
        d_time *= timeprops->TIME_SCALE;
        //chunk time
        int hours, minutes;
        double seconds;
        timeprops->chunktime(&hours, &minutes, &seconds);

        printf("At the end of time step %d the time is %d:%02d:%g (hrs:min:sec),\n", timeprops->iter, hours, minutes,
               seconds);
        printf("\ttime step length is %g [sec], volume is %g [m^3],\n", d_time, statvolume);
        printf("\tmax height is %g [m], max velocity is %g [m/s],\n", hmax, vmax);
        printf("\tave velocity is %g [m/s], v* = %g,\n", vmean, vstar);
        printf("\ttotal number of elements %.0f\n", tempout[13]);
        printf("\txyminmax %.9e %.9e %.9e %.9e\n", xyminmax[0], xyminmax[1], xyminmax[2], xyminmax[3]);
        printf("\n");
    }

    return;
}

void out_final_stats(TimeProps* timeprops, StatProps* statprops)
{
//round the run time to the nearest second and chunk it
    int walltime = (int) (time(NULL) - timeprops->starttime + 0.5);
    int wallhours = walltime / 3600;
    int wallminutes = (walltime % 3600) / 60;
    int wallseconds = walltime % 60;
    int cputime = (int) (clock() / CLOCKS_PER_SEC + 0.5);
    int cpuhours = cputime / 3600;
    int cpuminutes = (cputime % 3600) / 60;
    int cpuseconds = cputime % 60;

    char filename[256];
    sprintf(filename, "finalstats.%06d", statprops->runid);
//printf("runid=%d\n",statprops->runid);
    FILE* fp = fopen(filename, "w");
//print    runid  xC yC   rC  area maxh    wall time       cpu time
    fprintf(fp, "%6d   %E %E   %E   %E   %E   %E   %d:%02d:%02d   %d:%02d:%02d\n", statprops->runid, statprops->xcen,
            statprops->ycen, statprops->rmean, statprops->area, statprops->hmax, statprops->timereached, wallhours,
            wallminutes, wallseconds, cpuhours, cpuminutes, cpuseconds);
    fflush(fp);
    fclose(fp);

    return;
}

void InsanityCheck(ElementsHashTable* El_Table, int nump, int myid, TimeProps *timeprops_ptr)
{
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);

            if((Curr_El->refined_flag() == 0) && !((Curr_El->adapted_flag() == NOTRECADAPTED)
                    || (Curr_El->adapted_flag() == NEWFATHER) || (Curr_El->adapted_flag() == NEWSON)
                    || (Curr_El->adapted_flag() == BUFFER)))
            {
                cout<<"FUBAR 1 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";

                assert(0);
            }

            if((Curr_El->refined_flag() != 0) && ((Curr_El->adapted_flag() == NOTRECADAPTED)
                    || (Curr_El->adapted_flag() == NEWFATHER) || (Curr_El->adapted_flag() == NEWSON)
                    || (Curr_El->adapted_flag() == BUFFER)))
            {
                cout<<"FUBAR 2 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";
                assert(0);
            }

            if((Curr_El->refined_flag() == GHOST) && !((Curr_El->adapted_flag() <= -NOTRECADAPTED)
                    && (Curr_El->adapted_flag() >= -BUFFER)))
            {
                cout<<"FUBAR 3 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";
                assert(0);
            }

            if((Curr_El->refined_flag() != GHOST) && ((Curr_El->adapted_flag() <= -NOTRECADAPTED)
                    && (Curr_El->adapted_flag() >= -BUFFER)))
            {
                cout<<"FUBAR 4 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";
                assert(0);
            }

            if((Curr_El->refined_flag() > 0) && !((Curr_El->adapted_flag() == TOBEDELETED)
                    || (Curr_El->adapted_flag() == OLDFATHER) || (Curr_El->adapted_flag() == OLDSON)))
            {
                cout<<"FUBAR 5 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";
                assert(0);
            }

            if(!(Curr_El->refined_flag() > 0) && ((Curr_El->adapted_flag() == TOBEDELETED)
                    || (Curr_El->adapted_flag() == OLDFATHER) || (Curr_El->adapted_flag() == OLDSON)))
            {
                cout<<"FUBAR 6 in InsanityCheck()\nnump="<<nump<<" myid="<<myid<<" iter="<<timeprops_ptr->iter;
                cout<<" time="<<timeprops_ptr->timesec()<<"[sec]\nElement={"<<Curr_El->key();
                cout<<"} myprocess="<<Curr_El->myprocess()<<" refined="<<Curr_El->refined_flag()<<" adapted="<<Curr_El->adapted_flag()<<"\n";
                assert(0);
            }
        }
    }

    return;
}
