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
 * $Id: outsum.C 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"
#include "../header/stats.hpp"

#include <sstream>
#define ADAM_HEIGHT_FRAC 0.02

void output_summary(TimeProps* timeprops, StatProps* statprops, int savefileflag)
{
    //FILE *fp=fopen("output_summary.readme","a");
    ostringstream filename;

    filename<<statprops->output_prefix<<"output_summary."<<setw(6)<< setfill('0') <<internal<<statprops->runid<<std::ends;
    FILE *fp = fopen(filename.str().c_str(), "at");
    
    if(timeprops->ifstart())
        fprintf(fp, "Summary of Output from Titan2D\n\n");
    
    int hours, minutes;
    double seconds;
    timeprops->chunktime(&hours, &minutes, &seconds);
    
    fprintf(fp,
            "%d:%02d:%g (hrs:min:sec)   iter=%8d   Vave=%g [m/s]   hmax=%g [m]   Xcen=%g [m]   Ycen=%g [m]   <(X-<X>)^2>=%g [m^2]   <(Y-<Y>)^2>=%g [m^2]   area=%g [m^2]   stat volume=%g [m^3]   real volume=%g [m^3]   eroded volume=%g [m^3]   deposited volume=%g [m^3]   outflow volume=%g [m^3]   forceint=%g [m/s^2]   forcebed=%g [m/s^2]   slope=%g   savefile=%d\n",
            hours, minutes, seconds, timeprops->iter, statprops->vmean, statprops->hmax, statprops->xcen,
            statprops->ycen, statprops->xvar, statprops->yvar, statprops->area, statprops->statvolume,
            statprops->realvolume, statprops->erodedvol, statprops->depositedvol, statprops->outflowvol,
            statprops->forceint, statprops->forcebed, statprops->slopemean, savefileflag);
    fclose(fp);
    return;
}

void output_stoch_stats(MatProps* matprops, StatProps* statprops)
{
    ostringstream filename;

    filename<<statprops->output_prefix<<"statout_lhs."<<setw(6)<< setfill('0') <<internal<<statprops->runid<<std::ends;
    FILE *fp = fopen(filename.str().c_str(), "w");
    fprintf(fp, "%2d %6d %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n", statprops->lhs.refnum,
            statprops->lhs.runid, statprops->timereached, statprops->vstar * matprops->Vslump, statprops->hmax,
            statprops->xcen, statprops->ycen, statprops->xyminmax[1] - statprops->xyminmax[0],
            statprops->xyminmax[3] - statprops->xyminmax[2]);
    fflush(fp);
    fclose(fp);
    return;
}

void output_discharge(MatProps* matprops, TimeProps* timeprops, DischargePlanes* discharge, int myid)
{
    double *send, *receive, doubleswap;
    int iplane, num_planes = discharge->num_planes;
    
    if(num_planes > 0)
    {
        send = CAllocD1(num_planes);
        receive = CAllocD1(num_planes);
        for(iplane = 0; iplane < num_planes; iplane++)
            send[iplane] = discharge->planes[iplane][9];
#ifdef USE_MPI
        MPI_Reduce(send, receive, num_planes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else //USE_MPI
    for(int i=0;i<num_planes;++i)receive[i]=send[i];
#endif //USE_MPI
        if(myid == 0)
        {
            doubleswap = (matprops->scale.length) * (matprops->scale.length) * (matprops->scale.height);
            
            for(iplane = 0; iplane < num_planes; iplane++)
                receive[iplane] *= doubleswap;
            
            FILE* fp = fopen("discharge.out", "a");
            if(timeprops->iter == 0)
            {
                fprintf(fp, "number of discharge planes = %g\n", num_planes);
                for(iplane = 0; iplane < num_planes; iplane++)
                    fprintf(fp, "plane %d endpoints are (%16.10g,%16.10g)"
                            " (%16.10g,%16.10g)\n",
                            iplane + 1, discharge->planes[iplane][0] * (matprops->scale.length),
                            discharge->planes[iplane][2] * (matprops->scale.length),
                            discharge->planes[iplane][1] * (matprops->scale.length),
                            discharge->planes[iplane][3] * (matprops->scale.length));
            }
            
            fprintf(fp, "\n%16.10g  %16.10g", timeprops->timesec(), receive[0]);
            for(iplane = 1; iplane < num_planes; iplane++)
                fprintf(fp, "  %16.10g", receive[iplane]);
            
            fclose(fp);
        }
        
        CDeAllocD1(send);
        CDeAllocD1(receive);
    }
    
    return;
}

void output_globalquants(TimeProps* timeprops, StatProps* statprops, int myid) {

	if (myid == 0) {
		FILE* fp = fopen("TemporalSpatial.info", "a");
		fprintf(fp, "%16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g\n",
				statprops->force_gx, statprops->force_gy, statprops->force_bx,
				statprops->force_by, statprops->force_bcx, statprops->force_bcy,
				statprops->force_rx, statprops->force_ry, statprops->power_g,
				statprops->power_b, statprops->power_bc, statprops->power_r,
				statprops->Vol_, statprops->Area_, statprops->Velmean_,
				timeprops->cur_time * timeprops->TIME_SCALE);
		fclose(fp);
	}
}

void output_localquants(TimeProps* timeprops, LocalQuants* localq, int myid) {
	int iloc, num_locs = localq->no_locations;
	char filename[50];

	if (num_locs > 0 && myid == 0) {
		for (iloc = 0; iloc < num_locs; iloc++) {
			if (timeprops->iter == 0) {
				sprintf(filename, "Location-%04d.dat", iloc);
				FILE* fp = fopen(filename, "a");
				fprintf(fp,
						"Time-history of QoIs are requested at location %d:\n",
						iloc);
				fprintf(fp, "\t(UTM E, UTM N): %g, %g\n",
						localq->length_scale * localq->X[iloc],
						localq->length_scale * localq->Y[iloc]);
				fprintf(fp, "\tFlow height     Flow Velocity     S_gx     S_gy     S_bedx     S_bedy     S_bedcurvx     S_bedcurvy     S_resistx      S_resisty      Power_g      Power_bed      Power_bedcurv      Power_resist      Slope_x      Slope_y      Time\n\n");

				if (localq->temps[iloc].size() == 0)
					fclose(fp);
				else {
					fprintf(fp,
							"%16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g\n",
							localq->Height[iloc], localq->Velocity[iloc],
							localq->Fgx[iloc], localq->Fgy[iloc],
							localq->Fbx[iloc], localq->Fby[iloc],
							localq->Fbcx[iloc], localq->Fbcy[iloc],
							localq->Fix[iloc], localq->Fiy[iloc],
							localq->Pg[iloc], localq->Pb[iloc],
							localq->Pbc[iloc], localq->Pi[iloc],
							localq->zetax[iloc], localq->zetay[iloc],
							timeprops->cur_time * timeprops->TIME_SCALE);
					fclose(fp);
				}
			} else if (localq->temps[iloc].size() > 0) {
				sprintf(filename, "Location-%04d.dat", iloc);
				FILE* fp = fopen(filename, "a");
				fprintf(fp,
						"%16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g\n",
						localq->Height[iloc], localq->Velocity[iloc],
						localq->Fgx[iloc], localq->Fgy[iloc],
						localq->Fbx[iloc], localq->Fby[iloc],
						localq->Fbcx[iloc], localq->Fbcy[iloc],
						localq->Fix[iloc], localq->Fiy[iloc],
						localq->Pg[iloc], localq->Pb[iloc],
						localq->Pbc[iloc], localq->Pi[iloc],
						localq->zetax[iloc], localq->zetay[iloc],
						timeprops->cur_time * timeprops->TIME_SCALE);
				fclose(fp);
			}
		}
	}
}

void output_localquantsTimeIntegrals(TimeProps* timeprops, LocalQuants* localq, int myid) {

	int iloc, num_locs = localq->no_locations;

	if (num_locs > 0) {
		double tempin[12], tempout[12];
		for (iloc = 0; iloc < num_locs; iloc++) {

		    tempin[0] = localq->T_Fgx[iloc];
		    tempin[1] = localq->T_Fgy[iloc];
		    tempin[2] = localq->T_Fbx[iloc];
		    tempin[3] = localq->T_Fby[iloc];
		    tempin[4] = localq->T_Fbcx[iloc];
		    tempin[5] = localq->T_Fbcy[iloc];
		    tempin[6] = localq->T_Fix[iloc];
		    tempin[7] = localq->T_Fiy[iloc];
		    tempin[8] = localq->T_Pg[iloc];
		    tempin[9] = localq->T_Pb[iloc];
		    tempin[10] = localq->T_Pbc[iloc];
		    tempin[11] = localq->T_Pi[iloc];

#ifdef USE_MPI
		    	MPI_Reduce(tempin, tempout, 12, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else //USE_MPI
		    	for(int i=0;i<12;++i)tempout[i]=tempin[i];
#endif //USE_MPI

		    if (myid == 0) {
		    	FILE* fp = fopen("TimeIntegraredLocalQuants.info", "a");
		    	fprintf(fp,"%16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g, %16.10g\n",
					localq->Tst_scale * tempout[0],localq->Tst_scale * tempout[1],
					localq->Tst_scale * tempout[2],localq->Tst_scale * tempout[3],
					localq->Tst_scale * tempout[4],localq->Tst_scale * tempout[5],
					localq->Tst_scale * tempout[6],localq->Tst_scale * tempout[7],
					localq->Tp_scale * tempout[8],localq->Tp_scale * tempout[9],
					localq->Tp_scale * tempout[10],localq->Tp_scale * tempout[11]);
				fclose(fp);
		    }
		}
	}
}
/*************************************/
/* adam's function: written by keith */
/*************************************/

void OUTPUT_ADAM_STATS(ElementsHashTable* El_Table, MatProps* matprops_ptr, TimeProps* timeprops_ptr, StatProps* statprops_ptr)
{
    int myid, numprocs;
    IF_MPI(MPI_Status status);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    double velocity2 = 0.0;
    double vmax = 0, hmax = 0;
    double xy_cen[2] =
    { 0.0, 0.0 }, vh_cen[2] =
    { 0.0, 0.0 };
    double xyh_vmax[3] =
    { 0.0, 0.0, 0.0 };
    double xyv_hmax[3] =
    { 0.0, 0.0, 0.0 };
    double masscenterdist2 = 0.0, masscentermindist2 = HUGE_VAL, xycen[2] =
    { 0.0, 0.0 };
    double vmax_min_height = matprops_ptr->scale.max_negligible_height * 512.0 * ADAM_HEIGHT_FRAC;
    int i;
    struct
    { //for use with MPI_MAXLOC
        double val;
        int rank;
    } send, receive;
    
    xy_cen[0] = statprops_ptr->xcen / (matprops_ptr->scale.length);
    xy_cen[1] = statprops_ptr->ycen / (matprops_ptr->scale.length);
    
    double VxVy[2];
    
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
                double height=EmTemp->state_vars(0);
                EmTemp->eval_velocity(0.0, 0.0, VxVy);
                velocity2 = VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1];

                //get v and h at center of mass
                masscenterdist2 = (EmTemp->coord(0) - xy_cen[0]) * (EmTemp->coord(0) - xy_cen[0])
                        + (EmTemp->coord(1) - xy_cen[1]) * (EmTemp->coord(1) - xy_cen[1]);
                if(masscenterdist2 < masscentermindist2)
                {
                    masscentermindist2 = masscenterdist2;
                    vh_cen[0] = velocity2;
                    vh_cen[1] = height;
                    xycen[0] = EmTemp->coord(0);
                    xycen[1] = EmTemp->coord(1);
                }

                //eliminate fast moving very thin pile from consideration
                if(height >= vmax_min_height)
                {

                    if(velocity2 > vmax)
                    {
                        /* velocity2 is not a mistake... only need to take the root of 
                         the maximum value */
                        vmax = velocity2;

                        xyh_vmax[0] = EmTemp->coord(0);
                        xyh_vmax[1] = EmTemp->coord(1);
                        xyh_vmax[2] = height;
                    }
                }

                if(height > hmax)
                {

                    hmax = height;

                    xyv_hmax[0] = EmTemp->coord(0);
                    xyv_hmax[1] = EmTemp->coord(1);
                    xyv_hmax[2] = velocity2;

                }
            }
        }
    }
    
    vh_cen[0] = sqrt(vh_cen[0]);
    vmax = sqrt(vmax);
    xyv_hmax[2] = sqrt(xyv_hmax[2]);
    
    /* get the max value accross all processors */
#ifdef USE_MPI
    if(numprocs > 1)
    {
        send.rank = myid;
        
        //at center of mass
        send.val = masscentermindist2;
        MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        if(receive.rank != 0)
        { /* don't send location if it's already on the 
         root processor */
            if(receive.rank == myid)
                MPI_Send(vh_cen, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            else if(myid == 0)
                MPI_Recv(vh_cen, 2, MPI_DOUBLE, receive.rank, 0, MPI_COMM_WORLD, &status);
        }
        
        //at location of vmax
        send.val = vmax;
        MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        vmax = receive.val;
        
        if(receive.rank != 0)
        { /* don't send location if it's already on the 
         root processor */
            if(receive.rank == myid)
                MPI_Send(xyh_vmax, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            else if(myid == 0)
                MPI_Recv(xyh_vmax, 3, MPI_DOUBLE, receive.rank, 0, MPI_COMM_WORLD, &status);
        }
        
        //at location of hmax
        send.val = hmax;
        MPI_Allreduce(&send, &receive, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        hmax = receive.val;
        
        if(receive.rank != 0)
        { /* don't send location if it's already on the 
         root processor */
            if(receive.rank == myid)
                MPI_Send(xyv_hmax, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            else if(myid == 0)
                MPI_Recv(xyv_hmax, 3, MPI_DOUBLE, receive.rank, 0, MPI_COMM_WORLD, &status);
        }
    }
#endif //USE_MPI
    
    if(myid == 0)
    {
        FILE *fp = fopen("flow_dynamics.stats", "a");
        
        fprintf(fp,
                "%16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g,   %16.10g\n",
                timeprops_ptr->timesec(), //time in seconds
                
                statprops_ptr->vmean, //average velocity
                
                //x,y,v,h at center of mass
                statprops_ptr->xcen,
                statprops_ptr->ycen,
                vh_cen[0] * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)),
                vh_cen[1] * (matprops_ptr->scale.height),
                
                //x,y,v,h at location of vmax
                xyh_vmax[0] * matprops_ptr->scale.length,
                xyh_vmax[1] * matprops_ptr->scale.length,
                vmax * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)),
                xyh_vmax[2] * (matprops_ptr->scale.height),
                
                //x,y,v,h at location of hmax
                xyv_hmax[0] * matprops_ptr->scale.length,
                xyv_hmax[1] * matprops_ptr->scale.length,
                xyv_hmax[2] * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity)),
                hmax * (matprops_ptr->scale.height));
        
        fclose(fp);
    }
    return;
}
