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
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#if HAVE_HDF5_H
#include "../header/hd5calls.h"
#endif

//#define DEBUG
//#define LOAD_BAL_DEBUG  //turns on a whole mess of mpi barriers so it makes run time more sensitive to load imbalances i.e. more sensitive to the load balance weights, it just makes it easier to adjust the constants.
//#define PERFTEST
#define TARGETPROC  -1
#define TARGETPROCA -1

int REFINE_LEVEL = 3;

#include "../header/titan2d_utils.h"

#include "../header/titan_simulation.h"

TitanTimings titanTimings;
TitanTimings titanTimingsAlongSimulation;


#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "../header/titan_simulation.h"
#include <math.h>
#include "../header/constant.h"
#include "../header/properties.h"

cxxTitanSimulation::cxxTitanSimulation()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gis_format = -1;

    topomain = "";
    toposub = "";
    topomapset = "";
    topomap = "";
    topovector = "";

    region_limits_set = false;

    min_location_x = 0.0;
    max_location_x = 0.0;
    min_location_y = 0.0;
    max_location_y = 0.0;

    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanSimulation::~cxxTitanSimulation()
{

}
void cxxTitanSimulation::run()
{

}
void cxxTitanSimulation::input_summary()
{
    if(myid != 0)
    {
        MPI_Barrier (MPI_COMM_WORLD);
        return;
    }

    printf("GIS:\n");
    printf("\tgis_format: %d\n", gis_format);

    printf("\ttopomain: %s\n", topomain.c_str());
    printf("\ttoposub: %s\n", toposub.c_str());
    printf("\ttopomapset: %s\n", topomapset.c_str());
    printf("\ttopomap: %s\n", topomap.c_str());
    printf("\ttopovector: %s\n", topovector.c_str());

    printf("\tregion_limits_set %d\n", (int) region_limits_set);

    /*int i;
    printf("Piles:\n");
    printf("\tNumber of piles: %d\n", (int) piles.size());
    for(i = 0; i < piles.size(); i++)
    {
        printf("\tPile %d:\n", i);
        piles[i].print0();
    }
    printf("Flux sources:\n");
    printf("\tNumber of flux sources: %d\n", (int) flux_sources.size());
    for(i = 0; i < flux_sources.size(); i++)
    {
        printf("\tFlux_source %d:\n", i);
        flux_sources[i].print0();
    }
    printf("Discharge planes:\n");
    printf("\tNumber of discharge planes: %d\n", (int) discharge_planes.size());
    for(i = 0; i < discharge_planes.size(); i++)
    {
        printf("\tDischarge plane %d:\n", i);
        discharge_planes[i].print0();
    }*/

    MPI_Barrier (MPI_COMM_WORLD);
}

cxxTitanSinglePhase::cxxTitanSinglePhase() :
        cxxTitanSimulation()
{
    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanSinglePhase::~cxxTitanSinglePhase()
{

}


void cxxTitanSinglePhase::process_input(StatProps* statprops_ptr,
               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr)
{

    int i;
    int isrc;
    double doubleswap;

    /*************************************************************************/


    PileProps* pileprops_ptr=get_pileprops();

    int no_of_sources = fluxprops.no_of_sources;
    if(fluxprops.no_of_sources+pileprops_ptr->numpiles==0)
    {
        printf("ERROR: No material source was defined");
        exit(1);
    }

#ifdef STATISTICS_IS_NOT_TRANSFERED
    /*************************************************************************/
    //over ride regular input with statistical sample run data
    FILE* fp;
    int lhsref, lhsid, ifstatbed = 0;
    double statbed, statint;  //friction angles
    if((fp = fopen("statin.bed", "r")) != NULL)
    {
        ifstatbed = 1;
        fscanf(fp, "%d%d%lf%lf", &(statprops_ptr->lhs.refnum), &(statprops_ptr->lhs.runid), &statbed, &statint);
        fclose(fp);
        statbed *= PI / 180.0;
        statint *= PI / 180.0;
    }
    else if((fp = fopen("statin.vol", "r")) != NULL)
    {
        double volumescale;
        fscanf(fp, "%d%d%lf%lf", &(statprops_ptr->lhs.refnum), &(statprops_ptr->lhs.runid), &volumescale);
        fclose(fp);
        volumescale = pow(volumescale, 1 / 3.0);

        for(i = 0; i < numpiles; i++)
        {
            pileprops_ptr->pileheight[i] *= volumescale;
            pileprops_ptr->majorrad[i] *= volumescale;
            pileprops_ptr->minorrad[i] *= volumescale;
        }
    }

    //PCQ
    int isample = -1, ifpcqvolbed = 0;
    double pcqbedfrict;
    fp = fopen("sample.number", "r");
    if(fp != NULL)
    {
        fscanf(fp, "%d", &isample);
        fclose(fp);
        statprops_ptr->lhs.runid = isample;
        double volume;
        int intswap;
        char samplefilename[256];
        sprintf(samplefilename, "dirfluxsample.%06d", isample);
        fp = fopen(samplefilename, "r");
        if(fp != NULL)
        {
            assert(fluxprops.no_of_sources);
            fluxprops.no_of_sources = 1;

            char stringswap[4096];
            fgets(stringswap, 4096, fp);  //Nsample
            fscanf(fp, "isample=%d\n", &intswap);
            assert(intswap == isample);
            fgets(stringswap, 4096, fp);  //weight
            fgets(stringswap, 4096, fp);  //Nrand
            fgets(stringswap, 4096, fp);  //Nvol
            fgets(stringswap, 4096, fp);  //Ndir
            fgets(stringswap, 4096, fp);  //randvar
            fscanf(fp, "volume=%lf\n", &volume);
            double vel = sqrt(fluxprops.xVel[0] * fluxprops.xVel[0] + fluxprops.yVel[0] * fluxprops.yVel[0]);
            double vel_angle;
            fscanf(fp, "direction=%lf\n", &vel_angle);
            vel_angle *= PI / 180.;
            fluxprops.xVel[0] = vel * cos(vel_angle);
            fluxprops.yVel[0] = vel * sin(vel_angle);
            fluxprops.start_time[0] = 0.0;
            fscanf(fp, "flux duration=%lf\n", &(fluxprops.end_time[0]));
            fscanf(fp, "init center flux=%lf\n", &(fluxprops.influx[0]));
            printf("Vol=%g [m^3]\n", volume);
        }
    }
#endif
    statprops_ptr->runid = statprops_ptr->lhs.runid;

    /*************************************************************************/
    MatProps *matprops_ptr=get_matprops();

    /*************************************************************************/
    matprops_ptr->set_scale(length_scale, height_scale, gravity_scale,pileprops_ptr,&fluxprops);
    matprops_ptr->process_input();

    double TIME_SCALE = matprops_ptr->get_TIME_SCALE();
    //non-dimensionalize the inputs
    double VELOCITY_SCALE = matprops_ptr->get_VELOCITY_SCALE();


    pileprops_ptr->scale(matprops_ptr->LENGTH_SCALE,matprops_ptr->HEIGHT_SCALE,matprops_ptr->GRAVITY_SCALE);
    fluxprops.scale(matprops_ptr->LENGTH_SCALE,matprops_ptr->HEIGHT_SCALE,matprops_ptr->GRAVITY_SCALE);

    double smallestpileradius=min(pileprops_ptr->get_smallest_pile_radius(),fluxprops.get_smallest_source_radius());

    matprops_ptr->smallest_axis = 2.0 * smallestpileradius;

    //read in material map
    if(matprops_ptr->material_count > 1)
    {
        char *gis_matmap = (char *) malloc((topomap.size() + 5) * sizeof(char));
        strcpy(gis_matmap, topomap.c_str());
        strcat(gis_matmap + strlen(topomap.c_str()), "_Mat");

        if(Initialize_Raster_data(topomain.c_str(), toposub.c_str(), topomapset.c_str(), gis_matmap))
        {
            printf("Problem with GIS Material on processor %d\n", myid);
            exit(1);
        }
        free(gis_matmap);

        int nummat;
        Get_raster_categories(&nummat);
        if(nummat != matprops_ptr->material_count)
        {
            printf("frict.data has %d materials but material map has %d, aborting\n", matprops_ptr->material_count,
                   nummat);
            exit(1);
        }
    }
    /*************************************************************************/
    //time related info
    timeprops_ptr->inittime(maxiter, maxtime, timeoutput, timesave, TIME_SCALE);

    /*************************************************************************/

    int extramaps=0;
    if(use_gis_matmap)extramaps=1;

    /*************************************************************************/
    // read in GIS information
    mapnames_ptr->assign(topomain.c_str(), toposub.c_str(), topomapset.c_str(), topomap.c_str(), gis_format, extramaps);
    i = Initialize_GIS_data(topomain.c_str(), toposub.c_str(), topomapset.c_str(), topomap.c_str(), gis_format);
    if(i != 0)
    {
        printf("Problem with GIS on processor %d\n", myid);
        exit(1);
    }

    matprops_ptr->calc_Vslump(pileprops_ptr,&fluxprops);
    /*************************************************************************/
    //test point information
    statprops_ptr->hxyminmax=edge_height;
    if(statprops_ptr->hxyminmax == -1.0)
    {
        statprops_ptr->hxyminmax = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT * 10.0;

    }
    else if(statprops_ptr->hxyminmax <= 0.0)
    {
        printf("bogus edge height=%g read in from simulation.data\n \
            (edge height = -1 is the flag for using the default height)\n \
            Exitting!!\n",
               statprops_ptr->hxyminmax);
        exit(1);
    }
    statprops_ptr->hxyminmax /= matprops_ptr->HEIGHT_SCALE;

    statprops_ptr->heightifreach=test_height;
    char chargarb1[64], chargarb2[64];
    if(statprops_ptr->heightifreach != -2)
    {

        //default test height is 10 time the maximum negligible height
        if(statprops_ptr->heightifreach == -1)
            statprops_ptr->heightifreach = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT * 10.0;

        statprops_ptr->heightifreach /= matprops_ptr->HEIGHT_SCALE;

        statprops_ptr->xyifreach[0]=test_location_x/matprops_ptr->LENGTH_SCALE;
        statprops_ptr->xyifreach[1]=test_location_y/matprops_ptr->LENGTH_SCALE;
    }
    else
    {
        statprops_ptr->heightifreach = statprops_ptr->xyifreach[0] = statprops_ptr->xyifreach[1] =
        HUGE_VAL;
    }

    //to get rid on uninitiallized memory error in saverun() (restart.C)
    statprops_ptr->forceint = statprops_ptr->forcebed = 0.0;

    /*************************************************************************/
    //the discharge plane section starts here
    int iplane, num_planes;
    num_planes=discharge_planes.num_planes;
    discharge_planes.scale(matprops_ptr->LENGTH_SCALE);

    //the discharge plane section ends here
    /*************************************************************************/



    /* physically we don't know how to specify a changing internal friction
     angle, we don't know how the contents of an avalanche changes as it
     picks up new material.  We don't know how to keep track of the material
     it picks up.  So we say that the internal friction angle is a constant
     */
#ifdef STATISTICS_IS_NOT_TRANSFERED
    //replace frict.data values of bed friction with ones from a
    //statistics sample run file read in above
    if(ifstatbed)
    {
        matprops_ptr->intfrict = statint;
        for(int imat = 1; imat <= matprops_ptr->material_count; imat++)
            matprops_ptr->bedfrict[imat] = statbed;
    }

    if(ifpcqvolbed)
    { //printf("PCQ bedfrict yada\n");
        matprops_ptr->bedfrict[1] = pcqbedfrict * PI / 180.0;

        double doubleswap = sqrt(tan(matprops_ptr->bedfrict[1]));
        for(int imat = 2; imat <= matprops_ptr->material_count; imat++)
        {
            matprops_ptr->bedfrict[imat] = matprops_ptr->bedfrict[1];
        }
    }
#endif
    /*************************************************************************/
    //to read in outline parameters here when it has been added
    return;
}
void cxxTitanSinglePhase::run()
{
    MPI_Barrier (MPI_COMM_WORLD);

    int i; //-- counters

    HashTable* BT_Node_Ptr;
    ElementsHashTable* BT_Elem_Ptr;

    //-- MPI
    MPI_Status status;

    double start, end;
    double t_start, t_end;
    start = MPI_Wtime();

    /* create new MPI datastructures for class objects */
    MPI_New_Datatype();

    /* read original data from serial preprocessing
     code and then initialize element
     stiffness routines info */

    int xdmerr;

    StatProps statprops;
    TimeProps timeprops;
    timeprops.starttime = time(NULL);

    MapNames mapnames;
    OutLine outline;

    double end_time = 10000.0;
    /*
     * viz_flag is used to determine which viz output to use
     * nonzero 1st bit of viz_flag means output tecplotxxxx.tec
     * nonzero 2nd bit of viz_flag means output mshplotxxxx.tec (debug purposes)
     * nonzero 3rd bit of viz_flag means output Paraview/XDMF format
     * nonzero 4th bit of viz_flag means output grass_sites files

     order_flag == 1 means use first order method
     order_flag == 2 means use second order method
     */

    //savefileflag will be flipped so first savefile will end in 0
    int savefileflag = 1;
    int Init_Node_Num, Init_Elem_Num;
    double v_star; // v/v_slump
    double nz_star; /* temporary... used for negligible velocity as stopping
     criteria paper... plan to include in v_star implicitly
     later */

    process_input(&statprops, &timeprops,
              &mapnames, &outline);
    input_summary();

    MatProps *matprops_ptr=get_matprops();
    PileProps* pileprops_ptr=get_pileprops();


    if(!loadrun(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, matprops_ptr, &timeprops, &mapnames, &adapt, &order,
                &statprops, &discharge_planes, &outline))
    {
        Read_grid(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, matprops_ptr, &outline);

        setup_geoflow(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, matprops_ptr, &timeprops);

        move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

        AssertMeshErrorFree(BT_Elem_Ptr, BT_Node_Ptr, numprocs, myid, -1.0);

        //initialize pile height and if appropriate perform initial adaptation
        init_piles(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, adapt, matprops_ptr, &timeprops, &mapnames, pileprops_ptr,
                   &fluxprops, &statprops);
    }
    else
    {
        //update temporary arrays of elements/nodes pointers
        BT_Elem_Ptr->updateElements();
        BT_Elem_Ptr->updateLocalElements();
        BT_Elem_Ptr->updatePointersToNeighbours();
    }

    /* for debug only, to check if exactly what's loaded will be saved again
     by doing a diff on the files.
     saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr,
     matprops_ptr, &timeprops, &mapnames, adaptflag, order_flag,
     &statprops, &discharge_planes, &outline, &savefileflag);
     */
    if(myid == 0)
    {
        for(int imat = 1; imat <= matprops_ptr->material_count; imat++)
            printf("bed friction angle for \"%s\" is %g\n", matprops_ptr->matnames[imat].c_str(),
                   matprops_ptr->bedfrict[imat] * 180.0 / PI);

        printf("internal friction angle is %g, epsilon is %g \n method order = %i\n", matprops_ptr->intfrict * 180.0 / PI,
               matprops_ptr->epsilon, order);
        printf("REFINE_LEVEL=%d\n", REFINE_LEVEL);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    calc_stats(BT_Elem_Ptr, BT_Node_Ptr, myid, matprops_ptr, &timeprops, &statprops, &discharge_planes, 0.0);

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);

    move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);
    if(myid == 0)
        output_summary(&timeprops, &statprops, savefileflag);

    if(vizoutput & 1)
        tecplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

    if(vizoutput & 2)
        meshplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
    if(vizoutput & 4)
    xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,matprops_ptr,&mapnames,XDMF_NEW);
#endif

    if(vizoutput & 8)
    {
        if(myid == 0)
            grass_sites_header_output(&timeprops);
        grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, matprops_ptr, &timeprops);
    }

    /*
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     Time Stepping Loop

     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     */
    printf("Time Stepping Loop\n");
    long element_counter = 0; // for performance count elements/timestep/proc
    int ifstop = 0;
    double max_momentum = 100;  //nondimensional

    /* ifend(0.5*statprops.vmean) is a hack, the original intent (when we were
     intending to use vstar as a stopping criteria) whas to have the
     calculation when vstar dropped back down below 1, instead we're
     using the ifend() function to stop the simulation when the volume
     averaged velocity falls back down below 2 meters... this hack is only
     for the colima hazard map runs, otherwise pass ifend() a constant
     valued */

    titanTimingsAlongSimulation.totalTime = MPI_Wtime();
    while (!(timeprops.ifend(0)) && !ifstop)
    {
        /*
         *  mesh adaption routines
         */
        t_start = MPI_Wtime();
        double TARGET = .05;
        double UNREFINE_TARGET = .01;
        int h_count = 0;
        if(timeprops.iter < 50)
            matprops_ptr->frict_tiny = 0.1;
        else
            matprops_ptr->frict_tiny = 0.000000001;

        //check for changes in topography and update if necessary
        if(timeprops.iter == 200)
        {
            update_topo(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, matprops_ptr, &timeprops, &mapnames);
        }

        if((adapt != 0) && (timeprops.iter % 5 == 4))
        {
            AssertMeshErrorFree(BT_Elem_Ptr, BT_Node_Ptr, numprocs, myid, -2.0);

            H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_count, TARGET, matprops_ptr, &fluxprops, &timeprops, 5);

            move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

            unrefine(BT_Elem_Ptr, BT_Node_Ptr, UNREFINE_TARGET, myid, numprocs, &timeprops, matprops_ptr);

            //this move_data() here for debug... to make AssertMeshErrorFree() Work
            move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

            if((numprocs > 1) && (timeprops.iter % 10 == 9))
            {
                repartition2(BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

                //this move_data() here for debug... to make AssertMeshErrorFree() Work
                move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);
            }
            move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

            //update temporary arrays of elements/nodes pointers
            BT_Elem_Ptr->updateElements();
            BT_Elem_Ptr->updateLocalElements();
            BT_Elem_Ptr->updatePointersToNeighbours();
        }
        titanTimings.meshAdaptionTime += MPI_Wtime() - t_start;
        titanTimingsAlongSimulation.meshAdaptionTime += MPI_Wtime() - t_start;

        t_start = MPI_Wtime();
        step(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, matprops_ptr, &timeprops, pileprops_ptr, &fluxprops, &statprops,
             &order, &outline, &discharge_planes, adapt);
        titanTimings.stepTime += MPI_Wtime() - t_start;
        titanTimingsAlongSimulation.stepTime += MPI_Wtime() - t_start;

        t_start = MPI_Wtime();
        /*
         * save a restart file
         */
        if(timeprops.ifsave())
            saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr, matprops_ptr, &timeprops, &mapnames, adapt, order,
                    &statprops, &discharge_planes, &outline, &savefileflag);

        /*
         * output results to file
         */
        if(timeprops.ifoutput())
        {
            move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

            output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);

            if(myid == 0)
                output_summary(&timeprops, &statprops, savefileflag);

            if(vizoutput & 1)
                tecplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

            if(vizoutput & 2)
                meshplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
            if(vizoutput & 4)
            xdmerr=write_xdmf(BT_Elem_Ptr, BT_Node_Ptr, &timeprops,
                    matprops_ptr, &mapnames, XDMF_OLD);
#endif

            if(vizoutput & 8)
            {
                if(myid == 0)
                    grass_sites_header_output(&timeprops);
                grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, matprops_ptr, &timeprops);
            }
        }

#ifdef PERFTEST
        int countedvalue=timeprops.iter%2+1;
        int e_buckets=BT_Elem_Ptr->get_no_of_buckets();
        HashEntry* entryp;
        for(i=0; i<e_buckets; i++)
        {
            entryp = *(BT_Elem_Ptr->getbucketptr() + i);
            while(entryp)
            {
                Element * EmTemp = (Element*)entryp->value;
                assert(EmTemp);
                assert(EmTemp->get_counted()!=countedvalue);

                if((EmTemp->get_adapted_flag()>=NOTRECADAPTED)&&
                        (EmTemp->get_adapted_flag()<=BUFFER)
                )
                {
                    //if this element doesn't belong on this processor don't involve
                    element_counter++;
                    EmTemp->put_counted(countedvalue);
                }
                entryp = entryp->next;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        titanTimings.resultsOutputTime += MPI_Wtime() - t_start;
        titanTimingsAlongSimulation.resultsOutputTime += MPI_Wtime() - t_start;

        if(timeprops.iter % 200 == 0)
        {
            titanTimingsAlongSimulation.totalTime = MPI_Wtime() - titanTimingsAlongSimulation.totalTime;
            if(myid == 0)
                titanTimingsAlongSimulation.print("Timings over last 200 steps (seconds):");
            titanTimingsAlongSimulation.reset();
            titanTimingsAlongSimulation.totalTime = MPI_Wtime();
        }
    }

    move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

    /*
     * save a restart file
     */

    saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr, matprops_ptr, &timeprops, &mapnames, adapt, order,
            &statprops, &discharge_planes, &outline, &savefileflag);
    MPI_Barrier(MPI_COMM_WORLD);

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);
    MPI_Barrier(MPI_COMM_WORLD);

    if(myid == 0)
        output_summary(&timeprops, &statprops, savefileflag);

    if(vizoutput & 1)
        tecplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

    if(vizoutput & 2)
        meshplotter(BT_Elem_Ptr, BT_Node_Ptr, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
    if(vizoutput & 4)
    xdmerr=write_xdmf(BT_Elem_Ptr, BT_Node_Ptr, &timeprops,
            matprops_ptr, &mapnames, XDMF_CLOSE);
#endif

    if(vizoutput & 8)
    {
        if(myid == 0)
            grass_sites_header_output(&timeprops);
        grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, matprops_ptr, &timeprops);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // write out ending warning, maybe flow hasn't finished moving
    sim_end_warning(BT_Elem_Ptr, matprops_ptr, &timeprops, statprops.vstar);
    MPI_Barrier(MPI_COMM_WORLD);

    //write out the final pile statistics (and run time)
    if(myid == 0)
        out_final_stats(&timeprops, &statprops);

    MPI_Barrier(MPI_COMM_WORLD);

    //write out stochastic simulation statistics
    if(myid == 0)
        output_stoch_stats(matprops_ptr, &statprops);
    MPI_Barrier(MPI_COMM_WORLD);

    //output maximum flow depth a.k.a. flow outline
    if(numprocs > 1)
    {
        OutLine outline2;
        double dxy[2];
        dxy[0] = outline.dx;
        dxy[1] = outline.dy;
        outline2.init2(dxy, outline.xminmax, outline.yminmax);
        int NxNyout = outline.Nx * outline.Ny;

//TWO_PHASES is:MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#ifndef TWO_PHASES
        MPI_Reduce(*(outline.max_kinergy), *(outline2.max_kinergy), NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(*(outline.cum_kinergy), *(outline2.cum_kinergy), NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        if(myid == 0)
            outline2.output(matprops_ptr, &statprops);
    }
    else
        outline.output(matprops_ptr, &statprops);

#ifdef PERFTEST
    long m = element_counter, ii;
    MPI_Allreduce ( &element_counter, &ii, 1,
            MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
    end=MPI_Wtime();
    char perffilename[256];
    sprintf(perffilename,"perform%04d.%04d",numprocs,myid);
    FILE *fpperf=fopen(perffilename,"w");
    fprintf(fpperf,"%d Finished -- used %ld elements of %ld total in %e seconds, %e\n",myid,m,ii,end-start, ii/(end-start));
    fclose(fpperf);
#endif

    titanTimings.totalTime = MPI_Wtime() - start;
    if(myid == 0)
        titanTimings.print();

    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}

void cxxTitanSinglePhase::input_summary()
{
    if(myid != 0)
    {
        MPI_Barrier (MPI_COMM_WORLD);
        return;
    }
    
    printf("GIS:\n");
    printf("\tgis_format: %d\n", gis_format);
    
    printf("\ttopomain: %s\n", topomain.c_str());
    printf("\ttoposub: %s\n", toposub.c_str());
    printf("\ttopomapset: %s\n", topomapset.c_str());
    printf("\ttopomap: %s\n", topomap.c_str());
    printf("\ttopovector: %s\n", topovector.c_str());
    
    printf("\tregion_limits_set %d\n", (int) region_limits_set);
    
    int i;
    get_pileprops()->print0();

    fluxprops.print0();

    discharge_planes.print0();

    get_matprops()->print0();

    

    MPI_Barrier (MPI_COMM_WORLD);
}

cxxTitanTwoPhases::cxxTitanTwoPhases() :
        cxxTitanSinglePhase()
{
    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanTwoPhases::~cxxTitanTwoPhases()
{

}
