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
#include "../header/hadapt.h"

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
TitanProfiling titanProfiling;
int threads_number;


#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "../header/titan_simulation.h"
#include <math.h>
#include "../header/constant.h"
#include "../header/properties.h"

int NUM_STATE_VARS;
bool SHORTSPEED;

cxxTitanSimulation::cxxTitanSimulation():
        integrator(nullptr)
{
    elementType=ElementType::UnknownElementType;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanSimulation::~cxxTitanSimulation()
{
    FREE_VAR_IF_NOT_NULLPTR(integrator);
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
    NodeTable=NULL;
    ElemTable=NULL;

    NUM_STATE_VARS = 3;
    SHORTSPEED=false;
    get_outline()->elementType=ElementType::SinglePhase;
    elementType=ElementType::SinglePhase;
    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanSinglePhase::~cxxTitanSinglePhase()
{

}
void cxxTitanSinglePhase::set_short_speed(bool short_speed)
{
    SHORTSPEED=short_speed;
}

void cxxTitanSinglePhase::process_input()
{

    int i;
    int isrc;
    double doubleswap;

    MatProps *matprops_ptr=get_matprops();

    PileProps* pileprops_ptr=get_pileprops();
    FluxProps* fluxprops_ptr=get_fluxprops();

    StatProps* statprops_ptr=get_statprops();
    TimeProps* timeprops_ptr=get_timeprops();
    MapNames* mapnames_ptr=get_mapnames();
    OutLine* outline_ptr=get_outline();

    /*************************************************************************/

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
    matprops_ptr->set_scale(length_scale, height_scale, gravity_scale,pileprops_ptr,fluxprops_ptr);
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
        char *gis_matmap = (char *) malloc((mapnames_ptr->gis_map.size() + 5) * sizeof(char));
        strcpy(gis_matmap, mapnames_ptr->gis_map.c_str());
        strcat(gis_matmap + strlen(mapnames_ptr->gis_map.c_str()), "_Mat");

        if(Initialize_Raster_data(mapnames_ptr->gis_main.c_str(), mapnames_ptr->gis_sub.c_str(), mapnames_ptr->gis_mapset.c_str(), gis_matmap))
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
    timeprops_ptr->scale(TIME_SCALE);

    /*************************************************************************/

    int extramaps=0;
    if(use_gis_matmap)mapnames_ptr->extramaps=1;

    /*************************************************************************/
    // read in GIS information
    i = Initialize_GIS_data(mapnames_ptr->gis_main.c_str(), mapnames_ptr->gis_sub.c_str(), mapnames_ptr->gis_mapset.c_str(), mapnames_ptr->gis_map.c_str(), mapnames_ptr->gis_format);
    if(i != 0)
    {
        printf("Problem with GIS on processor %d\n", myid);
        exit(1);
    }

    matprops_ptr->calc_Vslump(pileprops_ptr,&fluxprops);
    /*************************************************************************/
    //test point information
    statprops_ptr->scale(matprops_ptr);

    /*************************************************************************/
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

    //-- MPI
    MPI_Status status;

    TIMING1_DEFINE(start);
    TIMING1_DEFINE(end);
    TIMING1_DEFINE(t_start);
    TIMING1_DEFINE(t_start2);
    TIMING1_DEFINE(t_end);

    TIMING1_START(start);

    /* create new MPI datastructures for class objects */
    MPI_New_Datatype();

    /* read original data from serial preprocessing
     code and then initialize element
     stiffness routines info */

    int xdmerr;



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

    MatProps *matprops_ptr=get_matprops();
    PileProps* pileprops_ptr=get_pileprops();
    StatProps* statprops_ptr=get_statprops();
    TimeProps* timeprops_ptr=get_timeprops();
    MapNames* mapnames_ptr=get_mapnames();
    OutLine* outline_ptr=get_outline();
    
    

    process_input();
    input_summary();



    int result=0;

    //check if restart is available
    //result=loadrun(myid, numprocs, &NodeTable,&ElemTable, matprops_ptr, &timeprops, &mapnames, &adapt, &order,
    //              &statprops, &discharge_planes, &outline);

    if(result==0)
        Read_grid(myid, numprocs, &NodeTable,&ElemTable, matprops_ptr, &outline);
    
    

    if(result==0)
    {
        setup_geoflow(ElemTable, NodeTable, myid, numprocs, matprops_ptr, &timeprops);

        move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

        AssertMeshErrorFree(ElemTable, NodeTable, numprocs, myid, -1.0);

        //initialize pile height and if appropriate perform initial adaptation
        init_piles();
    }
    else
    {
    	NodeTable->flushNodeTable();
        ElemTable->flushElemTable();
        //update temporary arrays of elements/nodes pointers
        ElemTable->updateLocalElements();
        ElemTable->updateNeighboursIndexes();
    }
    
    ElementsProperties ElemProp(ElemTable, NodeTable);

    HAdapt hadapt(ElemTable, NodeTable, &ElemProp,&timeprops,matprops_ptr,5);
    HAdaptUnrefine Unrefine(ElemTable, NodeTable,&timeprops,matprops_ptr);

    FREE_VAR_IF_NOT_NULLPTR(integrator);
    if(elementType == ElementType::TwoPhases)
    {
        if(order==1)
                    integrator=new Integrator(this);
    }
    if(elementType == ElementType::SinglePhase)
    {
        if(order==1)
            integrator=new Integrator_SinglePhase_CoulombMat_FirstOrder(this);
    }
    assert(integrator!=nullptr);

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
    calc_stats(elementType, ElemTable, NodeTable, myid, matprops_ptr, &timeprops, &statprops, &discharge_planes, 0.0);

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);

    move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);
    if(myid == 0)
        output_summary(&timeprops, &statprops, savefileflag);

    if(vizoutput & 1)
        tecplotter(elementType, ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

    if(vizoutput & 2)
        meshplotter(ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
    if(vizoutput & 4)
    {
        if(elementType == ElementType::TwoPhases)
        {
            xdmerr=write_xdmf_two_phases(ElemTable,NodeTable,&timeprops,matprops_ptr,&mapnames,XDMF_NEW);
        }
        if(elementType == ElementType::SinglePhase)
        {
            xdmerr=write_xdmf_single_phase(ElemTable,NodeTable,&timeprops,matprops_ptr,&mapnames,XDMF_NEW);
        }
    }

#endif

    if(vizoutput & 8)
    {
        if(myid == 0)
            grass_sites_header_output(&timeprops);
        grass_sites_proc_output(ElemTable, NodeTable, myid, matprops_ptr, &timeprops);
    }

    /*
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     Time Stepping Loop

     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     */

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

    ASSERT2(ElemTable->checkPointersToNeighbours("Prestep index check",false)==0);

    titanTimingsAlongSimulation.totalTime = MPI_Wtime();
    while (!(timeprops.ifend(0)) && !ifstop)
    {
        /*
         *  mesh adaption routines
         */
        TIMING1_START(t_start);
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
            update_topo(ElemTable, NodeTable, myid, numprocs, matprops_ptr, &timeprops, &mapnames);
        }

        if((adapt != 0) && (timeprops.iter % 5 == 4))
        {
            AssertMeshErrorFree(ElemTable, NodeTable, numprocs, myid, -2.0);

            TIMING1_START(t_start2);
            hadapt.adapt(h_count, TARGET);
            TIMING1_STOPADD(refinementTime, t_start2);
            
            
            TIMING1_START(t_start2);
            Unrefine.unrefine(UNREFINE_TARGET);

            //this move_data() here for debug... to make AssertMeshErrorFree() Work
            move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);
            TIMING1_STOPADD(unrefinementTime, t_start2);

            if((numprocs > 1) && (timeprops.iter % 10 == 9))
            {
                repartition2(ElemTable, NodeTable, &timeprops);

                //this move_data() here for debug... to make AssertMeshErrorFree() Work
                move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

                ElemTable->updateLocalElements();
                ElemTable->updateNeighboursIndexes();
            }


            //update temporary arrays of elements/nodes pointers
            NodeTable->flushNodeTable();
            ElemTable->flushElemTable();
            
            ASSERT2(ElemTable->checkPointersToNeighbours("After all adaptions",false)==0);
        }
        TIMING1_STOPADD(meshAdaptionTime, t_start);

        TIMING1_START(t_start);
        integrator->step();
        //step(matprops_ptr, &timeprops, pileprops_ptr, &fluxprops, &statprops,
        //     &order, &outline, &discharge_planes, adapt);
        TIMING1_STOPADD(stepTime, t_start);

        TIMING1_START(t_start);
        /*
         * save a restart file
         */
        if(timeprops.ifsave())
            saverun(&NodeTable, myid, numprocs, &ElemTable, matprops_ptr, &timeprops, &mapnames, adapt, order,
                    &statprops, &discharge_planes, &outline, &savefileflag);

        /*
         * output results to file
         */
        if(timeprops.ifoutput())
        {
            move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

            output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);

            if(myid == 0)
                output_summary(&timeprops, &statprops, savefileflag);

            if(vizoutput & 1)
                tecplotter(elementType, ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

            if(vizoutput & 2)
                meshplotter(ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
            if(vizoutput & 4)
            {
                if(elementType == ElementType::TwoPhases)
                {
                    xdmerr=write_xdmf_two_phases(ElemTable, NodeTable, &timeprops,
                                        matprops_ptr, &mapnames, XDMF_OLD);
                }
                if(elementType == ElementType::SinglePhase)
                {
                    xdmerr=write_xdmf_single_phase(ElemTable, NodeTable, &timeprops,
                                        matprops_ptr, &mapnames, XDMF_OLD);
                }
            }
#endif

            if(vizoutput & 8)
            {
                if(myid == 0)
                    grass_sites_header_output(&timeprops);
                grass_sites_proc_output(ElemTable, NodeTable, myid, matprops_ptr, &timeprops);
            }
        }

#ifdef PERFTEST
        int countedvalue=timeprops.iter%2+1;
        int e_buckets=ElemTable->get_no_of_buckets();
        HashEntry* entryp;
        for(i=0; i<e_buckets; i++)
        {
            entryp = *(ElemTable->getbucketptr() + i);
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
        TIMING1_STOPADD(resultsOutputTime, t_start);

        if(timeprops.iter % 200 == 0)
        {
            titanTimingsAlongSimulation.totalTime = MPI_Wtime() - titanTimingsAlongSimulation.totalTime;
            if(myid == 0)
                titanTimingsAlongSimulation.print("Timings over last 200 steps (seconds):");
            titanTimingsAlongSimulation.reset();
            titanTimingsAlongSimulation.totalTime = MPI_Wtime();
        }
    }

    move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

    /*
     * save a restart file
     */

    saverun(&NodeTable, myid, numprocs, &ElemTable, matprops_ptr, &timeprops, &mapnames, adapt, order,
            &statprops, &discharge_planes, &outline, &savefileflag);
    MPI_Barrier(MPI_COMM_WORLD);

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);
    MPI_Barrier(MPI_COMM_WORLD);

    if(myid == 0)
        output_summary(&timeprops, &statprops, savefileflag);

    if(vizoutput & 1)
        tecplotter(elementType, ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

    if(vizoutput & 2)
        meshplotter(ElemTable, NodeTable, matprops_ptr, &timeprops, &mapnames, statprops.vstar);

#if HAVE_LIBHDF5
    if(vizoutput & 4)
    {
        if(elementType == ElementType::TwoPhases)
        {
            xdmerr=write_xdmf_two_phases(ElemTable, NodeTable, &timeprops,
                        matprops_ptr, &mapnames, XDMF_CLOSE);
        }
        if(elementType == ElementType::SinglePhase)
        {
            xdmerr=write_xdmf_single_phase(ElemTable, NodeTable, &timeprops,
                        matprops_ptr, &mapnames, XDMF_CLOSE);
        }
    }
#endif

    if(vizoutput & 8)
    {
        if(myid == 0)
            grass_sites_header_output(&timeprops);
        grass_sites_proc_output(ElemTable, NodeTable, myid, matprops_ptr, &timeprops);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // write out ending warning, maybe flow hasn't finished moving
    sim_end_warning(elementType, ElemTable, matprops_ptr, &timeprops, statprops.vstar);
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

//TWO PHASES is:MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(*(outline.max_kinergy), *(outline2.max_kinergy), NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(*(outline.cum_kinergy), *(outline2.cum_kinergy), NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
    TIMING1_STOP(totalTime, start);
    if(myid == 0)
        titanTimings.print();
    if(myid == 0)
        titanProfiling.print();

    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}

void cxxTitanSinglePhase::input_summary()
{
    if(myid == 0)
    {
        get_mapnames()->print0();
        get_pileprops()->print0();
        get_fluxprops()->print0();
        get_discharge_planes()->print0();
        get_matprops()->print0();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}

cxxTitanTwoPhases::cxxTitanTwoPhases() :
        cxxTitanSinglePhase()
{
    NUM_STATE_VARS = 6;
    get_outline()->elementType=ElementType::TwoPhases;
    elementType=ElementType::TwoPhases;
    MPI_Barrier (MPI_COMM_WORLD);
}
cxxTitanTwoPhases::~cxxTitanTwoPhases()
{

}
