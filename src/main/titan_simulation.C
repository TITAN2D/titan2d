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
    matprops_ptr->set_scale(length_scale, height_scale, gravity_scale);
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

#ifdef TWO_PHASES
    /*************************************************************************/
    /* the non-dimensional velocity stopping criteria is an idea that
     didn't work for anything other than a slumping pile on a horizontal
     surface, it's only still here because I didn't want to bother
     with removing it.  --Keith Dalbey 2005.10.28

     kappa is a to be determined constant, calculation stops when
     v*=v_ave/v_slump<kappa (or perhaps v* < kappa/tan(intfrict)) */
    double kappa = 1.0;   //should eventually move to a header file
    double gravity = 9.8; //[m/s^2]
    matprops_ptr->Vslump = 1.0; //kappa*sqrt(gravity*max_init_height);
#endif
#ifndef TWO_PHASES
    /*************************************************************************
     * Vslump doesn't mean it is related to slumping of pile
     * It is simply the maximum hypothetical velocity from
     * free fall of the pile.
     ************************************************************************/
    double gravity = 9.8; //[m/s^2]
    double zmin, res;
    // get DEM resolution from GIS
    int ierr = Get_max_resolution(&res);
    if(ierr == 0)
    {
        // Get minimum finite elevation from GIS
        ierr = Get_elev_min(res, &zmin);
        if(ierr != 0)
            zmin = 0;
    }
    if(isnan (zmin) || (zmin < 0))
        zmin = 0;

    // search highest point amongst piles
    double zcen = 0;
    int j = 0;
    for(i = 0; i < pileprops_ptr->numpiles; i++)
    {
        double xcen = matprops_ptr->LENGTH_SCALE * pileprops_ptr->xCen[i];
        double ycen = matprops_ptr->LENGTH_SCALE * pileprops_ptr->yCen[i];
        double ztemp = 0;
        ierr = Get_elevation(res, xcen, ycen, &ztemp);
        if(ierr != 0)
            ztemp = 0;

        if((ztemp + pileprops_ptr->pileheight[i]) > (zcen + pileprops_ptr->pileheight[i]))
        {
            zcen = ztemp;
            j = i;
        }
    }

    // calculate Vslump
    double hscale = (zcen - zmin) + pileprops_ptr->pileheight[j];
    matprops_ptr->Vslump = sqrt(gravity * hscale);
#endif
    /*************************************************************************/
    //test point information
    statprops_ptr->hxyminmax=edge_height;
    if(statprops_ptr->hxyminmax == -1)
        statprops_ptr->hxyminmax = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT * 10.0;
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
int hpfem();
void cxxTitanSinglePhase::run()
{
    printf("cxxTitanSimulation::run %d\n", myid);
    MPI_Barrier (MPI_COMM_WORLD);
    printf("cxxTitanSimulation::run::let's go %d\n", myid);
#ifndef TWO_PHASES
    hpfem();
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    
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

    matprops.print0();

    

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
