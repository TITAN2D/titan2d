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

cxxTitanFluxSource::cxxTitanFluxSource()
{
    influx = 0.0;
    start_time = 0.0;
    end_time = 0.0;
    xcenter = 1.0;
    ycenter = 1.0;
    majradius = 1.0;
    minradius = 1.0;
    orientation = 0.0;
    Vmagnitude = 0.0;
    Vdirection = 0.0;
}
cxxTitanFluxSource::~cxxTitanFluxSource()
{
    
}
double cxxTitanFluxSource::get_effective_height()
{
    //approx: h=influx*t-0.5*a*t^2
    //if no s => t1=N*(2*h/g)^0.5  N is a empirical constant,
    //for cylindrical piles of aspect ratio (height/radius) of approx 1
    //2<=N<=3 (closer to 2) but there are 3 reasons we should increase N
    //(1) cylindrical pile does not collapse the whole way, shorter
    //distance means decreased acceleration means increased time, N
    //(2) paraboloid piles are closer to conical than cylinder so it
    //should collapse even less, so increase N
    //(3) "influx" is a constant source "velocity" not an initial
    //velocity which should increase h in "approx: h=..." equation, so
    //as a fudge factor increase N some more
    //calibrated on a single starting condition at tungaruhau says
    //N=3.21   N=X
    //anyway a=2*h/t1^2 = g/N^2
    //approx: v=influx-a*t2 at hmax v=0 => t2=influx/a = N^2*influx/g
    //t3=min(t2,end_time-start_time)
    //plug int first equation
    //approx hmax=influx*t3-0.5*a*t3^2
    //if t3==t2=> hmax= N^2/2*s^2/g
    //DEM: tungfla2
    //influx 12 m/s (vel 50 m/s at +35 degrees cc from +x direction
    //starts in filled crater which means this velocity points up hill
    //so pile is mostly stationary while flux source is active, 50 m/s
    //is just short of what is needed to top the crater)
    //end_time-start_time=20 gives actual hmax=75.6 m
    //g=9.8 m/s^2, N=3.21, t3=t2=12.62<20 s => computed hmax=75.7 m
    double X = 3.21;
    double g = 9.8;
    double a = g / X / X;
    double t3 = X * X * influx / g;
    if(t3 > (end_time - start_time))
        t3 = (end_time - start_time);
    return influx * t3 - 0.5 * a * t3 * t3;
}
cxxTitanFluxSource& cxxTitanFluxSource::operator=(const cxxTitanFluxSource& other)
{
    if(this != &other)
    {
        influx = other.influx;
        start_time = other.start_time;
        end_time = other.end_time;
        xcenter = other.xcenter;
        ycenter = other.ycenter;
        majradius = other.majradius;
        minradius = other.minradius;
        orientation = other.orientation;
        Vmagnitude = other.Vmagnitude;
        Vdirection = other.Vdirection;
    }
    return *this;
}
void cxxTitanFluxSource::print0()
{
    int i;
    //printf("Pile:\n");
    printf("\t\tExtrusion flux rate [m/s]:%f\n", influx);
    printf("\t\tActive Time [s], start, end: %f %f\n", start_time, end_time);
    
    printf("\t\tCenter of Initial Volume, xc, yc (UTM E, UTM N): %f %f\n", xcenter, ycenter);
    printf("\t\tMajor and Minor Extent, majorR, minorR (m, m): %f %f\n", majradius, minradius);
    printf("\t\tOrientation (angle [degrees] from X axis to major axis): %f\n", orientation);
    printf("\t\tInitial speed [m/s]: %f\n", Vmagnitude);
    printf("\t\tInitial direction ([degrees] from X axis): %f\n", Vdirection);
    printf("\t\tEffective Thickness, P (m):%f\n", get_effective_height());
}

cxxTitanDischargePlane::cxxTitanDischargePlane()
{
    x_a = 0.0;
    y_a = 0.0;
    x_b = 0.0;
    y_b = 0.0;
}

cxxTitanDischargePlane::cxxTitanDischargePlane(const double m_x_a, const double m_y_a, const double m_x_b,
                                               const double m_y_b)
{
    x_a = m_x_a;
    y_a = m_y_a;
    x_b = m_x_b;
    y_b = m_y_b;
}
cxxTitanDischargePlane::~cxxTitanDischargePlane()
{
    
}

cxxTitanDischargePlane& cxxTitanDischargePlane::operator=(const cxxTitanDischargePlane& other)
{
    
    if(this != &other)
    {
        x_a = other.x_a;
        y_a = other.y_a;
        x_b = other.x_b;
        y_b = other.y_b;
    }
    return *this;
}

void cxxTitanDischargePlane::print0()
{
    printf("\t\tPoint A (UTM E, UTM N): %f, %f\n", x_a, y_a);
    printf("\t\tPoint B (UTM E, UTM N): %f, %f\n", x_b, y_b);
}

MaterialMap::MaterialMap()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}
MaterialMap::~MaterialMap()
{
    
}
MaterialMap& MaterialMap::operator=(const MaterialMap& other)
{
    if(this != &other)
    {
        int i;
        name.resize(other.name.size());
        intfrict.resize(other.intfrict.size());
        bedfrict.resize(other.bedfrict.size());
        for(i = 0; i < other.name.size(); i++)
        {
            name[i] = other.name[i];
            intfrict[i] = other.intfrict[i];
            bedfrict[i] = other.bedfrict[i];
        }
    }
    return *this;
}
void MaterialMap::print0()
{
    if(myid != 0)
    {
        //MPI_Barrier(MPI_COMM_WORLD);
        return;
    }
    
    int i;
    printf("Material map:\n");
    for(i = 0; i < name.size(); i++)
    {
        printf("%d %s %f %f\n", i, name[i].c_str(), intfrict[i], bedfrict[i]);
    }
    
    //MPI_Barrier(MPI_COMM_WORLD);
}
int MaterialMap::get_material_count()
{
    return name.size();
}

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


void cxxTitanSinglePhase::process_input(MatProps* matprops_ptr, StatProps* statprops_ptr,
               TimeProps* timeprops_ptr, FluxProps* fluxprops, MapNames *mapnames_ptr, DISCHARGE* discharge_ptr, OutLine* outline_ptr)
{
    int isrc;
    /*************************************************************************/


#ifdef TWO_PHASES
    double maxphi = 0.;
    double minphi = HUGE_VAL;
    for(isrc = 0; isrc < numpiles; isrc++)
    {
        // search for min-max phi
        if(pileprops_ptr->vol_fract[isrc] > maxphi)
            maxphi = pileprops_ptr->vol_fract[isrc];
        if(pileprops_ptr->vol_fract[isrc] < minphi)
            minphi = pileprops_ptr->vol_fract[isrc];
    }
    // cut-off extremes
    assert(minphi <= maxphi);
    if(minphi < 0.2)
        matprops_ptr->flow_type = FLUID_FLOW;
    else if(maxphi > 0.9)
        matprops_ptr->flow_type = DRY_FLOW;
    else
        matprops_ptr->flow_type = TWOPHASE;
#endif

    int no_of_sources = flux_sources.size();
    if(no_of_sources>0)
    {
        double rotang = 0;
        double vel, vel_angle;
        fluxprops->allocsrcs(no_of_sources);

        for(isrc = 0; isrc < no_of_sources; isrc++)
        {
            fluxprops->influx[isrc]=flux_sources[isrc].influx;
            fluxprops->start_time[isrc]=flux_sources[isrc].start_time;
            fluxprops->end_time[isrc]=flux_sources[isrc].end_time;
            fluxprops->xCen[isrc]=flux_sources[isrc].xcenter;
            fluxprops->yCen[isrc]=flux_sources[isrc].ycenter;
            fluxprops->majorrad[isrc]=flux_sources[isrc].majradius;
            fluxprops->minorrad[isrc]=flux_sources[isrc].minradius;
            rotang=flux_sources[isrc].orientation;
            fluxprops->cosrot[isrc] = cos(rotang * PI / 180);
            fluxprops->sinrot[isrc] = sin(rotang * PI / 180);
            vel=flux_sources[isrc].Vmagnitude;
            vel_angle=flux_sources[isrc].Vdirection * PI / 180.;
            fluxprops->xVel[isrc] = vel * cos(vel_angle);
            fluxprops->yVel[isrc] = vel * sin(vel_angle);

        }

    }
    if(no_of_sources+pileprops_ptr->numpiles==0)
    {
        printf("ERROR: No material source was defined");
        exit(1);
    }
    int i;
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
            assert(fluxprops->no_of_sources);
            fluxprops->no_of_sources = 1;

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
            double vel = sqrt(fluxprops->xVel[0] * fluxprops->xVel[0] + fluxprops->yVel[0] * fluxprops->yVel[0]);
            double vel_angle;
            fscanf(fp, "direction=%lf\n", &vel_angle);
            vel_angle *= PI / 180.;
            fluxprops->xVel[0] = vel * cos(vel_angle);
            fluxprops->yVel[0] = vel * sin(vel_angle);
            fluxprops->start_time[0] = 0.0;
            fscanf(fp, "flux duration=%lf\n", &(fluxprops->end_time[0]));
            fscanf(fp, "init center flux=%lf\n", &(fluxprops->influx[0]));
            printf("Vol=%g [m^3]\n", volume);
        }
    }
#endif
    statprops_ptr->runid = statprops_ptr->lhs.runid;

    /*************************************************************************/
    //scaling info
    matprops_ptr->LENGTH_SCALE=length_scale;
    //all height scaling now based on cube root of predicted volume, see below
    matprops_ptr->HEIGHT_SCALE=height_scale;
    matprops_ptr->GRAVITY_SCALE=gravity_scale;

    double doubleswap;

    PileProps* pileprops_ptr=get_pileprops();

    //this is used in ../geoflow/stats.C ... might want to set
    //MAX_NEGLIGIBLE_HEIGHT to zero now that we have "good" thin
    //layer control, need to reevaluate this, we should also
    //reevaluate after we implement a "good" local stopping criteria
    matprops_ptr->MAX_NEGLIGIBLE_HEIGHT = matprops_ptr->HEIGHT_SCALE / 10000.0;

    if(height_scale==0.0)
    {
        double totalvolume = 0.0;

        if(pileprops_ptr->numpiles > 0)
            for(isrc = 0; isrc < pileprops_ptr->numpiles; isrc++)
                totalvolume += 0.5 * PI * pileprops_ptr->pileheight[isrc] * pileprops_ptr->majorrad[isrc]
                               * pileprops_ptr->minorrad[isrc];

        if(fluxprops->no_of_sources > 0)
            for(isrc = 0; isrc < fluxprops->no_of_sources; isrc++)
                totalvolume += 0.5 * PI * fluxprops->influx[isrc] * fluxprops->majorrad[isrc] * fluxprops->minorrad[isrc]
                               * 0.5 * (fluxprops->end_time[isrc] - //0.5 for linear decrease
                                       fluxprops->start_time[isrc]);

        doubleswap = pow(totalvolume, 1.0 / 3.0);

        if((matprops_ptr->GRAVITY_SCALE != 1.0) || (matprops_ptr->LENGTH_SCALE != 1.0))
            matprops_ptr->HEIGHT_SCALE = doubleswap;
        else
            matprops_ptr->HEIGHT_SCALE = 1.0;

        matprops_ptr->MAX_NEGLIGIBLE_HEIGHT = doubleswap / matprops_ptr->HEIGHT_SCALE / 10000.0;
    }

    matprops_ptr->epsilon = matprops_ptr->HEIGHT_SCALE / matprops_ptr->LENGTH_SCALE;

    double TIME_SCALE = sqrt(matprops_ptr->LENGTH_SCALE / matprops_ptr->GRAVITY_SCALE);

    //non-dimensionalize the inputs
    double VELOCITY_SCALE = sqrt(matprops_ptr->LENGTH_SCALE * matprops_ptr->GRAVITY_SCALE);

#ifdef TWO_PHASES
    double diameter = 0.005;
    double vterm = pow(diameter, 2.) * (matprops_ptr->den_solid - matprops_ptr->den_fluid) * matprops_ptr->GRAVITY_SCALE
            / (18. * matprops_ptr->viscosity);
    matprops_ptr->v_terminal = vterm;
#endif

    double smallestpileradius;

    pileprops_ptr->scale(matprops_ptr->LENGTH_SCALE,matprops_ptr->HEIGHT_SCALE,matprops_ptr->GRAVITY_SCALE);

    smallestpileradius=pileprops_ptr->get_smallest_pile_radius();

    for(isrc = 0; isrc < fluxprops->no_of_sources; isrc++)
    {
        fluxprops->influx[isrc] *= TIME_SCALE / (matprops_ptr->HEIGHT_SCALE);
        fluxprops->start_time[isrc] /= TIME_SCALE;
        fluxprops->end_time[isrc] /= TIME_SCALE;
        fluxprops->xCen[isrc] /= matprops_ptr->LENGTH_SCALE;
        fluxprops->yCen[isrc] /= matprops_ptr->LENGTH_SCALE;
        fluxprops->majorrad[isrc] /= matprops_ptr->LENGTH_SCALE;
        fluxprops->minorrad[isrc] /= matprops_ptr->LENGTH_SCALE;
        fluxprops->xVel[isrc] /= VELOCITY_SCALE;
        fluxprops->yVel[isrc] /= VELOCITY_SCALE;

        if(smallestpileradius > fluxprops->majorrad[isrc])
            smallestpileradius = fluxprops->majorrad[isrc];

        if(smallestpileradius > fluxprops->minorrad[isrc])
            smallestpileradius = fluxprops->minorrad[isrc];
    }

    matprops_ptr->smallest_axis = 2.0 * smallestpileradius;
    matprops_ptr->number_of_cells_across_axis=number_of_cells_across_axis;

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
    double **planes;
    num_planes=discharge_planes.size();

    if(num_planes > 0)
    {
        planes = CAllocD2(num_planes, 4);
        for(iplane = 0; iplane < num_planes; iplane++)
        {
            planes[iplane][0]=discharge_planes[iplane].x_a/matprops_ptr->LENGTH_SCALE;
            planes[iplane][1]=discharge_planes[iplane].y_a/matprops_ptr->LENGTH_SCALE;
            planes[iplane][2]=discharge_planes[iplane].x_b/matprops_ptr->LENGTH_SCALE;
            planes[iplane][3]=discharge_planes[iplane].y_b/matprops_ptr->LENGTH_SCALE;
        }
    }

    discharge_ptr->init(num_planes, planes);

    if(num_planes > 0)
        CDeAllocD2(planes);
    //the discharge plane section ends here
    /*************************************************************************/

    //read in material properties
    matprops_ptr->material_count=material_map.get_material_count();

    matprops_ptr->matnames = (char **) malloc((matprops_ptr->material_count + 1) * sizeof(char *));
    matprops_ptr->bedfrict = CAllocD1(matprops_ptr->material_count + 1);
    matprops_ptr->tanbedfrict = CAllocD1(matprops_ptr->material_count + 1);
    int imat;

    for(imat = 1; imat <= matprops_ptr->material_count; imat++)
    {
        matprops_ptr->matnames[imat] = allocstrcpy(material_map.name[imat-1].c_str());

        doubleswap=material_map.intfrict[imat-1];
        matprops_ptr->bedfrict[imat]=material_map.bedfrict[imat-1];
        matprops_ptr->bedfrict[imat] *= PI / 180.0;
        matprops_ptr->tanbedfrict[imat] = tan(matprops_ptr->intfrict);

        if(imat == 1)
        {
            matprops_ptr->intfrict = doubleswap * PI / 180.0;
            matprops_ptr->tanintfrict = tan(matprops_ptr->intfrict);
        }
    }

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
    }
    
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
