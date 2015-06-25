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
#include <mpi.h>
#include "../header/titan_simulation.h"
#include <math.h>

cxxTitanPile::cxxTitanPile()
{
    height=0.0;
    xcenter = 1.0;
    ycenter = 1.0;
    majradius = 1.0;
    minradius = 1.0;
    orientation = 0.0;
    Vmagnitude = 0.0;
    Vdirection = 0.0;
}
cxxTitanPile::~cxxTitanPile()
{

}
cxxTitanPile& cxxTitanPile::operator=(const cxxTitanPile& other){
    if (this != &other) {
        height = other.height;
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
double cxxTitanPile::get_volume(){
    return M_PI*height*majradius*minradius/2.0;
}
void cxxTitanPile::print0()
{
    int i;
    //printf("Pile:\n");
    printf("\t\tMaximum Initial Thickness, P (m):%f\n", height);
    printf("\t\tCenter of Initial Volume, xc, yc (UTM E, UTM N): %f %f\n",
            xcenter, ycenter);
    printf("\t\tMajor and Minor Extent, majorR, minorR (m, m): %f %f\n",
            majradius, minradius);
    printf("\t\tOrientation (angle [degrees] from X axis to major axis): %f\n",
            orientation);
    printf("\t\tInitial speed [m/s]: %f\n", Vmagnitude);
    printf("\t\tInitial direction ([degrees] from X axis): %f\n", Vdirection);
    printf("\t\tPile volume [m^3]: %f\n", get_volume());
}

cxxTitanFluxSource::cxxTitanFluxSource()
{
    influx=0.0;
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
double cxxTitanFluxSource::get_effective_height(){
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
    double a = g/X/X;
    double t3 = X*X*influx/g;
    if(t3 > (end_time-start_time))
        t3=(end_time-start_time);
    return influx*t3 - 0.5*a*t3*t3;
}
cxxTitanFluxSource& cxxTitanFluxSource::operator=(const cxxTitanFluxSource& other){
    if (this != &other) {
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
    printf("\t\tActive Time [s], start, end: %f %f\n", start_time,end_time);

    printf("\t\tCenter of Initial Volume, xc, yc (UTM E, UTM N): %f %f\n",
            xcenter, ycenter);
    printf("\t\tMajor and Minor Extent, majorR, minorR (m, m): %f %f\n",
            majradius, minradius);
    printf("\t\tOrientation (angle [degrees] from X axis to major axis): %f\n",
            orientation);
    printf("\t\tInitial speed [m/s]: %f\n", Vmagnitude);
    printf("\t\tInitial direction ([degrees] from X axis): %f\n", Vdirection);
    printf("\t\tEffective Thickness, P (m):%f\n", get_effective_height());
}

cxxTitanDischargePlane::cxxTitanDischargePlane(){
    x_a = 0.0;
    y_a = 0.0;
    x_b = 0.0;
    y_b = 0.0;
}

cxxTitanDischargePlane::cxxTitanDischargePlane(const double m_x_a, const double m_y_a, const double m_x_b, const double m_y_b){
    x_a = m_x_a;
    y_a = m_y_a;
    x_b = m_x_b;
    y_b = m_y_b;
}
cxxTitanDischargePlane::~cxxTitanDischargePlane(){

}

cxxTitanDischargePlane& cxxTitanDischargePlane::operator=(
        const cxxTitanDischargePlane& other) {

    if (this != &other) {
        x_a = other.x_a;
        y_a = other.y_a;
        x_b = other.x_b;
        y_b = other.y_b;
    }
    return *this;
}

void cxxTitanDischargePlane::print0(){
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
MaterialMap& MaterialMap::operator=(const MaterialMap& other){
    if (this != &other) {
        int i;
        name.resize(other.name.size());
        intfrict.resize(other.intfrict.size());
        bedfrict.resize(other.bedfrict.size());
        for(i=0;i<other.name.size();i++){
            name[i]=other.name[i];
            intfrict[i]=other.intfrict[i];
            bedfrict[i]=other.bedfrict[i];
        }
    }
    return *this;
}
void MaterialMap::print0()
{
    if(myid!=0){
        //MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    int i;
    printf("Material map:\n");
    for(i=0;i<name.size();i++){
        printf("%d %s %f %f\n",i,name[i].c_str(),intfrict[i],bedfrict[i]);
    }

    //MPI_Barrier(MPI_COMM_WORLD);
}
int MaterialMap::get_material_count(){
    return name.size();
}
cxxTitanSimulation::cxxTitanSimulation()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gis_format=-1;

    topomain="";
    toposub="";
    topomapset="";
    topomap="";
    topovector="";

    region_limits_set=false;

    min_location_x=0.0;
    max_location_x=0.0;
    min_location_y=0.0;
    max_location_y=0.0;

    MPI_Barrier(MPI_COMM_WORLD);
}
cxxTitanSimulation::~cxxTitanSimulation()
{

}
int hpfem();
void cxxTitanSimulation::run()
{
	printf("cxxTitanSimulation::run %d\n",myid);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("cxxTitanSimulation::run::let's go %d\n",myid);
	hpfem();
	MPI_Barrier(MPI_COMM_WORLD);

}

void cxxTitanSimulation::input_summary()
{
    if(myid!=0){
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    printf("GIS:\n");
    printf("\tgis_format: %d\n",gis_format);

    printf("\ttopomain: %s\n",topomain.c_str());
    printf("\ttoposub: %s\n",toposub.c_str());
    printf("\ttopomapset: %s\n",topomapset.c_str());
    printf("\ttopomap: %s\n",topomap.c_str());
    printf("\ttopovector: %s\n",topovector.c_str());

    printf("\tregion_limits_set %d\n",(int)region_limits_set);

    int i;
    printf("Piles:\n");
    printf("\tNumber of piles: %d\n",(int)piles.size());
    for(i=0;i<piles.size();i++){
        printf("\tPile %d:\n",i);
        piles[i].print0();
    }
    printf("Flux sources:\n");
    printf("\tNumber of flux sources: %d\n",(int)flux_sources.size());
    for(i=0;i<flux_sources.size();i++){
        printf("\tFlux_source %d:\n",i);
        flux_sources[i].print0();
    }
    printf("Discharge planes:\n");
    printf("\tNumber of discharge planes: %d\n",(int)discharge_planes.size());
    for(i=0;i<discharge_planes.size();i++){
        printf("\tDischarge plane %d:\n",i);
        discharge_planes[i].print0();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
