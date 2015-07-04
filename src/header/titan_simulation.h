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

#ifndef TITAN2D_SIMULATION_H
#define TITAN2D_SIMULATION_H

#include <string>
#include <vector>
#include "../gisapi/GisApi.h"

#include "properties.h"

/**
 * MaterialMap
 */
class MaterialMap
{
public:
    MaterialMap();
    ~MaterialMap();

    MaterialMap& operator=(const MaterialMap& other);
private:
    int myid;
    int numprocs;
public:
    std::vector<std::string> name;
    std::vector<double> intfrict;
    std::vector<double> bedfrict;

    void print0();
    int get_material_count();
};

/**
 * cxxTitanSimulation
 */
class cxxTitanSimulation
{
public:
    cxxTitanSimulation();
    virtual ~cxxTitanSimulation();

    int myid;
    int numprocs;

    static const int GDAL = GDAL;
    static const int GIS_GRASS = GIS_GRASS;

    //GIS
    //!>GIS data format, 1 if GIS 0 if gdal
    int gis_format;
    //If using GIS GRASS data format:
    //!>the full path of the GIS database
    std::string topomain;
    //!>the location
    std::string toposub;
    //!>the mapset,
    std::string topomapset;
    //!>the raster map name, used for Gdal too
    std::string topomap;
    //!>GIS Vector
    std::string topovector;

    bool region_limits_set;

    double min_location_x;
    double max_location_x;
    double min_location_y;
    double max_location_y;

    //! length scaling factor
    double length_scale;

    //! height scaling factor
    //legacy, not used, all height
    //scaling now based on cube root of predicted volume, see below
    double height_scale;

    //! gravity scaling factor
    double gravity_scale;

    //! the "maximum" number of cells across the smallest pile/flux-source minor axis
    int number_of_cells_across_axis;

    //! the maximum # of iterations (a.k.a. time steps) before the simulation ends
    int maxiter;
    //! the maximum amount of time (in seconds) before the simulation ends
    double maxtime;
    //! the amount of time (in seconds) between subsequent outputs (with one exception... when the simulation ends one final output is performed and that one will be less than "timeoutput" seconds after the previous one
    double timeoutput;
    //! the amount of time (in seconds) between subsequent saves (with one exception... when the simulation ends one final save is performed and that one will be less than "timeoutput" seconds after the previous one
    double timesave;

    //! adapt
    int adapt;

    //!use a GIS material map
    bool use_gis_matmap;
    /**
     * vizoutput is used to determine which viz output to use
     * nonzero 1st bit of viz_flag means output tecplotxxxx.tec
     * nonzero 2nd bit of viz_flag means output mshplotxxxx.tec (debug purposes)
     * nonzero 3rd bit of viz_flag means output Paraview/XDMF format
     * nonzero 4th bit of viz_flag means output grass_sites files
     */
    int vizoutput;
    //! order == 1 means use first order method
    //! order == 2 means use second order method
    int order;

    //! xyminmax holds the minimum and maximum x and y coordinates where the pile height is greater than hxyminmax
    double edge_height;
    //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
    double test_height;
    //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
    double test_location_x;
    double test_location_y;





    //!>Process input and initiate dependencies, replacing Read_data
    virtual void process_input(MatProps* matprops_ptr, StatProps* statprops_ptr,
                               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr)
    {}

    virtual void hpfem(){}
    virtual void run();
    virtual void input_summary();

    virtual PileProps* get_pileprops(){return NULL;}
};

/**
 * cxxTitanSinglePhase
 */
class cxxTitanSinglePhase:public cxxTitanSimulation
{
public:
    cxxTitanSinglePhase();
    ~cxxTitanSinglePhase();

    //!>Process input and initiate dependencies, replacing Read_data
    virtual void process_input(MatProps* matprops_ptr, StatProps* statprops_ptr,
                               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr);
    virtual void hpfem();
    virtual void run();
    virtual void input_summary();

    virtual PileProps* get_pileprops(){return &pileprops;}



    //!>Piles
    PileProps pileprops;

    //!>Flux sources
    FluxProps fluxprops;

    //!>Discharge planes
    DischargePlanes discharge_planes;
    //std::vector<cxxTitanDischargePlane> discharge_planes;


    MaterialMap material_map;
};

class cxxTitanTwoPhases:public cxxTitanSinglePhase
{
public:
    cxxTitanTwoPhases();
    ~cxxTitanTwoPhases();

    //!>Piles
    PilePropsTwoPhases pileprops;

    virtual PileProps* get_pileprops(){return &pileprops;}
};
#endif
