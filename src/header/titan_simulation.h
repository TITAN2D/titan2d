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

class ElementsHashTable;
class HashTable;

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
    virtual void process_input(StatProps* statprops_ptr,
                               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr)
    {}
    virtual void run();
    virtual void input_summary();

    virtual PileProps* get_pileprops()=0;
    virtual MatProps* get_matprops()=0;
    virtual FluxProps* get_fluxprops()=0;
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
    virtual void process_input(StatProps* statprops_ptr,
                               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr);

    virtual void run();
    virtual void input_summary();

    //!>Piles
    PileProps pileprops_single_phase;

    //!>Flux sources
    FluxProps fluxprops;

    //!>Discharge planes
    DischargePlanes discharge_planes;
    //std::vector<cxxTitanDischargePlane> discharge_planes;

    MatProps matprops_single_phase;

    virtual PileProps* get_pileprops(){return &pileprops_single_phase;}
    virtual MatProps* get_matprops(){return &matprops_single_phase;}
    virtual FluxProps* get_fluxprops(){return &fluxprops;}
protected:
    /** this function intializes the piles, by commenting/uncommenting define statements you can switch from
     * parabaloid to elliptical cylinder shaped piles, or even a hard coded pileshapes written to match particular
     * experiments.  Adaptive remeshing and pile reinitialization helps detect small piles and refine around pile
     * edges to obtain a more accurate initial solution and speed up the first few timesteps before adaptive
     * refinement and unrefinement would otherwise occur.
     */
    void init_piles(ElementsHashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, TimeProps* timeprops, MapNames* mapnames, StatProps* statprops);
public:
};

class cxxTitanTwoPhases:public cxxTitanSinglePhase
{
public:
    cxxTitanTwoPhases();
    ~cxxTitanTwoPhases();

    //!>Piles
    PilePropsTwoPhases pileprops_two_phases;

    MatPropsTwoPhases matprops_two_phases;

    virtual PileProps* get_pileprops(){return (PileProps*)&pileprops_two_phases;}
    virtual MatProps* get_matprops(){return (MatProps*)&matprops_two_phases;}
};
#endif
