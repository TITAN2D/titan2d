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
class NodeHashTable;

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

    //! length scaling factor
    double length_scale;

    //! height scaling factor
    //legacy, not used, all height
    //scaling now based on cube root of predicted volume, see below
    double height_scale;

    //! gravity scaling factor
    double gravity_scale;

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

    ElementType elementType;

    //!>Process input and initiate dependencies, replacing Read_data
    virtual void process_input(StatProps* statprops_ptr,
                               TimeProps* timeprops_ptr, MapNames *mapnames_ptr, OutLine* outline_ptr)
    {}
    virtual void run();
    virtual void input_summary();

    virtual PileProps* get_pileprops()=0;
    virtual MatProps* get_matprops()=0;
    virtual FluxProps* get_fluxprops()=0;
    virtual StatProps* get_statprops()=0;
    virtual TimeProps* get_timeprops()=0;
    virtual MapNames* get_mapnames()=0;
    virtual OutLine* get_outline()=0;
    virtual DischargePlanes* get_discharge_planes()=0;
    virtual NodeHashTable* get_HT_Node()=0;
    virtual ElementsHashTable* get_HT_Elem()=0;
};

/**
 * cxxTitanSinglePhase
 */
class cxxTitanSinglePhase:public cxxTitanSimulation
{
public:
    cxxTitanSinglePhase();
    ~cxxTitanSinglePhase();

    void set_short_speed(bool short_speed);

    //!>Process input and initiate dependencies, replacing Read_data
    virtual void process_input();

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

    StatProps statprops;
    TimeProps timeprops;
    MapNames mapnames;
    OutLine outline;

    NodeHashTable* NodeTable;
    ElementsHashTable* ElemTable;

    virtual PileProps* get_pileprops(){return &pileprops_single_phase;}
    virtual MatProps* get_matprops(){return &matprops_single_phase;}
    virtual FluxProps* get_fluxprops(){return &fluxprops;}
    virtual DischargePlanes* get_discharge_planes(){return &discharge_planes;}
    virtual StatProps* get_statprops(){return &statprops;}
    virtual TimeProps* get_timeprops(){return &timeprops;}
    virtual MapNames* get_mapnames(){return &mapnames;}
    virtual OutLine* get_outline(){return &outline;}

    virtual NodeHashTable* get_HT_Node(){return NodeTable;}
    virtual ElementsHashTable* get_HT_Elem(){return ElemTable;}

protected:
    /** this function intializes the piles, by commenting/uncommenting define statements you can switch from
     * parabaloid to elliptical cylinder shaped piles, or even a hard coded pileshapes written to match particular
     * experiments.  Adaptive remeshing and pile reinitialization helps detect small piles and refine around pile
     * edges to obtain a more accurate initial solution and speed up the first few timesteps before adaptive
     * refinement and unrefinement would otherwise occur.
     */
    void init_piles();

    //! this function implements 1 time step which consists of (by calling other functions) computing spatial derivatives of state variables, computing k active/passive and wave speeds and therefore timestep size, does a finite difference predictor step, followed by a finite volume corrector step, and lastly computing statistics from the current timestep's data.
    void step(MatProps* matprops_ptr,
              TimeProps* timeprops_ptr, PileProps *pileprops_ptr, FluxProps *fluxprops, StatProps* statprops_ptr,
              int* order_flag, OutLine* outline_ptr, DischargePlanes* discharge, int adaptflag);

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
