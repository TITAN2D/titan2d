/*******************************************************************
 * Copyright (C) 2015 University at Buffalo
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
#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "hpfem.h"

class cxxTitanSinglePhase;

/**
 * Base class for integrators
 *
 * these classes should implements 1 time step which consists of (by calling other functions) computing spatial
 * derivatives of state variables, computing k active/passive and wave speeds and therefore timestep size, does a
 * finite difference predictor step, followed by a finite volume corrector step, and lastly computing statistics
 * from the current timestep's data.
 *
 */
class Integrator:public EleNodeRef
{
public:
    Integrator(cxxTitanSinglePhase *_titanSimulation);
    virtual ~Integrator();

    //! Perform one integration step
    virtual void step();//=0;

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

protected:
    //!references to members of other classes
    int &order;
    ElementType &elementType;

    TimeProps *timeprops_ptr;
    MatProps *matprops_ptr;
    FluxProps* fluxprops_ptr;
    PileProps* pileprops_ptr;

    StatProps *statprops_ptr;
    OutLine *outline_ptr;
    DischargePlanes *discharge_ptr;

    int &adapt;

    //!its own members

    //!tiny flow
    double tiny;
    //!calculated current integration time step
    double dt;

    //!for comparison of magnitudes of forces in slumping piles
    double forceint;
    double forcebed;
    double eroded;
    double deposited;
    double realvolume;
};

/**
 * First order integrator for single phase and Coulomb material model
 */
class Integrator_SinglePhase_CoulombMat_FirstOrder:public Integrator
{
public:
    Integrator_SinglePhase_CoulombMat_FirstOrder(cxxTitanSinglePhase *_titanSimulation);

    //! Perform one integration step
    //virtual void step();
};

#endif
