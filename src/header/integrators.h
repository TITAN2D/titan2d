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
#include "properties.h"

class cxxTitanSimulation;

/**
 * Base class for integrators
 *
 * these classes should implements 1 time step which consists of (by calling other functions) computing spatial
 * derivatives of state variables, computing k active/passive and wave speeds and therefore timestep size, does a
 * finite difference predictor step, followed by a finite volume corrector step, and lastly computing statistics
 * from the current timestep's data.
 *
 */
class Integrator:public EleNodeRef,public TiScalableObject
{
public:
    friend class cxxTitanSimulation;

    Integrator(cxxTitanSimulation *_titanSimulation);
    virtual ~Integrator();

    //! Perform one integration step
    virtual void step();//=0;

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);

public:
    //!internal friction for Coulomb model, some other function uses it for other models after cleaning
    //!up it should be only in Coulomb models
    double int_frict;
    double frict_tiny;
    //! order == 1 means use first order method
    //! order == 2 means use second order method
    int order;

protected:
    /**
     * Predictor step for second order of nothing for first order
     * should be implemented in derivative class
     */
    virtual void predictor()=0;
    /**
     * Corrector step for second order of whole step for first order
     * should be implemented in derivative class
     */
    virtual void corrector()=0;

protected:
    //!references to members of other classes
    const ElementType &elementType;

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
 * Base class for SinglePhase Integrators
 */
class Integrator_SinglePhase:public Integrator
{
public:
    Integrator_SinglePhase(cxxTitanSimulation *_titanSimulation);

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

public:
    double threshold;
    double erosion_rate;
    int do_erosion;
protected:
    //!properly named references
    tivector<double> &h;
    tivector<double> &hVx;
    tivector<double> &hVy;

    tivector<double> &dh_dx;
    tivector<double> &dh_dy;
    tivector<double> &dhVx_dx;
    tivector<double> &dhVx_dy;
    tivector<double> &dhVy_dx;
    tivector<double> &dhVy_dy;

};


/**
 * First order integrator for single phase and Coulomb material model
 */
class Integrator_SinglePhase_Coulomb_FirstOrder:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Coulomb_FirstOrder(cxxTitanSimulation *_titanSimulation);

    virtual void print0(int spaces=0);
protected:

    /**
     * Predictor step for second order of nothing for first order
     *
     * half timestep update (finite difference predictor finite volume corrector)
     *
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();


};

/**
 * First order integrator for single phase and Vollmey material model
 */
class Integrator_SinglePhase_Vollmey_FirstOrder:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Vollmey_FirstOrder(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

public:
    double mu;
    double xi;
};

/**
 * First order integrator for single phase and Pouliquen material model
 */
class Integrator_SinglePhase_Pouliquen_FirstOrder:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Pouliquen_FirstOrder(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

public:
    double phi1;
    double phi2;
    double partdiam;
    double I_O;
};

/**
 * First order integrator for single phase and Vollmey material model
 */
class Integrator_SinglePhase_Maeno_FirstOrder:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Maeno_FirstOrder(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

public:
    double phis;
    double phi2;
    double partdiam;
    double I_not;
};

/**
 * Base class for Two Phases Integrators
 */
class Integrator_TwoPhases:public Integrator
{
public:
    Integrator_TwoPhases(cxxTitanSimulation *_titanSimulation);

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
    //!properly named references
    tivector<double> &h;
    tivector<double> &h_liq;
    tivector<double> &hVx_sol;
    tivector<double> &hVy_sol;
    tivector<double> &hVx_liq;
    tivector<double> &hVy_liq;

    /*tivector<double> &dh_dx;
    tivector<double> &dh_dy;
    tivector<double> &dhVx_dx;
    tivector<double> &dhVx_dy;
    tivector<double> &dhVy_dx;
    tivector<double> &dhVy_dy;*/
};

/**
 * class for Two Phases Coulomb Integrator
 */
class Integrator_TwoPhases_Coulomb:public Integrator_TwoPhases
{
public:
    Integrator_TwoPhases_Coulomb(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();
};
#endif
