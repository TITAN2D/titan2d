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
#include "stats.hpp"
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

    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

    //! Create integrater from hdf file content, will instantiate proper integrator class
    static Integrator* createIntegrator(const H5::CommonFG *parent, cxxTitanSimulation *_titanSimulation, const  string group_name="Integrator");

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

    //! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model)
    //! calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum
    //! allowable timestep for this iteration.
    //! should be implemented in derivative class
    virtual double get_coef_and_eigen(int ghost_flag)=0;

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords()=0;


protected:
    //!references to members of other classes
    const ElementType &elementType;

    ElementsProperties ElemProp;

    TimeProps *timeprops_ptr;
    MatProps *matprops_ptr;
    FluxProps* fluxprops_ptr;
    PileProps* pileprops_ptr;

    StatProps *statprops_ptr;
    OutLine *outline_ptr;
    DischargePlanes *discharge_ptr;
    LocalQuants *localquants_ptr;

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

    double force_gx;
    double force_gy;
    double force_bx;
    double force_by;
    double force_bcx;
    double force_bcy;
    double force_rx;
    double force_ry;
    double power_g;
    double power_b;
    double power_bc;
    double power_r;

    double Tforce_gx;
    double Tforce_bx;
    double Tforce_bcx;
    double Tforce_rx;
    double Tforce_gy;
    double Tforce_by;
    double Tforce_bcy;
    double Tforce_ry;
    double Tpower_g;
    double Tpower_b;
    double Tpower_bc;
    double Tpower_r;
};

/**
 * Base class for SinglePhase Integrators
 */
class Integrator_SinglePhase:public Integrator
{
public:
    Integrator_SinglePhase(cxxTitanSimulation *_titanSimulation);

    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");
protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

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
class Integrator_SinglePhase_Coulomb:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Coulomb(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);

    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

    int stopping_criteria;
    double thr;
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

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

    //! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model)
    //! calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum
    //! allowable timestep for this iteration.
    virtual double get_coef_and_eigen(int ghost_flag);


protected:
    //! calculation of k active passive
    virtual void gmfggetcoef_C(const double h,const double hVx,const double hVy,
            const double dh_dx,const double dhVx_dx,
            const double dh_dy,const double dhVy_dy,
            const double bedfrictang, const double intfrictang,
            double &Kactx, double &Kacty, const double tiny,
        const double epsilon)
    {
        //vel is used to determine if velocity gradients are converging or diverging
        double vel;

        //COEFFICIENTS
        double hSQ=h*h;
        double cosphiSQ = cos(intfrictang);
        double tandelSQ = tan(bedfrictang);
        cosphiSQ*=cosphiSQ;
        tandelSQ*=tandelSQ;



        if(h > tiny)
        {
             vel=dhVx_dx/h - hVx*dh_dx/hSQ+
                 dhVy_dy/h - hVy*dh_dy/hSQ;
             Kactx=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
                 sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;
             Kacty=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
                 sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;

             //if there is no yielding...
             if(fabs(hVx/h) < tiny && fabs(hVy/h) < tiny)
             {
                Kactx = 1.0;
                Kacty = 1.0;
             }
        }
        else
        {
        vel = 0.0;
        Kactx = 1.0;
        Kacty = 1.0;
        }
        Kactx = epsilon * Kactx;
        Kacty = epsilon * Kacty;
    }
    //! calculation of wave speeds (eigen vectors of the flux jacoboians)
    virtual void eigen_C( const double h, double &eigenvxmax, double &eigenvymax, double &evalue,
            const double tiny, double &kactx, const double gravity_z, const double *VxVy)
    {
        if (h > tiny)
        {
            //     iverson and denlinger
            if (kactx < 0.0)
            {
                //negative kactxy
                kactx = -kactx;
            }
            eigenvxmax = fabs(VxVy[0]) + sqrt(kactx * gravity_z * h);
            eigenvymax = fabs(VxVy[1]) + sqrt(kactx * gravity_z * h);

        }
        else
        {
            eigenvxmax = tiny;
            eigenvymax = tiny;
        }

        evalue = c_dmax1(eigenvxmax, eigenvymax);
    }
};

/**
 * First order integrator for single phase and Voellmy_Salm material model
 */
class Integrator_SinglePhase_Voellmy_Salm:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Voellmy_Salm(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

    //! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model)
    //! calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum
    //! allowable timestep for this iteration.
    virtual double get_coef_and_eigen(int ghost_flag);

public:
    double mu;
    double xi;
    double thr;

protected:
    //! calculation of k active passive
    virtual void gmfggetcoef_VS(double &Kactx, double &Kacty,const double epsilon)
    {
        Kactx = 1.0;
        Kacty = 1.0;

        Kactx = epsilon * Kactx;
        Kacty = epsilon * Kacty;
    }
    //! calculation of wave speeds (eigen vectors of the flux jacoboians)
    void eigen_VS( const double h, double &eigenvxmax, double &eigenvymax, double &evalue,
            const double tiny, double &kactx, const double gravity_z, const double *VxVy)
    {
        if (h > tiny)
        {
            //     iverson and denlinger
            if (kactx < 0.0)
            {
                //negative kactxy
                kactx = -kactx;
            }
            eigenvxmax = fabs(VxVy[0]) + sqrt(kactx * gravity_z * h);
            eigenvymax = fabs(VxVy[1]) + sqrt(kactx * gravity_z * h);

        }
        else
        {
            eigenvxmax = tiny;
            eigenvymax = tiny;
        }

        evalue = c_dmax1(eigenvxmax, eigenvymax);
    }

};

/**
 * First order integrator for single phase and Pouliquen_Forterre basal friction model
 */
class Integrator_SinglePhase_Pouliquen_Forterre:public Integrator_SinglePhase
{
public:
    Integrator_SinglePhase_Pouliquen_Forterre(cxxTitanSimulation *_titanSimulation);

    virtual bool scale();
    virtual bool unscale();
    virtual void print0(int spaces=0);
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

    //! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model)
    //! calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum
    //! allowable timestep for this iteration.
    virtual double get_coef_and_eigen(int ghost_flag);

public:
    double phi1;
    double phi2;
    double phi3;
    double Beta;
    double L_material;
    double thr;

protected:
    //! calculation of k active passive
    virtual void gmfggetcoef_PF(double &Kactx, double &Kacty,const double epsilon)
    {
        Kactx = 1.0;
        Kacty = 1.0;

        Kactx = epsilon * Kactx;
        Kacty = epsilon * Kacty;
    }
    //! calculation of wave speeds (eigen vectors of the flux jacoboians)
    void eigen_PF( const double h, double &eigenvxmax, double &eigenvymax, double &evalue,
            const double tiny, double &kactx, const double gravity_z, const double *VxVy)
    {
        if (h > tiny)
        {
            eigenvxmax = fabs(VxVy[0]) + sqrt(kactx * gravity_z * h);
            eigenvymax = fabs(VxVy[1]) + sqrt(kactx * gravity_z * h);

        }
        else
        {
            eigenvxmax = tiny;
            eigenvymax = tiny;
        }

        evalue = c_dmax1(eigenvxmax, eigenvymax);
    }
};

/**
 * Base class for Two Phases Integrators
 */
class Integrator_TwoPhases:public Integrator
{
public:
    Integrator_TwoPhases(cxxTitanSimulation *_titanSimulation);

    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

    //! this function computes k active/passive (which is necessary because of the use of the Coulomb friction model)
    //! calculates the wave speeds (eigen values of the flux jacobians) and based on them determines the maximum
    //! allowable timestep for this iteration.
    virtual double get_coef_and_eigen(int ghost_flag);

protected:
    //!properly named references
    tivector<double> &h;
    tivector<double> &h_liq;
    tivector<double> &hVx_sol;
    tivector<double> &hVy_sol;
    tivector<double> &hVx_liq;
    tivector<double> &hVy_liq;

    tivector<double> &dh_dx;
    tivector<double> &dh_dy;
    tivector<double> &dh_liq_dx;
    tivector<double> &dh_liq_dy;
    tivector<double> &dhVx_sol_dx;
    tivector<double> &dhVx_sol_dy;
    tivector<double> &dhVy_sol_dx;
    tivector<double> &dhVy_sol_dy;
    tivector<double> &dhVx_liq_dx;
    tivector<double> &dhVx_liq_dy;
    tivector<double> &dhVy_liq_dx;
    tivector<double> &dhVy_liq_dy;
protected:
    void gmfggetcoef2ph(const double h_liq,const double hVx_sol,const double hVy_sol,
            const double dh_dx_liq,const double dhVx_dx_sol,
            const double dh_dy_liq,const double dhVy_dy_sol,
            const double bedfrictang, const double intfrictang,
            double &Kactx, double &Kacty, const double tiny,
        const double epsilon)
    {
        //vel is used to determine if velocity gradients are converging or diverging
        double vel;

        //COEFFICIENTS
        double h_liqSQ=h_liq*h_liq;
        double cosphiSQ = cos(intfrictang);
        double tandelSQ = tan(bedfrictang);
        cosphiSQ*=cosphiSQ;
        tandelSQ*=tandelSQ;

        if(h_liq > tiny)
        {
             vel=dhVx_dx_sol/h_liq - hVx_sol*dh_dx_liq/h_liqSQ+
                 dhVy_dy_sol/h_liq - hVy_sol*dh_dy_liq/h_liqSQ;
             Kactx=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
                 sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;
             Kacty=(2.0/cosphiSQ)*(1.0-sgn_tiny(vel,tiny)*
                 sqrt(fabs(1.0-(1.0+tandelSQ)*cosphiSQ) )) -1.0;

             //if there is no yielding...
             if(fabs(hVx_sol/h_liq) < tiny && fabs(hVy_sol/h_liq) < tiny)
             {
                Kactx = 1.0;
                Kacty = 1.0;
             }
        }
        else
        {
        vel = 0.0;
        Kactx = 1.0;
        Kacty = 1.0;
        }
        Kactx = epsilon * Kactx;
        Kacty = epsilon * Kacty;
    }

    //@TODO is it really h2
    void eigen2ph( const double h_sol, const double h_liq, double &eigenvxmax, double &eigenvymax, double &evalue,
            const double tiny, double &kactx, const double gravity_z,
            const double *v_solid, const double *v_fluid, const int flowtype)

    {
        double sound_speed;
        if (h_sol > tiny) {
            //iverson and denlinger
            if (kactx < 0.0) {
                kactx = -kactx;
            }

            if (flowtype == 1)
                sound_speed = sqrt(h_sol * kactx * gravity_z);
            else if (flowtype == 2)
                sound_speed = sqrt(h_sol * gravity_z);
            else
                sound_speed = sqrt(h_liq * gravity_z * kactx+(h_sol - h_liq) * gravity_z);

            //x-direction
            eigenvxmax = c_dmax1(fabs(v_solid[0] + sound_speed), fabs(v_fluid[0] + sound_speed));

            //y-direction
            eigenvymax = c_dmax1(fabs(v_solid[1] + sound_speed), fabs(v_fluid[1] + sound_speed));
        }
        else {
            eigenvxmax = tiny;
            eigenvymax = tiny;
        }
        evalue = c_dmax1(eigenvxmax, eigenvymax);
    }
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
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="Integrator") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="Integrator");

protected:
    /**
     * Predictor step for second order of nothing for first order
     */
    virtual void predictor();
    /**
     * Corrector step for second order of whole step for first order
     */
    virtual void corrector();

    // This method is used when we need to get the records of flow characteristics such as forces.
    virtual void flowrecords();

    void calc_drag_force(const ti_ndx_t ndx, const double *vsolid, const double *vfluid,
            const double den_solid, const double den_fluid, const double vterminal,
            double *drag)
    {
        double temp, volf, denfrac;
        double delv[2];
        const double exponant = 3.0;

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!     model in pitman-le paper
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        volf = h_liq[ndx] / h[ndx];
        temp = h_liq[ndx] * pow(1. - volf, 1.0 - exponant) / vterminal;
        denfrac = den_fluid / den_solid;

        for (int i = 0; i < 4; ++i)
            drag[i] = 0.0;

        //fluid vel - solid vel
        if (h_liq[ndx] > tiny)
        {
            for (int i = 0; i < 2; ++i)
                delv[i] = vfluid[i] - vsolid[i];

            //compute individual drag-forces
            drag[0] = (1.0 - denfrac) * temp * delv[0];
            drag[1] = (1.0 - denfrac) * temp * delv[1];
            drag[2] = (1.0 - denfrac) * temp * delv[0] / denfrac;
            drag[3] = (1.0 - denfrac) * temp * delv[1] / denfrac;
        }
    }
};
#endif
