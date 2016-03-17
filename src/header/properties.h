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
 * $Id: properties.h 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifndef MAX_DEPTH_MAP
#define MAX_DEPTH_MAP
#endif

#ifndef PROPS
#define PROPS

#include "ticore.hpp"

#include<stdlib.h>
#include<stdio.h> //useful_lib.h has a function that return pointer to FILE
#include "../useful/useful_lib.h"
#include<time.h>  //for TimeProps
#include<math.h>
#include "../gisapi/GisApi.h"
#include "hd5calls.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>
using namespace std;

#include "constant.h"

class FluxProps;
class MatProps;
class NodeHashTable;
class ElementsHashTable;
class Element;
class DischargePlanes;
class TimeProps;
class StatProps;

class TiObject
{
public:
    virtual ~TiObject(){}
    virtual void print0(int spaces=0)=0;
};

//! all scaling
class TiScale:public TiObject
{
public:
    TiScale();
    //! length scaling factor
    double length;

    //! height scaling factor
    //legacy, not used, all height
    //scaling now based on cube root of predicted volume, see below
    double height;

    //! scaling value, ratio of HEIGHT_SCALE to LENGTH_SCALE
    double epsilon;

    //! gravity scaling factor
    double gravity;

    bool auto_calc_height_scale;

    //! cells with flow below this height are neglected for purposes of calculating statistics
    double max_negligible_height;

    //! to get the flow to start moving you have to decrease the friction angles, this variable is used to do that
    double frict_tiny;

    virtual void print0(int spaces=0);

    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="Scale") const
    {
        H5::Group group(parent->createGroup(group_name));
        TiH5_writeDoubleAttribute(group, length);
        TiH5_writeDoubleAttribute(group, height);
        TiH5_writeDoubleAttribute(group, epsilon);
        TiH5_writeDoubleAttribute(group, gravity);
        TiH5_writeBoolAttribute(group, auto_calc_height_scale);
        TiH5_writeDoubleAttribute(group, max_negligible_height);
        TiH5_writeDoubleAttribute(group, frict_tiny);
    }
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="Scale")
    {
        H5::Group group(parent->openGroup(group_name));
        TiH5_readDoubleAttribute(group, length);
        TiH5_readDoubleAttribute(group, height);
        TiH5_readDoubleAttribute(group, epsilon);
        TiH5_readDoubleAttribute(group, gravity);
        TiH5_readBoolAttribute(group, auto_calc_height_scale);
        TiH5_readDoubleAttribute(group, max_negligible_height);
        TiH5_readDoubleAttribute(group, frict_tiny);
    }
};

class TiScalableObject:public TiObject
{
public:
    TiScalableObject(const TiScale &scale__):scale_(scale__),scaled(false){}
    virtual ~TiScalableObject(){}

    //!scaling, return true if scaling is needed
    //!set scaled to true
    virtual bool scale()
    {
        if(scaled)
        {
            return false;
        }
        else
        {
            scaled=true;
            return true;
        }
    }
    //!unscaling return true if unscaling is needed
    //!set scaled to false
    virtual bool unscale()
    {
        if(scaled)
        {
            scaled=false;
            return true;
        }
        else
        {
            return false;
        }
    }
protected:
    const TiScale &scale_;
    bool scaled;
};

inline TiScale::TiScale()
{
    length=0.0;
    height=0.0;
    gravity=0.0;
    max_negligible_height=0.0;
    frict_tiny = 0.1;
    epsilon=1.0;
    auto_calc_height_scale=true;
}
inline void TiScale::print0(int spaces)
{
    printf("%*cScales:\n", spaces,' ');
    printf("%*c  Length: %.5f\n", spaces,' ',length);
    printf("%*c  Gravity: %.5f\n", spaces,' ',gravity);
    printf("%*c  Height: %.5f\n", spaces,' ',height);
    printf("%*c  Height auto-scaled: %s\n", spaces,' ',auto_calc_height_scale?"true":"false");
}


//! LHS stands for Latin Hypercube Sampling, it is a constrained sampling method whose convergence can be much faster than monte carlo 
struct LHS_Props
{
    
    //! the refinement number, the number of times each of the original intervals have been subdivided into 2 intervals in each direction 
    int refnum;

    //! the unique identifier for each sample run 
    int runid;

    //! this constructor initializes refnum and runid to -1
    LHS_Props()
    {
        refnum = runid = -1;
    }
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="LHS_Props") const;
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="LHS_Props");
};



//! the PileProps structure holds the pile properties read in in Read_data() so the pile can be placed at the proper locations shortly thereafter in init_piles() 
class PileProps
{
public:
    //!pile types on change don't forget to update H5::EnumType datatypePileType
    enum PileType {PARABALOID=1, CYLINDER=2,PLANE=3,CASITA=4,POPO=5,ID1=6,ID2=7 };
    //! the number of piles
    int numpiles;

    //!array holding pile height
    std::vector<double> pileheight;

    //!array holding x coordinate of pile center
    std::vector<double> xCen;

    //!array holding y coordinate of pile center
    std::vector<double> yCen;

    //!array holding the major (x before rotation) radius
    std::vector<double> majorrad;

    //!array holding the minor (y before rotation) radius
    std::vector<double> minorrad;

    //!array holding the cosine of the rotation angle
    std::vector<double> cosrot;

    //!array holding the sine of the rotation angle;
    std::vector<double> sinrot;

    //!array holding the initial x speed of the pile
    std::vector<double> initialVx;

    //!array holding the initial y speed of the pile
    std::vector<double> initialVy;

    std::vector<PileProps::PileType> pile_type;

protected:
    double height_scale;
    double length_scale;
    double velocity_scale;
public:
    //! this constuctor initializes the number of piles to zero.
    PileProps();
    virtual ~PileProps();
    //! function allocates space for the pile data
    virtual void allocpiles(int numpiles_in);

    virtual void addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                         double orientation, double Vmagnitude, double Vdirection, PileProps::PileType m_pile_type);
    virtual double get_volume(int i) const
    {
        return PI * pileheight[i] * majorrad[i] * minorrad[i] / 2.0;
    }
    virtual void scale(double length_scale, double height_scale, double gravity_scale);
    double get_smallest_pile_radius();

    /**
     * assign height to point of an elliptical (in (x,y)) shaped pile,
     * the pile can be either parabolic (in the z direction) or be
     * cylindrical (have uniform pile height)
     */
    virtual void set_element_height_to_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *EmTemp, MatProps* matprops);
    virtual double get_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *EmTemp, MatProps* matprops, double* m_xmom=NULL, double* ymom=NULL);

    virtual void print_pile(int i);
    virtual void print0();
    virtual PileProps::PileType get_default_piletype(){return CYLINDER;}

    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="PileProps") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="PileProps");
    //! Create PileProps from hdf file content, will instantiate proper PileProps class
    static PileProps* createPileProps(const H5::CommonFG *parent, const  string group_name="PileProps");

};

//! the PilePropsTwoPhases is PileProps for TwoPhases
class PilePropsTwoPhases: public PileProps
{
public:
    //! array of volume -fractions
    std::vector<double> vol_fract;

    PilePropsTwoPhases();
    virtual ~PilePropsTwoPhases();

    virtual void allocpiles(int numpiles_in);
#ifndef SWIG
    virtual void addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                         double orientation, double Vmagnitude, double Vdirection, PileProps::PileType m_pile_type);
#endif
    virtual void addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                         double orientation, double Vmagnitude, double Vdirection, PileProps::PileType m_pile_type,
                         double volfract);
    virtual void print_pile(int i);
    /**
     * assign height to point of an elliptical (in (x,y)) shaped pile,
     * the pile can be either parabolic (in the z direction) or be
     * cylindrical (have uniform pile height)
     */
    virtual void set_element_height_to_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *EmTemp, MatProps* matprops);
    virtual PileProps::PileType get_default_piletype()
    {
        return PARABALOID;
    }
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="PileProps") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="PileProps");
};

/*************************************************************************/
/* the gis map properties                                                */
/*************************************************************************/

//! this structure holds the path to and name of the GIS map and also a flag to say if there are any extra maps, such as a material properties map, associated with the DEM
class MapNames
{
public:
    //@TODO consolidate with GisApi.h as enums
    static const int GDAL = 2;
    static const int GIS_GRASS = 1;

    //! gis main directory
    std::string gis_main;

    //! gis sub directory
    std::string gis_sub;

    //! gis map set
    std::string gis_mapset;

    //! the actual gis map
    std::string gis_map;

    //! gis_vector
    std::string gis_vector;

    int gis_format;
    //! extra maps: 0=none  +1=2^0=<gis_map>_Mat  MATerial map  +2^(bit#-1)=<gis_map>_Xxx  not yet used
    int extramaps;

    bool region_limits_set;

    double min_location_x;
    double min_location_y;
    double max_location_x;
    double max_location_y;

    MapNames();
    ~MapNames();

    //! this function allocates space for and assigns the information about the GIS map
    void set(const int format, const std::string gis_main_in, const std::string gis_sub_in,
             const std::string gis_mapset_in, const std::string gis_map_in, const std::string gis_vector_in,
             const int extramaps_in);
    void set_region_limits(double min_location_x, double min_location_y, double max_location_x, double max_location_y);
    void print0();
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="MapNames") const;
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="MapNames");
};

/**************************************************************************/
/* the time properties                                                    */
/**************************************************************************/

//! this structure holds all the information about time and timestepping
class TimeProps
{
public:

    //! the maximum # of iterations (a.k.a. time steps) before the simulation ends
    //! if negative don't consider maxiter for stopping
    int maxiter;

    //! the current number of iterations
    int iter;

    //! the maximum amount of time (in seconds) before the simulation ends
    //! if negative don't consider maxtime for stopping
    double maxtime;

    //! the non-dimensional maxtime
    //! if negative don't consider ndmaxtime for stopping
    double ndmaxtime;

    //! the amount of time (in seconds) between subsequent outputs (with one exception... when the simulation ends one final output is performed and that one will be less than "timeoutput" seconds after the previous one
    double dTimeVizOutSaving;

    //!frequency of time series for visualization file savings in iterations, if negative don't consider this criteria for saving
    int dIterVizOutSaving;

    //! the amount of time (in seconds) between subsequent saves (with one exception... when the simulation ends one final save is performed and that one will be less than "timeoutput" seconds after the previous one
    double dTimeRestartSaving;

    //!frequency of restart file savings in iterations, if negative don't consider this criteria for saving
    int dIterRestartSaving;

    //! count of the number of times output has been done
    int iVizOutSaved;

    //! count of the number of times saving has been done
    int iRestartSaved;

    //! the non-dimensional time at which the next output should occur
    double ndNextVizOutput;

    //! the non-dimensional time at which the next save should occur
    double ndNextRestartSave;

    //! the value used to nondimensionalize time
    double TIME_SCALE;

    //! the non-dimensional current time
    double cur_time;

    //! the non-dimensional time step
    double dtime;

    //! velocity measure/Vslump used as a STOPPING CRITERIA,which is while it's stored under TimeProps (also members of MatProps are assigned once, are permanent for the run), see MatProps struct below for Vslump, see ../main/datread.C for initialization of Vslump
    double vstarmax;

    //! wallclock time shortly after titan starts running
    time_t starttime;

    TimeProps()
    {
        starttime = time(NULL);
        TIME_SCALE = 1.0;
        setTime(0, 0.0);
        setRestartOutput(-1,-1.0);
        setTimeSeriesOutput(-1,-1.0);
    }
    ~TimeProps()
    {

    }
    //! this function initializes the time properties at the beginning of the simulation
    void setTime(int maxiterin, double maxtimein)
    {
        maxiter = maxiterin;
        maxtime = maxtimein;
        ndmaxtime = maxtime / TIME_SCALE;
        iter = 0;
        cur_time = 0.0;
        dtime = 0.0;
        vstarmax = 0.0;
    }
    //! frequency of restart files output
    void setRestartOutput(int dIterOut, double dTimeOut)
    {
        dTimeRestartSaving = dTimeOut;
        dIterRestartSaving=dIterOut;
        ndNextRestartSave = dTimeRestartSaving / TIME_SCALE;
        iRestartSaved = 0;
    }
    //! frequency of time series output
    void setTimeSeriesOutput(int dIterOut, double dTimeOut)
    {
        dTimeVizOutSaving = dTimeOut;
        dIterVizOutSaving=dIterOut;
        ndNextVizOutput = dTimeVizOutSaving / TIME_SCALE;
        iVizOutSaved = 0;
    }

    void scale(double TIME_SCALEin)
    {
        TIME_SCALE = TIME_SCALEin;
        ndmaxtime = maxtime / TIME_SCALE;
        ndNextVizOutput = dTimeVizOutSaving / TIME_SCALE;
        ndNextRestartSave = dTimeRestartSaving / TIME_SCALE;
    }

    //! this function increments the time step, after, if neccessary, decreasing the time step to land evenly on the next time to output save or end the simulation
    void incrtime(double *dt)
    {
        // first reduce dt to hit output or end time "exactly"
        if(dTimeVizOutSaving>0.0)
            if(cur_time + *dt > ndNextVizOutput)
                *dt = ndNextVizOutput - cur_time;
        if(dTimeRestartSaving>0.0)
            if(cur_time + *dt > ndNextRestartSave)
                *dt = ndNextRestartSave - cur_time;
        if(maxtime>=0.0)
            if(cur_time + *dt > ndmaxtime)
                *dt = ndmaxtime - cur_time;
        dtime = *dt;
        // then increment time
        cur_time += *dt;
        iter++;
    }

    int ifstart()
    {
        return (iter == 0);
    } //! checks if it's before first time step
    int iffirst()
    {
        return (iter == 1);
    } //! checks if it's at first time step

    //! checks if the simulation should end due to running out of time, running out of timesteps or meeting a legacy "stopped" criteria
    int ifend(double vstar)
    {
        if(vstar > vstarmax)
            vstarmax = vstar;
        if(ndmaxtime>=0.0)
            if(cur_time >= ndmaxtime)
                return 1;
        if(maxiter>=0)
            if(iter > maxiter)
                return 1;
        if((vstarmax > 2.0) && !(vstar > 1.0))
            return 1;
        return 0;
    }

    //! checks if the simulation has passed 1/10th of the maximum time allowed
    int ifcheckstop()
    {
        return (cur_time > ndmaxtime / 10.0);
    }

    //! checks if the restart file should be saved now
    int ifTimeForRestartOutput()
    {
        if(dTimeRestartSaving>0.0)
        {
            if(cur_time >= ndNextRestartSave)
            {
                iRestartSaved++; //using isave eliminates roundoff
                ndNextRestartSave = ((iRestartSaved + 1) * dTimeRestartSaving) / TIME_SCALE;
                return 1;
            }
        }
        if(dIterRestartSaving>0)
        {
            if(iter%dIterRestartSaving==0)
                return 1;
        }
        return 0;
    }

    //! checks if the output files should be written now
    int ifTimeForTimeSeriesOutput()
    {
        if(dTimeVizOutSaving>0.0)
        {
            if(cur_time >= ndNextVizOutput)
            {
                iVizOutSaved++; //using ioutput eliminates roundoff
                ndNextVizOutput = ((iVizOutSaved + 1) * dTimeVizOutSaving) / TIME_SCALE;
                return (1);
            }
        }
        if(dIterVizOutSaving>0)
        {
            if(iter%dIterVizOutSaving==0)
                return 1;
        }
        return 0;
    }

    //! chunk simulated time into hours minutes and seconds
    void chunktime(int *hours, int *minutes, double *seconds)
    {
        double dimtime = cur_time * TIME_SCALE;
        *hours = ((int) dimtime) / 3600;
        *minutes = (((int) dimtime) % 3600) / 60;
        *seconds = dimtime - (double) (*hours * 3600 + *minutes * 60);
    }

    //! return the simulated time in seconds
    double timesec()
    {
        return (cur_time * TIME_SCALE);
    }
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="TimeProps") const
    {
        H5::Group group(parent->createGroup(group_name));
        TiH5_writeIntAttribute(group, maxiter);
        TiH5_writeIntAttribute(group, iter);
        TiH5_writeDoubleAttribute(group, maxtime);
        TiH5_writeDoubleAttribute(group, ndmaxtime);
        TiH5_writeDoubleAttribute(group, dTimeVizOutSaving);
        TiH5_writeIntAttribute(group, dIterVizOutSaving);
        TiH5_writeDoubleAttribute(group, dTimeRestartSaving);
        TiH5_writeIntAttribute(group, dIterRestartSaving);
        TiH5_writeIntAttribute(group, iVizOutSaved);
        TiH5_writeIntAttribute(group, iRestartSaved);
        TiH5_writeDoubleAttribute(group, ndNextVizOutput);
        TiH5_writeDoubleAttribute(group, ndNextRestartSave);
        TiH5_writeDoubleAttribute(group, TIME_SCALE);
        TiH5_writeDoubleAttribute(group, cur_time);
        TiH5_writeDoubleAttribute(group, dtime);
        TiH5_writeDoubleAttribute(group, vstarmax);
        //time_t starttime;
    }
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="TimeProps")
    {
        H5::Group group(parent->openGroup(group_name));
        TiH5_readIntAttribute(group, maxiter);
        TiH5_readIntAttribute(group, iter);
        TiH5_readDoubleAttribute(group, maxtime);
        TiH5_readDoubleAttribute(group, ndmaxtime);
        TiH5_readDoubleAttribute(group, dTimeVizOutSaving);
        TiH5_readIntAttribute(group, dIterRestartSaving);
        TiH5_readDoubleAttribute(group, dTimeRestartSaving);
        TiH5_readIntAttribute(group, dIterVizOutSaving);
        TiH5_readIntAttribute(group, iVizOutSaved);
        TiH5_readIntAttribute(group, iRestartSaved);
        TiH5_readDoubleAttribute(group, ndNextVizOutput);
        TiH5_readDoubleAttribute(group, ndNextRestartSave);
        TiH5_readDoubleAttribute(group, TIME_SCALE);
        TiH5_readDoubleAttribute(group, cur_time);
        TiH5_readDoubleAttribute(group, dtime);
        TiH5_readDoubleAttribute(group, vstarmax);
        //time_t starttime;
    }
};

/*****************************************************************************/
//! this struct holds constants for material properties as well as other constants note that the material id tag (used as the indice for material properties... matname, bedfrict) as returned by Get_raster_id() (a GIS function call) starts from 1 and not from 0 so arrays must be one element larger
/*****************************************************************************/
class MatProps
{
public:
    //! the "maximum" number of cells across the smallest pile/flux-source minor axis
    int number_of_cells_across_axis;

    //! the smallest minor radius
    double smallest_axis;

    //! the number of different materials on the map
    int material_count;

    //! the names of each material
    std::vector<std::string> matnames;

    //! phi_{int}, the internal friction angle (must be GREATER than the bedfriction angle)
    //double intfrict;

    //! tan(phi_{int}), tangent of the internal friction angle
    //double tanintfrict;

    //! phi_{bed}, the bed friction angle, must be LESS than the internal friction angle and should be greater than about 8 degrees, this minimum angle may change once Keith's local stopping criteria is enforced
    std::vector<double> bedfrict;

    //! tan(phi_{bed}), tangent of the bed friction angle
    std::vector<double> tanbedfrict;

    //! v_f, legacy not used
    double porosity;

    //! pore fluid viscosity, legacy not used
    double mu; //is it same as viscosity?

    //! density
    double rho;

    //! scaling value, ratio of HEIGHT_SCALE to LENGTH_SCALE
    //double epsilon;

    //! slope limiting stuff
    double gamma;

    TiScale &scale;

    //! Used in Bin Yu's legacy global stopping criteria, this never worked right except for initially cylindrical piles on a horizontal plane
    //! that were released and slumped... there Bin Yu's global stopping criteria worked great.
    double Vslump;

    //! to get the flow to start moving you have to decrease the friction angles, this variable is used to do that
    //double frict_tiny;

    //! this constructor allocates initial properties unfortunately the properties aren't known at the time this is called so dummy values are fed in instead for many if not all of these
    MatProps(TiScale &_scale):
        scale(_scale)
    {
        number_of_cells_across_axis = 0;
        smallest_axis = 0.0;
        material_count = 0;
        porosity = 1.0;

        mu = 0.0001;
        rho = 2200.0;

        //epsilon = 1.0;
        gamma = 1.0;
        Vslump = 0.0;
        //frict_tiny = 0.1;
        
        //something somewhere counting from 1, so add dummy values
        matnames.push_back("nothing");
        bedfrict.push_back(12.0);
        tanbedfrict.push_back(12.0);
    }
    //! this destructor deallocates the arrays of bed friction angles and their tangents
    virtual ~MatProps()
    {
    }
    virtual void process_input()
    {
        int imat;
        //intfrict = intfrict * PI / 180.0;
        //tanintfrict = tan(intfrict);

        tanbedfrict.resize(bedfrict.size());
        for(imat = 1; imat <= material_count; imat++)
        {
            bedfrict[imat] *= PI / 180.0;
            tanbedfrict[imat] = tan(bedfrict[imat]);
        }
    }

    virtual inline void set_scale(const PileProps *pileprops_ptr = NULL, const FluxProps *fluxprops_ptr = NULL);
    virtual inline void calc_Vslump(const PileProps *pileprops_ptr, const FluxProps *fluxprops_ptr);
    double get_TIME_SCALE()
    {
        return sqrt(scale.length / scale.gravity);
    }
    //non-dimensionalize the inputs
    double get_VELOCITY_SCALE()
    {
        return sqrt(scale.length * scale.gravity);
    }

    virtual void print0()
    {
        int i;
        printf("Material properties:\n");
        printf("\tNumber of cells across axis: %d\n", number_of_cells_across_axis);
        //printf("\tInternal friction: %f\n", intfrict * 180.0 / PI);
        printf("\tBed friction:\n");
        for(i = 1; i <= material_count; i++)
        {
            printf("\t\t%d %s %f\n", i, matnames[i].c_str(), bedfrict[i] * 180.0 / PI);
        }
    }
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="MatProps") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="MatProps");
    //! Create MatProps from hdf file content, will instantiate proper MatProps class
    static MatProps* createMatProps(const H5::CommonFG *parent,TiScale &_scale, const  string group_name="MatProps");


};

class MatPropsTwoPhases: public MatProps
{
public:
    MatPropsTwoPhases(TiScale &_scale) :
            MatProps(_scale)
    {
        mu = 0.1;
        rho = 2700;
        den_solid = 2700.0;
        den_fluid = 1200.0;
        viscosity = 0.1;
        v_terminal = 0.0;

        flow_type = 0;
    }
    virtual ~MatPropsTwoPhases()
    {
    }

    //! density
    double den_solid;

    //! fluid density
    double den_fluid;

    //! fluid viscosity
    double viscosity;

    //! terminal velocity of a single solid particle in fluid medium
    double v_terminal;

    //! Flow type flag
    int flow_type;

    virtual void process_input()
    {
        MatProps::process_input();

        double diameter = 0.005;
        v_terminal = pow(diameter, 2.) * (den_solid - den_fluid) * scale.gravity / (18. * viscosity);
    }
    virtual inline void set_scale(const PileProps *pileprops_ptr = NULL, const FluxProps *fluxprops_ptr = NULL);
    virtual inline void calc_Vslump(const PileProps *pileprops_ptr, const FluxProps *fluxprops_ptr);
    //! Dump object content to hdf5 file
    virtual void h5write(H5::CommonFG *parent, string group_name="MatProps") const;
    //! Load object content from hdf5 file
    virtual void h5read(const H5::CommonFG *parent, const  string group_name="MatProps");
};


/** the OutLine Structure holds the maximum throughout time flow depth at every spatial point
 *
 * first it collects max in pileheight_by_elm which is per element based (update_single_phase/update_two_phases)
 * when geometry is changed it converts pileheight_by_elm to pileheight_loc, which is a grid local to each thread (flush_stats)
 * at that step it also calculated the boundaries of each element (update_on_changed_geometry)
 * finally then the time came it combines all results from threads it single final grid pileheight (combine_results_from_threads)
 */
class OutLine
{
public:
    enum OutLineInitSize {AMR=1, DEM=2};

    //geometric id of temporary arrays
    int conformation;
    // the least squares interpolated height didn't work for some reason
    // that is didn't produce a meaningful pileheight contour map so it has
    // been commented out.  that is everything to deal with pileheight2
//! number of cells in the x direction on the map
    int Nx;

    //! number of cells in the y direction on the map 
    int Ny;

    //!2d elements accessed as a[ix+iy*stride]
    int stride;

    int size;


    //! length of a cell in the x direction
    double dx;

    //! length of a cell in the y direction
    double dy;

    //! min and max x coordinate on the map
    double xminmax[2];

    //! min and max y coordinate on the map
    double yminmax[2];

    //! 2-d array holding the maximum throughout time pileheight at every point
    //! pileheight[ithread][ix+iy*stride]
    double *pileheight;

    //! 2-d array holding the maximum throughout time kinetice energy at every point
    double *max_kinergy;

    //! 2-d array holding the maximum throughout time kinetice energy at every point
    double *max_dynamic_pressure;

    //! 2-d array holding the cummulative kinetic energy at every point
    double *cum_kinergy;

    ElementType elementType;

    //! do the calculations
    bool enabled;

    OutLineInitSize init_size;
    int max_linear_size;

    TiScale *scale;

    std::string output_prefix;

    //! this is the OutLine constructor it initializes the number of cells to zero
    OutLine();
    
    //! this is the OutLine it deallocates the 2 dimensional array holding maximum throughout time pileheight in every cell on the map
    ~OutLine();
    
    void setElemNodeTable(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable);

    //! this function initializes the OutLine map/2-dimensional array 
    void init(const double *dxy, int power, double *XRange, double *YRange);
    
    //! this function initializes the OutLine map/2-dimensional array
    void init_DEM_resolution(double resx, double resy, double *XRange, double *YRange);

    //! this function reinitializes the OutLine map/2-dimensional array during restart
    void init2(const double *dxy, double *XRange, double *YRange);

    //! this function updates the maximum throughout time pileheight in every cell covered by an arbitrary element
    void update();

    /*! this function outputs the maximum over time map of pileheights
     *  to the file pileheightrecord.xxxxxx
     */
    void output(MatProps* matprops_ptr, StatProps* statprops_ptr);
    
    //! this function reads in the previous map of maximum throughout time pileheight stored in the file pileheightrecord.xxxxxx during restart
    void reload(MatProps* matprops_ptr, StatProps* statprops_ptr);
    
    void combine_results_from_threads();

    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="OutLine");
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="OutLine");

protected:
    //!allocate internal arrays for new size
    void alloc_arrays();

    void update_on_changed_geometry();
    void flush_stats(bool zero_old_arrays=true);

    void update_single_phase();
    void update_two_phases();

    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;

    int myid;
    int numprocs;



    //! 2-d array holding the maximum throughout time pileheight at every point
    //! pileheight_loc[ithread][ix+iy*stride]
    double **pileheight_loc;
    vector<double, AlignmentAllocator<double> > pileheight_by_elm;

    //! 2-d array holding the maximum throughout time kinetice energy at every point
    double **max_kinergy_loc;
    vector<double, AlignmentAllocator<double> > max_kinergy_by_elm;

    //! 2-d array holding the maximum throughout time kinetice energy at every point
    double **max_dynamic_pressure_loc;
    vector<double, AlignmentAllocator<double> > max_dynamic_pressure_by_elm;

    //! 2-d array holding the cummulative kinetic energy at every point
    double **cum_kinergy_loc;
    vector<double, AlignmentAllocator<double> > cum_kinergy_by_elm;



    vector<int> el_x_start;
    vector<int> el_x_stop;
    vector<int> el_y_start;
    vector<int> el_y_stop;
};

/**********************************************************************/
//! this structure is for the calculation of volume that flows through user specified discharge planes.  The sign of the volume indicates which direction the flow went and follows the right hand rule convention.  velocity cross (point b-point a) is the sign of the flow through planes.  This means if you surround the only pile, specify the points defining the discharge planes in counter clockwise order, the flow "out of the box" will be positive.  if you specify the points in clockwise order flow "out of the box" will be negative.
class DischargePlanes
{
public:
    //! the number of discharge planes
    int num_planes;

    //! the discharge planes are lines (with planes normal to the surface intersecting the surface passing through the lines), this holds a lot of information associated with each planes, a lot of precomputed quantities to make updating the flux through the planes fast.
    //double **planes;
    std::vector<std::array<double,10> > planes;

    //! this constructor initializes the number of planes to zero
    DischargePlanes()
    {
        num_planes = 0;
        return;
    }
    
    //! this destructor deallocates the planes information
    ~DischargePlanes()
    {
        return;
    }
    void allocate(int m_num_planes)
    {
        num_planes = m_num_planes;
        planes.resize(num_planes);
    }
private:
    void calculateDerivativeProps(int i)
    {
        planes[i][4] = planes[i][1] - planes[i][0]; //xb-xa
        planes[i][5] = planes[i][3] - planes[i][2]; //yb-ya
        planes[i][6] = //(xb-xa)^2+(yb-ya)^2
                planes[i][4] * planes[i][4] + planes[i][5] * planes[i][5];
        planes[i][7] = //ya*(xb-xa)-xa*(yb-ya)
                planes[i][3] * planes[i][4] - planes[i][1] * planes[i][5];
        planes[i][8] = ((fabs(planes[i][4]) + fabs(planes[i][5])) * (fabs(planes[i][4]) + fabs(planes[i][5])))
                / planes[i][6];
        planes[i][9] = 0.0; //discharge through the planes[i]
    }
public:
    //! this function add plane and initializes the planes information (to zero flux through the planes) and precomputes a number of quantities to make updating the flux through the planes fast

    void addDischargePlane(const double m_x_a, const double m_y_a, const double m_x_b, const double m_y_b)
    {
        std::array<double,10> plane;

        plane[0] = m_x_a; //xa
        plane[1] = m_y_a; //xb
        plane[2] = m_x_b; //ya
        plane[3] = m_y_b; //yb

        //printf("plane %d: (%16.10g,%16.10g) (%16.10g,%16.10g)\n",iplane,planes[iplane][0],planes[iplane][2],planes[iplane][1],planes[iplane][3]);

        planes.push_back(plane);
        num_planes = planes.size();
        calculateDerivativeProps(num_planes - 1);
    }
    
    //reinitialized in load_run()
    void init(int num_planes_in, double **planes_in)
    {
        double planestemp[4];
        
        //printf("num_planes_in=%d   num_planes=%d\n",num_planes_in,num_planes);
        allocate(0);
        //printf("num_planes_in=%d   num_planes=%d\n",num_planes_in,num_planes);
        if(num_planes_in > 0)
        {
            for(int iplane = 0; iplane < num_planes_in; iplane++)
            {
                addDischargePlane(planes_in[iplane][0], planes_in[iplane][1], planes_in[iplane][2],
                                  planes_in[iplane][3]);
            }
        }
        return;
    }

    void scale(double length_scale)
    {
        for(int i = 0; i < num_planes; i++)
        {
            planes[i][0] /= length_scale;
            planes[i][1] /= length_scale;
            planes[i][2] /= length_scale;
            planes[i][3] /= length_scale;
            calculateDerivativeProps(i);
        }
    }

    void print_discharge_plane(int i)
    {
        printf("\tDischarge plane %d:\n", i);
        printf("\t\tPoint A (UTM E, UTM N): %f, %f\n", planes[i][0], planes[i][1]);
        printf("\t\tPoint B (UTM E, UTM N): %f, %f\n", planes[i][2], planes[i][3]);
    }

    void print0()
    {
        int i;
        if(num_planes > 0)
        {
            printf("Discharge planes:    (Number of discharge planes: %d)\n", num_planes);
            for(i = 0; i < num_planes; i++)
                print_discharge_plane(i);
        }
        else
        {
            printf("Discharge planes:    there is no discharge planes\n");
        }
    }
    
    //! this function updates the flux through the discharge planes
    void update(const double nodes[9][2], const double hVx, const double hVy, const double dt)
    {
        //FILE* fp=fopen("dischargedebug","a");
        
        double doubleswap1, doubleswap2, doubleswap3;
        double nearestpoint[2], intersectpoint[2][2];
        int iplane, iintersect, icorner1, icorner2;
        double dxnode[4], dynode[4];
        double dist2nearest2;
        double halfsidelength = //assumes grid is not slanted
                (fabs(nodes[6][0] - nodes[4][0]) + fabs(nodes[6][1] - nodes[4][1]) + fabs(nodes[7][0] - nodes[5][0])
                 + fabs(nodes[7][1] - nodes[5][1]))
                * 0.25;
        double halfsidelength2 = halfsidelength * halfsidelength;
        double err = halfsidelength / pow(2.0, 19.0);
        
        for(icorner1 = 0; icorner1 < 4; icorner1++)
        {
            icorner2 = (icorner1 + 1) % 4;
            dxnode[icorner1] = nodes[icorner2][0] - nodes[icorner1][0];
            dynode[icorner1] = nodes[icorner2][1] - nodes[icorner1][1];
        }
        
        for(iplane = 0; iplane < num_planes; iplane++)
        {
            /* nearest x & y along line (not necessaryly on line segment) a and b
             and the elements center node
             x=(-(yb-ya)*(ya*(xb-xa)-xa*(yb-ya))+(xb-xa)*(y8*(yb-ya)+x8*(xb-xa)))/
             ((xb-xa)^2+(yb-ya)^2);
             y=( (xb-xa)*(ya*(xb-xa)-xa*(yb-ya))+(yb-ya)*(y8*(yb-ya)+x8*(xb-xa)))/
             ((xb-xa)^2+(yb-ya)^2);
             */
            doubleswap1 = nodes[8][1] * planes[iplane][5] + nodes[8][0] * planes[iplane][4];
            //         =y8         *(yb-ya)          +x8         *(xb-xa)
            nearestpoint[0] = (-planes[iplane][5] * planes[iplane][7] + planes[iplane][4] * doubleswap1)
                    / planes[iplane][6];
            nearestpoint[1] = (+planes[iplane][4] * planes[iplane][7] + planes[iplane][5] * doubleswap1)
                    / planes[iplane][6];
            
            dist2nearest2 = (nodes[8][0] - nearestpoint[0]) * (nodes[8][0] - nearestpoint[0])
                    + (nodes[8][1] - nearestpoint[1]) * (nodes[8][1] - nearestpoint[1]);
            
            //check if line interesects with cell (not counting a single corner)
            if(dist2nearest2 < halfsidelength2 * planes[iplane][8])
            {
                /* fprintf(fp,"nodecoord=(%8f,%8f) nearestpoint=(%8f,%8f)",
                 nodes[8][0],nodes[8][1],nearestpoint[0],nearestpoint[1]);
                 */

                //it does intersect... if not terminated by endpoints
                iintersect = 0; //stop when you have 2 unique "intersections" 
                //(the end of a discharge plane line segment is considered to 
                //be an "intersection")
                for(icorner1 = 0; icorner1 < 4; icorner1++)
                {
                    doubleswap1 = nodes[icorner1][1] * dxnode[icorner1] - nodes[icorner1][0] * dynode[icorner1];
                    doubleswap2 = dxnode[icorner1] * planes[iplane][5] - dynode[icorner1] * planes[iplane][4];
                    
                    if((doubleswap2 < 0) || (doubleswap2 > 0))
                    {
                        
                        intersectpoint[iintersect][0] = (planes[iplane][4] * doubleswap1
                                - dxnode[icorner1] * planes[iplane][7])
                                                        / doubleswap2;
                        intersectpoint[iintersect][1] = (planes[iplane][5] * doubleswap1
                                - dynode[icorner1] * planes[iplane][7])
                                                        / doubleswap2;
                        
                        int ifprint = 0;
                        if(((intersectpoint[iintersect][0] < planes[iplane][0] - err) && (planes[iplane][0]
                                < planes[iplane][1]))
                           || ((intersectpoint[iintersect][0] > planes[iplane][0] + err) && (planes[iplane][0]
                                   > planes[iplane][1]))
                           || ((intersectpoint[iintersect][1] < planes[iplane][2] - err) && (planes[iplane][2]
                                   < planes[iplane][3]))
                           || ((intersectpoint[iintersect][1] > planes[iplane][2] + err) && (planes[iplane][2]
                                   > planes[iplane][3])))
                        {
                            intersectpoint[iintersect][0] = planes[iplane][0];
                            intersectpoint[iintersect][1] = planes[iplane][2];
                        }
                        else if(((intersectpoint[iintersect][0] < planes[iplane][1] - err)
                                && (planes[iplane][1] < planes[iplane][0]))
                                || ((intersectpoint[iintersect][0] > planes[iplane][1] + err) && (planes[iplane][1]
                                        > planes[iplane][0]))
                                || ((intersectpoint[iintersect][1] < planes[iplane][3] - err) && (planes[iplane][3]
                                        < planes[iplane][2]))
                                || ((intersectpoint[iintersect][1] > planes[iplane][3] + err) && (planes[iplane][3]
                                        > planes[iplane][2])))
                        {
                            intersectpoint[iintersect][0] = planes[iplane][1];
                            intersectpoint[iintersect][1] = planes[iplane][3];
                        }
                        
                        if(((((intersectpoint[iintersect][0] - intersectpoint[0][0]) * (intersectpoint[iintersect][0]
                                - intersectpoint[0][0])
                              + (intersectpoint[iintersect][1] - intersectpoint[0][1]) * (intersectpoint[iintersect][1]
                                      - intersectpoint[0][1]))
                             > (halfsidelength / 512.0))
                            || (iintersect == 0))
                           && (doubleswap2 != 0.0))
                            iintersect++;
                        
                        if(iintersect == 2)
                        {
                            
                            if((intersectpoint[1][0] - intersectpoint[0][0]) * (planes[iplane][1] - planes[iplane][0]) + (intersectpoint[1][1]
                                    - intersectpoint[0][1])
                                                                                                                         * (planes[iplane][3] - planes[iplane][2])
                               < 0.0)
                            {
                                doubleswap3 = intersectpoint[0][0];
                                intersectpoint[0][0] = intersectpoint[1][0];
                                intersectpoint[1][0] = doubleswap3;
                                doubleswap3 = intersectpoint[0][1];
                                intersectpoint[0][1] = intersectpoint[1][1];
                                intersectpoint[1][1] = doubleswap3;
                            }
                            
                            break; //we've got both intersection points
                        }
                    }
                } // for(icorner1=0;icorner1<4;icorner1++)	 
                
                if(iintersect == 1)
                {
                    //a plane end point is in this cell
                    
                    doubleswap1 = //dist from center node to 1st plane endpoint
                            (nodes[8][0] - planes[iplane][0]) * (nodes[8][0] - planes[iplane][0]) + (nodes[8][1]
                                    - planes[iplane][2])
                                                                                                    * (nodes[8][1] - planes[iplane][2]);
                    
                    doubleswap2 = //dist from center node to 2nd plane endpoint
                            (nodes[8][0] - planes[iplane][1]) * (nodes[8][0] - planes[iplane][1]) + (nodes[8][1]
                                    - planes[iplane][3])
                                                                                                    * (nodes[8][1] - planes[iplane][3]);
                    
                    if(doubleswap1 <= doubleswap2)
                    { //it's the 1st end point
                        intersectpoint[1][0] = planes[iplane][0];
                        intersectpoint[1][1] = planes[iplane][2];
                    }
                    else
                    { //it's the 2nd end point
                        intersectpoint[1][0] = planes[iplane][1];
                        intersectpoint[1][1] = planes[iplane][3];
                    }
                    
                    iintersect = 2; //we consider plane end point to be an "intersection"
                }
                
                if(iintersect == 2)
                {
                    
                    //discharge += dt*(hVx*dy-hVy*dx)
                    planes[iplane][9] += dt
                            * (hVx * (intersectpoint[1][1] - intersectpoint[0][1]) - hVy
                                    * (intersectpoint[1][0] - intersectpoint[0][0]));
                } // if(iintersect==2)
                
            } //if(halfsidelength*(absdy+absdx)>=absdx*absdx+absdy*absdy)
            
        } //for(iplane=0;iplane<num_planes;iplane++)
        
        return;
    }
    
    //! this function deallocates the array holding the information about the discharge planes
    void dealloc()
    {
        num_planes = 0;
        planes.resize(0);
        return;
    }
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="DischargePlanes") const;
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="DischargePlanes");
};

//! The FluxProps Structure holds all the data about extrusion flux sources (material flowing out of the ground) they can become active and later deactivate at any time during the simulation.  There must be at least 1 initial pile or one flux source that is active at time zero, otherwise the timestep will be set to zero and the simulation will never advance.
class FluxProps
{
public:
    //! number of extrusion flux sources
    int no_of_sources;

    //! array holding the influx rates for all extrusion flux sources, if calculation is not scaled this quantity has units of [m/s]
    std::vector<double> influx;

    //! array holding the activation time for all extrusion flux sources
    std::vector<double> start_time;

    //! array holding the deactivation time for all extrusion flux sources
    std::vector<double> end_time;

    //! array holding the x coordinate of the center of each the extrusion flux source
    std::vector<double> xCen;

    //! array holding the y coordinate of the center of each the extrusion flux source
    std::vector<double> yCen;

    //! array holding the radius along the major axis of each the extrusion flux source, the 2D pile shape is elliptical, in 3D it's a paraboloid
    std::vector<double> majorrad;

    //! array holding the radius along the minor axis of each the extrusion flux source, the 2D pile shape is elliptical, in 3D it's a paraboloid
    std::vector<double> minorrad;

    //! array holding the cosine of the rotation angle of each the extrusion flux source, the 2D pile shape is elliptical, and the ellipse's major axis does not have to be aligned with the x axis.
    std::vector<double> cosrot;

    //! array holding the sine of the rotation angle of each the extrusion flux source, the 2D pile shape is elliptical, and the ellipse's major axis does not have to be aligned with the x axis.
    std::vector<double> sinrot;

    //! array holding the x component of velocity of material extruding from each flux source
    std::vector<double> xVel;

    //! array holding the y component of velocity of material extruding from each flux source
    std::vector<double> yVel;

protected:
    double height_scale;
    double length_scale;
    double velocity_scale;
    double time_scale;

public:
    FluxProps()
    {
        no_of_sources = 0;
    }
    ~FluxProps()
    {
    }
    //! this function allocates space for all the extrusion rate fluxes, Dinesh Kumar added it, Keith modified it slightly
    void allocsrcs(int nsrcs)
    {
        no_of_sources = nsrcs;
        influx.resize(nsrcs);
        xCen.resize(nsrcs);
        yCen.resize(nsrcs);
        majorrad.resize(nsrcs);
        minorrad.resize(nsrcs);
        cosrot.resize(nsrcs);
        sinrot.resize(nsrcs);
        start_time.resize(nsrcs);
        end_time.resize(nsrcs);
        xVel.resize(nsrcs);
        yVel.resize(nsrcs);
    }
    //! add flux source
    virtual void addFluxSource(double m_influx, double m_start_time, double m_end_time, double xcenter, double ycenter,
                               double majradius, double minradius, double orientation, double Vmagnitude,
                               double Vdirection)
    {
        no_of_sources++;
        influx.push_back(m_influx);
        xCen.push_back(xcenter);
        yCen.push_back(ycenter);
        majorrad.push_back(majradius);
        minorrad.push_back(minradius);
        cosrot.push_back(cos(orientation * PI / 180.0));
        sinrot.push_back(sin(orientation * PI / 180.0));
        start_time.push_back(m_start_time);
        end_time.push_back(m_end_time);
        xVel.push_back(Vmagnitude * cos(Vdirection * PI / 180.0));
        yVel.push_back(Vmagnitude * sin(Vdirection * PI / 180.0));
    }
    virtual void scale(double m_length_scale, double m_height_scale, double m_gravity_scale)
    {
        length_scale = m_length_scale;
        height_scale = m_height_scale;
        //non-dimensionalize the inputs
        time_scale = sqrt(length_scale / m_gravity_scale);
        velocity_scale = sqrt(length_scale * m_gravity_scale);
        int isrc;
        for(isrc = 0; isrc < no_of_sources; isrc++)
        {
            influx[isrc] *= time_scale / height_scale;
            start_time[isrc] /= time_scale;
            end_time[isrc] /= time_scale;
            xCen[isrc] /= length_scale;
            yCen[isrc] /= length_scale;
            majorrad[isrc] /= length_scale;
            minorrad[isrc] /= length_scale;
            xVel[isrc] /= velocity_scale;
            yVel[isrc] /= velocity_scale;
        }
    }
    double get_smallest_source_radius()
    {
        double smallest_source_radius = HUGE_VAL;
        int isrc;
        for(isrc = 0; isrc < no_of_sources; isrc++)
        {
            if(smallest_source_radius > majorrad[isrc])
                smallest_source_radius = majorrad[isrc];

            if(smallest_source_radius > minorrad[isrc])
                smallest_source_radius = minorrad[isrc];
        }
        return smallest_source_radius;
    }
    double get_effective_height(int i)
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
        double t3 = X * X * influx[i] / g;
        if(t3 > (end_time[i] - start_time[i]))
            t3 = (end_time[i] - start_time[i]);
        return influx[i] * t3 - 0.5 * a * t3 * t3;
    }
    double get_volume(int isrc) const
    {
        return 0.5 * PI * influx[isrc] * majorrad[isrc] * minorrad[isrc] * 0.5 * (end_time[isrc] - //0.5 for linear decrease
                start_time[isrc]);
    }
    virtual void print_source(int i)
    {
        printf("\tFlux_source %d:\n", i);

        printf("\t\tExtrusion flux rate [m/s]:%f\n", influx[i] * height_scale/time_scale );
        printf("\t\tActive Time [s], start, end: %f %f\n", start_time[i], end_time[i]);

        printf("\t\tCenter of Initial Volume, xc, yc (UTM E, UTM N): %f %f\n", xCen[i] * length_scale, yCen[i] * length_scale);
        printf("\t\tMajor and Minor Extent, majorR, minorR (m, m): %f %f\n", majorrad[i] * length_scale, minorrad[i] * length_scale);
        double orientation = atan2(sinrot[i], cosrot[i]) * 180.0 / PI;
        printf("\t\tOrientation (angle [degrees] from X axis to major axis): %f\n", orientation);
        double Vmagnitude = sqrt(xVel[i] * xVel[i] + yVel[i] * yVel[i]);
        double Vdirection = atan2(yVel[i], xVel[i]) * 180.0 / PI;
        printf("\t\tInitial speed [m/s]: %f\n", Vmagnitude * velocity_scale);
        printf("\t\tInitial direction ([degrees] from X axis): %f\n", Vdirection);
        printf("\t\tEffective Thickness, P (m):%f\n", get_effective_height(i) * height_scale);
    }
    virtual void print0()
    {
        int i;
        if(no_of_sources > 0)
        {
            printf("Flux sources:    (Number of flux sources: %d)\n", no_of_sources);
            for(i = 0; i < no_of_sources; i++)
                print_source(i);
        }
        else
        {
            printf("Flux sources:    there is no flux sources\n");
        }
    }
    //! this function returns 1 if any flux sources become active during the current timestep, this is used to trigger "initial adaptation" of the flux source area, Keith wrote this function
    int IfAnyStart(TimeProps *timeprops_ptr)
    {

        for(int isrc = 0; isrc < no_of_sources; isrc++)
            if(((timeprops_ptr->cur_time - timeprops_ptr->dtime <= start_time[isrc]) && (start_time[isrc]
                    < timeprops_ptr->cur_time))
               || ((timeprops_ptr->iter == 0) && (start_time[isrc] == 0.0)))
                return (1);

        return (0);
    }

    //! this function returns the maximum of all currently active extrusion fluxes, Keith wrote this function
    double MaxInfluxNow(MatProps *matprops_ptr, TimeProps *timeprops_ptr)
    {
        double tempinflux, maxinflux = 0.0;

        for(int isrc = 0; isrc < no_of_sources; isrc++)
            if(((start_time[isrc] <= timeprops_ptr->cur_time) && (timeprops_ptr->cur_time <= end_time[isrc])) || ((timeprops_ptr
                    ->iter
                                                                                                                   == 0)
                                                                                                                  && (start_time[isrc] == 0)))
            {
                tempinflux = sqrt(
                        influx[isrc] * influx[isrc] + (xVel[isrc] * xVel[isrc] + yVel[isrc] * yVel[isrc])
                                / ((matprops_ptr->scale.epsilon) * (matprops_ptr->scale.epsilon)));
                if(tempinflux > maxinflux)
                    maxinflux = tempinflux;
            }

        return (maxinflux);
    }
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="FluxProps") const;
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="FluxProps");

};

inline void MatProps::set_scale(const PileProps *pileprops_ptr, const FluxProps *fluxprops_ptr)
{
    int isrc;
    double doubleswap;

    //this is used in ../geoflow/stats.C ... might want to set
    //MAX_NEGLIGIBLE_HEIGHT to zero now that we have "good" thin
    //layer control, need to reevaluate this, we should also
    //reevaluate after we implement a "good" local stopping criteria
    scale.max_negligible_height = scale.height / 10000.0;

    if(scale.height == 0.0 || scale.auto_calc_height_scale)
    {

        double totalvolume = 0.0;

        if(pileprops_ptr != NULL)
            for(isrc = 0; isrc < pileprops_ptr->numpiles; isrc++)
                totalvolume += pileprops_ptr->get_volume(isrc);
        if(fluxprops_ptr != NULL)
            for(isrc = 0; isrc < fluxprops_ptr->no_of_sources; isrc++)
                totalvolume += fluxprops_ptr->get_volume(isrc);

        doubleswap = pow(totalvolume, 1.0 / 3.0);

        if((scale.gravity != 1.0) || (scale.length != 1.0))
            scale.height = doubleswap;
        else
            scale.height = 1.0;

        if(scale.auto_calc_height_scale==true && totalvolume == 0.0)
        {
            printf("ERROR: totalvolume of all moving materials is zero, do you have piles or flux sources? Can not autocalculate scale.height!");
            assert(0);
        }

        scale.auto_calc_height_scale=true;

        scale.max_negligible_height = doubleswap / scale.height / 10000.0;
    }

    scale.epsilon = scale.height / scale.length;
}
inline void MatProps::calc_Vslump(const PileProps *pileprops_ptr, const FluxProps *fluxprops_ptr)
{
    /*************************************************************************
     * Vslump doesn't mean it is related to slumping of pile
     * It is simply the maximum hypothetical velocity from
     * free fall of the pile.
     ************************************************************************/
    /*int i;
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
        double xcen = scale.length * pileprops_ptr->xCen[i];
        double ycen = scale.length * pileprops_ptr->yCen[i];
        double ztemp = 0;
        ierr = Get_elevation(res, xcen, ycen, ztemp);
        if(ierr != 0)
            ztemp = 0;

        if((ztemp + pileprops_ptr->pileheight[i]) > (zcen + pileprops_ptr->pileheight[i]))
        {
            zcen = ztemp;
            j = i;
        }
    }

    // calculate Vslump
    if(pileprops_ptr->numpiles>0)
    {
        double hscale = (zcen - zmin) + pileprops_ptr->pileheight[j];
        Vslump = sqrt(gravity * hscale);
    }
    else
    {
        Vslump = 1.0;
    }*/
    Vslump = 1.0;
}
inline void MatPropsTwoPhases::set_scale(const PileProps *pileprops1_ptr,
                                         const FluxProps *fluxprops_ptr)
{
    MatProps::set_scale(pileprops1_ptr, fluxprops_ptr);

    PilePropsTwoPhases *pileprops_ptr = (PilePropsTwoPhases*) pileprops1_ptr;

    int isrc;
    double maxphi = 0.;
    double minphi = HUGE_VAL;
    for(isrc = 0; isrc < pileprops_ptr->numpiles; isrc++)
    {
        // search for min-max phi
        if(pileprops_ptr->vol_fract[isrc] > maxphi)
            maxphi = pileprops_ptr->vol_fract[isrc];
        if(pileprops_ptr->vol_fract[isrc] < minphi)
            minphi = pileprops_ptr->vol_fract[isrc];
    }
    // cut-off extremes
    //assert(minphi <= maxphi);
    if(minphi < 0.2)
        flow_type = FLUID_FLOW;
    else if(maxphi > 0.9)
        flow_type = DRY_FLOW;
    else
        flow_type = TWOPHASE;
}
inline void MatPropsTwoPhases::calc_Vslump(const PileProps *pileprops_ptr, const FluxProps *fluxprops_ptr)
{
    /*************************************************************************/
    /* the non-dimensional velocity stopping criteria is an idea that
     didn't work for anything other than a slumping pile on a horizontal
     surface, it's only still here because I didn't want to bother
     with removing it.  --Keith Dalbey 2005.10.28

     kappa is a to be determined constant, calculation stops when
     v*=v_ave/v_slump<kappa (or perhaps v* < kappa/tan(intfrict)) */
    double kappa = 1.0;   //should eventually move to a header file
    double gravity = 9.8; //[m/s^2]
    Vslump = 1.0; //kappa*sqrt(gravity*max_init_height);
}










#endif
