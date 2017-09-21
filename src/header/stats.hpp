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


#ifndef SRC_HEADER_STATS_HPP_
#define SRC_HEADER_STATS_HPP_

#include "hpfem.h"
#include "hd5calls.h"
//! the StatProps structure holds statistics about the flow
class StatProps: public EleNodeRef
{
public:
    //note all means are mass/volume averages

    //! job number for monte carlo or lhs simulations
    int runid;

    //! x coordinate of pile centroid
    double xcen;

    //! y coordinate of pile centroid
    double ycen;

    //! variance of location of pile material in the x direction
    double xvar;

    //! variance of location of pile material in the y direction
    double yvar;

    //! mean distance from the point (0,0), this was the mean pile starting point for the ordered reduced number of runs monte carlo method validation
    double rmean;

    //! area covered by pile of thickness greater than cutoffheight
    double area;

    //! mean speed
    double vmean;

    //! mean x velocity
    double vxmean;

    //! mean y velocity
    double vymean;

    //! mean slope in the direction of velocity a negative number   indicates the flow is heading uphill
    double slopemean;

    //! nondimensionalized mean speed
    double vstar;

    //! volume of ALL material, NOT used for other statistics
    double realvolume;

    //! STAT_VOL_FRAC*realvolume, all statistics computed in ../geoflow/stats.C are based on statvolume not realvolume
    double statvolume;

    //! volume of material that have flown off the map
    double outflowvol;

    //! volume of material that has been eroded
    double erodedvol;

    //! volume of material that is currently deposited
    double depositedvol;

    //! pile height of contour line that encloses statvolume
    double cutoffheight;

    //! an estimate of the radius of the pile
    double piler;

    //! current spatial maximum of pile height
    double hmax;

    //! current spatial maximum of speed
    double vmax;

    //! the integrated magnitude of acceleration due to internal friction force (based on realvolume not statvolume)
    double forceint;

    //! the same thing but for bed friction
    double forcebed;

    //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
    double heightifreach;

    //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
    double xyifreach[2];

    //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
    double timereached;

    //! xyminmax holds the minimum and maximum x and y coordinates where the pile height is greater than hxyminmax
    double xyminmax[4];

    //! xyminmax holds the minimum and maximum x and y coordinates where the pile height is greater than hxyminmax
    double hxyminmax;

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

    double Vol_;
    double Area_;
    double Velmean_;

    //! the latin hypercube sampling specific stats
    LHS_Props lhs;

    std::string output_prefix;

    //! the constructor initializes a few statistics
    StatProps(ElementsHashTable *_ElemTable, NodeHashTable* _NodeTable);
    //! the constructor from hd5
    StatProps(ElementsHashTable *_ElemTable, NodeHashTable* _NodeTable, const H5::CommonFG *parent, const  string group_name="StatProps");
    ~StatProps();
    void set(const double edge_height, const double test_height, const double test_location_x, const double test_location_y);
    void scale(const MatProps* matprops_ptr);

    /** This function calculates the vast majority of statistics used for output, including most of what appears
     *  in output_summary.######, the friction body forces however are not calculated in here, Keith wrote this
     *  to replace calc_volume()
     */
    void calc_stats(int myid, MatProps* matprops, TimeProps* timeprops,
            DischargePlanes* discharge, double d_time);
    //! Dump object content to hdf5 file
    void h5write(H5::CommonFG *parent, string group_name="StatProps") const;
    //! Load object content from hdf5 file
    void h5read(const H5::CommonFG *parent, const  string group_name="StatProps");
};



#endif /* SRC_HEADER_STATS_HPP_ */
