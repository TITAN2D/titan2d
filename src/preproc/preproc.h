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

// node.h for the definition of node class
#ifndef TITAN_PREPROC_H
#define TITAN_PREPROC_H

#include <string>


class NodePreproc;
class ElementPreproc;
class BoundaryPreproc;

#include "../header/titan_simulation.h"

/*
 * "You entered %d arguments, preprocess requires 4, 7, 8 or 11 arguments."
    "In order, they are:\n"
    "1) number of processes,\n"
    "2) the number of cells across smallest pile dimension,\n"
    "3) GIS data format\n"
    "If using GIS GRASS data format:\n"
    "---------------------------------\n"
    "4) the full path of the GIS database,\n"
    "5) the location,\n"
    "6) the mapset,\n"
    "7) the raster map name,\n\n"
    "If using Gdal data format:\n"
    "---------------------------------\n"
    "4) raster map subdirectory,\n\n"
    "Optionally\n"
    "5|8) the requested minimum x,\n"
    "6|9) the requested minimum y,\n"
    "7|10) the requested maximum x\n"
    "8|11) the requested maximum y.\n"
    "Please enter the correct number of arguments.\n"
 */
class TitanPreproc {
public:
    TitanPreproc(cxxTitanSimulation *tSim);
    ~TitanPreproc();

    //!>number of processes
    int NumProc;
    //!>the number of cells across smallest pile dimension

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

    MaterialMap material_map;




    //If using Gdal data format:
    //raster map subdirectory, mapname
    //Optionally,requested xmin, ymin, xmax, ymax
    bool region_limits_set;
    double min_location_x;
    double max_location_x;
    double min_location_y;
    double max_location_y;

    bool validate();
    void run();

    void createfunky(double limits[4], int *node_count, NodePreproc **node,
        int *element_count, ElementPreproc **element, int *force_count, int *constraint_count,
        BoundaryPreproc **boundary, int *material_count, char ***materialnames, double **lambda, double **mu);
    void Read_material_data(int *material_count, char ***materialnames, double **lambda, double **mu);


};
#endif
