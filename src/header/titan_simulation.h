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

#ifndef TITAN_SIMULATION_H
#define TITAN2D_SIMULATION_H

#include <string>
#include <vector>
#include "../gisapi/GisApi.h"

class MaterialMap {
public:
    MaterialMap();
    ~MaterialMap();
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

class cxxTitanSimulation {
public:
	cxxTitanSimulation();
	~cxxTitanSimulation();
	void run();
	void input_summary();

	int myid;
	int numprocs;

	static const int GDAL=GDAL;
	static const int GIS_GRASS=GIS_GRASS;

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

	bool region_limits_set;

    double min_location_x;
    double max_location_x;
    double min_location_y;
    double max_location_y;

    MaterialMap materialMap;
};
#endif
