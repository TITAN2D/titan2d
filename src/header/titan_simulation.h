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

/**
 * Information for Pile
 * Thickness of Initial Volume, h(x,y)
 * P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)
 */
class cxxTitanPile{
public:
    cxxTitanPile();
    ~cxxTitanPile();

    cxxTitanPile& operator=(const cxxTitanPile& other);

    //! Maximum Initial Thickness, P (m)
    double height;
    //! Center of Initial Volume, xc, yc (UTM E, UTM N)
    double xcenter;
    double ycenter;
    //! Major and Minor Extent, majorR, minorR (m, m)
    double majradius;
    double minradius;
    //! Orientation (angle [degrees] from X axis to major axis)
    double orientation;
    //! Initial speed [m/s]
    double Vmagnitude;
    //! Initial direction ([degrees] from X axis)
    double Vdirection;

    //! calculate volume m^3
    double get_volume();

    void print0();
};

/**
 * cxxTitanFluxSource
 */
class cxxTitanFluxSource{
public:
    cxxTitanFluxSource();
    ~cxxTitanFluxSource();

    cxxTitanFluxSource& operator=(const cxxTitanFluxSource& other);

    //! Extrusion flux rate [m/s]
    double influx;

    //! Active Time [s], start, end
    double start_time;
    double end_time;
    //! Center of Initial Volume, xc, yc (UTM E, UTM N)
    double xcenter;
    double ycenter;
    //! Major and Minor Extent, majorR, minorR (m, m)
    double majradius;
    double minradius;
    //! Orientation (angle [degrees] from X axis to major axis)
    double orientation;
    //! Initial speed [m/s]
    double Vmagnitude;
    //! Initial direction ([degrees] from X axis)
    double Vdirection;

    double get_effective_height();
    void print0();
};

/**
 * cxxTitanDischargePlane
 */
class cxxTitanDischargePlane{
public:
    cxxTitanDischargePlane();
    cxxTitanDischargePlane::cxxTitanDischargePlane(const double m_x_a, const double m_y_a, const double m_x_b, const double m_y_b);
    ~cxxTitanDischargePlane();

    cxxTitanDischargePlane& operator=(const cxxTitanDischargePlane& other);

    double x_a,y_a,x_b,y_b;
    void print0();
};

/**
 * MaterialMap
 */
class MaterialMap {
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
    //!>GIS Vector
    std::string topovector;

    //!>Piles
    std::vector<cxxTitanPile> piles;

    //!>Flux sources
    std::vector<cxxTitanFluxSource> flux_sources;

    //!>Discharge planes
    std::vector<cxxTitanDischargePlane> discharge_planes;

	bool region_limits_set;

    double min_location_x;
    double max_location_x;
    double min_location_y;
    double max_location_y;

    MaterialMap material_map;
};
#endif
