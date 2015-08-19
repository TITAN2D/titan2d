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

#include "../header/properties.h"
#include "../header/hpfem.h"

MapNames::MapNames()
{
    gis_main = "";
    gis_sub = "";
    gis_mapset = "";
    gis_map = "";
    gis_vector = "";
    gis_format = MapNames::GIS_GRASS;
    extramaps = 0;

    region_limits_set = false;

    min_location_x = 0.0;
    min_location_y = 0.0;
    max_location_x = 0.0;
    max_location_y = 0.0;
}
MapNames::~MapNames()
{
}
void MapNames::set(const int format, const std::string gis_main_in, const std::string gis_sub_in,
                   const std::string gis_mapset_in, const std::string gis_map_in, const std::string gis_vector_in,
                   const int extramaps_in)
{
    gis_main = gis_main_in;
    gis_sub = gis_sub_in;
    gis_mapset = gis_mapset_in;
    gis_map = gis_map_in;
    gis_vector = gis_vector_in;
    gis_format = format;
    extramaps = extramaps_in;
    return;
}
void MapNames::set_region_limits(double m_min_location_x, double m_max_location_x, double m_min_location_y,
                                 double m_max_location_y)
{
    region_limits_set = true;
    min_location_x = m_min_location_x;
    min_location_y = m_min_location_y;
    max_location_x = m_max_location_x;
    max_location_y = m_max_location_y;
}
void MapNames::print0()
{
    printf("GIS:\n");
    printf("\tgis_format: %d\n", gis_format);

    printf("\tgis_main: %s\n", gis_main.c_str());
    printf("\tgis_sub: %s\n", gis_sub.c_str());
    printf("\tgis_mapset: %s\n", gis_mapset.c_str());
    printf("\tgis_map: %s\n", gis_map.c_str());
    printf("\tgis_vector: %s\n", gis_vector.c_str());
    printf("\textramaps: %d\n", extramaps);
    printf("\tregion_limits_set %d\n", (int) region_limits_set);
    if(region_limits_set)
    {
        printf("\tregion limits %g %g %g %g\n", min_location_x, min_location_y, max_location_x, max_location_y);
    }
    return;
}
StatProps::StatProps()
{
    timereached = -1.0;
    xcen = ycen = xvar = yvar = rmean = area = vmean = vxmean = vymean = slopemean = vstar = 0.0;
    realvolume = statvolume = outflowvol = erodedvol = depositedvol = cutoffheight = 0.0;
    piler = hmax = vmax = forceint = forcebed = 0.0;
    heightifreach = xyifreach[0] = xyifreach[1] = timereached = 0.0;
    xyminmax[0] = xyminmax[1] = xyminmax[2] = xyminmax[3] = hxyminmax = 0.0;
    lhs.refnum = lhs.runid = -1;
    runid = -1;
}
StatProps::~StatProps()
{
}
void StatProps::set(const double edge_height, const double test_height, const double test_location_x, const double test_location_y)
{
    hxyminmax = edge_height;

    heightifreach=test_height;
    if(heightifreach < 0.0)
    {
        xyifreach[0]=test_location_x;
        xyifreach[1]=test_location_y;
    }
    else
    {
        heightifreach = xyifreach[0] = xyifreach[1] = HUGE_VAL;
    }

    //to get rid on uninitiallized memory error in saverun() (restart.C)
    forceint = forcebed = 0.0;
}
void StatProps::scale(const MatProps* matprops_ptr)
{
    if(hxyminmax < 0.0)
    {
        hxyminmax = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT * 10.0;

    }
    hxyminmax /= matprops_ptr->HEIGHT_SCALE;

    if(heightifreach > -1.9)
    {

        //default test height is 10 time the maximum negligible height
        if(heightifreach > -1.1 && heightifreach < -0.9)
            heightifreach = matprops_ptr->MAX_NEGLIGIBLE_HEIGHT * 10.0;

        heightifreach /= matprops_ptr->HEIGHT_SCALE;

        xyifreach[0] /= matprops_ptr->LENGTH_SCALE;
        xyifreach[1] /= matprops_ptr->LENGTH_SCALE;
    }
    else
    {
        heightifreach = xyifreach[0] = xyifreach[1] =
        HUGE_VAL;
    }

    //to get rid on uninitiallized memory error in saverun() (restart.C)
    forceint = forcebed = 0.0;
}
PileProps::PileProps()
{
    numpiles = 0;
    length_scale = 1.0;
    height_scale = 1.0;
    velocity_scale = 1.0;
}
PileProps::~PileProps()
{
}
//! function allocates space for the pile data
void PileProps::allocpiles(int numpiles_in)
{
    numpiles = numpiles_in;
    pileheight.resize(numpiles);

    xCen.resize(numpiles);
    yCen.resize(numpiles);
    majorrad.resize(numpiles);
    minorrad.resize(numpiles);
    cosrot.resize(numpiles);
    sinrot.resize(numpiles);
    initialVx.resize(numpiles);
    initialVy.resize(numpiles);
    pile_type.resize(numpiles);
}

void PileProps::addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                        double orientation, double Vmagnitude, double Vdirection, PileType m_pile_type)
{
    numpiles++;
    pileheight.push_back(hight);
    xCen.push_back(xcenter);
    yCen.push_back(ycenter);
    majorrad.push_back(majradius);
    minorrad.push_back(minradius);
    cosrot.push_back(cos(orientation * PI / 180.0));
    sinrot.push_back(sin(orientation * PI / 180.0));
    initialVx.push_back(Vmagnitude * cos(Vdirection * PI / 180.0));
    initialVy.push_back(Vmagnitude * sin(Vdirection * PI / 180.0));
    pile_type.push_back(m_pile_type);
}

void PileProps::scale(double m_length_scale, double m_height_scale, double m_gravity_scale)
{
    length_scale = m_length_scale;
    height_scale = m_height_scale;
    //non-dimensionalize the inputs
    velocity_scale = sqrt(length_scale * m_gravity_scale);
    int isrc;
    for(isrc = 0; isrc < numpiles; isrc++)
    {
        pileheight[isrc] /= height_scale;
        xCen[isrc] /= length_scale;
        yCen[isrc] /= length_scale;
        majorrad[isrc] /= length_scale;
        minorrad[isrc] /= length_scale;
        initialVx[isrc] /= velocity_scale;
        initialVy[isrc] /= velocity_scale;
    }
}
double PileProps::get_smallest_pile_radius()
{
    double smallestpileradius = HUGE_VAL;
    int isrc;
    for(isrc = 0; isrc < numpiles; isrc++)
    {
        if(smallestpileradius > majorrad[isrc])
            smallestpileradius = majorrad[isrc];

        if(smallestpileradius > minorrad[isrc])
            smallestpileradius = minorrad[isrc];
    }
    return smallestpileradius;
}
void PileProps::print_pile(int i)
{
    printf("\tPile %d\n", i);
    printf("\t\tMaximum Initial Thickness, P (m):%f\n", pileheight[i] * height_scale);
    printf("\t\tCenter of Initial Volume, xc, yc (UTM E, UTM N): %f %f\n", xCen[i] * length_scale,
           yCen[i] * length_scale);
    printf("\t\tMajor and Minor Extent, majorR, minorR (m, m): %f %f\n", majorrad[i] * length_scale,
           minorrad[i] * length_scale);
    double orientation = atan2(sinrot[i], cosrot[i]) * 180.0 / PI;
    printf("\t\tOrientation (angle [degrees] from X axis to major axis): %f\n", orientation);
    double Vmagnitude = sqrt(initialVx[i] * initialVx[i] + initialVy[i] * initialVy[i]);
    double Vdirection = atan2(initialVy[i], initialVx[i]) * 180.0 / PI;
    printf("\t\tInitial speed [m/s]: %f\n", Vmagnitude * velocity_scale);
    printf("\t\tInitial direction ([degrees] from X axis): %f\n", Vdirection);
    printf("\t\tPile type: %d\n", pile_type[i]);
    printf("\t\tPile volume [m^3]: %f\n", get_volume(i) * height_scale * length_scale * length_scale);
}
void PileProps::print0()
{
    int i;
    if(numpiles > 0)
    {
        printf("Piles:    (Number of piles: %d)\n", numpiles);
        for(i = 0; i < numpiles; i++)
            print_pile(i);
    }
    else
    {
        printf("Piles:    there is no piles\n");
    }
}
void PileProps::set_element_height_to_elliptical_pile_height(HashTable* HT_Node_Ptr, Element *m_EmTemp, MatProps* matprops)
{
    double pileheight;
    double xmom, ymom;
    pileheight=get_elliptical_pile_height(HT_Node_Ptr, m_EmTemp, matprops, &xmom,&ymom);

    ElementSinglePhase* EmTemp=(ElementSinglePhase*)m_EmTemp;
    EmTemp->put_height_mom(pileheight, xmom, ymom);
}
double PileProps::get_elliptical_pile_height(HashTable* HT_Node_Ptr, Element *EmTemp, MatProps* matprops, double* m_xmom,
                                         double* m_ymom)
{
    SFC_Key nodes[9];

    //get corner and edge nodes
    SFC_Key *node_key = EmTemp->getNode();
    for(int inode = 0; inode < 8; inode++)
        nodes[inode] = node_key[inode];

    //get center node
    nodes[8] = *(EmTemp->pass_key());

    double node_pile_height[9];
    double sum_node_pile_height[9];
    double sum_node_xmom[9];
    double sum_node_ymom[9];
    double height;

    for(int inode = 0; inode < 9; inode++)
    {

        //get pile height at each node...
        Node* ndtemp = (Node*) HT_Node_Ptr->lookup(nodes[inode]);
        double* ndcoord = ndtemp->get_coord();

        // for multiple piles which may overlap, the highest value is used..
        node_pile_height[inode] = 0.0;
        sum_node_pile_height[inode] = 0.0;
        sum_node_xmom[inode] = 0.0;
        sum_node_ymom[inode] = 0.0;

        //check each pile to see which has max height at this node
        for(int ipile = 0; ipile < numpiles; ipile++)
        {
            //get position relative to pile center
            double major = ndcoord[0] - xCen[ipile];
            double minor = ndcoord[1] - yCen[ipile];

            /* "undo" elliptical pile rotation ... from (x,y)->(major,minor)
             also make  nondimensional (by dividing by major and minor radius) */
            double doubleswap = (major * cosrot[ipile] + minor * sinrot[ipile])
                    / majorrad[ipile];

            minor = (-major * sinrot[ipile] + minor * cosrot[ipile]) / minorrad[ipile];
            major = doubleswap;

            /* calculate pile height based on non dimensional position relative to
             center of pile */

            if(pile_type[ipile] == PileProps::PARABALOID)
            {
                height = pileheight[ipile] * (1. - major * major - minor * minor);
            }
            else if(pile_type[ipile] == PileProps::CYLINDER)
            {
                if(major * major + minor * minor < 1.0)
                    height = pileheight[ipile];
                else
                    height = 0.0;
            }
            else
            {
                printf("Unknown type of pile\n");
                assert(0);
            }
            height = (height >= 0.0) ? height : 0.0;

            sum_node_pile_height[inode] += height;
            sum_node_xmom[inode] += height * (initialVx[ipile]);
            sum_node_ymom[inode] += height * (initialVy[ipile]);

            if(node_pile_height[inode] < height)
                node_pile_height[inode] = height;
        }
        if(sum_node_pile_height[inode] <= GEOFLOW_TINY)
            sum_node_xmom[inode] = sum_node_ymom[inode] = 0.0;
        else
        {
            sum_node_xmom[inode] *= height / sum_node_pile_height[inode];
            sum_node_ymom[inode] *= height / sum_node_pile_height[inode];
            //these are now the averaged momentums at each node
        }
    }

    /* The pile_height value assigned is an "area" weighted average over the
     element's 9 nodes.  The element is divided into 4 squares, and each
     corner of each of the 4 squares count once.  Because the center node
     is repeated 4 times it's weight is 4 times as much as the element's
     corner nodes which are not repeated; each edge node is repeated
     twice */
    double pileheight = ( //corner nodes
    node_pile_height[0] + node_pile_height[1] + node_pile_height[2] + node_pile_height[3] +
    //edge nodes
    2.0 * (node_pile_height[4] + node_pile_height[5] + node_pile_height[6] + node_pile_height[7]) +
    //center node
    4.0 * node_pile_height[8])
                        / 16.0;

    double xmom = ( //corner nodes
    sum_node_xmom[0] + sum_node_xmom[1] + sum_node_xmom[2] + sum_node_xmom[3] +
    //edge nodes
    2.0 * (sum_node_xmom[4] + sum_node_xmom[5] + sum_node_xmom[6] + sum_node_xmom[7]) +
    //center node
    4.0 * sum_node_xmom[8])
                  / 16.0;

    double ymom = ( //corner nodes
    sum_node_ymom[0] + sum_node_ymom[1] + sum_node_ymom[2] + sum_node_ymom[3] +
    //edge nodes
    2.0 * (sum_node_ymom[4] + sum_node_ymom[5] + sum_node_ymom[6] + sum_node_ymom[7]) +
    //center node
    4.0 * sum_node_ymom[8])
                  / 16.0;

    if(m_xmom!=NULL)
        *m_xmom=xmom;
    if(m_ymom!=NULL)
            *m_ymom=ymom;

    return pileheight;
}

PilePropsTwoPhases::PilePropsTwoPhases() :
        PileProps()
{
}
PilePropsTwoPhases::~PilePropsTwoPhases()
{
}

void PilePropsTwoPhases::allocpiles(int numpiles_in)
{
    PileProps::allocpiles(numpiles_in);
    vol_fract.resize(numpiles);
}
void PilePropsTwoPhases::addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                                 double orientation, double Vmagnitude, double Vdirection, PileProps::PileType m_pile_type)
{
    addPile(hight, xcenter, ycenter, majradius, minradius, orientation, Vmagnitude, Vdirection, m_pile_type, 1.0);
}
void PilePropsTwoPhases::addPile(double hight, double xcenter, double ycenter, double majradius, double minradius,
                                 double orientation, double Vmagnitude, double Vdirection, PileProps::PileType m_pile_type, double volfract)
{
    PileProps::addPile(hight, xcenter, ycenter, majradius, minradius, orientation, Vmagnitude, Vdirection, m_pile_type);
    vol_fract.push_back(volfract);
}
void PilePropsTwoPhases::print_pile(int i)
{
    PileProps::print_pile(i);
    printf("\t\tInitial solid-volume fraction,(0:1.): %f\n", vol_fract[i]);
}
void PilePropsTwoPhases::set_element_height_to_elliptical_pile_height(HashTable* HT_Node_Ptr, Element *m_EmTemp, MatProps* matprops)
{
    double pileheight;
    double xmom, ymom;
    pileheight=get_elliptical_pile_height(HT_Node_Ptr, m_EmTemp, matprops, &xmom,&ymom);

    ElementTwoPhases* EmTemp=(ElementTwoPhases*)m_EmTemp;
    int ipile;
    double vfract = 0.;
    for(ipile = 0; ipile < numpiles; ipile++)
    {
        if(vol_fract[ipile] > vfract)
        vfract = vol_fract[ipile];
    }
    EmTemp->put_height_mom(pileheight, vfract, xmom, ymom);
}

