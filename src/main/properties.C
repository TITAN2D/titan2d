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

#include <assert.h>
#include "../header/properties.h"
#include "../header/hpfem.h"


#include <sstream>

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

void MapNames::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));

    TiH5_writeStringAttribute(group, gis_main);
    TiH5_writeStringAttribute(group, gis_sub);
    TiH5_writeStringAttribute(group, gis_mapset);
    TiH5_writeStringAttribute(group, gis_map);
    TiH5_writeStringAttribute(group, gis_vector);
    TiH5_writeIntAttribute(group, gis_format);
    TiH5_writeIntAttribute(group, extramaps);
    TiH5_writeBoolAttribute(group, region_limits_set);
    TiH5_writeDoubleAttribute(group, min_location_x);
    TiH5_writeDoubleAttribute(group, min_location_y);
    TiH5_writeDoubleAttribute(group, max_location_x);
    TiH5_writeDoubleAttribute(group, max_location_y);
}
void MapNames::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readStringAttribute(group, gis_main);
    TiH5_readStringAttribute(group, gis_sub);
    TiH5_readStringAttribute(group, gis_mapset);
    TiH5_readStringAttribute(group, gis_map);
    TiH5_readStringAttribute(group, gis_vector);
    TiH5_readIntAttribute(group, gis_format);
    TiH5_readIntAttribute(group, extramaps);
    TiH5_readBoolAttribute(group, region_limits_set);
    TiH5_readDoubleAttribute(group, min_location_x);
    TiH5_readDoubleAttribute(group, min_location_y);
    TiH5_readDoubleAttribute(group, max_location_x);
    TiH5_readDoubleAttribute(group, max_location_y);
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
void PileProps::set_element_height_to_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *m_EmTemp, MatProps* matprops)
{
    double pileheight;
    double xmom, ymom;
    pileheight=get_elliptical_pile_height(HT_Node_Ptr, m_EmTemp, matprops, &xmom,&ymom);

    m_EmTemp->put_height_mom(pileheight, xmom, ymom);
}
double PileProps::get_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *EmTemp, MatProps* matprops, double* m_xmom,
                                         double* m_ymom)
{
    SFC_Key nodes[9];

    //get corner and edge nodes
    for(int inode = 0; inode < 8; inode++)
        nodes[inode] = EmTemp->node_key(inode);

    //get center node
    nodes[8] = EmTemp->key();

    double node_pile_height[9];
    double sum_node_pile_height[9];
    double sum_node_xmom[9];
    double sum_node_ymom[9];
    double height;

    for(int inode = 0; inode < 9; inode++)
    {

        //get pile height at each node...
        Node* ndtemp = (Node*) HT_Node_Ptr->lookup(nodes[inode]);

        // for multiple piles which may overlap, the highest value is used..
        node_pile_height[inode] = 0.0;
        sum_node_pile_height[inode] = 0.0;
        sum_node_xmom[inode] = 0.0;
        sum_node_ymom[inode] = 0.0;

        //check each pile to see which has max height at this node
        for(int ipile = 0; ipile < numpiles; ipile++)
        {
            //get position relative to pile center
            double major = ndtemp->coord(0) - xCen[ipile];
            double minor = ndtemp->coord(1) - yCen[ipile];

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
void PileProps::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));
    TiH5_writeStringAttribute__(group,"PileProps","Type");
    TiH5_writeIntAttribute(group, numpiles);
    TiH5_writeDoubleAttribute(group, height_scale);
    TiH5_writeDoubleAttribute(group, length_scale);
    TiH5_writeDoubleAttribute(group, velocity_scale);

    TiH5_writeDoubleVectorAttribute(group, pileheight, numpiles);
    TiH5_writeDoubleVectorAttribute(group, xCen, numpiles);
    TiH5_writeDoubleVectorAttribute(group, yCen, numpiles);
    TiH5_writeDoubleVectorAttribute(group, majorrad, numpiles);
    TiH5_writeDoubleVectorAttribute(group, minorrad, numpiles);
    TiH5_writeDoubleVectorAttribute(group, cosrot, numpiles);
    TiH5_writeDoubleVectorAttribute(group, sinrot, numpiles);
    TiH5_writeDoubleVectorAttribute(group, initialVx, numpiles);
    TiH5_writeDoubleVectorAttribute(group, initialVy, numpiles);

    TiH5_writeDataTypeVectorAttribute(group,pile_type,datatypePileType,numpiles);
}
void PileProps::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readIntAttribute(group, numpiles);
    TiH5_readDoubleAttribute(group, height_scale);
    TiH5_readDoubleAttribute(group, length_scale);
    TiH5_readDoubleAttribute(group, velocity_scale);

    TiH5_readDoubleVectorAttribute(group, pileheight, numpiles);
    TiH5_readDoubleVectorAttribute(group, xCen, numpiles);
    TiH5_readDoubleVectorAttribute(group, yCen, numpiles);
    TiH5_readDoubleVectorAttribute(group, majorrad, numpiles);
    TiH5_readDoubleVectorAttribute(group, minorrad, numpiles);
    TiH5_readDoubleVectorAttribute(group, cosrot, numpiles);
    TiH5_readDoubleVectorAttribute(group, sinrot, numpiles);
    TiH5_readDoubleVectorAttribute(group, initialVx, numpiles);
    TiH5_readDoubleVectorAttribute(group, initialVy, numpiles);

    TiH5_readDataTypeVectorAttribute(group,pile_type,datatypePileType,numpiles);
}
PileProps* PileProps::createPileProps(const H5::CommonFG *parent, const  string group_name)
{
    PileProps *pileProps=nullptr;
    string PilePropsType;
    H5::Group group(parent->openGroup(group_name));
    TiH5_readStringAttribute__(group,PilePropsType,"Type");

    if (PilePropsType == "PileProps")
    {
        pileProps = new PileProps();
    }
    else if (PilePropsType == "PilePropsTwoPhases")
    {
        pileProps = new PilePropsTwoPhases();
    }
    else
    {
        cout << "ERROR: Unknown type of PileProps:" << PilePropsType << "\n";
    }
    assert(pileProps != nullptr);
    pileProps->h5read(parent,group_name);

    return pileProps;
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
void PilePropsTwoPhases::set_element_height_to_elliptical_pile_height(NodeHashTable* HT_Node_Ptr, Element *m_EmTemp, MatProps* matprops)
{
    double pileheight;
    double xmom, ymom;
    pileheight=get_elliptical_pile_height(HT_Node_Ptr, m_EmTemp, matprops, &xmom,&ymom);

    int ipile;
    double vfract = 0.;
    for(ipile = 0; ipile < numpiles; ipile++)
    {
        if(vol_fract[ipile] > vfract)
        vfract = vol_fract[ipile];
    }
    m_EmTemp->put_height_mom(pileheight, vfract, xmom, ymom);
}
void PilePropsTwoPhases::h5write(H5::CommonFG *parent, string group_name) const
{
    PileProps::h5write(parent, group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"PilePropsTwoPhases","Type");
    TiH5_writeDoubleVectorAttribute(group, vol_fract, numpiles);
}
void PilePropsTwoPhases::h5read(const H5::CommonFG *parent, const  string group_name)
{
    PileProps::h5read(parent, group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_readDoubleVectorAttribute(group, vol_fract, numpiles);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FluxProps::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));
    TiH5_writeIntAttribute(group, no_of_sources);

    TiH5_writeDoubleVectorAttribute(group, influx, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, start_time, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, end_time, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, xCen, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, yCen, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, majorrad, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, minorrad, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, cosrot, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, sinrot, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, xVel, no_of_sources);
    TiH5_writeDoubleVectorAttribute(group, yVel, no_of_sources);

    TiH5_writeDoubleAttribute(group, height_scale);
    TiH5_writeDoubleAttribute(group,  length_scale);
    TiH5_writeDoubleAttribute(group,  velocity_scale);
    TiH5_writeDoubleAttribute(group,  time_scale);

}
void FluxProps::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readIntAttribute(group, no_of_sources);

    TiH5_readDoubleVectorAttribute(group, influx, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, start_time, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, end_time, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, xCen, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, yCen, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, majorrad, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, minorrad, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, cosrot, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, sinrot, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, xVel, no_of_sources);
    TiH5_readDoubleVectorAttribute(group, yVel, no_of_sources);

    TiH5_readDoubleAttribute(group, height_scale);
    TiH5_readDoubleAttribute(group,  length_scale);
    TiH5_readDoubleAttribute(group,  velocity_scale);
    TiH5_readDoubleAttribute(group,  time_scale);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DischargePlanes::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));
    TiH5_writeIntAttribute(group, num_planes);

    TiH5_writeVectorArrayAttribute(group, planes);
}
void DischargePlanes::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readIntAttribute(group, num_planes);

    TiH5_readVectorArrayAttribute(group, planes);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MatProps* MatProps::createMatProps(const H5::CommonFG *parent,TiScale &_scale, const  string group_name)
{
    MatProps *matProps=nullptr;
    string matPropsType;
    H5::Group group(parent->openGroup(group_name));
    TiH5_readStringAttribute__(group,matPropsType,"Type");

    if (matPropsType == "MatProps")
    {
        matProps = new MatProps(_scale);
    }
    else if (matPropsType == "MatPropsTwoPhases")
    {
        matProps = new MatPropsTwoPhases(_scale);
    }
    else
    {
        cout << "ERROR: Unknown type of MatProps:" << matPropsType << "\n";
    }
    assert(matProps != nullptr);
    matProps->h5read(parent,group_name);

    return matProps;
}
void MatProps::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));
    TiH5_writeStringAttribute__(group,"MatProps","Type");

    TiH5_writeIntAttribute(group, number_of_cells_across_axis);
    TiH5_writeDoubleAttribute(group, smallest_axis);
    TiH5_writeIntAttribute(group, material_count);
    TiH5_writeVectorStringAttribute(group, matnames);
    TiH5_writeDoubleVectorAttribute(group, bedfrict, matnames.size());
    TiH5_writeDoubleVectorAttribute(group, tanbedfrict, matnames.size());
    TiH5_writeDoubleAttribute(group, porosity);
    TiH5_writeDoubleAttribute(group, mu);
    TiH5_writeDoubleAttribute(group, rho);
    TiH5_writeDoubleAttribute(group, gamma);
    TiH5_writeDoubleAttribute(group, Vslump);

}
void MatProps::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readIntAttribute(group, number_of_cells_across_axis);
    TiH5_readDoubleAttribute(group, smallest_axis);
    TiH5_readIntAttribute(group, material_count);
    TiH5_readVectorStringAttribute(group, matnames);
    TiH5_readDoubleVectorAttribute(group, bedfrict, matnames.size());
    TiH5_readDoubleVectorAttribute(group, tanbedfrict, matnames.size());
    TiH5_readDoubleAttribute(group, porosity);
    TiH5_readDoubleAttribute(group, mu);
    TiH5_readDoubleAttribute(group, rho);
    TiH5_readDoubleAttribute(group, gamma);
    TiH5_readDoubleAttribute(group, Vslump);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MatPropsTwoPhases::h5write(H5::CommonFG *parent, string group_name) const
{
    MatProps::h5write(parent, group_name);
    H5::Group group(parent->openGroup(group_name));
    TiH5_writeStringAttribute__(group,"MatPropsTwoPhases","Type");

    TiH5_writeDoubleAttribute(group, den_solid);
    TiH5_writeDoubleAttribute(group, den_fluid);
    TiH5_writeDoubleAttribute(group, viscosity);
    TiH5_writeDoubleAttribute(group, v_terminal);
    TiH5_writeIntAttribute(group, flow_type);
}
void MatPropsTwoPhases::h5read(const H5::CommonFG *parent, const  string group_name)
{
    MatProps::h5read(parent, group_name);
    H5::Group group(parent->openGroup(group_name));

    TiH5_readDoubleAttribute(group, den_solid);
    TiH5_readDoubleAttribute(group, den_fluid);
    TiH5_readDoubleAttribute(group, viscosity);
    TiH5_readDoubleAttribute(group, v_terminal);
    TiH5_readIntAttribute(group, flow_type);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LHS_Props::h5write(H5::CommonFG *parent, string group_name) const
{
    H5::Group group(parent->createGroup(group_name));

    TiH5_writeIntAttribute(group, refnum);
    TiH5_writeIntAttribute(group, runid);

}
void LHS_Props::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readIntAttribute(group, refnum);
    TiH5_readIntAttribute(group, runid);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
OutLine::OutLine()
{
    ElemTable=nullptr;
    NodeTable=nullptr;
    conformation=-1;
    Nx = Ny = stride = size = 0;

    pileheight=nullptr;
    max_kinergy=nullptr;
    max_dynamic_pressure=nullptr;
    cum_kinergy=nullptr;

    pileheight_loc=nullptr;
    max_kinergy_loc=nullptr;
    max_dynamic_pressure_loc=nullptr;
    cum_kinergy_loc=nullptr;

    elementType=ElementType::SinglePhase;

    enabled=true;
    init_size=AMR;
    max_linear_size=1024;

    output_prefix="";

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    return;
}

//! this is the OutLine it deallocates the 2 dimensional array holding maximum throughout time pileheight in every cell on the map
OutLine::~OutLine()
{
    TI_FREE(pileheight);
    TI_FREE(max_kinergy);
    TI_FREE(cum_kinergy);
    if (pileheight_loc != nullptr)
    {
        for (int j = 0; j < threads_number; j++)
        {
            TI_FREE(pileheight_loc[j]);
        }
    }
    TI_FREE(pileheight_loc);

    if (max_kinergy_loc != nullptr)
    {
        for (int j = 0; j < threads_number; j++)
        {
            TI_FREE(max_kinergy_loc[j]);
        }
    }
    TI_FREE(max_kinergy_loc);
    if (max_dynamic_pressure_loc != nullptr)
    {
        for (int j = 0; j < threads_number; j++)
        {
            TI_FREE(max_dynamic_pressure_loc[j]);
        }
    }
    TI_FREE(max_dynamic_pressure_loc);

    if (cum_kinergy_loc != nullptr)
    {
        for (int j = 0; j < threads_number; j++)
        {
            TI_FREE(cum_kinergy_loc[j]);
        }
    }
    TI_FREE(cum_kinergy_loc);
    return;
}
void OutLine::setElemNodeTable(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable)
{
    ElemTable=_ElemTable;
    NodeTable=_NodeTable;
}
//! this function initializes the OutLine map/2-dimensional array
void OutLine::init(const double *dxy, int power, double *XRange, double *YRange)
{
    if(power < 0)
        power = 0;

    dx = dxy[0] / pow(2.0, power);
    dy = dxy[1] / pow(2.0, power);
    //printf("dx=%g dy=%g  XRange={%g,%g} YRange={%g,%g}\n",dx,dy,XRange[0],XRange[1],YRange[0],YRange[1]);

    xminmax[0] = XRange[0];
    xminmax[1] = XRange[1];
    yminmax[0] = YRange[0];
    yminmax[1] = YRange[1];

    Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
    Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer


    while (Nx > max_linear_size && Ny > max_linear_size)
    {
        dx *= 2.0;
        dy *= 2.0;

        Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
        Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer
    }
    //Nx * Ny should be less then 65536x65536
    assert((float)Nx * (float)Ny <65536.0*65536.0*0.5);
    stride=Nx;
    if(myid==0 && enabled)
    {
        printf("Outline init: Nx=%d Ny=%d dx=%.5g dy=%.5g\n", Nx, Ny, dx*scale->length, dy*scale->length);
        printf("Outline init: temporary arrays for %d threads\n", threads_number);
    }
    size=Ny*stride;

    if(enabled)
        alloc_arrays();
}
//! this function initializes the OutLine map/2-dimensional array
void OutLine::init_DEM_resolution(double resx, double resy, double *XRange, double *YRange)
{
    dx = resx; // OutLine.dx
    dy = resy; // OutLine.dy
    //printf("dx=%g dy=%g  XRange={%g,%g} YRange={%g,%g}\n",dx,dy,XRange[0],XRange[1],YRange[0],YRange[1]);

    xminmax[0] = XRange[0];
    xminmax[1] = XRange[1];
    yminmax[0] = YRange[0];
    yminmax[1] = YRange[1];

    Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
    Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer

    if(myid==0 && enabled)
        printf("Outline init: will try to use resolution of DEM: %.5g %.5g \n", dx*scale->length, dy*scale->length);

//    while (Nx > max_linear_size && Ny > max_linear_size)
//    {
//        if(myid==0 && enabled)
//            printf("Outline init: resolution is too high, scaling down...\n");
//        dx *= 2.0;
//        dy *= 2.0;
//
//        Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
//        Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer
//    }
    //Nx * Ny should be less then 65536x65536
    assert((float)Nx * (float)Ny <65536.0*65536.0*0.5);
    stride=Nx;
    if(myid==0 && enabled)
    {
//        printf("Outline init: Nx=%d Ny=%d dx=%.5g dy=%.5g\n", Nx, Ny, dx*scale->length, dy*scale->length);
        printf("Outline init: temporary arrays for %d threads\n", threads_number);
    }

    size=Ny*stride;

    if(enabled)
        alloc_arrays();
}
void OutLine::alloc_arrays()
{
    pileheight = TI_ALLOC(double, size);
    max_kinergy = TI_ALLOC(double, size);
    max_dynamic_pressure = TI_ALLOC(double, size);
    cum_kinergy = TI_ALLOC(double, size);

    pileheight_loc = TI_ALLOC(double*, threads_number);
    max_kinergy_loc = TI_ALLOC(double*, threads_number);
    max_dynamic_pressure_loc = TI_ALLOC(double*, threads_number);
    cum_kinergy_loc = TI_ALLOC(double*, threads_number);


    for(int j = 0; j < threads_number; j++)
    {
        pileheight_loc[j] = TI_ALLOC(double, size);
        max_kinergy_loc[j] = TI_ALLOC(double, size);
        max_dynamic_pressure_loc[j] = TI_ALLOC(double, size);
        cum_kinergy_loc[j] = TI_ALLOC(double, size);

        for(int i = 0; i < size; i++)
        {
            pileheight_loc[j][i] = 0.0;
            max_kinergy_loc[j][i] = 0.0;
            max_dynamic_pressure_loc[j][i] = 0.0;
            cum_kinergy_loc[j][i] = 0.0;
        }
    }


    for(int i = 0; i < size; i++)
    {
        pileheight[i] = 0.0;
        max_kinergy[i] = 0.0;
        max_dynamic_pressure[i] = 0.0;
        cum_kinergy[i] = 0.0;
    }
    return;
}

//! this function reinitializes the OutLine map/2-dimensional array during restart
void OutLine::init2(const double *dxy, double *XRange, double *YRange)
{
    int ix, iy;

    dx = dxy[0];
    dy = dxy[1];

    xminmax[0] = XRange[0];
    xminmax[1] = XRange[1];
    yminmax[0] = YRange[0];
    yminmax[1] = YRange[1];

    Nx = (int) ((XRange[1] - XRange[0]) / dx + 0.5); //round to nearest integer
    Ny = (int) ((YRange[1] - YRange[0]) / dy + 0.5); //round to nearest integer
    stride=Nx;

    int size=Ny*stride;
    pileheight = TI_ALLOC(double, size);
    max_kinergy = TI_ALLOC(double, size);
    max_dynamic_pressure = TI_ALLOC(double, size);
    cum_kinergy = TI_ALLOC(double, size);

    for(int i = 0; i < size; i++)
    {
        pileheight[i] = 0.0;
        max_kinergy[i] = 0.0;
        max_dynamic_pressure[i] = 0.0;
        cum_kinergy[i] = 0.0;
    }
    return;
}
void OutLine::update()
{
    if(ElemTable->conformation!=conformation)
    {
        update_on_changed_geometry();
    }
    if(elementType == ElementType::TwoPhases)
    {
        update_two_phases();
    }
    else if(elementType == ElementType::SinglePhase)
    {
        update_single_phase();
    }
    else
    {
        assert(0);
    }
}

void OutLine::update_on_changed_geometry()
{
    //we here because geometry was changed
    flush_stats(false);

    //define pointers and give hits to compilers
    const int N=ElemTable->size();


    const double * RESTRICT el_dx=&(ElemTable->dx_[0][0]);
    const double * RESTRICT el_dy=&(ElemTable->dx_[1][0]);

    const double * RESTRICT el_x=&(ElemTable->coord_[0][0]);
    const double * RESTRICT el_y=&(ElemTable->coord_[1][0]);

    TI_ASSUME_ALIGNED(el_dx);
    TI_ASSUME_ALIGNED(el_dy);
    TI_ASSUME_ALIGNED(el_x);
    TI_ASSUME_ALIGNED(el_y);

    el_x_start.resize(N);
    el_x_stop.resize(N);
    el_y_start.resize(N);
    el_y_stop.resize(N);

    int * RESTRICT ix_start=&(el_x_start[0]);
    int * RESTRICT ix_stop=&(el_x_stop[0]);
    int * RESTRICT iy_start=&(el_y_start[0]);
    int * RESTRICT iy_stop=&(el_y_stop[0]);

    TI_ASSUME_ALIGNED(ix_start);
    TI_ASSUME_ALIGNED(ix_stop);
    TI_ASSUME_ALIGNED(iy_start);
    TI_ASSUME_ALIGNED(iy_stop);


    pileheight_by_elm.resize(N);
    max_kinergy_by_elm.resize(N);
    max_dynamic_pressure_by_elm.resize(N);
    cum_kinergy_by_elm.resize(N);

    double * RESTRICT m_pileheight_by_elm=&(pileheight_by_elm[0]);
    double * RESTRICT m_max_kinergy_by_elm=&(max_kinergy_by_elm[0]);
    double * RESTRICT m_max_dynamic_pressure_by_elm=&(max_dynamic_pressure_by_elm[0]);
    double * RESTRICT m_cum_kinergy_by_elm=&(cum_kinergy_by_elm[0]);

    TI_ASSUME_ALIGNED(m_pileheight_by_elm);
    TI_ASSUME_ALIGNED(m_max_kinergy_by_elm);
    TI_ASSUME_ALIGNED(m_cum_kinergy_by_elm);


    #pragma omp parallel
    {

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < N; ndx++)
        {
            m_pileheight_by_elm[ndx]=0.0;
            m_max_kinergy_by_elm[ndx]=0.0;
            m_max_dynamic_pressure_by_elm[ndx]=0.0;
            m_cum_kinergy_by_elm[ndx]=0.0;
        }

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < N; ndx++)
        {
            double xstart=el_x[ndx] - 0.5 * el_dx[ndx];
            double xstop=el_x[ndx] + 0.5 * el_dx[ndx];
            double ystart=el_y[ndx] - 0.5 * el_dy[ndx];
            double ystop=el_y[ndx] + 0.5 * el_dy[ndx];

            ix_start[ndx] = (int) ((xstart - xminmax[0]) / dx + 0.5);
            ix_stop[ndx]  = (int) ((xstop - xminmax[0]) / dx + 0.5);
            iy_start[ndx] = (int) ((ystart - yminmax[0]) / dy + 0.5);
            iy_stop[ndx]  = (int) ((ystop - yminmax[0]) / dy + 0.5);

            if(ix_start[ndx] < 0)
                ix_start[ndx] = 0;
            if(ix_stop[ndx] == ix_start[ndx])
            {
                ix_start[ndx] = (int) ((xstart - xminmax[0]) / dx);
                ix_stop[ndx] = ix_start[ndx] + 1;
            }
            if(ix_stop[ndx] > Nx)
                ix_stop[ndx] = Nx;

            if(iy_start[ndx] < 0)
                iy_start[ndx] = 0;
            if(iy_stop[ndx] == iy_start[ndx])
            {
                iy_start[ndx] = (int) ((ystart - yminmax[0]) / dy);
                iy_stop[ndx] = iy_start[ndx] + 1;
            }
            if(iy_stop[ndx] > Ny)
                iy_stop[ndx] = Ny;

        }
    }
    conformation=ElemTable->conformation;
}
void OutLine::update_single_phase()
{
    const ti_ndx_t N=pileheight_by_elm.size();

    int * RESTRICT adapted_=&(ElemTable->adapted_[0]);

    double * RESTRICT h=&(ElemTable->state_vars_[0][0]);
    double * RESTRICT hVx=&(ElemTable->state_vars_[1][0]);
    double * RESTRICT hVy=&(ElemTable->state_vars_[2][0]);

    TI_ASSUME_ALIGNED(adapted_);
    TI_ASSUME_ALIGNED(h);
    TI_ASSUME_ALIGNED(hVx);
    TI_ASSUME_ALIGNED(hVy);

    double * RESTRICT m_pileheight_by_elm=&(pileheight_by_elm[0]);
    double * RESTRICT m_max_kinergy_by_elm=&(max_kinergy_by_elm[0]);
    double * RESTRICT m_max_dynamic_pressure_by_elm=&(max_dynamic_pressure_by_elm[0]);
    double * RESTRICT m_cum_kinergy_by_elm=&(cum_kinergy_by_elm[0]);

    TI_ASSUME_ALIGNED(m_pileheight_by_elm);
    TI_ASSUME_ALIGNED(m_max_kinergy_by_elm);
    TI_ASSUME_ALIGNED(m_max_dynamic_pressure_by_elm);
    TI_ASSUME_ALIGNED(m_cum_kinergy_by_elm);

    if(numprocs>1)
    {
        #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < N; ndx++)
        {
            if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

            //update the record of maximum pileheight in the area covered by this element
            double ke = 0.0;
            double dp = 0.0;
            if (h[ndx] > 1.0E-04){
            	ke = 0.5 * (hVx[ndx] * hVx[ndx] + hVy[ndx] * hVy[ndx]) / h[ndx];
            	dp = ke / h[ndx];
            }

            m_cum_kinergy_by_elm[ndx] += ke;
            m_pileheight_by_elm[ndx] = max(m_pileheight_by_elm[ndx], h[ndx]);
            m_max_dynamic_pressure_by_elm[ndx] = max(m_max_dynamic_pressure_by_elm[ndx],dp);
            m_max_kinergy_by_elm[ndx] = max(m_max_kinergy_by_elm[ndx],ke);
        }
    }
    else
    {
        #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < N; ndx++)
        {
            //update the record of maximum pileheight in the area covered by this element
            double ke = 0.0;
            double dp = 0.0;
            if (h[ndx] > 1.0E-04){
            	ke = 0.5 * (hVx[ndx] * hVx[ndx] + hVy[ndx] * hVy[ndx]) / h[ndx];
            	dp = ke / h[ndx];
            }

            m_cum_kinergy_by_elm[ndx] += ke;
            m_pileheight_by_elm[ndx] = max(m_pileheight_by_elm[ndx], h[ndx]);
            m_max_dynamic_pressure_by_elm[ndx] = max(m_max_dynamic_pressure_by_elm[ndx],dp);
            m_max_kinergy_by_elm[ndx] = max(m_max_kinergy_by_elm[ndx],ke);
        }
    }
}
void OutLine::flush_stats(bool zero_old_arrays)
{
    //update thread local grids with previously collected data
    //hits to compiler
    const ti_ndx_t N=pileheight_by_elm.size();
    double * RESTRICT m_pileheight_by_elm=&(pileheight_by_elm[0]);
    double * RESTRICT m_max_kinergy_by_elm=&(max_kinergy_by_elm[0]);
    double * RESTRICT m_max_dynamic_pressure_by_elm=&(max_kinergy_by_elm[0]);
    double * RESTRICT m_cum_kinergy_by_elm=&(cum_kinergy_by_elm[0]);

    TI_ASSUME_ALIGNED(m_pileheight_by_elm);
    TI_ASSUME_ALIGNED(m_max_kinergy_by_elm);
    TI_ASSUME_ALIGNED(m_max_dynamic_pressure_by_elm);
    TI_ASSUME_ALIGNED(m_cum_kinergy_by_elm);

    int * RESTRICT ix_start=&(el_x_start[0]);
    int * RESTRICT ix_stop=&(el_x_stop[0]);
    int * RESTRICT iy_start=&(el_y_start[0]);
    int * RESTRICT iy_stop=&(el_y_stop[0]);

    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();
        double *m_pileheight=pileheight_loc[ithread];
        double *m_max_kinergy=max_kinergy_loc[ithread];
        double *m_max_dynamic_pressure=max_dynamic_pressure_loc[ithread];
        double *m_cum_kinergy=cum_kinergy_loc[ithread];

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < N; ndx++)
        {
            for(int iy = iy_start[ndx]; iy < iy_stop[ndx]; iy++)
            {
                for(int ix = ix_start[ndx]; ix < ix_stop[ndx]; ix++)
                {
                    m_cum_kinergy[iy*stride+ix] += m_cum_kinergy_by_elm[ndx];
                    if(m_pileheight_by_elm[ndx] > m_pileheight[iy*stride+ix])
                        m_pileheight[iy*stride+ix] = m_pileheight_by_elm[ndx];
                    if( m_max_kinergy_by_elm[ndx]> m_max_kinergy[iy*stride+ix])
                        m_max_kinergy[iy*stride+ix] = m_max_kinergy_by_elm[ndx];
                    if( m_max_dynamic_pressure_by_elm[ndx]> m_max_dynamic_pressure[iy*stride+ix])
                        m_max_dynamic_pressure[iy*stride+ix] = m_max_dynamic_pressure_by_elm[ndx];
                }
            }
        }
        if(zero_old_arrays)
        {
            #pragma omp for schedule(static,TITAN2D_DINAMIC_CHUNK)
            for(ti_ndx_t ndx = 0; ndx < pileheight_by_elm.size(); ndx++)
            {
                m_pileheight_by_elm[ndx]=0.0;
                m_max_kinergy_by_elm[ndx]=0.0;
                m_max_dynamic_pressure_by_elm[ndx]=0.0;
                m_cum_kinergy_by_elm[ndx]=0.0;
            }
        }
    }
}
void OutLine::combine_results_from_threads()
{
    flush_stats();

    #pragma omp parallel
    {
        for (int j = 0; j < threads_number; ++j)
        {
            #pragma omp for schedule(static,TITAN2D_DINAMIC_CHUNK)
            for (int i = 0; i < size; ++i)
            {
                cum_kinergy[i] += cum_kinergy_loc[j][i];
                pileheight[i] = max(pileheight[i], pileheight_loc[j][i]);
                max_kinergy[i] = max(max_kinergy[i], max_kinergy_loc[j][i]);
                max_dynamic_pressure[i] = max(max_dynamic_pressure[i], max_dynamic_pressure_loc[j][i]);

                pileheight_loc[j][i] = 0.0;
                max_kinergy_loc[j][i] = 0.0;
                max_dynamic_pressure_loc[j][i] = 0.0;
                cum_kinergy_loc[j][i] = 0.0;
            }
        }
    }
}
void OutLine::update_two_phases()
{
    const ti_ndx_t N=pileheight_by_elm.size();

    int * RESTRICT adapted_=&(ElemTable->adapted_[0]);

    double * RESTRICT h=&(ElemTable->state_vars_[0][0]);
    double * RESTRICT hVx_sol=&(ElemTable->state_vars_[1][0]);
    double * RESTRICT hVy_sol=&(ElemTable->state_vars_[2][0]);
    double * RESTRICT hVx_liq=&(ElemTable->state_vars_[4][0]);
    double * RESTRICT hVy_liq=&(ElemTable->state_vars_[5][0]);

    TI_ASSUME_ALIGNED(adapted_);
    TI_ASSUME_ALIGNED(h);
    TI_ASSUME_ALIGNED(hVx_sol);
    TI_ASSUME_ALIGNED(hVy_sol);
    TI_ASSUME_ALIGNED(hVx_liq);
    TI_ASSUME_ALIGNED(hVy_liq);

    double * RESTRICT m_pileheight_by_elm=&(pileheight_by_elm[0]);
    double * RESTRICT m_max_kinergy_by_elm=&(max_kinergy_by_elm[0]);
    double * RESTRICT m_max_dynamic_pressure_by_elm=&(max_dynamic_pressure_by_elm[0]);
    double * RESTRICT m_cum_kinergy_by_elm=&(cum_kinergy_by_elm[0]);

    TI_ASSUME_ALIGNED(m_pileheight_by_elm);
    TI_ASSUME_ALIGNED(m_max_kinergy_by_elm);
    TI_ASSUME_ALIGNED(m_max_dynamic_pressure_by_elm);
    TI_ASSUME_ALIGNED(m_cum_kinergy_by_elm);


    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(ti_ndx_t ndx = 0; ndx < N; ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!

        //@TODO check ke for two phases
        double ke = 0.0;
        double dp = 0.0;
        if (h[ndx] > 1.0E-04){
        	ke = 0.5 * (hVx_sol[ndx] * hVx_sol[ndx] + hVy_sol[ndx] * hVy_sol[ndx] + hVx_liq[ndx] * hVx_liq[ndx] + hVy_liq[ndx] * hVy_liq[ndx]) / h[ndx];
        	dp = ke / h[ndx];
        }

        m_cum_kinergy_by_elm[ndx] += ke;
        if(h[ndx] > m_pileheight_by_elm[ndx])
            m_pileheight_by_elm[ndx] = h[ndx];
        if(ke > m_max_kinergy_by_elm[ndx])
            m_max_kinergy_by_elm[ndx] = ke;
        if(dp > m_max_dynamic_pressure_by_elm[ndx])
            m_max_dynamic_pressure_by_elm[ndx] = dp;
    }
}

/*! this function outputs the maximum over time map of pileheights
 *  to the file pileheightrecord.xxxxxx
 */
void OutLine::output(MatProps* matprops_ptr, StatProps* statprops_ptr)
{
    double ENERGY_SCALE = matprops_ptr->scale.length * matprops_ptr->scale.gravity * matprops_ptr->scale.height;

    // output max over time pile-height
    {
        int ix, iy;
        FILE *fp;
        ostringstream filename1;

        filename1<<output_prefix<<"pileheightrecord."<<setw(6)<< setfill('0') <<internal<<statprops_ptr->runid<<std::ends;
        fp = fopen(filename1.str().c_str(), "wt");

        fprintf(fp, "Nx=%d: X={%20.14g,%20.14g}\n"
                "Ny=%d: Y={%20.14g,%20.14g}\n"
                "Pileheight=\n",
                Nx, xminmax[0] * matprops_ptr->scale.length, xminmax[1] * matprops_ptr->scale.length, Ny,
                yminmax[0] * matprops_ptr->scale.length, yminmax[1] * matprops_ptr->scale.length);
        for(iy = 0; iy < Ny; iy++)
        {
            for(ix = 0; ix < Nx - 1; ix++)
                fprintf(fp, "%g ", pileheight[iy*stride+ix] * matprops_ptr->scale.height);
            fprintf(fp, "%g\n", pileheight[iy*stride+ix] * matprops_ptr->scale.height);
        }
        fclose(fp);

        /*ostringstream filename2;
        filename2<<output_prefix<<"pileheightrecord."<<setw(6)<< setfill('0') <<internal<<statprops_ptr->runid<<".mat"<<std::ends;
        fp = fopen(filename2.str().c_str(), "wb");

        const char *map_name="max_pile_hight";
        const char *map_data="map_data";
        const char *map_title="map_title";
        const char *map_northern_edge="map_northern_edge";
        const char *map_southern_edge="map_southern_edge";
        const char *map_eastern_edge="map_eastern_edge";
        const char *map_western_edge="map_western_edge";

        //header format: {MOPT,mrows,ncols,imagf,namlen}
        //rows-y,cols-x
        unsigned int map_data_header[5]={0,(unsigned int)Ny,(unsigned int)Nx,0,(unsigned int)strlen(map_data)+1};
        unsigned int map_title_header[5]={51,1,(unsigned int)strlen(map_name)+1,0,(unsigned int)strlen(map_title)+1};
        unsigned int map_northern_edge_header[5]={0,1,1,0,(unsigned int)strlen(map_northern_edge)+1};
        unsigned int map_southern_edge_header[5]={0,1,1,0,(unsigned int)strlen(map_southern_edge)+1};
        unsigned int map_eastern_edge_header[5]={0,1,1,0,(unsigned int)strlen(map_eastern_edge)+1};
        unsigned int map_western_edge_header[5]={0,1,1,0,(unsigned int)strlen(map_western_edge)+1};

        double western_edge=xminmax[0] * matprops_ptr->scale.length;
        double eastern_edge=xminmax[1] * matprops_ptr->scale.length;
        double southern_edge=yminmax[0] * matprops_ptr->scale.length;
        double northern_edge= yminmax[1] * matprops_ptr->scale.length;

        double *buffer=new double[Ny];

        fwrite(map_data_header,sizeof(unsigned int),5,fp);
        fwrite(map_data,sizeof(char),strlen(map_data)+1,fp);

        for(ix = 0; ix < Nx; ix++)
        {
            //for(iy = 0; iy < Ny; iy++)
            int iy2=0;
            //why it is backward?
            for(iy = Ny-1; iy >=0; iy--,iy2++)
                buffer[iy2]=pileheight[iy*stride+ix] * matprops_ptr->scale.height;
            fwrite(buffer,sizeof(double),Ny,fp);
        }

        fwrite(map_title_header,sizeof(unsigned int),5,fp);
        fwrite(map_title,sizeof(char),strlen(map_title)+1,fp);
        fwrite(map_name,sizeof(char),strlen(map_name)+1,fp);

        fwrite(map_northern_edge_header,sizeof(unsigned int),5,fp);
        fwrite(map_northern_edge,sizeof(char),strlen(map_northern_edge)+1,fp);
        fwrite(&northern_edge,sizeof(double),1,fp);

        fwrite(map_southern_edge_header,sizeof(unsigned int),5,fp);
        fwrite(map_southern_edge, sizeof(char), strlen(map_southern_edge) + 1, fp);
        fwrite(&southern_edge, sizeof(double), 1, fp);

        fwrite(map_eastern_edge_header, sizeof(unsigned int), 5, fp);
        fwrite(map_eastern_edge, sizeof(char), strlen(map_eastern_edge) + 1, fp);
        fwrite(&eastern_edge, sizeof(double), 1, fp);

        fwrite(map_western_edge_header, sizeof(unsigned int), 5, fp);
        fwrite(map_western_edge, sizeof(char), strlen(map_western_edge) + 1, fp);
        fwrite(&western_edge, sizeof(double), 1, fp);


        fclose(fp);
        delete [] buffer;*/
    }

    //output max over time kinetic energy
    if(elementType==ElementType::SinglePhase)
    {
        int ix, iy;
        ostringstream filename;

        filename<<output_prefix<<"maxkerecord."<<setw(6)<< setfill('0') <<internal<<statprops_ptr->runid<<std::ends;
        FILE *fp = fopen(filename.str().c_str(), "wt");

        fprintf(fp, "Nx=%d: X={%20.14g,%20.14g}\n"
                "Ny=%d: Y={%20.14g,%20.14g}\n"
                "KineticEnergy=\n",
                Nx, xminmax[0] * matprops_ptr->scale.length, xminmax[1] * matprops_ptr->scale.length, Ny,
                yminmax[0] * matprops_ptr->scale.length, yminmax[1] * matprops_ptr->scale.length);
        for(iy = 0; iy < Ny; iy++)
        {
            for(ix = 0; ix < Nx - 1; ix++)
                fprintf(fp, "%g ", max_kinergy[iy*stride+ix] * ENERGY_SCALE);
            fprintf(fp, "%g\n", max_kinergy[iy*stride+ix] * ENERGY_SCALE);
        }
        fclose(fp);
    }

    //dynamic_pressure
    if(elementType==ElementType::SinglePhase)
    {
        int ix, iy;
        ostringstream filename;

        filename<<output_prefix<<"max_dynamic_pressure_record."<<setw(6)<< setfill('0') <<internal<<statprops_ptr->runid<<std::ends;
        FILE *fp = fopen(filename.str().c_str(), "wt");

        //!todo dynamic_pressure scale?
        fprintf(fp, "Nx=%d: X={%20.14g,%20.14g}\n"
                "Ny=%d: Y={%20.14g,%20.14g}\n"
                "dynamic_pressure=\n",
                Nx, xminmax[0] * matprops_ptr->scale.length, xminmax[1] * matprops_ptr->scale.length, Ny,
                yminmax[0] * matprops_ptr->scale.length, yminmax[1] * matprops_ptr->scale.length);
        for(iy = 0; iy < Ny; iy++)
        {
            for(ix = 0; ix < Nx - 1; ix++)
                fprintf(fp, "%g ", max_dynamic_pressure[iy*stride+ix] * ENERGY_SCALE/matprops_ptr->scale.height);
            fprintf(fp, "%g\n", max_dynamic_pressure[iy*stride+ix] * ENERGY_SCALE/matprops_ptr->scale.height);
        }
        fclose(fp);
    }

    // output cummulative kinetic-energy
    if(elementType==ElementType::SinglePhase)
    {
        int ix, iy;
        ostringstream filename;

        filename<<output_prefix<<"cumkerecord."<<setw(6)<< setfill('0') <<internal<<statprops_ptr->runid<<std::ends;
        FILE *fp = fopen(filename.str().c_str(), "wt");

        fprintf(fp, "Nx=%d: X={%20.14g,%20.14g}\n"
                "Ny=%d: Y={%20.14g,%20.14g}\n"
                "KineticEnergy=\n",
                Nx, xminmax[0] * matprops_ptr->scale.length, xminmax[1] * matprops_ptr->scale.length, Ny,
                yminmax[0] * matprops_ptr->scale.length, yminmax[1] * matprops_ptr->scale.length);
        for(iy = 0; iy < Ny; iy++)
        {
            for(ix = 0; ix < Nx - 1; ix++)
                fprintf(fp, "%g ", cum_kinergy[iy*stride+ix] * ENERGY_SCALE);
            fprintf(fp, "%g\n", cum_kinergy[iy*stride+ix] * ENERGY_SCALE);
        }
        fclose(fp);
    }
    // output elevation data
    {
        int ix, iy;
        ostringstream filename;

        filename<<output_prefix<<"elevation.grid"<<std::ends;
        FILE *fp = fopen(filename.str().c_str(), "wt");

        fprintf(fp, "Nx=%d: X={%20.14g,%20.14g}\n"
                "Ny=%d: Y={%20.14g,%20.14g}\n"
                "Pileheight=\n",
                Nx, xminmax[0] * matprops_ptr->scale.length, xminmax[1] * matprops_ptr->scale.length, Ny,
                yminmax[0] * matprops_ptr->scale.length, yminmax[1] * matprops_ptr->scale.length);

        double yy, xx, res = dx + dy, elevation;
        int ierr;
        for(iy = 0; iy < Ny; iy++)
        {
            yy = ((iy + 0.5) * dy + yminmax[0]) * matprops_ptr->scale.length;
            for(ix = 0; ix < Nx - 1; ix++)
            {
                xx = ((ix + 0.5) * dx + xminmax[0]) * matprops_ptr->scale.length;
                ierr = Get_elevation(res, xx, yy, elevation);
                fprintf(fp, "%g ", elevation);
            }
            xx = ((ix + 0.5) * dx + xminmax[0]) * matprops_ptr->scale.length;
            ierr = Get_elevation(res, xx, yy, elevation);
            fprintf(fp,"%g\n",elevation);
        }
        fclose(fp);
    }
    return;
}

//! this function reads in the previous map of maximum throughout time pileheight stored in the file pileheightrecord.xxxxxx during restart
void OutLine::reload(MatProps* matprops_ptr, StatProps* statprops_ptr)
{
    int ix, iy;
    char filename[256];
    sprintf(filename, "pileheightrecord.%06d", statprops_ptr->runid);
    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
        printf("pileheightrecord.%06d can not be found.\nRestarting from zero instead\n", statprops_ptr->runid);
    else
    {

        int Nxtemp, Nytemp;
        double xminmaxtemp[2], yminmaxtemp[2];
        fscanf(fp, "Nx=%d: X={%lf,%lf}\nNy=%d: Y={%lf,%lf}\nPileheight=\n", &Nxtemp, (xminmaxtemp + 0),
               (xminmaxtemp + 1), &Nytemp, (yminmaxtemp + 0), (yminmaxtemp + 1));

        if((Nxtemp == Nx) && (fabs(xminmaxtemp[0] / matprops_ptr->scale.length - xminmax[0])
                <= fabs(xminmax[0]) / 10000000000.0)
           && (fabs(xminmaxtemp[1] / matprops_ptr->scale.length - xminmax[1]) <= fabs(xminmax[1]) / 10000000000.0)
           && (Nytemp == Ny)
           && (fabs(yminmaxtemp[0] / matprops_ptr->scale.length - yminmax[0]) <= fabs(yminmax[0]) / 10000000000.0)
           && (fabs(yminmaxtemp[1] / matprops_ptr->scale.length - yminmax[1]) <= fabs(yminmax[1]) / 10000000000.0))
        {

            for(iy = 0; iy < Ny; iy++)
                for(ix = 0; ix < Nx; ix++)
                {
                    fscanf(fp, "%lf", pileheight[iy] + ix);
                    pileheight[iy*stride+ix] /= matprops_ptr->scale.height;
                }
        }
        else
        {
            printf("the pileheightrecord.%06d that is present does not match the restart file.\nRestarting from zero instead\n",
                   statprops_ptr->runid);
            printf("Nx=%d Nxtemp=%d\n", Nx, Nxtemp);
            printf("Ny=%d Nytemp=%d\n", Ny, Nytemp);
            printf("xmin=%20.14g xmintemp=%20.14g  xmax=%20.14g, xmaxtemp=%20.14g\n",
                   xminmax[0] * matprops_ptr->scale.length, xminmaxtemp[0], xminmax[1] * matprops_ptr->scale.length,
                   xminmaxtemp[1]);
            printf("ymin=%20.14g ymintemp=%20.14g  ymax=%20.14g, ymaxtemp=%20.14g\n",
                   yminmax[0] * matprops_ptr->scale.length, yminmaxtemp[0], yminmax[1] * matprops_ptr->scale.length,
                   yminmaxtemp[1]);
            assert(0);
        }

        fclose(fp);
    }

    return;
}

void OutLine::h5write(H5::CommonFG *parent, string group_name)
{
    H5::Group group(parent->createGroup(group_name));

    if(enabled)
        combine_results_from_threads();

    TiH5_writeIntAttribute(group, conformation);
    TiH5_writeIntAttribute(group, Nx);
    TiH5_writeIntAttribute(group, Ny);
    TiH5_writeIntAttribute(group, stride);
    TiH5_writeIntAttribute(group, size);
    TiH5_writeDoubleAttribute(group, dx);
    TiH5_writeDoubleAttribute(group, dy);
    TiH5_writeDoubleAttribute(group, xminmax[2]);
    TiH5_writeDoubleAttribute(group, yminmax[2]);

    TiH5_writeBoolAttribute(group, enabled);
    TiH5_writeScalarDataTypeAttribute(group,init_size,datatypeOutLineInitSize);
    TiH5_writeIntAttribute(group, max_linear_size);
    TiH5_writeStringAttribute(group, output_prefix);
    if(enabled)
    {
        TiH5_writeArray2DDataSet(group, pileheight,Ny,stride);
        TiH5_writeArray2DDataSet(group, max_kinergy,Ny,stride);
        TiH5_writeArray2DDataSet(group, max_dynamic_pressure,Ny,stride);
        TiH5_writeArray2DDataSet(group, cum_kinergy,Ny,stride);
    }
    TiH5_writeScalarDataTypeAttribute(group, elementType, datatypeElementType);


}
void OutLine::h5read(const H5::CommonFG *parent, const  string group_name)
{
    H5::Group group(parent->openGroup(group_name));

    TiH5_readScalarDataTypeAttribute(group, elementType, datatypeElementType);

    TiH5_readIntAttribute(group, conformation);
    TiH5_readIntAttribute(group, Nx);
    TiH5_readIntAttribute(group, Ny);
    TiH5_readIntAttribute(group, stride);
    TiH5_readIntAttribute(group, size);
    TiH5_readDoubleAttribute(group, dx);
    TiH5_readDoubleAttribute(group, dy);
    TiH5_readDoubleAttribute(group, xminmax[2]);
    TiH5_readDoubleAttribute(group, yminmax[2]);

    TiH5_readBoolAttribute(group, enabled);
    TiH5_readScalarDataTypeAttribute(group,init_size,datatypeOutLineInitSize);
    TiH5_readIntAttribute(group, max_linear_size);
    TiH5_readStringAttribute(group, output_prefix);
    if(enabled)
    {
        alloc_arrays();

        TiH5_readArray2DDataSet(group, pileheight,Ny,stride);
        TiH5_readArray2DDataSet(group, max_kinergy,Ny,stride);
        TiH5_readArray2DDataSet(group, max_dynamic_pressure,Ny,stride);
        TiH5_readArray2DDataSet(group, cum_kinergy,Ny,stride);

        update_on_changed_geometry();
    }
}

