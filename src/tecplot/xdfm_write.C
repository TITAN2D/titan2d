/*
 * =====================================================================================
 *
 *       Filename:  xdfm_write.cc
 *
 *    Description:  XDMF writer generates paraview readable output
 *                  1. writes titan output to HDF5 file 
 *                  2. writes metadata to XML file
 *                  TODO 1: Need to re-refine UNREFINED cells 
 *
 *        Created:  05/15/2007 12:17:29 PM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */
#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <hdf5.h>
#include "../header/ticore/omp_mpi.hpp"

#include "../header/elenode.hpp"
#include <GisApi.h>
#include <hd5calls.h>

//just an innocent declaration. defined at end of the file
double interpolate_elv(ElementsHashTable *, NodeHashTable *, Element *, int);
void xdmf_fopen(char *);
void xdmf_fclose(ofstream &);

int write_xdmf_two_phases(ElementsHashTable *El_Table, NodeHashTable *NodeTable, TimeProps *timeprops_ptr, MatProps *matprops_ptr,
               MapNames *mapnames, const int mode, const char * output_prefix)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    ofstream xdmf;
    char xmlf_filename[256];
    sprintf(xmlf_filename, "%s_xdmf_p%04d.xmf", output_prefix,myid);

    if(mode == XDMF_ONLYCLOSE)
    {
        xdmf.open(xmlf_filename, ios::app);
        xdmf_fclose(xdmf);
        xdmf.close();
        return 0;
    }

    //if Need to have generic form, do vector of vectors
    vector<double> pheight, xmom, ymom, xcoord, ycoord, zcoord,tot_eros;
    vector<int> C1, C2, C3, C4;
    int num_nodes = 0, num_elm = 0;
    int i, j, k;
    double elevation;
    int conn[4];
    
    // scaling factor for the momentums
    double momentum_scale = matprops_ptr->scale.height
            * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity));
    
    /* scan HashTable and store coordinates and variables in vectors */
    Element *EmTemp = NULL;
    Node *NodeTemp = NULL;
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!(EmTemp->refined_flag()))
            {
                pheight.push_back(EmTemp->state_vars(0) * matprops_ptr->scale.height);
                xmom.push_back(EmTemp->state_vars(1) * momentum_scale);
                ymom.push_back(EmTemp->state_vars(2) * momentum_scale);
                tot_eros.push_back(EmTemp->TOT_EROS()*matprops_ptr->scale.height);
                num_elm++;
                for(j = 0; j < 4; j++)
                {
                    NodeTemp = (Node *) NodeTable->lookup(EmTemp->node_key(j));
                    xcoord.push_back(NodeTemp->coord(0) * matprops_ptr->scale.length);
                    ycoord.push_back(NodeTemp->coord(1) * matprops_ptr->scale.length);
                    conn[j] = num_nodes;
                    // Side-Corner node will have discontinuity in elevation
                    if(NodeTemp->info() != S_C_CON)
                    {
                        elevation = NodeTemp->elevation() * matprops_ptr->scale.length;
                        zcoord.push_back(elevation);
                    }
                    // hence <get a fake one> interpolate from neighbor
                    else
                    {
                        elevation = interpolate_elv(El_Table, NodeTable, EmTemp, j) * matprops_ptr->scale.length;
                        zcoord.push_back(elevation);
                    }
                    num_nodes++;
                }
            }
        }
    }
    

    /* 
     * Write HDF5 File 
     */
    char hdf5file[256];
    sprintf(hdf5file, "%s/xdmf_p%04d_i%08d.h5", output_prefix, myid, timeprops_ptr->iter);
    hid_t h5fid = GH5_fopen(hdf5file, 'n');
    
    //allocate memory for xyz points
    double *xyz = new double[num_nodes * 3];
    int ix = 0;
    for(i = 0; i < num_nodes; i++)
    {
        xyz[ix + 0] = xcoord[i];
        xyz[ix + 1] = ycoord[i];
        xyz[ix + 2] = zcoord[i];
        ix += 3;
    }
    
    //allocate memory for conncetions
    int *conns = new int[num_elm * 4];
    int icon = 0;
    for(i = 0; i < (4 * num_elm); i++)
    {
        conns[i] = i;
    }
    GH5_write_mesh_data(h5fid, num_elm, num_nodes, conns, xyz);
    
    /*  FREE ALLOCATED MEMORY */
    // xyz memory
    delete[] xyz;
    // conns memory
    delete[] conns;
    
    /* Write Variable Data */
    double *vars = new double[num_elm];
    
    //Pile Height
    copy(pheight.begin(), pheight.end(), vars);
    GH5_write_state_vars(h5fid, num_elm, vars, "PILE_HEIGHT");
    
    //X-Momentum
    copy(xmom.begin(), xmom.end(), vars);
    GH5_write_state_vars(h5fid, num_elm, vars, "XMOMENTUM");
    
    //Y-Momentum
    copy(ymom.begin(), ymom.end(), vars);
    GH5_write_state_vars(h5fid, num_elm, vars, "YMOMENTUM");

    //Total Erosion
    copy(tot_eros.begin(), tot_eros.end(), vars);
    GH5_write_state_vars(h5fid, num_elm, vars, "EROSION");
    
    delete[] vars;
    GH5_fclose(h5fid);
    
    /* generate XML file if required */
    if(mode == XDMF_NEW)
        xdmf_fopen(xmlf_filename);
    xdmf.open(xmlf_filename, ios::app);
    
    /*xmdf <<  add comments about when the data was generated */
    xdmf << "<Grid Name=\"" << mapnames->gis_map << "\">" << endl;
    xdmf << "<Time Value=\"" << timeprops_ptr->timesec() << "\" />" << endl;
    // connectivity data
    xdmf << "<Topology Type=\"QUADRILATERAL\" Dimensions=\"" << num_elm << "\" Order=\"0 3 2 1\">" << endl;
    xdmf << "<DataItem Name=\"Connections\" DataType=\"Int\" Precision=\"4\"" << endl;
    xdmf << "Dimensions=\"" << num_elm << " 4\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Mesh/Connections" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Topology>" << endl;
    // grid-points
    xdmf << "<Geometry Type=\"XYZ\">" << endl;
    xdmf << "<DataItem Name=\"Coordinates\" DataType=\"Float\" Precision=\"8\"" << endl;
    xdmf << "Dimensions=\"" << num_nodes << " 3\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Mesh/Points" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Geometry>" << endl;
    // pile-height
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Pile Height\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elm << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/PILE_HEIGHT" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    // x-momentum
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"X_Momentum\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elm << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/XMOMENTUM" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    // y-momentum
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Y_Momentum\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elm << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/YMOMENTUM" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    // Total-Erosion
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Total_Erosion\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elm << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/EROSION" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    xdmf << "</Grid>" << endl;
    
    if(mode == XDMF_CLOSE)
    {
        xdmf_fclose(xdmf);
    }
    xdmf.close();  // close till next output is written
    
    return 0;
}
int write_xdmf_single_phase(ElementsHashTable *El_Table, NodeHashTable *NodeTable, TimeProps *timeprops_ptr, MatProps *matprops_ptr,
               MapNames *mapnames, const int mode, const char * output_prefix)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    ofstream xdmf;
    char xdmf_filename[256];
    sprintf(xdmf_filename, "%s_xdmf_p%04d.xmf", output_prefix,myid);

    if(mode == XDMF_ONLYCLOSE)
    {
        xdmf.open(xdmf_filename, ios::app);
        xdmf_fclose(xdmf);
        xdmf.close();
        return 0;
    }

    int i, j, k;
    double elevation;
    
    // scaling factor for the momentums
    double momentum_scale = matprops_ptr->scale.height
            * sqrt(matprops_ptr->scale.length * (matprops_ptr->scale.gravity));
    static double time_prev;
    
    /* generate XML file if required */

    if(mode == XDMF_NEW)
        xdmf_fopen(xdmf_filename);
    xdmf.open(xdmf_filename, ios::app);
    
    /* save time-value to compare at end */
    if(mode != XDMF_CLOSE)
        time_prev = timeprops_ptr->timesec();
    
    /* don't rewrite the last time-step */
    /* just close and return */
    if((mode == XDMF_CLOSE) && (timeprops_ptr->timesec() - time_prev < 0.0001))
    {
        xdmf_fclose(xdmf);
        xdmf.close();
        return 0;
    }
    
    Element *EmTemp = NULL;
    Node *NodeTemp = NULL;
    // reset connection-ids
    //@NodesSingleLoop
    for(i = 0; i < NodeTable->elenode_.size(); i++)
    {
        if(NodeTable->status_[i]>=0)
            NodeTable->elenode_[i].connection_id(-1);
    }
    
    /* Add connection information to each node */
    int num_elem = 0;
    int num_node = 0;
    int id = 0;
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
            {
                num_elem++;
                for(j = 0; j < 4; j++)
                {
                    NodeTemp = (Node *) NodeTable->lookup(EmTemp->node_key(j));
                    if(NodeTemp->connection_id() < 0)
                    {
                        NodeTemp->connection_id(id++);
                        num_node++;
                    }
                }
            }
        }
    }
    
    /* allocate memory for connections and element verticies */
    double * xyz = new double[num_node * 3];
    int * connections = new int[num_elem * 4];
    
    /* and other output variables */
    double * pheight = new double[num_elem];
    double * xmom = new double[num_elem];
    double * ymom = new double[num_elem];
    
    /* scan HashTable and store coordinates and variables in vectors */
    int ielm = 0;
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int jelm = 0; jelm < bucket[ibuck].ndx.size(); jelm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[jelm]]);
            if(!(EmTemp->refined_flag()))
            {
                pheight[ielm] = EmTemp->state_vars(0) * matprops_ptr->scale.height;
                xmom[ielm] = EmTemp->state_vars(1) * momentum_scale;
                ymom[ielm] = EmTemp->state_vars(2) * momentum_scale;
                for(j = 0; j < 4; j++)
                {
                    NodeTemp = (Node *) NodeTable->lookup(EmTemp->node_key(j));
                    int inode = NodeTemp->connection_id();
                    xyz[3 * inode] = NodeTemp->coord(0) * matprops_ptr->scale.length;
                    xyz[3 * inode + 1] = NodeTemp->coord(1) * matprops_ptr->scale.length;
                    connections[4 * ielm + j] = inode;
                    // Side-Corner node will have discontinuity in elevation
                    if(NodeTemp->info() != S_C_CON)
                    {
                        elevation = NodeTemp->elevation() * matprops_ptr->scale.length;
                        xyz[3 * inode + 2] = elevation;
                    }
                    // hence <get a fake one> interpolate from neighbor
                    else
                    {
                        elevation = interpolate_elv(El_Table, NodeTable, EmTemp, j) * matprops_ptr->scale.length;
                        xyz[3 * inode + 2] = elevation;
                    }
                }
                ielm++;
            }
        }
    }
    
    /* 
     * Write HDF5 File 
     */
    char hdf5file[256];
    sprintf(hdf5file, "%s/xdmf_p%04d_i%08d.h5", output_prefix, myid, timeprops_ptr->iter);
    hid_t h5fid = GH5_fopen(hdf5file, 'n');
    
    // Write Mesh data
    GH5_write_mesh_data(h5fid, num_elem, num_node, connections, xyz);
    
    /*  FREE ALLOCATED MEMORY */
    // xyz memory
    delete[] xyz;
    // conns memory
    delete[] connections;
    
    /* Write Variable Data */
    //Pile Height
    GH5_write_state_vars(h5fid, num_elem, pheight, "PILE_HEIGHT");
    
    //X-Momentum
    GH5_write_state_vars(h5fid, num_elem, xmom, "XMOMENTUM");
    
    //Y-Momentum
    GH5_write_state_vars(h5fid, num_elem, ymom, "YMOMENTUM");
    
    GH5_fclose(h5fid);
    
    // free-memory
    delete[] pheight;
    delete[] xmom;
    delete[] ymom;
    
    /*xdmf <<  add comments about when the data was generated */
    xdmf << "<Grid Name=\"" << mapnames->gis_map << "\">" << endl;
    xdmf << "<Time Value=\"" << timeprops_ptr->timesec() << "\" />" << endl;
    // connectivity data
    xdmf << "<Topology Type=\"QUADRILATERAL\" Dimensions=\"" << num_elem << "\" Order=\"0 3 2 1\">" << endl;
    xdmf << "<DataItem Name=\"Connections\" DataType=\"Int\" Precision=\"4\"" << endl;
    xdmf << "Dimensions=\"" << num_elem << " 4\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Mesh/Connections" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Topology>" << endl;
    // grid-points
    xdmf << "<Geometry Type=\"XYZ\">" << endl;
    xdmf << "<DataItem Name=\"Coordinates\" DataType=\"Float\" Precision=\"8\"" << endl;
    xdmf << "Dimensions=\"" << num_node << " 3\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Mesh/Points" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Geometry>" << endl;
    // pile-height
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Pile Height\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elem << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/PILE_HEIGHT" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    // x-momentum
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"X_Momentum\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elem << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/XMOMENTUM" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    // y-momentum
    xdmf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Y_Momentum\">" << endl;
    xdmf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
    xdmf << "Dimensions=\"" << num_elem << " 1\" Format=\"HDF\">" << endl;
    xdmf << "\t\t" << hdf5file << ":/Properties/YMOMENTUM" << endl;
    xdmf << "</DataItem>" << endl;
    xdmf << "</Attribute>" << endl;
    xdmf << "</Grid>" << endl;
    
    if(mode == XDMF_CLOSE)
        xdmf_fclose(xdmf);
    
    xdmf.close();  // close till next output is written
    
    return 0;
}
/* open an XML file for Metadata and write header information */
void xdmf_fopen(char *filename)
{
    
    ofstream xmf(filename, ios::out);
    if(!xmf.is_open())
    {
        cout << "ERROR: unable to open file for writing," << endl << "check disc space or check permissions" << endl;
        exit(1);
    }
    xmf << "<?xml version=\"1.0\" ?>" << endl;
    xmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]> " << endl;
    xmf << "<Xdmf Version=\"2.0\">" << endl;
    xmf << "<Domain>" << endl;
    xmf << "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << endl;
    xmf.close();
    return;
}

/* close XML header file */
void xdmf_fclose(ofstream & xmlf)
{
    xmlf << "</Grid>" << endl;
    xmlf << "</Domain>" << endl;
    xmlf << "</Xdmf>" << endl;
}

/*
 * interpolate_elv is exactly the same code that is
 * being used in meshplotter(...), in tecplot.C
 */
double interpolate_elv(ElementsHashTable *El_Table, NodeHashTable *NodeTable, Element *EmTemp, int j)
{
    double elev;
    int neighside, mynode;
    if(j == EmTemp->which_son() + 1)
    {
        mynode = j - 1;
        neighside = j;
    }
    else if(j == EmTemp->which_son() - 1)
    {
        mynode = j + 1;
        neighside = j - 1;
        if(neighside < 0)
            neighside = 3;
    }
    else if(EmTemp->which_son() == 0)
    {
        mynode = 0;
        neighside = 2;
    }
    else if(EmTemp->which_son() == 3)
    {
        mynode = 3;
        neighside = 0;
    }
    else
    {
        mynode = 1;
        neighside = 1;
    }
    Node* NodeTemp2 = (Node*) NodeTable->lookup(EmTemp->node_key(mynode));
    elev = 0.5 * NodeTemp2->elevation();
    Element* EmTemp2 = (Element*) El_Table->lookup((EmTemp->neighbor(neighside)));
    NodeTemp2 = (Node*) NodeTable->lookup(EmTemp2->node_key(j));
    elev += 0.5 * NodeTemp2->elevation();
    return elev;
}
