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
 * $Id: node.h 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifndef NODE_H
#define NODE_H

#include "properties.h"
#include "constant.h"
#include "hashtab.h"
#include "struct.h"
#include "sfc.h"

class TimeProps;
class FluxProps;

class Node
{
    friend class NodeHashTable;
protected:
    //nodes constructors should not be called directly but by NodeHashTable::createAddNode methods
    //! this is the constructor that creates a node when the initial grid is read in
    Node(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr)
    {
        init(keyi, coordi, matprops_ptr);
    }
        //! this is the constructor that creates bubble and edge nodes for son Elements when the father Element is refined
    Node(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr)/*for refined*/
    {
        init(keyi, coordi, inf, ord, matprops_ptr);
    }


    
    //! this is the node constructor that is called in construct_el() in update_element_info.C
    Node(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada)
    {
        init(keyi, coordi, inf, ord, elev, yada);
    }
  
    //! this is the constructor that recreates/restores a node that was saved in a restart file.
    Node(FILE* fp, MatProps* matprops_ptr) //for restart
    {
        init(fp, matprops_ptr);
    }
    public:  
    ~Node();
public:
    void init(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr);

    //! this is the constructor that creates bubble and edge nodes for son Elements when the father Element is refined
    void init(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr);/*for refined*/
    
    //! this is the node constructor that is called in construct_el() in update_element_info.C
    void init(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada);

    //! this is the constructor that recreates/restores a node that was saved in a restart file.
    void init(FILE* fp, MatProps* matprops_ptr); //for restart
    
    //! this function writes all of one Node's data necessary for restart to a file in a single fwrite statement
    void save_node(FILE* fp); //for restart

    //! this function returns the node type, the options are listed in constant.h and include: NODEINIT, CORNER, BUBBLE, SIDE, CONSTRAINED, S_C_CON, S_S_CON, ASSIGNED,and UNASSIGNED.
    int info() const {return info_;}
    //! this function sets the node type, the options are listed in constant.h and include: NODEINIT, CORNER, BUBBLE, SIDE, CONSTRAINED, S_C_CON, S_S_CON, ASSIGNED,and UNASSIGNED.
    void info(int in){info_ = in;}

    //! this function returns the node key, a key is a single number that is 2 unsigned variables long and is used to access the pointer to a Node or Element through the HashTable
    const SFC_Key& key() const {return key_;}
    void set_key(const SFC_Key& new_key){key_=new_key;}

    //! this function returns the global x and y coordinates of the node, in the finite difference version of Titan this is not always reliable, use the coordinates of the element instead.  It is reliable in the Discontinuous Galerkin version of Titan however.
    double coord(int idim) const {return coord_[idim];}
    void coord(int idim, double new_crd){coord_[idim]=new_crd;}   

    //! this is legacy afeapi and is not used
    int order(){return order_;}
    //! this is legacy afeapi and is not used
    void order(int i){order_ = i;}

    //! this function sets the node information and order, node order is legacy afeapi but node information is currently used, this function is called in update_element_info.C, another distict function with a similar name refined_neighbor::set_parameters also existis and is used in updatenei.C, these should not be confused
    void set_parameters(int inf, int ord){info(inf);order(ord);}

    //! this function returns the id of a node, it is used in repartitioning, 
    int id() const {return id_;}
    //! this function sets the id of a node, it is used in repartitioning, 
    void id(int id_in){id_ = id_in;}


    //! this function returns the vector of fluxes stored in an edge node between elements 
    double flux(int idim) const {return flux_[idim];}
    void flux(int idim, double value){flux_[idim]=value;}

    //! this function zeros the flux used during refinement, this is only distinct from the regular flux if flux velocity is being zero'd because the experimental stopping criteria says it should be.  This feature is disabled by default.  Keith implemented it.
    void zero_flux(){for(int i = 0; i < NUM_STATE_VARS; i++){flux(i,0.0);refinementflux(i, 0.0);}}
    
    double refinementflux(int idim) const {return refinementflux_[idim];}
    void refinementflux(int idim, double value){refinementflux_[idim]=value;}

    //! this function returns the elevation of a node, in the finite difference version of Titan this is not always reliable, use the elevation of the element instead.  It is reliable in the Discontinuous Galerkin version of Titan however.
    double elevation() const {return elevation_;}
    double & elevation_ref() {return elevation_;}
    //! this function sets the elevation of a node
    void elevation(MatProps* matprops_ptr);
    void elevation(double new_elevation){elevation_=new_elevation;}

    //! this function stores the number of elements associated with this node
    void num_assoc_elem(int numin){num_assoc_elem_ = numin;}
    //! this function returns the number of elements associated with this node
    int num_assoc_elem() const {return num_assoc_elem_;}

    //! set connection id
    void connection_id(int id){connection_id_ = id;}
    //! get connection id 
    int connection_id() const {return connection_id_;}
    
protected:
    //! used in delete_unused_nodes_and_elements() function 
    int id_;

    //! the number of associated elements, it is used in extraneous node 
    //deletion and debugging function AssertMeshErrorFree()
    int num_assoc_elem_;

    //! says what type of node this is see the comments of Node::get_info()
    int info_;

    //! this is legacy afeapi and is not important though it would involve quite a bit of work to remove because it occurs frequently in Titan
    int order_;

    //! the global x and y coordinates of the node
    double coord_[DIMENSION];

    //! this is the node key, a key is a single number that is 2 unsigned variables long and is used to access the pointer to a Node or Element through the HashTable
    SFC_Key key_;

    //! this elevation should currently be the GIS elevation at the finest "scale"
    double elevation_;

    //! these are the so called "regular fluxes" that is the ones that are used to update the elements, assume that element normal is parallel to either the x or y axis, Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default
    double flux_[MAX_NUM_STATE_VARS];

    //! the "refinement flux" is necessary when using the stopping criteria to reset the "regular" fluxes to what they would be if velocity was zero in the cell(s) involved.  The refinement flux is what the flux would have been if it had not been reset, they are needed since refinement is based on fluxes (and also pileheight gradient but that's not relevant here) Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default.
    double refinementflux_[MAX_NUM_STATE_VARS];

    //! node number for connection data -- varies with adaptation
    int connection_id_;
};


#endif
