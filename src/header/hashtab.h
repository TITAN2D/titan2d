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
 * $Id: hashtab.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

/* Hash table */
/* Every table can process NBUCKETS entries */
/* All buckets with the same entry are linked together*/
/* The pointer to the first bucket is stored in the table */
/*---------------------------------------------------------*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#include "constant.h"


#define HASHTABLE_LOOKUP_LINSEARCH 8

class Element;
class Node;
class BC;
class MatProps;


class HashEntryLine
{
public:
    vector<SFC_Key> key; //!<key: object key word
    vector<ti_ndx_t> ndx; //index of record
    
    //! find position of m_key in key
    ti_ndx_t lookup_local_ndx(const SFC_Key& m_key)
    {
        int size = key.size();

        SFC_Key *keyArr = &(key[0]);
        int i;

        if(size == 0)
            return ti_ndx_doesnt_exist;
        if(m_key < keyArr[0])
            return ti_ndx_doesnt_exist;
        if(m_key > keyArr[size - 1])
            return ti_ndx_doesnt_exist;

        if(size < HASHTABLE_LOOKUP_LINSEARCH)
        {
            for(i = 0; i < size; i++)
            {
                if(m_key == keyArr[i])
                {
                    return i;
                }
            }
        }
        else
        {
            int i0, i1, i2;
            i0 = 0;
            i1 = size / 2;
            i2 = size - 1;
            while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
            {
                if(m_key > keyArr[i1])
                {
                    i0 = i1 + 1;
                    i1 = (i0 + i2) / 2;
                }
                else
                {
                    i2 = i1;
                    i1 = (i0 + i2) / 2;
                }
            }
            for(i = i0; i <= i2; i++)
            {
                if(m_key == keyArr[i])
                {
                    return i;
                }
            }
        }
        return ti_ndx_doesnt_exist;
    }
    
};

enum ContentStatus {CS_Permanent=0,CS_Removed=-1,CS_Added=1};

template<typename T>
class HashTable
{
public:
    double doublekeyrange[2];
    double hashconstant;
    double Xrange[2];
    double Yrange[2];
    double invdxrange;
    double invdyrange;
    
    int NBUCKETS;
    int ENTRIES;
    
    vector<HashEntryLine> bucket;
    
    tivector<SFC_Key> key_;
    tivector<ContentStatus> status_;
    tivector<T> elenode_;
    
public:
    HashTable(double *doublekeyrangein, int size, double XR[], double YR[],tisize_t reserved_size);
    ~HashTable(){}
    
    int hash(const SFC_Key& keyi) const
    {
        //Keith made this change 20061109; and made hash an inline function
        /* NBUCKETS*2 is NBUCKETS*integer integer is empirical could be 1
         return (((int) ((key[0]*doublekeyrange[1]+key[1])/
         (doublekeyrange[0]*doublekeyrange[1]+doublekeyrange[1])*
         NBUCKETS*2+0.5) )%NBUCKETS);
         */
    #if USE_ARRAY_SFC_KEY
        return (((int) ((key.key[0] * doublekeyrange[1] + key.key[1]) * hashconstant + 0.5)) % NBUCKETS);
    #else
        unsigned oldkey[KEYLENGTH];
        SET_OLDKEY(oldkey,keyi);
        return (((int) ((oldkey[0] * doublekeyrange[1] + oldkey[1]) * hashconstant + 0.5)) % NBUCKETS);
        //return (((int) (keyi * hashconstant + 0.5)) % NBUCKETS);
    #endif
    }
    
    ti_ndx_t lookup_ndx(const SFC_Key& keyi);
    T* lookup(const SFC_Key& keyi);
protected:
    ti_ndx_t add_ndx(const SFC_Key& keyi);
    T* add(const SFC_Key& keyi);
    
    void remove(const SFC_Key& keyi);
    
public:
    //plain getters and setters
    int get_no_of_buckets() const {return NBUCKETS;}
    double* get_Xrange() {return Xrange;}
    double* get_Yrange() {return Yrange;}
    double* get_doublekeyrange() {return doublekeyrange;}
    double get_invdxrange() const {return invdxrange;}
    double get_invdyrange() const {return invdyrange;}
    int get_no_of_entries() const {return ENTRIES;}

    void print0() const {}
    
    void print_out(int) const {}
    
    tisize_t size() const{return elenode_.size();}
protected:
    void flushTable();//actually delete, removed nodes and rearrange added (sort according to keys)
    //!temporary arrays used in flush()
    vector<ti_ndx_t> ndx_map;
    vector<ti_ndx_t> ndx_map_old;
    vector<SFC_Key> key_map;
};
template class HashTable<Node>;
template class HashTable<Element>;

////////////////////////////////////////////////////////////////////////////////
//
constexpr tisize_t elem_reserved_size=1000000;
constexpr tisize_t node_reserved_size=elem_reserved_size*9;

class NodeHashTable: public HashTable<Node>
{
public:
    NodeHashTable(double *doublekeyrangein, int size, double XR[], double YR[]);
    ~NodeHashTable();
    
    Node* createAddNode(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr);
    Node* createAddNode(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr);
    Node* createAddNode(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada);
    Node* createAddNode(FILE* fp, MatProps* matprops_ptr);
    
    void removeNode(Node* node);

    Node& node(const ti_ndx_t ndx){return elenode_[ndx];}
    void flushNodeTable();
private:
     Node* addNode(const SFC_Key& keyi);
public: 
    //! used in delete_unused_nodes_and_elements() function 
    tivector<int> id_;

    //! the number of associated elements, it is used in extraneous node 
    //deletion and debugging function AssertMeshErrorFree()
    tivector<int> num_assoc_elem_;

    //! says what type of node this is see the comments of Node::get_info()
    tivector<int> info_;

    //! this is legacy afeapi and is not important though it would involve quite a bit of work to remove because it occurs frequently in Titan
    tivector<int> order_;

    //! the global x and y coordinates of the node
    tivector<double> coord_[DIMENSION];

    //! this is the node key, a key is a single number that is 2 unsigned variables long and is used to access the pointer to a Node or Element through the HashTable
    //SFC_Key key_;

    //! this elevation should currently be the GIS elevation at the finest "scale"
    tivector<double> elevation_;

    //! these are the so called "regular fluxes" that is the ones that are used to update the elements, assume that element normal is parallel to either the x or y axis, Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default
    tivector<double> flux_[MAX_NUM_STATE_VARS];

    //! the "refinement flux" is necessary when using the stopping criteria to reset the "regular" fluxes to what they would be if velocity was zero in the cell(s) involved.  The refinement flux is what the flux would have been if it had not been reset, they are needed since refinement is based on fluxes (and also pileheight gradient but that's not relevant here) Keith is the one who introduced a distinction between regular and refinement fluxes for use with the stopping criteria, this distinction is disabled by default.
    tivector<double> refinementflux_[MAX_NUM_STATE_VARS];

    //! node number for connection data -- varies with adaptation
    tivector<int> connection_id_;
};
extern NodeHashTable *nodeHashTable;

////////////////////////////////////////////////////////////////////////////////
//! Hashtables for Elements
class ElementsHashTable: public HashTable<Element>
{
    
    //friend int hash(unsigned* keyi);
    friend class Element;
protected:
    NodeHashTable* NodeTable;

    vector<SFC_Key> ukeyElements;
    vector<Element*> elements;

    int NlocalElements;
    vector<SFC_Key> ukeyLocalElements;
    vector<Element*> localElements;
public:
    ElementsHashTable(double *doublekeyrangein, int size, double XR[], double YR[], NodeHashTable* nodeTable);
    ~ElementsHashTable();

    Element** getElements()
    {
        return &(elements[0]);
    }
    
    void updateElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckElementsPointers(const char *prefix);

    int getNumberOfLocalElements()
    {
        return NlocalElements;
    }
    Element** getLocalElementsValues()
    {
        return &(localElements[0]);
    }
    //only works for elements
    void updateLocalElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckLocalElementsPointers(const char *prefix);

    //! set neighboring elements and nodes pointers in elements
    void updatePointersToNeighbours();

    //! check that neighboring elements and nodes pointers are correct
    int checkPointersToNeighbours(const char *prefix);
        
    //! default Element generator constructor, does nothing except set stoppedflags=2, this should never be used
    Element* generateAddElement(const SFC_Key& keyi);

    //! constructor that creates an original element when funky is read in
    Element* generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int mat, int *elm_loc_in,
                             double pile_height, int myid, const SFC_Key& opposite_brother);

    //! constructor that creates a son element from its father during refinement
    Element* generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen, int elm_loc_in[],
                int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in, ElementsHashTable *El_Table,
                NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather, double Awetfather,
                double *drypoint_in);
    //! constructor that creates a father element from its four sons during unrefinement
    Element* generateAddElement(Element *sons[], NodeHashTable *NodeTable, ElementsHashTable *El_Table, MatProps *matprops_ptr);

    //! constructor that creates/restores a saved element during restart
    Element* generateAddElement(FILE* fp, NodeHashTable* NodeTable, MatProps* matprops_ptr, int myid);
    
    //!should not be used
    void removeElement(Element* elm);
    
    void flushElemTable();
    

    Element& elem(const ti_ndx_t ndx){return elenode_[ndx];}
private:
     Element* addElement(const SFC_Key& keyi);
public:
    //here goes element content storage probably should be separate class at the end
    
    //! Element type
    ElementType elementType_;
    
    //! myprocess is id of the process(or) that owns this element
    tivector<int> myprocess_;

    //! generation is how many times this element has been refined, currently -8<=generation<=3, a negative generation number means it has been unrefined beyond the orignal coarse mesh, a positive generation number means it has been refined (is smaller than the original element size)
    tivector<int> generation_;

    //! opposite_brother_flag indicate if we have the correct key for the non-neighbor brother (0:= don't have info, 1:= have info)
    tivector<int> opposite_brother_flag_;

    //! the material flag indicates which material should be used to set this element's bed friction, this is for when a GIS material map, specifying different materials in different spatial regions of the map, the GIS material map is a non standard grass map format that Laercio Namikawa developed, it's stored in the "cats" folder under a grass mapset directory
    tivector<int> material_; /*! ! ! THE MAT. FLAG ! ! !*/

    //! this is the load-balancing weight
    tivector<double> lb_weight_;

    //! this is the key for load-balancing, if there is no constrained node, it is the element key, otherwise it is a construct of the element "bunch", keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    tivector<SFC_Key> lb_key_;

    //! this is the element key, which has the same value as the key of the element's bubble node, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    //tivector<SFC_Key> key_;

    //! this array holds the first 8 (0->7) of this element's nodes' keys, the n9th (8 out of 0->8) node is the bubble node it's key is not stored separately since it has the same key as the element, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    tivector<SFC_Key> node_key_[8];

    //!same as node_key but pointers, can be out-dated
    tivector<Node*> node_keyPtr[8];

    //! this array holds the keys of this element's 8 neighbors (2 neigbors to a side if the neighbor is more refined than this element, otherwise the two neighbor keys for that side are identical in value), having 8 neighbors is an outcome of the 1 irregularity refinement rule, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    tivector<SFC_Key> neighbors_[8];

    //!same as neighbor but pointers, can be out-dated
    tivector<Element*> neighborPtr[8];

    //! the key of the father it is assigned in the refine() and unrefine_elements() functions
    tivector<SFC_Key> father_;

    //! this array holds the keys of this element's 4 sons, it is only used temporarily in the refinement process before the father (this element) is deleted, there's was an old comment "garantee ccw" associated with this variable, I don't know what it means, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    tivector<SFC_Key> son_[4];

    //! this array holds the process(or) id of this element's 8 neighbors, there can be 8 neighbors because of the 1 irregularity rule.  neigh_proc[4:7] != -2 only if it has 2 neighbors on that side, a value of -1 for neigh_proc means that this edge is a boundary of the computational domain. 
    tivector<int> neigh_proc_[8];

    //! this is legacy afeapi, all finite volume "elements"/cells are piece wise constant, but I believe this is actually used in the DG (Discontinuous Galerkin) version of titan
    tivector<int> order_[5];

    //! neigh_gen is an array that holds the "generation" (how refined it is) of this element's 8 neighbors, there can-be/are 2 neighbors to a side because of the 1 irregularity rule
    tivector<int> neigh_gen_[8];

    //! pointer to the boundary condition class, if this element is not a boundary element the pointer holds the NULL value
    tivector<BC*> bcptr_;

    //! the number of degrees of freedom, since Titan is a finite difference/volume code, ndof is afeapi legacy, but the DG (Discontinuous Galerkin) version of Titan actually uses this
    tivector<int> ndof_;

    //! this is legacy afeapi, it is not used, but do not remove it, it could cause problems if you do
    tivector<int> no_of_eqns_;

    //! this holds the "error" in the element's solution, which is useful in determining refinement, this may actually be afeapi legacy
    tivector<double> el_error_[EQUATIONS];

    //! this holds the element solution, I believe this is legacy afeapi
    tivector<double> el_solution_[EQUATIONS];

    //! refined is a flag that usually has the value 0, but will be 1 if the element has been refined this iteration (used to enforce the 1 irregularity rule), or have the value "GHOST" if it is a ghost cell, refined and ghost cells are not updated, see constant.h for the value of GHOST
    tivector<int> refined_;

    //! The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
    tivector<int> adapted_;

    //! which_son holds the value of which son this element is, which of the 4 squares that makes up the father elements square.
    tivector<int> which_son_;

    //! the new_old flag is used in mesh adaptation and repartitioning
    tivector<int> new_old_;

    //! this array holds the keys of this element's 4 brothers (an element is considered to be it's own brother), this information is used during mesh unrefinement (combining the 4 brothers to make their father), keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
    tivector<SFC_Key> brothers_[4];

    //! coord holds the coordinates of the elements cell center, these are the same as the coordinates of the element's bubble node's
    tivector<double> coord_[DIMENSION];

    //! elm_loc is used in unrefining beyond the original coarse mesh
    tivector<int> elm_loc_[DIMENSION];


    /* variables for hyperbolic geoflow problem */

    //! state_vars is an array that holds the current state variables: h, hVx, and hVy 
    tivector<double> state_vars_[MAX_NUM_STATE_VARS];
    //! these are the values of the state variables from before the predictor step
    tivector<double> prev_state_vars_[MAX_NUM_STATE_VARS];

    //! these are the spatial (x and y) derivatives of the state variables: (dh/dx, dhVx/dx, dhVy/dx, dh/dy, dhVx/dy, dhVy/dy)
    tivector<double> d_state_vars_[MAX_NUM_STATE_VARS * DIMENSION];

    //! the short speed is the speed computed as: shortspeed=|v|=|dhv/dh|=|v*dh/dh+h*dv/dh|=|v+h*dv/dh| which goes to |v| in the limit of h->0, this is a more accurate way to compute speed when the pile in this cell is short, hence the name "shortspeed" but it is not accurate when the pile is tall, that is when h*dv/dh is large, this is the value from the previous iteration (so there is lagging when using the shortspeed, but this should still be much more accurate than hV/h when h->0. Keith implemented this in late summer 2006, 
    tivector<double> shortspeed_;

    //! length of the element in the global x and y directions: dx and dy 
    tivector<double> dx_[DIMENSION];

    //! for structured grid, tells which side is the positive x direction
    tivector<int> positive_x_side_;

    //! maximum x and y direction wavespeeds for this element, wavespeeds are eigenvalues of the flux jacobians
    tivector<double> eigenvxymax_[DIMENSION];

    //! k active/passive in the x and y directions, k active/passive is part of the coulomb friction model for Granular Flows
    tivector<double> kactxy_[DIMENSION];

    //! terrain elevation at this elements center/bubble node 
    tivector<double> elevation_;

    //! terrain slope in the global x and y directions
    tivector<double> zeta_[DIMENSION];

    //! Curvature itself is the inverse of radius of curvature.  The exact value of curvature is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation. 
    tivector<double> curvature_[DIMENSION];

    //! the gravity vector in local x,y,z coordinates (z is normal to the terrain surface, the projections of the x and y local directions onto a horizontal plane are aligned with the global x and y directions)
    tivector<double> gravity_[3];

    //! the spatial (x and y) derivatives of the local z component of the gravity vector
    tivector<double> d_gravity_[DIMENSION];

    //! part of the new stopping criteria under development, has a 0 if flow is not stopped, has the value 1 if it should not be sliding but should be slumping, has the value 2 if it should neither be sliding or slumping (it should be completely stopped), I (Keith) am rather confident in the criteria used to set this the problem is determining what to do about it after you know the flow SHOULD be stopped
    tivector<int> stoppedflags_;

    //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus the effective bed friction angle holds either the value of the actual bed friction angle if it should not be stopped or the value of the internal friction angle if it should not be stopped
    tivector<double> effect_bedfrict_;

    //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus effect_tanbedfrict holds the value of the effective bed friction angle
    tivector<double> effect_tanbedfrict_;

    //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, if the effective bed friction angle equals the internal friction angle effect_kactxy takes on the value 1, k active/passive comes from using a Coulomb friction model for granular flows
    tivector<double> effect_kactxy_[2];

    //! extrusion flux rate for this timestep for this element, used when having material flow out of the ground, a volume per unit area influx rate source term
    tivector<double> Influx_[MAX_NUM_STATE_VARS];

    //! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
    tivector<int> ithelem_;

    //! the node number {0,1,..,7} of this element's "most wet node", this dictates both the orientation of the "dryline" and which side of it is wet.  The "dryline" is the line that divides a partially wetted element into a dry part and a wet part, for the sake of simplicity only 4 orientations are allowed, horizontal, vertical, parallel to either diagonal of the element.  If the iwetnode is an edge node of this element then the dryline is parallel to the edge the element is on, if the iwetnode is a corner node of this element then dryline is parallel to the diagonal of the element that the iwetnode is not on.  Which side of the dryline is wet is the same side in which the iwetnode resides (and is determined each timestep based soley on which of the elements neighbors currently have pile height greater than GEOFLOW_TINY)... as such iwetnode can be thought of as the MOST WET NODE  of the element.  Having iwetnode==8 indicates that the element is uniformly wet if this element's pile height is greater than GEOFLOW_TINY or is uniformly dry if pileheight is less than or equal to GEOFLOW_TINY.
    tivector<int> iwetnode_;

    //! Awet is the ratio of this element's wet area to total area (always between 0 and 1 inclusive) when taken together with iwetnode, this uniquely determines the exact placement of the "dryline" within the current element.  Awet is initially set by source placement to be either 0 (no material) or 1 (material) and is updated by the corrector part of the predictor-corrector method, the new value is determined by where the dry line has been convected to over this timestep. Keith wrote this May 2007.
    tivector<double> Awet_;

    //! center point of the "dryline", x and y coordinates value ranges between -0.5 and 0.5 with 0 being the center of the element, since the wet/dry interface is taken to be a non-deforming non rotating (within the timestep) "dryline" convecting a single point on the dryline (called the drypoint) is sufficient to determine the new placement of the dryline which allows us to update Awet... Keith wrote this May 2007.
    tivector<double> drypoint_[2];

    //! when an element edge is partially wet and partially dry... Swet is the fraction of a cell edge that is partially wet, because it can only be horizontal, vertical, or parallel to either diagonal, all of one element's partially wet sides are have the same fraction of wetness.  The state variables (used to compute the physical fluxes) at the element/cell edge are adjusted to be the weighted by wetness average over an element/cell edge.  As such physical fluxes through completely dry edges of partially wet elements/cells are zeroed, while physical fluxes through completely wet edges are left unchanged.  Because of the definition as "wetness weighted average" physical fluxes through a partially wet edge shared with a neighbor of the same generation is also left left unchanged but, when a partially wet edge is shared with two more refined neighbors the total mass and momentum at the edge is split between the two neighbors in proportion to how much of their boundary shared with this element is wet.  This "scaling" of the physical fluxes is the "adjustment of fluxes in partially wetted cells" facet of our multifaceted thin-layer problem mitigation approach.  And it has been shown to significantly reduce the area covered by a thin layer of material.  Keith wrote this May 2007.
    tivector<double> Swet_;
};

extern ElementsHashTable *elementsHashTable;

#endif
