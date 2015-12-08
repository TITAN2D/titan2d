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

#ifndef HADAPT_H
#define	HADAPT_H

#include <vector>
#include <array>

//!interfce for seed refinements finders
class SeedRefinementsFinder
{
public:
    SeedRefinementsFinder(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
    ~SeedRefinementsFinder(){};
	virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement)=0;
protected:
	ElementsHashTable* ElemTable;
	NodeHashTable* NodeTable;
	ElementsProperties* ElemProp;

    tivector<Element> &elements;
    tivector<ContentStatus> &status;
    tivector<int> &adapted;
    tivector<int> &generation;
};

//!PrimaryRefinementsFinder
class PrimaryRefinementsFinder:public SeedRefinementsFinder
{
public:
	PrimaryRefinementsFinder(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
	virtual ~PrimaryRefinementsFinder(){}
	virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
public:
	double geo_target;
private:
	tivector<double> &el_error;

	vector< vector<ti_ndx_t> > loc_SeedRefinement;
};

//!BuferFirstLayerRefinementsFinder
class BuferFirstLayerRefinementsFinder:public SeedRefinementsFinder
{
public:
    BuferFirstLayerRefinementsFinder(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
    virtual ~BuferFirstLayerRefinementsFinder(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    vector< vector<ti_ndx_t> > loc_SeedRefinement;
};


//!BuferNextLayerRefinementsFinder
class BuferNextLayerRefinementsFinder:public SeedRefinementsFinder
{
public:
    BuferNextLayerRefinementsFinder(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
    virtual ~BuferNextLayerRefinementsFinder(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    vector< vector<ti_ndx_t> > loc_SeedRefinement;
};


//! this is the normal grid adaptive refinement function it also refreshes the flux sources
class HAdapt
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    ElementsProperties* ElemProp;
    MatProps* matprops_ptr;
    TimeProps* timeprops_ptr;
    int num_buffer_layer;
public:
    HAdapt(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable, ElementsProperties* ElemProp,TimeProps* _timeprops, MatProps* _matprops, const int _num_buffer_layer);
    void adapt(int h_count, double target);
    
    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);


    

private:
    void refine2(SeedRefinementsFinder &seedRefinementsFinder);

    void findPrimaryRefinements(vector<ti_ndx_t> &primaryRefinement, const double geo_target);
    void findTriggeredRefinements(const vector<ti_ndx_t> &primaryRefinement, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);
    void findBuferFirstLayerRefinements(vector<ti_ndx_t> &primaryRefinement);
    void findBuferNextLayerRefinements(vector<ti_ndx_t> &primaryRefinement);

    void refineElements(const vector<ti_ndx_t> &allRefinement);

    void check_create_new_node(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info, int& RefNe, const int boundary);
    void check_create_new_node2(const int iElm, const int iNode, const int info, int& RefNe, const int boundary);

    void create_new_node3(const int which, const ti_ndx_t Node1, const ti_ndx_t Node2,
            SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info);

    void refinedNeighboursUpdate(const vector<ti_ndx_t> &allRefinement);
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const array<double,2> Node1, const array<double,2> Node2);
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const ti_ndx_t Node1, const ti_ndx_t Node2);
private:
    PrimaryRefinementsFinder primaryRefinementsFinder;
    BuferFirstLayerRefinementsFinder buferFirstLayerRefinementsFinder;
    BuferNextLayerRefinementsFinder buferNextLayerRefinementsFinder;

private:
    //temporary arrays used during refinement
    vector<int> set_for_refinement;
    vector<ti_ndx_t> seedRefinement;
    vector< vector<ti_ndx_t> > loc_SeedRefinement;
    vector<ti_ndx_t> allRefinement;

    //temporary arrays used during refinement (refineElements)
    vector<array<ti_ndx_t,16> > new_node_ndx;
    vector<array<SFC_Key,16> > new_node_key;
    vector<array<array<double,2>, 16> > new_node_coord;
    vector<array<bool, 16> > new_node_isnew;
    vector<array<ti_ndx_t,9> > node_ndx_ref;
    //!map for indexes of elements for refinement: refining_elem_map[allRefinement[i]]=i or -1 if not refined
    vector<int> refining_elem_map;

    vector<vector<int> > create_node_ielm;
    vector<vector<int> > create_node_iwhich;
    //vector<SFC_Key> new_key;
    //vector<array<double,2> > new_coord;

private:
    int myid;
    int numprocs;
    ElemPtrList TempList;
    vector<ti_ndx_t> tempList;

private:
    int which_neighbor(ti_ndx_t ndx,ti_ndx_t neigh_elm_ndx)
    {
        int i = 0;
        int which = -1;
        while (i < 4 && which == -1)
        {
            if(ElemTable->neighbors_[i][neigh_elm_ndx]==ElemTable->key_[ndx])
                which = i;
            i++;
        }
        ASSERT2(which != -1);
        return which;
    }

};

//! class for unrefinement
class HAdaptUnrefine
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    MatProps* matprops_ptr;
    TimeProps* timeprops_ptr;
public:
    HAdaptUnrefine(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable,TimeProps* _timeprops, MatProps* _matprops);


    //! this function loops through all the elements on this processor and (by calling other functions) checks which elements satisfy criteria for being okay to unrefine, if they can be it unrefines them.
    void unrefine(const double target);
private:
    void unrefine_neigh_update();
    void unrefine_interp_neigh_update();
    void delete_oldsons();

private:
    //temporary arrays used during refinement
    vector<ti_ndx_t> NewFatherList;
    vector<ti_ndx_t> OtherProcUpdate;
private:
    int myid;
    int numprocs;
};




#endif	/* HADAPT_H */

