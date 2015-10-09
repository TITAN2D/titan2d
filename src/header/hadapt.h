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

//!interfce for seed refinements finders
class SeedRefinementsFinder
{
public:
	virtual ~SeedRefinementsFinder(){};
	virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement)=0;
};

//!PrimaryRefinementsFinder
class PrimaryRefinementsFinder:public SeedRefinementsFinder
{
public:
	PrimaryRefinementsFinder(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable);
	virtual ~PrimaryRefinementsFinder(){}
	virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
public:
	double geo_target;
private:
	ElementsHashTable* ElemTable;
	NodeHashTable* NodeTable;
	tivector<Element> &elements;
	tivector<ContentStatus> &status;
	tivector<int> &adapted;
	tivector<int> &generation;
	tivector<double> &el_error;
};

//!BuferFirstLayerRefinementsFinder
class BuferFirstLayerRefinementsFinder:public SeedRefinementsFinder
{
public:
    BuferFirstLayerRefinementsFinder(ElementsHashTable* _ElemTable);
    virtual ~BuferFirstLayerRefinementsFinder(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    ElementsHashTable* ElemTable;
    tivector<Element> &elements;
    tivector<ContentStatus> &status;
};


//!BuferNextLayerRefinementsFinder
class BuferNextLayerRefinementsFinder:public SeedRefinementsFinder
{
public:
    BuferNextLayerRefinementsFinder(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable);
    virtual ~BuferNextLayerRefinementsFinder(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    tivector<Element> &elements;
    tivector<ContentStatus> &status;
};


//! this is the normal grid adaptive refinement function it also refreshes the flux sources
class HAdapt
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    MatProps* matprops_ptr;
    TimeProps* timeprops_ptr;
    int num_buffer_layer;
public:
    HAdapt(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable,TimeProps* _timeprops, MatProps* _matprops, const int _num_buffer_layer);
    void adapt(int h_count, double target);
    
    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);


    

private:
    void refine2(SeedRefinementsFinder &seedRefinementsFinder);

    void findPrimaryRefinements(vector<ti_ndx_t> &primaryRefinement, const double geo_target);
    void findTriggeredRefinements(const vector<ti_ndx_t> &primaryRefinement, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);
    void findBuferFirstLayerRefinements(vector<ti_ndx_t> &primaryRefinement);
    void findBuferNextLayerRefinements(vector<ti_ndx_t> &primaryRefinement);

    void refineElements(const vector<ti_ndx_t> &allRefinement);

    void create_new_node2(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], const int info, int& RefNe, const int boundary, const int order);
    void refinedNeighboursUpdate(const vector<ti_ndx_t> &allRefinement);

private:
    PrimaryRefinementsFinder primaryRefinementsFinder;
    BuferFirstLayerRefinementsFinder buferFirstLayerRefinementsFinder;
    BuferNextLayerRefinementsFinder buferNextLayerRefinementsFinder;

private:
    //temporary arrays used during refinement
    vector<int> set_for_refinement;
    vector<ti_ndx_t> seedRefinement;
    vector<ti_ndx_t> allRefinement;

private:
    int myid;
    int numprocs;
    ElemPtrList TempList;
    vector<ti_ndx_t> tempList;
};

#endif	/* HADAPT_H */

