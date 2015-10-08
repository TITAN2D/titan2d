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

//! this is the normal grid adaptive refinement function it also refreshes the flux sources
class HAdapt
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    MatProps* matprops_ptr;
    TimeProps* timeprops_ptr;
public:
    HAdapt(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable,TimeProps* _timeprops, MatProps* _matprops):TempList(_ElemTable, 384)
    {
        ElemTable=_ElemTable;
        NodeTable=_NodeTable;
        matprops_ptr=_matprops;
        timeprops_ptr=_timeprops;
    }
    void adapt(int h_count, double target,
             FluxProps *fluxprops_ptr, int num_buffer_layer);
    
    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);


    

private:
    void findPrimaryRefinements(vector<ti_ndx_t> &primaryRefinement, const double geo_target);
    void findTriggeredRefinements(const vector<ti_ndx_t> &primaryRefinement, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);

    void refineElements(const vector<ti_ndx_t> &allRefinement);

    void create_new_node2(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], const int info, int& RefNe, const int boundary, const int order);
    void refinedNeighboursUpdate(const vector<ti_ndx_t> &allRefinement);

private:
    int myid;
    int numprocs;
    ElemPtrList TempList;
    vector<ti_ndx_t> tempList;
};

#endif	/* HADAPT_H */

