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

//! this is the normal grid adaptive refinement function it also refreshes the flux sources
class HAdapt
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
public:
    HAdapt(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable):TempList(_ElemTable, 384)
    {
        ElemTable=_ElemTable;
        NodeTable=_NodeTable;
    }
    void adapt(int h_count, double target, MatProps* matprops_ptr,
             FluxProps *fluxprops_ptr, TimeProps* timeprops_ptr, int num_buffer_layer);
    
    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);
    void depchk2(ti_ndx_t ndx, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);
    void refine2(Element* EmTemp, MatProps* matprops_ptr);
    
private:
    ElemPtrList TempList;
};

#endif	/* HADAPT_H */

