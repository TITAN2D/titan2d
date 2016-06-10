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
 * $Id: refine2.C 127 2007-06-07 19:48:25Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/ticore/omp_mpi.hpp"

#include "../header/ticore.hpp"

#include "../header/hpfem.h"
#include "../header/hadapt.h"

extern void fhsfc2d_(double, unsigned, unsigned);
extern void hsfc2d(unsigned*, unsigned*, unsigned*);
extern void create_new_node(int, int, int, NodeHashTable*, Node*[], SFC_Key[], int, int*, int, int, MatProps*);


void HAdapt::check_create_new_node(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info, int& RefNe, const int boundary)
{
    double NewNodeCoord[2];
    double norm_coord[2];
    unsigned u_norm_coord[2];
    unsigned nkey = 2;
    SFC_Key key;
    unsigned oldkey[KEYLENGTH];
    
    ti_ndx_t ndx;
    static double XRange[2];
    static double YRange[2];
    int i;
    
    for(i = 0; i < 2; i++)
    {
        XRange[i] = NodeTable->get_Xrange()[i];
        YRange[i] = NodeTable->get_Yrange()[i];
    }
    
    for(i = 0; i < 2; i++)
        NewNodeCoord[i] = (NodeTable->coord_[i][ndxNodeTemp[Node1]] + NodeTable->coord_[i][ndxNodeTemp[Node2]]) * .5;
    
    norm_coord[0] = (NewNodeCoord[0] - XRange[0]) / (XRange[1] - XRange[0]);
    norm_coord[1] = (NewNodeCoord[1] - YRange[0]) / (YRange[1] - YRange[0]);
    
    fhsfc2d_(norm_coord, &nkey, oldkey);
    
    SET_NEWKEY(key,oldkey);
    NewNodeKey[which]=key;
    
    ndx = NodeTable->lookup_ndx(key);
    
    if(ti_ndx_negative(ndx))
    {
        ndx = NodeTable->createAddNode_ndx(key, NewNodeCoord, info, matprops_ptr);
    }
    else if(NodeTable->coord_[0][ndx] != NewNodeCoord[0] || NodeTable->coord_[1][ndx] != NewNodeCoord[1])
    {
        short same_key = 0;
        assert(same_key);
    }
    else
    {
        RefNe = 1;
    }
    NewNodeNdx[which]=ndx;
    
    if(RefNe || boundary)
        NodeTable->info_[ndx]=SIDE;
    
    return;
}
void HAdapt::check_create_new_node2(const int iElm, const int iNode, const int info, int& RefNe, const int boundary)
{
    ti_ndx_t ndx = new_node_ndx[iElm][iNode];

    if(new_node_isnew[iElm][iNode])
    {
        NodeTable->elenode_[ndx].init(new_node_key[iElm][iNode], &(new_node_coord[iElm][iNode][0]), info, -3, matprops_ptr);
        new_node_isnew[iElm][iNode]=false;
    }
    else if(NodeTable->coord_[0][ndx] != new_node_coord[iElm][iNode][0] || NodeTable->coord_[1][ndx] != new_node_coord[iElm][iNode][1])
    {
        assert(0);
    }
    else
    {
        RefNe = 1;
    }

    if(RefNe || boundary)
        NodeTable->info_[ndx]=SIDE;

    return;
}
void HAdapt::create_new_node3(const int which, const ti_ndx_t Node1, const ti_ndx_t Node2,
                     SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info)
{
    double NewNodeCoord[2];
    double norm_coord[2];
    unsigned u_norm_coord[2];
    unsigned nkey = 2;
    SFC_Key key;
    unsigned oldkey[KEYLENGTH];

    ti_ndx_t ndx;

    for(int i = 0; i < 2; i++)
        NewNodeCoord[i] = (NodeTable->coord_[i][Node1] + NodeTable->coord_[i][Node2]) * .5;

    norm_coord[0] = (NewNodeCoord[0] - NodeTable->Xrange[0]) / (NodeTable->Xrange[1] - NodeTable->Xrange[0]);
    norm_coord[1] = (NewNodeCoord[1] - NodeTable->Yrange[0]) / (NodeTable->Yrange[1] - NodeTable->Yrange[0]);

    fhsfc2d_(norm_coord, &nkey, oldkey);

    SET_NEWKEY(key,oldkey);
    NewNodeKey[which]=key;

    ndx=NodeTable->lookup_ndx(key);
    assert(ti_ndx_not_negative(ndx));
    NodeTable->elenode_[ndx].init(key, NewNodeCoord, info, -3, matprops_ptr);

    NewNodeNdx[which]=ndx;
    return;
}

// (#) - NewNodeKey - new node numbering
//  #  - old node numbering
//  E# - old neighbouring elements numbering
//
//              Side 2
//           E6        E2
//       3---(14)--6---(15)--2
//       |         |         |
//       |         |         |
// S E3 (11) (3)  (12) (2)  (13) E5 S
// i     |         |         |      i
// d     |         |         |      d
// e     7---(9)--8/E--(10)--5      e
//       |         |         |
// 3     |         |         |      1
//   E7 (6)  (0)  (7)  (1)  (8)  E1
//       |         |         |
//       |         |         |
//       0---(4)---4---(5)---1
//            E0       E4
//              Side 0

//max & min coordinates need to be passed later now only for L-shape!!!
//only 4 one step because of the info FLAG!!!
//if the new node is on INTERFACE flag will be -1

// (0) to (3) also which_sons value

void HAdapt::refineElements(const vector<ti_ndx_t> &allRefinement)
{
    PROFILING3_DEFINE(pt_start);
    PROFILING3_START(pt_start);

    //convinience references
    const int numElemToRefine=allRefinement.size();
    const ti_ndx_t *ElemToRefine=&(allRefinement[0]);


    refining_elem_map.assign(ElemTable->size(),-1);
    for(int iElm=0;iElm<numElemToRefine;++iElm)
        refining_elem_map[ElemToRefine[iElm]]=iElm;


	tivector<int> &adapted=ElemTable->adapted_;
	tivector<int> &refined=ElemTable->refined_;

	//find position of 16 nodes

	//resize and init new_node_ndx

	node_ndx_ref.resize(numElemToRefine);
    for(int iElm=0;iElm<numElemToRefine;++iElm)
        for(int k=0;k<9;++k)
            node_ndx_ref[iElm][k]=ti_ndx_doesnt_exist;

    new_node_ndx.resize(numElemToRefine);
	for(int iElm=0;iElm<numElemToRefine;++iElm)
	    for(int k=0;k<16;++k)
	        new_node_ndx[iElm][k]=ti_ndx_doesnt_exist;
    new_node_key.resize(numElemToRefine);
    for(int iElm=0;iElm<numElemToRefine;++iElm)
        for(int k=0;k<16;++k)
            new_node_key[iElm][k]=sfc_key_null;
    new_node_isnew.resize(numElemToRefine);
    for(int iElm=0;iElm<numElemToRefine;++iElm)
        for(int k=0;k<16;++k)
            new_node_isnew[iElm][k]=false;
    new_node_coord.resize(numElemToRefine);
#ifdef DEB3
    for(int iElm=0;iElm<numElemToRefine;++iElm)
            for(int k=0;k<16;++k)
                for(int i=0;i<2;++i)
                    new_node_coord[iElm][k][i]=i;
#endif
    for(int ithread=0;ithread<threads_number;++ithread)
    {
        create_node_ielm[ithread].resize(0);
        create_node_iwhich[ithread].resize(0);
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_init,pt_start);

	//find position of corners, sides and bubbles
    #pragma omp parallel for schedule(static)
	for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];

		for(int i = 0; i < 8; i++) //-- corners and sides
		{
		    node_ndx_ref[iElm][i] = ElemTable->node_key_ndx_[i][ndx];
			ASSERT2(ElemTable->node_key_ndx_[i][ndx]==NodeTable->lookup_ndx(ElemTable->node_key_[i][ndx]));
			ASSERT2(ti_ndx_not_negative(node_ndx_ref[iElm][i]));
		}

        //-- bubble
		node_ndx_ref[iElm][8] = ElemTable->node_bubble_ndx_[ndx];
		ASSERT2(ElemTable->node_bubble_ndx_[ndx]==NodeTable->lookup_ndx(ElemTable->key_[ndx]));
        ASSERT2(ti_ndx_not_negative(node_ndx_ref[iElm][8]));
    }
	PROFILING3_STOPADD_RESTART(HAdapt_refineElements_find_corners_sides_bubbles,pt_start);
	//SIDE 0
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iElm=0;iElm<numElemToRefine;++iElm)
        {
            ti_ndx_t ndx=ElemToRefine[iElm];
            if(ElemTable->neigh_proc_[0][ndx] == -1 || ElemTable->neigh_gen_[0][ndx] <= ElemTable->generation_[ndx])
            {
                //i.e. boundary of the computational domain or neighbor generation same or smaller then this one
                rE__find_new_side_node_or_set_for_alloc(
                        iElm, 4/*which*/, ndx, 0 /*neigh*/,
                        node_ndx_ref[iElm][0], node_ndx_ref[iElm][4],ithread);
                rE__find_new_side_node_or_set_for_alloc(
                        iElm, 5/*which*/, ndx, 4 /*neigh*/,
                        node_ndx_ref[iElm][1], node_ndx_ref[iElm][4],ithread);
            }
            else
            {
                //i.e. not boundary of the computational domain and neighbor generation higher then this one
                rE__find_new_side_node(iElm,4/*which*/,ndx,0/*neigh*/,6/*neigh_which*/);
                rE__find_new_side_node(iElm,5/*which*/,ndx,4/*neigh*/,6/*neigh_which*/);
            }
        }
	}
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side0_find_nodes,pt_start);

    //rE__alloc_new_nodes();
	PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side0_add_new_nodes,pt_start);
    //SIDE 1
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iElm=0;iElm<numElemToRefine;++iElm)
        {
            ti_ndx_t ndx=ElemToRefine[iElm];
            if(ElemTable->neigh_proc_[1][ndx] == -1 || ElemTable->neigh_gen_[1][ndx] <= ElemTable->generation_[ndx])
            {
                //i.e. boundary of the computational domain or neighbor generation same or smaller then this one
                rE__find_new_side_node_or_set_for_alloc(
                       iElm, 8/*which*/, ndx, 1 /*neigh*/,
                       node_ndx_ref[iElm][1], node_ndx_ref[iElm][5],ithread);
                rE__find_new_side_node_or_set_for_alloc(
                        iElm, 13/*which*/, ndx, 5 /*neigh*/,
                        node_ndx_ref[iElm][2], node_ndx_ref[iElm][5],ithread);
            }
            else
            {
                //i.e. not boundary of the computational domain and neighbor generation higher then this one
                rE__find_new_side_node(iElm,8/*which*/,ndx,1/*neigh*/,7/*neigh_which*/);
                rE__find_new_side_node(iElm,13/*which*/,ndx,5/*neigh*/,7/*neigh_which*/);
            }
        }
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side1_find_nodes,pt_start);

    rE__alloc_new_nodes();
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side1_add_new_nodes,pt_start);
	//SIDE 2
    #pragma omp parallel
	{
        int ithread=omp_get_thread_num();

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iElm=0;iElm<numElemToRefine;++iElm)
        {
            ti_ndx_t ndx=ElemToRefine[iElm];
            if(ElemTable->neigh_proc_[2][ndx] == -1 || ElemTable->neigh_gen_[2][ndx] <= ElemTable->generation_[ndx])
            {
                //i.e. boundary of the computational domain or neighbor generation same or smaller then this one
                rE__find_new_side_node_or_set_for_alloc_refcheck(
                        iElm, 14/*which*/, ndx, 6 /*neigh*/, 4/*neigh_which*/,
                        node_ndx_ref[iElm][3], node_ndx_ref[iElm][6],ithread);
                rE__find_new_side_node_or_set_for_alloc_refcheck(
                        iElm, 15/*which*/, ndx, 2 /*neigh*/, 5/*neigh_which*/,
                        node_ndx_ref[iElm][2], node_ndx_ref[iElm][6],ithread);
            }
            else
            {
                //i.e. not boundary of the computational domain and neighbor generation higher then this one
                rE__find_new_side_node(iElm,14/*which*/,ndx,6/*neigh*/,4/*neigh_which*/);
                rE__find_new_side_node(iElm,15/*which*/,ndx,2/*neigh*/,4/*neigh_which*/);
            }
        }
	}
	PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side2_find_nodes,pt_start);

	//rE__alloc_new_nodes();
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side2_add_new_nodes,pt_start);

    //SIDE 3
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();

        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iElm=0;iElm<numElemToRefine;++iElm)
        {
            ti_ndx_t ndx=ElemToRefine[iElm];
            if(ElemTable->neigh_proc_[3][ndx] == -1 || ElemTable->neigh_gen_[3][ndx] <= ElemTable->generation_[ndx])
            {
                //i.e. boundary of the computational domain or neighbor generation same or smaller then this one
                rE__find_new_side_node_or_set_for_alloc_refcheck(
                        iElm, 6/*which*/, ndx, 3 /*neigh*/, 8/*neigh_which*/,
                        node_ndx_ref[iElm][0], node_ndx_ref[iElm][7],ithread);
                rE__find_new_side_node_or_set_for_alloc_refcheck(
                        iElm, 11/*which*/, ndx, 3 /*neigh*/, 13/*neigh_which*/,
                        node_ndx_ref[iElm][3], node_ndx_ref[iElm][7],ithread);
            }
            else
            {
                //i.e. not boundary of the computational domain and neighbor generation higher then this one
                rE__find_new_side_node(iElm,6/*which*/,ndx,7/*neigh*/,5/*neigh_which*/);
                rE__find_new_side_node(iElm,11/*which*/,ndx,3/*neigh*/,5/*neigh_which*/);
            }
        }
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side3_find_nodes,pt_start);

    //rE__alloc_new_nodes();
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side3_add_new_nodes,pt_start);

	//INTERNAL SIDES and NEW BUBBLES
    //first find keys and coords for new nodes
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();
        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iElm=0;iElm<numElemToRefine;++iElm)
        {
            //++++++++++++++++INTERNAL SIDE NODES 7, 12, 9, 10
            //---Seventh new node---
            calc_coord_and_key(new_node_key[iElm][7],new_node_coord[iElm][7], node_ndx_ref[iElm][4], node_ndx_ref[iElm][8]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(7);
            //---Twelwth new node---
            calc_coord_and_key(new_node_key[iElm][12],new_node_coord[iElm][12], node_ndx_ref[iElm][6], node_ndx_ref[iElm][8]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(12);
            //---Ninth new node---
            calc_coord_and_key(new_node_key[iElm][9],new_node_coord[iElm][9], node_ndx_ref[iElm][7], node_ndx_ref[iElm][8]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(9);
            //---Tenth new node---
            calc_coord_and_key(new_node_key[iElm][10],new_node_coord[iElm][10], node_ndx_ref[iElm][5], node_ndx_ref[iElm][8]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(10);
            //+++++++++++++++++++THE NEW BUBBLES 0, 1, 2, 3
            //---0th new node---
            calc_coord_and_key(new_node_key[iElm][0],new_node_coord[iElm][0], new_node_coord[iElm][6], new_node_coord[iElm][7]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(0);
            //---1st new node---
            calc_coord_and_key(new_node_key[iElm][1],new_node_coord[iElm][1], new_node_coord[iElm][7], new_node_coord[iElm][8]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(1);
            //---2nd new node---
            calc_coord_and_key(new_node_key[iElm][2],new_node_coord[iElm][2], new_node_coord[iElm][12], new_node_coord[iElm][13]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(2);
            //---3rd new node---
            calc_coord_and_key(new_node_key[iElm][3],new_node_coord[iElm][3], new_node_coord[iElm][11], new_node_coord[iElm][12]);
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(3);
        }
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_int_nodes_calc_keys,pt_start);
    //allocate and insert new nodes
    rE__alloc_new_nodes();
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_int_nodes_alloc,pt_start);

	//SIDE 0
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        //@todo replace which with macro or enum
        //@todo reuse indexes for nodes and elements
        int which;
        ti_ndx_t n1_ndx, n2_ndx, n3_ndx, n4_ndx;

        ti_ndx_t *ndxNodeTemp=&(node_ndx_ref[iElm][0]);
        SFC_Key *NewNodeKey=&(new_node_key[iElm][0]);
        ti_ndx_t *NewNodeNdx=&(new_node_ndx[iElm][0]);

        int i;
        ti_ndx_t neigh_elm_ndx;
        int RefinedNeigh = 0;
        int info;
        int other_proc = 0;
        int boundary;

		//SIDE 0
		if(ElemTable->neigh_proc_[0][ndx] == -1)
			boundary = 1;
		else
			boundary = 0;

		if(boundary == 1 || ElemTable->neigh_gen_[0][ndx] <= ElemTable->generation_[ndx])
		{
		    //i.e. boundary of the computational domain or neighbor generation same or smaller then this one
			RefinedNeigh = 0;
			info = S_S_CON;
			if(ElemTable->neigh_proc_[0][ndx] != myid)
			{
				other_proc = 1;
				info = -1;
			}
			else
			    other_proc = 0;

			which = 4;
			//---Fourth new node---
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);

			//---Fourth old node---
			if(RefinedNeigh || boundary)
				NodeTable->info_[ndxNodeTemp[4]]=CORNER;
			else if(other_proc)
				NodeTable->info_[ndxNodeTemp[4]]=-1;
			else
				NodeTable->info_[ndxNodeTemp[4]]=S_C_CON;

			//---Fifth new node---
			which = 5;
			//n2 = (Node*) NodeTable->lookup(EmTemp->node_key(1));

			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);
		}
		else
		{
		    //i.e. not boundary of the computational domain and neighbor generation higher then this one
		    //
			//Keith Added this if
			if((ElemTable->neigh_proc_[0][ndx] != myid) || ((ElemTable->neigh_proc_[4][ndx] != myid)
					&& (ElemTable->neigh_proc_[4][ndx] != -2)))
				other_proc = 1;
			else
				other_proc = 0;

			// fourth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[0][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[0][ndx]));

			n1_ndx = NewNodeNdx[4];
			ASSERT2(NewNodeNdx[4] == NodeTable->lookup_ndx(NewNodeKey[4]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;
			//fourth old node
			NodeTable->info_[ndxNodeTemp[4]]=CORNER;
			if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
				NodeTable->info_[ndxNodeTemp[4]]=-1;

			// fifth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[4][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[4][ndx]));

			n1_ndx = NewNodeNdx[5];
			ASSERT2(NewNodeNdx[5] == NodeTable->lookup_ndx(NewNodeKey[5]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;
		}
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side0_init,pt_start);

    //SIDE1
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        //@todo replace which with macro or enum
        //@todo reuse indexes for nodes and elements
        int which;
        ti_ndx_t n1_ndx, n2_ndx, n3_ndx, n4_ndx;

        ti_ndx_t *ndxNodeTemp=&(node_ndx_ref[iElm][0]);
        SFC_Key *NewNodeKey=&(new_node_key[iElm][0]);
        ti_ndx_t *NewNodeNdx=&(new_node_ndx[iElm][0]);

        int i;
        ti_ndx_t neigh_elm_ndx;
        int RefinedNeigh = 0;
        int info;
        int other_proc = 0;
        int boundary;

	    //+++++++++++++++++++++++++++SIDE1

		if(ElemTable->neigh_proc_[1][ndx] == -1)
			boundary = 1;
		else
			boundary = 0;

		if(boundary == 1 || ElemTable->neigh_gen_[1][ndx] <= ElemTable->generation_[ndx])
		{
			RefinedNeigh = 0;
			info = S_S_CON;
			if(ElemTable->neigh_gen_[1][ndx] != myid)    // && *(EmTemp->get_neigh_proc()+1)>0)
			{
				other_proc = 1;
				info = -1;
			}
			else
				other_proc = 0;

			//---Eight new node---
			which = 8;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);

			//---Fifth old node---
			if(RefinedNeigh || boundary)
				NodeTable->info_[ndxNodeTemp[5]]=CORNER;
			else if(other_proc)
				NodeTable->info_[ndxNodeTemp[5]]=info;
			else
				NodeTable->info_[ndxNodeTemp[5]]=S_C_CON;

			//---Thirteenth new node---
			which = 13;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);
			//check_create_new_node(which, 2, 5, ndxNodeTemp, NewNodeKey, NewNodeNdx, info, RefinedNeigh, boundary);
		}
		else
		{
			//Keith Added this if
			if((ElemTable->neigh_proc_[1][ndx] != myid) || ((ElemTable->neigh_proc_[5][ndx] != myid)
					&& (ElemTable->neigh_proc_[5][ndx] != -2)))
				other_proc = 1;
			else
				other_proc = 0;

			// eighth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[1][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[1][ndx]));

			n1_ndx =NewNodeNdx[8];
			ASSERT2(NewNodeNdx[8] == NodeTable->lookup_ndx(NewNodeKey[8]));

			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;

			// fifth old node
			NodeTable->info_[ndxNodeTemp[5]]=CORNER;
			if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
				NodeTable->info_[ndxNodeTemp[5]]=-1;

			// thirteenth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[5][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[5][ndx]));

			n1_ndx = NewNodeNdx[13];
			ASSERT2(NewNodeNdx[13] == NodeTable->lookup_ndx(NewNodeKey[13]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;
		}
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side1_init,pt_start);

    //SIDE2
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        //@todo replace which with macro or enum
        //@todo reuse indexes for nodes and elements
        int which;
        ti_ndx_t n1_ndx, n2_ndx, n3_ndx, n4_ndx;

        ti_ndx_t *ndxNodeTemp=&(node_ndx_ref[iElm][0]);
        SFC_Key *NewNodeKey=&(new_node_key[iElm][0]);
        ti_ndx_t *NewNodeNdx=&(new_node_ndx[iElm][0]);

        int i;
        ti_ndx_t neigh_elm_ndx;
        int RefinedNeigh = 0;
        int info;
        int other_proc = 0;
        int boundary;

		//+++++++++++++++++++++++++++SIDE2
		if(ElemTable->neigh_proc_[2][ndx] == -1)
			boundary = 1;
		else
			boundary = 0;

		if(boundary == 1 || ElemTable->neigh_gen_[2][ndx] <= ElemTable->generation_[ndx])
		{
			info = S_S_CON;

			if(ElemTable->neigh_gen_[2][ndx] != myid)    // && *(EmTemp->get_neigh_proc()+2)>0)
			{
				other_proc = 1;
				info = -1;
			}
			else
				other_proc = 0;

			RefinedNeigh = 0;

			//---Fourteenth new node---
			which = 14;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);

			//---Sixth old node---
			if(RefinedNeigh || boundary)
				NodeTable->info_[ndxNodeTemp[6]]=CORNER;
			else if(other_proc)
				NodeTable->info_[ndxNodeTemp[6]]=info;
			else
				NodeTable->info_[ndxNodeTemp[6]]=S_C_CON;

			//---Fifteenth new node---
			which = 15;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);
		}
		else
		{
			//Keith Added this if
			if((ElemTable->neigh_proc_[2][ndx] != myid) || ((ElemTable->neigh_proc_[6][ndx] != myid)
					&& (ElemTable->neigh_proc_[6][ndx] != -2)))
				other_proc = 1;
			else
				other_proc = 0;

			// fourteenth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[6][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[6][ndx]));

			n1_ndx = NewNodeNdx[14];
			ASSERT2(NewNodeNdx[14] == NodeTable->lookup_ndx(NewNodeKey[14]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;

			// sixth old node
			NodeTable->info_[ndxNodeTemp[6]]=CORNER;
			if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
				NodeTable->info_[ndxNodeTemp[6]]=-1;

			// fifteenth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[2][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[2][ndx]));

            n1_ndx = NewNodeNdx[15];
            ASSERT2(NewNodeNdx[15] == NodeTable->lookup_ndx(NewNodeKey[15]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;
		}
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side2_init,pt_start);

    //SIDE 3
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        //@todo replace which with macro or enum
        //@todo reuse indexes for nodes and elements
        int which;
        ti_ndx_t n1_ndx, n2_ndx, n3_ndx, n4_ndx;

        ti_ndx_t *ndxNodeTemp=&(node_ndx_ref[iElm][0]);
        SFC_Key *NewNodeKey=&(new_node_key[iElm][0]);
        ti_ndx_t *NewNodeNdx=&(new_node_ndx[iElm][0]);

        int i;
        ti_ndx_t neigh_elm_ndx;
        int RefinedNeigh = 0;
        int info;
        int other_proc = 0;
        int boundary;
		//+++++++++++++++++++++++++++SIDE 3
		if(ElemTable->neigh_proc_[3][ndx] == -1)
			boundary = 1;
		else
			boundary = 0;

		if(boundary == 1 || ElemTable->neigh_gen_[3][ndx] <= ElemTable->generation_[ndx])
		{
			info = S_S_CON;

			if(ElemTable->neigh_gen_[3][ndx] != myid)  //&& *(EmTemp->get_neigh_proc()+3)>0)
			{
				other_proc = 1;
				info = -1;
			}
			else
				other_proc = 0;

			RefinedNeigh = 0;

			//---Sixth new node----
			which = 6;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);
			//check_create_new_node(which, 0, 7, ndxNodeTemp, NewNodeKey, NewNodeNdx, info, RefinedNeigh, boundary);

			//---Seventh old node---
			if(RefinedNeigh || boundary)
				NodeTable->info_[ndxNodeTemp[7]]=CORNER;
			else if(other_proc)
				NodeTable->info_[ndxNodeTemp[7]]=-1;
			else
				NodeTable->info_[ndxNodeTemp[7]]=S_C_CON;

			//---Eleventh new node---
			which = 11;
			check_create_new_node2(iElm, which, info, RefinedNeigh, boundary);
			//n1 = (Node*) NodeTable->lookup(EmTemp->node_key(7));
			//n2 = (Node*) NodeTable->lookup(EmTemp->node_key(3));

			//check_create_new_node(which, 3, 7, ndxNodeTemp, NewNodeKey, NewNodeNdx, info, RefinedNeigh, boundary);
		}
		else
		{
			//Keith Added this if
			if((ElemTable->neigh_proc_[3][ndx] != myid) || ((ElemTable->neigh_proc_[7][ndx] != myid)
					&& (ElemTable->neigh_proc_[7][ndx] != -2)))
				other_proc = 1;
			else
				other_proc = 0;

			// sixth new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[7][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[7][ndx]));

			n1_ndx = NewNodeNdx[6];
			ASSERT2(NewNodeNdx[6] == NodeTable->lookup_ndx(NewNodeKey[6]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;

			// seventh old node
			NodeTable->info_[ndxNodeTemp[7]]=CORNER;
			if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
				NodeTable->info_[ndxNodeTemp[7]]=-1;
			// eleventh new node
			neigh_elm_ndx = ElemTable->neighbor_ndx_[3][ndx];
			ASSERT2(neigh_elm_ndx == ElemTable->lookup_ndx(ElemTable->neighbors_[3][ndx]));

            n1_ndx = NewNodeNdx[11];
            ASSERT2(NewNodeNdx[11] == NodeTable->lookup_ndx(NewNodeKey[11]));
			if(ElemTable->refined_[neigh_elm_ndx] == 0 || ElemTable->refined_[neigh_elm_ndx] == GHOST)
				NodeTable->info_[n1_ndx]=SIDE;
			else
				NodeTable->info_[n1_ndx]=S_C_CON;
		}
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_side3_init,pt_start);

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        int which;

        //changing the old bubble
		NodeTable->info_[ElemTable->node_bubble_ndx_[iElm]]=CORNER;

		//++++++++++++++++INTERNAL SIDE NODES 7, 8OLD, 12, 9, 10
		//---Seventh new node---
		which=7;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), SIDE, -3, matprops_ptr);

		//---Twelwth new node---
		which=12;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), SIDE, -3, matprops_ptr);

		//---Ninth new node---
		which=9;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), SIDE, -3, matprops_ptr);

		//---Tenth new node---
		which=10;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), SIDE, -3, matprops_ptr);

		//+++++++++++++++++++THE NEW BUBBLES 0, 1, 2, 3
		//---0th new node---
		which=0;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), BUBBLE, -3, matprops_ptr);

		//---1st new node---
		which=1;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), BUBBLE, -3, matprops_ptr);

		//---2nd new node---
		which=2;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), BUBBLE, -3, matprops_ptr);

		//---3rd new node---
		which=3;
		NodeTable->elenode_[new_node_ndx[iElm][which]].init(new_node_key[iElm][which], &(new_node_coord[iElm][which][0]), BUBBLE, -3, matprops_ptr);
    }
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_int_nodes_init,pt_start);

    //first we will create 4 new elements and then init it as we will need indexes of brothers during initiation
    new_sons_ndx.resize(numElemToRefine);
    ElemTable->groupCreateAddNode(new_sons_ndx,new_node_key,new_node_coord,
            new_node_ndx,new_node_isnew);

    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_new_elm_aloc,pt_start);
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(int iElm=0;iElm<numElemToRefine;++iElm)
    {
        ti_ndx_t ndx=ElemToRefine[iElm];
        //@todo replace which with macro or enum
        //@todo reuse indexes for nodes and elements
        int which;
        ti_ndx_t n1_ndx, n2_ndx, n3_ndx, n4_ndx;

        SFC_Key *NewNodeKey=&(new_node_key[iElm][0]);
        ti_ndx_t *NewNodeNdx=&(new_node_ndx[iElm][0]);

        int i;
        ti_ndx_t neigh_elm_ndx;
        int RefinedNeigh = 0;
        int info;
        int other_proc = 0;
        int boundary;

        ti_ndx_t *ndxSons=&(new_sons_ndx[iElm][0]);




		SFC_Key nodes[9];
		ti_ndx_t nodes_ndx[9];
		SFC_Key neigh[8];
		ti_ndx_t neigh_ndx[8];
		int neigh_proc[8];
		int generation = ElemTable->generation_[ndx] + 1;
		int neigh_gen[4];
		int material = ElemTable->material_[ndx];

		double coord[DIMENSION];
		//---0th new element---

        ti_ndx_t ndxQuad9P;

		//the nodes

		nodes[0] = ElemTable->node_key_[0][ndx];
		nodes[1] = ElemTable->node_key_[4][ndx];
		nodes[2] = ElemTable->key_[ndx];
		nodes[3] = ElemTable->node_key_[7][ndx];
		nodes[4] = NewNodeKey[4];
		nodes[5] = NewNodeKey[7];
		nodes[6] = NewNodeKey[9];
		nodes[7] = NewNodeKey[6];
		nodes[8] = NewNodeKey[0];

		nodes_ndx[0] = ElemTable->node_key_ndx_[0][ndx];
		nodes_ndx[1] = ElemTable->node_key_ndx_[4][ndx];
		nodes_ndx[2] = ElemTable->node_bubble_ndx_[ndx];
		nodes_ndx[3] = ElemTable->node_key_ndx_[7][ndx];
        nodes_ndx[4] = NewNodeNdx[4];
        nodes_ndx[5] = NewNodeNdx[7];
        nodes_ndx[6] = NewNodeNdx[9];
        nodes_ndx[7] = NewNodeNdx[6];
        nodes_ndx[8] = NewNodeNdx[0];

		n1_ndx = nodes_ndx[8];
		for(i = 0; i < DIMENSION; i++)
			coord[i] = NodeTable->coord_[i][n1_ndx];
		//neighbors
		neigh[0] = neigh[4] = ElemTable->neighbors_[0][ndx]; //Why is this ok if not ok for 3 down
		neigh[1] = neigh[5] = NewNodeKey[1];
		neigh[2] = neigh[6] = NewNodeKey[3];
		if(ElemTable->neigh_proc_[7][ndx]!= -2)
			neigh[3] = neigh[7] = ElemTable->neighbors_[7][ndx]; //This should be okay no matter what
		else
			neigh[3] = neigh[7] = ElemTable->neighbors_[3][ndx]; //This is only ok if neigh_proc==-2

		neigh_ndx[0] = neigh_ndx[4] = ElemTable->neighbor_ndx_[0][ndx]; //Why is this ok if not ok for 3 down
		neigh_ndx[1] = neigh_ndx[5] = ndxSons[1];
		neigh_ndx[2] = neigh_ndx[6] = ndxSons[3];
        if(ElemTable->neigh_proc_[7][ndx]!= -2)
            neigh_ndx[3] = neigh_ndx[7] = ElemTable->neighbor_ndx_[7][ndx]; //This should be okay no matter what
        else
            neigh_ndx[3] = neigh_ndx[7] = ElemTable->neighbor_ndx_[3][ndx]; //This is only ok if neigh_proc==-2

		//process of the neighbors

		neigh_proc[0] = ElemTable->neigh_proc_[0][ndx];
		neigh_proc[1] = myid;
		neigh_proc[2] = myid;
		if(ElemTable->neigh_proc_[7][ndx] != -2)
			neigh_proc[3] = ElemTable->neigh_proc_[7][ndx]; //depending if the neighboring element is already refined
		else
			neigh_proc[3] = ElemTable->neigh_proc_[3][ndx];

		neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;

		neigh_gen[0] = ElemTable->neigh_gen_[0][ndx];
		neigh_gen[1] = generation;
		neigh_gen[2] = generation;
		neigh_gen[3] = ElemTable->neigh_gen_[3][ndx];


		double err = ElemTable->el_error_[0][ndx] * .5; //added by jp oct11
		double sol = ElemTable->el_solution_[0][ndx] * .5; //added by jp oct11
		// son 0 can use elm_loc
		int iwetnodefather = ElemTable->iwetnode_[ndx];
		double Awetfather = ElemTable->Awet_[ndx];
		double dpson[2];
		dpson[0] = ElemTable->drypoint_[0][ndx] * 2 + 0.5;
		dpson[1] = ElemTable->drypoint_[1][ndx] * 2 + 0.5;

        int elm_loc[2], my_elm_loc[2];
        elm_loc[0] = 2 * ElemTable->elm_loc_[0][ndx];
        elm_loc[1] = 2 * ElemTable->elm_loc_[1][ndx];

		//init new element
		ndxQuad9P = ndxSons[0];
		ElemTable->elenode_[ndxQuad9P].init(nodes, nodes_ndx, neigh, neigh_ndx, neigh_proc, generation, elm_loc, NULL, neigh_gen, material,
							 ndx, coord, ElemTable, NodeTable, myid, matprops_ptr, iwetnodefather, Awetfather,
							 dpson);
		ElemTable->which_son_[ndxQuad9P]=0;  //--by jp, 0 means son 0
		ElemTable->elenode_[ndxQuad9P].putel_sq(sol, err);  //added by jp oct11


		//---1st new element---

		//the nodes

		nodes[0] = ElemTable->node_key_[4][ndx];
		nodes[1] = ElemTable->node_key_[1][ndx];
		nodes[2] = ElemTable->node_key_[5][ndx];
		nodes[3] = ElemTable->key_[ndx];
		nodes[4] = NewNodeKey[5];
		nodes[5] = NewNodeKey[8];
		nodes[6] = NewNodeKey[10];
		nodes[7] = NewNodeKey[7];
		nodes[8] = NewNodeKey[1];

		nodes_ndx[0] = ElemTable->node_key_ndx_[4][ndx];
		nodes_ndx[1] = ElemTable->node_key_ndx_[1][ndx];
		nodes_ndx[2] = ElemTable->node_key_ndx_[5][ndx];
		nodes_ndx[3] = ElemTable->node_bubble_ndx_[ndx];
		nodes_ndx[4] = NewNodeNdx[5];
		nodes_ndx[5] = NewNodeNdx[8];
		nodes_ndx[6] = NewNodeNdx[10];
		nodes_ndx[7] = NewNodeNdx[7];
		nodes_ndx[8] = NewNodeNdx[1];

		n1_ndx = nodes_ndx[8];
		for(i = 0; i < DIMENSION; i++)
			coord[i] = NodeTable->coord_[i][n1_ndx];

		//neighbors
		if(ElemTable->neigh_proc_[4][ndx] != -2)
			neigh[0] = neigh[4] = ElemTable->neighbors_[4][ndx]; //this should be ok now matter what
		else
			neigh[0] = neigh[4] = ElemTable->neighbors_[0][ndx]; //this is only ok if neigh_proc==-2
		neigh[1] = neigh[5] = ElemTable->neighbors_[1][ndx];
		neigh[2] = neigh[6] = NewNodeKey[2];
		neigh[3] = neigh[7] = NewNodeKey[0];

		if(ElemTable->neigh_proc_[4][ndx] != -2)
		    neigh_ndx[0] = neigh_ndx[4] = ElemTable->neighbor_ndx_[4][ndx]; //this should be ok now matter what
        else
            neigh_ndx[0] = neigh_ndx[4] = ElemTable->neighbor_ndx_[0][ndx]; //this is only ok if neigh_proc==-2
		neigh_ndx[1] = neigh_ndx[5] = ElemTable->neighbor_ndx_[1][ndx];
		neigh_ndx[2] = neigh_ndx[6] = ndxSons[2];
        neigh_ndx[3] = neigh_ndx[7] = ndxSons[0];

		//process of the neighbors

		neigh_proc[0] = (ElemTable->neigh_proc_[4][ndx] != -2) ? ElemTable->neigh_proc_[4][ndx] : ElemTable->neigh_proc_[0][ndx];
		neigh_proc[1] = ElemTable->neigh_proc_[1][ndx];
		neigh_proc[2] = myid;
		neigh_proc[3] = myid;

		neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;

		neigh_gen[0] = ElemTable->neigh_gen_[0][ndx];
		neigh_gen[1] = ElemTable->neigh_gen_[1][ndx];
		neigh_gen[2] = generation;
		neigh_gen[3] = generation;

		my_elm_loc[0] = elm_loc[0] + 1;
		my_elm_loc[1] = elm_loc[1];
		dpson[0] = ElemTable->drypoint_[0][ndx] * 2 - 0.5;
		dpson[1] = ElemTable->drypoint_[1][ndx] * 2 + 0.5;

        //init new element
        ndxQuad9P = ndxSons[1];
        ElemTable->elenode_[ndxQuad9P].init(nodes, nodes_ndx, neigh, neigh_ndx, neigh_proc, generation, my_elm_loc, NULL, neigh_gen, material,
							 ndx, coord, ElemTable, NodeTable, myid, matprops_ptr, iwetnodefather, Awetfather,
							 dpson);
		ElemTable->which_son_[ndxQuad9P]=1;  //--by jp
		ElemTable->elenode_[ndxQuad9P].putel_sq(sol, err);  //added by jp oct11


		//---2nd new element---

		//the nodes
		nodes[0] = ElemTable->key_[ndx];
		nodes[1] = ElemTable->node_key_[5][ndx];
		nodes[2] = ElemTable->node_key_[2][ndx];
		nodes[3] = ElemTable->node_key_[6][ndx];
		nodes[4] = NewNodeKey[10];
		nodes[5] = NewNodeKey[13];
		nodes[6] = NewNodeKey[15];
		nodes[7] = NewNodeKey[12];
		nodes[8] = NewNodeKey[2];

		nodes_ndx[0] = ElemTable->node_bubble_ndx_[ndx];
        nodes_ndx[1] = ElemTable->node_key_ndx_[5][ndx];
        nodes_ndx[2] = ElemTable->node_key_ndx_[2][ndx];
        nodes_ndx[3] = ElemTable->node_key_ndx_[6][ndx];
        nodes_ndx[4] = NewNodeNdx[10];
        nodes_ndx[5] = NewNodeNdx[13];
        nodes_ndx[6] = NewNodeNdx[15];
        nodes_ndx[7] = NewNodeNdx[12];
        nodes_ndx[8] = NewNodeNdx[2];

		n1_ndx = nodes_ndx[8];
		for(i = 0; i < DIMENSION; i++)
			coord[i] = NodeTable->coord_[i][n1_ndx];

		//neighbors
		neigh[0] = neigh[4] = NewNodeKey[1];
		if(ElemTable->neigh_proc_[5][ndx] != -2)
			neigh[1] = neigh[5] = ElemTable->neighbors_[5][ndx]; //This should be ok no matter what
		else
			neigh[1] = neigh[5] = ElemTable->neighbors_[1][ndx]; //this is only ok is neigh_proc==-2
		neigh[2] = neigh[6] = ElemTable->neighbors_[2][ndx];
		neigh[3] = neigh[7] = NewNodeKey[3];

		neigh_ndx[0] = neigh_ndx[4] = ndxSons[1];
        if(ElemTable->neigh_proc_[5][ndx] != -2)
            neigh_ndx[1] = neigh_ndx[5] = ElemTable->neighbor_ndx_[5][ndx]; //This should be ok no matter what
        else
            neigh_ndx[1] = neigh_ndx[5] = ElemTable->neighbor_ndx_[1][ndx]; //this is only ok is neigh_proc==-2
        neigh_ndx[2] = neigh_ndx[6] = ElemTable->neighbor_ndx_[2][ndx];
        neigh_ndx[3] = neigh_ndx[7] = ndxSons[3];


		//process of the neighbors

		neigh_proc[0] = myid;
		neigh_proc[1] = (ElemTable->neigh_proc_[5][ndx] != -2) ? ElemTable->neigh_proc_[5][ndx] : ElemTable->neigh_proc_[1][ndx];
		neigh_proc[2] = ElemTable->neigh_proc_[2][ndx];
		neigh_proc[3] = myid;

		neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;

		neigh_gen[0] = generation;
		neigh_gen[1] = ElemTable->neigh_gen_[1][ndx];
		neigh_gen[2] = ElemTable->neigh_gen_[2][ndx];
		neigh_gen[3] = generation;

		my_elm_loc[0] = elm_loc[0] + 1;
		my_elm_loc[1] = elm_loc[1] + 1;
		dpson[0] = ElemTable->drypoint_[0][ndx] * 2 - 0.5;
		dpson[1] = ElemTable->drypoint_[1][ndx] * 2 - 0.5;

        //init new element
        ndxQuad9P = ndxSons[2];
        ElemTable->elenode_[ndxQuad9P].init(nodes, nodes_ndx, neigh, neigh_ndx, neigh_proc, generation, my_elm_loc, NULL, neigh_gen, material,
							 ndx, coord, ElemTable, NodeTable, myid, matprops_ptr, iwetnodefather, Awetfather,
							 dpson);
		ElemTable->which_son_[ndxQuad9P]=2;  //--by jp
		ElemTable->elenode_[ndxQuad9P].putel_sq(sol, err);  //added by jp oct11



		//---3rd new element---

		//the nodes
		nodes[0] = ElemTable->node_key_[7][ndx];
		nodes[1] = ElemTable->key_[ndx];
		nodes[2] = ElemTable->node_key_[6][ndx];
		nodes[3] = ElemTable->node_key_[3][ndx];
		nodes[4] = NewNodeKey[9];
		nodes[5] = NewNodeKey[12];
		nodes[6] = NewNodeKey[14];
		nodes[7] = NewNodeKey[11];
		nodes[8] = NewNodeKey[3];

		nodes_ndx[0] = ElemTable->node_key_ndx_[7][ndx];
        nodes_ndx[1] = ElemTable->node_bubble_ndx_[ndx];
        nodes_ndx[2] = ElemTable->node_key_ndx_[6][ndx];
        nodes_ndx[3] = ElemTable->node_key_ndx_[3][ndx];
        nodes_ndx[4] = NewNodeNdx[9];
        nodes_ndx[5] = NewNodeNdx[12];
        nodes_ndx[6] = NewNodeNdx[14];
        nodes_ndx[7] = NewNodeNdx[11];
        nodes_ndx[8] = NewNodeNdx[3];

		n1_ndx = nodes_ndx[8];
		for(i = 0; i < DIMENSION; i++)
			coord[i] = NodeTable->coord_[i][n1_ndx];

		//neighbors
		neigh[0] = neigh[4] = NewNodeKey[0];
		neigh[1] = neigh[5] = NewNodeKey[2];
		if(ElemTable->neigh_proc_[6][ndx] != -2)
			neigh[2] = neigh[6] = ElemTable->neighbors_[6][ndx];
		else
			neigh[2] = neigh[6] = ElemTable->neighbors_[2][ndx];
		neigh[3] = neigh[7] = ElemTable->neighbors_[3][ndx];

		neigh_ndx[0] = neigh_ndx[4] = ndxSons[0];
		neigh_ndx[1] = neigh_ndx[5] = ndxSons[2];
        if(ElemTable->neigh_proc_[6][ndx] != -2)
            neigh_ndx[2] = neigh_ndx[6] = ElemTable->neighbor_ndx_[6][ndx];
        else
            neigh_ndx[2] = neigh_ndx[6] = ElemTable->neighbor_ndx_[2][ndx];
        neigh_ndx[3] = neigh_ndx[7] = ElemTable->neighbor_ndx_[3][ndx];


		//process of the neighbors

		neigh_proc[0] = myid;
		neigh_proc[1] = myid;
		neigh_proc[2] = (ElemTable->neigh_proc_[6][ndx] != -2) ? ElemTable->neigh_proc_[6][ndx] : ElemTable->neigh_proc_[2][ndx];
		neigh_proc[3] = ElemTable->neigh_proc_[3][ndx];

		neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;

		neigh_gen[0] = generation;
		neigh_gen[1] = generation;
		neigh_gen[2] = ElemTable->neigh_gen_[2][ndx];
		neigh_gen[3] = ElemTable->neigh_gen_[3][ndx];

		my_elm_loc[0] = elm_loc[0];
		my_elm_loc[1] = elm_loc[1] + 1;
		dpson[0] = ElemTable->drypoint_[0][ndx] * 2 + 0.5;
		dpson[1] = ElemTable->drypoint_[1][ndx] * 2 - 0.5;

        //init new element
        ndxQuad9P = ndxSons[3];
        ElemTable->elenode_[ndxQuad9P].init(nodes, nodes_ndx, neigh, neigh_ndx, neigh_proc, generation, my_elm_loc, NULL, neigh_gen, material,
							 ndx, coord, ElemTable, NodeTable, myid, matprops_ptr, iwetnodefather, Awetfather,
							 dpson);
		ElemTable->which_son_[ndxQuad9P]=3;  //--by jp
		ElemTable->elenode_[ndxQuad9P].putel_sq(sol, err);  //added by jp oct11



		//---CHANGING THE FATHER---
		for(i = 0; i < 4; i++)
		{
		    ElemTable->son_[i][ndx]=NewNodeKey[i];
		    ElemTable->son_ndx_[i][ndx]=ndxSons[i];
		}
		// putting in brother info
		for(i = 0; i < 4; i++)
		{
			ElemTable->elenode_[ndxSons[i]].set_brothers(NewNodeKey);
			for(int j = 0; j < 4; j++)
			{
			    ElemTable->brothers_ndx_[j][ndxSons[i]]=ndxSons[j];
			}
		}

        adapted[ndx]=OLDFATHER;
        refined[ndx]=1;
	}
    PROFILING3_STOPADD_RESTART(HAdapt_refineElements_new_elm_init,pt_start);

	if(numprocs>1)ElemTable->update_neighbours_ndx_on_ghosts();
	PROFILING3_STOPADD_RESTART(HAdapt_refineElements_update_neighbours_ndx_on_ghosts,pt_start);

#ifdef DEB3
	//some checking
    /*for(ti_ndx_t ndx:allRefinement)
    {
        for(int i = 0; i < 4; ++i)
        {
            if(ElemTable->node_bubble_ndx_[ElemTable->son_ndx_[i][ndx]]!=NodeTable->lookup_ndx(ElemTable->key_[ElemTable->son_ndx_[i][ndx]]))
            {
                printf("RRR %d %d %d\n",i,ElemTable->node_bubble_ndx_[ElemTable->son_ndx_[i][ndx]],NodeTable->lookup_ndx(ElemTable->key_[ElemTable->son_ndx_[i][ndx]]));
            }
        }
    }*/
#endif
	return;
}
//()---new node numbering

//  3---(14)--6---(15)--2
//  |         |         |
//  |         |         |
// (11) (3)  (12) (2)  (13)
//  |         |         |
//  |         |         |
//  7---(9)---E---(10)--5
//  |         |         |
//  |         |         |
// (6)  (0)  (7)  (1)  (8)
//  |         |         |
//  |         |         |
//  0---(4)---4---(5)---1

//max & min coordinates need to be passed later now only for L-shape!!!
//only 4 one step because of the info FLAG!!!
//if the new node is on INTERFACE flag will be -1

void refine(Element* EmTemp, ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, MatProps* matprops_ptr)
{
    //printf("refining element %u %u \n",*(EmTemp->pass_key()), *(EmTemp->pass_key()+1));
    int which;
    Node *n1, *n2, *n3, *n4;
    
    Node* NodeTemp[9];
    SFC_Key NewNodeKey[16];
    Element* Quad9P;
    int numprocs, myid, i;
    Element* neigh_elm;
    //SFC_Key* neigh_node_key;
    int RefinedNeigh = 0;
    int info;
    int other_proc = 0;
    int boundary;
    int order[5];
    int NewOrder[4][5];
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int elm_loc[2], my_elm_loc[2];
    elm_loc[0] = 2 * EmTemp->elm_loc(0);
    elm_loc[1] = 2 * EmTemp->elm_loc(1);
    
    
    for(i = 0; i < 8; i++) //-- corners and sides
    {
        NodeTemp[i] = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(i));
        assert(NodeTemp[i]);
    }
    
    
    /*filling up the new order array
     str: side orders remain;
     newsides get the higher order of the already existing sides
     bubbles get the order of the old bubble
     */
    for(i = 0; i < 4; i++) //-- orders for the 4 sons
    {
        int a = i - 1;
        if(a == -1)
            a = 3;
        NewOrder[i][a] = order[a];
        NewOrder[i][i] = order[i];
        NewOrder[i][4] = order[4];
    }
    
    //for the new internal sides they get the order of the max side order
    for(i = 0; i < 4; i++)
    {
        int a = i + 1;
        int b = i + 2;
        int c = i - 1;
        if(i == 2)
            b = 0;
        if(i == 3)
        {
            a = 0;
            b = 1;
        } //in the case of the 3rd element
        if(c == -1)
            c = 3; //in the case of the 0th element
            
        NewOrder[i][a] = NewOrder[i][b] = order[0];
        
        for(int j = 1; j < 4; j++)
            if(order[j] > NewOrder[i][a])
                NewOrder[i][a] = NewOrder[i][b] = order[j];
    }
    
    NodeTemp[8] = (Node*) HT_Node_Ptr->lookup(EmTemp->key()); //-- bubble
            
    //SIDE 0
    if(EmTemp->neigh_proc(0) == -1)
        boundary = 1;
    else
        boundary = 0;
    
    if(boundary == 1 || EmTemp->neigh_gen(0) <= EmTemp->generation())
    {
        RefinedNeigh = 0;
        info = S_S_CON;
        if(EmTemp->neigh_proc(0) != myid)
        {
            other_proc = 1;
            info = -1;
        }
        
        which = 4;
        //---Fourth new node---
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(4));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(0));
        
        create_new_node(which, 0, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[0],
                        matprops_ptr);
        
        //---Fourth old node---
        if(RefinedNeigh || boundary)
            NodeTemp[4]->info(CORNER);
        else if(other_proc)
            NodeTemp[4]->info(-1);
        else
            NodeTemp[4]->info(S_C_CON);
        
        //---Fifth new node---
        which = 5;
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(1));
        
        create_new_node(which, 1, 4, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[0],
                        matprops_ptr);
    }
    else
    {
        
        //Keith Added this if
        if((EmTemp->neigh_proc(0) != myid) || ((EmTemp->neigh_proc(4) != myid)
                && (EmTemp->neigh_proc(4) != -2)))
            other_proc = 1;
        else
            other_proc = 0;
        
        // fourth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(0));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[4] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[4]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
        //fourth old node
        NodeTemp[4]->info(CORNER);
        if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
            NodeTemp[4]->info(-1);
        // fifth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(4));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[5] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[5]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
        
    }
    
//+++++++++++++++++++++++++++SIDE1
    
    if(EmTemp->neigh_proc(1) == -1)
        boundary = 1;
    else
        boundary = 0;
    
    if(boundary == 1 || EmTemp->neigh_gen(1) <= EmTemp->generation())
    {
        RefinedNeigh = 0;
        info = S_S_CON;
        if(EmTemp->neigh_proc(1) != myid)    // && *(EmTemp->get_neigh_proc()+1)>0)
        {
            other_proc = 1;
            info = -1;
        }
        else
            other_proc = 0;
        
        //---Eight new node---
        which = 8;
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(5));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(1));
        
        create_new_node(which, 1, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[1],
                        matprops_ptr);
        
        //---Fifth old node---
        if(RefinedNeigh || boundary)
            NodeTemp[5]->info(CORNER);
        else
        {
            if(other_proc)
                NodeTemp[5]->info(info);
            else
                NodeTemp[5]->info(S_C_CON);
        }
        
        //---Thirteenth new node---
        which = 13;
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(5));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(2));
        
        create_new_node(which, 2, 5, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[1],
                        matprops_ptr);
    }
    else
    {
        //Keith Added this if
        if((EmTemp->neigh_proc(1) != myid) || ((EmTemp->neigh_proc(5) != myid)
                && (EmTemp->neigh_proc(5) != -2)))
            other_proc = 1;
        else
            other_proc = 0;
        
        // eighth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(1));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[8] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
        // fifth old node
        NodeTemp[5]->info(CORNER);
        if(other_proc) //ERROR: other_proc is set based on side 0 neigbor not being more refined or never set, we never checked to see if the more refined neighbor was on another processor
            NodeTemp[5]->info(-1);
        // thirteenth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(5));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[13] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
    }
    
    //+++++++++++++++++++++++++++SIDE2
    
    if(EmTemp->neigh_proc(2) == -1)
        boundary = 1;
    else
        boundary = 0;
    
    if(boundary == 1 || EmTemp->neigh_gen(2) <= EmTemp->generation())
    {
        info = S_S_CON;
        
        if(EmTemp->neigh_proc(2) != myid)    // && *(EmTemp->get_neigh_proc()+2)>0)
        {
            other_proc = 1;
            info = -1;
        }
        else
            other_proc = 0;
        
        RefinedNeigh = 0;
        
        //---Fourteenth new node---
        which = 14;
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(3));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(6));
        
        create_new_node(which, 3, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[2],
                        matprops_ptr);
        
        //---Sixth old node---
        if(RefinedNeigh || boundary)
            NodeTemp[6]->info(CORNER);
        else if(other_proc)
            NodeTemp[6]->info(-1);
        else
            NodeTemp[6]->info(S_C_CON);
        
        //---Fifteenth new node---
        which = 15;
        // geoflow info
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(6));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(2));
        
        create_new_node(which, 2, 6, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[2],
                        matprops_ptr);
    }
    else
    {
        //Keith Added this if
        if((EmTemp->neigh_proc(2) != myid) || ((EmTemp->neigh_proc(6) != myid)
                && (EmTemp->neigh_proc(6) != -2)))
            other_proc = 1;
        else
            other_proc = 0;
        
        // fourteenth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(6));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[14] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[14]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
        // sixth old node
        NodeTemp[6]->info(CORNER);
        if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
            NodeTemp[6]->info(-1);
        // fifteenth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(2));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[15] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[15]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
    }
    
    //+++++++++++++++++++++++++++SIDE 3
    
    if(EmTemp->neigh_proc(3) == -1)
        boundary = 1;
    else
        boundary = 0;
    
    if(boundary == 1 || EmTemp->neigh_gen(3) <= EmTemp->generation())
    {
        info = S_S_CON;
        
        if(EmTemp->neigh_proc(3) != myid)  //&& *(EmTemp->get_neigh_proc()+3)>0)
        {
            other_proc = 1;
            info = -1;
        }
        else
            other_proc = 0;
        
        RefinedNeigh = 0;
        
        //---Sixth new node----
        which = 6;
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(7));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(0));
        
        create_new_node(which, 0, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[3],
                        matprops_ptr);
        
        //---Seventh old node---
        if(RefinedNeigh || boundary)
            NodeTemp[7]->info(CORNER);
        else if(other_proc)
            NodeTemp[7]->info(-1);
        else
            NodeTemp[7]->info(S_C_CON);
        
        //---Eleventh new node---
        which = 11;
        n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(7));
        n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(3));
        
        create_new_node(which, 3, 7, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, order[3],
                        matprops_ptr);
    }
    else
    {
        //Keith Added this if
        if((EmTemp->neigh_proc(3) != myid) || ((EmTemp->neigh_proc(7) != myid)
                && (EmTemp->neigh_proc(7) != -2)))
            other_proc = 1;
        else
            other_proc = 0;
        
        // sixth new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(7));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[6] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
        // seventh old node
        NodeTemp[7]->info(CORNER);
        if(other_proc) //ERROR: other_proc is never set, we never checked to see if the more refined neighbor was on another processor
            NodeTemp[7]->info(-1);
        // eleventh new node
        neigh_elm = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(3));
        i = 0;
        which = -1;
        while (i < 4 && which == -1)
        {
            if(neigh_elm->neighbor(i)==EmTemp->key())
                which = i;
            i++;
        }
        assert(which != -1);
        NewNodeKey[11] = neigh_elm->node_key(which + 4);
        n1 = (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
        if(neigh_elm->refined_flag() == 0 || neigh_elm->refined_flag() == GHOST)
            n1->info(SIDE);
        //else if(neigh_elm->get_refined_flag()==GHOST)
        //n1->putinfo(-1);
        else
            n1->info(S_C_CON);
    }
    //++++++++++++++++INTERNAL SIDE NODES 7, 8OLD, 12, 9, 10
    
    RefinedNeigh = 0;
    boundary = 0;
    info = SIDE;
    
    //---Seventh new node---
    
    which = 7;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(4));
    n2 = (Node*) HT_Node_Ptr->lookup(EmTemp->key());
    
    create_new_node(which, 4, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[0][1],
                    matprops_ptr);
    
    NodeTemp[8]->info(CORNER);    //changing the old bubble
            
    //---Twelwth new node---
    
    which = 12;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(6));
    
    create_new_node(which, 6, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[3][1],
                    matprops_ptr);
    
    //---Ninth new node---
    
    which = 9;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(7));
    
    create_new_node(which, 7, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[0][2],
                    matprops_ptr);
    
    //---Tenth new node---
    
    which = 10;
    // geoflow info
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(5));
    
    create_new_node(which, 5, 8, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[1][2],
                    matprops_ptr);
    
    //+++++++++++++++++++THE NEW BUBBLES 0, 1, 2, 3
    
    info = BUBBLE;
    
    //---0th new node---
    
    NodeTemp[0] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[6]);
    NodeTemp[1] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
    which = 0;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(0));
    n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(4));
    n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(7));
    
    create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[0][4],
                    matprops_ptr);
    
    //---1st new node---
    
    NodeTemp[0] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[7]);
    NodeTemp[1] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[8]);
    which = 1;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(1));
    n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(5));
    
    create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[1][4],
                    matprops_ptr);
    
    //---2nd new node---
    
    NodeTemp[0] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
    NodeTemp[1] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[13]);
    which = 2;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(2));
    n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(5));
    n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(6));
    
    create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[2][4],
                    matprops_ptr);
    
    //---3rd new node---
    
    NodeTemp[0] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[11]);
    NodeTemp[1] = (Node*) HT_Node_Ptr->lookup(NewNodeKey[12]);
    which = 3;
    n1 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(6));
    n3 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(3));
    n4 = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(7));
    
    create_new_node(which, 0, 1, HT_Node_Ptr, NodeTemp, NewNodeKey, info, &RefinedNeigh, boundary, NewOrder[3][4],
                    matprops_ptr);
    
    //---NEW ELEMENTS---
    
    SFC_Key nodes[9];
    SFC_Key neigh[8];
    int neigh_proc[8];
    int generation = EmTemp->generation() + 1;
    int neigh_gen[4];
    int material = EmTemp->material();
    
    double coord[DIMENSION];
    //---0th new element---
    
    //the nodes
    
    nodes[0] = EmTemp->node_key(0);
    nodes[1] = EmTemp->node_key(4);
    nodes[2] = EmTemp->key();
    nodes[3] = EmTemp->node_key(7);
    nodes[4] = NewNodeKey[4];
    nodes[5] = NewNodeKey[7];
    nodes[6] = NewNodeKey[9];
    nodes[7] = NewNodeKey[6];
    nodes[8] = NewNodeKey[0];

    n1 = (Node*) HT_Node_Ptr->lookup(nodes[8]);
    for(i = 0; i < DIMENSION; i++)
        coord[i] = n1->coord(i);
    //neighbors
    neigh[0] = neigh[4] = EmTemp->neighbor(0); //Why is this ok if not ok for 3 down
    neigh[1] = neigh[5] = NewNodeKey[1];
    neigh[2] = neigh[6] = NewNodeKey[3];
    if(EmTemp->neigh_proc(7) != -2)
        neigh[3] = neigh[7] = EmTemp->neighbor(7); //This should be okay no matter what
    else
        neigh[3] = neigh[7] = EmTemp->neighbor(3); //This is only ok if neigh_proc==-2
    
    //process of the neighbors
    
    neigh_proc[0] = EmTemp->neigh_proc(0);
    neigh_proc[1] = myid;
    neigh_proc[2] = myid;
    if(EmTemp->neigh_proc(7) != -2)
        neigh_proc[3] = EmTemp->neigh_proc(7); //depending if the neighboring element is already refined
    else
        neigh_proc[3] = EmTemp->neigh_proc(3);
    
    neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;
    
    neigh_gen[0] = EmTemp->neigh_gen(0);
    neigh_gen[1] = generation;
    neigh_gen[2] = generation;
    neigh_gen[3] = EmTemp->neigh_gen(3);
    
    double err = EmTemp->el_error(0) * .5; //added by jp oct11
    double sol = EmTemp->el_solution(0) * .5; //added by jp oct11
    // son 0 can use elm_loc
    int iwetnodefather = EmTemp->iwetnode();
    double Awetfather = EmTemp->Awet();
    double dpson[2];
    dpson[0] = EmTemp->drypoint(0) * 2 + 0.5;
    dpson[1] = EmTemp->drypoint(1) * 2 + 0.5;
    
    Element* old_elm = (Element*) HT_Elem_Ptr->lookup(nodes[8]);
    if(old_elm != NULL)
    {
        //old_elm->set_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
        //old_elm->void_bcptr();
        //HT_Elem_Ptr->removeElement(old_elm);
        old_elm->init(nodes, neigh, neigh_proc, generation, elm_loc, &NewOrder[0][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
        Quad9P = old_elm;
    }
    else{
        Quad9P = HT_Elem_Ptr->generateAddElement(nodes, neigh, neigh_proc, generation, elm_loc, &NewOrder[0][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
    }
    //double* state_vars = Quad9P->get_state_varsABCD();
    //printf("state_vars= %g   %g   %g\n",state_vars[0],state_vars[1],state_vars[2]);
    
    Quad9P->set_which_son(0);  //--by jp, 0 means son 0
            
    Quad9P->putel_sq(sol, err);  //added by jp oct11
    
    
    //---1st new element---
    
    //the nodes
    
    nodes[0] = EmTemp->node_key(4);
    nodes[1] = EmTemp->node_key(1);
    nodes[2] = EmTemp->node_key(5);
    nodes[3] = EmTemp->key();
    nodes[4] = NewNodeKey[5];
    nodes[5] = NewNodeKey[8];
    nodes[6] = NewNodeKey[10];
    nodes[7] = NewNodeKey[7];
    nodes[8] = NewNodeKey[1];

    n1 = (Node*) HT_Node_Ptr->lookup(nodes[8]);
    for(i = 0; i < DIMENSION; i++)
        coord[i] = n1->coord(i);
    
    //neighbors
    if(EmTemp->neigh_proc(4) != -2)
        neigh[0] = neigh[4] = EmTemp->neighbor(4); //this should be ok now matter what
    else
        neigh[0] = neigh[4] = EmTemp->neighbor(0); //this is only ok if neigh_proc==-2
    neigh[1] = neigh[5] = EmTemp->neighbor(1);
    neigh[2] = neigh[6] = NewNodeKey[2];
    neigh[3] = neigh[7] = NewNodeKey[0];
    
    //process of the neighbors
    
    neigh_proc[0] = (EmTemp->neigh_proc(4) != -2) ? EmTemp->neigh_proc(4) : EmTemp->neigh_proc(0);
    neigh_proc[1] = EmTemp->neigh_proc(1);
    neigh_proc[2] = myid;
    neigh_proc[3] = myid;
    
    neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;
    
    neigh_gen[0] = EmTemp->neigh_gen(0);
    neigh_gen[1] = EmTemp->neigh_gen(1);
    neigh_gen[2] = generation;
    neigh_gen[3] = generation;
    
    my_elm_loc[0] = elm_loc[0] + 1;
    my_elm_loc[1] = elm_loc[1];
    dpson[0] = EmTemp->drypoint(0) * 2 - 0.5;
    dpson[1] = EmTemp->drypoint(1) * 2 + 0.5;
    
    old_elm = (Element*) HT_Elem_Ptr->lookup(nodes[8]);
    if(old_elm != NULL)
    {
        //old_elm->set_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
        //old_elm->void_bcptr();
        //HT_Elem_Ptr->removeElement(old_elm);
        old_elm->init(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[1][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
        Quad9P = old_elm;
    }
    else{
        Quad9P = HT_Elem_Ptr->generateAddElement(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[1][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
    }
    //state_vars = Quad9P->get_state_varsABCD();
    //printf("state_vars= %g   %g   %g\n",state_vars[0],state_vars[1],state_vars[2]);
    
    Quad9P->set_which_son(1); //--by jp
            
    Quad9P->putel_sq(sol, err); //added by jp oct11
    
    
    //---2nd new element---
    
    //the nodes
    nodes[0] = EmTemp->key();
    nodes[1] = EmTemp->node_key(5);
    nodes[2] = EmTemp->node_key(2);
    nodes[3] = EmTemp->node_key(6);
    nodes[4] = NewNodeKey[10];
    nodes[5] = NewNodeKey[13];
    nodes[6] = NewNodeKey[15];
    nodes[7] = NewNodeKey[12];
    nodes[8] = NewNodeKey[2];

    n1 = (Node*) HT_Node_Ptr->lookup(nodes[8]);
    for(i = 0; i < DIMENSION; i++)
        coord[i] = n1->coord(i);
    
    //neighbors
    neigh[0] = neigh[4] = NewNodeKey[1];
    if(EmTemp->neigh_proc(5) != -2)
        neigh[1] = neigh[5] = EmTemp->neighbor(5); //This should be ok no matter what
    else
        neigh[1] = neigh[5] = EmTemp->neighbor(1); //this is only ok is neigh_proc==-2
    neigh[2] = neigh[6] = EmTemp->neighbor(2);
    neigh[3] = neigh[6] = NewNodeKey[3];
        
    
    //process of the neighbors
    
    neigh_proc[0] = myid;
    neigh_proc[1] = (EmTemp->neigh_proc(5) != -2) ? EmTemp->neigh_proc(5) : EmTemp->neigh_proc(1);
    neigh_proc[2] = EmTemp->neigh_proc(2);
    neigh_proc[3] = myid;
    
    neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;
    
    neigh_gen[0] = generation;
    neigh_gen[1] = EmTemp->neigh_gen(1);
    neigh_gen[2] = EmTemp->neigh_gen(2);
    neigh_gen[3] = generation;
    
    my_elm_loc[0] = elm_loc[0] + 1;
    my_elm_loc[1] = elm_loc[1] + 1;
    dpson[0] = EmTemp->drypoint(0) * 2 - 0.5;
    dpson[1] = EmTemp->drypoint(1) * 2 - 0.5;
    
    old_elm = (Element*) HT_Elem_Ptr->lookup(nodes[8]);
    if(old_elm != NULL)
    {
        //old_elm->set_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
        //old_elm->void_bcptr();
        //HT_Elem_Ptr->removeElement(old_elm);
        
        old_elm->init(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[2][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
        Quad9P = old_elm;
    }
    else{
        Quad9P = HT_Elem_Ptr->generateAddElement(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[2][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
    }
    //state_vars = Quad9P->get_state_varsABCD();
    //printf("state_vars= %g   %g   %g\n",state_vars[0],state_vars[1],state_vars[2]);
    
    Quad9P->set_which_son(2); //--by jp
            
    Quad9P->putel_sq(sol, err); //added by jp oct11
    
    
    
    //---3rd new element---
    
    //the nodes
    nodes[0] = EmTemp->node_key(7);
    nodes[1] = EmTemp->key();
    nodes[2] = EmTemp->node_key(6);
    nodes[3] = EmTemp->node_key(3);
    nodes[4] = NewNodeKey[9];
    nodes[5] = NewNodeKey[12];
    nodes[6] = NewNodeKey[14];
    nodes[7] = NewNodeKey[11];
    nodes[8] = NewNodeKey[3];
    
    n1 = (Node*) HT_Node_Ptr->lookup(nodes[8]);
    for(i = 0; i < DIMENSION; i++)
        coord[i] = n1->coord(i);
    
    //neighbors
    neigh[0] = neigh[4] = NewNodeKey[0];
    neigh[1] = neigh[5] = NewNodeKey[2];
    if(EmTemp->neigh_proc(6) != -2)
        neigh[2] = neigh[6] = EmTemp->neighbor(6);
    else
        neigh[2] = neigh[6] = EmTemp->neighbor(2);

    neigh[3] = neigh[6] = EmTemp->neighbor(3);

    
    //process of the neighbors
    
    neigh_proc[0] = myid;
    neigh_proc[1] = myid;
    neigh_proc[2] = (EmTemp->neigh_proc(6) != -2) ? EmTemp->neigh_proc(6) : EmTemp->neigh_proc(2);
    neigh_proc[3] = EmTemp->neigh_proc(3);
    
    neigh_proc[4] = neigh_proc[5] = neigh_proc[6] = neigh_proc[7] = -2;
    
    neigh_gen[0] = generation;
    neigh_gen[1] = generation;
    neigh_gen[2] = EmTemp->neigh_gen(2);
    neigh_gen[3] = EmTemp->neigh_gen(3);

    my_elm_loc[0] = elm_loc[0];
    my_elm_loc[1] = elm_loc[1] + 1;
    dpson[0] = EmTemp->drypoint(0) * 2 + 0.5;
    dpson[1] = EmTemp->drypoint(1) * 2 - 0.5;
    
    old_elm = (Element*) HT_Elem_Ptr->lookup(nodes[8]);
    if(old_elm != NULL)
    {
        //old_elm->set_adapted_flag(TOBEDELETED); //this line shouldn't be necessary just being redundantly careful
        //old_elm->void_bcptr();
        //HT_Elem_Ptr->removeElement(old_elm);
        old_elm->init(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[3][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
        Quad9P = old_elm;
    }
    else{
        Quad9P = HT_Elem_Ptr->generateAddElement(nodes, neigh, neigh_proc, generation, my_elm_loc, &NewOrder[3][0], neigh_gen, material,
                         EmTemp, coord, HT_Elem_Ptr, HT_Node_Ptr, myid, matprops_ptr, iwetnodefather, Awetfather,
                         dpson);
    }
    //state_vars = Quad9P->get_state_varsABCD();
    //printf("state_vars= %g   %g   %g\n\n",state_vars[0],state_vars[1],state_vars[2]);
    
    Quad9P->set_which_son(3); //--by jp
            
    Quad9P->putel_sq(sol, err); //added by jp oct11
    
    
    
    //---CHANGING THE FATHER---
    EmTemp->set_sons(NewNodeKey);
    // putting in brother info
    for(i = 0; i < 4; i++)
    {
        EmTemp = (Element*) HT_Elem_Ptr->lookup(NewNodeKey[i]);
        EmTemp->set_brothers(NewNodeKey);  //was  EmTemp->putbrothers(&NewNodeKey[i][0]);
    }
    
    return;
}

void create_new_node(int which, int Node1, int Node2, NodeHashTable* HT_Node_Ptr, Node* NodeTemp[],
                     SFC_Key NewNodeKey[], int info, int* RefNe, int boundary, int order,
                     MatProps* matprops_ptr)
{
    double NewNodeCoord[2];
    double norm_coord[2];
    unsigned u_norm_coord[2];
    unsigned nkey = 2;
    SFC_Key key;
    unsigned oldkey[KEYLENGTH];
    Node* NewNode;
    Node* p;
    static double XRange[2];
    static double YRange[2];
    int i;
    
    for(i = 0; i < 2; i++)
    {
        XRange[i] = *(HT_Node_Ptr->get_Xrange() + i);
        YRange[i] = *(HT_Node_Ptr->get_Yrange() + i);
    }
    
    for(i = 0; i < 2; i++)
        NewNodeCoord[i] = (NodeTemp[Node1]->coord(i) + NodeTemp[Node2]->coord(i)) * .5;
    
    norm_coord[0] = (NewNodeCoord[0] - XRange[0]) / (XRange[1] - XRange[0]);
    norm_coord[1] = (NewNodeCoord[1] - YRange[0]) / (YRange[1] - YRange[0]);
    
    fhsfc2d_(norm_coord, &nkey, oldkey);
    
    SET_NEWKEY(key,oldkey);
    NewNodeKey[which]=key;
    
    p = (Node*) HT_Node_Ptr->lookup(key);
    
    if(!p)

    {
        NewNode = HT_Node_Ptr->createAddNode(key, NewNodeCoord, info, order, matprops_ptr);
        p = NewNode;
        
    }
    
    else if(p->coord(0) != NewNodeCoord[0] || p->coord(1) != NewNodeCoord[1])
    {
        short same_key = 0;
        assert(same_key);
    }
    
    else
        *RefNe = 1;
    
    if(*RefNe || boundary)
        p->info(SIDE);
    
    return;
}

