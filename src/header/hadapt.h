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

//////////////////////////////////////////////////////////////////////////////


//!PrimaryRefinementsFinder
class PrimaryRefinementsFinderLevelSet:public SeedRefinementsFinder
{
public:
	PrimaryRefinementsFinderLevelSet(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
	virtual ~PrimaryRefinementsFinderLevelSet(){}
	virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
public:
	double geo_target;
private:
	tivector<double> &el_error;

	vector< vector<ti_ndx_t> > loc_SeedRefinement;
};

//!BuferFirstLayerRefinementsFinder
class BuferFirstLayerRefinementsFinderLevelSet:public SeedRefinementsFinder
{
public:
    BuferFirstLayerRefinementsFinderLevelSet(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
    virtual ~BuferFirstLayerRefinementsFinderLevelSet(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    vector< vector<ti_ndx_t> > loc_SeedRefinement;
};


//!BuferNextLayerRefinementsFinder
class BuferNextLayerRefinementsFinderLevelSet:public SeedRefinementsFinder
{
public:
    BuferNextLayerRefinementsFinderLevelSet(ElementsHashTable* _ElemTable,NodeHashTable* _NodeTable, ElementsProperties* _ElemProp);
    virtual ~BuferNextLayerRefinementsFinderLevelSet(){}
    virtual void findSeedRefinements(vector<ti_ndx_t> &seedRefinement);
private:
    vector< vector<ti_ndx_t> > loc_SeedRefinement;
};


/////////////////////////////////////////////////////////////////////////////

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
    HAdapt(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable, ElementsProperties* _ElemProp,TimeProps* _timeprops, MatProps* _matprops, const int _num_buffer_layer);
    void adapt(int h_count, double target);
    
    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);


    

private:
    void refine2(SeedRefinementsFinder &seedRefinementsFinder);

    void findPrimaryRefinements(vector<ti_ndx_t> &primaryRefinement, const double geo_target);  // What is this???
    void findTriggeredRefinements(const vector<ti_ndx_t> &primaryRefinement, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);
    void findBuferFirstLayerRefinements(vector<ti_ndx_t> &primaryRefinement);  // What is this???
    void findBuferNextLayerRefinements(vector<ti_ndx_t> &primaryRefinement);   // What is this???

    void refineElements(const vector<ti_ndx_t> &allRefinement);

    void check_create_new_node(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info, int& RefNe, const int boundary);
    void check_create_new_node2(const int iElm, const int iNode, const int info, int& RefNe, const int boundary);

    void create_new_node3(const int which, const ti_ndx_t Node1, const ti_ndx_t Node2,
            SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info);

    void refinedNeighboursUpdate(const vector<ti_ndx_t> &allRefinement);
    //void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const array<double,2> Node1, const array<double,2> Node2);
    //void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const ti_ndx_t Node1, const ti_ndx_t Node2);
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
    vector<array<ti_ndx_t,4> > new_sons_ndx;


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
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const ti_ndx_t Node1, const ti_ndx_t Node2)
    {
        double norm_coord[2];
        unsigned nkey = 2;
        unsigned oldkey[KEYLENGTH];

        ti_ndx_t ndx;

        for(int i = 0; i < 2; i++)
            coord[i] = (NodeTable->coord_[i][Node1] + NodeTable->coord_[i][Node2]) * .5;

        norm_coord[0] = (coord[0] - NodeTable->Xrange[0]) / (NodeTable->Xrange[1] - NodeTable->Xrange[0]);
        norm_coord[1] = (coord[1] - NodeTable->Yrange[0]) / (NodeTable->Yrange[1] - NodeTable->Yrange[0]);

        fhsfc2d_(norm_coord, &nkey, oldkey);

        SET_NEWKEY(key,oldkey);

        //ASSERT3(ti_ndx_negative(NodeTable->lookup_ndx_locked(key)));
        return;
    }
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const array<double,2> Node1, const array<double,2> Node2)
    {
        double norm_coord[2];
        unsigned nkey = 2;
        unsigned oldkey[KEYLENGTH];

        ti_ndx_t ndx;

        for(int i = 0; i < 2; i++)
            coord[i] = (Node1[i] + Node2[i]) * .5;

        norm_coord[0] = (coord[0] - NodeTable->Xrange[0]) / (NodeTable->Xrange[1] - NodeTable->Xrange[0]);
        norm_coord[1] = (coord[1] - NodeTable->Yrange[0]) / (NodeTable->Yrange[1] - NodeTable->Yrange[0]);

        fhsfc2d_(norm_coord, &nkey, oldkey);

        SET_NEWKEY(key,oldkey);

        ASSERT3(ti_ndx_negative(NodeTable->lookup_ndx_locked(key)));
        return;
    }
    /**
     * helper function for refineElements
     *
     * calculate new node's coordinates and keys and set for allocation, new node will locate half way
     * between node0_ndx and node1_ndx nodes
     *
     * in case of MPI parallel will use lookup_ndx if neighbor element is on other processor
     *
     * input:
     * iElm - index of element to refine within all elements for refinement (in allRefinement)
     * which - the node new local index in respect to element (in new_node_key, etc)
     * ndx  - index of element to refine (in ElemTable)
     * neigh - the neighbour element local index in respect to element
     * neigh_elm_ndx - index of neighbour which will share new node
     * node0_ndx,node1_ndx - indexes of nodes (in ElemTable) where
     * ithread - thread id which do work
     *
     * output:
     * new_node_key[iElm][which],
     * new_node_coord[iElm][which]
     * create_node_ielm[ithread].push_back(iElm)
     * create_node_iwhich[ithread].push_back(which)
     */

    void rE__find_new_side_node_or_set_for_alloc(const int iElm, const int which, const ti_ndx_t ndx,const int neigh,
                                                    const ti_ndx_t node0_ndx,const ti_ndx_t node1_ndx,const int ithread)
    {
        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];

        calc_coord_and_key(new_node_key[iElm][which],new_node_coord[iElm][which], node0_ndx, node1_ndx);
        if(ElemTable->myprocess_[neigh_elm_ndx]==myid)
        {
            ASSERT2(ti_ndx_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(which);
        }
        else
        {
            new_node_ndx[iElm][which]=NodeTable->lookup_ndx(new_node_key[iElm][which]);
            if(ti_ndx_negative(new_node_ndx[iElm][which]))
            {
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
        }
    }
    /**
     * same as rE__find_new_side_node_set_for_alloc but also will check if neighbour is also will be refined and
     * will get node from there
     */
    void rE__find_new_side_node_or_set_for_alloc_refcheck(const int iElm, const int which, const ti_ndx_t ndx,const int neigh, const int neigh_which,
                                                    const ti_ndx_t node0_ndx,const ti_ndx_t node1_ndx,const int ithread)
    {

        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];
        calc_coord_and_key(new_node_key[iElm][which],new_node_coord[iElm][which], node0_ndx, node1_ndx);
        if(ElemTable->myprocess_[neigh_elm_ndx]==myid)
        {
            if(refining_elem_map[neigh_elm_ndx]==-1 || ElemTable->neigh_gen_[neigh][ndx] < ElemTable->generation_[ndx])
            {
                //i.e. neigbour is not refined or neighour is of lower generation
                ASSERT2(ti_ndx_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
            else
            {
                ASSERT2(ti_ndx_not_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
                new_node_ndx[iElm][which]=new_node_ndx[refining_elem_map[neigh_elm_ndx]][neigh_which];
            }
        }
        else
        {
            new_node_ndx[iElm][which]=NodeTable->lookup_ndx(new_node_key[iElm][which]);
            if(ti_ndx_negative(new_node_ndx[iElm][which]))
            {
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
        }
    }
    /**
     * helper function for refineElements
     *
     * find new side node (4,5) in case of neighbor generation higher then this one
     *
     * input:
     * iElm - index of element to refine within all elements for refinement (in allRefinement)
     * which - the node new local index in respect to element (in new_node_key, etc)
     * ndx  - index of element to refine (in ElemTable)
     * neigh - the neighbour element local index in respect to element
     * neigh_which - local index in respect to neighboring element of the which node
     *
     * output:
     * new_node_key[iElm][which]
     * new_node_ndx[iElm][which]
     */
    void rE__find_new_side_node(const int iElm, const int which, const ti_ndx_t ndx,const int neigh, const int neigh_which)
    {
        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];
        ASSERT2(neigh_which == which_neighbor(ndx,neigh_elm_ndx) + 4);
        new_node_key[iElm][which] = ElemTable->node_key_[neigh_which][neigh_elm_ndx];
        new_node_ndx[iElm][which] = ElemTable->node_key_ndx_[neigh_which][neigh_elm_ndx];
        for(int j=0;j<DIMENSION;++j)
            new_node_coord[iElm][which][j]=NodeTable->coord_[j][new_node_ndx[iElm][which]];
    }

    /**
     * helper function for refineElements
     *
     * allocate new nodes, set indexes and coordinates
     */
    void rE__alloc_new_nodes()
    {
        ti_ndx_t number_of_new_elenodes=0;
        #pragma omp parallel
        {
            int ithread=omp_get_thread_num();
            ti_ndx_t start=0;
            if(ithread==0)
                for(int jthread=0;jthread<threads_number;++jthread)
                    number_of_new_elenodes+=create_node_ielm[jthread].size();
            else
                for(int jthread=0;jthread<ithread;++jthread)
                    start+=create_node_ielm[jthread].size();
            #pragma omp barrier
            if(ithread==0)
            {
                create_node_ielm[0].resize(number_of_new_elenodes);
                create_node_iwhich[0].resize(number_of_new_elenodes);
            }
            #pragma omp barrier
            if(ithread!=0)
            {
                for(ti_ndx_t i=0;i<create_node_ielm[ithread].size();++i)
                {
                    create_node_ielm[0][start+i]=create_node_ielm[ithread][i];
                    create_node_iwhich[0][start+i]=create_node_iwhich[ithread][i];
                }
                create_node_ielm[ithread].resize(0);
                create_node_iwhich[ithread].resize(0);
            }
        }
        NodeTable->groupCreateAddNode(create_node_ielm[0], create_node_iwhich[0],new_node_key,new_node_coord,new_node_ndx,new_node_isnew);

        create_node_ielm[0].resize(0);
        create_node_iwhich[0].resize(0);
    }
};


//! this is the normal grid adaptive refinement function it also refreshes the flux sources based on the "LEVEL SET" method
class HAdapt_LevelSet
{
public:
    ElementsHashTable* ElemTable;
    NodeHashTable* NodeTable;
    ElementsProperties* ElemProp;
    MatProps* matprops_ptr;
    TimeProps* timeprops_ptr;
    int num_buffer_layer;
public:
    HAdapt_LevelSet(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable, ElementsProperties* _ElemProp,TimeProps* _timeprops, MatProps* _matprops, const int _num_buffer_layer);
    void adapt(int h_count, double target);

//    void refinewrapper2(MatProps* matprops_ptr, ElemPtrList *RefinedList, Element *EmTemp);




private:
    void refine2(SeedRefinementsFinder &seedRefinementsFinder);

    void findPrimaryRefinements(vector<ti_ndx_t> &primaryRefinement, const double geo_target);  // What is this???
    void findTriggeredRefinements(const vector<ti_ndx_t> &primaryRefinement, vector<int> &set_for_refinement,vector<ti_ndx_t> &allRefinement);
    void findBuferFirstLayerRefinements(vector<ti_ndx_t> &primaryRefinement);  // What is this???
    void findBuferNextLayerRefinements(vector<ti_ndx_t> &primaryRefinement);   // What is this???

    void refineElements(const vector<ti_ndx_t> &allRefinement);

    void check_create_new_node(const int which, const int Node1, const int Node2,const ti_ndx_t * ndxNodeTemp,
                     SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info, int& RefNe, const int boundary);
    void check_create_new_node2(const int iElm, const int iNode, const int info, int& RefNe, const int boundary);

    void create_new_node3(const int which, const ti_ndx_t Node1, const ti_ndx_t Node2,
            SFC_Key NewNodeKey[], ti_ndx_t NewNodeNdx[], const int info);

    void refinedNeighboursUpdate(const vector<ti_ndx_t> &allRefinement);
    //void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const array<double,2> Node1, const array<double,2> Node2);
    //void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const ti_ndx_t Node1, const ti_ndx_t Node2);
private:
    PrimaryRefinementsFinderLevelSet primaryRefinementsFinder;
    BuferFirstLayerRefinementsFinderLevelSet buferFirstLayerRefinementsFinder;
    BuferNextLayerRefinementsFinderLevelSet buferNextLayerRefinementsFinder;

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
    vector<array<ti_ndx_t,4> > new_sons_ndx;


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
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const ti_ndx_t Node1, const ti_ndx_t Node2)
    {
        double norm_coord[2];
        unsigned nkey = 2;
        unsigned oldkey[KEYLENGTH];

        ti_ndx_t ndx;

        for(int i = 0; i < 2; i++)
            coord[i] = (NodeTable->coord_[i][Node1] + NodeTable->coord_[i][Node2]) * .5;

        norm_coord[0] = (coord[0] - NodeTable->Xrange[0]) / (NodeTable->Xrange[1] - NodeTable->Xrange[0]);
        norm_coord[1] = (coord[1] - NodeTable->Yrange[0]) / (NodeTable->Yrange[1] - NodeTable->Yrange[0]);

        fhsfc2d_(norm_coord, &nkey, oldkey);

        SET_NEWKEY(key,oldkey);

        //ASSERT3(ti_ndx_negative(NodeTable->lookup_ndx_locked(key)));
        return;
    }
    void calc_coord_and_key(SFC_Key &key, array<double,2> &coord, const array<double,2> Node1, const array<double,2> Node2)
    {
        double norm_coord[2];
        unsigned nkey = 2;
        unsigned oldkey[KEYLENGTH];

        ti_ndx_t ndx;

        for(int i = 0; i < 2; i++)
            coord[i] = (Node1[i] + Node2[i]) * .5;

        norm_coord[0] = (coord[0] - NodeTable->Xrange[0]) / (NodeTable->Xrange[1] - NodeTable->Xrange[0]);
        norm_coord[1] = (coord[1] - NodeTable->Yrange[0]) / (NodeTable->Yrange[1] - NodeTable->Yrange[0]);

        fhsfc2d_(norm_coord, &nkey, oldkey);

        SET_NEWKEY(key,oldkey);

        ASSERT3(ti_ndx_negative(NodeTable->lookup_ndx_locked(key)));
        return;
    }
    /**
     * helper function for refineElements
     *
     * calculate new node's coordinates and keys and set for allocation, new node will locate half way
     * between node0_ndx and node1_ndx nodes
     *
     * in case of MPI parallel will use lookup_ndx if neighbor element is on other processor
     *
     * input:
     * iElm - index of element to refine within all elements for refinement (in allRefinement)
     * which - the node new local index in respect to element (in new_node_key, etc)
     * ndx  - index of element to refine (in ElemTable)
     * neigh - the neighbour element local index in respect to element
     * neigh_elm_ndx - index of neighbour which will share new node
     * node0_ndx,node1_ndx - indexes of nodes (in ElemTable) where
     * ithread - thread id which do work
     *
     * output:
     * new_node_key[iElm][which],
     * new_node_coord[iElm][which]
     * create_node_ielm[ithread].push_back(iElm)
     * create_node_iwhich[ithread].push_back(which)
     */

    void rE__find_new_side_node_or_set_for_alloc(const int iElm, const int which, const ti_ndx_t ndx,const int neigh,
                                                    const ti_ndx_t node0_ndx,const ti_ndx_t node1_ndx,const int ithread)
    {
        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];

        calc_coord_and_key(new_node_key[iElm][which],new_node_coord[iElm][which], node0_ndx, node1_ndx);
        if(ElemTable->myprocess_[neigh_elm_ndx]==myid)
        {
            ASSERT2(ti_ndx_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
            create_node_ielm[ithread].push_back(iElm);
            create_node_iwhich[ithread].push_back(which);
        }
        else
        {
            new_node_ndx[iElm][which]=NodeTable->lookup_ndx(new_node_key[iElm][which]);
            if(ti_ndx_negative(new_node_ndx[iElm][which]))
            {
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
        }
    }
    /**
     * same as rE__find_new_side_node_set_for_alloc but also will check if neighbour is also will be refined and
     * will get node from there
     */
    void rE__find_new_side_node_or_set_for_alloc_refcheck(const int iElm, const int which, const ti_ndx_t ndx,const int neigh, const int neigh_which,
                                                    const ti_ndx_t node0_ndx,const ti_ndx_t node1_ndx,const int ithread)
    {

        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];
        calc_coord_and_key(new_node_key[iElm][which],new_node_coord[iElm][which], node0_ndx, node1_ndx);
        if(ElemTable->myprocess_[neigh_elm_ndx]==myid)
        {
            if(refining_elem_map[neigh_elm_ndx]==-1 || ElemTable->neigh_gen_[neigh][ndx] < ElemTable->generation_[ndx])
            {
                //i.e. neigbour is not refined or neighour is of lower generation
                ASSERT2(ti_ndx_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
            else
            {
                ASSERT2(ti_ndx_not_negative(NodeTable->lookup_ndx(new_node_key[iElm][which])));
                new_node_ndx[iElm][which]=new_node_ndx[refining_elem_map[neigh_elm_ndx]][neigh_which];
            }
        }
        else
        {
            new_node_ndx[iElm][which]=NodeTable->lookup_ndx(new_node_key[iElm][which]);
            if(ti_ndx_negative(new_node_ndx[iElm][which]))
            {
                create_node_ielm[ithread].push_back(iElm);
                create_node_iwhich[ithread].push_back(which);
            }
        }
    }
    /**
     * helper function for refineElements
     *
     * find new side node (4,5) in case of neighbor generation higher then this one
     *
     * input:
     * iElm - index of element to refine within all elements for refinement (in allRefinement)
     * which - the node new local index in respect to element (in new_node_key, etc)
     * ndx  - index of element to refine (in ElemTable)
     * neigh - the neighbour element local index in respect to element
     * neigh_which - local index in respect to neighboring element of the which node
     *
     * output:
     * new_node_key[iElm][which]
     * new_node_ndx[iElm][which]
     */
    void rE__find_new_side_node(const int iElm, const int which, const ti_ndx_t ndx,const int neigh, const int neigh_which)
    {
        const ti_ndx_t neigh_elm_ndx = ElemTable->neighbor_ndx_[neigh][ndx];
        ASSERT2(neigh_which == which_neighbor(ndx,neigh_elm_ndx) + 4);
        new_node_key[iElm][which] = ElemTable->node_key_[neigh_which][neigh_elm_ndx];
        new_node_ndx[iElm][which] = ElemTable->node_key_ndx_[neigh_which][neigh_elm_ndx];
        for(int j=0;j<DIMENSION;++j)
            new_node_coord[iElm][which][j]=NodeTable->coord_[j][new_node_ndx[iElm][which]];
    }

    /**
     * helper function for refineElements
     *
     * allocate new nodes, set indexes and coordinates
     */
    void rE__alloc_new_nodes()
    {
        ti_ndx_t number_of_new_elenodes=0;
        #pragma omp parallel
        {
            int ithread=omp_get_thread_num();
            ti_ndx_t start=0;
            if(ithread==0)
                for(int jthread=0;jthread<threads_number;++jthread)
                    number_of_new_elenodes+=create_node_ielm[jthread].size();
            else
                for(int jthread=0;jthread<ithread;++jthread)
                    start+=create_node_ielm[jthread].size();
            #pragma omp barrier
            if(ithread==0)
            {
                create_node_ielm[0].resize(number_of_new_elenodes);
                create_node_iwhich[0].resize(number_of_new_elenodes);
            }
            #pragma omp barrier
            if(ithread!=0)
            {
                for(ti_ndx_t i=0;i<create_node_ielm[ithread].size();++i)
                {
                    create_node_ielm[0][start+i]=create_node_ielm[ithread][i];
                    create_node_iwhich[0][start+i]=create_node_iwhich[ithread][i];
                }
                create_node_ielm[ithread].resize(0);
                create_node_iwhich[ithread].resize(0);
            }
        }
        NodeTable->groupCreateAddNode(create_node_ielm[0], create_node_iwhich[0],new_node_key,new_node_coord,new_node_ndx,new_node_isnew);

        create_node_ielm[0].resize(0);
        create_node_iwhich[0].resize(0);
    }
};

//! class for unrefinement
class HAdaptUnrefine:public EleNodeRef
{
public:
    ElementsProperties* ElemProp;
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

    //! it finds this element's brothers, the element should be 0th brother
    void find_brothers_to_unrefine(ti_ndx_t ndx,double target, int ithread);
    void find_brothers_to_unrefine__create_new_fathers(double target);


    //! it prevents refinement when one or more of the brothers does not belong to this processor
    int check_unrefinement(ti_ndx_t ndx, double target);

private:
    //temporary arrays used during refinement
    vector< vector< array<ti_ndx_t,4> > > brothers_to_unrefine_ndx;

    vector<ti_ndx_t> NewFatherList;
    vector<ti_ndx_t> OtherProcUpdate;

    vector< vector<ti_ndx_t> > nodesToDelete;
    vector< vector<ti_ndx_t> > elementsToDelete;
private:
    int myid;
    int numprocs;
};




#endif	/* HADAPT_H */

