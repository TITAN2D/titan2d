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
 * $Id: hashtab2.C 224 2011-12-04 20:49:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <assert.h>

#include "../header/boundary.h"
#include "../header/elenode.hpp"
#include "../header/titan2d_utils.h"

/*#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR*/
#include <mpi.h>
#include <limits.h>

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

ElementsHashTable *elementsHashTable;
NodeHashTable *nodeHashTable;

#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
template <typename T>
HashTable<T>::HashTable(double *doublekeyrangein, int size, double XR[], double YR[],tisize_t reserved_size)
: key_(reserved_size),status_(reserved_size),elenode_(reserved_size)
{
    int i;

    NBUCKETS = size;
    ENTRIES = 0;
    
    for(i = 0; i < KEYLENGTH; i++)
        doublekeyrange[i] = doublekeyrangein[i];
    
    hashconstant = 1.0 * NBUCKETS / (doublekeyrange[0] * doublekeyrange[1] + doublekeyrange[1]);
    
    bucket.resize(NBUCKETS);

    
    for(i = 0; i < 2; i++)
    {
        Xrange[i] = XR[i];
        Yrange[i] = YR[i];
    }
    
    invdxrange = 1.0 / (Xrange[1] - Xrange[0]);
    invdyrange = 1.0 / (Yrange[1] - Yrange[0]);
    
    //Keith made this change 20061109; and made hash an inline function
        /* NBUCKETS*2 is NBUCKETS*integer integer is empirical could be 1
         return (((int) ((key[0]*doublekeyrange[1]+key[1])/
         (doublekeyrange[0]*doublekeyrange[1]+doublekeyrange[1])*
         NBUCKETS*2+0.5) )%NBUCKETS);
         
        unsigned oldkey[KEYLENGTH];
        SET_OLDKEY(oldkey,key);
        return (((int) ((oldkey[0] * doublekeyrange[1] + oldkey[1]) * hashconstant + 0.5)) % NBUCKETS);*/    
}
template <typename T>
ti_ndx_t HashTable<T>::lookup_ndx(const SFC_Key& keyi)
{
    int entry = hash(keyi);
    int size = bucket[entry].key.size();
    
    SFC_Key *keyArr = &(bucket[entry].key[0]);
    int i;
    
    if(size == 0)
        return ti_ndx_doesnt_exist;
    if(keyi < keyArr[0])
        return ti_ndx_doesnt_exist;
    if(keyi > keyArr[size - 1])
        return ti_ndx_doesnt_exist;
    
    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(keyi == keyArr[i])
            {
                return bucket[entry].ndx[i];
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
            if(keyi > keyArr[i1])
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
            if(keyi == keyArr[i])
            {
                return bucket[entry].ndx[i];
            }
        }
    }
    return ti_ndx_doesnt_exist;
}
template <typename T>
T* HashTable<T>::lookup(const SFC_Key& keyi)
{
    ti_ndx_t ndx=lookup_ndx(keyi);
    /*if(status_[ndx]==CS_Removed)
    {
        assert(0);
        return nullptr;
    }*/
    if(ti_ndx_not_negative(ndx))
        return &(elenode_[ndx]);
    else
        return nullptr;
}

/*template <typename T>
T HashTable<T>::lookup(const SFC_Key& keyi)
{
    ti_ndx_t ndx=lookup_ndx(keyi);
    if(ti_ndx_not_negative(ndx))
        return elenode_[ndx];
    else
        return T;
}*/

template <typename T>
ti_ndx_t HashTable<T>::add_ndx(const SFC_Key& keyi)
{
    ti_ndx_t ndx=lookup_ndx(keyi);
    if(ti_ndx_not_negative(ndx))
        return ndx;
    
    int entry = hash(keyi);
    int entry_size = bucket[entry].key.size();

    //get space
    ndx=key_.push_back();
    elenode_.push_back();
    status_.push_back();
    
    //set values
    key_[ndx]=keyi;
    status_[ndx]=CS_Added;

    //place to hash table
    if(entry_size>0)
    {
        //this place is already occupied
        //find proper place to insert it
        int i;
        SFC_Key *keyArr = &(bucket[entry].key[0]);
        for(i=0;i<entry_size&&keyi>keyArr[i];++i){}

        bucket[entry].key.insert(bucket[entry].key.begin() + i, keyi);
        bucket[entry].ndx.insert(bucket[entry].ndx.begin() + i, ndx);
    }
    else
    {
        //will be first member of the bucket entry
        bucket[entry].key.push_back(keyi);
        bucket[entry].ndx.push_back(ndx);
    }
    ENTRIES+=1;

    elenode_[ndx].ndx(ndx);
    elenode_[ndx].set_key(keyi);
    return ndx;
}
template <typename T>
T* HashTable<T>::add(const SFC_Key& keyi)
{
    ti_ndx_t ndx=add_ndx(keyi);
    return &(elenode_[ndx]);
}

template <typename T>
void HashTable<T>::remove(const SFC_Key& keyi)
{
    int entry = hash(keyi);
    int entry_size = bucket[entry].key.size();
    ti_ndx_t bucket_entry_ndx=bucket[entry].lookup_local_ndx(keyi);
    
    if(ti_ndx_not_negative(bucket_entry_ndx))
    {
        //set status
        status_[bucket[entry].ndx[bucket_entry_ndx]]=CS_Removed;
        //delete        
        bucket[entry].key.erase(bucket[entry].key.begin() + bucket_entry_ndx);
        bucket[entry].ndx.erase(bucket[entry].ndx.begin() + bucket_entry_ndx);
    }
    ENTRIES-=1;
}
template <typename T>
void HashTable<T>::flushTable()
{
    int size_old=key_.size();
    ndx_map.resize(size_old);
    key_map.resize(size_old);
    ndx_map_old.resize(size_old);
    
    //key_.equate_reserved_size();
    //status_.equate_reserved_size();
    
    int j=0;
    for(int i=0;i<size_old;++i)
    {
        if(status_[i]<0)continue;
        key_map[j]=key_map[i];
        ndx_map[j]=i;
        ++j;
    }
    int size=j;
    ndx_map.resize(size);
    key_map.resize(size);
    
    sort(ndx_map.begin(), ndx_map.end(),
        [&](const int& a, const int& b) {
            return (key_map[a] < key_map[b]);
        }
    );
    for(int i=0;i<size_old;++i)
    {
        ndx_map_old[i]=ti_ndx_doesnt_exist;
    }
    for(int i=0;i<size;++i)
    {
        ndx_map_old[ndx_map[i]]=i;
    }
    
    //reoder content
    key_.reorder(&(ndx_map[0]), size);
    elenode_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<size;++i)
    {
        elenode_[i].ndx(i);
    }
    
    status_.resize(size,false);
    status_.set(CS_Permanent);
    //update baskets
    
    
    for(int entry=0;entry<bucket.size();++entry)
    {
        for(int i=0;i<bucket[entry].ndx.size();++i)
        {
            assert(ndx_map_old[bucket[entry].ndx[i]]!=ti_ndx_doesnt_exist);
            bucket[entry].ndx[i]=ndx_map_old[bucket[entry].ndx[i]];
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
NodeHashTable::NodeHashTable(double *doublekeyrangein, int size, double XR[], double YR[])
    :HashTable<Node>(doublekeyrangein, size, XR, YR,node_reserved_size)
{
    nodeHashTable=this;
}
NodeHashTable::~NodeHashTable()
{
    nodeHashTable=nullptr;
}

Node* NodeHashTable::addNode(const SFC_Key& keyi)
{
    Node* node = add(keyi);
    
    id_.push_back();
    num_assoc_elem_.push_back();
    info_.push_back();
    for(int i=0;i<DIMENSION;i++)
        coord_[i].push_back();
    elevation_.push_back();
    for(int i=0;i<NUM_STATE_VARS;i++)
    {
        flux_[i].push_back();
        refinementflux_[i].push_back();
    }
    connection_id_.push_back();
    return node;
}

ti_ndx_t NodeHashTable::addNode_ndx(const SFC_Key& keyi)
{
    ti_ndx_t ndx = add_ndx(keyi);
    
    id_.push_back();
    num_assoc_elem_.push_back();
    info_.push_back();
    for(int i=0;i<DIMENSION;i++)
        coord_[i].push_back();
    elevation_.push_back();
    for(int i=0;i<NUM_STATE_VARS;i++)
    {
        flux_[i].push_back();
        refinementflux_[i].push_back();
    }
    connection_id_.push_back();
    return ndx;
}

Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr)
{
    Node* node = addNode(keyi);
    node->init(keyi, coordi, matprops_ptr);
    return node;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr)
{
    Node* node = addNode(keyi);
    node->init(keyi, coordi, inf, ord, matprops_ptr);
    return node;
}
ti_ndx_t NodeHashTable::createAddNode_ndx(const SFC_Key& keyi, const double *coordi, const int inf, const MatProps *matprops_ptr)
{
    ti_ndx_t ndx=addNode_ndx(keyi);
    //@TODO remove ord
    elenode_[ndx].init(keyi, coordi, inf, -3, matprops_ptr);
    return ndx;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada)
{
    Node* node = addNode(keyi);
    node->init(keyi, coordi, inf, ord, elev, yada);
    return node;
}

Node* NodeHashTable::createAddNode(FILE* fp, MatProps* matprops_ptr) //for restart
{
    //not implemented
    assert(0);
    return nullptr;
    //Node* node = add(keyi);
    //node->init(node->key(), node);
    //return node;
}
void NodeHashTable::removeNode(Node* node)
{
	if(node->key()==65175920631581991ull)
		printf("delete something");
    remove(node->key());
    //delete node;
}
void NodeHashTable::flushNodeTable()
{
    double t_start = MPI_Wtime();
    flushTable();
    
    int size=ndx_map.size();
    
    id_.reorder(&(ndx_map[0]), size);
    num_assoc_elem_.reorder(&(ndx_map[0]), size);
    info_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;i++)
        coord_[i].reorder(&(ndx_map[0]), size);
    elevation_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<NUM_STATE_VARS;i++)
    {
        flux_[i].reorder(&(ndx_map[0]), size);
        refinementflux_[i].reorder(&(ndx_map[0]), size);
    }
    connection_id_.reorder(&(ndx_map[0]), size);
    
    titanTimings.flushNodeTableTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.flushNodeTableTime += MPI_Wtime() - t_start;
}


////////////////////////////////////////////////////////////////////////////////
ElementsHashTable::ElementsHashTable(double *doublekeyrangein, int size, double XR[], double YR[], NodeHashTable* nodeTable)
        :HashTable<Element>(doublekeyrangein, size, XR, YR,elem_reserved_size)
{
    NlocalElements = 0;
    NodeTable = nodeTable;
    elementsHashTable=this;
    
    if(NUM_STATE_VARS == 3)
        elementType_=ElementType::SinglePhase;
    else if(NUM_STATE_VARS == 6)
        elementType_=ElementType::TwoPhases;
    else
        elementType_=ElementType::UnknownElementType;
}

ElementsHashTable::~ElementsHashTable()              //evacuate the table
{
    elementsHashTable=nullptr;
}
void ElementsHashTable::updateLocalElements()
{
    int count = 0, NEntriesInBucket;
    
    ukeyLocalElements.resize(ENTRIES);
    localElements.resize(ENTRIES);
#ifdef HASH_PRESERVE_ORDER
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(elenode_[i].adapted_flag() > 0)
            {//if this element does not belong on this processor don't involve!!!
                ukeyLocalElements[count] = key_[i];
                localElements[count] = &(elenode_[i]);
                count++;
            }
        }
    }
#else
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0 && adapted_[i] > 0)
        {//if this element does not belong on this processor don't involve!!!
            ukeyLocalElements[count] = key_[i];
            localElements[count] = &(elenode_[i]);
            count++;
        }
    }
    /*for(int i = 0; i < bucket.size(); i++)
    {
        for(int j = 0; j < bucket[i].ndx.size(); j++)
        {
            if(elenode_[bucket[i].ndx[j]].adapted_flag() > 0)
            {//if this element does not belong on this processor don't involve!!!
                ukeyLocalElements[count] = bucket[i].key[j];
                localElements[count] = &(elenode_[bucket[i].ndx[j]]);
                count++;
            }
        }
    }*/
#endif
    NlocalElements = count;
    ukeyLocalElements.resize(NlocalElements);
    localElements.resize(NlocalElements);
}
int ElementsHashTable::ckeckLocalElementsPointers(const char *prefix)
{
    int count = 0, NEntriesInBucket, mismatch = 0;
#ifdef HASH_PRESERVE_ORDER
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(elenode_[i].adapted_flag() > 0)
            {//if this element does not belong on this processor don't involve!!!
                if(ukeyLocalElements[count] != key_[i] || localElements[count] != &(elenode_[i]))
                    mismatch++;
                localElements[count] = &(elenode_[i]);
                count++;
            }
                
        }
    }
#else
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(adapted_[i] > 0)
            {//if this element does not belong on this processor don't involve!!!
                if(ukeyLocalElements[count] != key_[i] || localElements[count] != &(elenode_[i]))
                    mismatch++;
                localElements[count] = &(elenode_[i]);
                count++;
            }
                
        }
    }
    /*for(int i = 0; i < bucket.size(); i++)
    {
        for(int j = 0; j < bucket[i].ndx.size(); j++)
        {
            if(elenode_[bucket[i].ndx[j]].adapted_flag() > 0)
            {//if this element does not belong on this processor don't involve!!!
                if(ukeyLocalElements[count] != bucket[i].key[j] || localElements[count] != &(elenode_[bucket[i].ndx[j]]))
                    mismatch++;
                count++;
            }
        }
    }*/
#endif
    if(mismatch > 0)
    {
        printf("%s WARNING: AllEntriesPointersLocal are out-dated. %d values pointers/keys do not match.\n", prefix,
               mismatch);
    }
    return mismatch;
}
void ElementsHashTable::updatePointersToNeighbours()
{
    for(ti_ndx_t ndx = 0; ndx < size(); ndx++)
    {
        if(status_[ndx]>=0)
        {
            //if(adapted_[ndx] > 0)
            {                
                for (int j = 0; j < 8; j++) {
                    ti_ndx_t neigh_ndx=lookup_ndx(neighbors_[j][ndx]);
                    neighbor_ndx_[j][ndx]=neigh_ndx;

                    if(ti_ndx_not_negative(neigh_ndx))
                        neighborPtr_[j][ndx] = &(elenode_[neigh_ndx]);
                    else
                        neighborPtr_[j][ndx] = nullptr;
                }
                
                for (int j = 0; j < 8; j++) {
                    ti_ndx_t node_ndx=NodeTable->lookup_ndx(node_key_[j][ndx]);
                    node_key_ndx_[j][ndx]=node_ndx;

                    if(ti_ndx_not_negative(node_ndx))
                        node_keyPtr_[j][ndx] = &(NodeTable->elenode_[node_ndx]);
                    else
                        node_keyPtr_[j][ndx] = nullptr;
                }
                node_bubble_ndx_[ndx] = NodeTable->lookup_ndx(key_[ndx]);
            }
        }
    }
    return;
}
int ElementsHashTable::checkPointersToNeighbours(const char *prefix)
{
    int count = 0;
    
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(elenode_[i].adapted_flag() > 0) //if this element does not belong on this processor don't involve!!!
                count += elenode_[i].check_neighbors_nodes_and_elements_pointers(this, NodeTable);
        }
    }
    
    if(count > 0)
        printf("%s WARNING: neighbors nodes and elements pointers mismatch to key. %d mismatched.\n", prefix, count);
    return count;
}
ti_ndx_t ElementsHashTable::addElement_ndx(const SFC_Key& keyi)
{
    ti_ndx_t ndx=add_ndx(keyi);
    
    myprocess_.push_back();
    generation_.push_back();
    opposite_brother_flag_.push_back();
    material_.push_back(); /*! ! ! THE MAT. FLAG ! ! !*/
    lb_weight_.push_back();
    lb_key_.push_back();
    for(int i=0;i<8;++i)node_key_[i].push_back();
    for(int i=0;i<8;++i)node_keyPtr_[i].push_back();
    for(int i=0;i<8;++i)node_key_ndx_[i].push_back();
    node_bubble_ndx_.push_back();
    for(int i=0;i<8;++i)neighbors_[i].push_back();
    for(int i=0;i<8;++i)neighborPtr_[i].push_back();
    for(int i=0;i<8;++i)neighbor_ndx_[i].push_back();
    father_.push_back();
    for(int i=0;i<4;++i)son_[i].push_back();
    for(int i=0;i<8;++i)neigh_proc_[i].push_back();
    for(int i=0;i<8;++i)neigh_gen_[i].push_back();
    bcptr_.push_back();
    ndof_.push_back();
    no_of_eqns_.push_back();
    for(int i=0;i<EQUATIONS;++i)el_error_[i].push_back();
    for(int i=0;i<EQUATIONS;++i)el_solution_[i].push_back();
    refined_.push_back();
    adapted_.push_back();
    which_son_.push_back();
    new_old_.push_back();
    for(int i=0;i<4;++i)brothers_[i].push_back();
    for(int i=0;i<DIMENSION;++i)coord_[i].push_back();
    for(int i=0;i<DIMENSION;++i)elm_loc_[i].push_back();
    for(int i=0;i<NUM_STATE_VARS;++i)state_vars_[i].push_back();
    for(int i=0;i<NUM_STATE_VARS;++i)prev_state_vars_[i].push_back();
    for(int i=0;i<NUM_STATE_VARS * DIMENSION;++i)d_state_vars_[i].push_back();
    shortspeed_.push_back();
    for(int i=0;i<DIMENSION;++i)dx_[i].push_back();
    positive_x_side_.push_back();
    for(int i=0;i<DIMENSION;++i)eigenvxymax_[i].push_back();
    for(int i=0;i<DIMENSION;++i)kactxy_[i].push_back();
    elevation_.push_back();
    for(int i=0;i<DIMENSION;++i)zeta_[i].push_back();
    for(int i=0;i<DIMENSION;++i)curvature_[i].push_back();
    for(int i=0;i<3;++i)gravity_[i].push_back();
    for(int i=0;i<DIMENSION;++i)d_gravity_[i].push_back();
    stoppedflags_.push_back();
    effect_bedfrict_.push_back();
    effect_tanbedfrict_.push_back();
    for(int i=0;i<2;++i)effect_kactxy_[i].push_back();
    for(int i=0;i<NUM_STATE_VARS;++i)Influx_[i].push_back();
    ithelem_.push_back();
    iwetnode_.push_back();
    Awet_.push_back();
    for(int i=0;i<2;++i)drypoint_[i].push_back();
    Swet_.push_back();
    
    return ndx;
}

Element* ElementsHashTable::addElement(const SFC_Key& keyi)
{
    ti_ndx_t ndx=addElement_ndx(keyi);
    return &(elenode_[ndx]);
}
Element* ElementsHashTable::generateAddElement(const SFC_Key& keyi)
{
    Element* elm=addElement(keyi);
    elm->init(keyi);
    return elm;
}
    
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC* b, int mat,
                                            int* elm_loc_in, double pile_height, int myid, const SFC_Key& opposite_brother)
{
    Element* elm=addElement(nodekeys[8]); //--using bubble key to represent the element
    elm->init(nodekeys, neigh, n_pro, b, mat, elm_loc_in, pile_height, myid, opposite_brother);
    return elm;    
}
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen,
                 int elm_loc_in[], int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
                 ElementsHashTable *El_Table, NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
                 double Awetfather, double *drypoint_in)
{
    
    Element* elm=addElement(nodekeys[8]); //--using bubble key to represent the element
    elm->init(nodekeys, neigh, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    return elm;
}
ti_ndx_t ElementsHashTable::generateAddElement_ndx(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen, int elm_loc_in[],
                int *ord, int gen_neigh[], int mat, ti_ndx_t fthTemp, double *coord_in, ElementsHashTable *El_Table,
                NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather, double Awetfather,
                double *drypoint_in)
{
    ti_ndx_t ndx=addElement_ndx(nodekeys[8]); //--using bubble key to represent the element
    elenode_[ndx].init(nodekeys, neigh, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    return ndx;
}
Element* ElementsHashTable::generateAddElement(Element* sons[], NodeHashTable* NodeTable, ElementsHashTable* El_Table, MatProps* matprops_ptr)
{
    Element* elm=addElement(sons[2]->node_key(0));
    elm->init(sons, NodeTable, El_Table, matprops_ptr);
    return elm;
}
Element* ElementsHashTable::generateAddElement(FILE* fp, NodeHashTable* NodeTable, MatProps* matprops_ptr, int myid)
{
    //not implemented
    assert(0);
    return nullptr;
}

void ElementsHashTable::removeElement(Element* elm)
{
    remove(elm->key());
    //delete elm;
}

void ElementsHashTable::flushElemTable()
{
    double t_start = MPI_Wtime();
    //return;
    //@todo proper index reordering
    flushTable();
    int size=ndx_map.size();
    myprocess_.reorder(&(ndx_map[0]), size);
    generation_.reorder(&(ndx_map[0]), size);
    opposite_brother_flag_.reorder(&(ndx_map[0]), size);
    material_.reorder(&(ndx_map[0]), size); /*! ! ! THE MAT. FLAG ! ! !*/
    lb_weight_.reorder(&(ndx_map[0]), size);
    lb_key_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)node_key_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)node_keyPtr_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)node_key_ndx_[i].reorder(&(ndx_map[0]), size);
    node_bubble_ndx_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neighbors_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neighborPtr_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neighbor_ndx_[i].reorder(&(ndx_map[0]), size);
    father_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<4;++i)son_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neigh_proc_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neigh_gen_[i].reorder(&(ndx_map[0]), size);
    bcptr_.reorder(&(ndx_map[0]), size);
    ndof_.reorder(&(ndx_map[0]), size);
    no_of_eqns_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<EQUATIONS;++i)el_error_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<EQUATIONS;++i)el_solution_[i].reorder(&(ndx_map[0]), size);
    refined_.reorder(&(ndx_map[0]), size);
    adapted_.reorder(&(ndx_map[0]), size);
    which_son_.reorder(&(ndx_map[0]), size);
    new_old_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<4;++i)brothers_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)coord_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)elm_loc_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<NUM_STATE_VARS;++i)state_vars_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<NUM_STATE_VARS;++i)prev_state_vars_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<NUM_STATE_VARS * DIMENSION;++i)d_state_vars_[i].reorder(&(ndx_map[0]), size);
    shortspeed_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)dx_[i].reorder(&(ndx_map[0]), size);
    positive_x_side_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)eigenvxymax_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)kactxy_[i].reorder(&(ndx_map[0]), size);
    elevation_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)zeta_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)curvature_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<3;++i)gravity_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<DIMENSION;++i)d_gravity_[i].reorder(&(ndx_map[0]), size);
    stoppedflags_.reorder(&(ndx_map[0]), size);
    effect_bedfrict_.reorder(&(ndx_map[0]), size);
    effect_tanbedfrict_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<2;++i)effect_kactxy_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<NUM_STATE_VARS;++i)Influx_[i].reorder(&(ndx_map[0]), size);
    ithelem_.reorder(&(ndx_map[0]), size);
    iwetnode_.reorder(&(ndx_map[0]), size);
    Awet_.reorder(&(ndx_map[0]), size);
    for(int i=0;i<2;++i)drypoint_[i].reorder(&(ndx_map[0]), size);
    Swet_.reorder(&(ndx_map[0]), size);
    titanTimings.flushElemTableTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.flushElemTableTime += MPI_Wtime() - t_start;
}
