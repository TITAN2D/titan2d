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
HashTable<T>::HashTable(tisize_t reserved_size)
: elenode_(reserved_size)
{
    all_elenodes_are_permanent=false;

}
template <typename T>
void HashTable<T>::init(double *doublekeyrangein, int size, double XR[], double YR[])
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
    
    all_elenodes_are_permanent=false;

#ifdef _OPENMP
    bucket_lock.resize(NBUCKETS);
    for(i=0;i<NBUCKETS;++i)
        omp_init_lock(&(bucket_lock[i]));
    omp_init_lock(&content_table_lock);
#endif

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
    //IF_OMP(omp_set_lock(&(bucket_lock[entry])));
    int size = bucket[entry].key.size();
    
    SFC_Key *keyArr = &(bucket[entry].key[0]);
    int i;
    
    if(size == 0)
    {
        //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }
    if(keyi < keyArr[0])
    {
        //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }
    if(keyi > keyArr[size - 1])
    {
        //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }
    
    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(keyi == keyArr[i])
            {
                //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
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
                //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
                return bucket[entry].ndx[i];
            }
        }
    }
    //IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
    return ti_ndx_doesnt_exist;
}
template <typename T>
ti_ndx_t HashTable<T>::lookup_ndx_locked(const SFC_Key& keyi)
{
    int entry = hash(keyi);
    IF_OMP(omp_set_lock(&(bucket_lock[entry])));
    int size = bucket[entry].key.size();

    SFC_Key *keyArr = &(bucket[entry].key[0]);
    int i;

    if(size == 0)
    {
        IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }
    if(keyi < keyArr[0])
    {
        IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }
    if(keyi > keyArr[size - 1])
    {
        IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
        return ti_ndx_doesnt_exist;
    }

    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(keyi == keyArr[i])
            {
                IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
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
                IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
                return bucket[entry].ndx[i];
            }
        }
    }
    IF_OMP(omp_unset_lock(&(bucket_lock[entry])));
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
    
    all_elenodes_are_permanent=false;
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
    return ndx;
}
template <typename T>
ti_ndx_t HashTable<T>::add_ndx_locked(const SFC_Key& keyi)
{
    ti_ndx_t ndx=lookup_ndx_locked(keyi);
    if(ti_ndx_not_negative(ndx))
        return ndx;

    all_elenodes_are_permanent=false;

    //get space
    IF_OMP(omp_set_lock(&content_table_lock));
    ndx=key_.push_back();
    elenode_.push_back();
    status_.push_back();
    IF_OMP(omp_unset_lock(&content_table_lock));

    int entry = hash(keyi);
    IF_OMP(omp_set_lock(&(bucket_lock[entry])));

    int entry_size = bucket[entry].key.size();

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

    IF_OMP(omp_unset_lock(&(bucket_lock[entry])));

    elenode_[ndx].ndx(ndx);
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
    all_elenodes_are_permanent=false;
}

template <typename T>
bool HashTable<T>::check_that_all_elenodes_are_permanent()
{
    bool results=true;
    for(int i = 0; i < status_.size(); i++)
    {
        results=results && (status_[i]==CS_Permanent);

    }
    return results;
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

    all_elenodes_are_permanent=true;
}
template <typename T>
void HashTable<T>::reserve_base(const tisize_t new_reserve_size)
{
    key_.reserve(new_reserve_size);
    elenode_.reserve(new_reserve_size);
    status_.reserve(new_reserve_size);
}
template <typename T>
void HashTable<T>::reserve_at_least_base(const tisize_t new_reserve_size)
{
    key_.reserve_at_least(new_reserve_size);
    elenode_.reserve_at_least(new_reserve_size);
    status_.reserve_at_least(new_reserve_size);
}

//explicit implementation
template class HashTable<Node>;
template class HashTable<Element>;
////////////////////////////////////////////////////////////////////////////////
NodeHashTable::NodeHashTable()
    :HashTable<Node>(node_reserved_size)
{
    nodeHashTable=this;
}
NodeHashTable::~NodeHashTable()
{
    nodeHashTable=nullptr;
}

void NodeHashTable::init(double *doublekeyrangein, int size, double XR[], double YR[])
{
    HashTable<Node>::init(doublekeyrangein, size, XR, YR);
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
ti_ndx_t NodeHashTable::addNode_ndx_locked(const SFC_Key& keyi)
{
    ti_ndx_t ndx = add_ndx_locked(keyi);

    IF_OMP(omp_set_lock(&content_table_lock));
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
    IF_OMP(omp_unset_lock(&content_table_lock));
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
ti_ndx_t NodeHashTable::createAddNode_ndx_locked(const SFC_Key& keyi, const double *coordi, const int inf, const MatProps *matprops_ptr)
{
    ti_ndx_t ndx=addNode_ndx_locked(keyi);
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
    remove(node->key());
}
void NodeHashTable::removeNode(const ti_ndx_t ndx)
{
    remove(key_[ndx]);
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
void NodeHashTable::reserve(const tisize_t new_reserve_size)
{
    reserve_base(new_reserve_size);
}
void NodeHashTable::reserve_at_least(const tisize_t new_reserve_size)
{
    reserve_at_least_base(new_reserve_size);
}

////////////////////////////////////////////////////////////////////////////////
ElementsHashTable::ElementsHashTable(NodeHashTable* nodeTable)
        :HashTable<Element>(elem_reserved_size)
{
    NlocalElements = 0;
    NodeTable = nodeTable;
    elementsHashTable=this;
    
    elementType_=ElementType::UnknownElementType;
    conformation=0;
}

ElementsHashTable::~ElementsHashTable()              //evacuate the table
{
    elementsHashTable=nullptr;
}
void ElementsHashTable::init(double *doublekeyrangein, int size, double XR[], double YR[])
{
    HashTable<Element>::init(doublekeyrangein, size, XR, YR);
}
void ElementsHashTable::set_element_type(const ElementType m_elementType)
{
    elementType_=m_elementType;

    if(elementType_==ElementType::SinglePhase)
    {
        assert(NUM_STATE_VARS == 3);
    }
    else if(elementType_==ElementType::TwoPhases)
    {
        assert(NUM_STATE_VARS == 6);
    }
    else
    {
        printf("Unknown type of element!\n");
        assert(0);
    }
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

void ElementsHashTable::update_neighbours_ndx_on_ghosts(const bool check_neigh_proc)
{
    bool update=false;
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    for(ti_ndx_t ndx = 0; ndx < size(); ndx++)
    {
        if(status_[ndx]>=0)
        {
            if(check_neigh_proc)
            {
                for(int ineigh = 0; ineigh < 8; ineigh++)
                {
                    if((neigh_proc_[ineigh][ndx] >= 0) && (neigh_proc_[ineigh][ndx] != myid))
                    {
                        update=true;
                        break;
                    }
                }
            }
            if(-6 < adapted_[ndx] && adapted_[ndx] < 0 )
            {
                update=true;
            }

            if(update)
            {
                for (int j = 0; j < 8; j++) {
                    neighbor_ndx_[j][ndx]=lookup_ndx(neighbors_[j][ndx]);
                }
                for (int j = 0; j < 4; j++) {
                    brothers_ndx_[j][ndx]=lookup_ndx(brothers_[j][ndx]);
                }
                for (int j = 0; j < 4; j++) {
                    son_ndx_[j][ndx]=lookup_ndx(son_[j][ndx]);
                }
                for (int j = 0; j < 8; j++) {
                    node_key_ndx_[j][ndx]=NodeTable->lookup_ndx(node_key_[j][ndx]);
                }
                node_bubble_ndx_[ndx] = NodeTable->lookup_ndx(key_[ndx]);



                //now the neighbours
                for (int i = 0; i < 8; i++)
                {
                    ti_ndx_t ndx2=neighbor_ndx_[i][ndx];
                    if(ti_ndx_not_negative(ndx2))
                    {
                        for (int j = 0; j < 8; j++) {
                            neighbor_ndx_[j][ndx2]=lookup_ndx(neighbors_[j][ndx2]);
                        }
                        for (int j = 0; j < 4; j++) {
                            brothers_ndx_[j][ndx2]=lookup_ndx(brothers_[j][ndx2]);
                        }
                        for (int j = 0; j < 4; j++) {
                            son_ndx_[j][ndx2]=lookup_ndx(son_[j][ndx2]);
                        }
                        for (int j = 0; j < 8; j++) {
                            node_key_ndx_[j][ndx2]=NodeTable->lookup_ndx(node_key_[j][ndx2]);
                        }
                        node_bubble_ndx_[ndx2] = NodeTable->lookup_ndx(key_[ndx2]);
                    }
                }
            }
        }
    }
    return;
}
void ElementsHashTable::updateBrothersIndexes(const bool onlyForNewElements)
{
    ContentStatus lowestType=CS_Permanent;
    if(onlyForNewElements)
        lowestType=CS_Added;

    for(ti_ndx_t ndx = 0; ndx < size(); ndx++)
    {
        if(status_[ndx]>=lowestType)
        {
            for (int j = 0; j < 4; j++) {
                brothers_ndx_[j][ndx]=lookup_ndx(brothers_[j][ndx]);
            }
        }
    }
    return;
}
void ElementsHashTable::updateNeighboursIndexes()
{
    for(ti_ndx_t ndx = 0; ndx < size(); ndx++)
    {
        if(status_[ndx]>=0)
        {
            //if(adapted_[ndx] > 0)
            {
                //elements
                for (int j = 0; j < 8; j++) {
                    ti_ndx_t neigh_ndx=lookup_ndx(neighbors_[j][ndx]);
                    neighbor_ndx_[j][ndx]=neigh_ndx;

                    if(ti_ndx_not_negative(neigh_ndx))
                        neighborPtr_[j][ndx] = &(elenode_[neigh_ndx]);
                    else
                        neighborPtr_[j][ndx] = nullptr;
                }
                for (int j = 0; j < 4; j++) {
                    brothers_ndx_[j][ndx]=lookup_ndx(brothers_[j][ndx]);
                    son_ndx_[j][ndx]=lookup_ndx(son_[j][ndx]);
                }
                father_ndx_[ndx]=lookup_ndx(father_[ndx]);

                //nodes
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
void ElementsHashTable::checkPointersToNeighbours(const ti_ndx_t ndx, int &count_elem_neigbours_ndx, int &count_elem_brothers_ndx, int &count_elem_sons_ndx, int &count_elem_father_ndx,
        int &count_node_ndx, const bool checkPointers, int &count_ptr, const bool checkBrothers, int &count)
{
    for (int j = 0; j < 8; j++) {
        ti_ndx_t neigh_ndx=lookup_ndx(neighbors_[j][ndx]);
        if(neighbor_ndx_[j][ndx]!=neigh_ndx)
        {
            if(count_elem_neigbours_ndx<10)
                printf("==neighbor== ndx: %d status:%d neighbor_ndx[%d]:%d lookup(neighbor):%d status:%d key:%" PRIu64 "\n",
                   ndx,status_[ndx],j,neighbor_ndx_[j][ndx],
                   lookup_ndx(neighbors_[j][ndx]),
                   status_[lookup_ndx(neighbors_[j][ndx])],
                   neighbors_[j][ndx]);
            ++count_elem_neigbours_ndx;
            ++count;
        }

        if(checkPointers)
        {
            if(ti_ndx_not_negative(neigh_ndx))
            {
                if(neighborPtr_[j][ndx] != &(elenode_[neigh_ndx]))
                {
                    ++count_ptr;
                    ++count;
                }
            }
            else
            {
                if(neighborPtr_[j][ndx] != nullptr)
                {
                    ++count_ptr;
                    ++count;
                }
            }
        }
    }
    for (int j = 0; j < 4; j++)
    {

        if(checkBrothers && brothers_ndx_[j][ndx]!=lookup_ndx(brothers_[j][ndx]))
        {
            if(count_elem_brothers_ndx<10)
                printf("==brothers== ndx: %d status:%d brothers_ndx[%d]:%d lookup(brother):%d status:%d key:%" PRIu64 "\n",
                   ndx,status_[ndx],j,brothers_ndx_[j][ndx],
                   lookup_ndx(brothers_[j][ndx]),
                   status_[lookup_ndx(brothers_[j][ndx])],
                   brothers_[j][ndx]);
            ++count_elem_brothers_ndx;
            ++count;
        }
        if(son_ndx_[j][ndx]!=lookup_ndx(son_[j][ndx]))
        {
            if(count_elem_sons_ndx<10)
                printf("==sons== ndx: %d status:%d son_ndx[%d]:%d lookup(son):%d status:%d key:%" PRIu64 "\n",
                   ndx,status_[ndx],j,son_ndx_[j][ndx],
                   lookup_ndx(son_[j][ndx]),
                   status_[lookup_ndx(son_[j][ndx])],
                   son_[j][ndx]);
            ++count_elem_sons_ndx;
            ++count;
        }
    }
    if(father_ndx_[ndx]!=lookup_ndx(father_[ndx]))
    {
        ++count_elem_father_ndx;
        ++count;
    }

    for (int j = 0; j < 8; j++) {
        ti_ndx_t node_ndx=NodeTable->lookup_ndx(node_key_[j][ndx]);
        if(node_key_ndx_[j][ndx]!=node_ndx)
        {
            ++count_node_ndx;
            ++count;
        }

        if(checkPointers)
        {
            if(ti_ndx_not_negative(node_ndx))
            {
                if(node_keyPtr_[j][ndx] != &(NodeTable->elenode_[node_ndx]))
                {
                    ++count_ptr;
                    ++count;
                }
            }
            else
            {
                if(node_keyPtr_[j][ndx] != nullptr)
                {
                    ++count_ptr;
                    ++count;
                }
            }
        }
    }
    if(node_bubble_ndx_[ndx] != NodeTable->lookup_ndx(key_[ndx]))
    {
        ++count_node_ndx;
        ++count;
    }
}
int ElementsHashTable::checkPointersToNeighbours(const char *prefix,const bool checkPointers,const bool checkNewElements, const bool checkBrothers)
{
    int count=0;
    int count_ptr = 0;
    int count_elem_ndx = 0;
    int count_elem_neigbours_ndx = 0;
    int count_elem_brothers_ndx = 0;
    int count_elem_sons_ndx = 0;
    int count_elem_father_ndx = 0;
    int count_node_ndx = 0;
    bool check_element;
    for(ti_ndx_t ndx = 0; ndx < size(); ndx++)
    {
        if(checkNewElements)
            check_element=status_[ndx]>=0;
        else
            check_element=status_[ndx]==CS_Permanent;


        if(check_element)
        {
            checkPointersToNeighbours(ndx,count_elem_neigbours_ndx, count_elem_brothers_ndx, count_elem_sons_ndx,count_elem_father_ndx,count_node_ndx,checkPointers,count_ptr,checkBrothers, count);
        }
    }

    if(count> 0)
    {
        count_elem_ndx=count_elem_neigbours_ndx+count_elem_brothers_ndx+count_elem_sons_ndx+count_elem_father_ndx;
        printf("%s WARNING: neighbors nodes and elements pointers mismatch to key (totally %d).\n", prefix,count);
        if(checkPointers)printf("%s WARNING: %d mismatched pointers.\n", prefix, count_ptr);
        printf("%s WARNING: %d mismatched indexes to elements.\n", prefix, count_elem_ndx);
        printf("%s WARNING: %d mismatched indexes to elements' neigbours.\n", prefix, count_elem_neigbours_ndx);
        printf("%s WARNING: %d mismatched indexes to elements' brothers.\n", prefix, count_elem_brothers_ndx);
        printf("%s WARNING: %d mismatched indexes to elements' sons.\n", prefix, count_elem_sons_ndx);
        printf("%s WARNING: %d mismatched indexes to elements' fathers.\n", prefix, count_elem_father_ndx);

        printf("%s WARNING: %d mismatched indexes to nodes.\n", prefix, count_node_ndx);
    }
    return count_ptr+count_elem_ndx;
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
    father_ndx_.push_back();
    for(int i=0;i<4;++i)son_[i].push_back();
    for(int i=0;i<4;++i)son_ndx_[i].push_back();
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
    for(int i=0;i<4;++i)brothers_ndx_[i].push_back();
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
ti_ndx_t ElementsHashTable::generateAddElement_ndx(const SFC_Key* nodekeys, const ti_ndx_t* nodes_ndx, const SFC_Key* neigh, const ti_ndx_t* neigh_ndx, int n_pro[], BC *b, int gen, int elm_loc_in[],
                int *ord, int gen_neigh[], int mat, ti_ndx_t fthTemp, double *coord_in, ElementsHashTable *El_Table,
                NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather, double Awetfather,
                double *drypoint_in)
{
    ti_ndx_t ndx=addElement_ndx(nodekeys[8]); //--using bubble key to represent the element
    elenode_[ndx].init(nodekeys, nodes_ndx, neigh, neigh_ndx, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    return ndx;
}
ti_ndx_t ElementsHashTable::generateAddElement_ndx(ti_ndx_t *sons_ndx, MatProps *matprops_ptr)
{
    ti_ndx_t ndx=addElement_ndx(node_key_[0][sons_ndx[2]]);
    elenode_[ndx].init(sons_ndx, NodeTable, this, matprops_ptr);
    return ndx;
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
    for(int i=0;i<8;++i)node_key_ndx_[i].reorder_ndx(&(ndx_map[0]),&(NodeTable->ndx_map_old[0]), size);
    node_bubble_ndx_.reorder_ndx(&(ndx_map[0]),&(NodeTable->ndx_map_old[0]), size);
    for(int i=0;i<8;++i)neighbors_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neighborPtr_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<8;++i)neighbor_ndx_[i].reorder_ndx(&(ndx_map[0]),&(ndx_map_old[0]), size);
    father_.reorder(&(ndx_map[0]), size);
    father_ndx_.reorder_ndx(&(ndx_map[0]),&(ndx_map_old[0]), size);
    for(int i=0;i<4;++i)son_[i].reorder(&(ndx_map[0]), size);
    for(int i=0;i<4;++i)son_ndx_[i].reorder_ndx(&(ndx_map[0]),&(ndx_map_old[0]), size);
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
    for(int i=0;i<4;++i)brothers_ndx_[i].reorder_ndx(&(ndx_map[0]),&(ndx_map_old[0]), size);
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

    for(ti_ndx_t ndx = 0; ndx < size; ndx++)
    {
        if(status_[ndx]>=0)
        {
            //elements
            for (int j = 0; j < 8; j++) {
                ti_ndx_t neigh_ndx=neighbor_ndx_[j][ndx];

                if(ti_ndx_not_negative(neigh_ndx))
                    neighborPtr_[j][ndx] = &(elenode_[neigh_ndx]);
                else
                    neighborPtr_[j][ndx] = nullptr;
            }
            //nodes
            for (int j = 0; j < 8; j++) {
                ti_ndx_t node_ndx=node_key_ndx_[j][ndx];

                if(ti_ndx_not_negative(node_ndx))
                    node_keyPtr_[j][ndx] = &(NodeTable->elenode_[node_ndx]);
                else
                    node_keyPtr_[j][ndx] = nullptr;
            }
        }
    }
    updateLocalElements();

    conformation++;
    titanTimings.flushElemTableTime += MPI_Wtime() - t_start;
    titanTimingsAlongSimulation.flushElemTableTime += MPI_Wtime() - t_start;
}
void ElementsHashTable::reserve(const tisize_t new_reserve_size)
{
    reserve_base(new_reserve_size);

    myprocess_.reserve(new_reserve_size);
    generation_.reserve(new_reserve_size);
    opposite_brother_flag_.reserve(new_reserve_size);
    material_.reserve(new_reserve_size); /*! ! ! THE MAT. FLAG ! ! !*/
    lb_weight_.reserve(new_reserve_size);
    lb_key_.reserve(new_reserve_size);
    for(int i=0;i<8;++i)node_key_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)node_keyPtr_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)node_key_ndx_[i].reserve(new_reserve_size);
    node_bubble_ndx_.reserve(new_reserve_size);
    for(int i=0;i<8;++i)neighbors_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)neighborPtr_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)neighbor_ndx_[i].reserve(new_reserve_size);
    father_.reserve(new_reserve_size);
    father_ndx_.reserve(new_reserve_size);
    for(int i=0;i<4;++i)son_[i].reserve(new_reserve_size);
    for(int i=0;i<4;++i)son_ndx_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)neigh_proc_[i].reserve(new_reserve_size);
    for(int i=0;i<8;++i)neigh_gen_[i].reserve(new_reserve_size);
    bcptr_.reserve(new_reserve_size);
    ndof_.reserve(new_reserve_size);
    no_of_eqns_.reserve(new_reserve_size);
    for(int i=0;i<EQUATIONS;++i)el_error_[i].reserve(new_reserve_size);
    for(int i=0;i<EQUATIONS;++i)el_solution_[i].reserve(new_reserve_size);
    refined_.reserve(new_reserve_size);
    adapted_.reserve(new_reserve_size);
    which_son_.reserve(new_reserve_size);
    new_old_.reserve(new_reserve_size);
    for(int i=0;i<4;++i)brothers_[i].reserve(new_reserve_size);
    for(int i=0;i<4;++i)brothers_ndx_[i].reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)coord_[i].reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)elm_loc_[i].reserve(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)state_vars_[i].reserve(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)prev_state_vars_[i].reserve(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS * DIMENSION;++i)d_state_vars_[i].reserve(new_reserve_size);
    shortspeed_.reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)dx_[i].reserve(new_reserve_size);
    positive_x_side_.reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)eigenvxymax_[i].reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)kactxy_[i].reserve(new_reserve_size);
    elevation_.reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)zeta_[i].reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)curvature_[i].reserve(new_reserve_size);
    for(int i=0;i<3;++i)gravity_[i].reserve(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)d_gravity_[i].reserve(new_reserve_size);
    stoppedflags_.reserve(new_reserve_size);
    effect_bedfrict_.reserve(new_reserve_size);
    effect_tanbedfrict_.reserve(new_reserve_size);
    for(int i=0;i<2;++i)effect_kactxy_[i].reserve(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)Influx_[i].reserve(new_reserve_size);
    ithelem_.reserve(new_reserve_size);
    iwetnode_.reserve(new_reserve_size);
    Awet_.reserve(new_reserve_size);
    for(int i=0;i<2;++i)drypoint_[i].reserve(new_reserve_size);
    Swet_.reserve(new_reserve_size);
}
void ElementsHashTable::reserve_at_least(const tisize_t new_reserve_size)
{
    reserve(new_reserve_size);
    return;
    reserve_at_least_base(new_reserve_size);

    myprocess_.reserve_at_least(new_reserve_size);
    generation_.reserve_at_least(new_reserve_size);
    opposite_brother_flag_.reserve_at_least(new_reserve_size);
    material_.reserve_at_least(new_reserve_size); /*! ! ! THE MAT. FLAG ! ! !*/
    lb_weight_.reserve_at_least(new_reserve_size);
    lb_key_.reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)node_key_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)node_keyPtr_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)node_key_ndx_[i].reserve_at_least(new_reserve_size);
    node_bubble_ndx_.reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)neighbors_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)neighborPtr_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)neighbor_ndx_[i].reserve_at_least(new_reserve_size);
    father_.reserve_at_least(new_reserve_size);
    father_ndx_.reserve_at_least(new_reserve_size);
    for(int i=0;i<4;++i)son_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<4;++i)son_ndx_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)neigh_proc_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<8;++i)neigh_gen_[i].reserve_at_least(new_reserve_size);
    bcptr_.reserve_at_least(new_reserve_size);
    ndof_.reserve_at_least(new_reserve_size);
    no_of_eqns_.reserve_at_least(new_reserve_size);
    for(int i=0;i<EQUATIONS;++i)el_error_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<EQUATIONS;++i)el_solution_[i].reserve_at_least(new_reserve_size);
    refined_.reserve_at_least(new_reserve_size);
    adapted_.reserve_at_least(new_reserve_size);
    which_son_.reserve_at_least(new_reserve_size);
    new_old_.reserve_at_least(new_reserve_size);
    for(int i=0;i<4;++i)brothers_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<4;++i)brothers_ndx_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)coord_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)elm_loc_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)state_vars_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)prev_state_vars_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS * DIMENSION;++i)d_state_vars_[i].reserve_at_least(new_reserve_size);
    shortspeed_.reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)dx_[i].reserve_at_least(new_reserve_size);
    positive_x_side_.reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)eigenvxymax_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)kactxy_[i].reserve_at_least(new_reserve_size);
    elevation_.reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)zeta_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)curvature_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<3;++i)gravity_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<DIMENSION;++i)d_gravity_[i].reserve_at_least(new_reserve_size);
    stoppedflags_.reserve_at_least(new_reserve_size);
    effect_bedfrict_.reserve_at_least(new_reserve_size);
    effect_tanbedfrict_.reserve_at_least(new_reserve_size);
    for(int i=0;i<2;++i)effect_kactxy_[i].reserve_at_least(new_reserve_size);
    for(int i=0;i<NUM_STATE_VARS;++i)Influx_[i].reserve_at_least(new_reserve_size);
    ithelem_.reserve_at_least(new_reserve_size);
    iwetnode_.reserve_at_least(new_reserve_size);
    Awet_.reserve_at_least(new_reserve_size);
    for(int i=0;i<2;++i)drypoint_[i].reserve_at_least(new_reserve_size);
    Swet_.reserve_at_least(new_reserve_size);
}

EleNodeRef::EleNodeRef(ElementsHashTable *_ElemTable, NodeHashTable* _NodeTable):
                ElemTable(_ElemTable),
                NodeTable(_NodeTable),

                elements_(ElemTable->elenode_),
                status_(ElemTable->status_),
                adapted_(ElemTable->adapted_),
                generation_(ElemTable->generation_),
                neigh_proc_(ElemTable->neigh_proc_),
                state_vars_(ElemTable->state_vars_),
                prev_state_vars_(ElemTable->prev_state_vars_),
                d_state_vars_(ElemTable->d_state_vars_),
                stoppedflags_(ElemTable->stoppedflags_),
                d_gravity_(ElemTable->d_gravity_),
                gravity_(ElemTable->gravity_),
                curvature_(ElemTable->curvature_),
                zeta_(ElemTable->zeta_),
                effect_bedfrict_(ElemTable->effect_bedfrict_),
                effect_kactxy_(ElemTable->effect_kactxy_),
                kactxy_(ElemTable->kactxy_),
                Influx_(ElemTable->Influx_),
                neighbor_ndx_(ElemTable->neighbor_ndx_),
                positive_x_side_(ElemTable->positive_x_side_),
                node_key_ndx_(ElemTable->node_key_ndx_),
                el_error_(ElemTable->el_error_),
                dx_(ElemTable->dx_),
                node_refinementflux_(_NodeTable->refinementflux_),
                node_flux_(_NodeTable->flux_),
                node_info_(_NodeTable->info_),
                material_(ElemTable->material_),
                myprocess_(ElemTable->myprocess_),
                iwetnode_(ElemTable->iwetnode_),
                Awet_(ElemTable->Awet_),
                drypoint_(ElemTable->drypoint_),
                Swet_(ElemTable->Swet_),
                eigenvxymax_(ElemTable->eigenvxymax_),
                coord_(ElemTable->coord_)

{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
}

