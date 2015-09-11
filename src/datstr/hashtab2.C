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
    
    hashconstant = 8.0 * NBUCKETS / (doublekeyrange[0] * doublekeyrange[1] + doublekeyrange[1]);
    
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
void HashTable<T>::check()
{
    for(int entry=0;entry<bucket.size();++entry)
    {
        for(int bucket_entry_ndx=0;bucket_entry_ndx<bucket[entry].ndx.size();++bucket_entry_ndx)
        {
            if(elenode_[bucket[entry].ndx[bucket_entry_ndx]].key()!=key_[bucket[entry].ndx[bucket_entry_ndx]])
            {
            	cout<<elenode_[bucket[entry].ndx[bucket_entry_ndx]].key()<<"\n";
            	cout<<bucket[entry].key[bucket_entry_ndx]<<"\n";
            	cout<<bucket[entry].ndx[bucket_entry_ndx]<<"\n";
                cout<<elenode_[bucket[entry].ndx[bucket_entry_ndx]].ndx()<<"\n";
            	cout<<key_[bucket[entry].ndx[bucket_entry_ndx]]<<"\n";
            	T* t=&(elenode_[bucket[entry].ndx[bucket_entry_ndx]]);
            	cout<<t->key()<<"\n";
            	assert(0);
            }
        }
    }
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
void HashTable<T>::flush()
{
    return;
    int size_old=key_.size();
    static vector<ti_ndx_t> ndx_map;
    static vector<ti_ndx_t> ndx_map_old;
    static vector<SFC_Key> key_map;
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


Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, matprops_ptr);
    
    int entry = hash(keyi);
    int entry_size = bucket[entry].key.size();
    ti_ndx_t bucket_entry_ndx=bucket[entry].lookup_local_ndx(keyi);
    
    assert(node->key()==keyi);
    assert(elenode_[node->ndx()].key()==keyi);
    assert(bucket[entry].key[bucket_entry_ndx]==keyi);
    assert(elenode_[bucket[entry].ndx[bucket_entry_ndx]].key()==keyi);
    assert(key_[node->ndx()]==keyi);
    
    return node;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, inf, ord, matprops_ptr);
    
    int entry = hash(keyi);
    int entry_size = bucket[entry].key.size();
    ti_ndx_t bucket_entry_ndx=bucket[entry].lookup_local_ndx(keyi);
    
    assert(node->key()==keyi);
    assert(elenode_[node->ndx()].key()==keyi);
    assert(bucket[entry].key[bucket_entry_ndx]==keyi);
    assert(elenode_[bucket[entry].ndx[bucket_entry_ndx]].key()==keyi);
    assert(key_[node->ndx()]==keyi);
    
    return node;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, inf, ord, elev, yada);
    
    int entry = hash(keyi);
    int entry_size = bucket[entry].key.size();
    ti_ndx_t bucket_entry_ndx=bucket[entry].lookup_local_ndx(keyi);
    
    assert(node->key()==keyi);
    assert(elenode_[node->ndx()].key()==keyi);
    assert(bucket[entry].key[bucket_entry_ndx]==keyi);
    assert(elenode_[bucket[entry].ndx[bucket_entry_ndx]].key()==keyi);
    assert(key_[node->ndx()]==keyi);
    
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
    for(int i = 0; i < bucket.size(); i++)
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
    }
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
    for(int i = 0; i < bucket.size(); i++)
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
    }
#endif
    if(mismatch > 0)
    {
        printf("%s WARNING: AllEntriesPointersLocal are out-dated. %d values pointers/keys do not match.\n", prefix,
               mismatch);
    }
    return mismatch;
}
void ElementsHashTable::updateElements()
{
    int count = 0, NEntriesInBucket;
    ukeyElements.resize(ENTRIES);
    elements.resize(ENTRIES);
    
#ifdef HASH_PRESERVE_ORDER
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            ukeyLocalElements[count] = key_[i];
            localElements[count] = &(elenode_[i]);
            count++;
        }
    }
#else
    for(int i = 0; i < bucket.size(); i++)
    {
        for(int j = 0; j < bucket[i].ndx.size(); j++)
        {
            ukeyElements[count] = bucket[i].key[j];
            elements[count] = &(elenode_[bucket[i].ndx[j]]);
            count++;
        }
    }
#endif
}
int ElementsHashTable::ckeckElementsPointers(const char *prefix)
{
    int i, j, count = 0, NEntriesInBucket, mismatch = 0;
    if(ENTRIES != ukeyElements.size())
    {
        printf("%s WARNING: AllEntriesPointers are out-dated, number of entries do not match.\n", prefix);
        return ukeyElements.size();
    }
#ifdef HASH_PRESERVE_ORDER
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(ukeyElements[count] != key_[i] || elements[count] != &(elenode_[i]))
                mismatch++;
            count++;
        }
    }
#else
    for(int i = 0; i < bucket.size(); i++)
    {
        for(int j = 0; j < bucket[i].ndx.size(); j++)
        {
            if(ukeyElements[count] != bucket[i].key[j] || elements[count] != &(elenode_[bucket[i].ndx[j]]))
                mismatch++;
            count++;
        }
    }
#endif
    if(mismatch > 0)
    {
        printf("%s WARNING: AllEntriesPointers are out-dated. %d values pointers/keys do not match.\n", prefix,
               mismatch);
    }
    return mismatch;
}
void ElementsHashTable::updatePointersToNeighbours()
{
    for(int i = 0; i < elenode_.size(); i++)
    {
        if(status_[i]>=0)
        {
            if(elenode_[i].adapted_flag() > 0)
            {
                elenode_[i].update_neighbors_nodes_and_elements_pointers(this, NodeTable);
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

Element* ElementsHashTable::generateAddElement(const SFC_Key& keyi)
{
    Element* elm=add(keyi);
    elm->init(keyi);
    return elm;
}
    
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC* b, int mat,
                                            int* elm_loc_in, double pile_height, int myid, const SFC_Key& opposite_brother)
{
    Element* elm=add(nodekeys[8]); //--using bubble key to represent the element
    elm->init(nodekeys, neigh, n_pro, b, mat, elm_loc_in, pile_height, myid, opposite_brother);
    return elm;    
}
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen,
                 int elm_loc_in[], int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
                 ElementsHashTable *El_Table, NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
                 double Awetfather, double *drypoint_in)
{
    
    Element* elm=add(nodekeys[8]); //--using bubble key to represent the element
    elm->init(nodekeys, neigh, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    return elm;
}
Element* ElementsHashTable::generateAddElement(Element* sons[], NodeHashTable* NodeTable, ElementsHashTable* El_Table, MatProps* matprops_ptr)
{
    Element* elm=add(sons[2]->node_key(0));
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

