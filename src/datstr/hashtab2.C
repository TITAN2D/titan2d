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
#include "../header/hashtab.h"
#include "../header/element2.h"
#include "../header/node.h"
/*#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR*/
#include <mpi.h>
#include <limits.h>

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

ElementsHashTable *elementsHashTable;
#include <algorithm>

HashTableBase::HashTableBase(double *doublekeyrangein, int size, double XR[], double YR[])
{
    int i;

    NBUCKETS = size;
    
    for(i = 0; i < KEYLENGTH; i++)
        doublekeyrange[i] = doublekeyrangein[i];
    
    hashconstant = 8.0 * NBUCKETS / (doublekeyrange[0] * doublekeyrange[1] + doublekeyrange[1]);
    
    bucket = new HashEntryPtr[NBUCKETS];
    
    for(i = 0; i < NBUCKETS; i++)
        *(bucket + i) = 0;
    
    for(i = 0; i < 2; i++)
    {
        Xrange[i] = XR[i];
        Yrange[i] = YR[i];
    }
    
    invdxrange = 1.0 / (Xrange[1] - Xrange[0]);
    invdyrange = 1.0 / (Yrange[1] - Yrange[0]);
    
    ukeyBucket = new vector<uint64_t> [NBUCKETS];
    hashEntryBucket = new vector<HashEntry*> [NBUCKETS];
    ENTRIES = 0;
}

HashTableBase::~HashTableBase()              //evacuate the table
{
    for(int i = 0; i < NBUCKETS; i++)
    {
        HashEntryPtr p = *(bucket + i);
        while (p)
        {
            HashEntryPtr p_next = p->next;
            delete p;
            p = p_next;
        }
    }
    delete[] bucket;
    delete[] ukeyBucket;
    delete[] hashEntryBucket;
    ENTRIES = 0;
}

HashEntryPtr HashTableBase::searchBucket(HashEntryPtr p, const SFC_Key& keyi)
{
    int i;
    uint64_t ukey = get_ukey_from_sfc_key(keyi);
    while (p)
    {
        if(p->ukey == ukey)
            return p;
        p = p->next;
    }
    return NULL;
}

HashEntryPtr HashTableBase::addElement(int entry, const SFC_Key& key)
{
    
    HashEntryPtr p = new HashEntry(key);
    int i = 0;
    
    if(*(bucket + entry)) //this place is already occupied
    {
        HashEntryPtr currentPtr = *(bucket + entry);
        /*while (currentPtr != 0 && (key[0] > currentPtr->key[0]))
        {
            p->pre = currentPtr;
            currentPtr = currentPtr->next;
            i++;
        }
        
        if(currentPtr != 0 && key[0] == currentPtr->key[0])
        {
            while (currentPtr != 0 && (key[1] > currentPtr->key[1]))
            {
                p->pre = currentPtr;
                currentPtr = currentPtr->next;
                i++;
            }
            
        }*/
        while (currentPtr != nullptr && (key > currentPtr->key))
        {
            p->pre = currentPtr;
            currentPtr = currentPtr->next;
            i++;
        }
        

        if(currentPtr)
            currentPtr->pre = p;
        p->next = currentPtr;
        currentPtr = p->pre;
        //assert(currentPtr);
        if(currentPtr)
            currentPtr->next = p;
        else
            bucket[entry] = p;
        
        ukeyBucket[entry].insert(ukeyBucket[entry].begin() + i, p->ukey);
        hashEntryBucket[entry].insert(hashEntryBucket[entry].begin() + i, p);
    }
    
    //  p->next = *(bucket+entry);        //add the bucket to the head
    else
    {
        bucket[entry] = p;                //else eliminate it
        ukeyBucket[entry].push_back(p->ukey);
        hashEntryBucket[entry].push_back(p);
    }
    ENTRIES++;
    return p;
}

void*
HashTableBase::lookup(const SFC_Key &key)
{
    int entry = hash(key);
    
    uint64_t ukey = get_ukey_from_sfc_key(key);
    int size = ukeyBucket[entry].size();
    uint64_t *ukeyArr = &(ukeyBucket[entry][0]);
    int i;
    
    if(size == 0)
        return NULL;
    if(ukey < ukeyArr[0])
        return NULL;
    if(ukey > ukeyArr[size - 1])
        return NULL;
    
    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(ukey == ukeyArr[i])
            {
                return hashEntryBucket[entry][i]->value;
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
            if(ukey > ukeyArr[i1])
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
            if(ukey == ukeyArr[i])
            {
                return hashEntryBucket[entry][i]->value;
            }
        }
    }
    return NULL;
}

void HashTableBase::add(const SFC_Key& key, void* value)
{
    void* v = lookup(key);
    if(v == NULL)
    {
        int entry = hash(key);
        HashEntryPtr p = addElement(entry, key);
        p->value = value;
    }
    return;
}

void HashTableBase::remove(const SFC_Key& key)
{
    /* if(key[0] == (unsigned) 270752286)
     printf(" removing an unknown object with key %u %u\n",key[0], key[1]);*/
    int entry = hash(key);
    
    HashEntryPtr p = searchBucket(*(bucket + entry), key);
    
    if(!p)
        return;
    
    if(*(bucket + entry) == p)
    {
        *(bucket + entry) = p->next;
        delete p;
        ukeyBucket[entry].erase(ukeyBucket[entry].begin());
        hashEntryBucket[entry].erase(hashEntryBucket[entry].begin());
    }
    
    else
    {
        int i;
        for(i = 0; i < ukeyBucket[entry].size(); ++i)
        {
            if(ukeyBucket[entry][i] == p->ukey)
            {
                ukeyBucket[entry].erase(ukeyBucket[entry].begin() + i);
                hashEntryBucket[entry].erase(hashEntryBucket[entry].begin() + i);
                break;
            }
        }
        if(!(p->next))
            delete p;
        else
        {
            (p->pre)->next = p->next;
            (p->next)->pre = p->pre;
            delete p;
        }
    }
    ENTRIES--;
}
void HashTableBase::print0()
{
    bool minmax_init=false;
    uint64_t ukey_min;
    uint64_t ukey_max;

    for(int i = 0; i < NBUCKETS; i++)
    {
        int NEntriesInBucket = ukeyBucket[i].size();
        for(int j = 0; j < NEntriesInBucket; j++)
        {
            HashEntry* ent = hashEntryBucket[i][j];
            if(minmax_init)
            {
                if(ent->ukey<ukey_min)ukey_min=ent->ukey;
                if(ent->ukey>ukey_max)ukey_max=ent->ukey;
            }
            else
            {
                ukey_min=ent->ukey;
                ukey_max=ent->ukey;
                minmax_init=true;
            }
        }
    }
    return;
}
/*
 int HashTableBase:: hash(unsigned* key)
 {
 unsigned numer = *key;
 float    igaze = (float)numer/Range;
 int      igazee = igaze*NBUCKETS+.5;

 if(igazee>=NBUCKETS) igazee = NBUCKETS-1;

 return (igazee);
 }
 */
////////////////////////////////////////////////////////////////////////////////
template <typename T>
HashTable<T>::HashTable(double *doublekeyrangein, int size, double XR[], double YR[])
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
    
}
template <typename T>
ti_ndx_t HashTable<T>::lookup_ndx(const SFC_Key& key)
{
    int entry = hash(key);
    int size = bucket[entry].key.size();
    
    SFC_Key *keyArr = &(bucket[entry].key[0]);
    int i;
    
    if(size == 0)
        return ti_ndx_doesnt_exist;
    if(key < keyArr[0])
        return ti_ndx_doesnt_exist;
    if(key > keyArr[size - 1])
        return ti_ndx_doesnt_exist;
    
    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(key == keyArr[i])
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
            if(key > keyArr[i1])
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
            if(key == keyArr[i])
            {
                return bucket[entry].ndx[i];
            }
        }
    }
    return ti_ndx_doesnt_exist;
}
template <typename T>
T* HashTable<T>::lookup(const SFC_Key& key)
{
    ti_ndx_t ndx=lookup_ndx(key);
    if(status_[ndx]==CS_Removed)
    {
        assert(0);
        return nullptr;
    }
    if(ti_ndx_not_negative(ndx))
        return &(elenode_[ndx]);
    else
        return nullptr;
}
template <typename T>
ti_ndx_t HashTable<T>::add_ndx(const SFC_Key& key)
{
    ti_ndx_t ndx=lookup_ndx(key);
    if(ti_ndx_not_negative(ndx))
        return ndx;
    
    int entry = hash(key);
    int entry_size = bucket[entry].key.size();

    //get space
    ndx=key_.push_back();
    elenode_.push_back();
    status_.push_back();
    
    //set values
    key_[ndx]=key;
    status_[ndx]=CS_Added;

    //place to hash table
    if(entry_size>0)
    {
        //this place is already occupied
        //find proper place to insert it
        int i;
        SFC_Key *keyArr = &(bucket[entry].key[0]);
        for(i=0;i<entry_size&&key>keyArr[i];++i){}

        bucket[entry].key.insert(bucket[entry].key.begin() + i, key);
        bucket[entry].ndx.insert(bucket[entry].ndx.begin() + i, ndx);
    }
    else
    {
        //will be first member of the bucket entry
        bucket[entry].key.push_back(key);
        bucket[entry].ndx.push_back(ndx);
    }
    ENTRIES+=1;
    return ndx;
}
template <typename T>
T* HashTable<T>::add(const SFC_Key& key)
{
    ti_ndx_t ndx=add_ndx(key);
    return &(elenode_[ndx]);
}
template <typename T>
void HashTable<T>::remove(const SFC_Key& key)
{
    int entry = hash(key);
    int entry_size = bucket[entry].key.size();
    ti_ndx_t bucket_entry_ndx=bucket[entry].lookup_local_ndx(key);
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
    :HashTable<Node>(doublekeyrangein, size, XR, YR)
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
    
}


Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, matprops_ptr);
    return node;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, inf, ord, matprops_ptr);
    return node;
}
Node* NodeHashTable::createAddNode(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada)
{
    Node* node = add(keyi);
    node->init(keyi, coordi, inf, ord, elev, yada);
    return node;
}

Node* NodeHashTable::createAddNode(FILE* fp, MatProps* matprops_ptr) //for restart
{
    //unimplemented
    assert(0);
    return nullptr;
    //Node* node = add(keyi);
    //node->init(node->key(), node);
    //return node;
}
void NodeHashTable::removeNode(Node* node)
{
    remove(node->key());
    //delete node;
}




ElementsHashTable::ElementsHashTable(double *doublekeyrangein, int size, double XR[], double YR[], NodeHashTable* nodeTable) :
        HashTableBase(doublekeyrangein, size, XR, YR)
{
    NlocalElements = 0;
    NodeTable = nodeTable;
    elementsHashTable=this;
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
    for(int i = 0; i < NBUCKETS; i++)
    {
        int NEntriesInBucket = ukeyBucket[i].size();
        for(int j = 0; j < NEntriesInBucket; j++)
        {
            Element* Curr_El = (Element*) hashEntryBucket[i][j]->value;
            if(Curr_El->adapted_flag() > 0)
            { //if this element does not belong on this processor don't involve!!!
                ukeyLocalElements[count] = ukeyBucket[i][j];
                localElements[count] = (Element*) hashEntryBucket[i][j]->value;
                count++;
            }
        }
    }
    NlocalElements = count;
    ukeyLocalElements.resize(NlocalElements);
    localElements.resize(NlocalElements);
}
int ElementsHashTable::ckeckLocalElementsPointers(const char *prefix)
{
    int count = 0, NEntriesInBucket, mismatch = 0;
    for(int i = 0; i < NBUCKETS; i++)
    {
        NEntriesInBucket = ukeyBucket[i].size();
        for(int j = 0; j < NEntriesInBucket; j++)
        {
            Element* Curr_El = (Element*) hashEntryBucket[i][j]->value;
            if(Curr_El->adapted_flag() > 0)
            { //if this element does not belong on this processor don't involve!!!
                if(ukeyLocalElements[count] != ukeyBucket[i][j] || localElements[count] != hashEntryBucket[i][j]->value)
                    mismatch++;
                count++;
            }
        }
    }
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
    for(int i = 0; i < NBUCKETS; i++)
    {
        NEntriesInBucket = ukeyBucket[i].size();
        for(int j = 0; j < NEntriesInBucket; j++)
        {
            ukeyElements[count] = ukeyBucket[i][j];
            elements[count] = (Element*) hashEntryBucket[i][j]->value;
            count++;
        }
    }
}
int ElementsHashTable::ckeckElementsPointers(const char *prefix)
{
    int i, j, count = 0, NEntriesInBucket, mismatch = 0;
    if(ENTRIES != ukeyElements.size())
    {
        printf("%s WARNING: AllEntriesPointers are out-dated, number of entries do not match.\n", prefix);
        return ukeyElements.size();
    }
    for(int i = 0; i < NBUCKETS; i++)
    {
        NEntriesInBucket = ukeyBucket[i].size();
        for(j = 0; j < NEntriesInBucket; j++)
        {
            if(ukeyElements[count] != ukeyBucket[i][j] || elements[count] != hashEntryBucket[i][j]->value)
                mismatch++;
            count++;
        }
    }
    if(mismatch > 0)
    {
        printf("%s WARNING: AllEntriesPointers are out-dated. %d values pointers/keys do not match.\n", prefix,
               mismatch);
    }
    return mismatch;
}
void ElementsHashTable::updatePointersToNeighbours()
{
    int i;
    HashEntryPtr currentPtr;
    Element* Curr_El;
    for(i = 0; i < NBUCKETS; i++)
    {
        if(*(bucket + i))
        {
            
            currentPtr = *(bucket + i);
            while (currentPtr)
            {
                Curr_El = (Element*) (currentPtr->value);
                if(Curr_El->adapted_flag() > 0) //if this element does not belong on this processor don't involve!!!
                    Curr_El->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
                currentPtr = currentPtr->next;
            }
        }
    }
    return;
}
int ElementsHashTable::checkPointersToNeighbours(const char *prefix)
{
    int i;
    int count = 0;
    HashEntryPtr currentPtr;
    Element* Curr_El;
    for(i = 0; i < NBUCKETS; i++)
    {
        if(*(bucket + i))
        {
            
            currentPtr = *(bucket + i);
            while (currentPtr)
            {
                Curr_El = (Element*) (currentPtr->value);
                if(Curr_El->adapted_flag() > 0) //if this element does not belong on this processor don't involve!!!
                    count += Curr_El->check_neighbors_nodes_and_elements_pointers(this, NodeTable);
                currentPtr = currentPtr->next;
            }
        }
    }
    if(count > 0)
        printf("%s WARNING: neighbors nodes and elements pointers mismatch to key. %d mismatched.\n", prefix, count);
    /*int i;
     int pointersUpToDate = El_Table->isSortedBucketsUpToDate();
     int nodesPointersUpToDate = NodeTable->isSortedBucketsUpToDate();
     int NElements = El_Table->getNumberOfSortedBuckets();
     int Nnodes = NodeTable->getNumberOfSortedBuckets();
     Element** ElArr = (Element**) El_Table->getSortedBuckets();

     if (pointersUpToDate == 0 || nodesPointersUpToDate==0) {
     for (i = 0; i < NElements; i++) {
     ElArr[i]->update_neighbors_pointers(El_Table, NodeTable);
     }
     }*/
    return count;
}

Element* ElementsHashTable::generateAddElement(const SFC_Key& key)
{
    Element* elm=new Element(key);
    add(key, elm);
    return elm;
    
    /*tipos_t pos=elements_.push_back();
    elements_[pos].init(key);
    add(key, elements_.array_+pos);
    return elements_.array_+pos;*/
}
    
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC* b, int mat,
                                            int* elm_loc_in, double pile_height, int myid, const SFC_Key& opposite_brother)
{
    Element* elm=new Element(nodekeys, neigh, n_pro, b, mat, elm_loc_in, pile_height, myid, opposite_brother);
    add(elm->key(), elm);
    return elm;
    
    /*tipos_t pos=elements_.push_back();
    elements_[pos].init(nodekeys, neigh, n_pro, b, mat, elm_loc_in, pile_height, myid, opposite_brother);
    add(elements_[pos].key(), elements_.array_+pos);
    return elements_.array_+pos;*/
    
}
Element* ElementsHashTable::generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen,
                 int elm_loc_in[], int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
                 ElementsHashTable *El_Table, NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
                 double Awetfather, double *drypoint_in)
{
    Element* elm=new Element(nodekeys, neigh, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    add(elm->key(), elm);
    return elm;
    
    /*tipos_t pos=elements_.push_back();
    elements_[pos].init(nodekeys, neigh, n_pro, b, gen,
                                  elm_loc_in, ord, gen_neigh, mat, fthTemp, coord_in,
                                  El_Table, NodeTable, myid, matprops_ptr, iwetnodefather,
                                  Awetfather, drypoint_in);
    add(elements_[pos].key(), elements_.array_+pos);
    return elements_.array_+pos;*/
}
Element* ElementsHashTable::generateAddElement(Element* sons[], NodeHashTable* NodeTable, ElementsHashTable* El_Table, MatProps* matprops_ptr)
{
    Element* elm=new Element(sons, NodeTable, El_Table, matprops_ptr);
    add(elm->key(), elm);
    return elm;
    
    /*tipos_t pos=elements_.push_back();
    elements_[pos].init(sons, NodeTable, El_Table, matprops_ptr);
    add(elements_[pos].key(), elements_.array_+pos);
    return elements_.array_+pos;*/
}
Element* ElementsHashTable::generateAddElement(FILE* fp, NodeHashTable* NodeTable, MatProps* matprops_ptr, int myid)
{
    Element* elm=new Element(fp, NodeTable, matprops_ptr, myid);
    add(elm->key(), elm);
    return elm;
    
    /*tipos_t pos=elements_.push_back();
    elements_[pos].init(fp, NodeTable, matprops_ptr, myid);
    add(elements_[pos].key(), elements_.array_+pos);
    return elements_.array_+pos;*/
}

void ElementsHashTable::removeElement(Element* elm)
{
    HashTableBase::remove(elm->key());
    delete elm;
}

void* ElementsHashTable::lookup(const SFC_Key &key)
{
    int entry = hash(key);
    
    uint64_t ukey = get_ukey_from_sfc_key(key);
    int size = ukeyBucket[entry].size();
    uint64_t *ukeyArr = &(ukeyBucket[entry][0]);
    int i;
    
    if(size == 0)
        return NULL;
    if(ukey < ukeyArr[0])
        return NULL;
    if(ukey > ukeyArr[size - 1])
        return NULL;
    
    if(size < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < size; i++)
        {
            if(ukey == ukeyArr[i])
            {
                return hashEntryBucket[entry][i]->value;
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
            if(ukey > ukeyArr[i1])
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
            if(ukey == ukeyArr[i])
            {
                return hashEntryBucket[entry][i]->value;
            }
        }
    }
    return NULL;
}

/*void ElementsHashTable::add(unsigned* key, void* value) {
 int i;
 HashTable::add(key, value);
 Element* Curr_El;
 Curr_El = (Element*) value;
 //if (Curr_El->get_adapted_flag() > 0)  //if this element does not belong on this processor don't involve!!!
 //assuming neigbours keys are already good
 Curr_El->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
 for(i=0;i<8;i++){
 if(Curr_El->getNeighborPtr(i)!=NULL)
 Curr_El->getNeighborPtr(i)->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
 }
 }

 void ElementsHashTable::remove(unsigned* key) {
 int i;
 unsigned* neighbors;

 Element* Curr_El,*neighborElement[8];

 Curr_El=(Element*)lookup(key);
 neighbors=Curr_El->get_neighbors();

 for(i=0;i<8;i++){
 neighborElement[i]=(Element*)lookup(neighbors+i*KEYLENGTH);
 }
 HashTable::remove(key);
 for(i=0;i<8;i++){
 if(neighborElement[i]!=NULL)
 neighborElement[i]->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
 }
 }
 // for debugging...
 void ElementsHashTable::remove(unsigned* key, int whatflag) {
 int i;
 unsigned* neighbors;

 Element* Curr_El,*neighborElement[8];
 Curr_El=(Element*)lookup(key);
 neighbors=Curr_El->get_neighbors();

 for(i=0;i<8;i++){
 neighborElement[i]=(Element*)lookup(neighbors+i*KEYLENGTH);
 }
 HashTable::remove(key,whatflag);
 for(i=0;i<8;i++){
 if(neighborElement[i]!=NULL)
 neighborElement[i]->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
 }
 }
 void ElementsHashTable::remove(unsigned* key, int whatflag, FILE *fp, int myid, int where){
 int i;
 unsigned* neighbors;

 Element* Curr_El,*neighborElement[8];
 Curr_El=(Element*)lookup(key);

 neighbors=Curr_El->get_neighbors();

 for(i=0;i<8;i++){
 neighborElement[i]=(Element*)lookup(neighbors+i*KEYLENGTH);
 }
 HashTable::remove(key,whatflag, fp, myid, where);
 for(i=0;i<8;i++){
 if(neighborElement[i]!=NULL)
 neighborElement[i]->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
 }
 }*/
