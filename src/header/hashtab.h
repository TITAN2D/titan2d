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
 * $Id: hashtab.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

/* Hash table */
/* Every table can process NBUCKETS entries */
/* All buckets with the same entry are linked together*/
/* The pointer to the first bucket is stored in the table */
/*---------------------------------------------------------*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "constant.h"
#include "sfc.h"
#include "tivector.h"

class Element;
class BC;
class MatProps;

struct HashEntry
{
    SFC_Key key;  //key: object key word
    uint64_t ukey; //ukey: single 64 bit key
    void* value;   //value: poiter to record
    HashEntry* pre;
    HashEntry* next;    //pre, next: objects with same entry will be stored in a two-way link
    
    HashEntry(const SFC_Key& keyi)
    {
        int i;
        key=keyi;
        ukey = get_ukey_from_sfc_key(key);
        value = NULL;
        pre = NULL;
        next = NULL;
    }

    ~HashEntry()
    {         //keep the follower when deleting an object
        if(next)
            next->pre = pre;
        if(pre)
            pre->next = next;
        
    }
    
};

typedef HashEntry* HashEntryPtr;

//! Hashtables store pointers to each Element or Node (of which HashTable is a friend class), these pointers can be accessed by giving the hashtable the "key" of the element number you want to "lookup."  The keys are ordered sequentially by a space filling curve that ensures that the pointers to elements (or nodes) that are located close to each other in physical space will usually be located close to each other in memory, which speeds up access time.  Each key is a single number that spans several unsigned variables (elements of an array).
class HashTable
{
    
    //friend int hash(unsigned* keyi);
    friend class Element;

protected:
    double doublekeyrange[2];
    double hashconstant;
    double Xrange[2];
    double Yrange[2];
    double invdxrange;
    double invdyrange;

    HashEntryPtr* bucket;
    int NBUCKETS;
    int ENTRIES;

    vector<uint64_t> *ukeyBucket;
    vector<HashEntry*> *hashEntryBucket;

    HashEntryPtr addElement(int entry, const SFC_Key& key);
    HashEntryPtr searchBucket(HashEntryPtr p, const SFC_Key& key);

public:
    HashTable(double *doublekeyrangein, int size, double XR[], double YR[]);
    virtual ~HashTable();

    int hash(const SFC_Key& key) const;
    void add(const SFC_Key& key, void* value);
    void* lookup(const SFC_Key& key);
    void remove(const SFC_Key& key);
    void print_out(int);
    int get_no_of_buckets();
//    void   get_element_stiffness(HashTable*);
    HashEntryPtr* getbucketptr();
    void* get_value();

    double* get_Xrange();
    double* get_Yrange();
    double* get_doublekeyrange();
    double get_invdxrange();
    double get_invdyrange();
    int get_nbuckets();
    int get_no_of_entries();

    void print0();
};

inline double* HashTable::get_doublekeyrange()
{
    return doublekeyrange;
}

inline HashEntryPtr* HashTable::getbucketptr()
{
    return bucket;
}

inline int HashTable::get_no_of_buckets()
{
    return NBUCKETS;
}

inline double HashTable::get_invdxrange()
{
    return invdxrange;
}

inline double HashTable::get_invdyrange()
{
    return invdyrange;
}

inline double* HashTable::get_Xrange()
{
    return Xrange;
}

inline double* HashTable::get_Yrange()
{
    return Yrange;
}


inline int HashTable::get_nbuckets()
{
    return NBUCKETS;
}

inline int HashTable::get_no_of_entries()
{
    return ENTRIES;
}

inline int HashTable::hash(const SFC_Key& key) const
{
    //Keith made this change 20061109; and made hash an inline function
    /* NBUCKETS*2 is NBUCKETS*integer integer is empirical could be 1
     return (((int) ((key[0]*doublekeyrange[1]+key[1])/
     (doublekeyrange[0]*doublekeyrange[1]+doublekeyrange[1])*
     NBUCKETS*2+0.5) )%NBUCKETS);
     */
    return (((int) ((key.key[0] * doublekeyrange[1] + key.key[1]) * hashconstant + 0.5)) % NBUCKETS);
}
//! Hashtables for Elements
class ElementsHashTable: public HashTable
{
    
    //friend int hash(unsigned* keyi);
    friend class Element;
protected:
    HashTable* NodeTable;

    vector<uint64_t> ukeyElements;
    vector<Element*> elements;

    int NlocalElements;
    vector<uint64_t> ukeyLocalElements;
    vector<Element*> localElements;
public:
    ElementsHashTable(double *doublekeyrangein, int size, double XR[], double YR[], HashTable* nodeTable);
    virtual ~ElementsHashTable();

    Element** getElements()
    {
        return &(elements[0]);
    }
    
    /*	virtual void add(unsigned* key, void* value);
     virtual void remove(unsigned* key);
     virtual void remove(unsigned* key, int whatflag);  //for debugging
     virtual void remove(unsigned* key, int whatflag, FILE *fp, int myid, int where);  //for debugging
     */
    void updateElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckElementsPointers(const char *prefix);

    int getNumberOfLocalElements()
    {
        return NlocalElements;
    }
    Element** getLocalElementsValues()
    {
        return &(localElements[0]);
    }
    //only works for elements
    void updateLocalElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckLocalElementsPointers(const char *prefix);

    //! set neighboring elements and nodes pointers in elements
    void updatePointersToNeighbours();

    //! check that neighboring elements and nodes pointers are correct
    int checkPointersToNeighbours(const char *prefix);
    
    //! default Element generator constructor, does nothing except set stoppedflags=2, this should never be used
    virtual Element* generateElement();

    //! constructor that creates an original element when funky is read in
    virtual Element* generateElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int mat, int *elm_loc_in,
                             double pile_height, int myid, const SFC_Key& opposite_brother);

    //! constructor that creates a son element from its father during refinement
    virtual Element* generateElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen, int elm_loc_in[],
                int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in, HashTable *El_Table,
                HashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather, double Awetfather,
                double *drypoint_in);
    //! constructor that creates a father element from its four sons during unrefinement
    virtual Element* generateElement(Element *sons[], HashTable *NodeTable, HashTable *El_Table, MatProps *matprops_ptr);

    //! constructor that creates/restores a saved element during restart
    virtual Element* generateElement(FILE* fp, HashTable* NodeTable, MatProps* matprops_ptr, int myid);
    
    
    //here goes element content storage probably should be separate class at the end
    
    tivector<SFC_Key> key_;

};

#endif
