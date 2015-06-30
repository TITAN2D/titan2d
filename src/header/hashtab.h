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

class Element;

struct HashEntry
{
    unsigned key[KEYLENGTH];  //key: object key word
    uint64_t ukey; //ukey: single 64 bit key
    void* value;   //value: poiter to record
    HashEntry* pre;
    HashEntry* next;    //pre, next: objects with same entry will be stored in a two-way link
    
    HashEntry(unsigned* keyi)
    {
        int i;
        for(i = 0; i < KEYLENGTH; i++)
            key[i] = *(keyi + i);
        ukey = (((uint64_t) keyi[0]) << 32) | keyi[1];
        next = NULL;
        pre = NULL;
    }
    
    /*HashEntry() {
     value = NULL;
     next = NULL;
     pre = NULL;
     }*/

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
    unsigned MinKey[2];
    unsigned MaxKey[2];
    unsigned Range;
    double doublekeyrange[2];
    double hashconstant;
    double Xrange[2];
    double Yrange[2];
    double invdxrange;
    double invdyrange;

    HashEntryPtr* bucket;
    int NBUCKETS;
    int PRIME;
    int ENTRIES;

    vector<uint64_t> *ukeyBucket;
    vector<HashEntry*> *hashEntryBucket;

    HashEntryPtr addElement(int entry, unsigned* key);
    HashEntryPtr searchBucket(HashEntryPtr p, unsigned* key);

public:
    HashTable(unsigned*, unsigned*, int, int);
    HashTable(double *doublekeyrangein, int, int, double* XR, double* YR, int ifrestart);
    virtual ~HashTable();

    int hash(unsigned* key);
    void add(unsigned* key, void* value);
    void* lookup(unsigned* key);
    void remove(unsigned* key);
    void remove(unsigned* key, int whatflag);  //for debugging
    void remove(unsigned* key, int whatflag, FILE *fp, int myid, int where);  //for debugging
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
    unsigned* get_MinKey();
    unsigned* get_MaxKey();
    int get_nbuckets();

    /*
     double* getXrange();
     double* getYrange();
     */
    int get_no_of_entries();
};

inline double* HashTable::get_doublekeyrange()
{
    return doublekeyrange;
}
;

inline HashEntryPtr* HashTable::getbucketptr()
{
    return bucket;
}

inline int HashTable::get_no_of_buckets()
{
    return NBUCKETS;
}
;

inline double HashTable::get_invdxrange()
{
    return invdxrange;
}
;

inline double HashTable::get_invdyrange()
{
    return invdyrange;
}
;

inline double* HashTable::get_Xrange()
{
    return Xrange;
}
;

inline double* HashTable::get_Yrange()
{
    return Yrange;
}
;

inline unsigned* HashTable::get_MinKey()
{
    return MinKey;
}
;

inline unsigned* HashTable::get_MaxKey()
{
    return MaxKey;
}
;

inline int HashTable::get_nbuckets()
{
    return NBUCKETS;
}
;
inline int HashTable::get_no_of_entries()
{
    return ENTRIES;
}
;

inline int HashTable::hash(unsigned* key)
{
    //Keith made this change 20061109; and made hash an inline function
    /* NBUCKETS*2 is NBUCKETS*integer integer is empirical could be 1
     return (((int) ((key[0]*doublekeyrange[1]+key[1])/
     (doublekeyrange[0]*doublekeyrange[1]+doublekeyrange[1])*
     NBUCKETS*2+0.5) )%NBUCKETS);
     */
    return (((int) ((key[0] * doublekeyrange[1] + key[1]) * hashconstant + 0.5)) % NBUCKETS);
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
    ElementsHashTable(unsigned*, unsigned*, int, int, HashTable*);
    ElementsHashTable(double *doublekeyrangein, int, int, double* XR, double* YR, int ifrestart, HashTable* nodeTable);
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
    
};

#endif
