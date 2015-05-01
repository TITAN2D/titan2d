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
#include "../header/hashtab.h"
#include "../header/element2.h"
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include <limits.h>

#define HASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

HashTable::HashTable(unsigned* min, unsigned* max, int size, int prime) {
	assert(0);  //avoid using this since it doesn't intialize doublekeyrange

	MinKey[0] = *min;
	MinKey[1] = *(min + 1);
	MaxKey[0] = *max;
	MaxKey[1] = *(max + 1);
	// extend the hashtable bounds a little bit to make it more efficient for adaptive meshes
	unsigned hashtable_extender = HASHTABLE_EXTENDER;
	if (MinKey[0] >= hashtable_extender)
		MinKey[0] -= hashtable_extender;
	else
		MinKey[0] = 0;
	unsigned umax = IScale;
	if ((hashtable_extender / 2 + MaxKey[0] / 2) <= (umax / 2))
		MaxKey[0] += hashtable_extender;
	else
		MaxKey[0] = umax;

	NBUCKETS = size;
	PRIME = prime;
	//Range  = *(MaxKey);
	Range = *(MaxKey) - *(MinKey); //Keith Made this change 20061109

	bucket = new HashEntryPtr[NBUCKETS];

	for (int i = 0; i < NBUCKETS; i++)
		*(bucket + i) = 0;

	/*  MaxMinX[0]=MaxMinY[0]=-1;
	 
	 MaxMinX[1]=MaxMinY[1]=1;*/

	ukeyBucket=new vector<uint64_t>[NBUCKETS];
	hashEntryBucket=new vector<HashEntry*> [NBUCKETS];

	ENTRIES=0;
}

HashTable::HashTable(double *doublekeyrangein, int size, int prime, double XR[], double YR[],
    int ifrestart) {
	int i;
	/*
	 MinKey[0] = *min;
	 MinKey[1] = *(min+1);
	 MaxKey[0] = *max;
	 MaxKey[1] = *(max+1);
	 */

	NBUCKETS = size;
	PRIME = prime;

	for (i = 0; i < KEYLENGTH; i++)
		doublekeyrange[i] = doublekeyrangein[i];

	hashconstant = 8.0 * NBUCKETS / (doublekeyrange[0] * doublekeyrange[1] + doublekeyrange[1]);

	bucket = new HashEntryPtr[NBUCKETS];

	for (i = 0; i < NBUCKETS; i++)
		*(bucket + i) = 0;

	for (i = 0; i < 2; i++) {
		Xrange[i] = XR[i];
		Yrange[i] = YR[i];
	}

	invdxrange = 1.0 / (Xrange[1] - Xrange[0]);
	invdyrange = 1.0 / (Yrange[1] - Yrange[0]);

	ukeyBucket=new vector<uint64_t>[NBUCKETS];
	hashEntryBucket=new vector<HashEntry*> [NBUCKETS];
	ENTRIES=0;
}

HashTable::~HashTable()              //evacuate the table
{
	for (int i = 0; i < NBUCKETS; i++) {
		HashEntryPtr p = *(bucket + i);
		while (p) {
			HashEntryPtr p_next = p->next;
			delete p;
			p = p_next;
		}
	}
	delete[] bucket;
	delete[] ukeyBucket;
	delete[] hashEntryBucket;
	ENTRIES=0;
}

HashEntryPtr HashTable::searchBucket(HashEntryPtr p, unsigned* keyi) {
	int i;
	uint64_t ukey=(((uint64_t) keyi[0]) << 32) | keyi[1];
	while(p) {
		if(p->ukey==ukey)
			return p;
		p = p->next;
	}
	return NULL;
}

HashEntryPtr HashTable::addElement(int entry, unsigned key[]) {

	HashEntryPtr p = new HashEntry(key);
	int i=0;

	if (*(bucket + entry)) //this place is already occupied
	{
		HashEntryPtr currentPtr = *(bucket + entry);
		while (currentPtr != 0 && (key[0] > currentPtr->key[0])) {
			p->pre = currentPtr;
			currentPtr = currentPtr->next;
			i++;
		}

		if (currentPtr != 0 && key[0] == currentPtr->key[0]) {
			while (currentPtr != 0 && (key[1] > currentPtr->key[1])) {
				p->pre = currentPtr;
				currentPtr = currentPtr->next;
				i++;
			}

		}

		if (currentPtr)
			currentPtr->pre = p;
		p->next = currentPtr;
		currentPtr = p->pre;
		//assert(currentPtr);
		if (currentPtr)
			currentPtr->next = p;
		else
			bucket[entry] = p;

		ukeyBucket[entry].insert(ukeyBucket[entry].begin()+i,p->ukey);
		hashEntryBucket[entry].insert(hashEntryBucket[entry].begin()+i,p);
	}

	//  p->next = *(bucket+entry);        //add the bucket to the head
	else{
		bucket[entry] = p;                //else eliminate it
		ukeyBucket[entry].push_back(p->ukey);
		hashEntryBucket[entry].push_back(p);
	}
	ENTRIES++;
	return p;
}

#define HASHTABLE_LOOKUP_LINSEARCH 8
void*
HashTable::lookup(unsigned* key) {
	int entry = hash(key);

	uint64_t ukey = (((uint64_t) key[0]) << 32) | key[1];
	int size = ukeyBucket[entry].size();
	uint64_t *ukeyArr = &(ukeyBucket[entry][0]);
	int i;

	if (size == 0)
		return NULL;
	if (ukey < ukeyArr[0])
		return NULL;
	if (ukey > ukeyArr[size - 1])
		return NULL;

	if (size < HASHTABLE_LOOKUP_LINSEARCH) {
		for (i = 0; i < size; i++) {
			if (ukey == ukeyArr[i]) {
				return hashEntryBucket[entry][i]->value;
			}
		}
	} else {
		int i0, i1, i2;
		i0 = 0;
		i1 = size / 2;
		i2 = size - 1;
		while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH) {
			if (ukey > ukeyArr[i1]) {
				i0 = i1 + 1;
				i1 = (i0 + i2) / 2;
			} else {
				i2 = i1;
				i1 = (i0 + i2) / 2;
			}
		}
		for (i = i0; i <= i2; i++) {
			if (ukey == ukeyArr[i]) {
				return hashEntryBucket[entry][i]->value;
			}
		}
	}
	return NULL;
}

void HashTable::add(unsigned* key, void* value) {
	void* v=lookup(key);
	if (v == NULL) {
		int entry = hash(key);
		HashEntryPtr p = addElement(entry, key);
		p->value = value;
	}
	return;
}

void HashTable::remove(unsigned* key) {
	/* if(key[0] == (unsigned) 270752286)
	 printf(" removing an unknown object with key %u %u\n",key[0], key[1]);*/
	int entry = hash(key);

	HashEntryPtr p = searchBucket(*(bucket + entry), key);

	if (!p)
		return;

	if (*(bucket + entry) == p) {
		*(bucket + entry) = p->next;
		delete p;
		ukeyBucket[entry].erase(ukeyBucket[entry].begin());
		hashEntryBucket[entry].erase(hashEntryBucket[entry].begin());
	}

	else {
		int i;
		for (i = 0; i < ukeyBucket[entry].size(); ++i) {
			if (ukeyBucket[entry][i] == p->ukey) {
				ukeyBucket[entry].erase(ukeyBucket[entry].begin() + i);
				hashEntryBucket[entry].erase(hashEntryBucket[entry].begin() + i);
				break;
			}
		}
		if (!(p->next))
			delete p;
		else {
			(p->pre)->next = p->next;
			(p->next)->pre = p->pre;
			delete p;
		}
	}
	ENTRIES--;
}
// for debugging...
void HashTable::remove(unsigned* key, int whatflag) {
	/*  if(whatflag == 1) {
	 printf(" removing an element with key %u %u\n",key[0], key[1]);
	 } 
	 else {
	 printf(" removing a node with key %u %u\n",key[0], key[1]);
	 }*/
	int entry = hash(key);

	HashEntryPtr p = searchBucket(*(bucket + entry), key);

	if (!p) {
		printf(" this object has already been deleted, entry type=%d  ****************************\n",
		    whatflag);
		return;
	}

	if (*(bucket + entry) == p) {
		*(bucket + entry) = p->next;
		delete p;
		ukeyBucket[entry].erase(ukeyBucket[entry].begin());
		hashEntryBucket[entry].erase(hashEntryBucket[entry].begin());
	}

	else {
		int i;
		for (i = 0; i < ukeyBucket[entry].size(); ++i) {
			if (ukeyBucket[entry][i] == p->ukey) {
				ukeyBucket[entry].erase(ukeyBucket[entry].begin() + i);
				hashEntryBucket[entry].erase(hashEntryBucket[entry].begin() + i);
				break;
			}
		}
		if (!(p->next))
			delete p;
		else {
			(p->pre)->next = p->next;
			(p->next)->pre = p->pre;
			delete p;
		}
	}
	ENTRIES--;
}

void HashTable::remove(unsigned* key, int whatflag, FILE *fp, int myid, int where) {
	if (fp == NULL)
		fp = stdout;

	/*  if(whatflag == 1) {
	 printf(" removing an element with key %u %u\n",key[0], key[1]);
	 } 
	 else {
	 printf(" removing a node with key %u %u\n",key[0], key[1]);
	 }*/
	int entry = hash(key);

	HashEntryPtr p = searchBucket(*(bucket + entry), key);

	if (!p) {
		fprintf(fp,
		    " this object has already been deleted, myid=%d entry type=%d key={%u,%u} where=%d ****************************\n",
		    myid, whatflag, key[0], key[1], where);
		return;
	}

	if (*(bucket + entry) == p) {
		*(bucket + entry) = p->next;
		delete p;
		ukeyBucket[entry].erase(ukeyBucket[entry].begin());
		hashEntryBucket[entry].erase(hashEntryBucket[entry].begin());
	}

	else {
		int i;
		for (i = 0; i < ukeyBucket[entry].size(); ++i) {
			if (ukeyBucket[entry][i] == p->ukey) {
				ukeyBucket[entry].erase(ukeyBucket[entry].begin() + i);
				hashEntryBucket[entry].erase(hashEntryBucket[entry].begin() + i);
				break;
			}
		}
		if (!(p->next))
			delete p;
		else {
			(p->pre)->next = p->next;
			(p->next)->pre = p->pre;
			delete p;
		}
	}
	ENTRIES--;
}

/*
 int HashTable:: hash(unsigned* key)
 {
 unsigned numer = *key;
 float    igaze = (float)numer/Range;
 int      igazee = igaze*NBUCKETS+.5;

 if(igazee>=NBUCKETS) igazee = NBUCKETS-1;

 return (igazee);
 }
 */
ElementsHashTable::ElementsHashTable(unsigned* min, unsigned* max, int size, int prime,HashTable* nodeTable)
	: HashTable(min, max, size, prime)
{
	NlocalElements=0;
	NodeTable=nodeTable;
}

ElementsHashTable::ElementsHashTable(double *doublekeyrangein, int size, int prime, double XR[], double YR[],
    int ifrestart,HashTable* nodeTable)
	: HashTable(doublekeyrangein, size, prime, XR, YR, ifrestart)
{
	NlocalElements=0;
	NodeTable=nodeTable;
}

ElementsHashTable::~ElementsHashTable()              //evacuate the table
{

}
void ElementsHashTable::updateLocalElements(){
	int i,j,count=0,NEntriesInBucket;

	ukeyLocalElements.resize(ENTRIES);
	localElements.resize(ENTRIES);
	for(int i = 0; i < NBUCKETS; i++){
		NEntriesInBucket=ukeyBucket[i].size();
		for(j=0;j<NEntriesInBucket;j++){
			Element* Curr_El=(Element*)hashEntryBucket[i][j]->value;
			if (Curr_El->get_adapted_flag() > 0) { //if this element does not belong on this processor don't involve!!!
				ukeyLocalElements[count]=ukeyBucket[i][j];
				localElements[count]=(Element*)hashEntryBucket[i][j]->value;
				count++;
			}
		}
	}
	NlocalElements=count;
	ukeyLocalElements.resize(NlocalElements);
	localElements.resize(NlocalElements);
}
int ElementsHashTable::ckeckLocalElementsPointers(const char *prefix){
	int i,j,count=0,NEntriesInBucket,mismatch=0;
	for(int i = 0; i < NBUCKETS; i++){
		NEntriesInBucket=ukeyBucket[i].size();
		for(j=0;j<NEntriesInBucket;j++){
			Element* Curr_El=(Element*)hashEntryBucket[i][j]->value;
			if (Curr_El->get_adapted_flag() > 0) { //if this element does not belong on this processor don't involve!!!
				if(ukeyLocalElements[count]!=ukeyBucket[i][j] || localElements[count]!=hashEntryBucket[i][j]->value)
					mismatch++;
				count++;
			}
		}
	}
	if(mismatch>0){
		printf("%s WARNING: AllEntriesPointersLocal are out-dated. %d values pointers/keys do not match.\n",prefix,mismatch);
	}
	return mismatch;
}
void ElementsHashTable::updateElements(){
	int i,j,count=0,NEntriesInBucket;
	ukeyElements.resize(ENTRIES);
	elements.resize(ENTRIES);
	for(int i = 0; i < NBUCKETS; i++){
		NEntriesInBucket=ukeyBucket[i].size();
		for(j=0;j<NEntriesInBucket;j++){
			ukeyElements[count]=ukeyBucket[i][j];
			elements[count]=(Element*)hashEntryBucket[i][j]->value;
			count++;
		}
	}
}
int ElementsHashTable::ckeckElementsPointers(const char *prefix){
	int i,j,count=0,NEntriesInBucket,mismatch=0;
	if(ENTRIES!=ukeyElements.size()){
		printf("%s WARNING: AllEntriesPointers are out-dated, number of entries do not match.\n",prefix);
		return ukeyElements.size();
	}
	for(int i = 0; i < NBUCKETS; i++){
		NEntriesInBucket=ukeyBucket[i].size();
		for(j=0;j<NEntriesInBucket;j++){
			if(ukeyElements[count]!=ukeyBucket[i][j] || elements[count]!=hashEntryBucket[i][j]->value)
				mismatch++;
			count++;
		}
	}
	if(mismatch>0){
		printf("%s WARNING: AllEntriesPointers are out-dated. %d values pointers/keys do not match.\n",prefix,mismatch);
	}
	return mismatch;
}
void ElementsHashTable::updatePointersToNeighbours() {
	int i;
	HashEntryPtr currentPtr;
	Element* Curr_El;
	for (i = 0; i < NBUCKETS; i++){
		if (*(bucket + i)) {

			currentPtr = *(bucket + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)  //if this element does not belong on this processor don't involve!!!
					Curr_El->update_neighbors_nodes_and_elements_pointers(this, NodeTable);
				currentPtr = currentPtr->next;
			}
		}
	}
	return;
}
int ElementsHashTable::checkPointersToNeighbours(const char *prefix) {
	int i;
	int count=0;
	HashEntryPtr currentPtr;
	Element* Curr_El;
	for (i = 0; i < NBUCKETS; i++){
		if (*(bucket + i)) {

			currentPtr = *(bucket + i);
			while (currentPtr) {
				Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0)  //if this element does not belong on this processor don't involve!!!
					count+=Curr_El->check_neighbors_nodes_and_elements_pointers(this, NodeTable);
				currentPtr = currentPtr->next;
			}
		}
	}
	if(count>0)
		printf("%s WARNING: neighbors nodes and elements pointers mismatch to key. %d mismatched.\n",prefix,count);
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
