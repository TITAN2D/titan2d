/*******************************************************************
 * Copyright (C) 2015 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Description: ticore test
 *
 *******************************************************************
 */

#ifndef TESTS_TICORE_TICORE_TISORT_HPP_
#define TESTS_TICORE_TICORE_TISORT_HPP_

#include <assert.h>
#include <stdint.h>

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/header/sfc.h"
#include "../../src/header/tivector.h"
#include "../../src/header/titan2d.h"
#include "../../src/header/titan2d_utils.h"

#include "../../src/header/ticore/tisort.hpp"

//!build index list
void build_ndx(const vector<SFC_Key> &key, vector<ti_ndx_t> &ndx)
{
    ndx.resize(key.size());
    for(ti_ndx_t i=i;i<key.size();++i)
        ndx[i]=i;
}

//!internals for reference sort
struct KeyNdxPair {
    SFC_Key key_;
    ti_ndx_t ndx_;

    bool operator<( const KeyNdxPair& rhs ) const
        { return key_ < rhs.key_; }
};
//!reference sort
void sort_ref(vector<SFC_Key> &key, vector<ti_ndx_t> &ndx)
{
    vector<KeyNdxPair> v;
    v.resize(key.size());
    for(int i=0;i<key.size();++i)
    {
        v[i].key_=key[i];
        v[i].ndx_=ndx[i];
    }
    sort(v.begin(), v.end());
    for(int i=0;i<key.size();++i)
    {
        key[i]=v[i].key_;
        ndx[i]=v[i].ndx_;
    }
}

//!return true if vectors are equal
template<typename T>
bool vectors_are_equal(const vector<T> &A, const vector<T> &B)
{
    if(A.size()!=B.size())return false;
    for(int i=0;i<A.size();++i)
    {
        if(A[i]!=B[i])return false;
    }
    return true;
}

class UniqRandGen
{
public:
    std::mt19937_64 gen;
    std::set<SFC_Key> already_generated;
    std::uniform_int_distribution<SFC_Key> dis;

    UniqRandGen()
        :dis(0,UINT64_MAX),gen(2015)
    {
    }
    void reset()
    {
        already_generated.clear();
    }
    SFC_Key rand()
    {
        while(1)
        {
            SFC_Key number=dis(gen);
            if(already_generated.count(number)==0)
            {
                return number;
            }
        }
    }
};

void fill_with_randomnumbers(vector<SFC_Key> &key, const int N, const bool sort_vector, UniqRandGen &rand)
{
    key.resize(N);

    for(ti_ndx_t i=0;i<N;++i)
    {
        key[i]=rand.rand();
    }

    if(sort_vector)sort(key.begin(), key.end() );
}


void MergeSort_BottomUp(vector<SFC_Key> &key, vector<ti_ndx_t> &ndx)
{
    static vector<SFC_Key> key_work;
    static vector<ti_ndx_t> ndx_work;
    key_work.resize(key.size());
    ndx_work.resize(ndx.size());
    MergeSort_BottomUp(key.size(), &(key[0]), &(ndx[0]), &(key_work[0]), &(ndx_work[0]));
}


inline void __tisort_sortedarrays_merge_test(vector<SFC_Key> &A,vector<SFC_Key> &B,vector<SFC_Key> &C, vector<ti_ndx_t> &Cndx)
{
    vector<SFC_Key> AB;
    vector<ti_ndx_t> ABndx;
    vector<SFC_Key> ABref;
    vector<ti_ndx_t> ABndxref;
    ti_ndx_t Asize,Bsize,Csize;

    Asize=A.size(),Bsize=B.size();Csize=Asize+Bsize;
    AB=A;
    AB.insert(AB.end(), B.begin(), B.end());
    build_ndx(AB,ABndx);
    C.resize(Csize);Cndx.resize(Csize);
    ABref=AB;
    ABndxref=ABndx;

    MergeSortedArraysToSortedArray(Asize, &(AB[0]), &(ABndx[0]),
                                   Bsize, &(AB[Asize]), &(ABndx[Asize]),
                                   Csize, &(C[0]), &(Cndx[0]));

    sort_ref(ABref,ABndxref);

    /*for(ti_ndx_t i=0;i<ABref.size();++i)
    {
        cout<<i<<" "<<C[i]<<" "<<ABref[i]<<" "<<Cndx[i]<<" "<<ABndxref[i]<<"\n";
    }*/


    assert(vectors_are_equal(C,ABref));
    assert(vectors_are_equal(Cndx,ABndxref));
}
inline void __tisort_sortedarrays_merge_test_omp(vector<SFC_Key> &A,vector<SFC_Key> &B,vector<SFC_Key> &C, vector<ti_ndx_t> &Cndx)
{
    vector<SFC_Key> AB;

    vector<ti_ndx_t> ABndx;

    vector<SFC_Key> ABref;
    vector<ti_ndx_t> ABndxref;
    ti_ndx_t Asize,Bsize,Csize;

    Asize=A.size(),Bsize=B.size();Csize=Asize+Bsize;
    AB=A;
    AB.insert(AB.end(), B.begin(), B.end());
    build_ndx(AB,ABndx);
    C.resize(Csize);Cndx.resize(Csize);
    ABref=AB;
    ABndxref=ABndx;

    #pragma omp parallel
    {
        MergeSortedArraysToSortedArray_omp(Asize, &(AB[0]), &(ABndx[0]),
                                       Bsize, &(AB[Asize]), &(ABndx[Asize]),
                                       Csize, &(C[0]), &(Cndx[0]));
    }
    sort_ref(ABref,ABndxref);

    assert(vectors_are_equal(C,ABref));
    assert(vectors_are_equal(Cndx,ABndxref));
}
inline void tisort_sortedarrays_merge_test()
{
    vector<SFC_Key> C;
    vector<ti_ndx_t> Cndx;

    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short01\n");
        vector<SFC_Key> A={1};
        vector<SFC_Key> B={5};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short03\n");
        vector<SFC_Key> A={1};
        vector<SFC_Key> B;
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }

    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short1\n");
        vector<SFC_Key> A={1,7};
        vector<SFC_Key> B={5};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
        /*for(int i=0;i<C.size();++i)
        {
            cout<<i<<" "<<C[i]<<"\n";
        }*/
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short2_Bempty1\n");
        vector<SFC_Key> A={0,1,2,3,4,5,6,7,8,9,10};
        vector<SFC_Key> B={};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short3\n");
        vector<SFC_Key> A={0,1,2,4,5,8,10};
        vector<SFC_Key> B={3,6,7,9};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short4_left\n");
        vector<SFC_Key> A={2,3,4,6,7,8,9,10};
        vector<SFC_Key> B={0,1,5};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short5_left\n");
        vector<SFC_Key> A={2,3,4,6,7,8};
        vector<SFC_Key> B={5,9,10};
        __tisort_sortedarrays_merge_test(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short60\n");
        vector<SFC_Key> A={3};
        vector<SFC_Key> B={2};
        __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);

    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short61\n");
        vector<SFC_Key> A={3};
        vector<SFC_Key> B;
        __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);

    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short6\n");
        vector<SFC_Key> A={0,1,3,4,5,7,8,10,11,12,15};
        vector<SFC_Key> B={2,6,9,13,14};
        __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);

    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short7\n");
        vector<SFC_Key> A = { 1, 3, 4, 5, 7, 8, 10, 11, 12, 15 };
        vector<SFC_Key> B = { 0, 6, 9, 13 };
        __tisort_sortedarrays_merge_test_omp(A, B, C, Cndx);
    }

    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short8\n");
        vector<SFC_Key> A = { 0, 1, 3, 4, 5, 7, 8, 10, 11, 12, };
        vector<SFC_Key> B = { 14, 15 };
        __tisort_sortedarrays_merge_test_omp(A, B, C, Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short9\n");
        vector<SFC_Key> A={1,2,3,4,5,7,8,10,11,12,15};
        vector<SFC_Key> B={0};
        __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short10\n");
        vector<SFC_Key> A={1,2,3,4,5,7,8,10,11,12,15};
        vector<SFC_Key> B={16};
        __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);
    }
    //random
    UniqRandGen rand;
    vector<int> test_sizesA={40,128,250,516,1234,5832,12111,24127};
    vector<int> test_sizesB={10, 28, 77,516, 234, 832, 8111, 8127};
    for(int j=0;j<test_sizesA.size();++j)
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/random %d %d\n",test_sizesA[j],test_sizesB[j]);
        for(int i=0;i<20;++i)
        {
            vector<SFC_Key> A,B;
            rand.reset();
            fill_with_randomnumbers(A,test_sizesA[j],true,rand);
            fill_with_randomnumbers(B,test_sizesB[j],true,rand);
            __tisort_sortedarrays_merge_test_omp(A,B,C,Cndx);
        }
    }
    cout<<"tisearch_test/tisort_sortedarrays_merge_test ... passed\n";
}
inline void __tisort_mergesort_test(vector<SFC_Key> &key)
{
    vector<SFC_Key> key_ref;
    vector<ti_ndx_t> ndx_ref;
    vector<ti_ndx_t> ndx;

    build_ndx(key,ndx);
    key_ref=key;ndx_ref=ndx;
    sort_ref(key_ref,ndx_ref);

    MergeSort_BottomUp(key,ndx);

    assert(vectors_are_equal(key,key_ref));
    assert(vectors_are_equal(ndx,ndx_ref));
}
inline void tisort_mergesort_test()
{
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short00\n");
        //simple test
        vector<SFC_Key> key={3};
        __tisort_mergesort_test(key);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short0\n");
        //simple test
        vector<SFC_Key> key={3,1};
        __tisort_mergesort_test(key);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short1\n");
        //simple test
        vector<SFC_Key> key={3,2,6,9,7,1,5};
        __tisort_mergesort_test(key);
    }
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/short2\n");
        //reverse
        vector<SFC_Key> key={9,8,7,6,5,4,3,2,1};
        __tisort_mergesort_test(key);
    }
    //random
    UniqRandGen rand;
    vector<int> test_sizes={40,128,250,516,1234,5832,12111,24127};
    for(int j=0;j<test_sizes.size();++j)
    {
        if(verbose)printf("tisearch_test/tisort_sortedarrays_merge_test/random %d\n",test_sizes[j]);
        for(int i=0;i<20;++i)
        {
            vector<SFC_Key> key;
            rand.reset();
            fill_with_randomnumbers(key,test_sizes[j],false,rand);

            __tisort_mergesort_test(key);
        }
    }
    cout<<"tisearch_test/tisort_mergesort_test ... passed\n";
}
inline void tisort_test()
{
    vector<SFC_Key> key_ref;
    vector<ti_ndx_t> ndx_ref;
    vector<ti_ndx_t> ndx;

    cout<<"tisort_test ...\n";


    tisort_sortedarrays_merge_test();
    tisort_mergesort_test();


    cout<<"tisearch_test passed\n";
}



#endif /* TESTS_TICORE_TICORE_TISORT_HPP_ */
