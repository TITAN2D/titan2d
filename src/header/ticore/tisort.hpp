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
 * Author:
 * Description:
 *
 *******************************************************************
 */


#ifndef SRC_HEADER_TICORE_TISORT_HPP_
#define SRC_HEADER_TICORE_TISORT_HPP_

#include "../constant.h"

#include "tisearch.hpp"

#include "omp_mpi.hpp"
#include <cmath>



/**
 * Merge two sorted A and B array to C, C will be sorted, size of array A must be bigger or equal to B
 * this is not a merge sort algorithm,
 * the implementation is somewhat close to insertion sort
 */
inline void MergeSortedArraysToSortedArray(const int Asize, const SFC_Key * RESTRICT Akey, const ti_ndx_t * RESTRICT Andx,
                                           const int Bsize, const SFC_Key * RESTRICT Bkey, const ti_ndx_t * RESTRICT Bndx,
                                           const int Csize, SFC_Key * RESTRICT Ckey, ti_ndx_t * RESTRICT Cndx)
{
    ASSERT2(Asize+Bsize==Csize);
    //assert(Asize>=Bsize);
    int iA=0;
    int iC=0;
    for(int iB=0;iB<Bsize;iB++)
    {
        //find position to place Bkey[iB] and meanwhile fill with Akey[iA]
        while(Bkey[iB]>Akey[iA]&&iA<Asize)
        {
            Ckey[iC]=Akey[iA];
            Cndx[iC]=Andx[iA];
            iA++;
            iC++;
        }
        Ckey[iC]=Bkey[iB];
        Cndx[iC]=Bndx[iB];
        iC++;
    }
    for(int i=0;i<Asize-iA;++i)
    {
        Ckey[iC+i]=Akey[iA+i];
        Cndx[iC+i]=Andx[iA+i];
    }
}


/**
 * Merge two sorted A and B array to C, C will be sorted
 * this is not a merge sort algorithm,
 * the implementation is somewhat close to insertion sort
 *
 * for omp parallel should be enclosed in external #pragma omp parallel block
 *
 */
inline void MergeSortedArraysToSortedArray_omp(const int Asize, const SFC_Key * RESTRICT Akey, const ti_ndx_t * RESTRICT Andx,
                                           const int Bsize, const SFC_Key * RESTRICT Bkey, const ti_ndx_t * RESTRICT Bndx,
                                           const int Csize, SFC_Key * RESTRICT Ckey, ti_ndx_t * RESTRICT Cndx)
{
    int ithread=omp_get_thread_num();

    assert(Asize+Bsize==Csize);
    assert(Asize>=Bsize);

    if(Asize<threads_number*2)
    {
        if(ithread==0)MergeSortedArraysToSortedArray(Asize, Akey, Andx,Bsize, Bkey, Bndx, Csize, Ckey, Cndx);
    }
    else
    {
        int thread_portion=Bsize/threads_number;
        if(thread_portion<2)thread_portion=2;

        //start of Ai
        int iA0=ithread*thread_portion;
        //end of Ai, dot inclusive
        int iA1=(ithread==threads_number-1)?Asize:(ithread+1)*thread_portion;

        int iB0;
        int iB1;

        if(ithread==0)iB0=0;
        else iB0=find_insertion_spot(Bsize,Bkey,Akey[iA0]);


        if(ithread==threads_number-1)iB1=Bsize;
        else iB1=find_insertion_spot(Bsize,Bkey,Akey[iA1]);

        int iC0=iA0+iB0;
        int Bi_size=iB1-iB0;
        int Ai_size=iA1-iA0;

        MergeSortedArraysToSortedArray(Ai_size, Akey+iA0, Andx+iA0,Bi_size,Bkey+iB0,Bndx+iB0,Ai_size+Bi_size,Ckey+iC0,Cndx+iC0);
    }
    #pragma omp barrier
}

/**
 * function for MergeSort_BottomUp
 * should not be used on itself
 */
inline void BottomUpMerge(ti_ndx_t iLeft, ti_ndx_t iRight, ti_ndx_t iEnd, SFC_Key *A, SFC_Key *B, ti_ndx_t *ndxA, ti_ndx_t *ndxB)
{
    ti_ndx_t i0 = iLeft;
    ti_ndx_t i1 = iRight;
    ti_ndx_t j;

    /* while there are elements in the left or right lists */
    for (j = iLeft; j < iEnd; j++)
    {
      /* if left list head exists and is <= existing right list head */
      if (i0 < iRight && (i1 >= iEnd || A[i0] <= A[i1]))
      {
          B[j] = A[i0];
          ndxB[j] = ndxA[i0];
          i0 = i0 + 1;
      }
      else
      {
          B[j] = A[i1];
          ndxB[j] = ndxA[i1];
          i1 = i1 + 1;
      }
    }
}
/**
 * Merge sort, bottom up, serial
 *  array A[] has the items to sort; array B[] is a work array
 *  */
inline void MergeSort_BottomUp_serial(ti_ndx_t n, SFC_Key *m_A, ti_ndx_t *m_ndxA, SFC_Key *tmp_B, ti_ndx_t *tmp_ndxB)
{
    if(n==1)return;

    ti_ndx_t width;
    SFC_Key *A=m_A;
    ti_ndx_t *ndxA=m_ndxA;
    SFC_Key *B=tmp_B;
    ti_ndx_t *ndxB=tmp_ndxB;

    SFC_Key *T;
    ti_ndx_t *ndxT;



    // each 1-element run in A is already "sorted".

    // Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted
    for (width = 1; width < n; width = 2 * width)
    {
        //array A is full of runs of length width
        #pragma omp for schedule(static)
        for (ti_ndx_t i = 0; i < n; i = i + 2 * width)
        {
            //merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[]
            //  or copy A[i:n-1] to B[] ( if(i+width >= n) )
            BottomUpMerge(i, min(i+width, n), min(i+2*width, n), A, B,ndxA,ndxB);
        }
        //now work array B is full of runs of length 2*width
        //swap arrays for next iteration
        T=A;
        ndxT=ndxA;
        A = B;
        ndxA = ndxB;
        B = T;
        ndxB = ndxT;
        // now array A is full of runs of length 2*width
    }
    if(m_A!=A)
    {
        //i.e. A is not m_A and the values should be copied
        #pragma omp for schedule(static)
        for(ti_ndx_t i=0; i<n; i++)
        {
            m_A[i] = A[i];
            m_ndxA[i] = ndxA[i];
        }
    }
}
/**
 * Merge sort, bottom up
 *  array A[] has the items to sort; array B[] is a work array
 *  */
inline void MergeSort_BottomUp(ti_ndx_t n, SFC_Key *m_A, ti_ndx_t *m_ndxA, SFC_Key *tmp_B, ti_ndx_t *tmp_ndxB)
{
    if(n==1)return;
    if(n<2*threads_number)return MergeSort_BottomUp_serial(n, m_A, m_ndxA, tmp_B, tmp_ndxB);
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();
        ti_ndx_t width;
        SFC_Key *A=m_A;
        ti_ndx_t *ndxA=m_ndxA;
        SFC_Key *B=tmp_B;
        ti_ndx_t *ndxB=tmp_ndxB;

        SFC_Key *T;
        ti_ndx_t *ndxT;

        int log2threads=log2(threads_number)+0.49;
        int width_switch0=log2(n)+0.49;
        int width_switch;

        if(width_switch0<log2threads)width_switch=0;
        else width_switch=1<<(width_switch0-log2threads);


        /*if(ithread==0)
        {
            printf("n=%d threads_number=%d log2threads=%d width_switch0=%d  width_switch=%d\n",n,threads_number,log2threads,width_switch0,width_switch);
            for (width = 1; width < n; width = 2 * width)
            {
                printf("\twidth %d\n",width);
            }
        }*/


        #pragma omp barrier

        // each 1-element run in A is already "sorted".

        // Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted
        for (width = 1; width < width_switch; width = 2 * width)
        {
            //array A is full of runs of length width
            #pragma omp for schedule(static)
            for (ti_ndx_t i = 0; i < n; i = i + 2 * width)
            {
                //merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[]
                //  or copy A[i:n-1] to B[] ( if(i+width >= n) )
                BottomUpMerge(i, min(i+width, n), min(i+2*width, n), A, B,ndxA,ndxB);
            }
            //now work array B is full of runs of length 2*width
            //swap arrays for next iteration
            T=A;
            ndxT=ndxA;
            A = B;
            ndxA = ndxB;
            B = T;
            ndxB = ndxT;
            // now array A is full of runs of length 2*width
        }
        #pragma omp barrier
        for (width = width_switch; width < n; width = 2 * width)
        {
            //array A is full of runs of length width
            //#pragma omp for schedule(static)
            //if(ithread==0)
            for (ti_ndx_t i = 0; i < n; i = i + 2 * width)
            {
                //merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[]
                //  or copy A[i:n-1] to B[] ( if(i+width >= n) )
                //BottomUpMerge(i, min(i+width, n), min(i+2*width, n), A, B,ndxA,ndxB);
                ti_ndx_t iLeft=i;
                ti_ndx_t iRight=min(i+width, n);
                ti_ndx_t iEnd=min(i+2*width, n);
                MergeSortedArraysToSortedArray_omp(iRight-iLeft,A+iLeft,ndxA+iLeft,iEnd-iRight,A+iRight,ndxA+iRight,iEnd-iLeft,B+iLeft,ndxB+iLeft);

            }
            //now work array B is full of runs of length 2*width
            //swap arrays for next iteration
            T=A;
            ndxT=ndxA;
            A = B;
            ndxA = ndxB;
            B = T;
            ndxB = ndxT;
            // now array A is full of runs of length 2*width
        }
        /*for (i = iswitch; i < n; i = i + 2 * width)
        {
            //merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[]
            //  or copy A[i:n-1] to B[] ( if(i+width >= n) )
            if(ithread==0)BottomUpMerge(i, min(i+width, n), min(i+2*width, n), A, B,ndxA,ndxB);
        }
        #pragma omp barrier*/
        #pragma omp barrier
        if(m_A!=A)
        {
            //i.e. A is not m_A and the values should be copied
            #pragma omp for schedule(static)
            for(ti_ndx_t i=0; i<n; i++)
            {
                m_A[i] = A[i];
                m_ndxA[i] = ndxA[i];
            }
        }
    }
}


#endif /* SRC_HEADER_TICORE_TISORT_HPP_ */
