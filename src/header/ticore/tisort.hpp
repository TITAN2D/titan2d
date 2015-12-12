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

inline int BottomUpMerge(ti_ndx_t iLeft, ti_ndx_t iRight, ti_ndx_t iEnd, SFC_Key *A, SFC_Key *B, ti_ndx_t *ndxA, ti_ndx_t *ndxB)
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
    return 0;
}

/**
 * Merge sort, bottom up
 *  array A[] has the items to sort; array B[] is a work array
 *  */
inline void BottomUpSort(ti_ndx_t n, SFC_Key *m_A, ti_ndx_t *m_ndxA, SFC_Key *tmp_B, ti_ndx_t *tmp_ndxB)
{
    ti_ndx_t i;
    ti_ndx_t width;
    SFC_Key *A=m_A;
    ti_ndx_t *ndxA=m_ndxA;
    SFC_Key *B=tmp_B;
    ti_ndx_t *ndxB=tmp_ndxB;

    SFC_Key *T;
    ti_ndx_t *ndxT;

    /* each 1-element run in A is already "sorted". */

    /* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted */
    for (width = 1; width < n; width = 2 * width)
    {
        /* array A is full of runs of length width */
        for (i = 0; i < n; i = i + 2 * width)
        {
            /* merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
            /*  or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
            BottomUpMerge(i, min(i+width, n), min(i+2*width, n), A, B,ndxA,ndxB);
        }

        /* now work array B is full of runs of length 2*width */
        /* swap arrays for next iteration */
        T=A;
        ndxT=ndxA;
        A = B;
        ndxA = ndxB;
        B = T;
        ndxB = ndxT;
        /* now array A is full of runs of length 2*width */
    }
    if(m_A!=A)
    {
        /*i.e. A is not m_A and the values should be copied*/
        //printf("m_A!=A\n");
        for(i=0; i<n; i++)
        {
            m_A[i] = A[i];
            m_ndxA[i] = ndxA[i];
        }
    }
}

/**
 * Merge two sorted A and B array to C, C will be sorted
 */
inline void MergeSortedArraysToSortedArray(const int Asize, const SFC_Key * RESTRICT Akey, const ti_ndx_t * RESTRICT Andx,
                                           const int Bsize, const SFC_Key * RESTRICT Bkey, const ti_ndx_t * RESTRICT Bndx,
                                           const int Csize, SFC_Key * RESTRICT Ckey, ti_ndx_t * RESTRICT Cndx)
{
    assert(Asize+Bsize==Csize);
    int iA=0;
    int iC=0;
    for(int iB=0;iB<Bsize;iB++)
    {
        //find position to place Bkey[iB] and meanwhile fill with Akey[iA]
        while(Bkey[iB]>Akey[iA]&&iA<Asize)
        {
            Ckey[iC]=Akey[iA];
            Cndx[iC]=Andx[iA];
            if(iC>0)
            {
                if(Ckey[iC-1]>Ckey[iC])
                {
                    assert(0);
                }
            }
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
#endif /* SRC_HEADER_TICORE_TISORT_HPP_ */
