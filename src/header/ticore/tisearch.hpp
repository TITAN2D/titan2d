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

#ifndef SRC_HEADER_TICORE_TISEARCH_HPP_
#define SRC_HEADER_TICORE_TISEARCH_HPP_

#include "../constant.h"
#include "../tivector.h"



/**
 * find an index i such that A[i-1]<b<A[i]
 */
template<typename T>
ti_ndx_t find_insertion_spot(const ti_ndx_t N, const T* RESTRICT A, const T b)
{
    //@TODO do binary tree search
    if(N==0)return 0;

    if(b<A[0])return 0;
    if(A[N-1]<b)return N;


    int i;

    if(N < HASHTABLE_LOOKUP_LINSEARCH)
    {
        for(i = 0; i < N; i++)
        {
            if(A[i-1]<b && b<A[i])return i;
        }
    }
    else
    {
        int i0, i1, i2;
        i0 = 0;
        i1 = N / 2;
        i2 = N - 1;
        while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
        {
            if(b>A[i1])
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
            if(A[i-1]<b && b<A[i])return i;
        }
    }
    return ti_ndx_unknown;
}


#endif /* SRC_HEADER_TICORE_TISEARCH_HPP_ */
