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

#ifndef TESTS_TICORE_TICORE_TISEARCH_HPP_
#define TESTS_TICORE_TICORE_TISEARCH_HPP_

#include <assert.h>

#include <iostream>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/header/sfc.h"
#include "../../src/header/ticore/tisearch.hpp"


inline void tisearch_test()
{
    cout<<"tisearch_test ...\n";
    {
        //some where in the middle
        vector<SFC_Key> vec={0,1,2,3,5,6,7};
        SFC_Key f=4;
        assert(find_insertion_spot<SFC_Key>(vec.size(),&(vec[0]),f)==4);
    }
    {
        //all the way to the left
        vector<SFC_Key> vec={1,2,3,5,6,7};
        SFC_Key f=0;
        assert(find_insertion_spot<SFC_Key>(vec.size(),&(vec[0]),f)==0);
    }
    {
        //all the way to the right
        vector<SFC_Key> vec={1,2,3,5,6,7};
        SFC_Key f=8;
        assert(find_insertion_spot<SFC_Key>(vec.size(),&(vec[0]),f)==6);
    }
    cout<<"tisearch_test passed\n";
}


#endif /* TESTS_TICORE_TICORE_TISEARCH_HPP_ */
