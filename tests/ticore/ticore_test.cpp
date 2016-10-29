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

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include <assert.h>

#include <iostream>
#include <vector>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../titest_common.hpp"

#include "../../src/header/tivector.h"

#include "ticore_tisearch.hpp"
#include "ticore_tisort.hpp"

int threads_number;
int verbose;
int main(int argc, char *argv[])
{
    process_cmdline(argc,argv);

    cout<<"tocore test\n";

    tisearch_test();
    tisort_test();

    cout<<"test passed\n";
    return 0;
}


