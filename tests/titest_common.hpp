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
 * Description: common routines for tests
 *
 *******************************************************************
 */

#ifndef TESTS_TITEST_COMMON_HPP_
#define TESTS_TITEST_COMMON_HPP_

#ifdef _OPENMP
#include <omp.h>
#endif


extern int threads_number;
extern int verbose;

void process_cmdline(int argc, char *argv[])
{
    int num_threads = 1;
    verbose=0;
    if(argc > 1)
    {
        for(int i = 0; i < argc; i++)
        {
            string arg(argv[i]);
            //strip off arguments which should not make to python
            if(arg == "-nt")
            {

                if((i + 1 >= argc) || (!isdigit(argv[i + 1][0])))
                {
                    printf("Error: Incorrect command line format\n");
                    break;
                }
                ++i;

                num_threads = atoi(argv[i]);

                //set openmp threads
#ifndef _OPENMP
                if(num_threads > 1)
                    printf("Warning: This version was compiled without openmp support!\n");

                num_threads = 1;
                printf("Set threads number to %d\n\n", num_threads);
#endif
                continue;
            }
            if(arg == "-v")
            {
                verbose=1;
            }
        }
    }

#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    threads_number=omp_get_max_threads();
#else
    threads_number = 1;
#endif
    printf("Max threads number is %d\n\n", threads_number);
}



#endif /* TESTS_TITEST_COMMON_HPP_ */
