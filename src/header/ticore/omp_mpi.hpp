/*******************************************************************
 * Copyright (C) 2016 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:
 * Description: mpi and openmp utils and helpers should be include instead of mpi.h and omp.h
 *
 *******************************************************************
 */

#ifndef SRC_HEADER_TICORE_OMP_MPI_HPP_
#define SRC_HEADER_TICORE_OMP_MPI_HPP_

#ifdef _OPENMP
#include <omp.h>

#define IF_OMP(statement) statement

#else

#define IF_OMP(statement)
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1

#endif

#ifdef USE_MPI

#include <mpi.h>
#define IF_MPI(statement) statement

#else

#define IF_MPI(statement)

#define MPI_Comm_size(comm_world, numprocs) *(numprocs)=1
#define MPI_Comm_rank(comm_world, myid) *(myid)=0

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <chrono>

inline double MPI_Wtime()
{
    int64_t t=std::chrono::duration_cast<std::chrono::microseconds>
            (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    return (double)t*0.000001;
}

#endif

#endif /* SRC_HEADER_TICORE_OMP_MPI_HPP_ */
