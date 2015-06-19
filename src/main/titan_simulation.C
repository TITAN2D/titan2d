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
 */

#include <stdio.h>
#include <mpi.h>
#include "../header/titan_simulation.h"

cxxTitanSimulation::cxxTitanSimulation()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}
cxxTitanSimulation::~cxxTitanSimulation()
{

}
int hpfem();
void cxxTitanSimulation::run()
{
	printf("cxxTitanSimulation::run\n");
	hpfem();

}
