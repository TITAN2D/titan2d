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

#ifndef TITAN_SIMULATION_H
#define TITAN2D_SIMULATION_H


class cxxTitanSimulation {
public:
	cxxTitanSimulation();
	~cxxTitanSimulation();
	void run();

	int myid;
	int numprocs;
};
#endif
