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

#ifndef TITAN2D_UTILS_H
#define TITAN2D_UTILS_H

#include "stdio.h"

class TitanTimings {
public:
	TitanTimings(){
		totalTime=0.0;
		meshAdaptionTime=0.0;
		stepTime=0.0;
		resultsOutputTime=0.0;
	}
	~TitanTimings(){}
	void print(){
		printf("Titan timings:\n");
		printf("  Total execution time:... %.3f\n",totalTime);
		printf("    Mesh adoption time:... %.3f (%.2f %%)\n",meshAdaptionTime,100.0*meshAdaptionTime/totalTime);
		printf("    Step time:............ %.3f (%.2f %%)\n",stepTime,100.0*stepTime/totalTime);
		printf("    Results dump time:.... %.3f (%.2f %%)\n",resultsOutputTime,100.0*resultsOutputTime/totalTime);

	}

	double totalTime;
	double meshAdaptionTime;
	double stepTime;
	double resultsOutputTime;
};

#endif
