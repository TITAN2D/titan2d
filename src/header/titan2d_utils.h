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

class TitanTimings
{
public:
    TitanTimings()
    {
        reset();
    }
    ~TitanTimings()
    {
    }
    void print(const char *title = NULL)
    {
        if(title == NULL)
        {
            printf("\nTitan timings, seconds:\n");
        }
        else
        {
            printf("%s\n", title);
        }
        printf("  Total execution time:... %.3f\n", totalTime);
        printf("    Mesh adoption time:... %.3f (%.2f %%)\n", meshAdaptionTime, 100.0 * meshAdaptionTime / totalTime);
        printf("    Step time:............ %.3f (%.2f %%)\n", stepTime, 100.0 * stepTime / totalTime);
        printf("      Predictor time:..... %.3f (%.2f %%, %.2f %%)\n", predictorStepTime,
               100.0 * predictorStepTime / stepTime, 100.0 * predictorStepTime / totalTime);
        printf("      Corrector time:..... %.3f (%.2f %%, %.2f %%)\n", correctorStepTime,
               100.0 * correctorStepTime / stepTime, 100.0 * correctorStepTime / totalTime);
        printf("      Outline time:....... %.3f (%.2f %%, %.2f %%)\n", outlineStepTime,
               100.0 * outlineStepTime / stepTime, 100.0 * outlineStepTime / totalTime);
        printf("      Slope Calc. time:... %.3f (%.2f %%, %.2f %%)\n", slopesCalcTime,
               100.0 * slopesCalcTime / stepTime, 100.0 * slopesCalcTime / totalTime);
        printf("    Results dump time:.... %.3f (%.2f %%)\n", resultsOutputTime, 100.0 * resultsOutputTime / totalTime);
        printf("    Flush ElemTable time:. %.3f (%.2f %%)\n", flushElemTableTime, 100.0 * flushElemTableTime / totalTime);
        printf("    Flush NodeTable time:. %.3f (%.2f %%)\n", flushNodeTableTime, 100.0 * flushNodeTableTime / totalTime);
        
        printf("\n");
    }
    void reset()
    {
        totalTime = 0.0;
        meshAdaptionTime = 0.0;
        stepTime = 0.0;
        predictorStepTime = 0.0;
        correctorStepTime = 0.0;
        outlineStepTime = 0.0;
        slopesCalcTime = 0.0;
        resultsOutputTime = 0.0;
        flushElemTableTime = 0.0;
        flushNodeTableTime = 0.0;
        
    }
    void devideBySteps(double steps)
    {
        totalTime /= steps;
        meshAdaptionTime /= steps;
        stepTime /= steps;
        predictorStepTime /= steps;
        correctorStepTime /= steps;
        outlineStepTime /= steps;
        slopesCalcTime /= steps;
        resultsOutputTime /= steps;
        flushElemTableTime /= steps;
        flushNodeTableTime /= steps;
    }
    
    double totalTime;
    double meshAdaptionTime;
    double stepTime;
    double predictorStepTime;
    double correctorStepTime;
    double outlineStepTime;
    double slopesCalcTime;
    double resultsOutputTime;
    double flushElemTableTime;
    double flushNodeTableTime;
    
};
extern TitanTimings titanTimings;
extern TitanTimings titanTimingsAlongSimulation;
#endif
