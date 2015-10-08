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
#define TIMING3_DEFINE(t_start) double t_start

#define TIMING3_START(t_start) t_start=MPI_Wtime()
#define TIMING3_STOP(where,t_start) titanTimingsAlongSimulation.where+=MPI_Wtime()-t_start;titanTimings.where+=MPI_Wtime()-t_start;

#define IF_DEF_TIMING3(statement) statement

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
        printf("  Total execution time:................... %.3f\n", totalTime);
        printf("    Mesh adoption time:................... %.3f (%.2f %%)\n", meshAdaptionTime, 100.0 * meshAdaptionTime / totalTime);
        printf("      Refinement time:.................... %.3f (%.2f %%, %.2f %%)\n", refinementTime,
               100.0 * refinementTime / meshAdaptionTime, 100.0 * refinementTime / totalTime);
        IF_DEF_TIMING3(
        printf("        primaryRefinementsSearch:......... %.3f (%.2f %%, %.2f %%)\n", primaryRefinementsSearch,
                       100.0 * primaryRefinementsSearch / refinementTime, 100.0 * primaryRefinementsSearch / totalTime);
        )
        IF_DEF_TIMING3(
		printf("        triggeredRefinementsSearch:....... %.3f (%.2f %%, %.2f %%)\n", triggeredRefinementsSearch,
					   100.0 * triggeredRefinementsSearch / refinementTime, 100.0 * triggeredRefinementsSearch / totalTime);
		)
		IF_DEF_TIMING3(
		printf("        refineElements:................... %.3f (%.2f %%, %.2f %%)\n", refineElements,
					   100.0 * refineElements / meshAdaptionTime, 100.0 * refineElements / totalTime);
		)
		IF_DEF_TIMING3(
		printf("        refinedElementsNeigboursUpdate:... %.3f (%.2f %%, %.2f %%)\n", refinedElementsNeigboursUpdate,
					   100.0 * refinedElementsNeigboursUpdate / refinementTime, 100.0 * refinedElementsNeigboursUpdate / totalTime);
		)
        IF_DEF_TIMING3(
		printf("        sourcesRefinements:............... %.3f (%.2f %%, %.2f %%)\n", sourcesRefinements,
					   100.0 * sourcesRefinements / refinementTime, 100.0 * sourcesRefinements / totalTime);
		)

        printf("      Unrefinement time:.................. %.3f (%.2f %%, %.2f %%)\n", unrefinementTime,
               100.0 * unrefinementTime / meshAdaptionTime, 100.0 * unrefinementTime / totalTime);
        printf("    Step time:............................ %.3f (%.2f %%)\n", stepTime, 100.0 * stepTime / totalTime);
        printf("      Predictor time:..................... %.3f (%.2f %%, %.2f %%)\n", predictorStepTime,
               100.0 * predictorStepTime / stepTime, 100.0 * predictorStepTime / totalTime);
        printf("      Corrector time:..................... %.3f (%.2f %%, %.2f %%)\n", correctorStepTime,
               100.0 * correctorStepTime / stepTime, 100.0 * correctorStepTime / totalTime);
        printf("      Outline time:....................... %.3f (%.2f %%, %.2f %%)\n", outlineStepTime,
               100.0 * outlineStepTime / stepTime, 100.0 * outlineStepTime / totalTime);
        printf("      Slope Calc. time:................... %.3f (%.2f %%, %.2f %%)\n", slopesCalcTime,
               100.0 * slopesCalcTime / stepTime, 100.0 * slopesCalcTime / totalTime);
        printf("    Results dump time:.................... %.3f (%.2f %%)\n", resultsOutputTime, 100.0 * resultsOutputTime / totalTime);
        printf("    Flush ElemTable time:................. %.3f (%.2f %%)\n", flushElemTableTime, 100.0 * flushElemTableTime / totalTime);
        printf("    Flush NodeTable time:................. %.3f (%.2f %%)\n", flushNodeTableTime, 100.0 * flushNodeTableTime / totalTime);
        
        printf("\n");
    }
    void reset()
    {
        totalTime = 0.0;
        meshAdaptionTime = 0.0;
        refinementTime = 0.0;
        IF_DEF_TIMING3(primaryRefinementsSearch = 0.0);
        IF_DEF_TIMING3(triggeredRefinementsSearch = 0.0);
        IF_DEF_TIMING3(refineElements = 0.0);
        IF_DEF_TIMING3(refinedElementsNeigboursUpdate = 0.0);
        IF_DEF_TIMING3(sourcesRefinements = 0.0);
        unrefinementTime = 0.0;
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
        refinementTime /= steps;
        IF_DEF_TIMING3(primaryRefinementsSearch /= steps);
        IF_DEF_TIMING3(triggeredRefinementsSearch /= steps);
        IF_DEF_TIMING3(refineElements /= steps);
        IF_DEF_TIMING3(refinedElementsNeigboursUpdate /= steps);
        IF_DEF_TIMING3(sourcesRefinements /= steps);
        unrefinementTime /= steps;
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
    double refinementTime;
    IF_DEF_TIMING3(double primaryRefinementsSearch);
    IF_DEF_TIMING3(double triggeredRefinementsSearch);
    IF_DEF_TIMING3(double refineElements);
    IF_DEF_TIMING3(double refinedElementsNeigboursUpdate);
    IF_DEF_TIMING3(double sourcesRefinements);
    double unrefinementTime;
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
