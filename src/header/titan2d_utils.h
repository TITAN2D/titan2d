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

#ifdef _OPENMP
#include <omp.h>
#endif


#define TIMING1_DEFINE(t_start) double t_start
#define TIMING1_START(t_start) t_start=MPI_Wtime()
#define TIMING1_STOP(where,t_start) titanTimings.where=MPI_Wtime()-t_start;
#define TIMING1_STOPADD(where,t_start) titanTimingsAlongSimulation.where+=MPI_Wtime()-t_start;titanTimings.where+=MPI_Wtime()-t_start;
#define IF_DEF_TIMING1(statement) statement

#define TIMING2_DEFINE(t_start) double t_start
#define TIMING2_START(t_start) t_start=MPI_Wtime()
#define TIMING2_STOPADD(where,t_start) titanTimingsAlongSimulation.where+=MPI_Wtime()-t_start;titanTimings.where+=MPI_Wtime()-t_start;
#define IF_DEF_TIMING2(statement) statement


#define TIMING3_DEFINE(t_start) double t_start
#define TIMING3_START(t_start) t_start=MPI_Wtime()
#define TIMING3_STOPADD(where,t_start) titanTimingsAlongSimulation.where+=MPI_Wtime()-t_start;titanTimings.where+=MPI_Wtime()-t_start;
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
        printf("  Total execution time:................... %10.3f\n", totalTime);
        printf("    Mesh adoption time:................... %10.3f (%5.2f %%)\n", meshAdaptionTime, 100.0 * meshAdaptionTime / totalTime);
        printf("      Refinement time:.................... %10.3f (%5.2f %%, %5.2f %%)\n", refinementTime,
               100.0 * refinementTime / totalTime, 100.0 * refinementTime / meshAdaptionTime);
        IF_DEF_TIMING3(
		printf("        elementWeightCalc:................ %10.3f (%5.2f %%, %5.2f %%)\n", elementWeightCalc,
					   100.0 * elementWeightCalc / totalTime, 100.0 * elementWeightCalc / refinementTime);
		)
        IF_DEF_TIMING3(
        printf("        seedRefinementsSearch:............ %10.3f (%5.2f %%, %5.2f %%)\n", seedRefinementsSearch,
                       100.0 * seedRefinementsSearch / totalTime, 100.0 * seedRefinementsSearch / refinementTime);
        )
        IF_DEF_TIMING3(
		printf("        triggeredRefinementsSearch:....... %10.3f (%5.2f %%, %5.2f %%)\n", triggeredRefinementsSearch,
					   100.0 * triggeredRefinementsSearch / totalTime, 100.0 * triggeredRefinementsSearch / refinementTime);
		)
		IF_DEF_TIMING3(
		printf("        refineElements:................... %10.3f (%5.2f %%, %5.2f %%)\n", refineElements,
					   100.0 * refineElements / totalTime, 100.0 * refineElements / meshAdaptionTime);
		)
		IF_DEF_TIMING3(
		printf("        refinedElementsNeigboursUpdate:... %10.3f (%5.2f %%, %5.2f %%)\n", refinedElementsNeigboursUpdate,
					   100.0 * refinedElementsNeigboursUpdate / totalTime, 100.0 * refinedElementsNeigboursUpdate / refinementTime);
		)
        IF_DEF_TIMING3(
		printf("        refinementsPostProc:.............. %10.3f (%5.2f %%, %5.2f %%)\n", refinementsPostProc,
					   100.0 * refinementsPostProc / totalTime, 100.0 * refinementsPostProc / refinementTime);
		)

        printf("      Unrefinement time:.................. %10.3f (%5.2f %%, %5.2f %%)\n", unrefinementTime,
               100.0 * unrefinementTime / totalTime, 100.0 * unrefinementTime / meshAdaptionTime);
        printf("    Step time:............................ %10.3f (%5.2f %%)\n", stepTime, 100.0 * stepTime / totalTime);
        printf("      Predictor time:..................... %10.3f (%5.2f %%, %5.2f %%)\n", predictorStepTime,
               100.0 * predictorStepTime / totalTime, 100.0 * predictorStepTime / stepTime);
        printf("      Corrector time:..................... %10.3f (%5.2f %%, %5.2f %%)\n", correctorStepTime,
               100.0 * correctorStepTime / totalTime, 100.0 * correctorStepTime / stepTime);
        printf("      Outline time:....................... %10.3f (%5.2f %%, %5.2f %%)\n", outlineStepTime,
               100.0 * outlineStepTime / totalTime, 100.0 * outlineStepTime / stepTime);
        printf("      Slope Calc. time:................... %10.3f (%5.2f %%, %5.2f %%)\n", slopesCalcTime,
               100.0 * slopesCalcTime / totalTime, 100.0 * slopesCalcTime / stepTime);
        printf("    Results dump time:.................... %10.3f (%5.2f %%)\n", resultsOutputTime, 100.0 * resultsOutputTime / totalTime);
        printf("    Flush ElemTable time:................. %10.3f (%5.2f %%)\n", flushElemTableTime, 100.0 * flushElemTableTime / totalTime);
        printf("    Flush NodeTable time:................. %10.3f (%5.2f %%)\n", flushNodeTableTime, 100.0 * flushNodeTableTime / totalTime);
        
        printf("\n");
    }
    void reset()
    {
        totalTime = 0.0;
        meshAdaptionTime = 0.0;
        refinementTime = 0.0;
        IF_DEF_TIMING3(elementWeightCalc = 0.0);
        IF_DEF_TIMING3(seedRefinementsSearch = 0.0);
        IF_DEF_TIMING3(triggeredRefinementsSearch = 0.0);
        IF_DEF_TIMING3(refineElements = 0.0);
        IF_DEF_TIMING3(refinedElementsNeigboursUpdate = 0.0);
        IF_DEF_TIMING3(refinementsPostProc = 0.0);
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
        IF_DEF_TIMING3(elementWeightCalc /= steps);
        IF_DEF_TIMING3(seedRefinementsSearch /= steps);
        IF_DEF_TIMING3(triggeredRefinementsSearch /= steps);
        IF_DEF_TIMING3(refineElements /= steps);
        IF_DEF_TIMING3(refinedElementsNeigboursUpdate /= steps);
        IF_DEF_TIMING3(refinementsPostProc /= steps);
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
    IF_DEF_TIMING3(double elementWeightCalc);
    IF_DEF_TIMING3(double seedRefinementsSearch);
    IF_DEF_TIMING3(double triggeredRefinementsSearch);
    IF_DEF_TIMING3(double refineElements);
    IF_DEF_TIMING3(double refinedElementsNeigboursUpdate);
    IF_DEF_TIMING3(double refinementsPostProc);
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









#ifndef _OPENMP
#define omp_get_thread_num() 0
#endif

template<typename T>
void merge_vectors_from_threads(vector<T> &where, const vector< vector<T> > &what)
{
    /*size_t  N=0;
    for(int ithread;ithread<threads_number;++ithread)
    {
        N+=what[ithread].size();
    }*/
    where.resize(0);
    for(int ithread=0;ithread<threads_number;++ithread)
    {
        where.insert(where.end(),what[ithread].begin(), what[ithread].end());
    }
}

/*void parallel_for_get_my_part(const int ithread, const tisize_t N, ti_ndx_t &ndx_start, ti_ndx_t &ndx_end)
{
    ndx_start=ithread*N/threads_number;
    ndx_end=(ithread==threads_number-1)?N:(ithread+1)*N/threads_number;
}*/



#ifdef DEB2
#define DEB1
#endif


#ifdef DEB3
#define DEB1
#define DEB2
#endif

#ifdef DEB1
#define IF_DEB1(statement) {statement}
#define ASSERT1(statement) assert(statement)
#else
#define IF_DEB1(statement) {}
#define ASSERT1(statement) {}
#endif

#ifdef DEB2
#define IF_DEB2(statement) {statement}
#define ASSERT2(statement) assert(statement)
#else
#define IF_DEB2(statement) {}
#define ASSERT2(statement) {}
#endif

#ifdef DEB3
#define IF_DEB3(statement) {statement}
#define ASSERT3(statement) assert(statement)
#else
#define IF_DEB3(statement) {}
#define ASSERT3(statement) {}
#endif




#endif
