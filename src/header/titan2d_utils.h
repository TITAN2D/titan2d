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

#include <stdio.h>
#include <time.h>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>

using namespace std;

#include "titan2d.h"
#include "ticore.hpp"
#include "constant.h"

#include "ticore/omp_mpi.hpp"
#include "tivector.h"

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
        printf("    Mesh adaptation time:................. %10.3f (%5.2f %%)\n", meshAdaptionTime, 100.0 * meshAdaptionTime / totalTime);
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

typedef std::chrono::high_resolution_clock Clock;

#define PROFILING1_DEFINE(t_start) Clock::time_point t_start
#define PROFILING1_START(t_start) t_start=Clock::now();
#define PROFILING1_STOP(where,t_start) titanProfiling.where=Clock::now()-t_start
#define PROFILING1_STOPADD(where,t_start) titanProfiling.where+=Clock::now()-t_start
#define PROFILING1_STOPADD_RESTART(where,t_start) titanProfiling.where+=Clock::now()-t_start;t_start=Clock::now();
#define IF_DEF_PROFILING1(statement) statement

#define PROFILE1_TIMINGS_PRINTING(variable) std::cout<<"  "<<#variable<<" "<<variable.count()<<" "<< std::chrono::duration<double, std::chrono::milliseconds::period>(variable).count() <<"\n"

#define PROFILING2_DEFINE(t_start) Clock::time_point t_start
#define PROFILING2_START(t_start) t_start=Clock::now();
#define PROFILING2_STOP(where,t_start) titanProfiling.where=Clock::now()-t_start
#define PROFILING2_STOPADD(where,t_start) titanProfiling.where+=Clock::now()-t_start
#define IF_DEF_PROFILING2(statement) statement

#define PROFILE2_TIMINGS_PRINTING(variable) std::cout<<"  "<<#variable<<" "<<variable.count()<<" "<< std::chrono::duration<double, std::chrono::milliseconds::period>(variable).count() <<"\n"

#define PROFILING3_DEFINE(t_start) Clock::time_point t_start
#define PROFILING3_START(t_start) t_start=Clock::now()
#define PROFILING3_STOP(where,t_start) titanProfiling.where=Clock::now()-t_start
#define PROFILING3_STOPADD(where,t_start) titanProfiling.where+=Clock::now()-t_start
#define PROFILING3_STOPADD_RESTART(where,t_start) titanProfiling.where+=Clock::now()-t_start;t_start=Clock::now();
#define IF_DEF_PROFILING3(statement) statement

#define PROFILE3_TIMINGS_PRINTING(variable) std::cout<<"  "<<#variable<<" "<<variable.count()<<" "<< std::chrono::duration<double, std::chrono::milliseconds::period>(variable).count() <<"\n"


class TitanProfiling
{
public:
    TitanProfiling()
    {
        zero=Clock::duration::zero();
        reset();
    }
    ~TitanProfiling()
    {
    }
    void print(const char *title = NULL)
    {
        if(title == NULL)
        {
            printf("\nTitan profiling, cycles, miliseconds:\n");
        }
        else
        {
            printf("%s\n", title);
        }
        PROFILE1_TIMINGS_PRINTING(tsim_iter);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_update_topo);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_adapt);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_step);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_saverestart);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_output);
        PROFILE1_TIMINGS_PRINTING(tsim_iter_post);

        PROFILE1_TIMINGS_PRINTING(PrimaryRefinementsFinder_findSeedRefinements_loop);
        PROFILE1_TIMINGS_PRINTING(PrimaryRefinementsFinder_findSeedRefinements_merge);
        PROFILE1_TIMINGS_PRINTING(BuferFirstLayerRefinementsFinder_findSeedRefinements_loop);
        PROFILE1_TIMINGS_PRINTING(BuferFirstLayerRefinementsFinder_findSeedRefinements_merge);
        PROFILE1_TIMINGS_PRINTING(BuferNextLayerRefinementsFinder_findSeedRefinements_loop);
        PROFILE1_TIMINGS_PRINTING(BuferNextLayerRefinementsFinder_findSeedRefinements_merge);

        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_prolog);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_element_weight);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_reset_adapt);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_allrefine);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_allrefine_buffer_handling);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_htflush2);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_remove_elements);
        PROFILE3_TIMINGS_PRINTING(HAdapt_adapt_epilog);

        PROFILE3_TIMINGS_PRINTING(step_other);
        PROFILE3_TIMINGS_PRINTING(step_slopesCalc);
        PROFILE3_TIMINGS_PRINTING(step_get_coef_and_eigen);
        PROFILE3_TIMINGS_PRINTING(step_adapt_fluxsrc_region);
        PROFILE3_TIMINGS_PRINTING(step_predict);
        PROFILE3_TIMINGS_PRINTING(step_apply_bc);
        PROFILE3_TIMINGS_PRINTING(step_corrector);
        PROFILE3_TIMINGS_PRINTING(step_calc_edge_states);
        PROFILE3_TIMINGS_PRINTING(step_outline);
        PROFILE3_TIMINGS_PRINTING(step_calc_wet_dry_orient);
        PROFILE3_TIMINGS_PRINTING(step_calc_stats);

        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_other);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_find_corners_sides_bubbles);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side0_find_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side0_add_new_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side1_find_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side1_add_new_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side2_find_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side2_add_new_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side3_find_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side3_add_new_nodes);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_int_nodes_calc_keys);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_int_nodes_alloc);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side0_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side1_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side2_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_side3_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_int_nodes_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_new_elm_aloc);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_new_elm_init);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refineElements_update_neighbours_ndx_on_ghosts);

        PROFILE3_TIMINGS_PRINTING(HAdapt_refinedNeighboursUpdate_mpi_prep);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refinedNeighboursUpdate_onproc_updates);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refinedNeighboursUpdate_mpi_update);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refinedNeighboursUpdate_removeElements);
        PROFILE3_TIMINGS_PRINTING(HAdapt_refinedNeighboursUpdate_update_neighbours_ndx_on_ghosts);

        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_init);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_find);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_updateBrothersIndexes);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_unrefine_neigh_update);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_unrefine_interp_neigh_update);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_delete_oldsons);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_update_neighbours_ndx_on_ghosts);
        PROFILE3_TIMINGS_PRINTING(HAdaptUnrefine_unrefine_calc_wet_dry_orient);

        PROFILE3_TIMINGS_PRINTING(flushTable_sort_prep);
        PROFILE3_TIMINGS_PRINTING(flushTable_sort);
        PROFILE3_TIMINGS_PRINTING(flushTable_sort2);
        PROFILE3_TIMINGS_PRINTING(flushTable_sort_post);
        PROFILE3_TIMINGS_PRINTING(flushTable_HashTable_reorder);
        PROFILE3_TIMINGS_PRINTING(flushTable_NodeHashTable_reorder);
        PROFILE3_TIMINGS_PRINTING(flushTable_ElementsHashTable_reorder);
        PROFILE3_TIMINGS_PRINTING(flushTable_ElementsHashTable_reorder2);
        PROFILE3_TIMINGS_PRINTING(flushTable_ElementsHashTable_updateLocalElements);

        //PROFILE1_TIMINGS_PRINTING();
        //PROFILE3_TIMINGS_PRINTING();
    }
    void reset()
    {
        IF_DEF_PROFILING1(tsim_iter=zero);
        IF_DEF_PROFILING1(tsim_iter_update_topo=zero);
        IF_DEF_PROFILING1(tsim_iter_adapt=zero);
        IF_DEF_PROFILING1(tsim_iter_step=zero);
        IF_DEF_PROFILING1(tsim_iter_saverestart=zero);
        IF_DEF_PROFILING1(tsim_iter_output=zero);
        IF_DEF_PROFILING1(tsim_iter_post=zero);

        IF_DEF_PROFILING1(PrimaryRefinementsFinder_findSeedRefinements_loop=zero);
        IF_DEF_PROFILING1(PrimaryRefinementsFinder_findSeedRefinements_merge=zero);
        IF_DEF_PROFILING1(BuferFirstLayerRefinementsFinder_findSeedRefinements_loop=zero);
        IF_DEF_PROFILING1(BuferFirstLayerRefinementsFinder_findSeedRefinements_merge=zero);
        IF_DEF_PROFILING1(BuferNextLayerRefinementsFinder_findSeedRefinements_loop=zero);
        IF_DEF_PROFILING1(BuferNextLayerRefinementsFinder_findSeedRefinements_merge=zero);

        IF_DEF_PROFILING3(HAdapt_adapt_prolog=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_element_weight=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_reset_adapt=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_allrefine=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_allrefine_buffer_handling=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_htflush2=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_remove_elements=zero);
        IF_DEF_PROFILING3(HAdapt_adapt_epilog=zero);

        IF_DEF_PROFILING3(step_other=zero);
        IF_DEF_PROFILING3(step_slopesCalc=zero);
        IF_DEF_PROFILING3(step_get_coef_and_eigen=zero);
        IF_DEF_PROFILING3(step_adapt_fluxsrc_region=zero);
        IF_DEF_PROFILING3(step_predict=zero);
        IF_DEF_PROFILING3(step_apply_bc=zero);
        IF_DEF_PROFILING3(step_corrector=zero);
        IF_DEF_PROFILING3(step_calc_edge_states=zero);
        IF_DEF_PROFILING3(step_outline=zero);
        IF_DEF_PROFILING3(step_calc_wet_dry_orient=zero);
        IF_DEF_PROFILING3(step_calc_stats=zero);

        IF_DEF_PROFILING3(HAdapt_refineElements_other = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_find_corners_sides_bubbles = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side0_find_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side0_add_new_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side1_find_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side1_add_new_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side2_find_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side2_add_new_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side3_find_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side3_add_new_nodes = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_int_nodes_calc_keys = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_int_nodes_alloc = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side0_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side1_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side2_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_side3_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_int_nodes_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_new_elm_init = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_new_elm_aloc = zero);
        IF_DEF_PROFILING3(HAdapt_refineElements_update_neighbours_ndx_on_ghosts = zero);

        IF_DEF_PROFILING3(HAdapt_refinedNeighboursUpdate_mpi_prep = zero);
        IF_DEF_PROFILING3(HAdapt_refinedNeighboursUpdate_onproc_updates = zero);
        IF_DEF_PROFILING3(HAdapt_refinedNeighboursUpdate_mpi_update = zero);
        IF_DEF_PROFILING3(HAdapt_refinedNeighboursUpdate_removeElements = zero);
        IF_DEF_PROFILING3(HAdapt_refinedNeighboursUpdate_update_neighbours_ndx_on_ghosts = zero);

        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_init = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_find = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_updateBrothersIndexes = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_unrefine_neigh_update = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_unrefine_interp_neigh_update = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_delete_oldsons = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_update_neighbours_ndx_on_ghosts = zero);
        IF_DEF_PROFILING3(HAdaptUnrefine_unrefine_calc_wet_dry_orient = zero);

        IF_DEF_PROFILING3(flushTable_sort_prep = zero);
        IF_DEF_PROFILING3(flushTable_sort = zero);
        IF_DEF_PROFILING3(flushTable_sort2 = zero);
        IF_DEF_PROFILING3(flushTable_sort_post = zero);
        IF_DEF_PROFILING3(flushTable_HashTable_reorder = zero);
        IF_DEF_PROFILING3(flushTable_NodeHashTable_reorder = zero);
        IF_DEF_PROFILING3(flushTable_ElementsHashTable_reorder = zero);
        IF_DEF_PROFILING3(flushTable_ElementsHashTable_reorder2 = zero);
        IF_DEF_PROFILING3(flushTable_ElementsHashTable_updateLocalElements = zero);
        //IF_DEF_PROFILING1(=zero);
        //IF_DEF_PROFILING3(=zero);
    }
public:
    Clock::duration zero;
    IF_DEF_PROFILING1(Clock::duration tsim_iter);
    IF_DEF_PROFILING1(Clock::duration tsim_iter_update_topo);
    IF_DEF_PROFILING1(Clock::duration tsim_iter_adapt);
    IF_DEF_PROFILING1(Clock::duration tsim_iter_step);
    IF_DEF_PROFILING1(Clock::duration tsim_iter_saverestart);
    IF_DEF_PROFILING1(Clock::duration tsim_iter_output);
    IF_DEF_PROFILING2(Clock::duration tsim_iter_post);
    //IF_DEF_PROFILING2(Clock::duration );

    IF_DEF_PROFILING1(Clock::duration PrimaryRefinementsFinder_findSeedRefinements_loop);
    IF_DEF_PROFILING1(Clock::duration PrimaryRefinementsFinder_findSeedRefinements_merge);
    IF_DEF_PROFILING1(Clock::duration BuferFirstLayerRefinementsFinder_findSeedRefinements_loop);
    IF_DEF_PROFILING1(Clock::duration BuferFirstLayerRefinementsFinder_findSeedRefinements_merge);
    IF_DEF_PROFILING1(Clock::duration BuferNextLayerRefinementsFinder_findSeedRefinements_loop);
    IF_DEF_PROFILING1(Clock::duration BuferNextLayerRefinementsFinder_findSeedRefinements_merge);

    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_prolog);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_element_weight);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_reset_adapt);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_allrefine);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_allrefine_buffer_handling);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_htflush2);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_remove_elements);
    IF_DEF_PROFILING3(Clock::duration HAdapt_adapt_epilog);

    IF_DEF_PROFILING3(Clock::duration step_other);
    IF_DEF_PROFILING3(Clock::duration step_slopesCalc);
    IF_DEF_PROFILING3(Clock::duration step_get_coef_and_eigen);
    IF_DEF_PROFILING3(Clock::duration step_adapt_fluxsrc_region);
    IF_DEF_PROFILING3(Clock::duration step_predict);
    IF_DEF_PROFILING3(Clock::duration step_apply_bc);
    IF_DEF_PROFILING3(Clock::duration step_corrector);
    IF_DEF_PROFILING3(Clock::duration step_calc_edge_states);
    IF_DEF_PROFILING3(Clock::duration step_outline);
    IF_DEF_PROFILING3(Clock::duration step_calc_wet_dry_orient);
    IF_DEF_PROFILING3(Clock::duration step_calc_stats);

    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_other);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_find_corners_sides_bubbles);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side0_find_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side0_add_new_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side1_find_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side1_add_new_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side2_find_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side2_add_new_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side3_find_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side3_add_new_nodes);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_int_nodes_calc_keys);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_int_nodes_alloc);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side0_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side1_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side2_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_side3_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_int_nodes_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_new_elm_init);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_new_elm_aloc);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refineElements_update_neighbours_ndx_on_ghosts);

    IF_DEF_PROFILING3(Clock::duration HAdapt_refinedNeighboursUpdate_mpi_prep);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refinedNeighboursUpdate_onproc_updates);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refinedNeighboursUpdate_mpi_update);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refinedNeighboursUpdate_removeElements);
    IF_DEF_PROFILING3(Clock::duration HAdapt_refinedNeighboursUpdate_update_neighbours_ndx_on_ghosts);

    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_init);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_find);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_updateBrothersIndexes);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_unrefine_neigh_update);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_unrefine_interp_neigh_update);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_delete_oldsons);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_update_neighbours_ndx_on_ghosts);
    IF_DEF_PROFILING1(Clock::duration HAdaptUnrefine_unrefine_calc_wet_dry_orient);

    IF_DEF_PROFILING1(Clock::duration flushTable_sort_prep);
    IF_DEF_PROFILING1(Clock::duration flushTable_sort);
    IF_DEF_PROFILING1(Clock::duration flushTable_sort2);
    IF_DEF_PROFILING1(Clock::duration flushTable_sort_post);
    IF_DEF_PROFILING1(Clock::duration flushTable_HashTable_reorder);
    IF_DEF_PROFILING1(Clock::duration flushTable_NodeHashTable_reorder);
    IF_DEF_PROFILING1(Clock::duration flushTable_ElementsHashTable_reorder);
    IF_DEF_PROFILING1(Clock::duration flushTable_ElementsHashTable_reorder2);
    IF_DEF_PROFILING1(Clock::duration flushTable_ElementsHashTable_updateLocalElements);
    //IF_DEF_PROFILING1(Clock::duration );
    //IF_DEF_PROFILING3(Clock::duration );
};

extern TitanProfiling titanProfiling;



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

template<typename T>
void merge_vectors_from_threads_to0(vector< vector<T> > &v)
{
    /*size_t  N=0;
    for(int ithread;ithread<threads_number;++ithread)
    {
        N+=what[ithread].size();
    }*/
    for(int ithread=1;ithread<v.size();++ithread)
    {
        v[0].insert(v[0].end(),v[ithread].begin(), v[ithread].end());
        v[ithread].resize(0);
    }
}
/**
 * move the content of vector< vector<T> > vec to vec[0]
 * usefull for omp
 */
template<typename T>
void merge_vectors_from_threads_to0_omp(vector< vector<T> > &vec)
{
    ti_ndx_t number_of_new_elenodes=0;
            #pragma omp parallel
            {
                int ithread=omp_get_thread_num();
                ti_ndx_t start=0;
                if(ithread==0)
                    for(int jthread=0;jthread<threads_number;++jthread)
                        number_of_new_elenodes+=vec[jthread].size();
                else
                    for(int jthread=0;jthread<ithread;++jthread)
                        start+=vec[jthread].size();
                #pragma omp barrier
                if(ithread==0)
                {
                    vec[0].resize(number_of_new_elenodes);
                }
                #pragma omp barrier
                if(ithread!=0)
                {
                    for(ti_ndx_t i=0;i<vec[ithread].size();++i)
                    {
                        vec[0][start+i]=vec[ithread][i];
                    }
                    vec[ithread].resize(0);
                }
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

//!string conversions

extern vector<string> coord_names;
extern vector<string> SinglePhaseVarNames;
extern vector<string> TwoPhasesVarNames;

//!file system
bool ti_dir_exists(const char *filename);
void ti_mkdir(const char *filename);

#endif
