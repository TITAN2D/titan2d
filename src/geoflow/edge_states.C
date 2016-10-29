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
 * $Id: edge_states.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"


/*! calc_edge_states() cycles through the element Hashtable (listing of all 
 *  elements) and for each element (that has not been refined this iteration 
 *  and is not a ghost_element) calls Element member function 
 *  Element::calc_edge_states() (which calculates the Riemann fluxes across 
 *  the element's boundaries), and adds local boundary-element outflow to 
 *  GIS map's cummulative outflow (defined as the mass flow off of the
 *  GIS map).  Also, the elements are checked for multiple pile-height values 
 */
void ElementsProperties::calc_edge_states(MatProps* matprops_ptr, TimeProps* timeprops_ptr,Integrator *integrator,
                      int myid, const int order, double &outflow)
{
    assert(ElemTable->all_elenodes_are_permanent);
    
    //-------------------go through all the elements of the subdomain and  
    //-------------------find the edge states
#ifdef BINARY_IDENTICAL_OMP
    vector<double> &localoutflow=dtmp;
    localoutflow.resize(elements_.size());
    
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        localoutflow[ndx]=0.0;
        if(adapted_[ndx] > 0)//if this element does not belong on this processor don't involve!!!
        {
                //if this element doesn't belong on this processor don't involve
                double pheight = state_vars_[0][ndx];
                elements_[ndx].calc_edge_states(ElemTable, NodeTable, matprops_ptr, integrator, myid, timeprops_ptr->dtime, order,
                                          &(localoutflow[ndx]));

                double pheight2 = state_vars_[0][ndx];
                if(pheight != pheight2)
                    printf("prolbem of changing height here,,,.....\n");
        }
    }

    outflow = 0.0;
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
        outflow+=localoutflow[ndx];
#else
    double localoutflow_sum=0.0;

    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK) reduction(+:localoutflow_sum)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] > 0)//if this element does not belong on this processor don't involve!!!
        {
                //if this element doesn't belong on this processor don't involve
                double pheight = state_vars_[0][ndx];
                double localoutflow;
                elements_[ndx].calc_edge_states(ElemTable, NodeTable, matprops_ptr, integrator, myid, timeprops_ptr->dtime, order,
                                          &localoutflow);
                localoutflow_sum+=localoutflow;
                double pheight2 = state_vars_[0][ndx];
                if(pheight != pheight2)
                    printf("prolbem of changing height here,,,.....\n");
        }
    }
    outflow=localoutflow_sum;
#endif
    return;
}
