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
 * $Id: setup_geoflow.C 224 2011-12-04 20:49:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"

void setup_geoflow(ElementsHashTable* El_Table, NodeHashTable* NodeTable, int myid, int nump, MatProps* matprops_ptr,
                   TimeProps *timeprops_ptr)
{
    int i;
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    /* zero out the fluxes for all of the nodes */
    //@NodesSingleLoop
    for(i = 0; i < NodeTable->elenode_.size(); i++)
    {
        if(NodeTable->status_[i]>=0)
            NodeTable->elenode_[i].zero_flux();
    }

    
    /* put the coord for the center node in the element */
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element *Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);
            int refined = Curr_El->refined_flag();
            if(Curr_El->adapted_flag() > 0) //if this is a refined element don't involve!!!
            {
                Curr_El->find_positive_x_side(NodeTable);
                Curr_El->calculate_dx(NodeTable);
                Curr_El->calc_topo_data(matprops_ptr);
                Curr_El->calc_gravity_vector(matprops_ptr);
            }
        }
    }
    
    /* transfer ghost elements to proper processors */
    move_data(nump, myid, El_Table, NodeTable, timeprops_ptr);
    
    /* calculate d_gravity array for each element */
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element *Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);
            int refined = Curr_El->refined_flag();
            if(Curr_El->adapted_flag() > 0) //if this is a refined element don't involve!!!
            {
                Curr_El->calc_d_gravity(El_Table);
            }
        }
    }
    return;
}
