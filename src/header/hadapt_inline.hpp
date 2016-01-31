/*******************************************************************
 * Copyright (C) 2015 University at Buffalo
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
#ifndef HADAPT_INLINE_HPP
#define	HADAPT_INLINE_HPP

#include "hadapt.h"


inline int HAdaptUnrefine::check_unrefinement(ti_ndx_t ndx, double target)
{
    int unrefine_flag = 1;

    //  put in good element check here!!!
    //if((if_pile_boundary(El_Table,target))||
    //   (if_source_boundary(El_Table)))
    //if((state_vars[0] >= target)||(Influx[0]>0.0))
    if(adapted_[ndx] != NOTRECADAPTED)
        //This rules out NEWFATHERs, NEWSONs, BUFFERs, GHOSTs, TOBEDELETEDs, and OLDFATERs
        //This is a redundant check but is is better to be safe than sorry
        return (0);

    for(int i = 0; i < 8; i++)
    {
        if(((neigh_proc_[i][ndx] != myprocess_[ndx]) && (neigh_proc_[i][ndx] >= 0) && (generation_[ndx] <= 0)) || (neigh_gen_[i][ndx] > generation_[ndx]))
            return (0);
    }

    return (1);
}


#endif	/* ELEMENTS_PROPERTIES_INLINE_HPP */

