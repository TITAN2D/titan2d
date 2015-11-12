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
 * $Id: slopes.C 128 2007-06-07 19:51:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/geoflow.h"

void slopes(ElementsHashTable* El_Table, NodeHashTable* NodeTable, MatProps* matprops_ptr)
{
    int i;
    //-------------------go through all the elements of the subdomain------------------------
    //-------------------and   --------------------------
    
    /* mdj 2007-02 */
    //#pragma omp parallel for private(currentPtr,Curr_El)
    for(ti_ndx_t ndx = 0; ndx < El_Table->size(); ndx++)
    {
        if(El_Table->adapted_[ndx] > 0)//if this element does not belong on this processor don't involve!!!
        {
            El_Table->elenode_[ndx].get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
        }
    }
    return;
}
