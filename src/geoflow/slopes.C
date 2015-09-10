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
    Element* Curr_El;
    int Nelms = El_Table->getNumberOfLocalElements();
    //if this element does not belong on this processor don't involve!!!
    Element** Elms = (Element**) El_Table->getLocalElementsValues();
//#pragma omp parallel for private(currentPtr,Curr_El)
    for(i = 0; i < Nelms; i++)
    {
        Elms[i]->get_slopes(El_Table, NodeTable, matprops_ptr->gamma);
    }
    return;
}
