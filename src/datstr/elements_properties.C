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
//#define DEBUG_SAVE_ELEM

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include <math.h>



ElementsProperties::ElementsProperties(ElementsHashTable *_ElemTable, NodeHashTable* _NodeTable):
        ElemTable(_ElemTable),
        NodeTable(_NodeTable),

        elements_(ElemTable->elenode_),
        status_(ElemTable->status_),
        adapted_(ElemTable->adapted_),
        generation_(ElemTable->generation_),
        neigh_proc_(ElemTable->neigh_proc_),
        state_vars_(ElemTable->state_vars_),
        neighbor_ndx_(ElemTable->neighbor_ndx_)
{

}
int ElementsProperties::if_pile_boundary(ti_ndx_t ndx, double contour_height)
{
    int ineigh;

    ASSERT3(state_vars_[0][ndx] >= 0.0);

    if(state_vars_[0][ndx] >= contour_height)
    {
        for(int ineigh = 0; ineigh < 8; ineigh++)
        {
            if(neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                if(state_vars_[0][neighbor_ndx_[ineigh][ndx]] < contour_height)
                    return 2; //inside of pileheight contour line
            }
        }
    }
    else
    {
        for(int ineigh = 0; ineigh < 8; ineigh++)
        {
            if(ElemTable->neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                ASSERT3(ElemTable->state_vars_[0][ElemTable->neighbor_ndx_[ineigh][ndx]] >= 0.0);
                if(ElemTable->state_vars_[0][ElemTable->neighbor_ndx_[ineigh][ndx]] >= contour_height)
                    return 1; //outside of pileheight contour line
            }
        }
    }
    return 0; //not on pileheight contour line
}
int ElementsProperties::if_source_boundary(ti_ndx_t ndx)
{

    int ineigh;
    Element* ElemNeigh;


    ASSERT3(ElemTable->Influx_[0][ndx] >= 0.0); //currently mass sinks are not allowed

    if(ElemTable->Influx_[0][ndx] > 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(ElemTable->neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                if(ElemTable->Influx_[0][ ElemTable->neighbor_ndx_[ineigh][ndx] ] <= 0.0)
                    return 2; //inside of line bounding area with a mass source
            }
    }

    else if(ElemTable->Influx_[0][ndx] == 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(ElemTable->neigh_proc_[ineigh][ndx] >= 0.0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                ASSERT3(ElemTable->Influx_[0][ ElemTable->neighbor_ndx_[ineigh][ndx] ] >= 0.0);
                if(ElemTable->Influx_[0][ ElemTable->neighbor_ndx_[ineigh][ndx] ] != 0.0)
                    return 1; //outside of line bounding area with a mass source/sink
            }
    }
    else if(ElemTable->Influx_[0][ndx] < 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(ElemTable->neigh_proc_[ineigh][ndx] >= 0.0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                if(ElemTable->Influx_[0][ ElemTable->neighbor_ndx_[ineigh][ndx] ] >= 0.0)
                    return -1; //inside of line bounding area with a mass sink
            }
    }

    return 0; //not on line bounding area with mass source/sink
}
int ElementsProperties::if_first_buffer_boundary(ti_ndx_t ndx, double contour_height) const
{

    int ineigh;
    int iffirstbuffer = 0;
    ti_ndx_t neig_ndx;

    if(ElemTable->adapted_[ndx] <= 0)
        return (ElemTable->adapted_[ndx] - 1);

    ASSERT3(ElemTable->state_vars_[0][ndx] >= 0.0);
    ASSERT3(ElemTable->Influx_[0][ndx] >= 0.0);

    if((ElemTable->state_vars_[0][ndx] < contour_height) && (ElemTable->Influx_[0][ndx] == 0.0))
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(ElemTable->neigh_proc_[ineigh][ndx]  >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                neig_ndx=ElemTable->neighbor_ndx_[ineigh][ndx];
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));

                if((ElemTable->state_vars_[0][neig_ndx] >= contour_height) || (ElemTable->Influx_[0][neig_ndx] > 0.0))
                {
                    iffirstbuffer = 1;
                    break;
                }
            }
    }
    else
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(ElemTable->neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                neig_ndx=ElemTable->neighbor_ndx_[ineigh][ndx];
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));

                if((ElemTable->state_vars_[0][neig_ndx] < contour_height) && (ElemTable->Influx_[0][neig_ndx] == 0.0))
                {
                    iffirstbuffer = 1;
                    break;
                }
            }
    }

    if(iffirstbuffer)
    {
        if((ElemTable->adapted_[ndx] >= NEWSON) || (ElemTable->generation_[ndx] == REFINE_LEVEL))
            return 2; //is a member of the buffer but doesn't need to be refined
        else
            return 1; //needs to be refined and some of its sons will be members
    }

    return 0;
}
