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
    EleNodeRef(_ElemTable,_NodeTable)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    //positive_x_side_

    positive_x_side_xm[0] = 2;
    positive_x_side_yp[0] = 1;
    positive_x_side_ym[0] = 3;

    positive_x_side_xm[1] = 3;
    positive_x_side_yp[1] = 2;
    positive_x_side_ym[1] = 0;

    positive_x_side_xm[2] = 0;
    positive_x_side_yp[2] = 3;
    positive_x_side_ym[2] = 1;

    positive_x_side_xm[3] = 1;
    positive_x_side_yp[3] = 0;
    positive_x_side_ym[3] = 2;
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
            if(neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                ASSERT3(ElemTable->state_vars_[0][ElemTable->neighbor_ndx_[ineigh][ndx]] >= 0.0);
                if(state_vars_[0][neighbor_ndx_[ineigh][ndx]] >= contour_height)
                    return 1; //outside of pileheight contour line
            }
        }
    }
    return 0; //not on pileheight contour line
}
int ElementsProperties::if_source_boundary(ti_ndx_t ndx)
{
    ASSERT3(Influx_[0][ndx] >= 0.0); //currently mass sinks are not allowed

    if(Influx_[0][ndx] > 0.0)
    {
        for(int ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc_[ineigh][ndx] >= 0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                if(Influx_[0][ neighbor_ndx_[ineigh][ndx] ] <= 0.0)
                    return 2; //inside of line bounding area with a mass source
            }
    }

    else if(Influx_[0][ndx] == 0.0)
    {
        for(int ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc_[ineigh][ndx] >= 0.0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                ASSERT3(ElemTable->Influx_[0][ ElemTable->neighbor_ndx_[ineigh][ndx] ] >= 0.0);
                if(Influx_[0][ neighbor_ndx_[ineigh][ndx] ] != 0.0)
                    return 1; //outside of line bounding area with a mass source/sink
            }
    }
    else if(Influx_[0][ndx] < 0.0)
    {
        for(int ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc_[ineigh][ndx] >= 0.0)
            {
                //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                if(Influx_[0][ neighbor_ndx_[ineigh][ndx] ] >= 0.0)
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
int ElementsProperties::if_next_buffer_boundary(ti_ndx_t ndx,  double contour_height)
{

    int ineigh;
    ti_ndx_t neigh_ndx;
    int ifnextbuffer = 0;

    if(adapted_[ndx] <= 0)
        //GHOST element or element that should be deleted soon
        return (adapted_[ndx] - 1);

    if((adapted_[ndx] != BUFFER) && //this element is not in the buffer
            ((Influx_[0][ndx] == 0.0)/*&&(state_vars[0]<contour_height)*/) //&& //this element is OUTSIDE the buffer layer "circle"
    //(state_vars[0]>=GEOFLOW_TINY)
    )
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
        {
            if(neigh_proc_[ineigh][ndx] >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ASSERT3(ti_ndx_not_negative(ElemTable->lookup_ndx(ElemTable->neighbors_[ineigh][ndx])));
                neigh_ndx=neighbor_ndx_[ineigh][ndx];

                if((abs(adapted_[neigh_ndx]) == BUFFER) && (state_vars_[0][ndx]
                        <= state_vars_[0][neigh_ndx]))
                { //this element is next to a member of the old buffer layer
                  //if((ElemNeigh->get_adapted_flag())==BUFFER){ //this element is next to a member of the old buffer layer
                    ifnextbuffer = 1; //which means this element is a member of the next outer boundary of the buffer layer
                    break;
                }
            }
        }
    }

    if(ifnextbuffer == 1)
    {
        if((adapted_[ndx] >= NEWSON) || (generation_[ndx] == REFINE_LEVEL))
            return 2; //is a member of the buffer but doesn't need to be refined
        else
            return 1; //needs to be refined and some of its sons will be members
    }

    return 0;
}

void ElementsProperties::calc_flux_balance(ti_ndx_t ndx)
{
    int i, j;
    double flux[3] ={ 0.0, 0.0, 0.0 };
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side_[ndx];
    xm = positive_x_side_xm[xp];
    yp = positive_x_side_yp[xp];
    ym = positive_x_side_ym[xp];
/*
    switch (positive_x_side_[ndx])
    {
        case 0:
            xm = 2;
            yp = 1;
            ym = 3;
            break;
        case 1:
            xm = 3;
            yp = 2;
            ym = 0;
            break;
        case 2:
            xm = 0;
            yp = 3;
            ym = 1;
            break;
        case 3:
            xm = 1;
            yp = 0;
            ym = 2;
            break;
    }*/
    /*Node *nd_xp, *nd_xn, *nd_yp, *nd_yn;
    nd_xp = (Node*) NodeTable->lookup(node_key(xp + 4));
    nd_xn = (Node*) NodeTable->lookup(node_key(xm + 4));
    nd_yp = (Node*) NodeTable->lookup(node_key(yp + 4));
    nd_yn = (Node*) NodeTable->lookup(node_key(ym + 4));*/
    ti_ndx_t nd_xp, nd_xn, nd_yp, nd_yn;
    nd_xp = node_key_ndx_[xp + 4][ndx];
    nd_xn = node_key_ndx_[xm + 4][ndx];
    nd_yp = node_key_ndx_[yp + 4][ndx];
    nd_yn = node_key_ndx_[ym + 4][ndx];
    ASSERT3(ti_ndx_not_negative(nd_xp));
    ASSERT3(ti_ndx_not_negative(nd_xn));
    ASSERT3(ti_ndx_not_negative(nd_yp));
    ASSERT3(ti_ndx_not_negative(nd_yn));

    for(j = 0; j < 3; j++)
        flux[j] = dabs(node_refinementflux_[j][nd_xp] - node_refinementflux_[j][nd_xn])
                + dabs(node_refinementflux_[j][nd_yp] - node_refinementflux_[j][nd_yn]);

    double el_error=0.0;
    for(j = 0; j < 3; j++)
        el_error+=flux[j];

    el_error= 2.0 * el_error * el_error / (dx_[0][ndx] + dx_[1][ndx]) + WEIGHT_ADJUSTER; //W_A is so that elements with pile height = 0 have some weight.
    el_error_[0][ndx]=el_error;
    return;
}

/*! element_weight() cycles through the element Hashtable (listing of all
 *  elements) and for each element (that has not been refined this iteration
 *  and is not a ghost_element) calls Element member function
 *  Element::calc_flux_balance() (which returns a double precision value
 *  representing the weight that an element is assigned based on the
 *  magnitude of its net mass/momentum fluxes). Note that this value is
 *  adjusted to give non-zero weight even to elements with zero pile-heights
 *  The cumulative weights (along with a count of the evaluated elements)
 *  are stored in sub_weight[]; based on this, the return value for this
 *  function is calculated and stored in global_weight[] (i.e. the sum of
 *  sub_weight[] from all processors).
 */
double ElementsProperties::element_weight()
{
    int i, j, k, counter;
    double tiny = GEOFLOW_TINY;
    int el_counter = 0;
    double evalue = 1;
    double sub_weight[2] =
    { 0, 0 }; // second number is to keep track of the number of objects

    //-------------------go through all the elements of the subdomain and
    //-------------------find the edge states

    int no_of_buckets = ElemTable->get_no_of_buckets();
    vector<HashEntryLine> &bucket=ElemTable->bucket;
    tivector<Element> &elenode_=ElemTable->elenode_;

    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(Curr_El->adapted_flag() > 0)
            {
                //if this element doesn't belong on this processor don't involve!!!
                calc_flux_balance(Curr_El->ndx());
                //sub_weight[0] += *(Curr_El->get_el_error())+1.;
                if(Curr_El->state_vars(0) > GEOFLOW_TINY)
                {
                    sub_weight[1] += 1;
                    sub_weight[0] += Curr_El->el_error(0);
                }
                else if(Curr_El->adapted_flag() == BUFFER)
                {
                    sub_weight[1] += 0.1;
                    sub_weight[0] += Curr_El->el_error(0) * 0.1;
                }
            }
        }
    }

    double global_weight[2];
    i = MPI_Allreduce(sub_weight, global_weight, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //global_weight[0] = (global_weight[0]-global_weight[1])/global_weight[1];
    if(global_weight[1] > 0.0)
        global_weight[0] = (global_weight[0]) / global_weight[1]; //to protect from division by zero
    else
        global_weight[0] = 1.0; //just to make it nonzero

    return global_weight[0];
}

