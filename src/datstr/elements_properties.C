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
    xm = (2 + xp) % 4;
    yp = (1 + xp) % 4;
    ym = (3 + xp) % 4;

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

void ElementsProperties::slopes(MatProps* matprops_ptr)
{
    assert(ElemTable->all_elenodes_are_permanent);
    double gamma=matprops_ptr->gamma;
    #pragma omp parallel for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] > 0)//if this element does not belong on this processor don't involve!!!
        {
            get_slopes(ndx,gamma);
        }
    }
    return;
}

void ElementsProperties::get_slopes(ti_ndx_t ndx, double gamma)
{
    int j = 0, bc = 0;
    /* check to see if this is a boundary */
    while (j < 4 && bc == 0)
    {
        if(neigh_proc_[j][ndx] == INIT)
            bc = 1;
        j++;
    }
    if(bc == 1)
    {
        for(j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
            d_state_vars_[j][ndx]=0.0;
        return;
    }

    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side_[ndx];
    xm = (2 + xp) % 4;
    yp = (1 + xp) % 4;
    ym = (3 + xp) % 4;

    /* x direction */
    ti_ndx_t ep = neighbor_ndx_[xp][ndx]; //(Element*) (ElemTable->lookup(&neighbor(xp)[0]));
    ti_ndx_t em = neighbor_ndx_[xm][ndx]; //(Element*) (ElemTable->lookup(&neighbor(xm)[0]));
    ti_ndx_t ep2 = ti_ndx_doesnt_exist;
    ti_ndx_t em2 = ti_ndx_doesnt_exist;
    //check if element has 2 neighbors on either side
    ti_ndx_t ndtemp = node_key_ndx_[xp + 4][ndx]; //(Node*) NodeTable->lookup(&node_key[xp + 4][0]);
    if(node_info_[ndtemp] == S_C_CON)
    {
        ep2 = neighbor_ndx_[xp + 4][ndx]; //(Element*) (ElemTable->lookup(&neighbor[xp + 4][0]));
        ASSERT3(neigh_proc_[xp + 4][ndx] >= 0 && ti_ndx_not_negative(ep2));
    }
    ndtemp = node_key_ndx_[xm + 4][ndx]; //(Node*) NodeTable->lookup(&node_key[xm + 4][0]);
    if(node_info_[ndtemp] == S_C_CON)
    {
        em2 = neighbor_ndx_[xm + 4][ndx]; //(Element*) (ElemTable->lookup(&neighbor[xm + 4][0]));
        ASSERT3(neigh_proc_[xm + 4][ndx] >= 0 && ti_ndx_not_negative(em2));
    }

    double dp, dm, dc, dxp, dxm;
    dxp = coord_[0][ep] - coord_[0][ndx];
    dxm = coord_[0][ndx] - coord_[0][em];
    for(j = 0; j < NUM_STATE_VARS; j++)
    {
        dp = (state_vars_[j][ep] - state_vars_[j][ndx]) / dxp;
        if(ti_ndx_not_negative(ep2))
            dp = .5 * (dp + (state_vars_[j][ep2] - state_vars_[j][ndx]) / dxp);
        dm = (state_vars_[j][ndx] - state_vars_[j][em]) / dxm;
        if(ti_ndx_not_negative(em2))
            dm = .5 * (dm + (state_vars_[j][ndx] - state_vars_[j][em2]) / dxm);

        dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
        //do slope limiting
        d_state_vars_[j][ndx]=0.5 * (c_sgn(dp) + c_sgn(dm)) * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
    }

    /* y direction */
    ep = neighbor_ndx_[yp][ndx];        //(Element*) (ElemTable->lookup(&neighbor(yp)[0]));
    em = neighbor_ndx_[ym][ndx];        //(Element*) (ElemTable->lookup(&neighbor(ym)[0]));
    ep2 = ti_ndx_doesnt_exist;
    em2 = ti_ndx_doesnt_exist;
    //check if element has 2 neighbors on either side
    ndtemp = node_key_ndx_[yp + 4][ndx];        //(Node*) NodeTable->lookup(&node_key[yp + 4][0]);
    if(node_info_[ndtemp] == S_C_CON)
    {
        ep2 = neighbor_ndx_[yp + 4][ndx];       //(Element*) (ElemTable->lookup(&neighbor[yp + 4][0]));
        ASSERT3(neigh_proc_[yp + 4][ndx] >= 0 && ti_ndx_not_negative(ep2));
    }
    ndtemp = node_key_ndx_[ym + 4][ndx];        //(Node*) NodeTable->lookup(&node_key[ym + 4][0]);
    if(node_info_[ndtemp] == S_C_CON)
    {
        em2 = neighbor_ndx_[ym + 4][ndx];       //(Element*) (ElemTable->lookup(&neighbor[ym + 4][0]));
        ASSERT3(neigh_proc_[ym + 4][ndx] >= 0 && ti_ndx_not_negative(em2));
    }

    dxp = coord_[1][ep] - coord_[1][ndx];
    dxm = coord_[1][ndx] - coord_[1][em];
    for(j = 0; j < NUM_STATE_VARS; j++)
    {
        dp = (state_vars_[j][ep] - state_vars_[j][ndx]) / dxp;
        if(ti_ndx_not_negative(ep2))
            dp = .5 * (dp + (state_vars_[j][ep2] - state_vars_[j][ndx]) / dxp);
        dm = (state_vars_[j][ndx] - state_vars_[j][em]) / dxm;
        if(ti_ndx_not_negative(em2))
            dm = .5 * (dm + (state_vars_[j][ndx] - state_vars_[j][em2]) / dxm);

        dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
        //do slope limiting
        d_state_vars_[j + NUM_STATE_VARS][ndx]=0.5 * (c_sgn(dp) + c_sgn(dm))
                                           * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
    }

    return;
}
void ElementsProperties::calc_wet_dry_orient()
{
    #pragma omp for schedule(static,TITAN2D_DINAMIC_CHUNK)
    for(ti_ndx_t ndx = 0; ndx < elements_.size(); ndx++)
    {
        if(adapted_[ndx] <= 0)continue;//if this element does not belong on this processor don't involve!!!
        //elements_[ndx].calc_wet_dry_orient(ElemTable);
        calc_wet_dry_orient(ndx);
    }
}
// the element member function calc_wet_dry_orient() calculates the orientation of the dryline (drylineorient),
// the wet length (Swet), the location of the drypoint, and the location of the wetpoint... it does NOT calculate
// the wet area (Awet)... these quantities are used in the wetted area adjustment of fluxes. calc_wet_dry_orient()
// is not coded for generic element orientation, positive_x_side must be side 1.  Keith wrote this may 2007
void ElementsProperties::calc_wet_dry_orient(ti_ndx_t ndx)
{

    int ifsidewet[4], numwetsides = 0;
    int ineigh;

    for(ineigh = 0; ineigh < 4; ineigh++)
    {
        if(neigh_proc_[ineigh][ndx] == -1)
        {
            //edge of map and cell has same wetness as the cell
            ifsidewet[ineigh] = (state_vars_[0][ndx] > GEOFLOW_TINY) ? 1 : 0;
        }
        else
        {
            if(state_vars_[0][neighbor_ndx_[ineigh][ndx]] > GEOFLOW_TINY)
            {
                //first neighbor on this side is wet
                ifsidewet[ineigh] = 1;
            }
            else if(neigh_proc_[ineigh + 4][ndx] == -2)
            {
                //only one neighbor on this side and it's not wet
                ifsidewet[ineigh] = 0;
            }
            else
            {
                //since first neighbor on this side is not wet, the edge has the wetness of the second neighbor on this side
                ifsidewet[ineigh] = (state_vars_[0][neighbor_ndx_[ineigh+4][ndx]] > GEOFLOW_TINY) ? 1 : 0;
            }
        }
        numwetsides += ifsidewet[ineigh];
    }

    if((ifsidewet[0] == ifsidewet[2]) && (ifsidewet[1] == ifsidewet[3]))
    {
        //if opposite sides of the element are the same (both wet or dry)
        iwetnode_[ndx]=8;
        drypoint_[0][ndx]=0.0;
        drypoint_[1][ndx]=0.0;
        if(state_vars_[0][ndx] > GEOFLOW_TINY)
        {
            Awet_[ndx]=1.0;
            Swet_[ndx]=1.0;
        }
        else
        {
            Awet_[ndx]=0.0;
            Swet_[ndx]=0.0;
        }
        //?1.0:0.0;
    }
    else if(numwetsides == 2)
    {
        //having exactly 2 adjacent wet edges means it has a diagonal orientation

        Swet_[ndx]=sqrt(2.0 * ((Awet_[ndx] > 0.5) ? 1.0 - Awet_[ndx] : Awet_[ndx])); //edge length of small triangle
        drypoint_[0][ndx]=0.5 * (1.0 - Swet_[ndx]);
        drypoint_[1][ndx]=drypoint_[0][ndx];
        if(Awet_[ndx] > 0.5)
            Swet_[ndx]=1.0 - Swet_[ndx]; //the small triangle is dry not wet

        if(ifsidewet[3] && ifsidewet[0])
        {
            iwetnode_[ndx]=0;
            if(Awet_[ndx] <= 0.5){
                drypoint_[0][ndx] = -drypoint_[0][ndx];
                drypoint_[1][ndx] =  drypoint_[0][ndx];
            }

        }
        else if(ifsidewet[0] && ifsidewet[1])
        {
            iwetnode_[ndx]=1;
            if(Awet_[ndx] <= 0.5)
                drypoint_[1][ndx] =  -drypoint_[1][ndx];
            else
                drypoint_[0][ndx] = -drypoint_[0][ndx];
        }
        else if(ifsidewet[1] && ifsidewet[2])
        {
            iwetnode_[ndx]=2;
            if(Awet_[ndx] > 0.5)
            {
                drypoint_[0][ndx] = -drypoint_[0][ndx];
                drypoint_[1][ndx] = drypoint_[0][ndx];
            }
        }
        else if(ifsidewet[2] && ifsidewet[3])
        {
            iwetnode_[ndx]=3;
            if(Awet_[ndx] > 0.5)
                drypoint_[1][ndx] =  -drypoint_[1][ndx];
            else
                drypoint_[0][ndx] = -drypoint_[0][ndx];
        }
    }
    else
    {
        //numwetsides is 1 or 3 i.e. it's a vertical or horizontal orientation
        if(numwetsides == 1)
        {
            //find the one wet side
            for(ineigh = 0; ineigh < 4; ineigh++)
                if(ifsidewet[ineigh])
                    break;
        }
        else
        {
            //find the one dry side
            for(ineigh = 0; ineigh < 4; ineigh++)
                if(!ifsidewet[ineigh])
                    break;
            //find the wet side opposite the one dry side
            ineigh = (ineigh + 2) % 4;
        }
        assert((-1 < ineigh) && (ineigh < 4));
        Swet_[ndx]=Awet_[ndx];

        iwetnode_[ndx]=ineigh + 4;
        switch (iwetnode_[ndx])
        {
            case 4:
                drypoint_[0][ndx] = 0.0;
                drypoint_[1][ndx] = -0.5 + Swet_[ndx];
                break;
            case 5:
                drypoint_[0][ndx] = +0.5 - Swet_[ndx];
                drypoint_[1][ndx] = 0.0;
                break;
            case 6:
                drypoint_[0][ndx] = 0.0;
                drypoint_[1][ndx] = +0.5 - Swet_[ndx];
                break;
            case 7:
                drypoint_[0][ndx] = -0.5 + Swet_[ndx];
                drypoint_[1][ndx] = 0.0;
                break;
            default:
                assert(0);
        }
    }

    if(iwetnode_[ndx] == 8)
        Awet_[ndx]=(state_vars_[0][ndx] > GEOFLOW_TINY) ? 1.0 : 0.0;

    return;
}

