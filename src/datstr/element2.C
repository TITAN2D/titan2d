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
 * $Id: element2.C 232 2012-03-27 00:33:41Z dkumar $ 
 */
//#define DEBUG_SAVE_ELEM

#define DISABLE_DRY_FLUX_ZEROING   //this disables the zeroing of fluxes through dry sides facet of thin layer control

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include <math.h>

//stopping criteria define statements have been moved to geoflow.h

//#define PRINT_GIS_ERRORS

/*  original element   */
Element::Element(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC* b, int mat,
                 int* elm_loc_in, double pile_height, int myid, const SFC_Key& opposite_brother)
{
    if(NUM_STATE_VARS==3)elementType(ElementType::SinglePhase);
    else if(NUM_STATE_VARS==6)elementType(ElementType::TwoPhases);
    else elementType(ElementType::UnknownElementType);

    int ikey;
    int i, j;
    set_adapted_flag(NOTRECADAPTED);

    for(i = 0; i < NUM_STATE_VARS; i++)
        prev_state_vars[i] = 0.;
    
    for(i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
        d_state_vars[i] = 0.;

    set_father(sfc_key_zero);
    for(i = 0; i < 4; i++){
        set_brother(i, sfc_key_zero);
        set_son(i, sfc_key_zero);
    }
    set_lb_key(sfc_key_zero);
    
    set_lb_weight(1.0);
    set_new_old(OLD);
    generation(0); //--first generation 
    set_material(mat);
    for(i = 0; i < EQUATIONS; i++)
        set_el_error(i, 0.0);
    
    set_key(nodekeys[8]); //--using bubble key to represent the element
    
    for(i = 0; i < 8; i++)
        set_node_key(i, nodekeys[i]);
    
    for(i = 0; i < 4; i++)
    {
        set_neigh_proc(i, n_pro[i]);
        set_neigh_proc(i + 4, -2); //-- -2 means regular element
        if(neigh_proc(i) != -1){
            set_neighbor(i, neigh[i]);
            set_neighbor(i + 4, neigh[i]);
        }
        else
        {
            set_neighbor(i, sfc_key_zero);
            set_neighbor(i + 4, sfc_key_zero);
        }
    }
    
    bcptr(b);
    
    for(i = 0; i < 5; i++)
        set_order(i, POWER);   //--used in initial uniform mesh
                
    for(i = 0; i < 8; i++)
        get_neigh_gen(i, 0);
    
    set_no_of_eqns(EQUATIONS);
    
    int help = 0;
    for(i = 0; i < 4; i++)
        help += order(i) * (no_of_eqns());
    help += pow((float) (order(4) - 1), 2) * (no_of_eqns());
    set_ndof(help);
    
    set_refined_flag(0);
    
    
    set_myprocess(myid);
    set_elm_loc(0, elm_loc_in[0]);
    set_elm_loc(1, elm_loc_in[1]);
    calc_which_son();
    
    set_brother(which_son(), key());
    set_brother((which_son() + 2) % 4, opposite_brother);

    switch (which_son())
    {
        case 0:
            set_brother(1, neighbor(1));
            set_brother(3, neighbor(2));
            break;
        case 1:
            set_brother(0, neighbor(3));
            set_brother(2, neighbor(2));
            break;
        case 2:
            set_brother(1, neighbor(0));
            set_brother(3, neighbor(3));
            break;
        case 3:
            set_brother(0, neighbor(0));
            set_brother(2, neighbor(1));
            break;
    }
    set_opposite_brother_flag(1);
    
    set_new_old(OLD);
    if(elementType() == ElementType::SinglePhase)
    {
        state_vars[0] = pile_height;
        state_vars[1] = 0;
        if(state_vars[0] != 0)
            state_vars[1] = 0.0001;
        state_vars[2] = 0;
    }

    if(elementType() == ElementType::TwoPhases)
    {
        state_vars[0] = pile_height;
        state_vars[1] = pile_height;
        state_vars[2] = state_vars[4] = 0.;
        if(state_vars[0] != 0)
            state_vars[2] = state_vars[4] = 0.0001;
        state_vars[3] = state_vars[5] = 0;
    }
    set_shortspeed(0.0);
    
    set_iwetnode(8);
    drypoint(0, 0.0);
    drypoint(1, 0.0);
    set_Awet((pile_height > GEOFLOW_TINY) ? 1.0 : 0.0);
    set_Swet(Awet());
    

    if(elementType() == ElementType::TwoPhases)
    {
        prev_state_vars[0] = pile_height;
        prev_state_vars[1] = pile_height;
        prev_state_vars[2] = prev_state_vars[4] = 0.;
        if(prev_state_vars[0] != 0)
            prev_state_vars[2] = prev_state_vars[4] = 0.0001;
        prev_state_vars[3] = prev_state_vars[5] = 0;
    }
    if(elementType() == ElementType::SinglePhase)
    {
        prev_state_vars[0] = pile_height;
        prev_state_vars[1] = 0;
        if(prev_state_vars[0] != 0)
            prev_state_vars[1] = 0.0001;
        prev_state_vars[2] = 0;
    }

    for(i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
        d_state_vars[i] = 0;
    
    for(i = 0; i < NUM_STATE_VARS; i++)
        Influx(i, 0.0);

    //printf("creating an original element\n");
    
    set_stoppedflags(2);

    if(elementType() == ElementType::TwoPhases)
    {
        // initialize kactxy
        kactxy(0,0.0);
        kactxy(1,0.0);
        effect_kactxy(0, 0.0);
        effect_kactxy(1, 0.0);
        set_effect_bedfrict(0.0);
        set_effect_tanbedfrict(0.0);
    }
}

//used for refinement
Element::Element(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen,
                 int elm_loc_in[], int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
                 HashTable *El_Table, HashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather,
                 double Awetfather, double *drypoint_in)
{
    if(NUM_STATE_VARS == 3)
        elementType(ElementType::SinglePhase);
    else if(NUM_STATE_VARS == 6)
        elementType(ElementType::TwoPhases);
    else
        elementType(ElementType::UnknownElementType);

    int i;
    int ikey;
            
    set_adapted_flag(NEWSON);
    
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        prev_state_vars[i] = 0.;
        Influx(i, 0.);
    }
    for(i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
        d_state_vars[i] = 0.;
    
    set_father(sfc_key_zero);
    for(i = 0; i < 4; i++){
        set_brother(i, sfc_key_zero);
        set_son(i, sfc_key_zero);
    }
    set_lb_key(sfc_key_zero);
    set_lb_weight(1.0);
    set_myprocess(myid);
    generation(gen); //--first generation
    set_opposite_brother_flag(1);
    set_material(mat);
    for(i = 0; i < EQUATIONS; i++)
        set_el_error(i, 0.0);
    
    set_key(nodekeys[8]); //--using buble key to represent the element
    
    set_elm_loc(0, elm_loc_in[0]);
    set_elm_loc(1, elm_loc_in[1]);
    
    for(i = 0; i < 8; i++)
        set_node_key(i, nodekeys[i]);
    
    for(i = 0; i < 4; i++)
    {
        set_neigh_proc(i, n_pro[i]);
        set_neigh_proc(i + 4, -2); //-- -2 means regular element
        if(neigh_proc(i) != -1){
            set_neighbor(i, neigh[i]);
            set_neighbor(i + 4, neigh[i]);
        }
        else
        {
            set_neighbor(i, sfc_key_zero);
            set_neighbor(i + 4, sfc_key_zero);
        }
        
        get_neigh_gen(i, gen_neigh[i]);
        get_neigh_gen(i + 4, gen_neigh[i]);
    }
    
    bcptr(b);
    
    for(i = 0; i < 5; i++)
    {
        set_order(i, *(ord + i));   //--used in initial uniform mesh
        //cout<<"In the constructor the order"<<order[i]<<"\n\n"<<flush;
    }
    
    set_no_of_eqns(EQUATIONS);
    
    int help = 0;
    for(i = 0; i < 4; i++)
        help += order(i) * (no_of_eqns());
    help += pow((float) (order(4) - 1), 2) * (no_of_eqns());
    set_ndof(help);
    
    set_refined_flag(0);
    
    
    set_new_old(NEW);
    //geoflow stuff
    dx(0,.5 * fthTemp->dx(0));  //assume constant refinement
    dx(1,.5 * fthTemp->dx(1));
    
    set_iwetnode(iwetnodefather);
    
    drypoint(0, drypoint_in[0]);
    drypoint(1, drypoint_in[1]);
    
    double myfractionoffather;
    if((Awetfather == 0.0) || (Awetfather == 1.0))
    {
        set_Awet(Awetfather);
        myfractionoffather = 1.0;
    }
    else
    {
        set_Awet(convect_dryline(dx(0),dx(1), 0.0)); //dx is a dummy stand in for convection speed... value doesn't matter because it's being multiplied by a timestep of zero
        myfractionoffather = Awet() / Awetfather;
    }
    set_Swet(1.0);
    
    double dxx = coord_in[0] - fthTemp->coord(0);
    double dyy = coord_in[1] - fthTemp->coord(1);
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        // state_vars[i] = fthTemp->state_vars[i]+fthTemp->d_state_vars[i]*dxx + fthTemp->d_state_vars[i+NUM_STATE_VARS]*dyy;
        state_vars[i] = fthTemp->state_vars[i] * myfractionoffather;
        prev_state_vars[i] = fthTemp->prev_state_vars[i] * myfractionoffather;
        set_shortspeed(fthTemp->shortspeed());
    }
    
    if(state_vars[0] < 0.)
        state_vars[0] = 0.;
    
    find_positive_x_side(NodeTable);
    
    calc_topo_data(matprops_ptr);
    calc_gravity_vector(matprops_ptr);
    
    set_coord(0, coord_in[0]);
    set_coord(1, coord_in[1]);
    
    set_stoppedflags(fthTemp->stoppedflags());
    
    return;
}
/*********************************
 making a father element from its sons
 *****************************************/
Element::Element(Element* sons[], HashTable* NodeTable, HashTable* El_Table, MatProps* matprops_ptr)
{
    if(NUM_STATE_VARS == 3)
        elementType(ElementType::SinglePhase);
    else if(NUM_STATE_VARS == 6)
        elementType(ElementType::TwoPhases);
    else
        elementType(ElementType::UnknownElementType);

            
    set_adapted_flag(NEWFATHER);
    
    for(int i = 0; i < NUM_STATE_VARS; i++)
    {
        prev_state_vars[i] = 0.;
        Influx(i, 0.);
    }
    for(int i = 0; i < DIMENSION * NUM_STATE_VARS; i++)
        d_state_vars[i] = 0.;
    
    set_father(sfc_key_zero);
    for(int i = 0; i < 4; i++){
        set_brother(i, sfc_key_zero);
        set_son(i, sfc_key_zero);
    }

    
    int i, j, ison, isonneigh, ineigh;
    
    set_key(sons[2]->node_key(0));
    
    for(ison = 0; ison < 4; ison++)
    {
        sons[ison]->set_adapted_flag(OLDSON);
        set_son(ison, sons[ison]->key());
        sons[ison]->set_father(key());
    }
    
    set_lb_key(sfc_key_zero);
    set_lb_weight(1.0);
    set_new_old(NEW);
    set_opposite_brother_flag(0);
    set_stoppedflags(2);
    for(i = 0; i < EQUATIONS; i++)
        set_el_error(i, 0.0);
    
    for(ison = 0; ison < 4; ison++)
    {
        if(sons[ison]->stoppedflags() < stoppedflags())
            set_stoppedflags(sons[ison]->stoppedflags());
    }
    
    set_father(sfc_key_zero);
    set_node_key(0, sons[0]->node_key(0));
    set_node_key(1, sons[1]->node_key(1));
    set_node_key(2, sons[2]->node_key(2));
    set_node_key(3, sons[3]->node_key(3));
    set_node_key(4, sons[0]->node_key(1));
    set_node_key(5, sons[1]->node_key(2));
    set_node_key(6, sons[2]->node_key(3));
    set_node_key(7, sons[3]->node_key(0));

    for(int idim = 0; idim < DIMENSION; idim++)
    {
        set_elm_loc(idim,sons[0]->elm_loc(idim) / 2);
    }
    set_myprocess(sons[0]->myprocess());
    generation(sons[0]->generation() - 1);
    set_material(sons[0]->material());
    set_no_of_eqns(EQUATIONS);
    
    calc_which_son();
    bcptr(sons[0]->bcptr());
    //order information -- keep the highest order
    set_order(0, sons[0]->order(0));
    if(order(0) < sons[1]->order(0))
        set_order(0, sons[1]->order(0));
    set_order(1, sons[1]->order(1));
    if(order(1) < sons[2]->order(1))
        set_order(1, sons[2]->order(1));
    set_order(2, sons[2]->order(2));
    if(order(2) < sons[3]->order(2))
        set_order(2, sons[3]->order(2));
    set_order(3, sons[3]->order(3));
    if(order(3) < sons[0]->order(3))
        set_order(3, sons[0]->order(3));
    set_order(4, sons[0]->order(4));
    for(i = 1; i < 4; i++)
        if(order(4) < sons[i]->order(4))
            set_order(4, sons[i]->order(4));
    
    set_ndof(0);
    for(i = 0; i < 4; i++)
        set_ndof(ndof() + order(i) * (no_of_eqns()));
    set_ndof(ndof() + pow((float) (order(4) - 1), 2) * (no_of_eqns()));
    
    set_refined_flag(1); // not an active element yet!!!
            
    // neighbor information
    for(ison = 0; ison < 4; ison++)
    {
        isonneigh = ison;
        ineigh = isonneigh;
        get_neigh_gen(ineigh, sons[ison]->neigh_gen(isonneigh));
        set_neighbor(ineigh, sons[ison]->neighbor(isonneigh));
        set_neigh_proc(ineigh, sons[ison]->neigh_proc(isonneigh));
        
        isonneigh = (ison + 3) % 4;
        ineigh = isonneigh + 4;
        get_neigh_gen(ineigh, sons[ison]->neigh_gen(isonneigh));
        set_neighbor(ineigh, sons[ison]->neighbor(isonneigh));
        if((sons[ison]->neigh_gen(isonneigh) == generation()) || (sons[ison]->neigh_proc(isonneigh)== -1))
            set_neigh_proc(ineigh, -2);
        else
            set_neigh_proc(ineigh, sons[ison]->neigh_proc(isonneigh));
    }
    
    /* brother information -- requires that atleast one of this
     element's neighboring brothers is on this process in 
     order to get information on the brother that is not a neighbor */
    Element* EmTemp;
    switch (which_son())
    {
        case 0:
            set_brother(0, key());
            if(neigh_proc(1) == -1)
            {
                set_brother(1, sfc_key_zero);
            }
            else if(neigh_gen(1) == generation())
            {
                set_brother(1, neighbor(1));
            }
            else if(neigh_gen(1) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(1));
                assert(EmTemp);
                set_brother(1, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            if(neigh_proc(2) == -1)
            {
                set_brother(3, sfc_key_zero);
            }
            else if(neigh_gen(2) == generation())
            {
                set_brother(3, neighbor(2));
            }
            else if(neigh_gen(2) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(2));
                assert(EmTemp);
                set_brother(3, EmTemp->father());
            }
            else
                assert(0);
            break;
        case 1:
            set_brother(1, key());
            if(neigh_proc(3) == -1)
            {
                set_brother(0, sfc_key_zero);
            }
            else if(neigh_gen(3) == generation())
            {
                set_brother(0, neighbor(3));
            }
            else if(neigh_gen(3) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(3));
                assert(EmTemp);
                set_brother(0, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            if(neigh_proc(2) == -1)
            {
                set_brother(2, sfc_key_zero);
            }
            else if(neigh_gen(2) == generation())
            {
                set_brother(2, neighbor(2));
            }
            else if(neigh_gen(2) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(2));
                assert(EmTemp);
                set_brother(2, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            break;
        case 2:
            set_brother(2, key());
            if(neigh_proc(0) == -1)
            {
                set_brother(1, sfc_key_zero);
            }
            else if(neigh_gen(0) == generation())
            {
                set_brother(1, neighbor(0));
            }
            else if(neigh_gen(0) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(0));
                assert(EmTemp);
                set_brother(1, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            if(neigh_proc(3) == -1)
            {
                set_brother(3, sfc_key_zero);
            }
            else if(neigh_gen(3) == generation())
            {
                set_brother(3, neighbor(3));
            }
            else if(neigh_gen(3) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(3));
                assert(EmTemp);
                set_brother(3, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            break;
        case 3:
            set_brother(3, key());
            if(neigh_proc(0) == -1)
            {
                set_brother(0, sfc_key_zero);
            }
            else if(neigh_gen(0) == generation())
            {
                set_brother(0, neighbor(0));
            }
            else if(neigh_gen(0) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(0));
                assert(EmTemp);
                set_brother(0, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            if(neigh_proc(1) == -1)
            {
                set_brother(2, sfc_key_zero);
            }
            else if(neigh_gen(1) == generation())
            {
                set_brother(2, neighbor(1));
            }
            else if(neigh_gen(1) == generation() + 1)
            {
                EmTemp = (Element*) El_Table->lookup(neighbor(1));
                assert(EmTemp);
                set_brother(2, EmTemp->father());
            }
            else
            {
                assert(0);
            }
            break;
    }
    
    find_positive_x_side(NodeTable);  //also inserts the coordinates
    calculate_dx(NodeTable);
    find_opposite_brother(El_Table);
    
    calc_topo_data(matprops_ptr);
    calc_gravity_vector(matprops_ptr);
    calc_d_gravity(El_Table);
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        state_vars[i] = 0.;
        prev_state_vars[i] = 0.;
        for(j = 0; j < 4; j++)
        {
            state_vars[i] += *(sons[j]->get_state_vars() + i) * 0.25;
            prev_state_vars[i] += *(sons[j]->get_prev_state_vars() + i) * 0.25;
        }
    }
    set_Awet(0.0);
    for(int ison = 0; ison < 4; ison++)
        set_Awet(Awet() + sons[ison]->Awet());
    set_Awet(Awet() * 0.25);
    
    //uninitialized flag values... will fix shortly 
    drypoint(0, 0.0);
    drypoint(1, 0.0);
    set_iwetnode(8);
    set_Swet(1.0);
    
    //calculate the shortspeed
    set_shortspeed(0.0);
    for(j = 0; j < 4; j++)
        set_shortspeed(shortspeed() + *(sons[j]->get_state_vars()) * sons[j]->shortspeed());
    if(state_vars[0] > 0.0)
       set_shortspeed(shortspeed()/(4.0 * state_vars[0]));
    
    return;
}

const SFC_Key& Element::father() const
{
    switch (which_son())
    {
        case 0:
            return node_key(2);
            break;
        case 1:
            return node_key(3);
            break;
        case 2:
            return node_key(0);
            break;
        case 3:
            return node_key(1);
            break;
    }
    cout<<"my key is "<<key()<<" in getfather on proc "<<myprocess()<<endl;
    assert(0); // 0 <= which_son <= 3 !!!
}

int Element::which_neighbor(const SFC_Key &FindNeigh)
{
    int i;
    for(i = 0; i < 8; i++)
        if((neighbor(i)==FindNeigh) && (neigh_proc(i) >= 0))
            return i;
    
    assert(i < 8);
    
    return i;
}

void Element::change_neighbor(const SFC_Key * newneighbs, int which_side, int proc, int reg)
{
    int j;
    switch (reg)
    {
        case 1:
            j = 0;
        case 3:
            assert(which_side < 4);
            set_neighbor(which_side, newneighbs[0]);
            set_neighbor(which_side + 4, newneighbs[1]);
            set_neigh_proc(which_side + 4, proc); //assuming no element movement
            get_neigh_gen(which_side, neigh_gen(which_side) + 1);
            get_neigh_gen(which_side + 4, neigh_gen(which_side) + 1);
            break;
        case 4:
            j = 0;
        case 2:
            j = 0;
        case 5:
            set_neighbor(which_side, newneighbs[0]);
            get_neigh_gen(which_side, neigh_gen(which_side) + 1);
            break;
            
        case 6:
            set_neighbor(which_side, newneighbs[0]);
            set_neighbor(which_side + 4, newneighbs[0]);
            get_neigh_gen(which_side, neigh_gen(which_side) + 1);
            get_neigh_gen(which_side + 4, neigh_gen(which_side) + 1);
            break;
            
            /*Andrew's section called from update_interproc*/
        case 10: //the refined element and old neighbor have the same gen.     
            assert(which_side < 4);
            set_neighbor(which_side, newneighbs[0]);
            set_neighbor(which_side + 4, newneighbs[1]);

            set_neigh_proc(which_side + 4, proc);
            
            get_neigh_gen(which_side, neigh_gen(which_side) + 1);
            get_neigh_gen(which_side + 4, neigh_gen(which_side) + 1);
            break;
            
        case 11:
            set_neighbor(which_side, newneighbs[0]);
            set_neighbor(which_side + 4, newneighbs[0]);
            get_neigh_gen(which_side, neigh_gen(which_side) + 1);
            get_neigh_gen(which_side + 4, neigh_gen(which_side) + 1);
            break;
            
    }
}

void Element::update_ndof()
{
    int help = 0;
    for(int i = 0; i < 4; i++)
        help += order(i) * (no_of_eqns());
    
    help += pow((float) (order(4) - 1), 2) * (no_of_eqns());
    set_ndof(help);
}

void Element::get_nelb_icon(HashTable* NodeTable, HashTable* HT_Elem_Ptr, int* Nelb, int* icon)

//for ONE step H-refinement (icon)

{
    int i;
    int ifg = 2;
    int Nc = ndof();
    double bc_value[4]; //--for poisson equ
    Node* NodePtr;
    Element* ElemPtr;
    
    for(i = 0; i < 4; i++)
    {
        Nelb[i] = 0;
        icon[i] = 0;
        //bc_value[i] = .0;
    } //-- -1 may be better
    
    /*modifid for elasticity 03.08*/
    
    for(i = 0; i < 4; i++)
    {
        if(neigh_proc(i) == -1)
        {
            
            if(bcptr() == nullptr)
                Nelb[i] = 2; //the element has no bc at all
            else
            {
                if(bcptr()->type[i] == 0)
                    Nelb[i] = 2; //no bc at that side
                            
                else if(bcptr()->type[i] == 2)
                    Nelb[i] = 1; //stress applied
                            
                else if(bcptr()->type[i] == 1 && bcptr()->value[i][0][0] == UN_CONSTRAINED)
                    Nelb[i] = 4; //y constrined
                            
                else if(bcptr()->type[i] == 1 && bcptr()->value[i][0][1] == UN_CONSTRAINED)
                    Nelb[i] = 3; //x constrined
                            
                else if(bcptr()->type[i] == 1)
                    Nelb[i] = 5; //x, y constrained
            }
        }
        
        bc_value[i] = 0;
    }
    
    if(generation()) //filling up icon array
    {
        ElemPtr = (Element*) HT_Elem_Ptr->lookup(father()); //if the father is not there it is not a constrained element
                
        if(ElemPtr)
        {
            int j = 0; //<---indicates which son it is
            while (ElemPtr->son(j) != key())
            {
                j++;
                if(j == 4)
                {
                    cerr << "error in get_el_stiffness\n\n" << flush;
                    exit(0);
                }
            }
            
            int a = j - 1;
            if(a == -1)
                a = 3;
            NodePtr = (Node*) (NodeTable->lookup(node_key(a)));
            assert(NodePtr);
            
            if(NodePtr->getinfo() == S_C_CON)
                icon[a] = -1;
            
            a = j + 1;
            if(a == 4)
                a = 0;
            NodePtr = (Node*) (NodeTable->lookup(node_key(a)));
            assert(NodePtr);
            
            if(NodePtr->getinfo() == S_C_CON)
                icon[j] = 1; //-- j was changed to a
                        
        }
    }
    
}

Element::~Element()
{
    if(bcptr())
        delete_bcptr();
    /*  if(key[0] == (unsigned) 2501998986) {
     int mmmyid;
     MPI_Comm_rank(MPI_COMM_WORLD, &mmmyid);
     printf("?????????????????????????????????????????????????????? \n");
     printf("?????????????????????????????????????????????????????? \n");
     printf("deleting element %u %u on %d $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n",key[0], key[1], mmmyid);
     printf("?????????????????????????????????????????????????????? \n");
     printf("?????????????????????????????????????????????????????? \n");
     }*/
}

/* routine also puts in the coords of the center node in the elm */
void Element::find_positive_x_side(HashTable* nodetable)
{
    int i, j, side;
    double xmax;
    Node* nodeptr;
    
    nodeptr = (Node*) (nodetable->lookup(key()));
    xmax = nodeptr->coord[0];
    set_coord(0, xmax);
    set_coord(1, nodeptr->get_coord(1));
    
    for(i = 4; i < 8; i++)
    {
        nodeptr = (Node*) (nodetable->lookup(node_key(i)));
        double xcoord = nodeptr->coord[0];
        if(xcoord > xmax)
        {
            xmax = xcoord;
            side = i - 4;
        }
    }
    
    set_positive_x_side(side);
    
    return;
}

void Element::get_slopes(HashTable* El_Table, HashTable* NodeTable, double gamma)
{
    int j = 0, bc = 0;
    /* check to see if this is a boundary */
    while (j < 4 && bc == 0)
    {
        if(neigh_proc(j) == INIT)
            bc = 1;
        j++;
    }
    if(bc == 1)
    {
        for(j = 0; j < NUM_STATE_VARS * DIMENSION; j++)
            d_state_vars[j] = 0;
        return;
    }
    
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side();
    switch (positive_x_side())
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
    }
    /* x direction */
    Element *ep = neighborPtr[xp]; //(Element*) (El_Table->lookup(&neighbor(xp)[0]));
    Element *em = neighborPtr[xm]; //(Element*) (El_Table->lookup(&neighbor(xm)[0]));
    Element *ep2 = NULL;
    Element *em2 = NULL;
    //check if element has 2 neighbors on either side
    Node* ndtemp = node_keyPtr[xp + 4]; //(Node*) NodeTable->lookup(&node_key[xp + 4][0]);
    if(ndtemp->info == S_C_CON)
    {
        ep2 = neighborPtr[xp + 4]; //(Element*) (El_Table->lookup(&neighbor[xp + 4][0]));
        assert(neigh_proc(xp + 4) >= 0 && ep2);
    }
    ndtemp = node_keyPtr[xm + 4]; //(Node*) NodeTable->lookup(&node_key[xm + 4][0]);
    if(ndtemp->info == S_C_CON)
    {
        em2 = neighborPtr[xm + 4]; //(Element*) (El_Table->lookup(&neighbor[xm + 4][0]));
        assert(neigh_proc(xm + 4) >= 0 && em2);
    }
    
    double dp, dm, dc, dxp, dxm;
    dxp = ep->coord(0) - coord(0);
    dxm = coord(0) - em->coord(0);
    for(j = 0; j < NUM_STATE_VARS; j++)
    {
        dp = (ep->state_vars[j] - state_vars[j]) / dxp;
        if(ep2 != NULL)
            dp = .5 * (dp + (ep2->state_vars[j] - state_vars[j]) / dxp);
        dm = (state_vars[j] - em->state_vars[j]) / dxm;
        if(em2 != NULL)
            dm = .5 * (dm + (state_vars[j] - em2->state_vars[j]) / dxm);
        
        dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average
        //do slope limiting
        d_state_vars[j] = .5 * (c_sgn(dp) + c_sgn(dm)) * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
    }
    
    /* y direction */
    ep = neighborPtr[yp];		//(Element*) (El_Table->lookup(&neighbor(yp)[0]));
    em = neighborPtr[ym];		//(Element*) (El_Table->lookup(&neighbor(ym)[0]));
    ep2 = NULL;
    em2 = NULL;
    //check if element has 2 neighbors on either side
    ndtemp = node_keyPtr[yp + 4];		//(Node*) NodeTable->lookup(&node_key[yp + 4][0]);
    if(ndtemp->info == S_C_CON)
    {
        ep2 = neighborPtr[yp + 4];		//(Element*) (El_Table->lookup(&neighbor[yp + 4][0]));
        assert(neigh_proc(yp + 4) >= 0 && ep2);
    }
    ndtemp = node_keyPtr[ym + 4];		//(Node*) NodeTable->lookup(&node_key[ym + 4][0]);
    if(ndtemp->info == S_C_CON)
    {
        em2 = neighborPtr[ym + 4];		//(Element*) (El_Table->lookup(&neighbor[ym + 4][0]));
        if(!(neigh_proc(ym + 4) >= 0 && em2))
        {
            printf("ym=%d neigh_proc[ym+4]=%d em2=%p\n", ym, neigh_proc(ym + 4), em2);
        }

        assert(neigh_proc(ym + 4) >= 0 && em2);
    }
    
    dxp = ep->coord(1) - coord(1);
    dxm = coord(1) - em->coord(1);
    for(j = 0; j < NUM_STATE_VARS; j++)
    {
        dp = (ep->state_vars[j] - state_vars[j]) / dxp;
        if(ep2 != NULL)
            dp = .5 * (dp + (ep2->state_vars[j] - state_vars[j]) / dxp);
        dm = (state_vars[j] - em->state_vars[j]) / dxm;
        if(em2 != NULL)
            dm = .5 * (dm + (state_vars[j] - em2->state_vars[j]) / dxm);
        
        dc = (dp * dxm + dm * dxp) / (dxm + dxp);  // weighted average 
        //do slope limiting
        d_state_vars[j + NUM_STATE_VARS] = .5 * (c_sgn(dp) + c_sgn(dm))
                                           * c_dmin1(gamma * dabs(dp), gamma * dabs(dm), dabs(dc));
    }
    
    return;
}

void Element::calculate_dx(HashTable* NodeTable)
{
    int i, j;
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side();
    switch (positive_x_side())
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
    }
    
    Node *np, *nm;
    
    np = (Node*) NodeTable->lookup(node_key(xp + 4));
    nm = (Node*) NodeTable->lookup(node_key(xm + 4));
    
    dx(0,(np->coord[0] - nm->coord[0])); /*(zeta[0]*zeta[0]+1)*/
    if(dx(0) == 0)
        printf("np %p,nm %p,dx, np_coord, nm_coord %e %e\n", np, nm, np->coord[0], nm->coord[0]);
    
    np = (Node*) NodeTable->lookup(node_key(yp + 4));
    nm = (Node*) NodeTable->lookup(node_key(ym + 4));
    
    dx(1,(np->coord[1] - nm->coord[1]) /*(zeta[1]*zeta[1]+1)*/);
    
    if(dx(1) == 0)
        printf("dy, np_coord, nm_coord %e %e\n", np->coord[1], nm->coord[1]);
    
    return;
}

void Element::insert_coord(HashTable* NodeTable)
{
    
    Node* node = (Node*) (NodeTable->lookup(key()));
    int i;
    for(i = 0; i < DIMENSION; i++)
        set_coord(i, node->get_coord(i));
    
    return;
}

double max(double x, double y)
{
    if(x >= y)
    {
        return (x);
    }
    else
    {
        return (y);
    }
}

double min(double x, double y)
{
    if(x >= y)
    {
        return (y);
    }
    else
    {
        return (x);
    }
}

// the element member function calc_wet_dry_orient() calculates the orientation of the dryline (drylineorient), the wet length (Swet), the location of the drypoint, and the location of the wetpoint... it does NOT calculate the wet area (Awet)... these quantities are used in the wetted area adjustment of fluxes. calc_wet_dry_orient() is not coded for generic element orientation, positive_x_side must be side 1.  Keith wrote this may 2007
void Element::calc_wet_dry_orient(HashTable *El_Table)
{
    
    int ifsidewet[4], numwetsides = 0;
    int ineigh;
    Element *EmTemp;
    
    for(ineigh = 0; ineigh < 4; ineigh++)
    {
        if(neigh_proc(ineigh) == -1)
            //edge of map and cell has same wetness as the cell 
            ifsidewet[ineigh] = (state_vars[0] > GEOFLOW_TINY) ? 1 : 0;
        else
        {
            EmTemp = (Element *) El_Table->lookup(neighbor(ineigh));
            if(*(EmTemp->get_state_vars() + 0) > GEOFLOW_TINY)
                //first neighbor on this side is wet
                ifsidewet[ineigh] = 1;
            else if(neigh_proc(ineigh + 4) == -2)
                //only one neighbor on this side and it's not wet
                ifsidewet[ineigh] = 0;
            else
            {
                //since first neighbor on this side is not wet, the edge has the wetness of the second neighbor on this side 
                EmTemp = (Element *) El_Table->lookup(neighbor(ineigh + 4));
                ifsidewet[ineigh] = (*(EmTemp->get_state_vars() + 0) > GEOFLOW_TINY) ? 1 : 0;
            }
        }
        numwetsides += ifsidewet[ineigh];
    }
    
    if((ifsidewet[0] == ifsidewet[2]) && (ifsidewet[1] == ifsidewet[3]))
    {
        //if opposite sides of the element are the same (both wet or dry)
        set_iwetnode(8);
        drypoint(0, 0.0);
        drypoint(1, 0.0);
        if(state_vars[0] > GEOFLOW_TINY)
        {
            set_Awet(1.0);
            set_Swet(1.0);
        }
        else
        {
            set_Awet(0.0);
            set_Swet(0.0);
        }
        //?1.0:0.0;    
    }
    else if(numwetsides == 2)
    {
        //having exactly 2 adjacent wet edges means it has a diagonal orientation
        
        set_Swet(sqrt(2.0 * ((Awet() > 0.5) ? 1.0 - Awet() : Awet()))); //edge length of small triangle
        drypoint(0, 0.5 * (1.0 - Swet()));
        drypoint(1,drypoint(0));
        if(Awet() > 0.5)
            set_Swet(1.0 - Swet()); //the small triangle is dry not wet
            
        if(ifsidewet[3] && ifsidewet[0])
        {
            set_iwetnode(0);
            if(Awet() <= 0.5){
                drypoint(0, -drypoint(0));
                drypoint(1,drypoint(0));
            }
            
        }
        else if(ifsidewet[0] && ifsidewet[1])
        {
            set_iwetnode(1);
            if(Awet() <= 0.5)
                drypoint(1, -drypoint(1));
            else
                drypoint(0, -drypoint(0));
        }
        else if(ifsidewet[1] && ifsidewet[2])
        {
            set_iwetnode(2);
            if(Awet() > 0.5)
            {
                drypoint(0, -drypoint(0));
                drypoint(1,drypoint(0));
            }
        }
        else if(ifsidewet[2] && ifsidewet[3])
        {
            set_iwetnode(3);
            if(Awet() > 0.5)
                drypoint(1, -drypoint(1));
            else
                drypoint(0, -drypoint(0));
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
        set_Swet(Awet());
        
        set_iwetnode(ineigh + 4);
        switch (iwetnode())
        {
            case 4:
                drypoint(0, 0.0);
                drypoint(1, -0.5 + Swet());
                break;
            case 5:
                drypoint(0, +0.5 - Swet());      
                drypoint(1, 0.0);
                break;
            case 6:
                drypoint(0, 0.0);
                drypoint(1, +0.5 - Swet());
                break;
            case 7:
                drypoint(0, -0.5 + Swet());
                drypoint(1, 0.0);
                break;
            default:
                assert(0);
        }
    }
    
    if(iwetnode() == 8)
        set_Awet((state_vars[0] > GEOFLOW_TINY) ? 1.0 : 0.0);
    
    return;
}

// The Element member function calc_elem_edge_wet_fraction() returns the "how much of this is wet" fraction of side that this element shares with its ineigh-th neighboring element . This fraction is used to determine wether or not to "zero" the state variables used to compute physical fluxes through a "dry" side, as part of "thin layer control".  Keith wrote this function may 2007
double Element::calc_elem_edge_wet_fraction(int ineigh, int ifusewholeside)
{
    /* double Element::calc_elem_edge_wet_fraction(unsigned *neighkey) {
     
     int ineigh=which_neighbor(neighkey);
     assert((-1<ineigh)&&(ineigh<8));
     */

    if(Awet() == 0.0)
        return 0.0;
    //return -1.0;
    
    if(Awet() == 1.0)
        return 1.0;
    
    if(iwetnode() == 8)
    {
        cout<<"calc_elem_edge_wet_fraction(): key={"<<key()<<"} adapted="<<adapted_flag()<<"\n";
        printf("  iwetnode=%d, Awet=%g, Swet=%g, drypoint={%g,%g}\n", iwetnode(), Awet(), Swet(), drypoint(0), drypoint(1));
        assert(iwetnode() != 8);
    }
    
    if(!((0.0 < Awet()) && (Awet() < 1.0)))
    {
        printf("Awet=%g\n", Awet());
        fflush(stdout);
        assert((0.0 < Awet()) && (Awet() < 1.0));
    }
    int ineighm4 = ineigh % 4;
    
    if((neigh_gen(ineighm4 + 4) == -2) || ifusewholeside)
    {
        //there is only one neighbor on this side
        switch (iwetnode())
        {
            case 0:
                switch (ineighm4)
                {
                    case 3:
                    case 0:
                        if(Awet() > 0.5)
                            return 1.0;
                        else
                            return Swet();
                    case 2:
                    case 1:
                        if(Awet() > 0.5)
                            return Swet();
                        else
                            return 0.0;
                    default:
                        assert(0);
                }
            case 1:
                switch (ineighm4)
                {
                    case 0:
                    case 1:
                        if(Awet() > 0.5)
                            return 1.0;
                        else
                            return Swet();
                    case 3:
                    case 2:
                        if(Awet() > 0.5)
                            return Swet();
                        else
                            return 0.0;
                    default:
                        assert(0);
                }
            case 2:
                switch (ineighm4)
                {
                    case 1:
                    case 2:
                        if(Awet() > 0.5)
                            return 1.0;
                        else
                            return Swet();
                    case 0:
                    case 3:
                        if(Awet() > 0.5)
                            return Swet();
                        else
                            return 0.0;
                    default:
                        assert(0);
                }
            case 3:
                switch (ineighm4)
                {
                    case 2:
                    case 3:
                        if(Awet() > 0.5)
                            return 1.0;
                        else
                            return Swet();
                    case 1:
                    case 0:
                        if(Awet() > 0.5)
                            return Swet();
                        else
                            return 0.0;
                    default:
                        assert(0);
                }
            case 4:
                switch (ineighm4)
                {
                    case 0:
                        return 1.0;
                    case 3:
                    case 1:
                        return Swet();
                    case 2:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 5:
                switch (ineighm4)
                {
                    case 1:
                        return 1.0;
                    case 0:
                    case 2:
                        return Swet();
                    case 3:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 6:
                switch (ineighm4)
                {
                    case 2:
                        return 1.0;
                    case 1:
                    case 3:
                        return Swet();
                    case 0:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 7:
                switch (ineighm4)
                {
                    case 3:
                        return 1.0;
                    case 2:
                    case 0:
                        return Swet();
                    case 1:
                        return 0.0;
                    default:
                        assert(0);
                }
            default:
                assert(0);
        }
    }
    else
    {
        //there is are 2 more refined neighbors on this side
        //therefore need to "double" the wetness (possibly 
        //minus 0.5) for each
        switch (iwetnode())
        {
            case 0:
                switch (ineigh)
                {
                    case 7:
                    case 0:
                        if(Awet() > 0.125)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 3:
                    case 4:
                        if(Awet() <= 0.125)
                            return 0.0;
                        else if(Awet() > 0.5)
                            return 1.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    case 6:
                    case 1:
                        if(Awet() <= 0.5)
                            return 0.0;
                        else if(Awet() > 0.875)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 2:
                    case 5:
                        if(Awet() <= 0.875)
                            return 0.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    default:
                        assert(0);
                }
            case 1:
                switch (ineigh)
                {
                    case 4:
                    case 1:
                        if(Awet() > 0.125)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 0:
                    case 5:
                        if(Awet() <= 0.125)
                            return 0.0;
                        else if(Awet() > 0.5)
                            return 1.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    case 7:
                    case 2:
                        if(Awet() <= 0.5)
                            return 0.0;
                        else if(Awet() > 0.875)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 3:
                    case 6:
                        if(Awet() <= 0.875)
                            return 0.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    default:
                        assert(0);
                }
            case 2:
                switch (ineigh)
                {
                    case 5:
                    case 2:
                        if(Awet() > 0.125)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 1:
                    case 6:
                        if(Awet() <= 0.125)
                            return 0.0;
                        else if(Awet() > 0.5)
                            return 1.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    case 4:
                    case 3:
                        if(Awet() <= 0.5)
                            return 0.0;
                        else if(Awet() > 0.875)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 0:
                    case 7:
                        if(Awet() <= 0.875)
                            return 0.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    default:
                        assert(0);
                }
            case 3:
                switch (ineigh)
                {
                    case 5:
                    case 2:
                        if(Awet() > 0.125)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 1:
                    case 6:
                        if(Awet() <= 0.125)
                            return 0.0;
                        else if(Awet() > 0.5)
                            return 1.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    case 4:
                    case 3:
                        if(Awet() <= 0.5)
                            return 0.0;
                        else if(Awet() > 0.875)
                            return 1.0;
                        else
                            return 2.0 * Swet();
                    case 0:
                    case 7:
                        if(Awet() <= 0.875)
                            return 0.0;
                        else
                            return 2.0 * (Swet() - 0.5);
                    default:
                        assert(0);
                }
            case 4:
                switch (ineigh)
                {
                    case 0:
                    case 4:
                        return 1.0;
                    case 7:
                    case 1:
                        if(Awet() > 0.5)
                            return 1;
                        else
                            return 2.0 * Swet();
                    case 3:
                    case 5:
                        if(Awet() > 0.5)
                            return 2.0 * (Swet() - 0.5);
                        else
                            return 0.0;
                    case 6:
                    case 2:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 5:
                switch (ineigh)
                {
                    case 1:
                    case 5:
                        return 1.0;
                    case 4:
                    case 2:
                        if(Awet() > 0.5)
                            return 1;
                        else
                            return 2.0 * Swet();
                    case 0:
                    case 6:
                        if(Awet() > 0.5)
                            return 2.0 * (Swet() - 0.5);
                        else
                            return 0.0;
                    case 7:
                    case 3:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 6:
                switch (ineigh)
                {
                    case 2:
                    case 6:
                        return 1.0;
                    case 5:
                    case 3:
                        if(Awet() > 0.5)
                            return 1;
                        else
                            return 2.0 * Swet();
                    case 1:
                    case 7:
                        if(Awet() > 0.5)
                            return 2.0 * (Swet() - 0.5);
                        else
                            return 0.0;
                    case 4:
                    case 0:
                        return 0.0;
                    default:
                        assert(0);
                }
            case 7:
                switch (ineigh)
                {
                    case 3:
                    case 7:
                        return 1.0;
                    case 6:
                    case 0:
                        if(Awet() > 0.5)
                            return 1;
                        else
                            return 2.0 * Swet();
                    case 2:
                    case 4:
                        if(Awet() > 0.5)
                            return 2.0 * (Swet() - 0.5);
                        else
                            return 0.0;
                    case 5:
                    case 1:
                        return 0.0;
                    default:
                        assert(0);
                }
            default:
                assert(0);
        }
    }
    
    assert(0);
    return 0.0;
}

// this function relaxes the zeroing of fluxes through cell edges that are completely dry at the beginning of the timestep, as indicated by calc_elem_edge_wet_fraction(), but will be at least partly wet by the end of the timestep... Keith wrote this function June 2007
double Element::calc_elem_edge_wetness_factor(int ineigh, double dt)
{
    
#ifdef DISABLE_DRY_FLUX_ZEROING
    //this disables the zeroing of fluxes through dry sides facet of thin layer control
    return 1.0;
#endif
    
    if(!(dt > 0.0))
        //is this a prank call?!! why are you bothering me with a timestep of zero?
        return 0.0;
    
    //handle completely wet or completely dry cells as a special case
    if((iwetnode() == 8) || !(state_vars[0] > GEOFLOW_TINY))
        return Awet();  //one or zero;
        
    double wetnessfactor = calc_elem_edge_wet_fraction(ineigh % 4, 1);
    if(wetnessfactor > 0.0)
        return 1.0; //for simpleness ... may change later but will need to change most of this function's logic
        
    //the logic below this point is only to compute the time averaged wetness of a cell side that is completely dry at the beginning of the timestep assuming the cell is at least partially wet (pileheight>GEOFLOW_TINY) at the beginning of the timestep... if a edge is partially or completely wet, a wetness factor of 1 has already been returned.
    
    double speed; //this is the rarefaction speed, and is positive when moving in the direction starting from the most wet node heading to the most dry node
    double dtnotwet = 0.0; //the ammount of time until the completely dry edge will be partially wet.
    
    double a; //speed of sound
    double VxVy[2]; //will be rarefaction velocity non-dimensionalized by cell size
    if(elementType() == ElementType::TwoPhases)
    {
        VxVy[0] = state_vars[2] / state_vars[1];
        if(VxVy[0] != 0.0)
        {
            a = sqrt(effect_kactxy(0) * gravity(2) * state_vars[1] + (state_vars[0] - state_vars[1]) * gravity(2));
            VxVy[0] *= (1.0 + 2.0 * a / fabs(VxVy[0])) / dx(0);
        }

        VxVy[1] = state_vars[3] / state_vars[1];
        if(VxVy[1] != 0.0)
        {
            a = sqrt(effect_kactxy(1) * gravity(2) * state_vars[1] + (state_vars[0] - state_vars[1]) * gravity(2));
            VxVy[1] *= (1.0 + 2.0 * a / fabs(VxVy[1])) / dx(1);
        }
    }
    if(elementType() == ElementType::SinglePhase)
    {
        VxVy[0] = state_vars[1] / state_vars[0];
        if(VxVy[0] != 0.0)
        {
            a = sqrt(effect_kactxy(0) * gravity(2) * state_vars[0]);
            VxVy[0] *= (1.0 + 2.0 * a / fabs(VxVy[0])) / dx(0);
        }

        VxVy[1] = state_vars[2] / state_vars[0];
        if(VxVy[1] != 0.0)
        {
            a = sqrt(effect_kactxy(1) * gravity(2) * state_vars[0]);
            VxVy[1] *= (1.0 + 2.0 * a / fabs(VxVy[1])) / dx(1);
        }
    }
    
    double doubleswap = 1.0 / sqrt(2.0);
    switch (iwetnode())
    {
        case 0:
            speed = VxVy[0] * doubleswap + VxVy[1] * doubleswap;
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.0 - drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 1:
            speed = VxVy[0] * -doubleswap + VxVy[1] * doubleswap;
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.0 + drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 2:
            speed = VxVy[0] * doubleswap + VxVy[1] * -doubleswap;
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.0 + drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 3:
            speed = VxVy[0] * -doubleswap + VxVy[1] * -doubleswap;
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.0 - drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 4:
            speed = VxVy[1];
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.5 - drypoint(1)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 5:
            speed = -VxVy[0];
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.5 + drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 6:
            speed = -VxVy[1];
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.5 + drypoint(1)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        case 7:
            speed = VxVy[0];
            if(speed <= 0.0)
                return 0.0;
            dtnotwet = (0.5 - drypoint(0)) / speed;
            if(dtnotwet >= dt)
                return 0.0;
            return 1.0 - dtnotwet / dt;
        default:
            assert(0);
    }
    
    return 0.0;
}

// The Element member function convect_dryline() calculates the coordinates of the "drypoint" in the element's local coordinate system.  This is used to determine the location of the wet-dry front (or dryline) inside this element, which in turn is used (in conjunction with the location of "iwetnode" - which indicates which side of the dryline is wet) to determine the fraction of its total area that is wet (Awet).  Awet is then returned by the function.  Keith wrote this function may 2007
double Element::convect_dryline(const double Vx, const double Vy, const double dt)
{
    
    //if dt>0.0 is used to disable Awet==0 during refinement, when split of a father's conservative variables have not yet been set, and need to be set based on the ratio's of Awet between the new sons... this is also why the function returns Awet rather than not returning any number
    if((state_vars[0] <= GEOFLOW_TINY) && (dt > 0.0))
    {
        set_Awet(0.0);
        return Awet();
    }
    
    if(iwetnode() == 8)
    {
        set_Awet(1.0);
        return Awet();
    }
    
    drypoint(0, drypoint(0) + Vx * dt / dx(0));
    drypoint(1, drypoint(1) + Vy * dt / dx(1));
    
    switch (iwetnode())
    {
        case 0: //diagonal: \
    drypoint[0]=0.5*(drypoint[0]+drypoint[1]);
            if(drypoint(0) < -0.5)
                set_Awet(0.0);
            else if(drypoint(0) > 0.5)
                set_Awet(1.0);
            else if(drypoint(0) < 0.0)
                set_Awet(2 * (0.5 + drypoint(0)) * (0.5 + drypoint(0)));
            else
                set_Awet(1.0 - 2.0 * (0.5 - drypoint(0)) * (0.5 - drypoint(0)));
            return Awet();
        case 1: //diagonal: /
            drypoint(0, 0.5 * (drypoint(0) - drypoint(1)));
            if(drypoint(0) > 0.5)
                set_Awet(0.0);
            else if(drypoint(0) < -0.5)
                set_Awet(1.0);
            else if(drypoint(0) > 0.0)
                set_Awet(2.0 * (0.5 - drypoint(0)) * (0.5 - drypoint(0)));
            else
                set_Awet(1.0 - 2 * (0.5 + drypoint(0)) * (0.5 + drypoint(0)));
            return Awet();
        case 2: //diagonal: \
            drypoint(0, 0.5*(drypoint(0)+drypoint(1)));
            if(drypoint(0) > 0.5)
                set_Awet(0.0);
            else if(drypoint(0) < -0.5)
                set_Awet(1.0);
            else if(drypoint(0) > 0.0)
                set_Awet(2.0 * (0.5 - drypoint(0)) * (0.5 - drypoint(0)));
            else
                set_Awet(1.0 - 2 * (0.5 + drypoint(0)) * (0.5 + drypoint(0)));
            return Awet();
        case 3: //diagonal: /
            drypoint(0, 0.5 * (drypoint(0) - drypoint(1)));
            if(drypoint(0) < -0.5)
                set_Awet(0.0);
            else if(drypoint(0) > 0.5)
                set_Awet(1.0);
            else if(drypoint(0) < 0.0)
                set_Awet(2 * (0.5 + drypoint(0)) * (0.5 + drypoint(0)));
            else
                set_Awet(1.0 - 2.0 * (0.5 - drypoint(0)) * (0.5 - drypoint(0)));
            return Awet();
        case 4: //horizontal: -
            if(drypoint(1) < -0.5)
                set_Awet(0.0);
            else if(drypoint(1) > 0.5)
                set_Awet(1.0);
            else
                set_Awet(0.5 + drypoint(1));
            return Awet();
        case 5: //vertical: |
            if(drypoint(0) > 0.5)
                set_Awet(0.0);
            else if(drypoint(0) < -0.5)
                set_Awet(1.0);
            else
                set_Awet(0.5 - drypoint(0));
            return Awet();
        case 6: //horizontal: -
            if(drypoint(1) > 0.5)
                set_Awet(0.0);
            else if(drypoint(1) < -0.5)
                set_Awet(1.0);
            else
                set_Awet(0.5 - drypoint(1));
            return Awet();
        case 7: //vertical: |
            if(drypoint(0) < -0.5)
                set_Awet(0.0);
            else if(drypoint(0) > 0.5)
                set_Awet(1.0);
            else
                set_Awet(0.5 + drypoint(0));
            return Awet();
        default:
            assert(0);
    }
    
    return Awet();
}

//x direction flux in current cell
void Element::xdirflux(MatProps* matprops_ptr2, double dz, double wetnessfactor, double hfv[3][MAX_NUM_STATE_VARS],
                       double hrfv[3][MAX_NUM_STATE_VARS])
{
    if(elementType() == ElementType::TwoPhases)
    {
        MatPropsTwoPhases* matprops_ptr=(MatPropsTwoPhases*)matprops_ptr2;
        int i, j;
        double a, Vel[4]; // Vel[0:1]: solid-vel, Vel[2:3]: fluid-vel
        double volf = 0.;
        double epsilon = matprops_ptr->epsilon;
        double den_frac = matprops_ptr->den_fluid / matprops_ptr->den_solid;

        //the "update flux" values (hfv) are the fluxes used to update the solution,
        // they may or may not be "reset" from their standard values based on whether
        // or not the stopping criteria is triggering a change intended to cause the flow to stop.
        if(state_vars[0] < GEOFLOW_TINY)
        {
            for(i = 0; i < 3; i++)
                for(j = 0; j < NUM_STATE_VARS; j++)
                    hfv[i][j] = 0.0;
        }
    #ifdef STOPCRIT_CHANGE_FLUX
        else if(stoppedflags()==2)
        {   
            //state variables
            hfv[0][0]=state_vars[0]+d_state_vars[0]*dz;
            hfv[0][1]=state_vars[1]+d_state_vars[1]*dz;
            for (i=2; i<NUM_STATE_VARS; i++)
            hfv[0][i]=0.0;
            if((Awet()>0.0)&&(Awet()<1.0))
            {
                hfv[0][0]*=wetnessfactor;
                hfv[0][1]*=wetnessfactor;
            }

            for (i=0; i<4; i++)
            Vel[i] = 0.;

            // a^2 = k_ap*h*ph*g(3) + h*(1-phi)*g(3)
            a=sqrt(effect_kactxy(0)*gravity(2)*hfv[0][1]
                    +(hfv[0][0]-hfv[0][1])*gravity(2));

            //fluxes
            for (i=0; i<NUM_STATE_VARS; i++)
            hfv[1][i] = 0.;

            //wave speeds
            hfv[2][0]=Vel[0]-a;
            hfv[2][1]=Vel[0];
            hfv[2][2]=Vel[0]+a;
            hfv[2][3]=Vel[2]-a;
            hfv[2][4]=Vel[2];
            hfv[2][5]=Vel[2]+a;
        }
    #endif
        else
        {
            //state variables
            for(i = 0; i < NUM_STATE_VARS; i++)
                hfv[0][i] = state_vars[i] + d_state_vars[i] * dz;

            if((0.0 < Awet()) && (Awet() < 1.0))
                for(i = 0; i < NUM_STATE_VARS; i++)
                    hfv[0][i] *= wetnessfactor;
            Vel[1] = Vel[3] = 0.; // not really, but don't need it here
            Vel[0] = hfv[0][2] / hfv[0][1];
            Vel[2] = hfv[0][4] / hfv[0][0];

            // a^2 = k_ap*h*ph*g(3) + h*(1-phi)*g(3)
            double temp = effect_kactxy(0) * hfv[0][1] * gravity(2);
            a = sqrt(temp + (hfv[0][0] - hfv[0][1]) * gravity(2));

            // get du/dy
            double dudy = (d_state_vars[NUM_STATE_VARS + 2] - d_state_vars[NUM_STATE_VARS + 1] * Vel[0]) / state_vars[1];
            double alphaxy = -c_sgn(dudy) * sin(matprops_ptr->intfrict) * effect_kactxy(0);
            double temp2 = alphaxy * hfv[0][0] * hfv[0][1] * gravity(2);
            if(hfv[0][0] > GEOFLOW_TINY)
                volf = hfv[0][1] / hfv[0][0];
            //fluxes
            hfv[1][0] = hfv[0][2] + hfv[0][4] * (1. - volf);
            hfv[1][1] = hfv[0][2];
            hfv[1][2] = hfv[0][2] * Vel[0] + 0.5 * (1. - den_frac) * temp * hfv[0][0];
            hfv[1][3] = hfv[0][3] * Vel[0] + 0.5 * (1. - den_frac) * temp2;
            hfv[1][4] = hfv[0][4] * Vel[2] + 0.5 * epsilon * hfv[0][0] * hfv[0][0] * gravity(2);
            hfv[1][5] = hfv[0][5] * Vel[2];

            //wave speeds
            hfv[2][0] = Vel[0] - a;
            hfv[2][1] = Vel[0];
            hfv[2][2] = Vel[0] + a;
            hfv[2][3] = Vel[2] - a;
            hfv[2][4] = Vel[2];
            hfv[2][5] = Vel[2] + a;
        }

        //the "refinement flux" values (hrfv) are what the flux would have
        //been if it had not been reset due to being "stopped,"
        //they are needed since refinement is based on fluxes
        //(and also pileheight gradient but that's not relevant here)
    #if defined STOPCRIT_CHANGE_FLUX || defined STOPCRIT_CHANGE_BED
        if(state_vars[0] < GEOFLOW_TINY)
        {
            for (i=0; i<3; i++)
            for (j=0; j<NUM_STATE_VARS; j++)
            hrfv[i][j] = 0.;
        }
        else
        {
            for (i=0; i<NUM_STATE_VARS; i++)
            hrfv[0][i]=state_vars[i]+d_state_vars[i]*dz;

            if((0.0<Awet())&&(Awet()<1.0))
            for (i=0; i<NUM_STATE_VARS; i++)
            hrfv[0][i]*=wetnessfactor;
    
            Vel[1]= Vel[3] = 0.; // not really, but don't need them here
            Vel[0]=hfv[0][2]/hfv[0][1];
            Vel[2]=hfv[0][4]/hfv[0][0];
            double temp = effect_kactxy(0)*hrfv[0][1]*gravity(2);
            a=sqrt(temp + (hrfv[0][0]-hrfv[0][1])*gravity(2));

            // get du/dy
            double dudy=(d_state_vars[NUM_STATE_VARS+2]-
                    d_state_vars[NUM_STATE_VARS+1]*Vel[0])/state_vars[1];
            double alphaxy=-c_sgn(dudy)*sin(matprops_ptr->intfrictang)*effect_kactxy(0);
            double temp2=alphaxy*hrfv[0][0]*hrfv[0][1]*gravity(2);
            if ( hrfv[0][0] > GEOFLOW_TINY )
            volf = hrfv[0][1]/hrfv[0][0];
            //fluxes
            hrfv[1][0]=hrfv[0][2]+hrfv[0][4]*(1.-volf);
            hrfv[1][1]=hrfv[0][2];
            hrfv[1][2]=hrfv[0][2]*Vel[0] + 0.5*(1.-den_frac)*temp*hrfv[0][0];
            hrfv[1][3]=hrfv[0][3]*Vel[0] + 0.5*(1.-den_frac)*temp2;
            hrfv[1][4]=hrfv[0][4]*Vel[2] + 0.5*epsilon*hrfv[0][0]*hrfv[0][0]*gravity(2);
            hrfv[1][5]=hrfv[0][5]*Vel[2];

            //wave speeds
            hrfv[2][0]=Vel[0]-a;
            hrfv[2][1]=Vel[0];
            hrfv[2][2]=Vel[0]+a;
            hrfv[2][0]=Vel[2]-a;
            hrfv[2][1]=Vel[2];
            hrfv[2][2]=Vel[2]+a;
        }
    #else
        for(i = 0; i < 3; i++)
            for(j = 0; j < NUM_STATE_VARS; j++)
                hrfv[i][j] = hfv[i][j];
    #endif
    }
    if(elementType() == ElementType::SinglePhase)
    {
        double speed, speed2, a, VxVy[2];

        //the "update flux" values (hfv) are the fluxes used to update the solution, they may or may not be "reset" from their standard values based on whether or not the stopping criteria is triggering a change intended to cause the flow to stop.
        if(state_vars[0] < GEOFLOW_TINY)
        {
            //printf("xdirflux case 1 ");
            hfv[0][0] = hfv[0][1] = hfv[0][2] = 0.0; //state variables
            hfv[1][0] = hfv[1][1] = hfv[1][2] = 0.0; //fluxes
            hfv[2][0] = hfv[2][1] = hfv[2][2] = 0.0; //wave speeds
        }
    #ifdef STOPCRIT_CHANGE_FLUX
        else if(stoppedflags()==2)
        {
            //printf("xdirflux case 2 ");
            //state variables
            hfv[0][0]=state_vars[0]+d_state_vars[0]*dz;
            hfv[0][1]=hfv[0][2]=0.0;
            if((0.0<Awet())&&(Awet()<1.0)) hfv[0][0]*=wetnessfactor;

            speed=0.0;
            a=sqrt(effect_kactxy(0)*gravity(2)*hfv[0][0]);

            //fluxes
            hfv[1][0]=speed*hfv[0][0];
            hfv[1][1]=speed*hfv[0][1]+0.5*hfv[0][0]*a*a;
            hfv[1][2]=speed*hfv[0][2];

            //wave speeds
            hfv[2][0]=speed-a;
            hfv[2][1]=speed;
            hfv[2][2]=speed+a;
        }
    #endif
        else
        {
            //printf("xdirflux case 3 ");
            //state variables
            hfv[0][0] = state_vars[0] + d_state_vars[0] * dz;
            hfv[0][1] = state_vars[1] + d_state_vars[1] * dz;
            hfv[0][2] = state_vars[2] + d_state_vars[2] * dz;

            if((0.0 < Awet()) && (Awet() < 1.0))
            {
                hfv[0][0] *= wetnessfactor;
                hfv[0][1] *= wetnessfactor;
                hfv[0][2] *= wetnessfactor;
            }

            //eval_velocity(dz,0.0,VxVy);
            //speed=VxVy[0];
            speed = speed2 = hfv[0][1] / hfv[0][0];

            a = sqrt(effect_kactxy(0) * gravity(2) * hfv[0][0]);

            //fluxes
            hfv[1][0] = speed * hfv[0][0];
            hfv[1][1] = speed * hfv[0][1] + 0.5 * hfv[0][0] * a * a;
            hfv[1][2] = speed * hfv[0][2];

            //wave speeds
            hfv[2][0] = speed2 - a;
            hfv[2][1] = speed2;
            hfv[2][2] = speed2 + a;
        }
        //the "refinement flux" values (hrfv) are what the flux would have been if it had not been reset due to being "stopped," they are needed since refinement is based on fluxes (and also pileheight gradient but that's not relevant here)
    #if defined STOPCRIT_CHANGE_FLUX || defined STOPCRIT_CHANGE_BED
        if(state_vars[0] < GEOFLOW_TINY)
        {
            hrfv[0][0]=hrfv[0][1]=hrfv[0][2]=0.0; //state variables
            hrfv[1][0]=hrfv[1][1]=hrfv[1][2]=0.0;//fluxes
            hrfv[2][0]=hrfv[2][1]=hrfv[2][2]=0.0;//wave speeds
        }
        else
        {
            hrfv[0][0]=state_vars[0]+d_state_vars[NUM_STATE_VARS+0]*dz;
            hrfv[0][1]=state_vars[1]+d_state_vars[NUM_STATE_VARS+1]*dz;
            hrfv[0][2]=state_vars[2]+d_state_vars[NUM_STATE_VARS+2]*dz;

            if((0.0<Awet())&&(Awet()<1.0))
            {
                hrfv[0][0]*=wetnessfactor;
                hrfv[0][1]*=wetnessfactor;
                hrfv[0][2]*=wetnessfactor;
            }

            //eval_velocity(dz,0.0,VxVy);
            //speed=VxVy[0];
            speed=speed2=hfv[0][1]/hfv[0][0];

            a=sqrt(kactxy(1)*gravity(2)*hrfv[0][0]);

            //fluxes
            hrfv[1][0]=speed*hrfv[0][0];
            hrfv[1][1]=speed*hrfv[0][1]+0.5*hrfv[0][0]*a*a;
            hrfv[1][2]=speed*hrfv[0][2];

            //wave speeds
            hrfv[2][0]=speed2-a;
            hrfv[2][1]=speed2;
            hrfv[2][2]=speed2+a;
        }
    #else
        hrfv[0][0] = hfv[0][0];
        hrfv[0][1] = hfv[0][1];
        hrfv[0][2] = hfv[0][2];

        hrfv[1][0] = hfv[1][0];
        hrfv[1][1] = hfv[1][1];
        hrfv[1][2] = hfv[1][2];

        hrfv[2][0] = hfv[2][0];
        hrfv[2][1] = hfv[2][1];
        hrfv[2][2] = hfv[2][2];
    #endif
    }
    return;
}

//y direction flux in current cell
void Element::ydirflux(MatProps* matprops_ptr2, double dz, double wetnessfactor, double hfv[3][MAX_NUM_STATE_VARS],
                       double hrfv[3][MAX_NUM_STATE_VARS])
{
    if(elementType() == ElementType::TwoPhases)
    {
        MatPropsTwoPhases* matprops_ptr=(MatPropsTwoPhases*)matprops_ptr2;
        int i, j;
        double Vel[4], a;
        double volf = 0.;
        double epsilon = matprops_ptr->epsilon;
        double den_frac = matprops_ptr->den_fluid / matprops_ptr->den_solid;
        
        //the "update flux" values (hfv) are the fluxes used to update the solution,
        // they may or may not be "reset" from their standard values based on whether
        // or not the stopping criteria is triggering a change intended to cause the flow to stop.
        if(state_vars[0] < GEOFLOW_TINY)
        {
            for(i = 0; i < 3; i++)
                for(j = 0; j < NUM_STATE_VARS; j++)
                    hfv[i][j] = 0.0; //state variables
        }
    #ifdef STOPCRIT_CHANGE_FLUX
        else if(stoppedflags()==2)
        {
            //state variables
            hfv[0][0]=state_vars[0]+d_state_vars[NUM_STATE_VARS+0]*dz;
            hfv[0][1]=state_vars[1]+d_state_vars[NUM_STATE_VARS+1]*dz;

            for (i=2; i<NUM_STATE_VARS; i++)
            hfv[0][i]=0.;
            if((0.0<Awet())&&(Awet()<1.0)) hfv[0][0]*=wetnessfactor;
            if((0.0<Awet())&&(Awet()<1.0)) hfv[0][1]*=wetnessfactor;

            double temp= effect_kactxy(1)*hfv[0][1]*gravity(2);
            a=sqrt(temp + (hfv[0][0]-hfv[0][1])*gravity(2));
            Vel[0]=Vel[1]=Vel[2]=Vel[3]=0.;

            //fluxes
            for (i=0; i<NUM_STATE_VARS; i++)
            hfv[1][i]=0.;
    
            //wave speeds
            hfv[2][0]=Vel[1]-a;
            hfv[2][1]=Vel[1];
            hfv[2][2]=Vel[1]+a;
            hfv[2][3]=Vel[3]-a;
            hfv[2][4]=Vel[3];
            hfv[2][5]=Vel[3]+a;
        }
    #endif
        else
        {
            //state variables
            for(i = 0; i < NUM_STATE_VARS; i++)
                hfv[0][i] = state_vars[i] + d_state_vars[NUM_STATE_VARS + i] * dz;

            if((0.0 < Awet()) && (Awet() < 1.0))
                for(i = 0; i < NUM_STATE_VARS; i++)
                    hfv[0][i] *= wetnessfactor;

            // a = speed of sound through the medium
            double temp = effect_kactxy(1) * hfv[0][1] * gravity(2);
            a = sqrt(temp + (hfv[0][0] - hfv[0][1]) * gravity(2));

            // velocities
            Vel[0] = Vel[2] = 0.; // don't need them here
            Vel[1] = hfv[0][3] / hfv[0][1];
            Vel[3] = hfv[0][5] / hfv[0][0];

            // hydostatic terms
            double dvdx = (d_state_vars[3] - d_state_vars[1] * Vel[1]) / state_vars[1];
            double alphayx = -c_sgn(dvdx) * sin(matprops_ptr->intfrict) * effect_kactxy(1);
            double temp2 = alphayx * hfv[0][0] * hfv[0][1] * gravity(2);
            if(hfv[0][0] > GEOFLOW_TINY)
                volf = hfv[0][1] / hfv[0][0];

            //fluxes
            hfv[1][0] = hfv[0][3] + hfv[0][5] * (1. - volf);
            hfv[1][1] = hfv[0][3];
            hfv[1][2] = hfv[0][2] * Vel[1] + 0.5 * (1. - den_frac) * temp2;
            hfv[1][3] = hfv[0][3] * Vel[1] + 0.5 * (1. - den_frac) * temp * hfv[0][0];
            hfv[1][4] = hfv[0][4] * Vel[3];
            hfv[1][5] = hfv[0][5] * Vel[3] + 0.5 * epsilon * hfv[0][0] * hfv[0][0] * gravity(2);

            //wave speeds
            hfv[2][0] = Vel[1] - a;
            hfv[2][1] = Vel[1];
            hfv[2][2] = Vel[1] + a;
            hfv[2][3] = Vel[3] - a;
            hfv[2][4] = Vel[3];
            hfv[2][5] = Vel[3] + a;
        }
        
        //the "refinement flux" values (hrfv) are what the flux would have been if it
        // had not been reset due to being "stopped," they are needed since refinement
        // is based on fluxes (and also pileheight gradient but that's not relevant here)
    #if defined STOPCRIT_CHANGE_FLUX || defined STOPCRIT_CHANGE_BED
        if(state_vars[0] < GEOFLOW_TINY)
        for (i=0; i<3; i++)
        for (j=0; j<NUM_STATE_VARS; j++)
        hrfv[i][j]=0.;
        else
        {
            for (i=0; i<NUM_STATE_VARS; i++)
            hrfv[0][i]=state_vars[i]+d_state_vars[NUM_STATE_VARS+i]*dz;
    
            if((0.0<Awet())&&(Awet()<1.0))
            for (i=0; i<NUM_STATE_VARS; i++)
            hrfv[0][i]*=wetnessfactor;
    
            double temp=effect_kactxy(1)*hrfv[0][1]*gravity(2);
            a=sqrt(temp + (hrfv[0][0]-hrfv[0][1])*gravity(2));
    
            // velocities
            Vel[0]=Vel[2]=0.;
            Vel[1]=hrfv[0][3]/hrfv[0][1];
            Vel[3]=hrfv[0][5]/hrfv[0][0];
    
            // hydostatic terms
            double dvdx=(d_state_vars[3]-d_state_vars[1]*Vel[1])/state_vars[1];
            double alphayx=-c_sgn(dudy)*sin(matprops_ptr->intfrictang)*effect_kactxy(1);
            double temp2=alphayx*hrfv[0][0]*hrfv[0][1]*gravity(2);
            if ( hrfv[0][0] > GEOFLOW_TINY )
            volf = hrfv[0][1]/hrfv[0][0];
            //fluxes
            hrfv[1][0]=hrfv[0][3]+hrfv[0][5]*(1.-volf);
            hrfv[1][1]=hrfv[0][3];
            hrfv[1][2]=hrfv[0][2]*Vel[1] + 0.5*(1.-den_frac)*temp2;
            hrfv[1][3]=hrfv[0][3]*Vel[1] + 0.5*(1.-den_frac)*temp*hrfv[0][0];
            hrfv[1][4]=hrfv[0][4]*Vel[3];
            hrfv[1][5]=hrfv[0][5]*Vel[3] + 0.5*epsilon*hrfv[0][0]*hrfv[0][0]*gravity(2);
    
            //wave speeds
            hrfv[2][0]=Vel[1]-a;
            hrfv[2][1]=Vel[1];
            hrfv[2][2]=Vel[1]+a;
            hrfv[2][3]=Vel[3]-a;
            hrfv[2][4]=Vel[3];
            hrfv[2][5]=Vel[3]+a;
        }
    #else
        for(i = 0; i < 3; i++)
            for(j = 0; j < NUM_STATE_VARS; j++)
                hrfv[i][j] = hfv[i][j];
    #endif
    }
    if(elementType() == ElementType::SinglePhase)
    {
        double speed, speed2, a, VxVy[2];
        
        //the "update flux" values (hfv) are the fluxes used to update the solution, they may or may not be "reset" from their standard values based on whether or not the stopping criteria is triggering a change intended to cause the flow to stop.
        if(state_vars[0] < GEOFLOW_TINY)
        {
            //    printf("choice 1 ");
            hfv[0][0] = hfv[0][1] = hfv[0][2] = 0.0; //state variables
            hfv[1][0] = hfv[1][1] = hfv[1][2] = 0.0; //fluxes
            hfv[2][0] = hfv[2][1] = hfv[2][2] = 0.0; //wave speeds

        }
    #ifdef STOPCRIT_CHANGE_FLUX
        else if(stoppedflags()==2)
        {
            //    printf("choice 2 ");

            //state variables
            hfv[0][0]=state_vars[0]+d_state_vars[NUM_STATE_VARS+0]*dz;
            hfv[0][1]=hfv[0][2]=0.0;
            if((0.0<Awet())&&(Awet()<1.0)) hfv[0][0]*=wetnessfactor;

            speed=0.0;
            a=sqrt(effect_kactxy(1)*gravity(2)*hfv[0][0]);

            //fluxes
            hfv[1][0]=speed*hfv[0][0];
            hfv[1][1]=speed*hfv[0][1];
            hfv[1][2]=speed*hfv[0][2]+0.5*hfv[0][0]*a*a;

            //wave speeds
            hfv[2][0]=speed-a;
            hfv[2][1]=speed;
            hfv[2][2]=speed+a;
        }
    #endif
        else
        {
            //    printf("choice 3 ");
            //state variables
            hfv[0][0] = state_vars[0] + d_state_vars[NUM_STATE_VARS + 0] * dz;
            hfv[0][1] = state_vars[1] + d_state_vars[NUM_STATE_VARS + 1] * dz;
            hfv[0][2] = state_vars[2] + d_state_vars[NUM_STATE_VARS + 2] * dz;

            if((0.0 < Awet()) && (Awet() < 1.0))
            {
                hfv[0][0] *= wetnessfactor;
                hfv[0][1] *= wetnessfactor;
                hfv[0][2] *= wetnessfactor;
            }

            //eval_velocity(0.0,dz,VxVy);
            //speed=VxVy[1];
            speed = speed2 = hfv[0][2] / hfv[0][0];

            a = sqrt(effect_kactxy(1) * gravity(2) * hfv[0][0]);

            //fluxes
            hfv[1][0] = speed * hfv[0][0];
            hfv[1][1] = speed * hfv[0][1];
            hfv[1][2] = speed * hfv[0][2] + 0.5 * hfv[0][0] * a * a;

            //wave speeds
            hfv[2][0] = speed2 - a;
            hfv[2][1] = speed2;
            hfv[2][2] = speed2 + a;
        }
        
        //the "refinement flux" values (hrfv) are what the flux would have been if it had not been reset due to being "stopped," they are needed since refinement is based on fluxes (and also pileheight gradient but that's not relevant here)
    #if defined STOPCRIT_CHANGE_FLUX || defined STOPCRIT_CHANGE_BED
        if(state_vars[0] < GEOFLOW_TINY)
        {
            hrfv[0][0]=hrfv[0][1]=hrfv[0][2]=0.0; //state variables
            hrfv[1][0]=hrfv[1][1]=hrfv[1][2]=0.0;//fluxes
            hrfv[2][0]=hrfv[2][1]=hrfv[2][2]=0.0;//wave speeds
        }
        else
        {   
            hrfv[0][0]=state_vars[0]+d_state_vars[NUM_STATE_VARS+0]*dz;
            hrfv[0][1]=state_vars[1]+d_state_vars[NUM_STATE_VARS+1]*dz;
            hrfv[0][2]=state_vars[2]+d_state_vars[NUM_STATE_VARS+2]*dz;

            if((0.0<Awet())&&(Awet()<1.0))
            {
                hrfv[0][0]*=wetnessfactor;
                hrfv[0][1]*=wetnessfactor;
                hrfv[0][2]*=wetnessfactor;
            }

            //eval_velocity(0.0,dz,VxVy);
            //speed=VxVy[1];
            speed=speed2=hfv[0][2]/hfv[0][0];

            a=sqrt(kactxy(1)*gravity(2)*hrfv[0][0]);

            //fluxes
            hrfv[1][0]=speed*hrfv[0][0];
            hrfv[1][1]=speed*hrfv[0][1];
            hrfv[1][2]=speed*hrfv[0][2]+0.5*hrfv[0][0]*a*a;

            //wave speeds
            hrfv[2][0]=speed2-a;
            hrfv[2][1]=speed2;
            hrfv[2][2]=speed2+a;
        }
    #else
        hrfv[0][0] = hfv[0][0];
        hrfv[0][1] = hfv[0][1];
        hrfv[0][2] = hfv[0][2];

        hrfv[1][0] = hfv[1][0];
        hrfv[1][1] = hfv[1][1];
        hrfv[1][2] = hfv[1][2];

        hrfv[2][0] = hfv[2][0];
        hrfv[2][1] = hfv[2][1];
        hrfv[2][2] = hfv[2][2];
    #endif
    }
    
    return;
}

//note z is not "z" but either x or y
void Element::zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int order_flag, int dir,
                       double hfv[3][MAX_NUM_STATE_VARS], double hrfv[3][MAX_NUM_STATE_VARS], Element *EmNeigh, double dt)
{
    double dz = 0.0;
    
    int ineigh = which_neighbor(EmNeigh->key());
    if(!((-1 < ineigh) && (ineigh < 8)))
    {
        printf("zdirflux: ineigh=%d, dir=%d\n", ineigh, dir);
        printf("this element******************************\n");
        ElemBackgroundCheck2(El_Table, NodeTable, this, stdout);
        fflush(stdout);
        printf("Neigh element******************************\n");
        ElemBackgroundCheck2(El_Table, NodeTable, EmNeigh, stdout);
        fflush(stdout);
        assert((-1 < ineigh) && (ineigh < 8));
    }
    
    double wetnessfactor = calc_elem_edge_wetness_factor(ineigh, dt);
    //double thissideSwet=(calc_elem_edge_wet_fraction(ineigh%4,1)>0.0)?1.0:0.0;
    
    if(order_flag == 2)
        dz = (1.0 + dir % 2 - dir) * 0.5 * dx(dir % 2); //+ or - 1/2 dx or dy
             
    if(dir % 2 == 0)
        xdirflux(matprops_ptr, dz, wetnessfactor, hfv, hrfv);
    else if(dir % 2 == 1)
        ydirflux(matprops_ptr, dz, wetnessfactor, hfv, hrfv);
    else
    {
        printf("zdirflux: direction %d not known\n", dir);
        exit(1);
    }
    return;
}

//need move this to step.C
void riemannflux(const ElementType elementType,double hfvl[3][MAX_NUM_STATE_VARS], double hfvr[3][MAX_NUM_STATE_VARS], double flux[MAX_NUM_STATE_VARS])
{
    //hfv: h=state variable, f=flux, v=wave speeds
    //l="left" (the minus side), r="right" (the plus side)
    
    int ivar;
    //this is the hll riemann flux
    if((hfvl[0][0] == 0.0) && (hfvr[0][0] == 0.0))
        for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            flux[ivar] = 0.0;
    else
    {
        double sl, sr;
        if(elementType == ElementType::TwoPhases)
        {
            if(hfvl[0][0] == 0.0)
            {
                sl = min(0, min(2.0 * hfvr[2][0] - hfvr[2][1], 2.0 * hfvr[2][3] - hfvr[2][4]));
                sr = max(0, max(2.0 * hfvr[2][2] - hfvr[2][1], 2.0 * hfvr[2][5] - hfvr[2][4]));
            }
            else if(hfvr[0][0] == 0.0)
            {
                sl = min(0, min(2.0 * hfvl[2][0] - hfvl[2][1], 2.0 * hfvl[2][3] - hfvl[2][4]));
                sr = max(0, max(2.0 * hfvl[2][2] - hfvl[2][1], 2.0 * hfvl[2][5] - hfvl[2][4]));
            }
            else
            {
                sl = min(0, min(min(hfvl[2][0], hfvl[2][3]), min(hfvr[2][0], hfvr[2][3])));
                sr = max(0, max(max(hfvl[2][2], hfvl[2][5]), max(hfvr[2][2], hfvr[2][5])));
            }
        }
        if(elementType == ElementType::SinglePhase)
        {
            if(hfvl[0][0] == 0.0)
            {
                sl = min(0, 2.0 * hfvr[2][0] - hfvr[2][1]);
                sr = max(0, 2.0 * hfvr[2][2] - hfvr[2][1]);
                //printf("hfvr=\n%12.6g   %12.6g   %12.6g\n%12.6g   %12.6g   %12.6g\n%12.6g   %12.6g   %12.6g\nsl=%12.6g   sr=%12.6g\n",hfvr[0][0],hfvr[1][0],hfvr[2][0],hfvr[0][1],hfvr[1][1],hfvr[2][1],hfvr[0][2],hfvr[1][2],hfvr[2][2],sl,sr);
            }
            else if(hfvr[0][0] == 0.0)
            {
                sl = min(0, 2.0 * hfvl[2][0] - hfvl[2][1]);
                sr = max(0, 2.0 * hfvl[2][2] - hfvl[2][1]);
                //printf("hfvl=\n%12.6g   %12.6g   %12.6g\n%12.6g   %12.6g   %12.6g\n%12.6g   %12.6g   %12.6g\nsl=%12.6g   sr=%12.6g\n",hfvl[0][0],hfvl[1][0],hfvl[2][0],hfvl[0][1],hfvl[1][1],hfvl[2][1],hfvl[0][2],hfvl[1][2],hfvl[2][2],sl,sr);

            }
            else
            {
                sl = min(0, min(hfvl[2][0], hfvr[2][0]));
                sr = max(0, max(hfvl[2][2], hfvr[2][2]));
            }
        }
        
        if(sl >= 0.0)
            for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                flux[ivar] = hfvl[1][ivar];
        else if(sr <= 0.0)
            for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                flux[ivar] = hfvr[1][ivar];
        else
            for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                flux[ivar] = (sr * hfvl[1][ivar] - sl * hfvr[1][ivar] + sl * sr * (hfvr[0][ivar] - hfvl[0][ivar]))
                        / (sr - sl);
    }
    
    return;
}

void Element::calc_edge_states(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int myid, double dt,
                               int* order_flag, double *outflow)
{
    Node *np, *np1, *np2, *nm, *nm1, *nm2;
    Element *elm1, *elm2;
    int side, zp, zm;
    int zp2, zm2; //positive_z_side_2 minus_z_side_2
    int zelmpos = -100, zelmpos_2 = -100;
    int ivar;
    double hfv[3][MAX_NUM_STATE_VARS], hfv1[3][MAX_NUM_STATE_VARS], hfv2[3][MAX_NUM_STATE_VARS]; //update flux
    double hrfv[3][MAX_NUM_STATE_VARS], hrfv1[3][MAX_NUM_STATE_VARS], hrfv2[3][MAX_NUM_STATE_VARS]; //refinement flux
            
    //ghost elements don't have nodes so you have to make temp storage for flux
    double ghostflux[MAX_NUM_STATE_VARS]; //, (*fluxptr)[MAX_NUM_STATE_VARS];
    
    //double min(double x,double y), max(double x,double y);
    
    *outflow = 0.0;
    
    for(side = 0; side < 2; side++)
    {
        zp = (positive_x_side() + side) % 4;
        zm = (zp + 2) % 4;
        
        np = node_keyPtr[zp + 4]; //(Node*) NodeTable->lookup(&node_key[zp + 4][0]);
        
        if(neigh_proc(zp) == -1)
        {
            nm = node_keyPtr[zm + 4]; //(Node*) NodeTable->lookup(&node_key[zm + 4][0]);
            *outflow += (nm->flux[0]) * dx(!side);
            
            //outflow boundary conditions
            for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
            {
                np->flux[ivar] = nm->flux[ivar];
                np->refinementflux[ivar] = nm->refinementflux[ivar];
            }
        }
        else if(neigh_proc(zp) != myid)
        {
            np = node_keyPtr[zp + 4]; //(Node*) NodeTable->lookup(&node_key[zp + 4][0]);
            elm1 = neighborPtr[zp]; //(Element*) El_Table->lookup(&neighbor(zp)[0]);
            assert(elm1);
            
            zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
            elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this, dt);
            
            riemannflux(elm1->elementType(),hfv, hfv1, np->flux);
            riemannflux(elm1->elementType(),hrfv, hrfv1, np->refinementflux);
            
            elm2 = neighborPtr[zp + 4]; //(Element*) El_Table->lookup(&neighbor[zp + 4][0]);
            assert(elm2);
            zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
            elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this, dt);
            
            //note a rectangular domain ensures that neigh_proc[zm+4]!=-1
            if(neigh_proc(zp + 4) == myid)
            {
                zm2 = elm2->which_neighbor(key()) % 4;
                nm2 = elm2->node_keyPtr[zm2 + 4]; //(Node*) NodeTable->lookup(&elm2->node_key[zm2 + 4][0]);
                
                riemannflux(elm2->elementType(),hfv, hfv2, nm2->flux);
                riemannflux(elm2->elementType(),hrfv, hrfv2, nm2->refinementflux);
                
                for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                {
                    np->flux[ivar] = 0.5 * (np->flux[ivar] + nm2->flux[ivar]);
                    np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + nm2->refinementflux[ivar]);
                }
            }
            else
            {
                riemannflux(elm2->elementType(),hfv, hfv2, ghostflux);
                for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    np->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);
                
                riemannflux(elm2->elementType(),hrfv, hrfv2, ghostflux);
                for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    np->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
            }
            
        }
        else
        {
            
            np = node_keyPtr[zp + 4]; //(Node*) NodeTable->lookup(&node_key[zp + 4][0]);
            elm1 = neighborPtr[zp]; //(Element*) El_Table->lookup(&neighbor(zp)[0]);
            assert(elm1);
            
            zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm1, dt);
            elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv1, hrfv1, this, dt);
            
            riemannflux(elm1->elementType(),hfv, hfv1, np->flux);
            riemannflux(elm1->elementType(),hrfv, hrfv1, np->refinementflux);
            
            /* CASE I
             ------------------- -------------------               
             |                   |                   |               
             |                   |                   |               
             |                   |                   |               
             |                   |                   |               
             |                   |                   |       
             |      this       np|         elm1      |           
             |        h          |         hp        |        
             |    kactxy_gz      |    kactxy_gz_n    |        
             |                   |                   |   
             |                   |                   |       
             |                   |                   |       
             ------------------- -------------------
             */

            /*
             Case II

             --------- ----------------------------         
             |         |         |                  |       
             |positive_z_side--->|<----zelmpos      |       
             |         | this  np|                  |       
             |         |   h     |                  |       
             |         |         |                  |
             ---------|---------|nm1   elm1        |
             |         |         |      hp          |
             |         |         |                  |
             |         | elm2 np2|                  |
             |         |         |                  |
             positive_z_side_2-->|                  |
             --------- ----------------------------

             */

            if(np->info == S_S_CON)
            {
                nm1 = NULL;
                np2 = NULL;
                
                zelmpos = elm1->which_neighbor(key());
                assert(zelmpos > -1);
                nm1 = elm1->node_keyPtr[zelmpos % 4 + 4]; //(Node*) NodeTable->lookup(&elm1->node_key[zelmpos % 4 + 4][0]);
                
                elm2 = elm1->neighborPtr[(zelmpos + 4) % 8]; //(Element*) El_Table->lookup(&elm1->neighbor[(zelmpos + 4) % 8][0]);
                assert(elm2);
                
                elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, elm2, dt);
                elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, elm1, dt);
                
                if(elm1->neigh_proc((zelmpos + 4) % 8) == myid)
                {
                    zp2 = elm2->which_neighbor(elm1->key()) % 4;
                    np2 = elm2->node_keyPtr[zp2 + 4]; //(Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);
                    riemannflux(elm2->elementType(),hfv2, hfv1, np2->flux);
                    riemannflux(elm2->elementType(),hrfv2, hrfv1, np2->refinementflux);
                    
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    {
                        nm1->flux[ivar] = 0.5 * (np->flux[ivar] + np2->flux[ivar]);
                        nm1->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + np2->refinementflux[ivar]);
                    }
                }
                else
                {
                    
                    riemannflux(elm2->elementType(),hfv2, hfv1, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                        nm1->flux[ivar] = 0.5 * (np->flux[ivar] + ghostflux[ivar]);
                    
                    riemannflux(elm2->elementType(),hrfv2, hrfv1, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                        nm1->refinementflux[ivar] = 0.5 * (np->refinementflux[ivar] + ghostflux[ivar]);
                }
            }
            
            /*  Case III

             ------------------- --------- ---------               
             |                   |         |         |
             |positive_z_side--->|<----zelmpos_2     |               
             |                   |nm2 elm2 |         |               
             |                   |     hp2 |         |               
             |                   |         |         |       
             |       this        |---------|---------               
             |        h          |         |         |        
             |                   |         |         |        
             |                   |nm1 elm1 |         |        
             |                   |     hp1 |         |       
             |                   |<----zelmpos       |       
             ------------------- --------- ---------


             */

            else if(np->info == S_C_CON)
            {
                
                nm1 = NULL;
                nm2 = NULL;
                
                zelmpos = elm1->which_neighbor(key()) % 4;
                nm1 = elm1->node_keyPtr[zelmpos + 4]; //(Node*) NodeTable->lookup(&elm1->node_key[zelmpos + 4][0]);
                
                elm2 = neighborPtr[zp + 4]; //(Element*) (El_Table->lookup(&neighbor[zp + 4][0]));
                assert(elm2);
                
                zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv, hrfv, elm2, dt);
                elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv2, hrfv2, this, dt);
                
                if(neigh_proc(zp + 4) == myid)
                {
                    zelmpos_2 = elm2->which_neighbor(key()) % 4;
                    nm2 = elm2->node_keyPtr[zelmpos_2 + 4]; //(Node*) NodeTable->lookup(&elm2->node_key[zelmpos_2 + 4][0]);
                    riemannflux(elm2->elementType(),hfv, hfv2, nm2->flux);
                    riemannflux(elm2->elementType(),hrfv, hrfv2, nm2->refinementflux);
                    
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    {
                        nm1->flux[ivar] = np->flux[ivar];
                        np->flux[ivar] = 0.5 * (nm1->flux[ivar] + nm2->flux[ivar]);
                        
                        nm1->refinementflux[ivar] = np->refinementflux[ivar];
                        np->refinementflux[ivar] = 0.5 * (nm1->refinementflux[ivar] + nm2->refinementflux[ivar]);
                    }
                }
                else
                {
                    riemannflux(elm2->elementType(),hfv, hfv2, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    {
                        nm1->flux[ivar] = np->flux[ivar];
                        np->flux[ivar] = 0.5 * (nm1->flux[ivar] + ghostflux[ivar]);
                    }
                    
                    riemannflux(elm2->elementType(),hrfv, hrfv2, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    {
                        nm1->refinementflux[ivar] = np->refinementflux[ivar];
                        np->refinementflux[ivar] = 0.5 * (nm1->refinementflux[ivar] + ghostflux[ivar]);
                    }
                }
                
            }
            
        }
        
        if(neigh_proc(zm) != myid)
        {
            
            if(neigh_proc(zm) == -1)
            {
                np = node_keyPtr[zp + 4]; //(Node*) NodeTable->lookup(&node_key[zp + 4][0]);
                nm = node_keyPtr[zm + 4]; //(Node*) NodeTable->lookup(&node_key[zm + 4][0]);
                *outflow -= (np->flux[0]) * dx(!side);
                //outflow boundary conditions
                for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                {
                    nm->flux[ivar] = np->flux[ivar];
                    nm->refinementflux[ivar] = np->refinementflux[ivar];
                }
            }
            else
            {
                /* if an interface is on the x-minus or y-minus side, need to 
                 calculate those edgestates in this element */
                // x-minus side
                /*
                 
                 interface   
                 |
                 |
                 left cells-GHOST CELLS   |     
                 (no need to        |
                 calculate         |
                 fluxes )         |
                 v
                 
                 
                 Case I

                 ------------------- -------------------               
                 |                   |                   |               
                 |                   |<--zm              |               
                 |                   |                   |               
                 |                   |                   |               
                 |                   |                   |       
                 |       elm1        |nm      this       |           
                 |        h          |         hp        |        
                 |                   |                   |        
                 |                   |                   |   
                 |                   |                   |       
                 |                   |                   |       
                 ------------------- -------------------


                 Case II  

                 --------- -----------------------------         
                 |         |         |                   |       
                 |         |         |<--zm(z-minus side)|       
                 |         |  elm1   |                   |       
                 |         |   h     |                   |       
                 |         |         |                   |
                 ---------|---------|nm    this         |
                 |         |         |       hp          |
                 |         |         |                   |
                 |         | elm2    |                   |
                 |         |   h2    |                   |
                 |         |         |                   |
                 --------- -----------------------------

                 */

                nm = node_keyPtr[zm + 4]; //(Node*) NodeTable->lookup(&node_key[zm + 4][0]);
                elm1 = neighborPtr[zm]; //(Element*) El_Table->lookup(&neighbor(zm)[0]);
                assert(elm1);
                
                zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm1, dt);
                elm1->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv1, hrfv1, this, dt);
                riemannflux(elm1->elementType(),hfv1, hfv, nm->flux);
                riemannflux(elm1->elementType(),hrfv1, hrfv, nm->refinementflux);
                
                elm2 = neighborPtr[zm + 4]; //(Element*) El_Table->lookup(&neighbor[zm + 4][0]);
                assert(elm2);
                
                zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side + 2, hfv, hrfv, elm2, dt);
                elm2->zdirflux(El_Table, NodeTable, matprops_ptr, *order_flag, side, hfv2, hrfv2, this, dt);
                
                //note a rectangular domain ensures that neigh_proc[zm+4]!=-1
                if(neigh_proc(zm + 4) == myid)
                {
                    zp2 = elm2->which_neighbor(key()) % 4;
                    np2 = elm2->node_keyPtr[zp2 + 4]; //(Node*) NodeTable->lookup(&elm2->node_key[zp2 + 4][0]);
                    
                    riemannflux(elm2->elementType(),hfv2, hfv, np2->flux);
                    riemannflux(elm2->elementType(),hrfv2, hrfv, np2->refinementflux);
                    
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                    {
                        nm->flux[ivar] = 0.5 * (nm->flux[ivar] + np2->flux[ivar]);
                        nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + np2->refinementflux[ivar]);
                    }
                }
                else
                {
                    riemannflux(elm2->elementType(),hfv2, hfv, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                        nm->flux[ivar] = 0.5 * (nm->flux[ivar] + ghostflux[ivar]);
                    
                    riemannflux(elm2->elementType(),hrfv2, hrfv, ghostflux);
                    for(ivar = 0; ivar < NUM_STATE_VARS; ivar++)
                        nm->refinementflux[ivar] = 0.5 * (nm->refinementflux[ivar] + ghostflux[ivar]);
                }
                
            }
        }
        
    }
    
    return;
}
#if 0
void Element::correct(HashTable* NodeTable, HashTable* El_Table, double dt, MatProps* matprops_ptr,
                      FluxProps *fluxprops, TimeProps *timeprops, double *forceint, double *forcebed, double *eroded,
                      double *deposited)
{
    int i, j, k;
    double dtdx = dt / dx(0);
    double dtdy = dt / dx(1);
    
    double tiny = GEOFLOW_TINY;
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side();
    switch (positive_x_side())
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
    }
    
    double fluxxp[3], fluxyp[3], fluxxm[3], fluxym[3];
    
    Node* nxp = (Node*) NodeTable->lookup(node_key(xp + 4));
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        fluxxp[i] = nxp->flux[i];
    }
    Node* nyp = (Node*) NodeTable->lookup(node_key(yp + 4));
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        fluxyp[i] = nyp->flux[i];
    }
    Node* nxm = (Node*) NodeTable->lookup(node_key(xm + 4));
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        fluxxm[i] = nxm->flux[i];
    }
    Node* nym = (Node*) NodeTable->lookup(node_key(ym + 4));
    for(i = 0; i < NUM_STATE_VARS; i++)
    {
        fluxym[i] = nym->flux[i];
    }
    
#ifdef SUNOS
    /* the values being passed to correct are for a SINGLE element, NOT a
     region, as such the only change that having variable bedfriction 
     requires is to pass the bedfriction angle for the current element 
     rather than the only bedfriction 

     I wonder if this is legacy code, it seems odd that it is only called
     for the SUN Operating System zee ../geoflow/correct.f */

#ifdef DO_EROSION
    int do_erosion=1;
#else
    int do_erosion = 0;
#endif
    
#ifdef STOPCRIT_CHANGE_SOURCE
    int IF_STOPPED=stoppedflags();
#else
    int IF_STOPPED = !(!stoppedflags());
#endif
    
    double VxVy[2];
    if(state_vars[0] > GEOFLOW_TINY)
    {
        VxVy[0] = state_vars[1] / state_vars[0];
        VxVy[1] = state_vars[2] / state_vars[0];
    }
    else
    {
        VxVy[0] = VxVy[1] = 0.0;
    }
    
    correct_(state_vars, prev_state_vars, fluxxp, fluxyp, fluxxm, fluxym, &tiny, &dtdx, &dtdy, &dt, d_state_vars,
             (d_state_vars + NUM_STATE_VARS), &(zeta_ref(0)), &(zeta_ref(1)), curvature, &(matprops_ptr->intfrict),
             effect_bedfrict_ptr(), gravity, effect_kactxy, d_gravity, &(matprops_ptr->frict_tiny), forceint, forcebed,
             &do_erosion, eroded, VxVy, //eval_velocity(0.0,0.0,VxVy),
             &IF_STOPPED, Influx);
    
    *forceint *= dx(0) * dx(1);
    *forcebed *= dx(0) * dx(1);
    *eroded *= dx(0) * dx(1);
    
    //convect the dry line
    //eval_velocity(0.0,0.0,VxVy);
    if(state_vars[0] > GEOFLOW_TINY)
    {
        VxVy[0] = state_vars[1] / state_vars[0];
        VxVy[1] = state_vars[2] / state_vars[0];
    }
    else
    {
        VxVy[0] = VxVy[1] = 0.0;
    }
    convect_dryline(VxVy[0], VxVy[1], dt);
    //convect_dryline(VxVy,0.5*dt);
#endif 
    
    if(stoppedflags() == 2)
        *deposited = state_vars[0] * dx(0) * dx(1);
    else
        *deposited = 0.0;
    
    if(stoppedflags())
        *eroded = 0.0;
    
    calc_shortspeed(1.0 / dt);
    
    return;
}
#endif

void Element::calc_shortspeed(double inv_dt)
{
    
    /****************************************************************/
    /* calculate the short cell magnitude of velocity as            */
    /* dhv/dh=v*dh/dh+h*dv/dh=>lim h->0=>v                          */
    /* v=dhv/dh=(Dhv/Dt)/(Dh/Dt)   (material derivative)            */
    /* ||v||^2= (dhvx/dh)^2+(dhvy/dh)^2                             */
    /* ||v||^2*(Dh/Dt)^2-(Dhvx/Dt)^2+(Dhvy/Dt)^2 = 0 = f            */
    /* ||v||^2*(dh/dt+||v||((hv/||hv||).grad(h)))  -                */
    /*    (dhvx/dt+||v||((hv/||hv||).grad(hvx)))^2 +                */
    /*    (dhvy/dt+||v||((hv/||hv||).grad(hvy)))^2 = 0 = f          */
    /* Newton iterate to find ||v||, starting at ||hv||/h           */
    /* ||v||_(i+1)=||v||_i-((d(f^2)/d||v||)/(d^2(f^2)/d||v||^2))_i  */
    /****************************************************************/

    set_shortspeed(0.0);
    
    if(SHORTSPEED)
    {
        if(state_vars[0] > GEOFLOW_TINY)
        {
            double Vmag;
            if(elementType() == ElementType::TwoPhases)
            {
                Vmag=sqrt(state_vars[2]*state_vars[2]+
                    state_vars[3]*state_vars[3]);
            }
            if(elementType() == ElementType::SinglePhase)
            {
                Vmag = pow(state_vars[1] * state_vars[1] + state_vars[2] * state_vars[2], 0.5);
            }

            if(!(Vmag > 0.0))
                return;

            double invnormhv = 1.0 / Vmag;

            if(elementType() == ElementType::TwoPhases)
            {
                Vmag/=state_vars[1];
            }
            if(elementType() == ElementType::SinglePhase)
            {
                Vmag /= state_vars[0];
            }

            double Vmag0 = Vmag; //tall cell speed
            assert(Vmag0 > 0.0);

            double doubleswap_h;
            double doubleswap_hvx;
            double doubleswap_hvy;
            if(elementType() == ElementType::TwoPhases)
            {
                doubleswap_h =(state_vars[2]*d_state_vars[1]+
                    state_vars[3]*d_state_vars[NUM_STATE_VARS+1])*invnormhv;
                doubleswap_hvx=(state_vars[2]*d_state_vars[2]+
                    state_vars[3]*d_state_vars[NUM_STATE_VARS+2])*invnormhv;
                doubleswap_hvy=(state_vars[2]*d_state_vars[3]+
                    state_vars[3]*d_state_vars[NUM_STATE_VARS+3])*invnormhv;
            }
            if(elementType() == ElementType::SinglePhase)
            {
                doubleswap_h = (state_vars[1] * d_state_vars[0] + state_vars[2] * d_state_vars[3]) * invnormhv;
                doubleswap_hvx = (state_vars[1] * d_state_vars[1] + state_vars[2] * d_state_vars[4]) * invnormhv;
                doubleswap_hvy = (state_vars[1] * d_state_vars[2] + state_vars[2] * d_state_vars[5]) * invnormhv;
            }

            double doubleswap_h_2, doubleswap_hvx_2, doubleswap_hvy_2;
            
            double f, df, d2f, df2, d2f2, dVmag, absdVmag0, VmagOld = -1.0, dVmag2 = -1.0;
            double toler = 1.0 / ((double) (1024 * 1024 * 1024)); //toler~10^-9
            int inewt, yada = 0;
            
            for(inewt = 0; inewt < 15; inewt++)
            { //should only need about 5 newton iterations
                if(elementType() == ElementType::TwoPhases)
                {
                    doubleswap_h_2 =((state_vars[1]-prev_state_vars[1])*inv_dt+Vmag*doubleswap_h );
                    doubleswap_hvx_2=((state_vars[2]-prev_state_vars[2])*inv_dt+Vmag*doubleswap_hvx);
                    doubleswap_hvy_2=((state_vars[3]-prev_state_vars[3])*inv_dt+Vmag*doubleswap_hvy);
                }
                if(elementType() == ElementType::SinglePhase)
                {
                    doubleswap_h_2 = ((state_vars[0] - prev_state_vars[0]) * inv_dt + Vmag * doubleswap_h);
                    doubleswap_hvx_2 = ((state_vars[1] - prev_state_vars[1]) * inv_dt + Vmag * doubleswap_hvx);
                    doubleswap_hvy_2 = ((state_vars[2] - prev_state_vars[2]) * inv_dt + Vmag * doubleswap_hvy);
                }

                /* f could be greater or less than zero but f2 is always non-negative,
                 note however that minimum/optimal value of f2 could be > 0 so solving
                 for f2==0 isn't full proof, instead we want minium f2 which means solve
                 for df2==0, but df2==0 could be min or max. but in newton if you take
                 the absolute value d2f2 in the denominator you cause f2 to ALWAYS
                 decrease, and thus you are guarenteed to find a local mimimum, recall
                 we do want a minimum because f2 is guaranteed >= 0 */

                f = Vmag * Vmag * doubleswap_h_2 * doubleswap_h_2 - doubleswap_hvx_2 * doubleswap_hvx_2
                    - doubleswap_hvy_2 * doubleswap_hvy_2;

                df = 2.0 * (Vmag * doubleswap_h_2 * doubleswap_h_2 + Vmag * Vmag * doubleswap_h_2 * doubleswap_h
                        - doubleswap_hvx_2 * doubleswap_hvx - doubleswap_hvy_2 * doubleswap_hvy);
                d2f = 2.0 * (doubleswap_h_2 * doubleswap_h_2 + 4.0 * Vmag * doubleswap_h_2 * doubleswap_h
                             + Vmag * Vmag * doubleswap_h * doubleswap_h
                             - doubleswap_hvx * doubleswap_hvx - doubleswap_hvy * doubleswap_hvy);
                //f2=f*f;
                df2 = 2.0 * f * df; //solving for df2==0 is either min or max of f2
                d2f2 = 2.0 * (f * d2f + df * df);
                dVmag = (df2) ? (-df2 / fabs(d2f2)) : 0.0; //fabs(d2f2) will cause f2 to always decrease and since f2 can't be less than zero we always want it to decrease, also need to prevent division of zero by zero = nan

                if(!(shortspeed() >= 0.0))
                {
                    printf("inewt=%d Vmag0=%g Vmag=%g doubleswap_h_2=%g doubleswap_hvx_2=%g doubleswap_hvy_2=%g f=%g df=%g d2f=%g df2=%g d2f2=%g\n",
                           inewt, Vmag0, Vmag, doubleswap_h_2, doubleswap_hvx_2, doubleswap_hvy_2, f, df, d2f, df2, d2f2);
                    assert(0);
                }

                VmagOld = Vmag;
                Vmag += dVmag;

                if(Vmag < 0.0)
                    Vmag = 0.0; //safety in case it finds wrong root but this
                //shouldn't be a problem equation looks like it should be smooth in
                //the region were interested in

                dVmag2 = Vmag - VmagOld;

                if(fabs(dVmag2) < 0.0)
                    printf("VmagOld=%g dVmag2=%g yada=%d\n", VmagOld, dVmag2, yada);

                yada = yada + yada;

                if(inewt == 0)
                    absdVmag0 = fabs(dVmag);
                else if((absdVmag0 * toler >= fabs(dVmag2)) || (Vmag * toler >= fabs(dVmag2)))
                    break;
            }
            
            assert(inewt >= 0.0);
            
            set_shortspeed(Vmag);
        }
    }
    
    if(!(shortspeed() >= 0.0))
        printf("shortspeed=%g\n", shortspeed());
    
    assert(shortspeed() >= 0.0);
    
    return;
}
double* Element::eval_velocity(double xoffset, double yoffset, double Vel[])
{
    if(!(shortspeed() >= 0.0))
        printf("shortspeed=%g\n", shortspeed());

    assert(shortspeed() >= 0.0);
    
    double temp_state_vars[MAX_NUM_STATE_VARS];

    for(int ivar = 0; ivar < NUM_STATE_VARS; ivar++)
        temp_state_vars[ivar] = state_vars[ivar] + d_state_vars[ivar] * xoffset +  //distfromcenter[0]+
                                d_state_vars[NUM_STATE_VARS + ivar] * yoffset;  //distfromcenter[1];
                                
    if(elementType() == ElementType::TwoPhases)
    {
        if(!(temp_state_vars[0] > 0))
        {
            for(int i = 0; i < 4; i++)
                Vel[i] = 0;
            return Vel;
        }

        if(SHORTSPEED)
        {
            double doubleswap = (temp_state_vars[2]) * (temp_state_vars[2])
                    + (temp_state_vars[3]) * (temp_state_vars[3]);

            if(!(doubleswap > (temp_state_vars[1] * temp_state_vars[1] *
            GEOFLOW_TINY
                               * GEOFLOW_TINY)))
            {
                for(int i = 0; i < 4; i++)
                    Vel[i] = 0.;
                return Vel;
            }

            if((temp_state_vars[1] < GEOFLOW_SHORT) && (doubleswap
                    > shortspeed() * shortspeed() * temp_state_vars[1] * temp_state_vars[1]))
            {
                doubleswap = sqrt(doubleswap);
                Vel[0] = shortspeed() * temp_state_vars[2] / doubleswap;
                Vel[1] = shortspeed() * temp_state_vars[3] / doubleswap;
                Vel[2] = 0.;
                Vel[3] = 0.;
            }
            else
            {
                Vel[0] = temp_state_vars[2] / temp_state_vars[1];
                Vel[1] = temp_state_vars[3] / temp_state_vars[1];
                Vel[2] = temp_state_vars[4] / temp_state_vars[0];
                Vel[3] = temp_state_vars[5] / temp_state_vars[0];
            }
        }
        else
        {
            Vel[0] = temp_state_vars[2] / temp_state_vars[1];
            Vel[1] = temp_state_vars[3] / temp_state_vars[1];
            Vel[2] = temp_state_vars[4] / temp_state_vars[0];
            Vel[3] = temp_state_vars[5] / temp_state_vars[0];
        }
        return Vel;
    }
    if(elementType() == ElementType::SinglePhase)
    {
        if(!(temp_state_vars[0] > 0))
        {
            Vel[0] = Vel[1] = 0.0;
            return Vel;
        }

        if(SHORTSPEED)
        {
            double doubleswap = (temp_state_vars[1]) * (temp_state_vars[1])
                    + (temp_state_vars[2]) * (temp_state_vars[2]);

            if(!(doubleswap > (temp_state_vars[0] * temp_state_vars[0] *
            GEOFLOW_TINY
                               * GEOFLOW_TINY)))
            {
                Vel[0] = Vel[1] = 0.0;
                return Vel;
            }

            if((temp_state_vars[0] < GEOFLOW_SHORT) && (doubleswap
                    > shortspeed() * shortspeed() * temp_state_vars[0] * temp_state_vars[0]))
            {
                doubleswap = sqrt(doubleswap);
                Vel[0] = shortspeed() * temp_state_vars[1] / doubleswap;
                Vel[1] = shortspeed() * temp_state_vars[2] / doubleswap;
            }
            else
            {
                Vel[0] = temp_state_vars[1] / temp_state_vars[0];
                Vel[1] = temp_state_vars[2] / temp_state_vars[0];
            }
        }
        else
        {
            Vel[0] = temp_state_vars[1] / temp_state_vars[0];
            Vel[1] = temp_state_vars[2] / temp_state_vars[0];
        }
        assert((Vel[0] * Vel[0] * Vel[1] * Vel[1]) >= 0.0);

        return Vel;
    }
    return NULL;
}
void Element::calc_gravity_vector(MatProps* matprops_ptr)
{
    double max_slope = sqrt(zeta(0) * zeta(0) + zeta(1) * zeta(1));
    double max_angle = atan(max_slope);
    
    double down_slope_gravity = 9.8 * sin(max_angle);
    if(dabs(down_slope_gravity) > GEOFLOW_TINY)
    {
        gravity(0, -down_slope_gravity * zeta(0) / max_slope);
        gravity(1, -down_slope_gravity * zeta(1) / max_slope);
        //   gravity[0] = -down_slope_gravity*cos(atan(1.0));
        //   gravity[1] = -down_slope_gravity*sin(atan(1.0));
        
        gravity(2, 9.8 * cos(max_angle));
    }
    else
    {
        gravity(0, 0.0);
        gravity(1, 0.0);
        gravity(2, 9.8);
    }
    
    for(int i = 0; i < 3; i++)
        gravity(i, gravity(i) / matprops_ptr->GRAVITY_SCALE);
    
    return;
}

int Element::determine_refinement(double target)
{
    int flag = 0, i;
    
    /*  for(i=0;i<NUM_STATE_VARS*DIMENSION;i++)
     if(dabs(d_state_vars[i]) > target)
     flag = 1;
     */
    if(state_vars[0] > target)
        flag = 1;
    return flag;
}

void Element::calc_d_gravity(HashTable* El_Table)
{
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side();
    switch (positive_x_side())
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
    }
    /* x direction */
    Element* ep = (Element*) (El_Table->lookup(neighbor(xp)));
    Element* em = (Element*) (El_Table->lookup(neighbor(xm)));
    int j;
    if(ep != NULL && em != NULL)
    {
        double dp, dm, dxp, dxm;
        dxp = ep->coord(0) - coord(0);
        dp = (ep->gravity(2) - gravity(2)) / dxp;
        dxm = coord(0) - em->coord(0);
        dm = (gravity(2) - em->gravity(2)) / dxm;
        
        d_gravity(0, (dp * dxm + dm * dxp) / (dxm + dxp));  // weighted average 
    }
    else if(em != NULL)
    {
        double dm, dxm;
        dxm = coord(0) - em->coord(0);
        d_gravity(0, (gravity(2) - em->gravity(2)) / dxm);
    }
    else if(ep != NULL)
    {
        double dp, dxp;
        dxp = ep->coord(0) - coord(0);
        d_gravity(0, (ep->gravity(2) - gravity(2)) / dxp);
    }
    else
        //no neighbors on either side -- assume that the ground is flat
        d_gravity(0, 0.0);
    
    /* y direction */
    ep = (Element*) (El_Table->lookup(neighbor(yp)));
    em = (Element*) (El_Table->lookup(neighbor(ym)));
    if(ep != NULL && em != NULL)
    {
        double dp, dm, dxp, dxm;
        dxp = ep->coord(1) - coord(1);
        dp = (ep->gravity(2) - gravity(2)) / dxp;
        dxm = coord(1) - em->coord(1);
        dm = (gravity(2) - em->gravity(2)) / dxm;
        
        d_gravity(1, (dp * dxm + dm * dxp) / (dxm + dxp));  // weighted average 
    }
    else if(em != NULL)
    {
        double dm, dxm;
        dxm = coord(1) - em->coord(1);
        d_gravity(1, (gravity(2) - em->gravity(2)) / dxm);
    }
    else if(ep != NULL)
    {
        double dp, dxp;
        dxp = ep->coord(1) - coord(1);
        d_gravity(1, (ep->gravity(2) - gravity(2)) / dxp);
    }
    else
        //no neighbors on either side -- assume that the ground is flat
        d_gravity(1, 0.0);
    
    return;
}

void Element::calc_topo_data(MatProps* matprops_ptr)
{
    double resolution = (dx(0)/*/(zeta[0]*zeta[0]+1)*/+ dx(1)/*/(zeta[1]*zeta[1]+1)*/) * (matprops_ptr->LENGTH_SCALE)
            / 2.0;  // element "size"
    double xcoord = coord(0) * (matprops_ptr->LENGTH_SCALE);
    double ycoord = coord(1) * (matprops_ptr->LENGTH_SCALE);
    //double eldif = elevation;
    int i = Get_elevation(resolution, xcoord, ycoord, elevation_ref());
#ifdef PRINT_GIS_ERRORS
    if(i != 0)
    {   
        printf("Error in Get_elevation(error code %d)\n", i);
        exit(1);
    }
#endif
    set_elevation(elevation() / matprops_ptr->LENGTH_SCALE);
    //eldif=(elevation-eldif)*matprops_ptr->LENGTH_SCALE;
    //if(fabs(eldif)>1.0) printf("calc_topo_data() after-before=%g\n",eldif);
    i = Get_slope(resolution, xcoord, ycoord, zeta_ref(0), zeta_ref(1));
#ifdef PRINT_GIS_ERRORS
    if(i != 0)
    {   
        printf("Error in Get_slope(error code %d)\n", i);
        exit(1);
    }
#endif
    i = Get_curvature(resolution, xcoord, ycoord, curvature_ref(0), curvature_ref(1));
#ifdef PRINT_GIS_ERRORS
    if(i != 0)
    {   
        printf("Error in Get_curvature(error code %d)\n", i);
        exit(1);
    }
#endif
    curvature(0, curvature(0) * (matprops_ptr->LENGTH_SCALE));
    curvature(1, curvature(1) * (matprops_ptr->LENGTH_SCALE));
    
    if(matprops_ptr->material_count == 1)  //only one material so don't need map  
        set_material(1);  //GIS material id tag/index starts from 1
    else //more than one material so need to get material from map
    {
        int m;
        Get_raster_id(resolution, xcoord, ycoord, &m);
        set_material(m);
        if((material() < 1) || (material() > 1000))
        {
            fprintf(stderr, "ERROR garbage material %d read from material map\n", material());
            exit(1);
        }
    }
    
    return;
}

void Element::calc_flux_balance(HashTable* NodeTable)
{
    int i, j;
    double flux[3] =
    { 0, 0, 0 };
    int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
    xp = positive_x_side();
    switch (positive_x_side())
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
    }
    Node *nd_xp, *nd_xn, *nd_yp, *nd_yn;
    nd_xp = (Node*) NodeTable->lookup(node_key(xp + 4));
    nd_xn = (Node*) NodeTable->lookup(node_key(xm + 4));
    nd_yp = (Node*) NodeTable->lookup(node_key(yp + 4));
    nd_yn = (Node*) NodeTable->lookup(node_key(ym + 4));
    for(j = 0; j < 3; j++)
        flux[j] = dabs(nd_xp->refinementflux[j] - nd_xn->refinementflux[j])
                + dabs(nd_yp->refinementflux[j] - nd_yn->refinementflux[j]);
    
    set_el_error(0, 0.0);
    for(j = 0; j < 3; j++)
        set_el_error(0, el_error(0)+flux[j]);
    
    set_el_error(0, 2. * el_error(0) * el_error(0) / (dx(0) + dx(1)) + WEIGHT_ADJUSTER); //W_A is so that elements with pile height = 0 have some weight.
            
    return;
}

// load-balancing stuff ...

void Element::put_coord(double* coord_in)
{
    int i;
    for(i = 0; i < KEYLENGTH; i++)
        set_coord(i, coord_in[i]);
    return;
}

/*
 using elm_loc, which_son is calculated
 */
void Element::calc_which_son()
{
    if(elm_loc(0) % 2 == 0)
    {
        if(elm_loc(1) % 2 == 0)
            set_which_son(0);
        else
            set_which_son(3);
    }
    else
    {
        if(elm_loc(1) % 2 == 0)
            set_which_son(1);
        else
            set_which_son(2);
    }
    
}

//should be full proof way to get key of opposite brother, 
//as long as you know your own coord, dx, and which_son
void Element::find_opposite_brother(HashTable* El_Table)
{
    
    if(opposite_brother_flag() == 1)
        return;
    
    set_brother((which_son() + 2) % 4, sfc_key_zero);

    SFC_Key nullkey = 0;

    if(!((brother((which_son() + 1) % 4)==sfc_key_null) && (brother((which_son() + 3) % 4)==sfc_key_null)))
    {
        //use space filling curve to compute the key of opposite
        //brother from it's bubble node coordinates
        double bro_norm_coord[2];
        unsigned oldkey[KEYLENGTH];
        unsigned nkey = KEYLENGTH;
        
        if((which_son() == 0) || (which_son() == 3))
            bro_norm_coord[0] = El_Table->get_invdxrange() * (coord(0) + dx(0) - *(El_Table->get_Xrange() + 0));
        else
            bro_norm_coord[0] = El_Table->get_invdxrange() * (coord(0) - dx(0) - *(El_Table->get_Xrange() + 0));
        
        if((which_son() == 0) || (which_son() == 1))
            bro_norm_coord[1] = El_Table->get_invdyrange() * (coord(1) + dx(1) - *(El_Table->get_Yrange() + 0));
        else
            bro_norm_coord[1] = El_Table->get_invdyrange() * (coord(1) - dx(1) - *(El_Table->get_Yrange() + 0));
        
        fhsfc2d_(bro_norm_coord, &nkey, oldkey);
        set_brother((which_son() + 2) % 4,sfc_key_from_oldkey(oldkey));
        
        set_opposite_brother_flag(1);
    }
    
    return;
    
}
#ifdef OLDCODE
void Element::find_opposite_brother(HashTable* El_Table)
{   
    /* brother information -- requires that atleast one of this
     element's neighboring brothers is on this process in 
     order to get information onthe brother that is not a neighbor */
    Element* EmTemp;
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int i;
    if(opposite_brother_flag == 1)
    return;
    switch(which_son)
    {   
        case 0:
        if(neigh_proc[2] != -1 && neigh_proc[1] != -1)
        {   
            if(neigh_proc[2] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(2));
                if(*(EmTemp->get_neigh_gen()+1) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[2][i] = *(EmTemp->get_neighbors()+KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[2][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }
            }
            else if(neigh_proc[1] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(5));
                if(*(EmTemp->get_neigh_gen()+2) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[2][i] = *(EmTemp->get_neighbors()+2*KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+6*KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[2][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }

            }
        }
        break;
        case 1:
        if(neigh_proc[2] != -1 && neigh_proc[3] != -1)
        {   
            if(neigh_proc[2] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(6));
                if(*(EmTemp->get_neigh_gen()+3) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[3][i] = *(EmTemp->get_neighbors()+3*KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+7*KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[3][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }
            }
            else if(neigh_proc[3] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(3));
                if(*(EmTemp->get_neigh_gen()+2) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[3][i] = *(EmTemp->get_neighbors()+2*KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+2*KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[3][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }

            }
        }
        break;
        case 2:
        if(neigh_proc[0] != -1 && neigh_proc[3] != -1)
        {   
            if(neigh_proc[0] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(0));
                if(*(EmTemp->get_neigh_gen()+3) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[0][i] = *(EmTemp->get_neighbors()+3*KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+3*KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[0][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }
            }
            else if(neigh_proc[3] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(3));
                if(*(EmTemp->get_neigh_gen()) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[0][i] = *(EmTemp->get_neighbors()+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors());
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[0][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }

            }
        }
        break;
        case 3:
        if(neigh_proc[4] != -1 && neigh_proc[1] != -1)
        {   
            if(neigh_proc[4] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(4));
                if(*(EmTemp->get_neigh_gen()+1) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[1][i] = *(EmTemp->get_neighbors()+KEYLENGTH+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors()+KEYLENGTH);
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[1][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }
            }
            else if(neigh_proc[1] == myid)
            {   
                EmTemp = (Element*) El_Table->lookup(neighbor(1));
                if(*(EmTemp->get_neigh_gen()) == generation)
                {   
                    for(i=0;i<KEYLENGTH;i++)
                    brothers[1][i] = *(EmTemp->get_neighbors()+i);
                    opposite_brother_flag = 1;
                }
                else
                {   
                    EmTemp = (Element*) El_Table->lookup(EmTemp->get_neighbors());
                    if(EmTemp != NULL && EmTemp->get_refined_flag() != GHOST && EmTemp->get_gen() == (generation+1))
                    {   
                        unsigned* elm_father = EmTemp->getfather();
                        for(i=0;i<KEYLENGTH;i++)
                        brothers[1][i] = elm_father[i];
                        opposite_brother_flag = 1;
                    }
                }
            }
        }
        break;
    }

    return;
}
#endif

/* the _only_ place calc_stop_crit() is called is in get_coef_and_eigen.C 
 immediately after k_active/passive is calculated AND in init_piles.C when
 computing the initial "deposited" volume.

 */
void Element::calc_stop_crit(MatProps *matprops_ptr)
{
    
    double stopcrit;
    
    effect_kactxy(0, kactxy(0));
    effect_kactxy(1, kactxy(1));
    set_effect_bedfrict(matprops_ptr->bedfrict[material()]);
    set_effect_tanbedfrict(matprops_ptr->tanbedfrict[material()]);
    
#ifdef STOPCRIT_CHANGE_BED
    if(stoppedflags()==2)
    {   
        effect_kactxy(0,matprops_ptr->epsilon);
        effect_kactxy(1,matprops_ptr->epsilon);
        set_effect_bedfrict(matprops_ptr->intfrict);
        set_effect_tanbedfrict(matprops_ptr->tanintfrict);
    }
#endif
    
    set_stoppedflags(0);
    if(elementType() == ElementType::TwoPhases)
    {
        return;
    }
    if(state_vars[0] <= GEOFLOW_TINY)
        stopcrit = HUGE_VAL;
    else
    {
        double dirx, diry, Vtemp, bedslope, slopetemp;
        double VxVy[2];
        eval_velocity(0.0, 0.0, VxVy);
        
        Vtemp = sqrt(VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1]);
        
        //these are the element force balance test, will the pile slide
        if(Vtemp > 0)
        {
            dirx = VxVy[0] / Vtemp;
            diry = VxVy[1] / Vtemp;
            
            bedslope = dirx * zeta(0) + diry * zeta(1);
            slopetemp = bedslope
                    + (dirx * effect_kactxy(0) * d_state_vars[0] + diry * effect_kactxy(1)
                                                                   * d_state_vars[NUM_STATE_VARS]);
            
            stopcrit = (slopetemp + effect_tanbedfrict()) / Vtemp * NUM_FREEFALLS_2_STOP
                       * sqrt(2.0 * 9.8 / matprops_ptr->GRAVITY_SCALE * state_vars[0] / (1 + bedslope * bedslope));
        }
        else
        {
            set_stoppedflags(1); //erosion off
            bedslope = -sqrt(zeta(0) * zeta(0) + zeta(1) * zeta(1));
            if(bedslope < 0)
            {
                slopetemp = bedslope
                        + (zeta(0) / bedslope * effect_kactxy(0) * d_state_vars[0] + zeta(1) / bedslope
                                * effect_kactxy(1) * d_state_vars[NUM_STATE_VARS]);
                
                stopcrit = sign(slopetemp + effect_tanbedfrict()) * HUGE_VAL;
            }
            else
                stopcrit = HUGE_VAL;
        }
        
        //this is the internal friction test, will the pile slump
        if(stopcrit >= 1.0)
        {
            set_stoppedflags(1);
            slopetemp = -sqrt(
                    (zeta(0) + effect_kactxy(0) * d_state_vars[0]) * (zeta(0) + effect_kactxy(0) * d_state_vars[0])
                    * (zeta(1) + effect_kactxy(1) * d_state_vars[NUM_STATE_VARS])
                    * (zeta(1) + effect_kactxy(1) * d_state_vars[NUM_STATE_VARS]));
            if(slopetemp + tan(matprops_ptr->intfrict) > 0)
                set_stoppedflags(2);
        }
    }
    
    //stoppedflags=2;
    
    /* 
     do you choose to set bed friction to internal friction in stopped cells?
     because friction at bed can only support up to the bed friction angle for
     an angle of repose 
     */
#ifdef STOPCRIT_CHANGE_BED
    if(stoppedflags()==2)
    {   
        effect_kactxy(0,matprops_ptr->epsilon);
        effect_kactxy(1,matprops_ptr->epsilon);
        set_effect_bedfrict(matprops_ptr->intfrict);
        set_effect_tanbedfrict(matprops_ptr->tanintfrict);
    }
    else
    {   
        effect_kactxy(0,kactxy(0));
        effect_kactxy(1,kactxy(1));
        set_effect_bedfrict(matprops_ptr->bedfrict[material()]);
        set_effect_tanbedfrict(matprops_ptr->tanbedfrict[material()]);
    }
#endif
    
    return;
    
}

int Element::if_pile_boundary(HashTable *ElemTable, double contour_height)
{
    
    int ineigh;
    Element* ElemNeigh;
    
    assert(state_vars[0] >= 0.0);
    
    if(state_vars[0] >= contour_height)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(ElemNeigh == NULL)
                {
                    cout<<"ElemNeigh==NULL ineigh="<<ineigh<<"\n";
                    cout<<" mykey   ={"<<key()<<"} myprocess ="<<myprocess()<<" generation="<<generation()<<" refined="<<refined_flag()<<" adapted="<<adapted_flag()<<"\n";
                    cout<<" neighbor={"<<neighbor(ineigh)<<"} neigh_proc="<<neigh_proc(ineigh)<<" neigh_gen ="<<neigh_gen(ineigh) <<"\n\n";
                    cout.flush();
                }
                assert(ElemNeigh);
                if(*(ElemNeigh->get_state_vars() + 0) < contour_height)
                    return (2); //inside of pileheight contour line
            }
        //else if(neigh_proc[ineigh%4]==-1) return(2); //pileheight on boundary of domain
    }
    else
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(ElemNeigh == NULL)
                {
                    cout<<"ElemNeigh==NULL ineigh="<<ineigh<<"\n";
                    cout<<" mykey   ={"<<key()<<"} myprocess ="<<myprocess()<<" generation="<<generation()<<" refined="<<refined_flag()<<" adapted="<<adapted_flag()<<"\n";
                    cout<<" neighbor={"<<neighbor(ineigh)<<"} neigh_proc="<<neigh_proc(ineigh)<<" neigh_gen ="<<neigh_gen(ineigh) <<"\n\n";
                    cout.flush();
                }
                assert(ElemNeigh);
                assert(*(ElemNeigh->get_state_vars() + 0) >= 0.0);
                if(*(ElemNeigh->get_state_vars() + 0) >= contour_height)
                    return (1); //outside of pileheight contour line
            }
    }
    
    return (0); //not on pileheight contour line
}

int Element::if_source_boundary(HashTable *ElemTable)
{
    
    int ineigh;
    Element* ElemNeigh;
    
    if(!(Influx(0) >= 0.0))
    {
        printf("if_source_boundary() Influx[0]=%g\n", Influx(0));
        fflush(stdout);
    }
    
    assert(Influx(0) >= 0.0); //currently mass sinks are not allowed
    
    if(Influx(0) > 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(ElemNeigh == NULL)
                {
                    cout<<"ElemNeigh==NULL ineigh="<<ineigh<<"\n";
                    cout<<" mykey   ={"<<key()<<"} myprocess ="<<myprocess()<<" generation="<<generation()<<" refined="<<refined_flag()<<" adapted="<<adapted_flag()<<"\n";
                    cout<<" neighbor={"<<neighbor(ineigh)<<"} neigh_proc="<<neigh_proc(ineigh)<<" neigh_gen ="<<neigh_gen(ineigh) <<"\n\n";
                    cout.flush();
                }
                assert(ElemNeigh);
                if(ElemNeigh->Influx(0) <= 0.0)
                    return (2); //inside of line bounding area with a mass source 
            }
        //else if(neigh_proc[ineigh%4]==-1) return(2); //mass source on boundary of domain
    }
    
    else if(Influx(0) == 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0.0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(ElemNeigh == NULL)
                {
                    cout<<"ElemNeigh==NULL ineigh="<<ineigh<<"\n";
                    cout<<" mykey   ={"<<key()<<"} myprocess ="<<myprocess()<<" generation="<<generation()<<" refined="<<refined_flag()<<" adapted="<<adapted_flag()<<"\n";
                    cout<<" neighbor={"<<neighbor(ineigh)<<"} neigh_proc="<<neigh_proc(ineigh)<<" neigh_gen ="<<neigh_gen(ineigh) <<"\n\n";
                    cout.flush();
                }
                assert(ElemNeigh);
                assert(ElemNeigh->Influx(0) >= 0.0);
                if(ElemNeigh->Influx(0) != 0.0)
                    return (1); //outside of line bounding area with a mass source/sink 
            }
    }
    else if(Influx(0) < 0.0)
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0.0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(ElemNeigh == NULL)
                {
                    cout<<"ElemNeigh==NULL ineigh="<<ineigh<<"\n";
                    cout<<" mykey   ={"<<key()<<"} myprocess ="<<myprocess()<<" generation="<<generation()<<" refined="<<refined_flag()<<" adapted="<<adapted_flag()<<"\n";
                    cout<<" neighbor={"<<neighbor(ineigh)<<"} neigh_proc="<<neigh_proc(ineigh)<<" neigh_gen ="<<neigh_gen(ineigh) <<"\n\n";
                    cout.flush();
                }
                assert(ElemNeigh);
                if(ElemNeigh->Influx(0) >= 0.0)
                    return (-1); //inside of line bounding area with a mass sink 
            }
        //else if(neigh_proc[ineigh%4]==-1) return(-1); //mass sink on boundary of domain
    }
    
    return (0); //not on line bounding area with mass source/sink
}

int Element::if_first_buffer_boundary(HashTable *ElemTable, double contour_height)
{
    
    int ineigh;
    Element* ElemNeigh;
    int iffirstbuffer = 0;
    
    if(adapted_flag() <= 0)
        return (adapted_flag() - 1);
    
    assert(state_vars[0] >= 0.0);
    assert(Influx(0) >= 0.0);
    if((state_vars[0] < contour_height) && (Influx(0) == 0.0))
    {
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                /*
                 if(ElemNeigh==NULL){
                 printf("ElemNeigh==NULL\n mykey   ={%u,%u} myprocess =%d generation=%d refined=%d adapted=%d\n",
                 key[0],key[1],myprocess,generation,refined,adapted);
                 printf(" neighbor={%u,%u} neigh_proc=%d neigh_gen =%d\n\n",
                 neighbor(ineigh)[0],neighbor(ineigh)[1],neigh_proc[ineigh],neigh_gen[ineigh]);
                 fflush(stdout);
                 }
                 */
                assert(ElemNeigh);
                /*
                 assert((ElemNeigh->get_adapted_flag()!=0)&&
                 ((abs(ElemNeigh->get_adapted_flag())<=BUFFER)||
                 (ElemNeigh->get_adapted_flag()==OLDFATHER)
                 )
                 );
                 
                 if(!(*(ElemNeigh->get_influx())>=0.0)){
                 printf("Neigh Influx<0\n");
                 printf("Elem     ={%u,%u} myprocess=%d generation=%d refined=%d adapted=%d Influx[0]=%g\n",
                 key[0],key[1],myprocess,generation,refined,adapted,Influx[0]);
                 printf("Elemneigh={%u,%u} myprocess=%d generation=%d refined=%d adapted=%d Influx[0]=%g\n",
                 neighbor(ineigh)[0],neighbor(ineigh)[1],ElemNeigh->get_myprocess(),ElemNeigh->get_gen(),
                 ElemNeigh->get_refined_flag(),ElemNeigh->get_adapted_flag(),ElemNeigh->get_influx());
                 fflush(stdout);
                 }
                 assert(*(ElemNeigh->get_influx())>=0.0);
                 assert(*(ElemNeigh->get_state_vars())>=0.0);
                 */
                if((*(ElemNeigh->get_state_vars() + 0) >= contour_height) || (ElemNeigh->Influx(0) > 0.0))
                {
                    iffirstbuffer = 1;
                    break;
                }
            }
    }
    else
    {/* if((state_vars[0]>=contour_height)||
     (Influx[0]>0)){*/
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                assert(ElemNeigh);
                if((*(ElemNeigh->get_state_vars() + 0) < contour_height) && (ElemNeigh->Influx(0) == 0.0))
                {
                    iffirstbuffer = 1;
                    break;
                }
            }
    }
    
    if(iffirstbuffer)
    {
        if((adapted_flag() >= NEWSON) || (generation() == REFINE_LEVEL))
            return (2); //is a member of the buffer but doesn't need to be refined
        else
            return (1); //needs to be refined and some of its sons will be members
    }
    
    return (0);
}

int Element::if_next_buffer_boundary(HashTable *ElemTable, HashTable *NodeTable, double contour_height)
{
    
    int ineigh;
    Element* ElemNeigh;
    int ifnextbuffer;
    ifnextbuffer = 0;
    if(adapted_flag() <= 0)
        //GHOST element or element that should be deleted soon
        return (adapted_flag() - 1);
    
    if((adapted_flag() != BUFFER) && //this element is not in the buffer
    ((Influx(0) == 0.0)/*&&(state_vars[0]<contour_height)*/) //&& //this element is OUTSIDE the buffer layer "circle"
    //(state_vars[0]>=GEOFLOW_TINY)
    )
        for(ineigh = 0; ineigh < 8; ineigh++)
            if(neigh_proc(ineigh) >= 0)
            { //don't check outside map boundary or duplicate neighbor
                ElemNeigh = (Element*) ElemTable->lookup(neighbor(ineigh));
                if(!ElemNeigh)
                {
                    cout<<"Elem={"<<key()<<"} missing neighbor ineigh="<<ineigh<<" {"<<neighbor(ineigh)<<"}\n";
                    ElemBackgroundCheck(ElemTable, NodeTable, key(), stdout);
                    assert(ElemNeigh);
                }
                //assert(ElemNeigh->get_adapted_flag()!=0);
                
                if((abs(ElemNeigh->adapted_flag()) == BUFFER) && (state_vars[0]
                        <= *(ElemNeigh->get_state_vars() + 0)))
                //if(abs(ElemNeigh->get_adapted_flag())==BUFFER)
                { //this element is next to a member of the old buffer layer
                  //if((ElemNeigh->get_adapted_flag())==BUFFER){ //this element is next to a member of the old buffer layer
                    ifnextbuffer = 1; //which means this element is a member of the next outer boundary of the buffer layer
                    break;
                }
            }
    
    if(ifnextbuffer == 1)
    {
        if((adapted_flag() >= NEWSON) || (generation() == REFINE_LEVEL))
            return (2); //is a member of the buffer but doesn't need to be refined
        else
            return (1); //needs to be refined and some of its sons will be members
    }
    
    return (0);
}

void Element::save_elem(FILE* fp, FILE *fptxt)
{
    
    FourBytes temp4;
    EightBytes temp8;
    unsigned writespace[138];
    
    int Itemp = 0, itemp, jtemp;
    for(itemp = 0; itemp < 2; itemp++)
    {
        temp4.i = elm_loc(itemp);
        writespace[Itemp++] = temp4.u;
    }
    assert(Itemp == 2);
    
#ifdef DEBUG_SAVE_ELEM
    //FILE *fpdb=fopen("save_elem.debug","w");
    FILE *fpdb=fptxt;
    fprintf(fpdb,"\n\nelm_loc=%d %d\n",elm_loc(0),elm_loc(1));
#endif
    
    temp4.i = generation();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 3);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"generation=%d\n",generation());
#endif
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"order={ ");
#endif
    for(itemp = 0; itemp < 5; itemp++)
    {
        temp4.i = order(itemp);
        writespace[Itemp++] = temp4.u;
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"%d ",order(itemp));
#endif
    }
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\n");
#endif
    assert(Itemp == 8);
    
    temp4.i = opposite_brother_flag();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 9);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"opposite_brother_flag=%d\n",opposite_brother_flag());
#endif
    
    temp4.i = new_old();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 10);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"new_old=%d\n",new_old());
#endif
    
    temp4.i = material();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 11);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"material=%d\n",material());
#endif
    
    temp4.i = ndof();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 12);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"ndof=%d\n",ndof());
#endif
    
    temp4.i = refined_flag();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 13);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"refined=%d\n",refined_flag());
#endif
    
    temp4.i = adapted_flag();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 14);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"adapted=%d\n",adapted_flag());
#endif
    
    temp8.d = lb_weight();
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    assert(Itemp == 16);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"lb_weight=%g\n",lb_weight());
#endif
    sfc_key_write_to_space(lb_key(),writespace,Itemp);
    assert(Itemp == 18);
/*#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"lb_key=%u %u\n",lb_key[0],lb_key[1]);
#endif*/
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"neigh_proc={ ");
#endif
    for(itemp = 0; itemp < 8; itemp++)
    {
        temp4.i = neigh_proc(itemp);
        writespace[Itemp++] = temp4.u;
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"%d ",neigh_proc(itemp));
#endif
    }
    assert(Itemp == 26);
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nneigh_gen={ ");
#endif
    for(itemp = 0; itemp < 8; itemp++)
    {
        temp4.i = neigh_gen(itemp);
        writespace[Itemp++] = temp4.u;
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"%d ",neigh_gen(itemp));
#endif 
    }
    assert(Itemp == 34);
    
    sfc_key_write_to_space(key(),writespace,Itemp);
    assert(Itemp == 36);
/*#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nkey=%u %u\nnode_key={ ",key[0],key[1]);
#endif  */
    
    for(itemp = 0; itemp < 8; itemp++)
    {
        sfc_key_write_to_space(node_key(itemp),writespace,Itemp);
/*#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"(%u %u) ",node_key[itemp][0],node_key[itemp][1]);
#endif */
    }
    assert(Itemp == 52);
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nneighbor={ ");
#endif 
    for(itemp = 0; itemp < 8; itemp++)
    {
        sfc_key_write_to_space(neighbor(itemp),writespace,Itemp);
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"(%u %u) ",neighbor(itemp)[0],neighbor(itemp)[1]);
#endif 
    }
    assert(Itemp == 68);
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nbrothers={ ");
#endif 
    for(itemp = 0; itemp < 4; itemp++)
    {
        sfc_key_write_to_space(brother(itemp),writespace,Itemp);
#ifdef DEBUG_SAVE_ELEM
        //fprintf(fpdb,"(%u %u) ",brother(itemp)[0],brother(itemp)[1]);
#endif 
    }
    assert(Itemp == 76);
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nson={ ");
#endif 
    for(itemp = 0; itemp < 4; itemp++)
    {
        sfc_key_write_to_space(son(itemp),writespace,Itemp);
#ifdef DEBUG_SAVE_ELEM
        //fprintf(fpdb,"(%u %u) ",son(itemp)[0],son(itemp)[1]);
#endif 
    }
    assert(Itemp == 84);

    //int ItempComp=84;
    
    for(itemp = 0; itemp < NUM_STATE_VARS; itemp++)
    {
        temp8.d = state_vars[itemp];
        writespace[Itemp++] = temp8.u[0];
        writespace[Itemp++] = temp8.u[1];
    }

    //ItempComp+=2*NUM_STATE_VARS;
    //assert(Itemp == 90);
    
    //don't need prev_state_vars or d_state_vars do need shortspeed
    temp8.d = shortspeed();
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    //assert(Itemp == 92);
    
#ifdef DEBUG_SAVE_ELEM
    fprintf(fpdb,"}\nstate_vars={ %g %g %g }\nshortspeed=%g\n",
            state_vars[0],state_vars[1],state_vars[2],shortspeed());
#endif 
    
    //bob  
    
    temp8.i[0] = iwetnode();
    writespace[Itemp++] = temp8.u[0];
    //assert(Itemp == 93);
    
    temp8.d = Awet();
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    //assert(Itemp == 95);
    
    temp8.d = Swet();
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    //assert(Itemp == 97);
    
    temp8.d = drypoint(0);
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    //assert(Itemp == 99);
    
    temp8.d = drypoint(1);
    writespace[Itemp++] = temp8.u[0];
    writespace[Itemp++] = temp8.u[1];
    //assert(Itemp == 101);
    
    //boundary conditions start here
    if(bcptr() == nullptr)
    {
        writespace[Itemp++] = 0;
        //assert(Itemp==118);
        //assert(Itemp == 102);
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"num_extra=0\n");
#endif
    }
    else
    {
        writespace[Itemp++] = 20;
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"num_extra=20\nbcptr->type={ ");
#endif
        for(itemp = 0; itemp < 4; itemp++)
        {
            temp4.i = bcptr()->type[itemp];
            writespace[Itemp++] = temp4.u;
#ifdef DEBUG_SAVE_ELEM
            fprintf(fpdb,"%d ",bcptr()->type[itemp]);
#endif
        }
        //assert(Itemp==122);
        //assert(Itemp == 106);
#ifdef DEBUG_SAVE_ELEM
        fprintf(fpdb,"}\n");
#endif
        for(itemp = 0; itemp < 4; itemp++)
        {
            temp4.f = bcptr()->value[itemp][0][0];
            writespace[Itemp++] = temp4.u;
            temp4.f = bcptr()->value[itemp][0][1];
            writespace[Itemp++] = temp4.u;
            temp4.f = bcptr()->value[itemp][1][0];
            writespace[Itemp++] = temp4.u;
            temp4.f = bcptr()->value[itemp][1][1];
            writespace[Itemp++] = temp4.u;
#ifdef DEBUG_SAVE_ELEM
            fprintf(fpdb,"bcptr->value={ %f %f %f %f }\n",
                    bcptr()->value[itemp][0][0],bcptr()->value[itemp][0][1],
                    bcptr()->value[itemp][1][0],bcptr()->value[itemp][1][1]);
#endif      
        }
        //assert(Itemp==138);
        //assert(Itemp == 122);
    }
#ifdef DEBUG_SAVE_ELEM
    //fclose(fpdb);
#endif
    
    fwrite(writespace, sizeof(unsigned), Itemp, fp);
    
    return;
}

Element::Element(FILE* fp, HashTable* NodeTable, MatProps* matprops_ptr, int myid)
{
    set_father(sfc_key_zero);
    for(int i = 0; i < 4; i++){
        set_brother(i, sfc_key_zero);
        set_son(i, sfc_key_zero);
    }
    set_lb_key(sfc_key_zero);
    
    for(int i = 0; i < NUM_STATE_VARS; i++)
        Influx(i, 0.0);
    set_myprocess(myid);
    set_no_of_eqns(EQUATIONS);
    //refined=0;
    
    FourBytes temp4;
    EightBytes temp8;
    unsigned readspace[102];
    
    fread(readspace, sizeof(unsigned), 102, fp);
    
    //read the element here
    int Itemp = 0, itemp, jtemp;
    for(itemp = 0; itemp < 2; itemp++)
    {
        temp4.u = readspace[Itemp++];
        set_elm_loc(itemp, temp4.i);
    }
    assert(Itemp == 2);
    
    temp4.u = readspace[Itemp++];
    generation(temp4.i);
    assert(Itemp == 3);
    
    for(itemp = 0; itemp < 5; itemp++)
    {
        temp4.u = readspace[Itemp++];
        set_order(itemp, temp4.i);
    }
    assert(Itemp == 8);
    
    temp4.u = readspace[Itemp++];
    set_opposite_brother_flag(temp4.i);
    assert(Itemp == 9);
    
    temp4.u = readspace[Itemp++];
    set_new_old(temp4.i);
    assert(Itemp == 10);
    
    temp4.u = readspace[Itemp++];
    set_material(temp4.i);
    assert(Itemp == 11);
    
    temp4.u = readspace[Itemp++];
    set_ndof(temp4.i);
    assert(Itemp == 12);
    
    temp4.u = readspace[Itemp++];
    set_refined_flag(temp4.i);
    assert(Itemp == 13);
    
    temp4.u = readspace[Itemp++];
    set_adapted_flag(temp4.i);
    assert(Itemp == 14);
    
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    set_lb_weight(temp8.d);
    assert(Itemp == 16);
    
    set_lb_key(sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 18);
    
    for(itemp = 0; itemp < 8; itemp++)
    {
        temp4.u = readspace[Itemp++];
        set_neigh_proc(itemp, temp4.i);
    }
    assert(Itemp == 26);
    
    for(itemp = 0; itemp < 8; itemp++)
    {
        temp4.u = readspace[Itemp++];
        get_neigh_gen(itemp, temp4.i);
    }
    assert(Itemp == 34);
    
    set_key(sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 36);
    
    for(itemp = 0; itemp < 8; itemp++)
        set_node_key(itemp,sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 52);
    
    for(itemp = 0; itemp < 8; itemp++)
        set_neighbor(itemp,sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 68);
    
    for(itemp = 0; itemp < 4; itemp++)
        set_brother(itemp,sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 76);
    
    for(itemp = 0; itemp < 4; itemp++)
        set_son(itemp,sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 84);
    
    for(itemp = 0; itemp < NUM_STATE_VARS; itemp++)
    {
        temp8.u[0] = readspace[Itemp++];
        temp8.u[1] = readspace[Itemp++];
        state_vars[itemp] = temp8.d;
    }
    assert(Itemp == 90);
    
    //don't need prev_state_vars or d_state_vars do need shortspeed
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    set_shortspeed(temp8.d);
    assert(Itemp == 92);
    
    for(int i = 0; i < NUM_STATE_VARS; i++)
    {
        prev_state_vars[i] = 0.;
        d_state_vars[i] = 0.;
        d_state_vars[NUM_STATE_VARS + i] = 0;
    }
    
    temp8.u[0] = readspace[Itemp++];
    set_iwetnode(temp8.i[0]);
    assert(Itemp == 93);
    
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    set_Awet(temp8.d);
    assert(Itemp == 95);
    
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    set_Swet(temp8.d);
    assert(Itemp == 97);
    
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    drypoint(0, temp8.d);
    assert(Itemp == 99);
    
    temp8.u[0] = readspace[Itemp++];
    temp8.u[1] = readspace[Itemp++];
    drypoint(1, temp8.d);
    assert(Itemp == 101);
    
    if(readspace[Itemp] > 0)
    {
        int num_extra = readspace[Itemp];
        fread(readspace, sizeof(unsigned), num_extra, fp);
        Itemp = 0;
        
        bcptr(new BC);
        
        //boundary conditions start here
        for(itemp = 0; itemp < 4; itemp++)
        {
            temp4.u = readspace[Itemp++];
            bcptr()->type[itemp] = temp4.i;
        }
        assert(Itemp == 4);
        
        for(itemp = 0; itemp < 4; itemp++)
        {
            temp4.u = readspace[Itemp++];
            bcptr()->value[itemp][0][0] = temp4.f;
            
            temp4.u = readspace[Itemp++];
            bcptr()->value[itemp][0][1] = temp4.f;
            
            temp4.u = readspace[Itemp++];
            bcptr()->value[itemp][1][0] = temp4.f;
            
            temp4.u = readspace[Itemp++];
            bcptr()->value[itemp][1][1] = temp4.f;
        }
        assert(Itemp == 20);
    }
    else
        void_bcptr();
    
    find_positive_x_side(NodeTable);
    calculate_dx(NodeTable);
    calc_topo_data(matprops_ptr);
    calc_gravity_vector(matprops_ptr);
    calc_which_son();
    return;
}

/* //for debugging purposes only, had to trick ddd into working right
 int Element::get_adapted_flag() {return adapted;} 

 void Element::put_adapted_flag(int new_adapted_status) {adapted = new_adapted_status;}
 */

//#define DEBUGLIST
#ifdef DEBUGLIST
ElemPtrList::ElemPtrList()
{   
    init(1024);
    return;
}

//! this constructor allocates space for user specified initial-size, the size_increment equals the initial size.
ElemPtrList::ElemPtrList(int initial_size)
{   
    if(initial_size==0) initial_size=1024;
    init(initial_size);
    return;
}

//! the destructor frees the list space so the programmer never has to worry about it
ElemPtrList:: ~ElemPtrList()
{   
    //printf("list_space=%d, num_elem=%d, inewstart=%d\n",list_space,num_elem,inewstart);
    free(list);
    return;
}

//! add an element pointer to the list, it will increase the size of the list by the size_increment if necessary, the size_increment is the initial size.
void ElemPtrList::add(Element* EmTemp)
{   
    if(num_elem==list_space-1)
    {   
        list_space+=size_increment;
        list=(Element**) realloc(list,list_space*sizeof(Element *));
    }

    list[num_elem]=EmTemp;
    num_elem++;
    return;
}

void ElemPtrList::init(int initial_size)
{   
    list_space=size_increment=initial_size;
    num_elem=inewstart=0;
    list=(Element **) malloc(list_space*sizeof(Element*));
    for(int i=0;i<list_space;i++) list[i]=NULL;
}

Element* ElemPtrList::get(int i)
{   return (((i>=0)&&(i<num_elem))?list[i]:NULL);};
unsigned* ElemPtrList::get_key(int i)
{   return (((i>=0)&&(i<num_elem))?list[i]->pass_key():NULL);};
int ElemPtrList::get_inewstart()
{   return inewstart;};
void ElemPtrList::set_inewstart(int inewstart_in)
{   inewstart=inewstart_in; return;};
int ElemPtrList::get_num_elem()
{   return num_elem;};

void ElemPtrList::trashlist()
{   
    for(int i=0;i<num_elem;i++) list[i]=NULL;
    num_elem=inewstart=0;
    return;
};
#endif
