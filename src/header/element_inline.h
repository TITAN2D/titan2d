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
#ifndef ELEMENT_INLINE_H
#define	ELEMENT_INLINE_H

#include "hashtab.h"





/*************************************************************************/
//getters and setters

inline const SFC_Key& Element::key() const {
    return elementsHashTable->key_[ndx_];
}

inline void Element::set_key(const SFC_Key& new_key) {
    elementsHashTable->key_[ndx_]=new_key;
}

//! returns the integer material flag for this element, needed for use of a material map which allows bedfriction to vary with physical position
inline int Element::material() const {
    return material_;
}

inline void Element::set_material(int m) {
    material_ = m;
}

//! returns the address of the first of 8 (nodes 0-7) node keys in an array, the node keys are used to access the nodes through the node hashtable
inline const SFC_Key& Element::node_key(const int i) const {
    return node_key_[i];
}

inline void Element::set_node_key(const int i, const SFC_Key& new_key) {
    node_key_[i] = new_key;
}

//! returns the pointers to the first of 8 (nodes 0-7) nodes, careful pointers can be outdated

inline Node** Element::getNodesPtrs() {
    return &(node_keyPtr[0]);
}

//! returns the pointers to the i-th of 8 (nodes 0-7) nodes, careful pointers can be outdated
inline Node* Element::getNodePtr(int i) {
    return node_keyPtr[i];
}

//! returns the pointers to the i-th of 8 (elements 0-7) elements , careful pointers can be outdated

inline Element* Element::getNeighborPtr(int i) {
    return neighborPtr[i];
}

//! not used in finite difference/volume version of titan, legacy, returns number of degrees of freedom, used is global stiffness matrices
inline int Element::ndof() const {
    return ndof_;
}

inline void Element::set_ndof(const int new_ndof) {
    ndof_ = new_ndof;
}

inline int Element::no_of_eqns() const {
    return no_of_eqns_;
}

inline void Element::set_no_of_eqns(const int new_no_of_eqns) {
    no_of_eqns_ = new_no_of_eqns;
}



//! returns this elements generation, that is how many times it's been refined -8<=generation<=+3, negative means courser than original mesh

inline int Element::generation() const {
    return generation_;
}
//! set the generation (number of times it's been refined -8<=gen<=+3) of this "element"/cell

inline void Element::generation(const int g) {
    generation_ = g;
}

//! this function returns the keys of an element's 4 brothers (an element is considered to be it's own brother) this is used during unrefinement to combine 4 brothers into their father element

inline const SFC_Key& Element::brother(const int i) const {
    return brothers_[i];
}

inline void Element::set_brother(const int i, const SFC_Key& new_key) {
    brothers_[i] = new_key;
}

//! returns the processors for the i-th neighbours of this element

inline const int& Element::neigh_proc(const int i) const {
    return neigh_proc_[i];
}
//! this function stores the processor id "proc" of neighbor "i" in the 8 element array of neighbor processors, use this function instead of putassoc.

inline void Element::set_neigh_proc(const int i, const int& proc) {
    neigh_proc_[i] = proc;
}

//! afeapi legacy not used in the finite difference/volume version of Titan, but it is used in the discontinuous galerkin version (a separate more accurate less stable implementation with a lot of things in common with the finite difference/volume code)

inline int Element::order(const int i) const {
    return order_[i];
}
//! afeapi legacy not used in the finite difference/volume version of Titan, but it is used in the discontinuous galerkin version (a separate more accurate less stable implementation with a lot of things in common with the finite difference/volume code)

inline void Element::set_order(const int i, const int ord) {
    order_[i] = ord;
}


//!only used in unrefinement
inline const SFC_Key& Element::father_by_ref() const {
    return father_;
}
//! store the father's key in the "father" variable, the "father's" key is zero until an element has been unrefined (and has not yet been deleted) it is only used in unrefinement. The getfather() member function computes the father key from "which_son" and it's nodes and is totally unrelated to the father variable.

inline void Element::set_father(const SFC_Key &fatherin) {
    father_ = fatherin;
}

//! return the element keys of this element's 4 sons, used during refinement

inline const SFC_Key& Element::son(const int i) const {
    return son_[i];
}

inline void Element::set_son(const int i, const SFC_Key& new_key) {
    son_[i] = new_key;
}

//! return the element's solution
inline double Element::el_solution(int i) const {
    return el_solution_[i];
}

inline void Element::set_el_solution(int i, double m_el_solution) {
    el_solution_[i] = m_el_solution;
}

//! returns the element's error

inline double Element::el_error(int i) const {
    return el_error_[i];
}

inline void Element::set_el_error(int i, double m_el_error) {
    el_error_[i] = m_el_error;
}




//! returns the key for this element's 8 neighbors

inline const SFC_Key& Element::neighbor(const int i) const {
    return neighbors_[i];
}
//! this function stores the key "n" of neighbor "i" in the array of the 8 keys of the neighbor keys

inline void Element::set_neighbor(const int i, const SFC_Key &new_key) {
    neighbors_[i] = new_key;
}

//! returns the pointer to this element's array of boundary conditions, not really all that important in titan since any flow that goes beyond the boundary of the GIS map leaves the computational domain.

inline BC* Element::bcptr() {
    return bcptr_;
}

inline void Element::bcptr(BC* new_bcptr) {
    bcptr_ = new_bcptr;
}
//! this function sets the pointer to an element's boundary conditions to NULL

inline void Element::void_bcptr() {
    bcptr_ = nullptr;
}

inline void Element::delete_bcptr() {
    if(bcptr_!=nullptr)
        delete bcptr_;
    void_bcptr();
}

//! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag().  refined can be permanently set to GHOST (defined in constant.h) or zero or temporarily set to 1 (with in the refinement and unrefinement routines), Keith believes it's not being unset (set from 1 to 0) when it should be after the refinement is done.  Keith believes the problem is located within H_adapt() or a function called from within it, recurse down.
//get_refined_flag
inline int Element::refined_flag() const {
    return refined_;
}
//! set this element's refined flag to i, can set it to normal (hasn't just been refined and isn't a ghost cell), "temporarily" set to "refined" (has just been refined so don't refine again), or say that it's a GHOST cell, see constant.h, (which means you don't update it, instead you get new values from the processor that owns it and you don't refine it.) refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag(). 
//put_refined_flag
inline void Element::set_refined_flag(const int i) {
    refined_ = i;
}

//! magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell.  This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
inline int Element::adapted_flag() const {
    return adapted_;
}
//! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag(). The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. These values are defined in constant.h.  The NEWSON value has allowed Keith to provide one time only immunity from unrefinement to recently refined elements, after which the "adapted" flag is resent to NOTRECADAPTED.
inline void Element::set_adapted_flag(const int new_adapted_status) {
    adapted_ = new_adapted_status;
}

//! this function returns an the generation of i-th this element's neighbors

inline int Element::neigh_gen(const int i) const {
    return neigh_gen_[i];
}
//! this function sets the ith neighbor's generation to "gen"

inline void Element::get_neigh_gen(const int i, const int gen) {
    neigh_gen_[i] = gen;
}

//! returns the which_son flag, which tells the portion of the father element that this element is physically located in

inline int Element::which_son() const {
    return which_son_;
}
//! this function sets the which_son flag when a father element is refined into its 4 sons, the which_son flag tells the portion of the father element that this element is physically located in

inline void Element::set_which_son(const int i) {
    which_son_ = i;
}

//! this function returns the vlaue of the new_old flag which is used during mesh adaptation and repartitioning

inline int Element::new_old() const {
    return new_old_;
}
//! this function sets the new or old flag, it is initialized in htflush.C and reset during repartitioning (repartition_BSFC.C and BSFC_update_and_send_elements.C)

inline void Element::set_new_old(const int i) {
    new_old_ = i;
}

//! this function returns the Load Balancing weight of an element which is used in repartitioning

inline double Element::lb_weight() const {
    return lb_weight_;
}
//! this function stores an element's load balancing weight, which is used during repartitioning

inline void Element::set_lb_weight(double dd_in) {
    lb_weight_ = dd_in;
}


//! this function returns the load balancing key, which is used during repartitioning

inline const SFC_Key& Element::lb_key() const {
    return lb_key_;
}
//! this function sets the load balancing key, which is used during repartitioning

inline void Element::set_lb_key(const SFC_Key& new_key) {
    lb_key_ = new_key;
}

//! this function copies the elmenent key to the load balancing key

inline void Element::copy_key_to_lb_key() {
    set_lb_key(key());
}


//! this function returns the process(or) id of an element, it says which processor owns the element

inline int Element::myprocess() const {
    return myprocess_;
}
//! this function sets the process(or) id of an element, it says which processor owns the element.

inline void Element::set_myprocess(const int in_proc) {
    myprocess_ = in_proc;
}

//! this function returns the opposite_brother_flag, I (Keith) am not entirely sure what this flag is for, but I know that it is used in repartioning, see BSFC_combine_elements, I think it says if an element has an opposite brother, that is, can it be combined with it's brothers to form their father

inline int Element::opposite_brother_flag() const {
    return opposite_brother_flag_;
}

inline void Element::set_opposite_brother_flag(int new_opposite_brother_flag) {
    opposite_brother_flag_ = new_opposite_brother_flag;
}


/* geoflow functions */


//! this function returns the vector of state variables
inline double Element::state_vars(int idim) const {return state_vars_[idim];}
inline void Element::state_vars(int idim, double value) {state_vars_[idim]=value;}

//!statevars by name
inline double Element::h() const {return state_vars(0);}
//!hVx for single phase or hVx_sol for two phases
inline double Element::hVx() const {
    if(elementType() == ElementType::TwoPhases)return state_vars(2);
    else return state_vars(1);
}
//!hVy for single phase or hVy_sol for two phases 
inline double Element::hVy() const {
    if(elementType() == ElementType::TwoPhases)return state_vars(3);
    else return state_vars(2);
}
inline double Element::h2() const {return state_vars(1);}

//! this function returns the vector of x and y derivatives of state variables, all the x derivatives come first as a group followed by the y derivatives as a group
inline double Element::d_state_vars(int idim) const {return d_state_vars_[idim];}
inline void Element::d_state_vars(int idim, double value) {d_state_vars_[idim]=value;}

//d_state_vars by name dh/dx, dhVx/dx, dhVy/dx, dh/dy, dhVx/dy, dhVy/dyd
inline double Element::dh_dx() const {return d_state_vars(0);}
inline double Element::dh_dy() const {return d_state_vars(NUM_STATE_VARS);}
inline double Element::dhVx_dx() const {return d_state_vars(1);}
inline double Element::dhVx_dy() const {return d_state_vars(NUM_STATE_VARS+1);}
inline double Element::dhVy_dx() const {return d_state_vars(2);}
inline double Element::dhVy_dy() const {return d_state_vars(NUM_STATE_VARS+2);}

inline double Element::dh_dx_liq() const {return d_state_vars(1);}
inline double Element::dh_dy_liq() const {return d_state_vars(NUM_STATE_VARS+1);}
inline double Element::dhVx_dx_sol() const {return d_state_vars(2);}
inline double Element::dhVx_dy_sol() const {return d_state_vars(NUM_STATE_VARS+2);}
inline double Element::dhVy_dx_sol() const {return d_state_vars(3);}
inline double Element::dhVy_dy_sol() const {return d_state_vars(NUM_STATE_VARS+3);}


//! this function returns a vector containing the previous state variables, previous mean beginning of timestep before the finite difference predictor halfstep
inline double Element::prev_state_vars(int idim) const {return prev_state_vars_[idim];}
inline void Element::prev_state_vars(int idim, double value) {prev_state_vars_[idim]=value;}

//! updates prev_states variables to current states, for first order-calculations
inline void Element::update_prev_state_vars(){for (int i = 0; i < NUM_STATE_VARS; i++)prev_state_vars(i, state_vars(i));}

//! this function returns the length of an element in the x and y directions
inline double Element::dx(int idim) const {return dx_[idim];}
inline void Element::dx(int idim, double value){dx_[idim]=value;}

//! this function returns which side of the element is facing the positive x direction

inline int Element::positive_x_side() const {
    return positive_x_side_;
}

inline void Element::set_positive_x_side(const int new_positive_x_side) {
    positive_x_side_ = new_positive_x_side;
}



//! this function calculates the maximum x and y direction wavespeeds which are the eigenvalues of the flux jacobian
inline double Element::eigenvxymax(int idim) const {return eigenvxymax_[idim];}
inline double& Element::eigenvxymax_ref(int idim) {return eigenvxymax_[idim];}
inline void Element::eigenvxymax(int idim,double value) {eigenvxymax_[idim]=value;}

inline double Element::shortspeed() {
    return shortspeed_;
}

inline void Element::set_shortspeed(double shortspeedin) {
    shortspeed_ = shortspeedin;
}


//! this function returns the already calculated value(s) of k active passive, which comes from using th Coulomb friction model of granular flows (this is problem specific to titan and thus does not appear in the standard afeapi code)
inline double& Element::kactxy_ref(const int idim){return kactxy_[idim];}
inline double Element::kactxy(const int idim) const {return kactxy_[idim];}

inline void Element::kactxy(const int idim, double value) {kactxy_[idim]=value;}
//! interface to change value of earth-pressure coefficients
inline void Element::set_kactxy(double kap[DIMENSION]) {for (int i = 0; i < DIMENSION; i++)kactxy_[i] = kap[i];}

//! this function returns the precomputed elevation
inline double& Element::elevation_ref() {
    return elevation_;
};

inline double Element::elevation() const {
    return elevation_;
};

inline void Element::set_elevation(double new_elevation) {
    elevation_ = new_elevation;
};

//! this function returns the x and y slopes of the terrain elevation
inline double Element::zeta(int idim) const {return zeta_[idim];}
inline double& Element::zeta_ref(int idim) {return zeta_[idim];}
inline void Element::zeta(int idim, double value) {zeta_[idim]=value;}

//! returns the already computed gravity vector in local coordinates, the local z direction is normal to the terrain surface and the projection of the local x and y components into the horizontal plane are aligned with global x (UTM E) and y (UTM N) directions.
inline double Element::gravity(int idim) const {return gravity_[idim];}
inline void Element::gravity(int idim, double value) {gravity_[idim]=value;}

//! this function returns the precomputed derivatives of the z component of gravity, this is a purely terrain geometry dependant derivative, that is little diffent than curvature
inline double Element::d_gravity(int idim) const {return d_gravity_[idim];}
inline void Element::d_gravity(int idim, double value) {d_gravity_[idim]=value;}

//! this function returns the precomputed local terrain curvature.  Curvature itself is the inverse of radius of curvature.  The exact value of curvature  is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation. 
inline double Element::curvature(int idim) const {return curvature_[idim];}
inline double& Element::curvature_ref(int idim) {return curvature_[idim];}
inline void Element::curvature(int idim, double value) {curvature_[idim]=value;}


//! this function returns the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
inline int Element::elm_loc(int idim) const {return elm_loc_[idim];}
//! this function sets the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
inline void Element::set_elm_loc(const int idim, const int int_in){elm_loc_[idim] = int_in;}
//void put_elm_loc(int* int_in){for(i=0;i<DIMENSION;i++)elm_loc[i] = int_in[i];}


//! this function returns the precomputed and scaled coordinates of this element (which would be the same as its bubble node's coordinates)
inline double Element::coord(int idim) const {return coord_[idim];}
inline void Element::set_coord(int idim, double new_crd){coord_[idim]=new_crd;}


//! this function returns the value of "stoppedflags"
inline int Element::stoppedflags() {
    return stoppedflags_;
}
//! this function is used to assign a value to stopped flags, for when you don't want to compute the criteria to decide whether it's stopped or not, useful during developement

inline void Element::set_stoppedflags(int stoppedflagsin) {
    stoppedflags_ = stoppedflagsin;
}

//! this function returns the stored value of the extrusion (out of the ground) fluxes in this element
inline double Element::Influx(int idim) const {return Influx_[idim];}
inline void Element::Influx(int idim, double value) {Influx_[idim]=value;}
//! this function zeros the extrusion (out of the ground) fluxes in this element
inline void Element::zero_influx(){for (int i = 0; i < NUM_STATE_VARS; i++)Influx(i, 0.0);}



//! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong
inline int Element::ithelem() {return ithelem_;}

//! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.

inline void Element::set_ithelem(int i) {ithelem_ = i;}


//! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus the effective bed friction angle holds either the value of the actual bed friction angle if it should not be stopped or the value of the internal friction angle if it should not be stopped
inline double& Element::effect_bedfrict_ref() {return effect_bedfrict_;}
inline double* Element::effect_bedfrict_ptr() {return &effect_bedfrict_;}
inline double Element::effect_bedfrict() {return effect_bedfrict_;}
inline void Element::set_effect_bedfrict(double new_effect_bedfrict) {effect_bedfrict_ = new_effect_bedfrict;}

inline double Element::effect_tanbedfrict() {return effect_tanbedfrict_;}
inline void Element::set_effect_tanbedfrict(double new_effect_tanbedfrict) {effect_tanbedfrict_ = new_effect_tanbedfrict;}

//! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, if the effective bed friction angle equals the internal friction angle effect_kactxy takes on the value 1, k active/passive comes from using a Coulomb friction model for granular flows
inline double Element::effect_kactxy(int idim) const {return effect_kactxy_[idim];}
inline void Element::effect_kactxy(int idim, double value) {effect_kactxy_[idim]=value;}

//! this inline member function returns the stored value of Awet, Awet is the fraction of an element's area that is wet (has material), 0.0<=Awet<=1.0, where there is no flow (pileheight < GEOFLOW_TINY) Awet=0, in the interior of the Flow Awet=1.0, at the boundary of the flow, elements will be PARTIALLY WET (i.e. where the element SHOULD be separated into a dry part and wet part), Awet is the fraction that should be wet, Awet is updated during the corrector part of the (finite difference)predictor-(finite volume)corrector update.  Fluxes are adjusted to acount for how wet/dry an edge of an element is. Keith wrote this may 2007
inline double Element::Awet() const {return Awet_;}
//! this inline member function assigns a value to Awet, Awet is the fraction of an element's area that is wet (has material), 0.0<=Awet<=1.0, where there is no flow (pileheight < GEOFLOW_TINY) Awet=0, in the interior of the Flow Awet=1.0, at the boundary of the flow, elements will be PARTIALLY WET (i.e. where the element SHOULD be separated into a dry part and wet part), Awet is the fraction that should be wet, Awet is updated during the corrector part of the (finite difference)predictor-(finite volume)corrector update.  Fluxes are adjusted to acount for how wet/dry an edge of an element is. Keith wrote this may 2007
inline void Element::set_Awet(double Awet_in) {Awet_ = Awet_in;}

//! this inline member function returns the stored value of Swet.  Swet is the fraction of the element's partially wet sides that are wet (i.e. have material).  Where there is no flow (pileheight < GEOFLOW_TINY), Swet=0.  In the interior of a flow, Swet=1.0.  At the flow boundary, elements will be PARTIALLY WET, 0.0<=Swet<=1.0.  Due to symmetry, any element can only have 0 or 2 partially wet sides, each of which (for normalized elements) will have the same fraction that is wet, Swet.  Swet for each partially wet cell is updated every time-step when calc_wet_dry_orient() is called in step.C.  Fluxes are adjusted to account for how wet/dry an edge of an element is.  Keith wrote this function may 2007
inline double Element::Swet() const {return Swet_;}
//! this inline member function assigns a value to Swet.  Swet is the fraction of the element's partially wet sides that are wet (i.e. have material).  Where there is no flow (pileheight < GEOFLOW_TINY), Swet=0.  In the interior of a flow, Swet=1.0.  At the flow boundary, elements will be PARTIALLY WET, 0.0<=Swet<=1.0.  Due to symmetry, any element can only have 0 or 2 partially wet sides, each of which (for normalized elements) will have the same fraction that is wet, Swet.  Swet for each partially wet cell is updated every time-step when calc_wet_dry_orient() is called in step.C.  Fluxes are adjusted to account for how wet/dry an edge of an element is.  Keith wrote this function may 2007
inline void Element::set_Swet(double Swet_in) {Swet_ = Swet_in;}

//! this inline member function returns the value of iwetnode.  iwetnode is an integer that defines which of an element's 9 nodes is its "wettest" node (wet elements are those containing material).  In the interior of a flow, iwetnode=8 (the center node), indicating a fully wet element.  Outside of a flow (where an element and all it's neighbors have pileheight < GEOFLOW_TINY), iwetnode is also 8.  Along a flow boundary, partially wet elements with 1,2, or 3 wet sides can have an iwetnode other than 8.   iwetnode is used to determine which side of the dryline in a partially wet element has material.  Keith wrote this function may 2007
inline int Element::iwetnode() const {return iwetnode_;}
//! this inline member function sets the value of iwetnode. iwetnode is an integer that defines which of an element's 9 nodes is its "wettest" node (wet elements are those containing material).  In the interior of a flow, iwetnode=8 (the center node), indicating a fully wet element.  Outside of a flow (where an element and all it's neighbors have pileheight < GEOFLOW_TINY), iwetnode is also 8.  Along a flow boundary, partially wet elements with 1,2, or 3 wet sides can have an iwetnode other than 8.   iwetnode is used to determine which side of the dryline in a partially wet element has material.  Keith wrote this function may 2007
inline void Element::set_iwetnode(int iwetnode_in) {iwetnode_ = iwetnode_in;}

//! this inline member function returns the array "drypoint".  drypoint[0] is the local x-coordinate, and drypoint[1] the local y-coordinate of its namesake, which is used to specify the position of the flow-front (or dryline) inside a given element.  The position of the dryline along with iwetnode is used to determine Awet, i.e. which fraction of a partially wet element is wet (contains material).  Keith wrote this function may 2007
inline double Element::drypoint(int idim) const {return drypoint_[idim];}
//! this inline member function sets the values of the array "drypoint".  drypoint[0] is the local x-coordinate, and drypoint[1] the local y-coordinate of its namesake, which is used to specify the position of the flow-front (or dryline) inside a given element.  The position of the dryline along with iwetnode is used to determine Awet, i.e. which fraction of a partially wet element is wet (contains material).  Keith wrote this function may 2007  
inline void Element::drypoint(int idim, double value) {drypoint_[idim]=value;}


//! sgn of double @TODO replace this sign function]
inline double Element::sgn(double a) {
    return (a < 0.0 ? -1.0 : 1.0);
}

inline const ElementType& Element::elementType() const {
    return elementsHashTable->elementType_;
}

inline void Element::elementType(const ElementType& new_element_type) {
    elementsHashTable->elementType_ = new_element_type;
}

//index in storage
inline ti_ndx_t Element::ndx() const {return ndx_;}
inline void Element::ndx(ti_ndx_t new_ndx) {ndx_ = new_ndx;}

/*************************************************************************/
inline void Element::put_height_mom(double pile_height, double volf, double xmom, double ymom) {
    prev_state_vars(0, pile_height);
    state_vars(0, pile_height);
    prev_state_vars(1, pile_height * volf);
    state_vars(1, pile_height * volf);
    prev_state_vars(2, xmom);
    state_vars(2, xmom);
    prev_state_vars(3, ymom);
    state_vars(3, ymom);
    if (pile_height > GEOFLOW_TINY) {
        set_shortspeed(sqrt(xmom * xmom + ymom * ymom) / (pile_height * volf));
        set_Awet(1.0);
    } else {
        set_shortspeed(0.0);
        set_Awet(0.0);
    }
    return;
};

inline void Element::put_height_mom(double pile_height, double xmom, double ymom) {
    prev_state_vars(0, pile_height);
    state_vars(0, pile_height);
    prev_state_vars(1, xmom);
    state_vars(1, xmom);
    prev_state_vars(2, ymom);
    state_vars(2, ymom);
    if (pile_height > GEOFLOW_TINY) {
        set_shortspeed(sqrt(xmom * xmom + ymom * ymom) / pile_height);
        set_Awet(1.0);
    } else {
        set_shortspeed(0.0);
        set_Awet(0.0);
    }

    return;
};

inline void Element::put_height(double pile_height) {
    if (elementType() == ElementType::TwoPhases) {
        put_height_mom(pile_height, 1., 0., 0.);
    }
    if (elementType() == ElementType::SinglePhase) {
        put_height_mom(pile_height, 0.0, 0.0);
    }
    return;
};

inline void Element::set_sons(const SFC_Key* s) {
    for (int i = 0; i < 4; i++)
        set_son(i, s[i]);

    set_refined_flag(1);
    set_adapted_flag(OLDFATHER);
}

inline void Element::set_brothers(const SFC_Key* s) {
    for (int i = 0; i < 4; i++)
        set_brother(i, s[i]);
}

inline void Element::putel_sq(double solsq, double errsq) {
    set_el_solution(0, solsq);
    set_el_error(0, errsq);
}


inline void Element::update_neighbors_nodes_and_elements_pointers(ElementsHashTable* El_Table, NodeHashTable* NodeTable) {
    int i;
    if (El_Table != NULL) {
        for (i = 0; i < 8; i++) {
            neighborPtr[i] = (Element*) El_Table->lookup(neighbor(i));
        }
    }
    if (NodeTable != NULL) {
        for (i = 0; i < 8; i++) {
            node_keyPtr[i] = (Node*) NodeTable->lookup(node_key(i));
        }
    }
    return;
}

inline int Element::check_neighbors_nodes_and_elements_pointers(ElementsHashTable* El_Table, NodeHashTable* NodeTable) {
    int i;
    int count = 0;
    if (El_Table != NULL) {
        for (i = 0; i < 8; i++) {
            if (neighborPtr[i] != (Element*) El_Table->lookup(neighbor(i)))
                count++;
        }
    }
    if (NodeTable != NULL) {
        for (i = 0; i < 8; i++) {
            if (node_keyPtr[i] != (Node*) NodeTable->lookup(node_key(i)))
                count++;
        }
    }
    return count;
}


#endif	/* ELEMENT_INLINE_H */

