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
 * $Id: element2.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifndef ELEMENT_H
#define ELEMENT_H
#include <math.h>
#include "boundary.h"
#include "hashtab.h"
#include "node.h"
#include "struct.h"

#include <fstream>
#include <iostream>
using namespace std;

//#define USE_FATHER

//! The Element class is a data structure designed to hold all the information need for an h (cell edge length) p (polynomial order) adaptive finite element.  Titan doesn't use p adaptation because it is a finite difference/volume code, hence many of the members are legacy from afeapi (adaptive finite element application programmers interface) which serves as the core of titan.  There is a seperate Discontinuous Galerkin Method (finite elements + finite volumes) version of titan and the polynomial information is not legacy there.  However in this version of Titan elements function simply as finite volume cells.
class Element{

  friend class HashTable;

  friend void AssertMeshErrorFree(HashTable *El_Table, HashTable* NodeTable,
				  int numprocs, int myid, double loc);

  friend void ElemBackgroundCheck(HashTable* El_Table, HashTable* NodeTable,
				  unsigned *debugkey, FILE *fp);

  friend void ElemBackgroundCheck2(HashTable* El_Table, HashTable* NodeTable,
				   void *EmDebug, FILE *fp);

  friend void NodeBackgroundCheck(HashTable* El_Table, HashTable* NodeTable,
				  unsigned *debugkey, FILE *fp);

  friend void delete_oldsons(HashTable* El_Table, HashTable* NodeTable,
			     int myid, void *EmFather);

  friend void refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, 
				  int numprocs, int myid, void* RefinedList,
				  TimeProps* timeprops_ptr);

  friend void unrefine_neigh_update(HashTable* El_Table, HashTable* NodeTable,
				    int myid, void* NewFatherList);
  
  friend void  unrefine_interp_neigh_update(HashTable* El_Table, 
					    HashTable* NodeTable, int nump, 
					    int myid, void* OtherProcUpdate);

  friend void BSFC_combine_elements(int side, Element *EmTemp, 
				    HashTable *HT_Elem_Ptr, 
				    HashTable *HT_Node_Ptr, 
				    int destination_proc);
 
  friend void Pack_element(void *sendel, ElemPack* elem,
			   HashTable* HT_Node_Ptr, int destination_proc);

  friend void destroy_element(void  *r_element, HashTable* HT_Elem_Ptr, 
			      HashTable* HT_Node_Ptr, int target_pro, ELinkPtr* EL_head);
  
  friend void create_element(ElemPack* elem2, HashTable* HT_Elem_Ptr, 
			     HashTable* HT_Node_Ptr, int myid, double* e_error);

  friend void construct_el(Element* newelement, ElemPack* elem2,
			   HashTable* HT_Node_Ptr, int myid, double* e_error);

 public:

  //! default constructor, does nothing except set stoppedflags=2, this should never be used
  Element() {
    counted=0;
    father[0]=father[1]=0; //initialize the father key to zero
    for (int i=0; i<NUM_STATE_VARS; i++)
    {
      state_vars[i] = -1;
      Influx[i] =0.;
    }
    adapted=TOBEDELETED;
    refined=1;
    Awet=0.0;
    Swet=1.0;
    drypoint[0]=drypoint[1]=0.0;
    iwetnode=8;

    stoppedflags=2; //material in all elements start from rest
    //do_erosion=-1;
  };  
  
  //! constructor that creates an original element when funky is read in
  Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], 
	  int n_pro[], BC *b, int mat, int *elm_loc_in, 
	  double pile_height, int myid, unsigned *opposite_brother);

  //! constructor that creates a son element from its father during refinement
  Element(unsigned nodekeys[][KEYLENGTH], unsigned neigh[][KEYLENGTH], 
	  int n_pro[], BC *b, int gen, int elm_loc_in[], int *ord, 
	  int gen_neigh[], int mat, Element *fthTemp, double *coord_in,
	  HashTable *El_Table, HashTable *NodeTable, int myid, 
	  MatProps *matprops_ptr,
	  int iwetnodefather, double Awetfather, double *drypoint_in);

  //! constructor that creates a father element from its four sons during unrefinement
  Element(Element *sons[], HashTable *NodeTable, HashTable *El_Table, 
	  MatProps *matprops_ptr);

  //! constructor that creates/restores a saved element during restart
  Element(FILE* fp, HashTable* NodeTable, MatProps* matprops_ptr, int myid);

  //! destructor that does nothing except delete boundary condition pointer
  ~Element();

  //! this member function saves a single element to a file with a single fwrite call, this allows the element to be recreated/restored upon restart of a simulation
  void      save_elem(FILE* fp, FILE* fptxt); //for restart
  
  //! returns address of element (same as bubble node, node 8 out of 0->8) hashtable key
  unsigned* pass_key();

  //! returns the integer material flag for this element, needed for use of a material map which allows bedfriction to vary with physical position
  int       get_material();

  //! legacy afeapi function prototype, this function does not exist in the finite difference/volume version of Titan
  void      get_stiffness(HashTable*, HashTable*, double*, double*, Element*);

  //! returns the address of the first of 8 (nodes 0-7) node keys in an array, the node keys are used to access the nodes through the node hashtable
  unsigned* getNode();

  //! returns the array of 8 processors for the 8 neigbors of this element
  int*      getassoc();

  //! not used in finite difference/volume version of titan, legacy, returns number of degrees of freedom, used is global stiffness matrices
  int       get_no_of_dof();

  //! set the generation (number of times it's been refined -8<=gen<=+3) of this "element"/cell
  void      put_gen(int);

  //! store the keys for the four son "elements" in the father element, used temporarily during refinement
  void      putson(unsigned*);

  //! when a father element is refined into 4 son elements, the 4 son elements are "brothers" (they can be recombined into the father), this function stores the keys of all four brothers in one of them, it should be called 4 times one for each brother
  void      putbrothers(unsigned*);

  //! this function returns the keys of an element's 4 brothers (an element is considered to be it's own brother) this is used during unrefinement to combine 4 brothers into their father element
  unsigned* get_brothers();

  //! this function stores the processor id "a" of neighbor "i" in the 8 element array of neighbor processors, this functionality is duplicated by put_neigh_proc which is the preferred function to use (don't use this one it's legacy)
  void      putassoc(int a, int i);

  //! this function stores the key "n" of neighbor "i" in the array of the 8 keys of the neighbor keys
  void      putneighbor(unsigned *n, int i);

  //! this function stores the processor id "proc" of neighbor "i" in the 8 element array of neighbor processors, use this function instead of putassoc.
  void      put_neigh_proc(int i, int proc);

  //! afeapi legacy not used in the finite difference/volume version of Titan, but it is used in the discontinuous galerkin version (a separate more accurate less stable implementation with a lot of things in common with the finite difference/volume code)
  void      put_order(int, int);

  //! afeapi legacy not used in the finite difference/volume version of Titan, but it is used in the discontinuous galerkin version (a separate more accurate less stable implementation with a lot of things in common with the finite difference/volume code)
  int*      get_order();

  //! find and return what the key of this element's father element would be, very simple since the bubble node has the same key as the element, so all this function does is find which of its corner nodes will be the father element's bubble node, which it knows since it knows which_son it is.  
  unsigned* getfather();

  //! store the father's key in the "father" variable, the "father's" key is zero until an element has been unrefined (and has not yet been deleted) it is only used in unrefinement. The getfather() member function computes the father key from "which_son" and it's nodes and is totally unrelated to the father variable.
  void put_father(unsigned fatherin[KEYLENGTH]);

  //! return the element keys of this element's 4 sons, used during refinement
  unsigned* getson();

  //! stores the ?square? of the "solution" and solution error, used durring refinement
  void      putel_sq(double solsq, double ellsq);

  //! return the element's solution vector
  double*   get_el_solution();

  //! returns the element's error vector
  double*   get_el_error();

  //! returns the array of keys for this element's 8 neighbors
  unsigned* get_neighbors();

  //! returns the array of processor ids for this element's 8 neighbors
  int*      get_neigh_proc();

  //! returns the pointer to this element's array of boundary conditions, not really all that important in titan since any flow that goes beyond the boundary of the GIS map leaves the computational domain.
  BC*       get_bcptr();

  //! returns this elements generation, that is how many times it's been refined -8<=generation<=+3, negative means courser than original mesh
  int       get_gen();

  //! compare the FindNeigh key against the keys of this element's 8 neighbors to determine which if any neighbor FindNeigh is
  int       which_neighbor(unsigned *FindNeigh);

  //! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag().  refined can be permanently set to GHOST (defined in constant.h) or zero or temporarily set to 1 (with in the refinement and unrefinement routines), Keith believes it's not being unset (set from 1 to 0) when it should be after the refinement is done.  Keith believes the problem is located within H_adapt() or a function called from within it, recurse down.
  int       get_refined_flag(); 

  //! magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell.  This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
  int       get_adapted_flag();


  //! call this function after this element's neighbor(s) have been refined, proc is processor id for neighbor[which_side+4]
  void      change_neighbor(unsigned *newneighbs, int which_side, int proc, int reg);

  //! set this element's refined flag to i, can set it to normal (hasn't just been refined and isn't a ghost cell), "temporarily" set to "refined" (has just been refined so don't refine again), or say that it's a GHOST cell, see constant.h, (which means you don't update it, instead you get new values from the processor that owns it and you don't refine it.) refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag(). 
  void      put_refined_flag(int i);

  //! refined, get_refined_flag(), put_refined_flag() are the partly replaced predecessors of adapted, get_adapted_flag(), and put_adapted_flag(). The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. These values are defined in constant.h.  The NEWSON value has allowed Keith to provide one time only immunity from unrefinement to recently refined elements, after which the "adapted" flag is resent to NOTRECADAPTED.
  void      put_adapted_flag(int new_adapted_status);


  //! this function is only called in htflush.C to initialize it to zero, I (Keith) think it is afeapi legacy that isn't being used anymore but I'm not sure about that 
  void      put_send_flag(int, int);

  //! this function is only called in htflush.C to initialize it to zero, I (Keith) think it is afeapi legacy that isn't being used anymore but I'm not sure about that 
  void      put_recv_flag(int, int); 

  //! this function isn't being called anywhere, which means it is afeapi legacy
  int       get_send_flag(int); 

  //! this function isn't being called anywhere, which means it is afeapi legacy
  int       get_recv_flag(int);

  //! this function returns an array holding the generation of all 8 of this element's neighbors
  int*      get_neigh_gen();

  //! this function sets the ith neighbor's generation to "gen"
  void      put_neigh_gen(int i, int gen);

  //! this function sets the which_son flag when a father element is refined into its 4 sons, the which_son flag tells the portion of the father element that this element is physically located in
  void      put_which_son(int);

  //! returns the which_son flag, which tells the portion of the father element that this element is physically located in
  int       get_which_son();

  //! this function calculates the which_son flag for the original (or restored in case of a restart) element.  It also calculates which son of the grandfather element the father is durring unrefinement.
  void      calc_which_son();

  //! this function sets the new or old flag, it is initialized in htflush.C and reset during repartitioning (repartition_BSFC.C and BSFC_update_and_send_elements.C)
  void      put_new_old(int i);

  //! this function returns the vlaue of the new_old flag which is used during mesh adaptation and repartitioning
  int       get_new_old();
  
  //! this function is legacy afeapi code, the function is defined in element2.C but it is never called anywhere in the finite difference/volume version of titan because it's finite element (including Discontinuous Galerkin) specific code
  void      update_ndof();

  //! this function is called during repartitioning when one of an element's neighbors is sent to another processor
  void      change_neighbor_process(int which, int newp);
  
  //! the function returns the vector of element "error", element error is used to say when a function should be refined
  double*   get_el_err();

  //! this function is afeapi legacy it is not called anywhere in the finite difference/volume version of titan
  void      get_nelb_icon(HashTable*, HashTable*,int*,int*);

  //! this function sets the pointer to an element's boundary conditions to NULL
  void      void_bcptr();

  //! this function returns the Load Balancing weight of an element which is used in repartitioning
  double    get_lb_weight();

  //! this function stores an element's load balancing weight, which is used during repartitioning
  void      put_lb_weight(double dd_in);

  //! this function returns the load balancing key, which is used during repartitioning
  unsigned* get_lb_key();

  //! this function sets the load balancing key, which is used during repartitioning
  void      put_lb_key(unsigned* in_key);

  //! this function copies the elmenent key to the load balancing key
  void      copy_key_to_lb_key();

  //! this function sets the process(or) id of an element, it says which processor owns the element.
  void      put_myprocess(int in_proc);

  //! this function returns the process(or) id of an element, it says which processor owns the element
  int       get_myprocess();

  //! this function returns the opposite_brother_flag, I (Keith) am not entirely sure what this flag is for, but I know that it is used in repartioning, see BSFC_combine_elements, I think it says if an element has an opposite brother, that is, can it be combined with it's brothers to form their father
  int       get_opposite_brother_flag();

  //! this function computes searches for an element's brother, i.e. the brother (son of the same father) that is located diagonally from it, to get the brother information requires that atleast one of this element's neighboring brothers is on this process in order to get information onthe brother that is not a neighbor
  void      find_opposite_brother(HashTable*);
  //void      get_icon(HashTable*, HashTable*, int[4]);
  //void      get_boundary(int[4], double[4]);
  /* geoflow functions */

  //! this function initializes pileheight, momentums and shortspeed (also known as the L'Hosptial speed see calc_shortspeed for an explanation),this function is called in init_piles.C 
  void      put_height_mom(double pile_height, double vfract, double xmom, double ymom);

  //! this function assigns a specified value to the pileheight and zeros to the momentums and shortspeed
  void      put_height(double pile_height);

  //! this function returns the vector of state variables
  double*   get_state_vars();

  //! this function returns the vector of x and y derivatives of state variables, all the x derivatives come first as a group followed by the y derivatives as a group
  double*   get_d_state_vars();

  //! this function returns the x and y slopes of the terrain elevation
  double*   get_zeta();

  //! this function returns the length of an element in the x and y directions
  double*   get_dx();

  //! this function computes which side of the element is facing the positive x direction
  void      find_positive_x_side(HashTable*);

  //! this function returns which side of the element is facing the positive x direction
  int       get_positive_x_side();

  //! this function computes the x and y derivatives of the state variables
  void      get_slopes(HashTable*, HashTable*, double);

  //! this function returns a vector containing the previous state variables, previous mean beginning of timestep before the finite difference predictor halfstep
  double*   get_prev_state_vars();

  //! updates prev_states variables to current states, for first order-calculations
  void      update_prev_state_vars();

  //! this function calculates the lengths of the element in the (global) x and y directions
  void      calculate_dx(HashTable* NodeTable);

  //! this function assigns the element's coordinates to be its bubble node's coordinates
  void      insert_coord(HashTable* NodeTable);

  //! this function, based on the dir flag, chooses between calling xdirflux and ydirflux, which respectively, calculate either the x or y direction analytical cell center fluxes (or the fluxes at the the boundary if 2nd order flux option is checked on the gui). Keith wrote this.
  void      zdirflux(HashTable* El_Table, HashTable* NodeTable, MatProps* matprops_ptr, int order_flag, int dir, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS], Element* EmNeigh, double dt);

  //! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) x direction fluxes. Keith wrote this
  void      xdirflux(MatProps* matprops_ptr, double dz, double thissideSwet, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

  //! this function calculates the analytical cell center (or cell boundary if 2nd order flux flag is checked on the gui) y direction fluxes. Keith wrote this
  void      ydirflux(MatProps* matprops_ptr, double dz, double thissideSwet, double hfv[3][NUM_STATE_VARS], double hrfv[3][NUM_STATE_VARS]);

  //! this function (indirectly) calculates the fluxes that will be used to perform the finite volume corrector step and stores them in element edge nodes, indirectly because it calls other functions to calculate the analytical fluxes and then calls another function to compute the riemann fluxes from the analytical fluxes. Talk to me (Keith) before you modify this, as I am fairly certain that it is now completely bug free and parts of it can be slightly confusing.
  void      calc_edge_states(HashTable* El_Table, HashTable* NodeTable,
			     MatProps* matprops_ptr, int myid, double dt, 
			     int* order_flag, double *outflow);
  
  /*! this function calculates the maximum x and y direction wavespeeds 
   *  which are the eigenvalues of the flux jacobian
   */
  double*   get_eigenvxymax();

  /*! 
   * this function calculates the shortspeed,also known as the L'Hospital 
   * (pronounced Loo-pee-tal, you can look up L'Hospital's rule in almost 
   * any calculus book if you so desire). here is a brief explanation of 
   * shortspeed: shortspeed=|v|=|dhv/dh|=|v*dh/dh+h*dv/dh|=|v+h*dv/dh| which 
   * goes to |v| in the limit of h->0, this is a more accurate way to compute 
   * speed when the pile in this cell is short, hence the name "shortspeed" 
   * but it is not accurate when the pile is tall, that is when h*dv/dh is 
   * large, Keith implemented this in late summer 2006 
   */
  void      calc_shortspeed(double inv_dt);

  //! this function returns the already computed shortspeed
  double    get_shortspeed();

  //! this function assigns the value passed in to shortspeed
  void      put_shortspeed(double shortspeedin);

  /*! this function computes the velocity, either V=hV/h or shortspeed in 
   *  the direction of hV/h, if the pile is short, that is h is less than 
   *  the defined (nondimensional) value of GEOFLOW_SHORT, see geoflow.h, 
   *  it chooses the speed to be min(|hV/h|,shortspeed) if h is greater than 
   *  GEOFLOW_SHORT it chooses hV/h regardless of which one is smaller. 
   */
  void   eval_velocity(double xoffset, double yoffset, double Vel[]);

  //! this function is legacy afeapi code, it is never called in the finite difference/volume version of titan
  double*   get_coefABCD() {return coefABCD;};

  //! this function returns the already calculated value(s) of k active passive, which comes from using th Coulomb friction model of granular flows (this is problem specific to titan and thus does not appear in the standard afeapi code)
  double*   get_kactxy();

  //! returns the already computed gravity vector in local coordinates, the local z direction is normal to the terrain surface and the projection of the local x and y components into the horizontal plane are aligned with global x (UTM E) and y (UTM N) directions.
  double*   get_gravity();

  //! this function is titan legacy code it is defined in Element2.C but is not called anywhere
  int       determine_refinement(double);

  //! this function returns the precomputed elevation
  double    get_elevation();

  //! this function returns the precomputed derivatives of the z component of gravity, this is a purely terrain geometry dependant derivative, that is little diffent than curvature
  double*   get_d_gravity();

  //! this function returns the precomputed local terrain curvature.  Curvature itself is the inverse of radius of curvature.  The exact value of curvature  is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation. 
  double*   get_curvature();

  //! the only place this function is called is in move_data.C, I believe it is legacy titan code and could probably be removed
  void      put_lam(double lam_in);

  //! this function is never called in the finite difference/volume version of titan, it is legacy (I believe titan rather than afeapi legacy code) that probably should be removed
  double    get_lam();

  //! this function is called in element_weight.C, it is used in computing the load balancing weight
  void      calc_flux_balance(HashTable *NodeTable);

  //! this function calculates topographic data, it calls GIS commands to compute elevation, slopes, and curvatures from the GIS map and scales them appropriately 
  void      calc_topo_data(MatProps *matprops_ptr);

  //! this function calculates the (global) x and y derivatives of the local z component of gravity as an approximation of the local derivatives, it wouldn't be that difficult to correct incorporating the terrain slopes in the calculation it is calculated in the creation of a father element, after mesh refinement and, during a restart.
  void      calc_d_gravity(HashTable *El_Table);

  //! this function calculates the gravity vector in local coordinates
  void      calc_gravity_vector(MatProps *matprops_ptr);
  
   //! this function is defined in unrefine.C, it is also called in that file, it finds this element's brothers
  int find_brothers(HashTable* El_Table, HashTable* NodeTable, 
		    double target, int myid, MatProps* matprops_ptr,
		    void* NewFatherList, void* OtherProcUpdate);
  /*
  //! this function is defined in unrefine.C, it is also called in that file, it finds this element's brothers
  int find_brothers(HashTable* El_Table, HashTable* NodeTable, 
		    double target, int myid, MatProps* matprops_ptr,
		    Element **NewFatherList, int* NumNewFathers,
		    Element **OtherProcUpdate, int *NumOtherProcUpdate); 
  */

  //! this function is defined in unrefine.C, it is also called in that file and no where else, it prevents refinement when one or more of the brothers does not belong to this processor
  int check_unrefinement(HashTable *El_Table, double target);

  //! this function updates this elements neighbor info when one of its neighbors has been unrefined
  void change_neigh_info(unsigned *fth_key, unsigned *ng_key, int neworder, int ng_gen, int fth_proc);

  //! this function returns the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
  int* get_elm_loc();
  
  //! this function sets the elm_loc variable, which is used in unrefinement beyond the initial coarse grid
  void put_elm_loc(int* int_in);

  //! this function returns the precomputed and scaled coordinates of this element (which would be the same as its bubble node's coordinates)
  double* get_coord();

  //! this function stores the coordinates of this element (which would be the same as its bubble node's coordinates)
  void put_coord(double* coord_in);

  //! this function is part of the experimental _LOCAL_ (not Bin Yu's) stopping criteria which has not yet been validated, I (Keith) have faith in the criteria, but enforcing the stopped after it has been decided that it needs to stop still needs some work. the only place this function is called is in get_coef_and_eigen.C immediately after k_active/passive and in init_piles.C when computing the initial volume that "should already be" deposited.
  void calc_stop_crit(MatProps*);

  //! this function is used to assign a value to stopped flags, for when you don't want to compute the criteria to decide whether it's stopped or not, useful during developement
  void put_stoppedflags(int stoppedflagsin);

  //! interface to change value of earth-pressure cofficients
  void put_kactxy (double kap[]){ kactxy[0]=kap[0]; kactxy[1]=kap[1]; }

  //! this function returns the value of "stoppedflags"
  int  get_stoppedflags();
  
  //! this function zeros the extrusion (out of the ground) fluxes in this element
  void zero_influx();

  //! this function returns the stored value of the extrusion (out of the ground) fluxes in this element
  double *get_influx();

  //! this function calculates the extrusion (out of the ground) fluxes for this elements
  void calc_flux(HashTable *NodeTable,FluxProps *fluxprops, 
		 TimeProps *timeprops);

  //! this function returns 2 if this element contains pileheight>=contour_height and has a neighbor who contains pileheight<contour_height.  It returns 1 if this element contains pileheight<contour_height and has a neighbor who contains pileheight>=contour_height.  It returns 0 otherwise. The intended use if if(EmTemp->if_pile_boundary(ElemTable,contour_height)) but I (Keith) added the distinction bewteen 1 and 2 to allow future developers to distinguish between the inside and outside of a pileheight contour line, as this functionality could be useful in the future.
  int if_pile_boundary(HashTable *ElemTable, double contour_height);

  //! this function returns 2 if this element has Influx[0]>0 and has a neighbor who has Influx[0]<=0.  It returns 1 if this element has Influx[0]==0 and has a neighbor who has Influx[0]!=0.  It returns -1 if this element has Influx[0]<0 and a neighbor with Influx[0]>=0. It returns 0 otherwise. Influx[0] is a pileheight per unit time source term.  Currently Influx[0] is restricted to be non-negative (a source or no source with sinks not allowed), but I (Keith) have added the extra functionality because it may be useful at a future date. The intended use if if(EmTemp->if_source_boundary(ElemTable)), but the distinction between 1 and 2 allows futuredevelopers to distinguish between the strictly inside and strictly outside of an area with a flux source term.
  int if_source_boundary(HashTable *ElemTable);

  //! the buffer layer is a layer of refined cells on the outside of the pile, i.e. ((pileheight<contour_height)&&(Influx[0]==0)) and adjacent to the pile.  It is "N" elements wide, and the "N" element width is increased one element at a time.  This function returns 2 if this element a member of the innermost boundary of the buffer and does not need to be adapted.  It returns 1 if this elment needs to be refined and some of its sons will be members of the innermost boundary of the buffer layer 
  int if_first_buffer_boundary(HashTable *ElemTable, double contour_height);

  //! the buffer layer is a layer of refined cells on the outside of the pile, i.e. ((pileheight<contour_height)&&(Influx[0]==0)) and adjacent to the pile.  It is "N" elements wide, and the "N" element width is increased one element at a time.  This function returns 2 if this element a member of the boundary of the buffer that is one element wider than the current buffer and does not need to be adapted.  It returns 1 if this elment needs to be refined and some of its sons will be in the next buffer boundary
  int if_next_buffer_boundary(HashTable *ElemTable, HashTable *NodeTable,
			      double contour_height);
  
  //! for debugging only
  int get_counted();

  //! for debugging only
  void put_counted(int countedvalue);

  //! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
  int get_ithelem();

  //! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
  void put_ithelem(int i);

  //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus the effective bed friction angle holds either the value of the actual bed friction angle if it should not be stopped or the value of the internal friction angle if it should not be stopped
  double get_effect_bedfrict();

  //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, if the effective bed friction angle equals the internal friction angle effect_kactxy takes on the value 1, k active/passive comes from using a Coulomb friction model for granular flows
  double* get_effect_kactxy();


  //! this inline member function returns the stored value of Awet, Awet is the fraction of an element's area that is wet (has material), 0.0<=Awet<=1.0, where there is no flow (pileheight < GEOFLOW_TINY) Awet=0, in the interior of the Flow Awet=1.0, at the boundary of the flow, elements will be PARTIALLY WET (i.e. where the element SHOULD be separated into a dry part and wet part), Awet is the fraction that should be wet, Awet is updated during the corrector part of the (finite difference)predictor-(finite volume)corrector update.  Fluxes are adjusted to acount for how wet/dry an edge of an element is. Keith wrote this may 2007
  double get_Awet();
  
  //! this inline member function assigns a value to Awet, Awet is the fraction of an element's area that is wet (has material), 0.0<=Awet<=1.0, where there is no flow (pileheight < GEOFLOW_TINY) Awet=0, in the interior of the Flow Awet=1.0, at the boundary of the flow, elements will be PARTIALLY WET (i.e. where the element SHOULD be separated into a dry part and wet part), Awet is the fraction that should be wet, Awet is updated during the corrector part of the (finite difference)predictor-(finite volume)corrector update.  Fluxes are adjusted to acount for how wet/dry an edge of an element is. Keith wrote this may 2007
  void   put_Awet(double Awet_in);

  //! this inline member function returns the stored value of Swet.  Swet is the fraction of the element's partially wet sides that are wet (i.e. have material).  Where there is no flow (pileheight < GEOFLOW_TINY), Swet=0.  In the interior of a flow, Swet=1.0.  At the flow boundary, elements will be PARTIALLY WET, 0.0<=Swet<=1.0.  Due to symmetry, any element can only have 0 or 2 partially wet sides, each of which (for normalized elements) will have the same fraction that is wet, Swet.  Swet for each partially wet cell is updated every time-step when calc_wet_dry_orient() is called in step.C.  Fluxes are adjusted to account for how wet/dry an edge of an element is.  Keith wrote this function may 2007
  double get_Swet();

  //! this inline member function assigns a value to Swet.  Swet is the fraction of the element's partially wet sides that are wet (i.e. have material).  Where there is no flow (pileheight < GEOFLOW_TINY), Swet=0.  In the interior of a flow, Swet=1.0.  At the flow boundary, elements will be PARTIALLY WET, 0.0<=Swet<=1.0.  Due to symmetry, any element can only have 0 or 2 partially wet sides, each of which (for normalized elements) will have the same fraction that is wet, Swet.  Swet for each partially wet cell is updated every time-step when calc_wet_dry_orient() is called in step.C.  Fluxes are adjusted to account for how wet/dry an edge of an element is.  Keith wrote this function may 2007
  void   put_Swet(double Swet_in);

  //! this inline member function returns the value of iwetnode.  iwetnode is an integer that defines which of an element's 9 nodes is its "wettest" node (wet elements are those containing material).  In the interior of a flow, iwetnode=8 (the center node), indicating a fully wet element.  Outside of a flow (where an element and all it's neighbors have pileheight < GEOFLOW_TINY), iwetnode is also 8.  Along a flow boundary, partially wet elements with 1,2, or 3 wet sides can have an iwetnode other than 8.   iwetnode is used to determine which side of the dryline in a partially wet element has material.  Keith wrote this function may 2007
  int  get_iwetnode();

  //! this inline member function sets the value of iwetnode. iwetnode is an integer that defines which of an element's 9 nodes is its "wettest" node (wet elements are those containing material).  In the interior of a flow, iwetnode=8 (the center node), indicating a fully wet element.  Outside of a flow (where an element and all it's neighbors have pileheight < GEOFLOW_TINY), iwetnode is also 8.  Along a flow boundary, partially wet elements with 1,2, or 3 wet sides can have an iwetnode other than 8.   iwetnode is used to determine which side of the dryline in a partially wet element has material.  Keith wrote this function may 2007
  void put_iwetnode(int iwetnode_in);

  //! this inline member function returns the array "drypoint".  drypoint[0] is the local x-coordinate, and drypoint[1] the local y-coordinate of its namesake, which is used to specify the position of the flow-front (or dryline) inside a given element.  The position of the dryline along with iwetnode is used to determine Awet, i.e. which fraction of a partially wet element is wet (contains material).  Keith wrote this function may 2007
  double* get_drypoint();

  //! this inline member function sets the values of the array "drypoint".  drypoint[0] is the local x-coordinate, and drypoint[1] the local y-coordinate of its namesake, which is used to specify the position of the flow-front (or dryline) inside a given element.  The position of the dryline along with iwetnode is used to determine Awet, i.e. which fraction of a partially wet element is wet (contains material).  Keith wrote this function may 2007
  void    put_drypoint(double *drypoint_in);

  //! the element member function calc_wet_dry_orient() determines the orientation of the dryline and which side of it is wet, the wet fraction (Swet) of a partially wet edge, the location of the drypoint, it does NOT calculate the wet area (Awet)... these quantities are used in the adjustment of fluxes in partially wet elements. calc_wet_dry_orient() is not coded for generic element orientation, i.e. the positive_x_side must be side 1.  Keith wrote this may 2007
  void calc_wet_dry_orient(HashTable *El_Table);

  //! The Element member function calc_elem_edge_wet_fraction() returns the "how much of this is wet" fraction of side that this element shares with its ineigh-th neighboring element . This fraction is used to determine wether or not to "zero" the state variables used to compute physical fluxes through a "dry" side.  Keith wrote this function may 2007
  double calc_elem_edge_wet_fraction(int ineigh, int ifusewholeside);

  //! this function relaxes the zeroing of fluxes through cell edges that are completely dry at the beginning of the timestep, as indicated by calc_elem_edge_wet_fraction(), but will be at least partly wet by the end of the timestep... Keith wrote this function June 2007
  double calc_elem_edge_wetness_factor(int ineigh, double dt);

  //! The Element member function convect_dryline() calculates the coordinates of the "drypoint" in the element's local coordinate system.  This is used to determine the location of the wet-dry front (or dryline) inside this element, which in turn is used (in conjunction with the location of "iwetnode" - which indicates which side of the dryline is wet) to determine the fraction of its total area that is wet (Awet).  Awet is then returned by the function.  Keith wrote this function may 2007
  double convect_dryline(double VxVy[2], double dt);

  //! sgn of double
  double sgn (double a) { return ( a < 0 ? -1. : 1); }

 private:
  //! myprocess is id of the process(or) that owns this element
  int myprocess;

  //! generation is how many times this element has been refined, currently -8<=generation<=3, a negative generation number means it has been unrefined beyond the orignal coarse mesh, a positive generation number means it has been refined (is smaller than the original element size)
  int generation;

  //! opposite_brother_flag indicate if we have the correct key for the non-neighbor brother (0:= don't have info, 1:= have info)
  int opposite_brother_flag;  

  //! the material flag indicates which material should be used to set this element's bed friction, this is for when a GIS material map, specifying different materials in different spatial regions of the map, the GIS material map is a non standard grass map format that Laercio Namikawa developed, it's stored in the "cats" folder under a grass mapset directory
  int material;/*! ! ! THE MAT. FLAG ! ! !*/

  //! this is the load-balancing weight
  double lb_weight; 

  //! this is the key for load-balancing, if there is no constrained node, it is the element key, otherwise it is a construct of the element "bunch", keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned lb_key[KEYLENGTH];  

  //! this is the element key, which has the same value as the key of the element's bubble node, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned key[KEYLENGTH];

  //! this array holds the first 8 (0->7) of this element's nodes' keys, the n9th (8 out of 0->8) node is the bubble node it's key is not stored separately since it has the same key as the element, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned node_key[8][KEYLENGTH];

  //! this array holds the keys of this element's 8 neighbors (2 neigbors to a side if the neighbor is more refined than this element, otherwise the two neighbor keys for that side are identical in value), having 8 neighbors is an outcome of the 1 irregularity refinement rule, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned neighbor[8][KEYLENGTH];

  //! the key of the father it is assigned in the refine() and unrefine_elements() functions
  unsigned father[KEYLENGTH];

  //! this array holds the keys of this element's 4 sons, it is only used temporarily in the refinement process before the father (this element) is deleted, there's was an old comment "garantee ccw" associated with this variable, I don't know what it means, keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned son[4][KEYLENGTH];

  //! this array holds the process(or) id of this element's 8 neighbors, there can be 8 neighbors because of the 1 irregularity rule.  neigh_proc[4:7] != -2 only if it has 2 neighbors on that side, a value of -1 for neigh_proc means that this edge is a boundary of the computational domain. 
  int neigh_proc[8];

  //! this is legacy afeapi, all finite volume "elements"/cells are piece wise constant, but I believe this is actually used in the DG (Discontinuous Galerkin) version of titan
  int order[5];

  //! neigh_gen is an array that holds the "generation" (how refined it is) of this element's 8 neighbors, there can-be/are 2 neighbors to a side because of the 1 irregularity rule
  int neigh_gen[8];

  //! pointer to the boundary condition class, if this element is not a boundary element the pointer holds the NULL value
  BC* bcptr;

  //! the number of degrees of freedom, since Titan is a finite difference/volume code, ndof is afeapi legacy, but the DG (Discontinuous Galerkin) version of Titan actually uses this
  int ndof;

  //! this is legacy afeapi, it is not used, but do not remove it, it could cause problems if you do
  int no_of_eqns;

  //! this holds the "error" in the element's solution, which is useful in determining refinement, this may actually be afeapi legacy
  double el_error[EQUATIONS];

  //! this holds the element solution, I believe this is legacy afeapi
  double el_solution[EQUATIONS];

  //! refined is a flag that usually has the value 0, but will be 1 if the element has been refined this iteration (used to enforce the 1 irregularity rule), or have the value "GHOST" if it is a ghost cell, refined and ghost cells are not updated, see constant.h for the value of GHOST
  int    refined; 

//! The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell. This allowed Keith to implement one time only immunity to unrefinement for recently refined (NEWSON) elements, which allowed him to protect a refined layer of buffer cells around piles.  Keith has partially replaced refined, get_refined_flag() and put_refined_flag() with adapted, get_adapted_flag() and put_adapted_flag(), but has left the if statements in the code responsible for refinement and unrefinement untouched because he encountered a bug, that he narrowed to within H_adapt() or a function called from within H_adapt(), recurse down, but has not pinpointed.  Keith believes the bug is related to the refined flag being inappropriately set to 1, or not unset to zero when it should be.
  int    adapted;

  //! which_son holds the value of which son this element is, which of the 4 squares that makes up the father elements square.
  int    which_son;

  //! the new_old flag is used in mesh adaptation and repartitioning
  int    new_old;

  //! this array holds the keys of this element's 4 brothers (an element is considered to be it's own brother), this information is used during mesh unrefinement (combining the 4 brothers to make their father), keys are used to access elements or nodes through the appropriate hashtables, each key is a single number that fills 2 unsigned variables
  unsigned brothers[4][KEYLENGTH];

  //! coord holds the coordinates of the elements cell center, these are the same as the coordinates of the element's bubble node's
  double coord[DIMENSION];

  //! elm_loc is used in unrefining beyond the original coarse mesh
  int elm_loc[2];

  //! this is afeapi legacy
  int send[8];

  //! this is afeapi legacy
  int recv[8];

  /* variables for hyperbolic geoflow problem */

  //! state_vars is an array that holds the current state variables: h, hVx, and hVy 
  double state_vars[NUM_STATE_VARS]; 
  //! these are the values of the state variables from before the predictor step
  double prev_state_vars[NUM_STATE_VARS]; 

  //! these are the spatial (x and y) derivatives of the state variables: (dh/dx, dhVx/dx, dhVy/dx, dh/dy, dhVx/dy, dhVy/dy)
  double d_state_vars[NUM_STATE_VARS*DIMENSION]; 

  //! the short speed is the speed computed as: shortspeed=|v|=|dhv/dh|=|v*dh/dh+h*dv/dh|=|v+h*dv/dh| which goes to |v| in the limit of h->0, this is a more accurate way to compute speed when the pile in this cell is short, hence the name "shortspeed" but it is not accurate when the pile is tall, that is when h*dv/dh is large, this is the value from the previous iteration (so there is lagging when using the shortspeed, but this should still be much more accurate than hV/h when h->0. Keith implemented this in late summer 2006, 
  double shortspeed; 

  //! length of the element in the global x and y directions: dx and dy 
  double dx[DIMENSION]; 

  //! for structured grid, tells which side is the positive x direction
  int    positive_x_side; 

  //! maximum x and y direction wavespeeds for this element, wavespeeds are eigenvalues of the flux jacobians
  double eigenvxymax[DIMENSION];

  //! this is legacy afeapi not used in titan
  double coefABCD[4];

  //! k active/passive in the x and y directions, k active/passive is part of the coulomb friction model for Granular Flows
  double kactxy[DIMENSION];

  //! terrain elevation at this elements center/bubble node 
  double elevation; 

  //! terrain slope in the global x and y directions
  double zeta[DIMENSION];

  //! Curvature itself is the inverse of radius of curvature.  The exact value of curvature is the spatial second derivative of the normal coordinate of the surface along directions tangent to the surface at that point (local x and y).  However I believe that Laercio Namikawa implemented it approximately, i.e. as the global x and y second derivatives of terrain elevation. 
  double curvature[DIMENSION]; 

  //! the gravity vector in local x,y,z coordinates (z is normal to the terrain surface, the projections of the x and y local directions onto a horizontal plane are aligned with the global x and y directions)
  double gravity[3];

  //! the spatial (x and y) derivatives of the local z component of the gravity vector
  double d_gravity[DIMENSION]; 

  //! legacy titan not really used, lam :=p_{bed}/(rho*g_z*h)
  double lam; 

  //! part of the new stopping criteria under development, has a 0 if flow is not stopped, has the value 1 if it should not be sliding but should be slumping, has the value 2 if it should neither be sliding or slumping (it should be completely stopped), I (Keith) am rather confident in the criteria used to set this the problem is determining what to do about it after you know the flow SHOULD be stopped
  int stoppedflags;

  //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus the effective bed friction angle holds either the value of the actual bed friction angle if it should not be stopped or the value of the internal friction angle if it should not be stopped
  double effect_bedfrict;

  //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, thus effect_tanbedfrict holds the value of the effective bed friction angle
  double effect_tanbedfrict;

  //! one option for what to do when you know the flow should be stopped is to reset the bed friction angle to take on the value of the internal friction angle, if the effective bed friction angle equals the internal friction angle effect_kactxy takes on the value 1, k active/passive comes from using a Coulomb friction model for granular flows
  double effect_kactxy[2];

  //! extrusion flux rate for this timestep for this element, used when having material flow out of the ground, a volume per unit area influx rate source term
  double Influx[NUM_STATE_VARS];

  int counted;
  
  //! when sorted by keys this element is the ithelem element on this processor, ithelem is just storage for a value you have to assign before using, if you do not compute it before you use it will be wrong.
  int ithelem;

  //! the node number {0,1,..,7} of this element's "most wet node", this dictates both the orientation of the "dryline" and which side of it is wet.  The "dryline" is the line that divides a partially wetted element into a dry part and a wet part, for the sake of simplicity only 4 orientations are allowed, horizontal, vertical, parallel to either diagonal of the element.  If the iwetnode is an edge node of this element then the dryline is parallel to the edge the element is on, if the iwetnode is a corner node of this element then dryline is parallel to the diagonal of the element that the iwetnode is not on.  Which side of the dryline is wet is the same side in which the iwetnode resides (and is determined each timestep based soley on which of the elements neighbors currently have pile height greater than GEOFLOW_TINY)... as such iwetnode can be thought of as the MOST WET NODE  of the element.  Having iwetnode==8 indicates that the element is uniformly wet if this element's pile height is greater than GEOFLOW_TINY or is uniformly dry if pileheight is less than or equal to GEOFLOW_TINY.
  int iwetnode;

  //! Awet is the ratio of this element's wet area to total area (always between 0 and 1 inclusive) when taken together with iwetnode, this uniquely determines the exact placement of the "dryline" within the current element.  Awet is initially set by source placement to be either 0 (no material) or 1 (material) and is updated by the corrector part of the predictor-corrector method, the new value is determined by where the dry line has been convected to over this timestep. Keith wrote this May 2007.
  double Awet;
  
  //! center point of the "dryline", x and y coordinates value ranges between -0.5 and 0.5 with 0 being the center of the element, since the wet/dry interface is taken to be a non-deforming non rotating (within the timestep) "dryline" convecting a single point on the dryline (called the drypoint) is sufficient to determine the new placement of the dryline which allows us to update Awet... Keith wrote this May 2007.
  double drypoint[2];

  //! when an element edge is partially wet and partially dry... Swet is the fraction of a cell edge that is partially wet, because it can only be horizontal, vertical, or parallel to either diagonal, all of one element's partially wet sides are have the same fraction of wetness.  The state variables (used to compute the physical fluxes) at the element/cell edge are adjusted to be the weighted by wetness average over an element/cell edge.  As such physical fluxes through completely dry edges of partially wet elements/cells are zeroed, while physical fluxes through completely wet edges are left unchanged.  Because of the definition as "wetness weighted average" physical fluxes through a partially wet edge shared with a neighbor of the same generation is also left left unchanged but, when a partially wet edge is shared with two more refined neighbors the total mass and momentum at the edge is split between the two neighbors in proportion to how much of their boundary shared with this element is wet.  This "scaling" of the physical fluxes is the "adjustment of fluxes in partially wetted cells" facet of our multifaceted thin-layer problem mitigation approach.  And it has been shown to significantly reduce the area covered by a thin layer of material.  Keith wrote this May 2007.
  double Swet;
};

inline int Element::get_ithelem() {return ithelem;};

inline void Element::put_ithelem(int i) {ithelem=i;};

inline unsigned* Element::pass_key(){return key;};

inline int Element::get_material(){return material;};

inline unsigned* Element::get_brothers() {return &brothers[0][0];};

inline void Element::void_bcptr() {bcptr = NULL;};

inline double Element::get_lb_weight() {return lb_weight;}; 

inline void Element::put_lb_weight(double dd_in) {lb_weight = dd_in;}; 

inline unsigned* Element::get_lb_key() {return lb_key;}; 

inline void Element::put_myprocess(int in_proc) {myprocess = in_proc;}; 

inline int Element::get_myprocess() {return myprocess;}; 

inline int Element::get_opposite_brother_flag() {return opposite_brother_flag;};
  
inline void Element::put_height_mom(double pile_height,double volf,double xmom,double ymom) 
{
  prev_state_vars[0]=state_vars[0]=pile_height; 
  prev_state_vars[1]=state_vars[1]=pile_height*volf; 
  prev_state_vars[2]=state_vars[2]=xmom; 
  prev_state_vars[3]=state_vars[3]=ymom; 
  if(pile_height>GEOFLOW_TINY){
    shortspeed=sqrt(xmom*xmom+ymom*ymom)/(pile_height*volf);
    Awet=1.0;
  }
  else
  {
    shortspeed=0.0;
    Awet=0.0;
  }  
  return;
}

inline void Element::put_height(double pile_height) 
{
  put_height_mom(pile_height, 1., 0., 0.); 
  return;
}

inline double* Element::get_state_vars() {return state_vars;};

inline double* Element::get_d_state_vars() {return d_state_vars;};

inline double* Element::get_dx() {return dx;};

inline int Element::get_positive_x_side() {return positive_x_side;};

inline double* Element::get_prev_state_vars() {return prev_state_vars;};

inline void Element::update_prev_state_vars()
{
  for (int i=0; i<NUM_STATE_VARS; i++)
    prev_state_vars[i]=state_vars[i];
}

inline double* Element::get_eigenvxymax() {return eigenvxymax;};

inline double Element::get_shortspeed() {return shortspeed;};

inline void Element::put_shortspeed(double shortspeedin) {shortspeed=shortspeedin; return;};

inline double* Element::get_kactxy() {return kactxy;};  

inline double* Element::get_gravity() {return gravity;};

inline double Element::get_elevation() {return elevation;};  

inline double* Element::get_d_gravity() {return d_gravity;};  

inline double* Element::get_curvature() {return curvature;};

inline void Element::put_lam(double lam_in) {lam = lam_in;};

inline double Element::get_lam() {return lam;}; 

inline int* Element::get_elm_loc() {return elm_loc;};  

inline void Element::put_elm_loc(int* int_in) {elm_loc[0] = int_in[0]; elm_loc[1] = int_in[1];};  

inline double* Element::get_coord() {return coord;};  

inline void Element::put_stoppedflags(int stoppedflagsin) {stoppedflags=stoppedflagsin;};
  
inline int Element::get_stoppedflags() {return stoppedflags;};


//above this line Keith made inline 20061128
/* REALLY? member functions defined in class body are 
 * automatically inlined by the complier. When they are
 * not, "inline" keyword is not going to help the cause.
 * This was a huge waste of time
 */

inline int* Element::getassoc(){ return neigh_proc;}

inline unsigned* Element:: getNode(){ return &(node_key[0][0]);}

inline int Element::get_no_of_dof(){return ndof;}

inline void Element::put_gen(int g){generation = g;}

inline void Element::put_neigh_proc(int i, int proc) { neigh_proc[i] = proc;}

inline void Element::put_order(int i, int ord) {order[i] = ord;}

inline unsigned* Element::getson(){return &(son[0][0]);}

inline int* Element::get_order(){return order;}

inline void Element::putson(unsigned* s)
{

  for(int i=0; i<4; i++)
    for(int j=0; j<KEYLENGTH;j++)  
      son[i][j] = *(s+i*KEYLENGTH+j);

  refined=1;
  adapted=OLDFATHER;
}

inline void Element::putbrothers(unsigned* s)
{

  for(int i=0; i<4; i++)
    for(int j=0; j<KEYLENGTH;j++)  
      brothers[i][j] = *(s+i*KEYLENGTH+j);

  return;
}

inline void Element::putneighbor(unsigned* n, int i)
{
  int j;
  for(j=0; j<KEYLENGTH;j++)  neighbor[i][j] = *(n+j);
}

inline void Element::putassoc(int a, int i){neigh_proc[i] = a;}

inline void Element::putel_sq(double solsq, double errsq)
{  
  el_solution[0] = solsq;
  el_error[0] = errsq;
}

inline double* Element::get_el_solution()
{return el_solution;}

inline double* Element::get_el_error()
{return el_error;}

inline unsigned* Element::get_neighbors()
{return &neighbor[0][0];}

inline int* Element::get_neigh_proc()
{return neigh_proc;}

inline BC* Element::get_bcptr()
{return bcptr;}

inline int Element::get_gen()
{return generation;}

inline int Element::get_refined_flag()
{return refined;} 

inline void Element::put_refined_flag(int i)
{refined = i;}

inline int Element::get_adapted_flag() {return adapted;} 

inline void Element::put_adapted_flag(int new_adapted_status) {adapted = new_adapted_status;}



inline void Element::put_send_flag(int i, int j) {send[i] = j;}

inline int Element::get_send_flag(int i) {return send[i];}

inline void Element::put_recv_flag(int i, int j) {recv[i] = j;}

inline int Element::get_recv_flag(int i) {return recv[i];}

inline int* Element::get_neigh_gen() {return neigh_gen;}

inline void Element::put_neigh_gen(int i, int gen) {neigh_gen[i] = gen;}

inline void Element::put_which_son(int i){ which_son = i;}

inline int  Element::get_which_son(){return which_son;}

inline int  Element::get_new_old(){return new_old;}

inline void Element::put_new_old(int i) {new_old = i;}

inline void Element::change_neighbor_process(int which, int newp){neigh_proc[which]=newp;}

inline void Element::zero_influx() 
{
  for (int i=0; i<NUM_STATE_VARS; i++)
    Influx[i]=0.;
}

inline double* Element::get_influx() {return Influx;};

inline void Element::put_father(unsigned *fatherin){
  for(int ikey=0;ikey<KEYLENGTH;ikey++) father[ikey]=fatherin[ikey];};

inline int Element::get_counted() {return counted;};

inline void Element::put_counted(int countedvalue) {counted=countedvalue;};

inline double  Element::get_effect_bedfrict(){return effect_bedfrict;};
inline double* Element::get_effect_kactxy(){return effect_kactxy;};

inline double  Element::get_Awet(){return Awet;};
inline void    Element::put_Awet(double Awet_in){Awet=Awet_in; return;};
inline double  Element::get_Swet(){return Swet;};
inline void    Element::put_Swet(double Swet_in){Swet=Swet_in; return;};
inline int     Element::get_iwetnode(){return iwetnode;};
inline void    Element::put_iwetnode(int iwetnode_in){iwetnode=iwetnode_in; return;};
inline double* Element::get_drypoint(){return drypoint;};
inline void    Element::put_drypoint(double *drypoint_in){drypoint[0]=drypoint_in[0];drypoint[1]=drypoint_in[1]; return;};


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

//! The ElemPtrList class is basically just a "smart array" of pointers to Elements, by smart I mean it keeps track of its size and number of Elements in the list and expands/reallocates itself whenever you add an element ptr to the list when you've run out of space, it also keeps a record of the index of the first "new" element pointer you've added in the current series, which is useful for the intended purpose... ElemList was designed for use in refinement and unrefinement to replace fixed sized arrays (length=297200) of pointers to Elements.  The reason for this upgrade was it was causing valgrind to issue all kinds of warnings about the "client switching stacks" and "invalid write/read of size blah blah blah" because the stacksize was too large.  My 20061121 rewrite of hadapt.C and unrefine.C to make them "fast" caused this problem because I added a second (large) fixed sized array to both of them so I could reduce the number of hashtable scans by only revisiting the "new" additions to the array of pointers of Elements. --Keith wrote this on 20061124, i.e. the day after Thanksgiving, and I'm very thankful for having the inspiration to figure out the cause of valgrid warning
class ElemPtrList{
 public:

  /*
  ElemPtrList();
  ElemPtrList(int initial_size);
  ~ElemPtrList();
  void add(Element* EmTemp);
  */

  //! this constructor allocates space for an array of the default initial size (1024 Element pointers), the size_increment equal the initial size
  ElemPtrList(){
    init(1024);
    return;
  };

  //! this constructor allocates space for user specified initial-size, the size_increment equals the initial size.
  ElemPtrList(int initial_size){  
    if(initial_size==0) initial_size=1024;
    init(initial_size);
    return;
  };

  //! the destructor frees the list space so the programmer never has to worry about it
  ~ElemPtrList(){
    //printf("list_space=%d, num_elem=%d, inewstart=%d\n",list_space,num_elem,inewstart);
    free(list);
    return;
  }

  //! add an element pointer to the list, it will increase the size of the list by the size_increment if necessary, the size_increment is the initial size.
  void add(Element* EmTemp);

  //! returns the ith Element pointer stored in the list
  Element*  get(int i);

  //! returns the key of the ith Element whose pointer is stored in the list
  unsigned* get_key(int i);

  //! returns the number of elements in the list.
  int       get_num_elem();

  //! marks the "starting position" of a new "series" of Element pointers stored
  void      set_inewstart(int inewstart_in);
  
  //! returns the "starting position" of the new "series" of Element pointers stored in the list
  int       get_inewstart();

  //! zeros the list (does not deallocate space)
  void      trashlist();

 private:
  //! actually creates the list, is called by the constructors
  //void      init(int initial_size);
  void      init(int initial_size){
    list_space=size_increment=initial_size;    
    num_elem=inewstart=0;
    list=(Element **) malloc(list_space*sizeof(Element*));
    for(int i=0;i<list_space;i++) list[i]=NULL;
  };
  
  //! number of elements whose pointers are stored in the list
  int       num_elem;
  
  //! when the list runs out of space increase it by this much
  int       size_increment;

  //! the current size of the list (the ammount of memory allocated not number of nonzero entries)
  int       list_space;

  //! a convienient way to mark the start of a new series in the list
  int       inewstart;
  
  //! the "array" holding the list of element pointers, space is allocated by the constructor, increased automatically whenever needed, and freed automatically by the destructor
  Element** list;
};

inline Element*  ElemPtrList::get(int i){return (((i>=0)&&(i<num_elem))?list[i]:NULL);};
inline unsigned* ElemPtrList::get_key(int i){return (((i>=0)&&(i<num_elem))?list[i]->pass_key():NULL);};
inline int       ElemPtrList::get_inewstart(){return inewstart;};
inline void      ElemPtrList::set_inewstart(int inewstart_in){inewstart=inewstart_in; return;};
inline int       ElemPtrList::get_num_elem(){return num_elem;};

inline void      ElemPtrList::trashlist(){
  for(int i=0;i<num_elem;i++) list[i]=NULL; 
  num_elem=inewstart=0; 
  return;
};

inline void ElemPtrList::add(Element* EmTemp){
  if(num_elem==list_space-1){
    list_space+=size_increment;
    list=(Element **) realloc(list,list_space*sizeof(Element *));
  }

  list[num_elem]=EmTemp;
  num_elem++;
  return;
};
  


#endif


