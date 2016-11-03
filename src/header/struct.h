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
 * $Id: struct.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef STRUCT_H
#define STRUCT_H

//! ElemPack is a smaller (memory spacewise) version of Element that can be sent from one processor to another via MPI calls
struct ElemPack{
  //see ../repartition/new_datatype.C blockcounts[3]={58,25*KEYLENGTH,85}
  int        myprocess;                                              //  1
  int        generation;                                             //  2
  int        material;/*flag added by andrew*/                       //  3
  int        neigh_proc[8];                                          // 11
  int        order[5];                                               // 16
  int        neigh_gen[8];                                           // 24
  int        ndof;                                                   // 25
  int        no_of_eqns;                                             // 26
  int        refined;                                                // 27
  int        adapted;                                                // 28
  int        which_son;                                              // 29
  int        new_old;                                                // 30
  int        n_order[9];                                             // 39
  int        n_info[9];                                              // 48
  int        bc; /*flag indicating if there is or there is no bc*/   // 49
  int        bc_type[4];                                             // 53
  int        positive_x_side;                                        // 54
  int        elm_loc[2];                                             // 56
  int        opposite_brother_flag;                                  // 57
  int        iwetnode;                                               // 58

  unsigned   key[KEYLENGTH];/*contains the 9th node key*/   //  1
  unsigned   node_key[8][KEYLENGTH];                        //  9
  unsigned   neighbor[8][KEYLENGTH];                        // 17
  unsigned   son[4][KEYLENGTH];                             // 21
  unsigned   brothers[4][KEYLENGTH];                        // 25 * KEYLENGTH


  /*NODE DEF*/
  double elevation;                                   //  1
  double n_coord[9][2];                               // 19 
  double el_error[EQUATIONS];                         // 21
  double el_solution[EQUATIONS];                      // 23  
  double bc_value[4][2][2];                           // 39
  double state_vars[NUM_STATE_VARS];                  // 42
  double prev_state_vars[NUM_STATE_VARS];             // 45 
  double d_state_vars[NUM_STATE_VARS*DIMENSION];      // 51
  double shortspeed;                                  // 52
  double dx[DIMENSION];                               // 54 
  double eigenvxymax[DIMENSION];                      // 56
  double kactxy[DIMENSION];                           // 58
  double zeta[DIMENSION];                             // 60
  double curvature[DIMENSION];                        // 62
  double gravity[3];                                  // 65
  double d_gravity[DIMENSION];                        // 67
  double lam;                                         // 68
  double lb_weight;                                   // 69
  double node_elevation[9];                           // 78
  double Influx[NUM_STATE_VARS];                      // 81
  double Awet;                                        // 82
  double Swet;                                        // 83
  double drypoint[DIMENSION];                         // 85

  /*BC DEF*/



};

//                             \|||/    
//                             (o o)   
//---------Elementlink------oo0-(_)-0oo----------STARTS HERE-------


struct ElementLink{
  int             target_proc;
  int             new_proc;
  unsigned        elkey[KEYLENGTH]; 
  unsigned        targetkey[KEYLENGTH];

  ElementLink*    pre;
  ElementLink*    next;
  
  ElementLink(unsigned* keyi,unsigned* key2, int tp, int np)
    {
      target_proc = tp;
      new_proc    = np;
      int i;
      for(i=0;i<KEYLENGTH; i++)
	{
	  elkey[i] = *(keyi+i);
	  targetkey[i]= *(key2+i);
	}
      next   = NULL;
      pre    = NULL;
    }
  ElementLink(){
    next   = NULL;
    pre    = NULL;
  }
  ~ElementLink(){
    if(next)
      next->pre = pre;
    if(pre)
      pre->next = next;
  }


};

typedef ElementLink* ELinkPtr;


//                              \|||/    
//                              (o o)   
//--------MPI datatype-------oo0-(_)-0oo----------STARTS HERE-------


struct NeighborPack{
 
  int             target_proc;
  int             new_proc;
  unsigned        elkey[KEYLENGTH]; 
  unsigned        targetkey[KEYLENGTH];


};

struct Neigh_Sol_Pack{


int      nside;
int      norder[5];
unsigned key[KEYLENGTH];
double   solu[2][121];
double   Xnod[18];

};

typedef NeighborPack* NePtr;

#endif
