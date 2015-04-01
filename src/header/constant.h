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
 * $Id: constant.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef CONSTANTS
#define CONSTANTS
const int  KEYLENGTH = 2;

//const   NODE_TABLE_SIZE=6000;
//const   EL_TABLE_SIZE=2000;


const int  MAX_PROCS = 2056;

const int  ZERO      = '0';
const int  NONZERO   = 111;
//const   NBUCKETS  = 2000;        //every hashtable posseses 5000 entries.
const int   PRIME     = 2017;        //used for creating hash key
const int   DIMENSION = 2;
const int   EQUATIONS = 2;         /* be careful EQUATIONS and NUM_STATE_VARS
				      are completely different things and are 
				      not interchangeable, EQUATIONS appears 
				      to not be used except for the el_error
				      and el_solution variables (members of 
				      element) I believe they are only for 
				      FEM which we aren't really using in the
				      Finite difference/volume version of 
				      titan.  Just leave it at 2 and otherwise
				      ignore it. */
const int   POWER     = 1; 
const int   MAX_ORDER = 6;

const double PI   = 3.1415926;

const int   BCTYPE    = 3;

const float  C = 1.5;

const int   NODEINIT    = 0x0000;
const int   CORNER      = 0x0002;//--corner dof
const int   BUBBLE      = 0x0006;
const int   SIDE        = 0x0004;//--side dof
const int   CONSTRAINED = 0x0001;
const int   S_C_CON     = 0x0007;//--side dof
const int   S_S_CON     = 0x0005;//--no dof
const int   ASSIGNED    = 1;
const int   UNASSIGNED  = 0;

const double   UN_CONSTRAINED =-999999.0;/*for the unconstrained dof*/

const int   ON = 1;
const int   OFF = 0;

const int   NEW = 1;
const int   OLD = 0;

const int   INIT = -1;

const float  MAX_X = 1.0;
const float  MAX_Y = 1.0;
const float  MIN_X = -1.0;
const float  MIN_Y = -1.0;

const float  SMALL = 1.0e-5;
const float  XBC = -1.0;
const float  YBC = -1.0;
const float  LOAD_BALANCE_TOLERANCE = 1.0001;

extern void fhsfc3d(double*, unsigned*, unsigned*);
extern void fhsfc2d_(double*, unsigned*, unsigned*);

/* geoflow data */
const int NUM_STATE_VARS = 3;
#define GEOFLOW_TINY 0.0001
#define GEOFLOW_SHORT 0.01

const int GHOST         = -9876; //"refined" GHOST CELL FLAG

//The magnitude of the "adapted" flag indicates whether the cell is NEWSON, NEWFATHER, NOTRECADAPTED, or TOBEDELETED.  A postive value indicates it's on this processor, a negative sign indicates a GHOST cell.

#define NEWBUFFER      5  //NEWBUFFER elements have at least one current buffer element as a neighbor, a temporary marking needed for building the buffer layer, to be remarked as BUFFER elements
#define BUFFER         4  //Do not refine or unrefine elements in the buffer layer, has to be a separate flag from NEWSON because of triggered refinement
#define NEWSON         3  //Do not refine or unrefine new son elements
#define NEWFATHER      2  //Do not unrefine new father elements 
#define NOTRECADAPTED  1  //you can do anything you want to this element
#define TOBEDELETED    0  //This element has been marked for deletion, do ot involve
#define OLDFATHER     -6  //This is a temporary marking that says neighbors need to be updated with NEWSON information and then this element should be remarked as TOBEDELETED
#define OLDSON        -7  //The plan is make this analogous to OLDFATHER to allow unrefinement NEXT TO but not across interprocessor boundaries, this has not been implemented yet (date of this comment: 20061026)


#endif
