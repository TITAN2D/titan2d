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
 * $Id: boundary.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef BC_H
#define BC_H


/*
type==1 essential (dof):

value=UN_CONSTRINED that dof is not constriend
value=number the dof is constrined to number

type==2 natural (force)
value[0]=x-component of the load vector
value[1]=y-component of the load vector

type==3 both

value [4]       [2]             [2]
      side   0-natural        0-x dir.
            1-essential       1-y dir.

*/


//! the BC structure contains members: "type[4]" identifying type of boundary condition as essential=1, natural=2, or both=3; and "value[4][2][2]": identifying the element-side, type (0=natural, 1=essential), and load component (0=x, 1=y) comprising the boundary condition.
struct BC
{
  int type[4];
  float value[4][2][2]; 

  //! this function is a constuctor initializes the type and value of the boundary condition to zero
  BC()
    {
      type[0] = 0;
      type[1] = 0;
      type[2] = 0;
      type[3] = 0;
      for(int i=0;i<4;i++)
	for(int j = 0; j<2; j++) 
	  for(int k=0; k<2; k++)
	    value[i][j][k] = 0.0;
    }
};
#endif
