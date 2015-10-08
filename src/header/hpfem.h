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
 * $Id: hpfem.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

//modifed 2003/08/29 11:29:35 by kdalbey
#define  SUNOS  //definition for gmake architecture
#define TOPO_DATA

/*#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR*/
#include <mpi.h>

#ifdef CRAY
#include <fortran.h>
#endif 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "titan2d.h"

#include "properties.h"
#include "elenode.hpp"
#include "recv.h"
#include "extfun.h"
#include "geoflow.h"
#include "scale.h"
#include "../gisapi/GisApi.h"
#include "../useful/useful_lib.h"
#include "FileFormat.h"
#include "flux_srcs.h"

#include "titan2d_utils.h"

//! construct_el is a friend function of the Element class that fills an element with information it receives in a variable of the ElemPack class from an MPI call
void construct_el(Element* newelement, ElemPack* elem2, NodeHashTable* HT_Node_Ptr, int myid, double* e_error);

#undef CRAY

