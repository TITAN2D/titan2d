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

#ifdef CRAY
#include <fortran.h>
#endif 

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include "properties.h"
#include "node.h"
#include "element2.h"
#include "hashtab.h"
#include "recv.h"
#include "extfun.h"
#include "geoflow.h"
#include "scale.h"
#include "../gisapi/GisApi.h"
#include "../useful/useful_lib.h"
#include "FileFormat.h"
#include "flux_srcs.h"


#undef CRAY

