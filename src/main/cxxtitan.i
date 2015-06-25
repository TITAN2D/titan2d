/*******************************************************************
 * Copyright (C) 2003-2015 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: NA Simakov
 * Description: python api
 *
 *******************************************************************
 */

%module cxxtitan

%{
#include "../header/titan_simulation.h"
#include "../preproc/preproc.h"
#include "../vectordatapreproc/vectordatpreproc.h"
%}

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectord) vector<double>;
   %template(vectors) vector<string>;
   %template(vector_cxxTitanPile) vector<cxxTitanPile>;
   %template(vector_cxxTitanFluxSource) vector<cxxTitanFluxSource>;
   %template(vector_cxxTitanDischargePlane) vector<cxxTitanDischargePlane>;
};

%include "../header/titan_simulation.h"
%include "../preproc/preproc.h"
%include "../vectordatapreproc/vectordatpreproc.h"

