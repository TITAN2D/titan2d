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
%begin %{
#define __STDC_FORMAT_MACROS
%}
%{
#include "../header/constant.h"
#include "../header/properties.h"
#include "../header/titan_simulation.h"
#include "../preproc/preproc.h"
#include "../vectordatapreproc/vectordatpreproc.h"
#include "../header/integrators.h"
#include "../header/stats.hpp"
%}

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectord) vector<double>;
   %template(vectors) vector<string>;
};
%include "../header/constant.h"
%include "../header/properties.h"
%include "../header/integrators.h"
%include "../header/titan_simulation.h"
%include "../preproc/preproc.h"
%include "../vectordatapreproc/vectordatpreproc.h"
%include "../header/stats.hpp"

