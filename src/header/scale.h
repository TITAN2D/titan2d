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
 * $Id: scale.h 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifndef SCALE_H
#define SCALE_H
struct ScaleValues {
  double lengthscale;
  double heightscale;
  double gravityscale;

  ScaleValues() {
    ifstream inDatafile("scale.data", ios::in);
    if(inDatafile.fail())
      {
	// assume no scaling then...
	lengthscale = 1;
	heightscale = 1;
	gravityscale = 1;
      }
    inDatafile>>lengthscale;
    if(lengthscale < GEOFLOW_TINY)
      lengthscale = 1;
    inDatafile>>heightscale;
    if(heightscale < GEOFLOW_TINY)
      heightscale = 1;
    inDatafile>>heightscale;
    if(gravityscale < GEOFLOW_TINY)
      gravityscale = 1;
    
    inDatafile.close();
  }
};

#endif
