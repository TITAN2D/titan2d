/***************************************************************************
                          cgis_trigrid.cpp  -  description
                             -------------------
    begin                : Mon Feb 21 2005
    copyright            : (C) 2005 by Amrita Chanda
    email                : achanda@buffalo.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "cgis_trigrid.h"
#include "../gisapi/GisApi.h"

#include <iostream>
#include <stdlib.h>
using namespace std;

CGIS_TriGrid::CGIS_TriGrid(char *GIS_Dbase, char *GIS_Maplocation, char *GIS_Mapset, char *GIS_Rasterimage){

  int err = Initialize_GIS_data(GIS_Dbase, GIS_Maplocation, GIS_Mapset, GIS_Rasterimage);

  if(err) {
    cout<<"\n\nERROR: Unable to initialize GIS for "<<GIS_Dbase<<" : "<<GIS_Maplocation<<" : "<<GIS_Mapset<<" : "<<GIS_Rasterimage<<endl;
    exit(1);
  }

  err = Get_max_resolution(&_res);
  if(err){
    cout<<"\n\nERROR: Unable to get max resolution from GIS";
    exit(1);
  }

  err = Get_window(&_xmin, &_xmax, &_ymin, &_ymax);
  if(err){
    cout<<"\n\nERROR: Unable to get window extents from GIS";
    exit(1);
  }
  
  err = Get_number_of_rows(&_num_of_rows);
  if(err){
    cout<<"\n\nERROR: Unable to get number of rows from GIS";
    exit(1);
  }

  err = Get_number_of_columns(&_num_of_cols);
  if(err){
    cout<<"\n\nERROR: Unable to get number of cols from GIS";
    exit(1);
  }

  int number_of_locations =  _num_of_rows * _num_of_cols;
  
  _elevations = (double *)malloc( number_of_locations * sizeof(double));
  
  for(int i=0; i<_num_of_rows; ++i){
    for(int j=0; j<_num_of_cols; ++j){
      //reading elevations row by row (x increments first before y)
      if(Get_elevation(_res, _xmin + j*_res, _ymin + i*_res, &_elevations[i*_num_of_cols + j])){
      	cout<<"\n\nERROR: Unable to retrieve elevations from GIS";
      	exit(1);
      }  
    }//for(j)
  }//for(i)

  if(Delete_GIS_data()){
    cout<<"\n\nERROR: Unable to delete GIS data";
    exit(1);
  }

}


CGIS_TriGrid::~CGIS_TriGrid(){
  delete _elevations;  
}
