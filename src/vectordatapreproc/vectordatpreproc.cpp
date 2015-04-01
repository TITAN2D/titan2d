/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Mon Feb 21 17:41:53 EST 2005
    copyright            : (C) |YEAR| by Amrita Chanda
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
/*
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
*/

#include <iostream>
#include <stdlib.h>
using namespace std;

#include "cgis_vectordata.h"
#include "../gisapi/GisApi.h"

int main(int argc, char *argv[])
{
  
  if(argc < 6){
    cout<<"ERROR: Usage ./vectordatapreprocessor <GIS Dbase> <GIS Maplocation> <GIS Mapset> <GIS Rasterimage> <GIS Vectorfile>\n";
    return EXIT_FAILURE;
  }
  

  
  //Reads in the DEM elevations and stores them in a 2D array
  //CGIS_TriGrid* _triGrid = new CGIS_TriGrid(argv[1], argv[2], argv[3], argv[4]);
  //  CGIS_TriGrid* _triGrid = new CGIS_TriGrid("/home/achanda/Research/CheckedOutFromFire/example/grass.data/grass5", "Colima",
  //					    "ColimaR", "ColimaR");

  //Get res, minx, maxx, miny, maxy, numofRows, numOfCols from GIS
  double res, xmin, xmax, ymin, ymax;
  int rows, cols;

  int err = Initialize_GIS_data(argv[1], argv[2], argv[3], argv[4]);
  //int err = Initialize_GIS_data("/home/achanda/Research/CheckedOutFromFire/example/grass.data/grass5", "Colima",
  //				"ColimaR", "ColimaR");
  if(err)
    {
      
      cout<<"ERROR: Unable to initialize GIS data file\n";
      return EXIT_FAILURE; 

    }
  
  err = Get_max_resolution(&res);
  if(err){
   
    cout<<"\n\nERROR: Unable to get max resolution from GIS";
    return EXIT_FAILURE; 

  }

  err = Get_window(&xmin, &xmax, &ymin, &ymax);
  if(err){
   
    cout<<"\n\nERROR: Unable to get window extents from GIS";
    return EXIT_FAILURE; 

  }

  err = Get_number_of_rows(&rows);
  if(err){
    
    cout<<"\n\nERROR: Unable to get number of rows from GIS";
    return EXIT_FAILURE; 
    
  }

  err = Get_number_of_columns(&cols);
  if(err){
    
    cout<<"\n\nERROR: Unable to get number of columns from GIS";
    return EXIT_FAILURE; 
    
  }

  if(Delete_GIS_data()){

    cout<<"\n\nERROR: Unable to delete GIS data";
    return EXIT_FAILURE; 

  }

  //Reads in the XYZ values of points on all lines and creates a CPolyLine obj for each line
  CGIS_VectorData* _vecData = new  CGIS_VectorData(res, xmin, xmax, ymin, ymax,  rows, cols, argv[1], argv[2], argv[3], argv[5]);

  //  CGIS_VectorData* _vecData = new  CGIS_VectorData(res, xmin, xmax, ymin, ymax,  rows, cols, 
  //						   "/home/achanda/Research/CheckedOutFromFire/example/grass.data/grass5", "Colima",
  //						   "ColimaR", "Flow1991_lin.spr");

  if(_vecData->ReadVectorData())
    { 
      //reads XY values of all lines in vectorfile
      cout<<"ERROR: Unable to read vector data file\n";
      return EXIT_FAILURE;
    }

  // _vecData->ProjectOn2DGrid(minx, maxx, miny, maxy, numOfRows, numOfCols);

  if(_vecData->StoreInFile())
    {
      //stores XY values of all 2D projected lines in vector data file
      cout<<"ERROR: Unable to save vector data file\n";
      return EXIT_FAILURE;
    }  
  /*
  CGIS_VectorData* _vecData2 = new  CGIS_VectorData();

  if(_vecData2->ReadFromFile())
    {
      //Read XY values from vectordata file
      cout<<"ERROR: Unable to read vector data file\n";
      return EXIT_FAILURE;
    }
  */
      
  return EXIT_SUCCESS;
}
