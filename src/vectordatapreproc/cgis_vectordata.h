/***************************************************************************
                          cgis_vectordata.h  -  description
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

#ifndef CGIS_VECTORDATA_H
#define CGIS_VECTORDATA_H

#include <assert.h>
#include <vector>
using namespace std;

#include "cpolyline.h"

/**This class holds a vector of PolyLine objects. 
  *@author Amrita Chanda
  */

class CGIS_VectorData {
public: 
  CGIS_VectorData();
  CGIS_VectorData(double res, double xmin, double xmax, double ymin, double ymax, int rows, int cols,
		  char *GIS_Dbase, char *GIS_Maplocation, char *GIS_Mapset, char *GIS_VectorFile);
  ~CGIS_VectorData();

  int ReadVectorData();

  int ProjectOn2DGrid();

  int StoreInFile();

  int ReadFromFile();

private:
  vector<CPolyLine *> _polylines;
  double _res, _xmin, _xmax, _ymin, _ymax;
  int _npolylines, _rows, _cols;
  char _gisDBase[256], _gisMaplocation[256], _gisMapset[256], _gisVectorFile[256]; 
};

#endif
