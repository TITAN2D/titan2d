/***************************************************************************
                          cpolyline.h  -  description
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

#ifndef CPOLYLINE_H
#define CPOLYLINE_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
using namespace std;

/**This class represents a line made up of many segments. It hold a list of xy points on the line.
  *@author Amrita Chanda
  */


class CPolyLine {
public: 
  CPolyLine(int type, string lbl);
  CPolyLine();
  ~CPolyLine();

  void SetLinesXY(int npts, double *x, double *y);

  int Get_ith_xy(int idx, float *x, float *y);

  int WritePolyLine(ofstream *stream);

  int ReadPolyLine(ifstream *stream);

  int _numOfPts;
  int _type;
  float **_xylst;
  string _label;
    
};

#endif
