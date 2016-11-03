/***************************************************************************
                          cpolyline.cpp  -  description
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

#include <iostream>
#include <stdlib.h>
#include <list>
#include <cassert>
using namespace std;

#include <string.h>

#include "cpolyline.h"

CPolyLine::CPolyLine(int type, string lbl){
  _type = type;
  _label = lbl;
}

CPolyLine::CPolyLine(){
}


CPolyLine::~CPolyLine(){
}

/** No descriptions */
void CPolyLine::SetLinesXY(int npts, double *x, double *y){
  _numOfPts = npts;
  int k;

  _xylst = (float**)malloc(_numOfPts * sizeof(float*));

  for( k = 0; k < _numOfPts; ++k)
  {
    _xylst[k] = (float*)malloc(2 * sizeof(float));
    
    _xylst[k][0] = (float)x[k];
    _xylst[k][1] = (float)y[k];
 
  } 

}

int CPolyLine::Get_ith_xy(int idx, float *x, float *y){
  
  if(idx > _numOfPts-1 || idx < 0)
    return 1;
  
  float *xy = (float*)_xylst[idx];

  *x = xy[0];
  *y = xy[1];

  return 0;
}

int CPolyLine::WritePolyLine(ofstream *stream)
{
  float *xy;

  for(int i=0; i<_numOfPts; ++i)
    {
      xy = (float*)_xylst[i];
      
      if(xy == NULL)
	{

	  cout<<"ERROR: Unable to get "<<i<<"-th XY values from Polyline\n";
	  return 1;

	}

      *stream<<xy[0]<<","<<xy[1]<<endl;
    }

  return 0;
}


int CPolyLine::ReadPolyLine(ifstream *stream)
{
  list<double> xlst;
  list<double> ylst;
  
  char line[256], type[100];

  (*stream)>>_type;

  (*stream)>>_label;

  (*stream)>>line;

  while(strcmp(line, "END") != 0){

    char *xval, *yval;

    xval = strtok(line,",");
    yval = strtok(NULL," ");

    assert(xval != NULL);
    assert(yval != NULL);

    xlst.push_back(atof(xval));
    ylst.push_back(atof(yval));
  
    (*stream)>>line;

  }//while()

  _numOfPts = xlst.size();

  _xylst = (float**)malloc(_numOfPts * sizeof(float*));

  
  list<double>::iterator itr_x;
  list<double>::iterator itr_y;
  int k=0;
  for(itr_x=xlst.begin(), itr_y=ylst.begin(); itr_x!=xlst.end(); ++itr_x, ++itr_y, ++k)
    {

      _xylst[k] = (float*)malloc(2 * sizeof(float));
      
      _xylst[k][0] = (float)*itr_x;
      _xylst[k][1] = (float)*itr_y;

    }
  return 0;

}
