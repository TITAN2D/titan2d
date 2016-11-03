/***************************************************************************
                          cgis_vectordata.cpp  -  description
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
#include <fstream>
using namespace std;

#include <string.h>

#include "cgis_vectordata.h"
#include "../gisapi/GisApi.h"

CGIS_VectorData::CGIS_VectorData()
{}

CGIS_VectorData::CGIS_VectorData(double res, double xmin, double xmax, double ymin, double ymax, int rows, int cols,
				 char *GIS_Dbase, char *GIS_Maplocation, char *GIS_Mapset, char *GIS_VectorFile){
  _res = res;
  strcpy(_gisDBase, GIS_Dbase);
  strcpy(_gisMaplocation, GIS_Maplocation);
  strcpy(_gisMapset, GIS_Mapset);
  strcpy(_gisVectorFile, GIS_VectorFile);
}//Constructor


int CGIS_VectorData::ReadVectorData()
{
  int status_vec;
  int nlines;
  int ltype;
  int npts;
  int i, j;
  double* line_x;
  double* line_y;
  string lbl;

  //Initialize the GIS vector data
  status_vec = Initialize_Vector_data(_gisDBase, _gisMaplocation, _gisMapset, _gisVectorFile);
  
  if(status_vec == 0)
    {     
      if(Get_vector_n_lines(&nlines))
	{
	  
	  cout<<"Error retrieving number of lines from Vector Data"<<endl;
	  return 1;

	}

      _npolylines = nlines;
      
      for (i = 0; i < _npolylines; i++)
	{	 
	  if ( Get_vector_line_type(i, &ltype) )
	    {
	      
	      cout<<"Error retrieving line type from Vector Data"<<endl;
	      return 1;
	      
	    }
	  
	  if (Get_vector_line_size(i, &npts) )
	    {
	      
	      cout<<"Error retrieving line size from Vector Data"<<endl;
	      return 1;
	      
	    }
	  	  
	  if (Get_vector_line_label(i, &lbl) )
	    {
	      
	      cout<<"Error retrieving line label from Vector Data"<<endl;
	      return 1;
	      
	    }
	  
	  //	  cout<<"in my zone .."<<lbl<<endl;

	  if (npts > 0)
	    {
	      CPolyLine *pl = new CPolyLine(ltype, lbl);
	      _polylines.push_back(pl);
	      
	      line_x = new double[npts];
	      line_y = new double[npts];
	      
	      if ( Get_vector_line(i, line_x, line_y) )
		{
		  
		  cout<<"Error retrieving XY Vector Data"<<endl;
		  delete [] line_x;
		  delete [] line_y;
		  return 1;
		  
		}
	      
	      pl->SetLinesXY(npts, line_x, line_y);
	      
	      delete [] line_x;
	      delete [] line_y;
	      
	    }
	  
	}//for(i)
      
    }else{
    
    cout<<"\n\nError Initializing Vector Data : "<<_gisDBase<<":"<<_gisMaplocation<<":"<<_gisMapset<<":"<<_gisVectorFile<<endl;
    return 1;
    
  }//if-else(status_vec == 0)
  
  if(Delete_Vector_data()){
    cout<<"\n\nERROR: Unable to delete Vector data";
    return 1;
  }
  
  cout<<"DONE WITH CGIS_VectorData ReadVectorData()!!"<<endl;
  
  return 0;
  
}//ReadVectorData()


CGIS_VectorData::~CGIS_VectorData(){
}


/** This method takes each polyline and projects it on a 2D field.
    Intersections of each segment making up a polyline, with X & Y axes
    and the diagonals in the 2D field are noted.  */
int CGIS_VectorData::ProjectOn2DGrid(){
  /*
  list<float *>::iterator itr;
  int j;

  for(int i=0; i<_polylines.size(); ++i){
    j=0;
    itr = _polylines.begin();
    
    while(itr != _polylines.end()){
      grid->GetElevation(((CPolyLine*)itr)->_xyzList[j]);
      ++j;
    }// while()
    
  }//for(i<_polylines.size())
  */

   return 0;               
}

/** Stores XY values of 2D projected polylines in file "VectorDataOutput.dat" */
int CGIS_VectorData::StoreInFile(){

  ofstream vecoutfile("VectorDataOutput.data", ofstream::out);
   
  vecoutfile<<_npolylines<<endl;
                                                                                                                  
  for(int j=0; j<_npolylines; ++j){
    
    if(_polylines[j] == NULL)
      {
	
	cout<<"ERROR: Unable to access"<< j<<"-th Polyline\n";
	return 1; 
	
      }
    
    vecoutfile<<_polylines[j]->_type<<endl;
    vecoutfile<<_polylines[j]->_label<<endl;
    
    _polylines[j]->WritePolyLine(&vecoutfile);
    
    vecoutfile<<"END"<<endl;
    
  }
  
  vecoutfile.close();
  
  return 0;
}


int CGIS_VectorData::ReadFromFile()
{
  ifstream vecinfile("VectorDataOutput.data", ifstream::in);
  vecinfile>>_npolylines;
  
  //  cout<<"\n\n npolylines="<<_npolylines<<endl;

  for(int j=0; j<_npolylines; ++j){

    CPolyLine *pl = new CPolyLine();
    _polylines.push_back(pl);

    //    cout<<"#############"<<endl;

    pl->ReadPolyLine(&vecinfile);

  }

  vecinfile.close();
  
  return 0;

}
