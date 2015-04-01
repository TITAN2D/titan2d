/***************************************************************************
                          cgis_trigrid.h  -  description
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

#ifndef CGIS_TRIGRID_H
#define CGIS_TRIGRID_H

/**This class creates an array of elevations. 
  It answers the queries: Which triangle(s) contain point (x,y)?
  and which triangles fall within the bounding area (x1,y1) and (x2,y2)?
  *@author Amrita Chanda
  */

class CGIS_TriGrid {
public: 
	CGIS_TriGrid(char *GIS_Dbase, char *GIS_Maplocation, char *GIS_Mapset, char *GIS_Rasterimage);
	~CGIS_TriGrid();

  double *_elevations;
  double _xmax, _xmin, _ymax, _ymin;
  int _num_of_rows;
  int _num_of_cols;
  double _res;

};

#endif
