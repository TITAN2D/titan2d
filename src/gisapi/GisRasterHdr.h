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
 * $Id: GisRasterHdr.h 12 2003-11-07 17:58:49Z kdalbey $ 
 */

#ifndef GisRasterHdr_H
#define GisRasterHdr_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class GisRasterHdr
{
public:

	GisRasterHdr(const string& name);
	
	virtual ~GisRasterHdr(){} 

	bool isCompressed()
	{ return (_compressed == 1); }

	int Rows()
	{ return _rows; }

	int Cols()
	{ return _cols; }

	double XRes()
	{ return _ewresol;}

	double YRes()
	{ return _nsresol;}

	double North()
	{ return _north;}

	double South()
	{ return _south;}

	double East()
	{ return _east;}

	double West()
	{ return _west;}

	bool good()
	{ return _status; }

	int cellFormat()
	{ return _formatId; }

	void print();

protected:

	int _projId; //0=XY,1=UTM,2=SP(spherical?),3=LL(LatLong),99=other
	int _zoneId; //UTM zone: from (-180 + (zone-1)*6) to (-180 + zone*6)
	double _north;
	double _south;
	double _east;
	double _west;
	int _cols;
	int _rows;
	double _ewresol; //E-W direction resolution - Cell size
	double _nsresol; //N-S direction resolution - Cell size
	int _formatId;   //  0 means 1 byte, 1 means 2 bytes, etc. -1 indicates a floating point
	int _compressed; // 0 = uncompressed, 1 = compressed
	bool _status; // true = OK, 1false = some error

private:
	
// No copy allowed
	GisRasterHdr(const GisRasterHdr&);
	GisRasterHdr& operator=(const GisRasterHdr&);
};

#endif
