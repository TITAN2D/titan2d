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
 * $Id: GisLines.h 120 2007-06-07 19:21:02Z dkumar $ 
 */

#ifndef GisLines_H
#define GisLines_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class GisLabels;

class GisLine
{
public:

	GisLine();
	GisLine(vector<double>& xS, vector<double>& yS);
	
	virtual ~GisLine(){} 

	int numberOfPoints()
	{ return _nPoints; }

	void getPoints(double* xCoords, double* yCoords);

	void insertPoints (double x, double y);

	void setIndex (int lineIndex)
	{ _lineIndex = lineIndex; }
	
	void setLabel (string lineLabel)
	{ _lineLabel = lineLabel; }

	int getIndex ()
	{ return _lineIndex; }

	string getLabel ()
	{ return _lineLabel; }
	
	bool findPoints(double x, double y, double scale);

	void print();

	vector<double> _xS;
	vector<double> _yS;
	GisLine(const GisLine& a)
	{
		_xS = a._xS;
		_yS = a._yS;
		_lineIndex = 0;
		_nPoints = (int)_xS.size();
	}

	GisLine& operator=(const GisLine& a)
	{
		if ( this != &a )
		{	
			_xS = a._xS; _yS = a._yS;
			_lineIndex = 0;
			_nPoints = (int)_xS.size();
		}
		return *this;
	}

protected:
	int _nPoints;
	int _lineIndex;
	string _lineLabel;
private:

};

class GisLines
{
public:

	GisLines(const string& name);
	
	virtual ~GisLines(){} 

	bool good()
	{ return _status; }

	int numberOfLines()
	{ return _nLines; }

	void insertLine (const GisLine& readLine)
	{ _lines.push_back(readLine); }

	int getLineSize(int nLine, int* lineSize);

	int getIndex(int nLine, int* lineIndex);

	int getLabel(int nLine, string* lineStr);

	int getLine(int nLine, double* x, double* y);

	void setIndices(GisLabels& labels, double scale);

	void print();

protected:

	bool _status;

	int _nLines;
	vector<GisLine> _lines;

private:
	
// No copy allowed
	GisLines(const GisLines&);
	GisLines& operator=(const GisLines&);
};

#endif
