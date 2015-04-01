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
 * $Id: GisLines.C 120 2007-06-07 19:21:02Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif

#include "GisLines.h"
#include "GisLabels.h"
#include "GisSPRFile.h"
#include <math.h>

GisLine::GisLine ()
{
	_nPoints = 0;
	_lineIndex = 0;
}

GisLine::GisLine(vector<double>& xS, vector<double>& yS)
{
	_xS = xS;
	_yS = yS;

	_lineIndex = 0;
	_nPoints = (int)_xS.size();
}

void
GisLine::getPoints(double* xCoords, double* yCoords)
{
	int k;
	for( k = 0; k < _nPoints; k++)
	{
	  xCoords[k] = _xS[k];
	  yCoords[k] = _yS[k];
	}
}

// Perpendicular distance from location pin to segment (first -> last),  
// intersection at pinter 
double perpendicularDistance(double xi, double yi, double xf, double yf, double x, double y) 
{ 
	double	d12; 
 
	double dx = xf-xi; 
	double dy = yf-yi; 
	double a2 = (y-yi)*dx-(x-xi)*dy; 
	 
    if(dx==0. && dy==0.) 
	{ 
		d12= sqrt(((x-xi)*(x-xi))+((y-yi)*(y-yi))); 
		d12 *= d12; 
	} 
	else 
		d12= a2*a2/(dx*dx+dy*dy); 
 
	return (sqrt(d12)); 
} 


bool
GisLine::findPoints(double x, double y, double scale)
{
	bool pointFound = false;
	int k;
	double linePrecision = 0.0002 * scale; //0.2mm of the scale, Change if you must

	if ( _nPoints > 1 )
	{
		for( k = 0; k < _nPoints; k++)
		{
			if ( fabs(_xS[k] - x) < linePrecision && fabs(_yS[k] - y) < linePrecision )
			{
				pointFound = true;
				break;
			}
		}
		for( k = 0; k < _nPoints - 1; k++)
		{
			if ( _xS[k] < _xS[k+1] )
			{
				if ( x < _xS[k] || x > _xS[k+1] )
					continue;
			}
			else if ( _xS[k] > _xS[k+1] )
			{
				if ( x > _xS[k] || x < _xS[k+1] )
					continue;
			}
			if ( _yS[k] < _yS[k+1] )
			{
				if ( y < _yS[k] || y > _yS[k+1] )
					continue;
			}
			else if ( _yS[k] > _yS[k+1] )
			{
				if ( y > _yS[k] || y < _yS[k+1] )
					continue;
			}
			double dist=perpendicularDistance(_xS[k],_yS[k],_xS[k+1],_yS[k+1],x,y);

			if ( dist < linePrecision )
			{
				pointFound = true;
				break;
			}
		}
	}
	return pointFound;
}

void
GisLine::insertPoints (double x, double y)
{
	_xS.push_back(x);
	_yS.push_back(y);
	_nPoints = (int)_xS.size();
}

void
GisLine::print()
{
	vector<double>::iterator i;
	vector<double>::iterator j;

	cout << "Index:" << _lineIndex << " with " << _nPoints << " points." << endl;
	cout.precision(12);

	for( i  = _xS.begin(), j  = _yS.begin();
		 i != _xS.end(),   j != _yS.end();
		 i++, j++)
        cout << "(" << (*i) << " " << (*j) << ") ";
}

GisLines::GisLines ( const string& name )
{
	GisSPRFile sprFile( name );
	if (sprFile.good())
	{
		vector<double> fXs;
		vector<double> fYs;
		if ( sprFile.readFirstLine(fXs, fYs) )
		{
			GisLine firstLine(fXs, fYs);
			this->insertLine(firstLine);

			_nLines = (int)_lines.size();

			for (;;)
			{
				vector<double> nXs;
				vector<double> nYs;

				if ( sprFile.readNextLine(nXs, nYs) )
				{
					GisLine nextLine(nXs, nYs);
					this->insertLine(nextLine);

					_nLines = (int)_lines.size();
				}
				else
					break;
			}
		}
	}
}

int
GisLines::getLine(int nLine, double* x, double* y)
{
	if ( nLine < _nLines )
	{
		_lines[nLine].getPoints(x, y);
		return 0;
	}
	return -4;
}

int
GisLines::getIndex(int nLine, int* lineIndex)
{
	if ( nLine < _nLines )
	{
		*lineIndex = _lines[nLine].getIndex();
		return 0;
	}
	return -4;
}

int
GisLines::getLabel(int nLine, string* lineStr)
{
	if ( nLine < _nLines )
	{
	  string lbl = _lines[nLine].getLabel();
	  *lineStr = lbl;
	  return 0;
	}
	return -4;
}


int
GisLines::getLineSize(int nLine, int* lineSize)
{
	if ( nLine < _nLines )
	{
		*lineSize = _lines[nLine].numberOfPoints();
		return 0;
	}
	return -4;
}

void
GisLines::setIndices(GisLabels& labels, double scale)
{
  //DEUBG
 	int nlabels = labels.numberOfLabels();
	if ( nlabels > 0 )
	{
	  string labelStr;
	  double x, y;
	  vector<GisLine>::iterator i;
	  for( i  = _lines.begin(); i != _lines.end(); i++)
	    {
	      int j;
	      for ( j = 0; j < nlabels; j++ )
		{
		  if ( labels.getLabel(j, labelStr, x, y) )
		    {
		      if ( (*i).findPoints(x,y,scale) )
			{
			  //cout << "found: " << j << endl;
			  (*i).setIndex(labels.getIndex(labelStr));
			  (*i).setLabel(labelStr);
			  
			  //cout<<"label found "<<labelStr<<endl;
			  break;
			}
		    }
		}
	    }
	}
}

void
GisLines::print()
{
	vector<GisLine>::iterator i;
	for( i  = _lines.begin(); i != _lines.end(); i++)
        (*i).print();
}

