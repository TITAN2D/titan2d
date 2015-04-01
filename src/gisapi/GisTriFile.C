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
 * $Id: GisTriFile.C,v 1.2 2005/02/02 00:09:54 namikawa Exp $  
 */ 
 
#ifdef HAVE_CONFIG_H 
# include <config.h> 
#endif 
  
#ifdef WIN32 
#pragma warning ( disable: 4786 ) 
#endif 
 
#include <math.h> 
#include "GisTriFile.h" 
 
double 
GisTriOut::xMin() 
{ 
	double min = x0_; 
	if (min > x1_) 
		min = x1_; 
	if (min > x2_) 
		min = x2_; 
	return min; 
} 
 
double 
GisTriOut::yMin() 
{ 
	double min = y0_; 
	if (min > y1_) 
		min = y1_; 
	if (min > y2_) 
		min = y2_; 
	return min; 
} 
 
double 
GisTriOut::xMax() 
{ 
	double max = x0_; 
	if (max < x1_) 
		max = x1_; 
	if (max < x2_) 
		max = x2_; 
	return max; 
} 
 
double 
GisTriOut::yMax() 
{ 
	double max = y0_; 
	if (max < y1_) 
		max = y1_; 
	if (max < y2_) 
		max = y2_; 
	return max; 
} 
 
bool 
GisTriOut::contains(double x, double y) 
{ 
	double	totalArea, triangleArea; 
 
//	Calculate the base triangle area 
	triangleArea = fabs( ((x1_ - x0_) * (y2_ - y0_)) - 
	 	                 ((x2_ - x0_) * (y1_ - y0_)) ); 
	triangleArea *= 1.00001; 
 
	totalArea = fabs( ((x0_ - x) * (y1_ - y)) - 
		              ((x1_ - x) * (y0_ - y)) ); 
	if (totalArea > triangleArea) 
		return false; 
 
	totalArea += fabs( ((x1_ - x) * (y2_ - y) ) - 
		               ((x2_ - x) * (y1_ - y) ) ); 
	if (totalArea > triangleArea)  
		return false; 
 
	totalArea += fabs( ((x0_ - x) * (y2_ - y) ) - 
		               ((x2_ - x) * (y0_ - y) )); 
	if (totalArea > triangleArea)  
		return false; 
 
	return true; 
} 
 
void 
GisTriOut::print() 
{ 
	cout.precision(12); 
	cout << " " << (int)x0_ << " " << (int)y0_ << " " << (int)x1_ << " " << (int)y1_ << " " << (int)x2_ << " " << (int)y2_ << endl; 
} 
 
GisTriFile::GisTriFile ( const string& name, const char* mode ): 
	GisBinFile ( name, mode ) 
{ 
	this->setWordSize(4); 
} 
 
int 
GisTriFile::versionNumber() 
{ 
	if ( this->good() ) 
	{ 
		int intVal; 
		this->read(&intVal); 
		return intVal; 
	} 
	return -1; 
} 
 
bool 
GisTriFile::readTimeStepInfo() 
{ 
	if ( this->good() ) 
	{ 
		float floatVal; 
		this->read(&numTri_); 
		this->read(&timeStep_); 
		this->read(&simTime_); 
		this->read(&pileMax_); 
		this->read(&pileMin_); 
		this->read(&xMomMax_); 
		this->read(&xMomMin_); 
		this->read(&yMomMax_); 
		this->read(&yMomMin_); 
		this->read(&floatVal); xMax_ = floatVal; 
		this->read(&floatVal); xMin_ = floatVal; 
		this->read(&floatVal); yMax_ = floatVal; 
		this->read(&floatVal); yMin_ = floatVal; 
		this->read(&elevMax_); 
		this->read(&elevMin_); 
		this->read(&maxVel_); 
		return true; 
	} 
	return false; 
} 
 
bool 
GisTriFile::readTriData(GisTriOut& triOut) 
{ 
	if ( this->good() ) 
	{ 
		float floatVal; 
		this->read(&triOut.key1_); 
		this->read(&triOut.key2_); 
		this->read(&triOut.gen_); 
 
		this->read(&floatVal); triOut.x0_ = floatVal; 
		this->read(&floatVal); triOut.y0_ = floatVal; 
		this->read(&triOut.z0_); 
		this->read(&floatVal); triOut.x1_ = floatVal; 
		this->read(&floatVal); triOut.y1_ = floatVal; 
		this->read(&triOut.z1_); 
		this->read(&floatVal); triOut.x2_ = floatVal; 
		this->read(&floatVal); triOut.y2_ = floatVal; 
		this->read(&triOut.z2_); 
 
		this->read(&triOut.pHeight_); 
		this->read(&triOut.xmom_); 
		this->read(&triOut.ymom_); 
 
		this->read(&triOut.istriedge_); 
		return true; 
	} 
	return false; 
} 
