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
 * $Id: GisLabels.h 110 2005-02-02 00:09:54Z namikawa $ 
 */

#ifndef GisLabels_H
#define GisLabels_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

class GisLabels
{
public:

	GisLabels(const string& name);
	
	virtual ~GisLabels(){} 

	bool good()
	{ return _status; }

	int numberOfLabels()
	{ return _nlabels; }

	bool getLabel(int nlabel, string& labelString, double& x, double& y);

	int getIndex(string& labelStr);

	void print();

protected:

	bool _status;

	int _nlabels;
	vector<string> _labelStrings;
	vector<double> _labelXs;
	vector<double> _labelYs;
	map<string,int> _labelIndxs;


private:
	
// No copy allowed
	GisLabels(const GisLabels&);
	GisLabels& operator=(const GisLabels&);
};

#endif
