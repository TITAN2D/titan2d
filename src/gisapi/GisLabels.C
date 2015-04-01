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
 * $Id: GisLabels.C 110 2005-02-02 00:09:54Z namikawa $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif

#include "GisLabels.h"
#include "GisSPRFile.h"

GisLabels::GisLabels ( const string& name )
{
	_status = false;
	_nlabels = 0;
	GisSPRFile sprFile( name );
	if (sprFile.good())
	{
		if ( sprFile.readLabels(_labelXs, _labelYs, _labelStrings) )
		{
			_nlabels = _labelStrings.size();
			_status = true;
		}

		vector<string>::iterator i;
		map<string,int>::iterator j;
		for( i = _labelStrings.begin(); i != _labelStrings.end(); i++ )
		{
			j = _labelIndxs.find(*i);
			if (j == _labelIndxs.end() )
				_labelIndxs[(*i)] = _labelIndxs.size();
		}
	}
}

bool
GisLabels::getLabel(int nlabel, string& labelString, double& x, double& y)
{
	if ( nlabel < _nlabels )
	{
		labelString = _labelStrings[nlabel];
		x = _labelXs[nlabel];
		y = _labelYs[nlabel];
		return true;
	}
	return false;
}

int
GisLabels::getIndex(string& labelStr)
{
	map<string,int>::iterator i;
	i = _labelIndxs.find(labelStr);
	if (i != _labelIndxs.end() )
		return ( (*i).second ) + 1;

	return 0;
}

void
GisLabels::print()
{
	vector<string>::iterator i;
	vector<double>::iterator j;
	vector<double>::iterator k;
	for( i  = _labelStrings.begin(), j  = _labelXs.begin(), k  = _labelYs.begin();
		 i != _labelStrings.end(),   j != _labelXs.end(),   k != _labelYs.end();
		 i++, j++, k++)
        cout << (*i) << " " << (*j) << " " << (*k) << endl;

	map<string,int>::iterator m;
	for( m = _labelIndxs.begin(); m != _labelIndxs.end(); m++)
        cout << (*m).first << " " << (*m).second << endl;

}

