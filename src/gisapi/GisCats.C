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
 * $Id: GisCats.C 109 2005-02-02 00:04:41Z namikawa $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif

#include "GisCats.h"
#include "GisAscFile.h"



GisCats::GisCats ( const string& name )
{
	_status = false;
	GisAscFile catsFile( name );
	if (catsFile.good())
	{
		char charText[200];
		catsFile.getLine(charText, 200, '#');	//#
		catsFile.getAscInt( _ncats );	// number of categories
		catsFile.getLine(charText, 200);	//categories
		catsFile.getLine(charText, 200);	//Geology
		catsFile.getLine(charText, 200);	//Empty line
		catsFile.getLine(charText, 200);	//0.00 0.00 0.00 0.00

		while (	catsFile.good() )
		{
			string catName;
			int catNumber;
			catsFile.getAscInt( catNumber );	//category number
			catsFile.getLine(charText, 20, ':');
			catsFile.getString(catName);
			_catnames[catNumber] = catName;
			if (catNumber == _ncats)
				break;
		}
		//this->print();
		_status = true;
	}
}

void
GisCats::print()
{
	map<int,string>::iterator i;

	for( i = _catnames.begin(); i != _catnames.end(); i++)
        cout << (*i).first << " " << (*i).second << endl;
}

