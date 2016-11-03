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
 * $Id: GisAscFile.C 109 2005-02-02 00:04:41Z namikawa $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 

#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif
#include "GisAscFile.h"

// -- Constructor

GisAscFile::GisAscFile ( const string& name, const char* mode ):
	mode_ ( mode )
{
	const char* filename = name.c_str();
	long rw_flag = 0;
	if ( mode_ == string("r") )
		//file_.open ( filename, ios_base::in );
		file_.open ( filename, ios::in );
	else
		//file_.open ( filename, ios_base::out );
		file_.open ( filename, ios::out );
//	if ( ! file_.good() )
}

void
GisAscFile::rewind()
{
	if ( file_.good() )
	{
		if ( mode_ == string("r") ) //valid only for reading ??
			file_.seekg(0, ios::beg );
	}
}

bool
GisAscFile::findString(const string& toFindString)
{
	bool stringFound = false;
	if ( file_.good() )
	{
		string inString;
		for (;;)
		{
			file_ >> inString;
			if ( inString == string(toFindString) )
			{
				stringFound = true;
				break;
			}
		}
	}
	return stringFound;
}

