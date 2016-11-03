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
 * $Id: GisSPRFile.C 110 2005-02-02 00:09:54Z namikawa $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 

#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif
#include "GisSPRFile.h"

GisSPRFile::GisSPRFile ( const string& name, const char* mode ):
	GisAscFile ( name, mode )
{
}

bool
GisSPRFile::gotoLINESSection()
{
	bool stringFound = false;
	this->rewind();
	while ( this->good() )
	{
		string inString;
		string outString;
		this->getLine(inString);
		istringstream inStrStream(inString);
		inStrStream >> outString;
		if ( outString == string("LINES") )
		{
			stringFound = true;
			break;
		}
	}
	return stringFound;
}

bool
GisSPRFile::gotoPOINTSSection()
{
	bool stringFound = false;
	this->rewind();
	while ( this->good() )
	{
		string inString;
		string outString;
		this->getLine(inString);
		istringstream inStrStream(inString);
		inStrStream >> outString;
		if ( outString == string("POINTS") )
		{
			stringFound = true;
			break;
		}
	}
	return stringFound;
}

bool
GisSPRFile::gotoSection(string& sectionName)
{
	bool stringFound = false;
	while ( this->good() )
	{
		string inString;
		string outString;
		this->getLine(inString);
		istringstream inStrStream(inString);
		inStrStream >> outString;
		if ( outString == sectionName )
		{
			stringFound = true;
			break;
		}
	}
	return stringFound;
}

bool
GisSPRFile::readINFOSection()
{
	string infoStr("INFO");
        if ( this->gotoSection (infoStr) )
	{
		bool infoSection = false;
		while ( this->good() )
		{
			string inString;
			string outString;

			this->getLine(inString);
			istringstream inStrStream(inString);
			inStrStream >> outString;

			if ( outString.find ("//") == 0)
				continue;

			if ( outString == string("INFO_END") )
			{
				infoSection = true;
				break;
			}

			if ( outString == string("SEPARATOR") )
			{
				inStrStream >>  _sepStr;
				continue;
			}
		}
		return infoSection;
	}
	return false;
}

bool
GisSPRFile::readFirstLine(vector<double>& x, vector<double>& y)
{
	bool lineFound = false;
	if ( this->gotoLINESSection() )
	{
		if ( this->readINFOSection() )
		{
			double xString, yString;
			string inString, outString;
			for (;;)
			{
				this->getLine(inString);
				istringstream inStrStream(inString);
				inStrStream >> outString;
				if ( outString == string("END") )
					break;

				inStrStream.seekg(0, ios::beg );
				inStrStream >> xString >> yString;

				x.push_back(xString);
				y.push_back(yString);
			}
			if ( x.size() > 0 )
				lineFound = true;
		}
	}
	return lineFound;
}

bool
GisSPRFile::readNextLine(vector<double>& x, vector<double>& y)
{
	bool lineFound = false;
	double xString, yString;
	string inString, outString;
	for (;;)
	{
		this->getLine(inString);
		istringstream inStrStream(inString);
		inStrStream >> outString;
		if ( outString == string("END") )
			break;

		inStrStream.seekg(0, ios::beg );
		inStrStream >> xString >> yString;

		x.push_back(xString);
		y.push_back(yString);
	}
	if ( x.size() > 0 )
		lineFound = true;

	return lineFound;
}

bool
GisSPRFile::readLabels(vector<double>& x, vector<double>& y, vector<string>& labelsStr)
{
	bool labelFound = false;
	if ( this->gotoPOINTSSection() )
	{
		if ( this->readINFOSection() )
		{
			double xString, yString;
			string labelStr;
			string nextString1, nextString2;
			string inString, outString;

			for (;;)
			{
				this->getLine(inString);
				istringstream inStrStream(inString);
				inStrStream >> outString;
				if ( outString == string("END") )
					break;

				inStrStream.seekg(0, ios::beg );
				if (_sepStr.size() > 0)
					inStrStream >> xString >> nextString1 >> yString >> nextString2 >> labelStr;
				else
					inStrStream >> xString >> yString >> labelStr;
				x.push_back(xString);
				y.push_back(yString);
				labelsStr.push_back(labelStr);
				labelFound = true;
			}
		}
	}
	return labelFound;
}
