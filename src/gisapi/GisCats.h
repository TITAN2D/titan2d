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
 * $Id: GisCats.h 14 2003-11-21 16:21:49Z kdalbey $ 
 */

#ifndef GisCats_H
#define GisCats_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

class GisCats
{
public:

	GisCats(const string& name);
	
	virtual ~GisCats(){} 

	bool good()
	{ return _status; }

	int mumberOfCats()
	{ return _ncats; }

	const char* category(int index)
	{ return _catnames[index].c_str(); }

	void print();

protected:

	bool _status;

	int _ncats;
	map<int,string> _catnames;

private:
	
// No copy allowed
	GisCats(const GisCats&);
	GisCats& operator=(const GisCats&);
};

#endif
