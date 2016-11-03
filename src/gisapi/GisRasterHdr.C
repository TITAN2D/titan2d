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
 * $Id: GisRasterHdr.C 206 2009-01-26 17:32:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <string.h>

#include "GisRasterHdr.h"
#include "GisAscFile.h"

GisRasterHdr::GisRasterHdr ( const string& name )
{
	_status = false;
	GisAscFile headerFile( name );
	if (headerFile.good())
	{
		char charText[20];
		headerFile.getLine(charText, 20, ':');
		if ( strcmp(charText, "proj") == 0)
		{ //At least test first line!!!
			headerFile.getAscInt ( _projId );
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscInt ( _zoneId );
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_north);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_south);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_east);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_west);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscInt (_cols);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscInt (_rows);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_ewresol);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscDouble (_nsresol);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscInt (_formatId);
			headerFile.getLine(charText, 20, ':');
			headerFile.getAscInt (_compressed);
			_status = true;
		}
	}
}

void
GisRasterHdr::print()
{
	cout << "projId " << _projId << "\n";
	cout << "zoneId " << _zoneId << "\n";
	cout << "north " << _north << "\n";
	cout << "south " << _south << "\n";
	cout << "east " << _east << "\n";
	cout << "west " << _west << "\n";
	cout << "cols " << _cols << "\n";
	cout << "rows " << _rows << "\n";
	cout << "ewres " << _ewresol << "\n";
	cout << "nsres " << _nsresol << "\n";
	cout << "formatId " << _formatId << "\n";
	cout << "compressed " << _compressed << "\n";
}

