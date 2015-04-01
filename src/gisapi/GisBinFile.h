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
 * $Id: GisBinFile.h 211 2009-06-16 20:02:10Z dkumar $ 
 */

#ifndef GisBinFile_H
#define GisBinFile_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

bool GisUncompress (int inSize, unsigned char* inChars, int outSize, unsigned char* outChars );

bool Gis_is_little_endian();

class GisBinFile
{
public:
	GisBinFile( const string& name, const char* mode = "r");

	virtual ~GisBinFile(){};

	void setEndian (const char* mode = "little");

	void setWordSize (int wordSize)
	{ wordSize_ = wordSize; }

	void setDataSize (int dataSize)
	{ dataSize_ = dataSize; }

	void setIsInteger (bool trueFalse)
	{ isInteger_ = trueFalse; }

	bool readRow(int row, float* floatValues);
	bool readRow(int row, char* charValues);

	bool readCompressdRow(int row, float* floatValues);
	bool readCompressdRow(int row, char* charValues);

	bool read (float* floatValue);	//IEEE 32 bits
	
	bool read (off_t* offsetVals);	// read row pointers

	bool read (int* intValue);		// 32 bits

	bool read (short* shortValue); //16 bits

	bool read4Bytes(char* char4Values);
	
	bool readNChar (char* charValues, int nValues);

	bool good()
	{ return file_.good(); }

	int nRows ()
	{ return nRows_; }

	int nCols ()
	{ return nCols_; }

	void nRows (int nrows)
	{ nRows_ = nrows; }

	void nCols (int ncols)
	{ nCols_ = ncols; }

	bool isCompressed()
	{ return ( compressed_ == 1); }

	void isCompressed( bool compFlag )
	{ compressed_ = compFlag; }

	void gotoPos(long newPos);

protected:
// -- File pointer
	fstream file_;
	string fileName_;
	string mode_;
	bool swappMode_;
	int wordSize_;
	int dataSize_;
	bool isInteger_;
	char tempStore_[4];
	int nRows_;
	int nCols_;
	bool compressed_;

private:
// No copy allowed
	GisBinFile(const GisBinFile&);
	GisBinFile& operator=(const GisBinFile&);
};
#endif
