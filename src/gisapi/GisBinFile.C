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
 * $Id: GisBinFile.C 214 2009-07-14 22:28:13Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 #ifdef WIN32
#pragma warning ( disable: 4786 )
#endif
#include <cstdlib>
using namespace std;

#include "GisBinFile.h"
#include "zlib.h"
#include <string.h>

// -- Constructor

GisBinFile::GisBinFile ( const string& name, const char* mode ):
  fileName_ ( name ), mode_ ( mode ), swappMode_ (false),
  wordSize_ (2), isInteger_ (true), compressed_ (false)
{
  const char* filename = fileName_.c_str();
  if ( mode_ == string("r") )
    file_.open ( filename, ios::in | ios::binary );
  else
    file_.open ( filename, ios::out | ios::binary );
  //	if ( ! file_.good() )
}

void
GisBinFile::setEndian (const char* mode)
{
  swappMode_ = true;
  if ( mode == string("little") )
    {
      if ( Gis_is_little_endian() ) // Is hardware little endian?
	swappMode_ = false;
    }
  else if ( ! Gis_is_little_endian() )// mode not Big and system
    swappMode_ = false;				// Big Endian
}

bool
GisBinFile::read(short* shortValue)
{
  if ( file_.good() )
    {
      if ( ( wordSize_ == 2 ) && isInteger_ )
	{
	  if ( file_.read((char*)shortValue,wordSize_) )
	    {
	      if (swappMode_)
		*shortValue = ( ((*shortValue >> 8) & 0x00FF) +
				((*shortValue << 8) & 0xFF00) );
	    }
	  return true;
	}
    }
  return false;
}

bool
GisBinFile::read4Bytes(char* char4Values)
{
  if (swappMode_)
    {
      if ( file_.get(tempStore_[3]) )
	if ( file_.get(tempStore_[2]) )
	  if ( file_.get(tempStore_[1]) )
	    if ( file_.get(tempStore_[0]) )
	      {
		memcpy (char4Values, tempStore_, 4);
		return true;
	      }
    }
  else if ( file_.read(char4Values,4) )
    return true;
  return false;
}

bool
GisBinFile::read(int* intValue)
{
  if ( file_.good() )
    if ( wordSize_ == 4 )
      return this->read4Bytes((char*) intValue);
  return false;
}

bool
GisBinFile::read(float* floatValue)
{
  if ( file_.good() )
    if ( wordSize_ == 4 )
      return this->read4Bytes((char*) floatValue);
  return false;
}

bool
GisBinFile::read (off_t *offsetVal )
{
  char *buf, *b;
  unsigned char c;
  off_t v;
  int i,sz;
  sz = sizeof(off_t);

  if ( file_.good() )
    {
      buf=(char *) malloc(wordSize_);
      file_.read (buf,wordSize_);
      v = 0;
      b = buf;
      for ( i=0; i< wordSize_; i++ )
	{
	  c = *b++;
	  if ( wordSize_ > sz && i < wordSize_ - sz && c != 0 )
	    return false;
	  v <<= 8;
	  v += c;
	}
      *offsetVal = v;
      free (buf);
      return true;
    }
  return false;
}

bool
GisBinFile::readNChar (char* charValues, int nValues)
{
  if ( file_.good() )
    if ( file_.read ( charValues, nValues ) )
      return true;
  return false;
}

void
GisBinFile::gotoPos(long newPos)
{
  file_.seekg( newPos );
}

bool
GisBinFile::readRow(int row, float* floatValues)
{
  if ( file_.good() )
    {
      this->gotoPos( 0L );
      if (this->isCompressed())
	return this->readCompressdRow(row, floatValues);
      else
	{
	  char nBytes = 4;
	  int totBytes = this->nCols()*nBytes;
	  this->gotoPos( row*totBytes );
	  unsigned char* readValues = new unsigned char[totBytes];
	  this->readNChar ((char *)readValues, totBytes);
	  
	  int i = 0;
	  int j = 0;
	  while ( i < totBytes )
	    {
	      char fvalchar[4];
	      for ( int k = 0; k < 4 ; k++)
		{
		  if (swappMode_)
		    fvalchar[3-k] = readValues[i++];
		  else
		    fvalchar[k] = readValues[i++];
		}
	      floatValues[j++] = *((float*)&fvalchar[0]);
	    }
	  delete[] readValues;
	  return true;
	}
    }
  return false;
}

bool
GisBinFile::readRow(int row, char* charValues)
{
  if ( file_.good() )
    {
      this->gotoPos( 0L );
      if (this->isCompressed())
	return this->readCompressdRow(row, charValues);
      else
	{
	  int totBytes = this->nCols()*dataSize_;
	  this->gotoPos( row*totBytes );
	  this->readNChar (charValues, totBytes);
	  return true;
	}
    }
  return false;
}

bool
GisBinFile::readCompressdRow(int row, float* floatValues)
{
  char nBytes;
  if ( file_.get(nBytes) )
    {
      if ( (int)nBytes == 4  || (int) nBytes == 8 )
	{
	  setWordSize ((int)nBytes);
	  int totBytes = nCols()*dataSize_;
	  unsigned char* expandedValues = new unsigned char[totBytes];
	  gotoPos( row*nBytes + 1L );
	  int rowSize;
	  if ( (int) nBytes == 4 )
	    {
	      int rowPtr, nextRowPtr;
	      read (&rowPtr);
	      read (&nextRowPtr);
	      gotoPos( rowPtr );
	      rowSize = nextRowPtr - rowPtr -1;
	    }
	  else
	    {
	      off_t rowPtr, nextRowPtr;
	      read (&rowPtr);
	      read (&nextRowPtr);
	      gotoPos ( rowPtr );
	      rowSize = nextRowPtr - rowPtr -1;
	    }
	  char compressFlag;
	  file_.get(compressFlag);
	  if ( compressFlag == 0x31 )
	    {
	      unsigned char* charValues = new unsigned char[rowSize];
	      this->readNChar ((char *)charValues, rowSize);
	      if ( ! GisUncompress ( rowSize, charValues, 
				     totBytes, expandedValues ) )
		{
		  delete[] charValues;
		  return false;
		}
	      delete[] charValues;
	    }
	  else
	    this->readNChar ((char *)expandedValues, totBytes);
	  int i = 0;
	  int j = 0;
	  while ( i < totBytes )
	    {
	      char fvalchar[4];
	      for ( int k = 0; k < 4; k++)
		{
		  if (swappMode_)
		    fvalchar[3-k] = expandedValues[i++];
		  else
		    fvalchar[k] = expandedValues[i++];
		}
	      floatValues[j++] = *((float*)&fvalchar[0]);
	    }
	  delete[] expandedValues;
	  return true;
	}
    }
  return false;
}

bool
GisBinFile::readCompressdRow(int row, char* charValues)
{
  char nBytes;
  if ( file_.get(nBytes) )
    {
      if ( (int)nBytes == 4 || (int) nBytes == 8 )
	{
	  setWordSize ((int)nBytes);
	  int totBytes = nCols()*dataSize_;
	  gotoPos( row*nBytes + 1L );
	  int rowSize;
	  if ( (int) nBytes == 4 )
	    {
	      int rowPtr, nextRowPtr;
	      read ( &rowPtr );
	      read ( &nextRowPtr );
	      gotoPos ( rowPtr );
	      rowSize = nextRowPtr - rowPtr -1;
	    }
	  else
	    {
	      off_t rowPtr, nextRowPtr;
	      read (&rowPtr);
	      read (&nextRowPtr);
	      gotoPos( rowPtr );
	      rowSize = nextRowPtr - rowPtr;
	    }
	  char compressFlag;
	  file_.get(compressFlag);
	  //			if ( rowSize < totBytes )
	  if  ( ( compressFlag == 0x01 ) && ( (rowSize-1) < totBytes ) )
	    {
	      unsigned int i = 0;
	      while ( i < totBytes )
		{
		  char charCount;
		  file_.get(charCount);
		  char charValue;
		  file_.get(charValue);
		  unsigned int j = 0;
		  for ( j = i; j < i+(unsigned char)charCount; j++ )
		    charValues[j] = charValue;
		  i += (unsigned char)charCount;
		}
	    }
	  else
	    this->readNChar (charValues, totBytes);
	  return true;
	}
    }
  return false;
}

bool
GisUncompress ( int inSize, unsigned char* inChars, int outSize, unsigned char* outChars )
{
  z_stream streamctrl;
  streamctrl.zalloc = (alloc_func)0;
  streamctrl.zfree  = (free_func)0;
  streamctrl.opaque = (voidpf)0;
  streamctrl.avail_in  = inSize;
  streamctrl.next_in   = inChars;
  int decompSize = outSize;
  streamctrl.avail_out = decompSize;
  streamctrl.next_out  = outChars;
  
  int err = inflateInit (&streamctrl); //zlib init
  if (err == Z_OK)
    {
      err = inflate (&streamctrl, Z_FINISH);
      int availBytes = decompSize - streamctrl.avail_out;
      if (!(err == Z_STREAM_END || err == Z_OK))
	{	//Some error
	  if (!(err == Z_BUF_ERROR && availBytes == decompSize))
	    inflateEnd (&streamctrl);
	}
      else
	{ //No error
	  inflateEnd (&streamctrl);
	  return true;
	}
    }
  return false;
}

bool Gis_is_little_endian()
{
  union
  {
    int anyInt;
    char anyChar[sizeof(int)];
  } testUnion;
  testUnion.anyInt = 1;
  if (testUnion.anyChar[0] == 1)
    return true;
  else
    return false;
}

