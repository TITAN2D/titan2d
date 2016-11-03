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
 * $Id: GisGrid.C,v 1.2 2005/02/02 00:09:54 namikawa Exp $  
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif

#include <queue>
#include <math.h>
#include "GisGrid.h"
#include "GisTriFile.h"

GisGrid::GisGrid ()
{
  resX_ = resY_ = 0.;
  xMin_ = yMin_ = xMax_ = yMax_ = 0.;
  nCols_ = nRows_ = 0;
  noData_ = GISGRID_MAXFLOAT;
}

void
GisGrid::initGrid ()
{
  int gridSize = nCols_ * nRows_;
  if (gridSize > 0)
    gisGrid_.assign (gridSize, 0.);
}

void
GisGrid::setNumRowsCols (int nRows, int nCols)
{
  if (nRows > 0 && nCols > 0)
    {
      nRows_ = nRows;
      resY_ = this->deltaY () / nRows_;

      nCols_ = nCols;
      resX_ = this->deltaX () / nCols_;

      this->initGrid ();
    }
}
void
GisGrid::setNumRows (int nRows)
{
  if (nRows > 0)
    {
      nRows_ = nRows;
      resY_ = this->deltaY () / nRows_;
      this->initGrid ();
    }
}

void
GisGrid::setNumCols (int nCols)
{
  if (nCols > 0)
    {
      nCols_ = nCols;
      resX_ = this->deltaX () / nCols_;
      this->initGrid ();
    }
}

void
GisGrid::setRes (float res)
{
  if (res > 0.)
    {
      resX_ = res;
      nCols_ = (int) (this->deltaX () / resX_);
      resY_ = res;
      nRows_ = (int) (this->deltaY () / resY_);
      this->initGrid ();
    }
}

void
GisGrid::setResX (float res)
{
  if (res > 0.)
    {
      resX_ = res;
      nCols_ = (int) (this->deltaX () / resX_);
      this->initGrid ();
    }
}

void
GisGrid::setResY (float res)
{
  if (res > 0.)
    {
      resY_ = res;
      nRows_ = (int) (this->deltaY () / resY_);
      this->initGrid ();
    }
}

double
GisGrid::get (int row, int col)
{
  if ((row >= 0 && row < nRows_) && (col >= 0 && col < nCols_))
    {
      int index = (row * nCols_) + col;
      return gisGrid_[index];
    }
  return noData_;
}

void
GisGrid::set (int row, int col, float floatVal)
{
  if ((row >= 0 && row < nRows_) && (col >= 0 && col < nCols_))
    {
      int index = (row * nCols_) + col;
      gisGrid_[index] = floatVal;
    }
}

int
GisGrid::getRow (double y)
{
  int row = 0;
  if (resY_ > 0.)
    {
      row = (int) ((yMax_ - y) / resY_);

      if (row >= 0 && row < nRows_)
        return row;
    }
  return 0;
}

int
GisGrid::getCol (double x)
{
  int col = 0;
  if (resX_ > 0.)
    {
      col = (int) ((x - xMin_) / resX_);

      if (col >= 0 && col < nCols_)
        return col;
    }

  return 0;
}

double
GisGrid::getY (int row)
{
  double y;
  if (resY_ > 0.)
    y = yMax_ - (double) row *resY_;
  if (y >= yMin_ && y <= yMax_)
    return y;

  return yMin_;
}

double
GisGrid::getX (int col)
{
  double x;
  if (resX_ > 0.)
    x = xMin_ + (double) col *resX_;
  if (x >= xMin_ && x <= xMax_)
    return x;

  return xMin_;
}

void
GisGrid::setMax (GisTriOut & tri, int simVar)
{
  int fCol = this->getCol (tri.xMin ());
  int lCol = this->getCol (tri.xMax ());
  int fRow = this->getRow (tri.yMax ());
  int lRow = this->getRow (tri.yMin ());

  int i;
  for (i = fCol; i <= lCol; i++)
    {
      double x = this->getX (i);
      int j;
      for (j = fRow; j <= lRow; j++)
        {
          double y = this->getY (j);
          if (tri.contains (x, y))
            {
              if (simVar == 1)
                {
                  if (this->get (j, i) < tri.pHeight_)
                    this->set (j, i, tri.pHeight_);
                }
              else if (simVar == 2 && tri.pHeight_ > 0.2)
                {
                  float v =
                    (sqrt (tri.xmom_ * tri.xmom_ + tri.ymom_ * tri.ymom_)) /
                    tri.pHeight_;
                  if (this->get (j, i) < v)
                    this->set (j, i, v);
                }
            }
        }
    }
}

void
GisGrid::smooth3 ()
{
  queue < float >outBuf;
  int i;
  float fVal;
  for (i = 0; i < nRows_; i++)
    {
      int k1 = nCols_ * (i - 1);
      int k2 = nCols_ * i;
      int k3 = nCols_ * (i + 1);
      if (i == 0)
        k1 = k2;
      else if (i == nRows_ - 1)
        k3 = k2;

      int j;
      for (j = 0; j < nCols_; j++)
        {
          int m1 = j - 1;
          int m2 = j;
          int m3 = j + 1;
          if (j == 0)
            m1 = m2;
          else if (j == nCols_ - 1)
            m3 = m2;

          fVal = (gisGrid_[k1 + m1] + gisGrid_[k1 + m2] + gisGrid_[k1 + m3] +
                  gisGrid_[k2 + m1] + gisGrid_[k2 + m2] + gisGrid_[k2 + m3] +
                  gisGrid_[k3 + m1] + gisGrid_[k3 + m2] + gisGrid_[k3 +
                                                                   m3]) / 9.;
          outBuf.push (fVal);
          if (outBuf.size () == nCols_ + 2)
            {
              fVal = outBuf.front ();
              outBuf.pop ();
              gisGrid_[k1 + m1] = fVal;
            }
          if (k2 + m2 == gisGrid_.size () - 1)
            {
              int n;
              for (n = k1 + m1; n <= k2 + m2; n++)
                {
                  fVal = outBuf.front ();
                  outBuf.pop ();
                  gisGrid_[k1 + m1] = fVal;
                }
            }

        }

    }
}

void
GisGrid::print ()
{
  int col = 0;
  vector < float >::iterator i;
  for (i = gisGrid_.begin (); i != gisGrid_.end (); i++)
    {

      if ((*i) < 1.)
        continue;

      cout << (*i) << " ";
      col++;
      if (col >= nCols_)
        {
          col = 0;
          cout << endl;
        }
    }
}
