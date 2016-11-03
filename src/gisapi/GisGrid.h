#ifndef GisGrid_H
#define GisGrid_H

#define GISGRID_MAXFLOAT 1.0e31
#define GISGRID_BIGFLOAT 1.0e32

#include <vector>
using namespace std;

class GisTriOut;

class GisGrid
{
public:

  GisGrid ();

  virtual ~ GisGrid ()
  {
  }

  void setNumRowsCols (int nRows, int nCols);

  void setNumRows (int nRows);

  void setNumCols (int nCols);

  void setBox (double xMin, double yMin, double xMax, double yMax)
  {
    xMin_ = xMin;
    yMin_ = yMin;
    xMax_ = xMax;
    yMax_ = yMax;
  }

  void setRes (float res);

  void setResX (float res);

  double getResX ()
  {
    return resX_;
  }

  void setResY (float res);

  double getResY ()
  {
    return resY_;
  }

  double deltaX ()
  {
    return (xMax_ - xMin_);
  }

  double deltaY ()
  {
    return (yMax_ - yMin_);
  }

  double noDataValue ()
  {
    return noData_;
  }

  void noDataValue (double doubleVal)
  {
    noData_ = doubleVal;
  }

  void initGrid ();

  void setMax (GisTriOut & tri, int simVar);
  //    simVar: 1 - Pile height; 2 - Velocity";

  double get (int row, int col);

  int getRow (double y);

  int getCol (double x);

  double getY (int row);

  double getX (int col);

  void set (int row, int col, float floatVal);

  void smooth3 ();

  void print ();

  // protected: 

  double resX_;
  double resY_;
  double xMin_;
  double yMin_;
  double xMax_;
  double yMax_;
  double noData_;
  int nCols_;
  int nRows_;

  vector < float >gisGrid_;

private:

};

#endif
