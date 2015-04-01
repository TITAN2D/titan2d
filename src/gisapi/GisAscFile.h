#ifndef GisAscFile_H
#define GisAscFile_H
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class GisAscFile
{
public:

  GisAscFile (const string & name, const char *mode = "r");

  virtual ~ GisAscFile ()
  {
  }

  bool good ()
  {
    return file_.good ();
  }

  void getLine (char *outChar, int nChar, char termChar) //'\n' = 0x0A
  {
    file_.getline (outChar, nChar, termChar);
  }

  void getLine (char *outChar, int nChar)
  {
    file_.getline (outChar,nChar,'\n');
  }

  void getLine (string & outString)
  {
    getline (file_, outString);
  }

  void getChar (char *outChar)
  {
    file_.read (outChar, 1);
  }

  void getAscInt (int &intValue)
  {
    file_ >> intValue;
  }

  void getAscDouble (double &doubleValue)
  {
    file_ >> doubleValue;
  }

  void getString (string & rString)
  {
    file_ >> rString;
  }

  void rewind ();

  string mode ()
  {
    return mode_;
  }

  bool findString (const string & toFindString);

protected:

// -- File pointer

  fstream file_;
  string mode_;

private:

// No copy allowed
  GisAscFile (const GisAscFile &);
  GisAscFile & operator= (const GisAscFile &);
};

#endif
