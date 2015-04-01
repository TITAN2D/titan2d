#ifndef GisColors_H
#define GisColors_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class GisColors
{
public:

  GisColors (const string & name);

    virtual ~ GisColors ()
  {
  }

  bool good ()
  {
    return _status;
  }

  void getColor (int index, unsigned char &red, unsigned char &green,
                 unsigned char &blue);

  void print ();

protected:

  bool _status;
  unsigned int _red[256];
  unsigned int _green[256];
  unsigned int _blue[256];

private:

// No copy allowed
  GisColors (const GisColors &);
  GisColors & operator= (const GisColors &);
};

#endif
