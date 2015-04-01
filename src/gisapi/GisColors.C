#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif

#include "GisColors.h"
#include "GisAscFile.h"



GisColors::GisColors (const string & name)
{
  _status = false;
  GisAscFile colsFile (name);
  if (colsFile.good ())
    {
      char charText[200];
      int fcolor, lcolor;
      int i;

      for (i = 0; i <= 255; i++)
        _red[i] = _green[i] = _blue[i] = 0;

      colsFile.getLine (charText, 200, '%');    // % 0 255
      colsFile.getAscInt (fcolor);      // first color
      colsFile.getAscInt (lcolor);      // last color

      fcolor = (fcolor < 0) ? 0 : fcolor;
      lcolor = (lcolor > 255) ? 255 : lcolor;
      for (i = fcolor; i <= lcolor; i++)
        {
          int colNumber, red, green, blue;
          colsFile.getAscInt (colNumber);       //color number
          if (i != colNumber)
            {
              _status = false;
              return;
            }
          colsFile.getLine (charText, 20, ':');
          colsFile.getAscInt (red);
          colsFile.getChar (&charText[0]);
          if (charText[0] == 0x0a)
            {
              _green[i] = _blue[i] = _red[i] = red;
              continue;
            }
          colsFile.getAscInt (green);
          colsFile.getLine (charText, 20, ':');
          colsFile.getAscInt (blue);
          _green[i] = green;
          _blue[i] = blue;
          _red[i] = red;
        }
      _status = true;
    }
}

void
GisColors::getColor (int index, unsigned char &red, unsigned char &green,
                     unsigned char &blue)
{
  if (index >= 0 && index <= 255)
    {
      green = _green[index];
      blue = _blue[index];
      red = _red[index];
    }
  else
    {
      green = blue = red = 0;
    }
}

void
GisColors::print ()
{
//cout << colNumber << " " << red << " " << green << " " << blue << endl;
//
//      for( i = _catnames.begin(); i != _catnames.end(); i++)
//        cout << (*i).first << " " << (*i).second << endl;
}
