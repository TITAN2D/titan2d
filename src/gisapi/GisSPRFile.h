#ifndef GisSPRFile_H
#define GisSPRFile_H

#include <sstream>
#include <vector>
#include "GisAscFile.h"
#include "GisLines.h"

class GisSPRFile :	public GisAscFile{
public:

	GisSPRFile(const string& name, const char* mode = "r");
	
	virtual ~GisSPRFile(){} 

	bool gotoPOINTSSection();
	bool gotoLINESSection();

	bool readLabels(vector<double>& x, vector<double>& y, vector<string>& labelStr);

	bool readFirstLine(vector<double>& x, vector<double>& y);

	bool readNextLine(vector<double>& x, vector<double>& y);

	bool readINFOSection();

	bool gotoSection(string& sectionName);


protected:

	string _sepStr;

private:
	
// No copy allowed
	GisSPRFile(const GisSPRFile&);
	GisSPRFile& operator=(const GisSPRFile&);
};

#endif
