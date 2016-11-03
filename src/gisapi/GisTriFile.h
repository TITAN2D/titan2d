#ifndef GisTriFile_H 
#define GisTriFile_H 
 
#include <vector> 
#include "GisBinFile.h" 
 
class GisTriOut 
{ 
public: 
  
  GisTriOut(){} 
  
  virtual ~GisTriOut(){}  
  
  double x0_; double y0_; float z0_; 
  double x1_; double y1_; float z1_; 
  double x2_; double y2_; float z2_; 
  float pHeight_; 
  float xmom_; 
  float ymom_; 
  
  int key1_; int key2_; 
  int gen_; 
  int istriedge_; 
  
  double xMin(); 
  double xMax(); 
  double yMin(); 
  double yMax(); 
  
  bool contains(double x, double y); 
  void print(); 
  
 protected: 
  
 private: 
  
}; 

class GisTriFile :	public GisBinFile 
{ 
 public: 
  
  GisTriFile(const string& name, const char* mode = "r"); 
  
  virtual ~GisTriFile(){}  
  
  bool readTimeStepInfo(); 
  
  int versionNumber(); 
  //	return -1, if error. 
  
  bool readTriData(GisTriOut& triOut); 
  
  int numTri_; 
  int timeStep_; 
  float simTime_; 
  float pileMin_; 
  float pileMax_;  
  float xMomMin_; 
  float xMomMax_; 
  float yMomMin_; 
  float yMomMax_; 
  double xMin_; 
  double xMax_; 
  double yMin_; 
  double yMax_; 
  float elevMin_; 
  float elevMax_; 
  float maxVel_; 
  
 protected: 
  
 private: 
  
  // No copy allowed 
  GisTriFile(const GisTriFile&); 
  GisTriFile& operator=(const GisTriFile&); 
}; 

#endif 
