#define GW 100
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdio.h>
struct cellinfo{
	double CellHt[2][2];
  	double CellX[2][2];
	double CellY[2][2];
	double CellCenX;	
    	double CellCenY;	
    	double CellCenPh;	 
	short int touched;
};


class Contour{
 	public:
	int GWa;
	int xcells ;
	int ycells ;

  	// Starting Location of Contour
	int startXIndex;	// X grid location where curve started
	int startYIndex;	// Y grid location where curve started
	int startEdge;
	double startX;
	double startY;
	//currently loaded cell
	cellinfo cell;
	int touchpad[GW][GW];
	double 	SubTriX[4],SubTriY[4],SubTriPh[4];

	static int InI[]; //indexing I
	static int InJ[]; //indexing J
	double PhArray[GW][GW];
	node EvGrid[GW][GW];
	double heightlvl[4];
	int XCurrent;
	int YCurrent;

	double MaxPh;
	double MinPh;	
	FILE *cellfile;

// 
	void TraceCurveFromHere(int);

// Loads in Pile Height Array, Linked Lists and Contour Level array 
	void LoadData(double PhArray[][GW],node EvGrid[][GW],double heightlvl[4], FILE* cellfile);
// main contouring routine : all other methods are supplimentary
	void contfun(void);
// Loads Cell at location specified by xIncrement and yIncrement	
	void LoadCellAt(int,int);
// Checks the intersection of loaded cell with contourlevel number supplied
	int CheckIntersection(int);
// Forms triangle using edge number supplied as argument with center of loaded cell
	void FormTriangleWithEdge(int);
// should redirect to next cell and next edge directly
	int FollowCurve(int ,int & ,int &,int &,int &);
// returns the number of next triangle in cell using ring 0-1-2-3-0
	int GetNextTriangle(int);
// returns the number of next edge in cell using ring 0-1-2-0
	int GetNextEdge(int);
// resets touchpad
	void FlushTouchPad(void);
// Constructor
	Contour(int);
};	
// Indexing used in contouring
	int Contour::InI[]={0,1,1,0,0};
	int Contour::InJ[]={0,0,1,1,0};
//	double Contour::MaxPh=0.0;
//	double Contour::MinPh=1.1e10;

Contour::Contour(int a){
	GWa=a;
}

