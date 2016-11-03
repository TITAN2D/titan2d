//=============================================================================
// New algorithm being developed for tracing the contour 
// Brief description:
//	1) loop to find first cell qualified for drawing contour,
// 	   then step out of x,y loop
//	2) with first cell loop in triangles to find qualified triangle 
//	3) find line segment in this triangle , mark I,J and edge in global
//=============================================================================
#include "contour.h"
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

//=============================================================================
// Load Data function: gulps in all the info in contour type object
//=============================================================================
void Contour::LoadData(double PhArray[][GW], node  EvGrid[][GW],double  heightlvl[4], FILE* cellfile){

  for(int i=0;i<GW;i++){
    for(int j=0;j<GW;j++){
      this->PhArray[i][j]=PhArray[i][j];
      this->EvGrid[i][j]=EvGrid[i][j];
    }
  }
  for(int ii=0;ii<4;ii++) this->heightlvl[ii]=heightlvl[ii];
  this->cellfile=cellfile;

}

void Contour::FlushTouchPad(void){
	for(int i=0;i<GW;i++){
    		for(int j=0;j<GW;j++){
			touchpad[i][j]=0;
    		}
  	}

}

void Contour::TraceCurveFromHere(int ContNm){
    int firstPassFlag=1;
    int thisTriangle=0;	// current triangle
    int thisEdge=0;	// current edge	
    int startEdge,termEdge; // starting edge and termination edge of line segment in thisEdge
    int thisXmarch,thisYmarch;
    thisXmarch=startXIndex;
    thisYmarch=startYIndex;
    int CellChangeFlag=0;

    //starting location for this contour
    int XStart,YStart;
    int StartTriangle,StartEdge;
    int itercounter=0;
    int backToStartFlag=0;
    if(firstPassFlag==1){
	int triCount=0;
	while(triCount<3&&backToStartFlag==0){
		CellChangeFlag=0;
      		FormTriangleWithEdge(thisTriangle);
      		int LineVCount=0;
      		double  ConLineX[2];
      		double  ConLineY[2];
      		int edgeCount=0;
		do{
			double Xe1,Xe2,Ye1,Ye2;		
			double Phe1,Phe2;
			double m,n,X_int,Y_int;
			Xe1=SubTriX[thisEdge];
			Ye1=SubTriY[thisEdge];
			Xe2=SubTriX[thisEdge+1];
			Ye2=SubTriY[thisEdge+1];
			Phe1=SubTriPh[thisEdge];
			Phe2=SubTriPh[thisEdge+1];
			if(X_int==startX&&Y_int==startY&&firstPassFlag==0){
				backToStartFlag=1;
				fprintf(cellfile, "\n%f\t%f\t%f\n%f\t%f\t%lf\n",
				float(X_int),float(Y_int),float(Phe1),
				float(startX),float(startY),float(Phe2));
				fflush(cellfile);
				break;
			}
			if(MAX(Phe2,Phe1)>heightlvl[ContNm] && MIN(Phe2,Phe1)<=heightlvl[ContNm]){
				m=(fabs(heightlvl[ContNm]-Phe1)/fabs(Phe1-Phe2));
				n=1.0-m;
				X_int=(m*Xe2+n*Xe1)/(m+n);
				Y_int=(m*Ye2+n*Ye1)/(m+n);
				ConLineX[LineVCount]=X_int;
				ConLineY[LineVCount]=Y_int;
				LineVCount++;
				if(LineVCount==2){
					termEdge=thisEdge;
					fprintf(cellfile, "\n%f\t%f\t%f",
					float(ConLineX[0]),float(ConLineY[0]),float(heightlvl[ContNm]));
					fflush(cellfile);
					ConLineX[0]=ConLineX[1];	
					ConLineY[0]=ConLineY[1];
					int boundry_flag=0;
					if(startXIndex==0) boundry_flag=1;
					if(startYIndex==0) boundry_flag=2;
					CellChangeFlag=FollowCurve(ContNm,thisTriangle, startEdge,termEdge,thisEdge);
					triCount= -1;
					thisTriangle=GetNextTriangle(thisTriangle);
					thisTriangle=GetNextTriangle(thisTriangle);
					thisTriangle=GetNextTriangle(thisTriangle);
					firstPassFlag=0;
					break;
				}//IF 
				else{
					if(firstPassFlag==1){
						startX=X_int;
						startY=Y_int;
					}
				        startEdge=thisEdge;
					edgeCount++;
					thisEdge=GetNextEdge(thisEdge);
				}
			}// If maxmin
			else{
			        edgeCount++;
				thisEdge=GetNextEdge(thisEdge);
			}
		}while(edgeCount<3);
		triCount++;
		thisTriangle=GetNextTriangle(thisTriangle);
		itercounter++;
	}//tri counter
    }//first Pss Flag
   
 // }//for loop contours
  //fclose(cellfile);
 

} 

//=============================================================================
// Loads cell at x and y location index in Cell" " datavariable of class
//=============================================================================

void Contour::LoadCellAt(int xmarch,int ymarch){
  MaxPh=0e-10;
  MinPh=9.9e100;

  int CurNum=xmarch*xcells+ymarch; 
  // Pile Height at 4 corners of cell
  cell.CellHt[0][0]= PhArray[xmarch][ymarch];		 		
  cell.CellHt[1][0]= PhArray[xmarch+1][ymarch];
  cell.CellHt[1][1]= PhArray[xmarch+1][ymarch+1];
  cell.CellHt[0][1]= PhArray[xmarch][ymarch+1];		 		

  //X locations at 4 corners
  cell.CellX[0][0]= EvGrid[xmarch][ymarch].Xdata;
  cell.CellX[1][0]= EvGrid[xmarch+1][ymarch].Xdata;
  cell.CellX[1][1]= EvGrid[xmarch+1][ymarch+1].Xdata;
  cell.CellX[0][1]= EvGrid[xmarch][ymarch+1].Xdata;

  //Y  locations of 4 corners
  cell.CellY[0][0]= EvGrid[xmarch][ymarch].Ydata;
  cell.CellY[1][0]= EvGrid[xmarch+1][ymarch].Ydata;
  cell.CellY[1][1]= EvGrid[xmarch+1][ymarch+1].Ydata;
  cell.CellY[0][1]= EvGrid[xmarch][ymarch+1].Ydata;

  // finding maximum and minimum pile height for current cell
  for(int vertex=0;vertex<4;vertex++){
    if(MaxPh<=cell.CellHt[InI[vertex]][InJ[vertex]]) 
      MaxPh=cell.CellHt[InI[vertex]][InJ[vertex]];
    if(MinPh>=cell.CellHt[InI[vertex]][InJ[vertex]]) 
      MinPh=cell.CellHt[InI[vertex]][InJ[vertex]];
  }

  double sumX=0.0,sumY=0.0,sumPh=0.0;
  //adding X, Y values of corners and dividing by 4 to get center
  for(int vertex=0;vertex<4;vertex++){
    sumX=sumX+cell.CellX[InI[vertex]][InJ[vertex]];
    sumY=sumY+cell.CellY[InI[vertex]][InJ[vertex]];
    sumPh=sumPh+cell.CellHt[InI[vertex]][InJ[vertex]];
  }
  cell.CellCenX=sumX/4.0;
  cell.CellCenY=sumY/4.0;
  cell.CellCenPh=sumPh/4.0;
  XCurrent=xmarch;
  YCurrent=ymarch;
}

//=============================================================================
// returnds 1 if contour of value heightlvl[ContNm] passes thr loaded cell 
//=============================================================================
int Contour::CheckIntersection(int ContNm){
  if((MaxPh>=heightlvl[ContNm]&&MinPh<heightlvl[ContNm])){
    return 1;
  }
  else return 0;
}


//=============================================================================
// Forms triangle using edge number supplied
//=============================================================================
void Contour::FormTriangleWithEdge(int edge){
  SubTriX[0]=cell.CellX[InI[edge]][InJ[edge]];
  SubTriX[1]=cell.CellX[InI[edge+1]][InJ[edge+1]];
  SubTriX[2]=cell.CellCenX;
  SubTriX[3]=cell.CellX[InI[edge]][InJ[edge]];

  SubTriY[0]=cell.CellY[InI[edge]][InJ[edge]];
  SubTriY[1]=cell.CellY[InI[edge+1]][InJ[edge+1]];
  SubTriY[2]=cell.CellCenY;
  SubTriY[3]=cell.CellY[InI[edge]][InJ[edge]];

  SubTriPh[0]=cell.CellHt[InI[edge]][InJ[edge]];
  SubTriPh[1]=cell.CellHt[InI[edge+1]][InJ[edge+1]];
  SubTriPh[2]=cell.CellCenPh;
  SubTriPh[3]=cell.CellHt[InI[edge]][InJ[edge]];
}

int Contour::FollowCurve(int ContNm,int & thisTriangle, int & startEdge,int & termEdge,int & thisEdge){

touchpad[XCurrent][YCurrent]=1;
switch(thisTriangle){
  case 0:
    switch(startEdge){
    	case 0: if(termEdge==1){ 
      			LoadCellAt(XCurrent,YCurrent);
	      		thisTriangle=1;
      			thisEdge=2;
      			return 0;
    		}
    		else if(termEdge==2){
      			LoadCellAt(XCurrent,YCurrent);
      			thisTriangle=3;
      			thisEdge=1;
      			return 0;
    		}
      		break;
    	case 1:if(termEdge==0){
			if(YCurrent==0){
				LoadCellAt(XCurrent-1,YCurrent);
				touchpad[XCurrent][YCurrent]=1;
				int interFlag=CheckIntersection(ContNm);
				while(interFlag==0){
					XCurrent--;
					touchpad[XCurrent][YCurrent]=1;
					LoadCellAt(XCurrent,YCurrent);
					interFlag=CheckIntersection(ContNm);
				}
	    			thisTriangle=0;
    				thisEdge=0;
			}
			else{
      				LoadCellAt(XCurrent,YCurrent-1);
	      			thisTriangle=2;
	      			thisEdge=0;
      				return 1;
			}
    		}
    		else if(termEdge==2){
      			LoadCellAt(XCurrent,YCurrent);
      			thisTriangle=3;
      			thisEdge=1;
      			return 0;
    		}
      		break;
    	case 2 :if(termEdge==0){
			if(YCurrent==0){
				LoadCellAt(XCurrent-1,YCurrent);
				touchpad[XCurrent][YCurrent]=1;
				int interFlag=CheckIntersection(ContNm);
				while(interFlag==0){
					XCurrent--;
					touchpad[XCurrent][YCurrent]=1;
					LoadCellAt(XCurrent,YCurrent);
					interFlag=CheckIntersection(ContNm);
				}
	    			thisTriangle=3;
    				thisEdge=0;
			}
			else{
      				LoadCellAt(XCurrent,YCurrent-1);
	      			thisTriangle=2;
	      			thisEdge=0;
	      			return 1;
			}
    		}
    		else if(termEdge==1){
    	  		LoadCellAt(XCurrent,YCurrent);
      			thisTriangle=1;
      			thisEdge=2;
      			return 0;
    		}
      		break;
    	};
    	break;
  case 1:switch(startEdge){
  		case 0: if(termEdge==1){ 
  	  			LoadCellAt(XCurrent,YCurrent);
    				thisTriangle=2;
    				thisEdge=2;
    				return 0;
  			}
  			else if(termEdge==2){
	    			LoadCellAt(XCurrent,YCurrent);
	    			thisTriangle=0;
	    			thisEdge=1;
    				return 0;
	  		}
	    		break;
		  case 1:if(termEdge==0){
			    	LoadCellAt(XCurrent+1,YCurrent);
				thisTriangle=3;
    				thisEdge=0;
    				return 1;
  			}
  			else if(termEdge==2){
    				LoadCellAt(XCurrent,YCurrent);
    				thisTriangle=0;
    				thisEdge=1;
    				return 0;
  			}
    			break;
  		case 2 :if(termEdge==0){
    				LoadCellAt(XCurrent+1,YCurrent);
    				thisTriangle=3;
   	 			thisEdge=0;
    				return 1;
  			}
  			else if(termEdge==1){
    				LoadCellAt(XCurrent,YCurrent);
   	 			thisTriangle=2;
    				thisEdge=2;
    				return 0;
  			}
   			break;
  		};
    		break;
 case 2: switch(startEdge){
	 case 0: if(termEdge==1){ 
		   LoadCellAt(XCurrent,YCurrent);
		   thisTriangle=3;
		   thisEdge=2;
		   return 0;
		 }
		 else if(termEdge==2){
		    LoadCellAt(XCurrent,YCurrent);
		    thisTriangle=1;
		    thisEdge=1;
		    return 0;
  		}
  		break;
  	case 1:if(termEdge==0){
	 	 	LoadCellAt(XCurrent,YCurrent+1);
			thisTriangle=0;
			thisEdge=0;
    			return 1;
  		}
  		else if(termEdge==2){
    			LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=1;
    			thisEdge=1;
    			return 0;
  		}
    		break;
  	case 2 :if(termEdge==0){
    			LoadCellAt(XCurrent,YCurrent+1);
    			thisTriangle=0;
    			thisEdge=0;
    			return 1;
  		}
  		else if(termEdge==1){
    			LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=3;
    			thisEdge=2;
    			return 0;
  		}
    		break;
	 };
    	 break;
  case 3:switch(startEdge){
  	case 0: if(termEdge==1){ 
    			LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=0;
    			thisEdge=2;
    			return 0;
  		}
  		else if(termEdge==2){
    			LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=2;
    			thisEdge=1;
    			return 0;
  		}
    		break;
  	case 1:if(termEdge==0){
			if(XCurrent==0){
  		  		LoadCellAt(XCurrent,YCurrent+1);
				touchpad[XCurrent][YCurrent]=1;
				int interFlag=CheckIntersection(ContNm);	
				while(interFlag==0){
					YCurrent++;
					touchpad[XCurrent][YCurrent]=1;
					LoadCellAt(XCurrent,YCurrent);
					interFlag=CheckIntersection(ContNm);
				}
	    			thisTriangle=3;
    				thisEdge=0;
			}
			else{
  		  		LoadCellAt(XCurrent-1,YCurrent);
	    			thisTriangle=1;
    				thisEdge=0;
			}
    			return 1;
  		}
  		else if(termEdge==2){
    			LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=2;
    			thisEdge=1;
    			return 0;
  		}
    		break;
  	case 2 :if(termEdge==0){
	                if(XCurrent==0){
  		  		LoadCellAt(XCurrent,YCurrent+1);
	    			thisTriangle=3;
    				thisEdge=0;
                        }
                        else{
			  LoadCellAt(XCurrent-1,YCurrent);
			  thisTriangle=1;
			  thisEdge=0;
			}
    			return 1;
 	 	}
  		else if(termEdge==1){
   	 		LoadCellAt(XCurrent,YCurrent);
    			thisTriangle=0;
    			thisEdge=2;
    			return 0;
  		}
    		break;
  	};
    	break;
  }
	
}

int Contour::GetNextTriangle(int thisTria){
  switch(thisTria){
  case 0:
    return 1;
    break;
  case 1:
    return 2;
    break;
  case 2:
    return 3;
    break;
  case 3:
    return 0;
    break;
  };
}

int Contour::GetNextEdge(int thisEdge){
  switch(thisEdge){
  case 0:
    return 1;
    break;
  case 1:
    return 2;
    break;
  case 2:
    return 0;
    break;
  };
}



//=============================================================================
// Contouring function : most imp function
//=============================================================================

void Contour::contfun(void){
  // iterating between 4 contour levels
	for(int ContNm=0;ContNm<4;ContNm++){
		FlushTouchPad();
    	// location where tracing of contour starts
    		startXIndex=-1;
    		startYIndex=-1;
    		startEdge=0;
    		MaxPh=0e-10;
    		MinPh=9.9e100;
    		int xmarch,ymarch;
    		while(startXIndex<0&&startYIndex<0){ 
			//iterating from X min to Xmax
		      	for(xmarch=0;xmarch<GW-1;xmarch++){
				// iterating from ymin to ymax
				for(ymarch=0;ymarch<GW-1;ymarch++){
		  			// counting the current cell number
		  			//int CurNum=xmarch*xcells+ymarch; 
		 			LoadCellAt(xmarch,ymarch);
		  			int interFlag=CheckIntersection(ContNm);	
		  			if(interFlag&&touchpad[xmarch][ymarch]==0){
						if(xmarch==0){
						  cout<< "contour "<<ContNm<<" intersecting x=0 at ("<<xmarch<<","<<ymarch<<")"<<endl;
						  continue;
						}
		    				startXIndex= xmarch;
		    				startYIndex= ymarch;
						TraceCurveFromHere(ContNm);
						touchpad[xmarch][ymarch]=1;
		  			}
					else continue;
						
				}//for loop ymarch
		      	}// for loop xmarch
		      	break;
	    	}//while loop
        // finding center of the cell and averaging pile height at it
	}
}
