#define  GW 100
//============================================================================
//  	returns class curve containting locations having pile height 
//	within 'error'(say 0.01) limits from 'entry'(say 5)
//============================================================================
curve getnearest(curve & megalist,double entry,double error){
	NodePtr temp1,finallist;
	curve finalcurve;
	temp1=megalist.head;
	while(temp1->next!=NULL){
		if(fabs(temp1->ph-entry)<=error){
			finalcurve.arrange(temp1->Xdata,temp1->Ydata,
			          temp1->Zdata,temp1->ph,temp1->totht);
		}
		temp1=temp1->next;
	}
	return finalcurve;
}

//============================================================================
//	returns distance between points (X1,Y1) and (X2,Y2)
//============================================================================
double GetDistance(double X1,double Y1,double X2,double Y2){				//returns distance betwen two points

	double distance;

	distance=sqrt(pow((X2-X1),2)+pow((Y2-Y1),2));
	return(distance);
}


//============================================================================
//	Determines pile height at even contour grid locations by
//      averaging points within proximity radius=10(hard coded)
//      from the grid node coordinates
//============================================================================
void SetPile(curve & megalist,node EvGrid[][GW],int XC,int YC){

	NodePtr temp1;
	double XBoxMax,XBoxMin;
	double YBoxMax,YBoxMin;
	if(XC==0){
		if(YC==0){
			 XBoxMin=EvGrid[XC][YC].Xdata;
			 XBoxMax=EvGrid[XC+1][YC].Xdata;
			 YBoxMin=EvGrid[XC][YC].Ydata;
			 YBoxMax=EvGrid[XC][YC+1].Ydata;
		}
		else if(YC==GW-1){
			 XBoxMin=EvGrid[XC][YC-1].Xdata;
			 XBoxMax=EvGrid[XC+1][YC-1].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC+1][YC].Ydata;
		}
		else{
		 XBoxMin=EvGrid[XC][YC].Xdata;
			 XBoxMax=EvGrid[XC+1][YC].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC][YC+1].Ydata;
		}
	}
	else if(XC==GW-1){

		if(YC==0){
			 XBoxMin=EvGrid[XC-1][YC].Xdata;
			 XBoxMax=EvGrid[XC][YC].Xdata;
			 YBoxMin=EvGrid[XC][YC].Ydata;
			 YBoxMax=EvGrid[XC][YC+1].Ydata;
		}
		else if(YC==GW-1){
			 XBoxMin=EvGrid[XC-1][YC-1].Xdata;
			 XBoxMax=EvGrid[XC][YC-1].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC][YC].Ydata;
		}
		else{
			 XBoxMin=EvGrid[XC-1][YC].Xdata;
			 XBoxMax=EvGrid[XC][YC].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC][YC+1].Ydata;
		}
	}
	else{
		if(YC==0){
			 XBoxMin=EvGrid[XC-1][YC].Xdata;
			 XBoxMax=EvGrid[XC+1][YC].Xdata;
			 YBoxMin=EvGrid[XC-1][YC].Ydata;
			 YBoxMax=EvGrid[XC+1][YC+1].Ydata;
		}
		else if(YC==GW-1){
			 XBoxMin=EvGrid[XC-1][YC-1].Xdata;
			 XBoxMax=EvGrid[XC+1][YC-1].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC][YC].Ydata;
		}
		else{
			 XBoxMin=EvGrid[XC-1][YC].Xdata;
			 XBoxMax=EvGrid[XC+1][YC].Xdata;
			 YBoxMin=EvGrid[XC][YC-1].Ydata;
			 YBoxMax=EvGrid[XC][YC+1].Ydata;
		}

	}
			
	double X1,Y1,Ph=0.0;
	X1=EvGrid[XC][YC].Xdata;
	Y1=EvGrid[XC][YC].Ydata;

	temp1= megalist.head;
	int count=0;
	double Xcurr,Ycurr,distance;
	while(temp1->Xdata<XBoxMin) temp1=temp1->next;
	while(temp1->Xdata<=XBoxMax){
		while(temp1->Ydata<YBoxMin || temp1->Ydata>YBoxMax){
				temp1=temp1->next;
				if(temp1->Xdata>=XBoxMax||
					temp1->Xdata<=XBoxMin ) break;

		}
		if(temp1->Ydata>=YBoxMin&&temp1->Ydata<=YBoxMax){
			Xcurr=temp1->Xdata;
			Ycurr=temp1->Ydata;

			distance=GetDistance(X1,Y1,Xcurr,Ycurr);
			if(distance<=10 ){
				count++;
				//double weight=10.0-distance;
				double weight=1.0;
				Ph=Ph+weight*temp1->ph;
			}
		}
	temp1=temp1->next;
	}
	//cout<<endl<<count;
	Ph=Ph/count;
	if(count==0)Ph=0.0;
	EvGrid[XC][YC].ph=Ph;
}

