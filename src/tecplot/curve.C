#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "curvedef.h"
#include "extrafun.C"
//=============================================================================
//
//=============================================================================
void curve::merge(curve import){

	NodePtr InNode,temp1;
	int InLength;

	InNode=import.head;
	InLength=import.length;

	double	InX,InY,InZ,InPh,InTotHt;
	for(int i=0;i<InLength;i++){
		InX=	InNode->Xdata;
		InY=	InNode->Ydata;
		InZ=	InNode->Zdata;
		InPh=	InNode->ph;
		InTotHt= InNode->totht;
		
		arrange(InX,InY,InZ,InPh,InTotHt);
		InNode=InNode->next;
	}	
}
//=============================================================================
//  returns 
//=============================================================================

void curve::getextents(double* Xarray ,double*  Yarray,double* Zarray, 
							double* Pharray){
	updatelength();
	for(int i=0;i<2;i++){
		Xarray[i]=xbounds[i];
		Yarray[i]=ybounds[i];
		Zarray[i]=zbounds[i];
		Pharray[i]=phbounds[i];
	}
}

//=============================================================================
//
//=============================================================================

void curve::arrange(double Xentry,double Yentry,double Zentry,double pht,
							double totheight){
	NodePtr temp1;
	int location=1;
	temp1=head;
	if(head==NULL) added(Xentry,Yentry,Zentry,pht,totheight);
	else{
		while(temp1->Xdata<=Xentry){
			if(temp1->next==NULL){
				location ++;
				break;
			}
			if(temp1->Xdata==Xentry){
			    while(temp1->Ydata<=Yentry &&temp1->Xdata==Xentry){
				if(temp1->next==NULL) break;
				else {
					temp1=temp1->next;	
					location++;
				}
			   }
			   break;
			}
			else{
				location++;
				temp1=temp1->next;
			}
		}
		push(Xentry,Yentry,Zentry,pht,totheight,location);
	}
}

//=============================================================================
//
//=============================================================================
void curve::reorder(void){

	NodePtr carrier,temp1,temp2;
	temp1=head;
	getinfo();
	cout<<endl<<endl;
	int base=0;
	int offset=1;
	int occ=1;
	double Xnum=0;	
	double Ynum=0;	
	double Znum=0;
	double Ph=0;
	double TotHt=0;
	
	if(temp1==NULL)	{
			cout<<endl<<"No memeber to reorder"<<endl;
	}
	else{
		while(temp1!=NULL){
			Xnum=temp1->Xdata;
			Ynum=temp1->Ydata;
			Znum=temp1->Zdata;
			Ph=temp1->ph;
			TotHt=temp1->totht;
			temp2=temp1->next;		// b 0  off 1  occ 1 
			while(temp2!=NULL){
			    if(temp2->Xdata==Xnum){
				push(temp2->Xdata,temp2->Ydata,temp2->Zdata,
					temp2->ph,temp2->totht,base+occ);
				occ++;
				temp2=temp2->next;
				deloc(base+offset+occ);
				updatelength();
				if(temp2==NULL) break;
			    }
			    else{
				temp2=temp2->next;
				offset++;
			    }
		        }	
		      	if(temp1->next==NULL) break;
	 		temp1=temp1->next;
			base=base+occ;
			occ=1;
			offset=1;
		}
	}
}

//=============================================================================
//
//=============================================================================
void curve::push(double Xentry,double Yentry,double Zentry,double Ph,
								double TotHt){
	NodePtr temp1;
	temp1=new node;
	temp1->Xdata=Xentry;
	temp1->Ydata=Yentry;
	temp1->Zdata=Zentry;
	temp1->ph=Ph;
	temp1->totht=TotHt;
	temp1->next=NULL;
	if(head==NULL) {
		head=temp1;
	}
	else{
		temp1->next=head;
		head=temp1;
	}
	updatelength();
}

//=============================================================================
//
//=============================================================================
void curve::push(double Xentry,double Yentry,double Zentry,double Ph,
							double TotHt,int loc){

	updatelength();
	int temp_len;
	int prloc=0;
	getlength(temp_len);
	
	NodePtr temp1;
	temp1=new node;
	temp1->Xdata=Xentry;
	temp1->Ydata=Yentry;
	temp1->Zdata=Zentry;
	temp1->ph=Ph;
	temp1->totht=TotHt;
	temp1->next=NULL;

	if(temp_len+1>=loc &&loc!=0){
	      NodePtr temp2;
	      temp2=head;
	      if(loc==1) push(Xentry,Yentry,Zentry,Ph,TotHt);
	      else if(loc==temp_len+1){
		 for(temp_len=1;temp_len<loc-1;temp_len++) temp2=temp2->next;
		 temp2->next=temp1;
	      }
	      else{
	       	 for(temp_len=1;temp_len<loc-1;temp_len++) temp2=temp2->next;
			temp1->next=temp2->next;
			temp2->next=temp1;
		}
	}
	else cout<<endl<<"ERROR: location must be NON Zero less than "<<temp_len+1;

}

//=============================================================================
//
//=============================================================================
void curve::deloc(int loc){

	updatelength();
	int temp_len;
	int prloc=0;
	getlength(temp_len);
	
	NodePtr temp1,temp2;
	temp1=new node;

	if(temp_len>=loc && loc!=0){
		temp1=head;
		if(loc==1) head=head->next;
		else if (loc==2) head->next=head->next->next;
		else if (loc==temp_len){
			for( int count=1;count<loc-1;count++) temp1=temp1->next;
			temp1->next=NULL;
		}
		else{
			for(int count=1;count<loc-1;count++){
				temp1=temp1->next;
				temp2=temp1->next;
			}
			temp1->next=temp2->next;
		}
	}
	else 
	   cout<<endl<<"ERROR: location must be NON ZERO less than "<<temp_len;

	getlength(temp_len);

}
//=============================================================================
//
//=============================================================================
void curve::remove(double Xentry,double Yentry,double Zentry){

	NodePtr temp1,temp2;
	temp1=head;
	int temp_len;

	if(temp1==NULL)cout<<endl<<"List has no entry to delete";
	if((temp1->Xdata==Xentry&&temp1->Ydata==Yentry)||temp1!=NULL){
		while(temp1->Xdata==Xentry&&temp1->Ydata==Yentry){
			head=temp1->next;
			temp1=head;
			updatelength();
		}
	}

	while(temp1!=NULL){
		if(temp1->next!=NULL) temp2=temp1->next;
		else break;
		while(temp2->Xdata==Xentry&&temp2->Ydata==Yentry){
			temp1->next=temp2->next;
			if(temp1->next!=NULL) temp2=temp1->next;
			else break;
			updatelength();
		}
		temp1=temp1->next;	
	}

	updatelength();
	getinfo(temp_len);
	cout<<endl<<"Number of members = "<<temp_len;
}

//=============================================================================
//
//=============================================================================
void curve::updatelength(void){		//internal function to update length

	length=0;
	node* temp1;
	temp1=head;
	double xmin=9.9e9,ymin=9.9e9,zmin=9.9e9,phmin=9.9e9;
	double xmax= -9.9e9,ymax= -9.9e9,zmax= -9.9e9,phmax= -9.9e9;
	while(temp1->next!=NULL){
		if(temp1->Xdata<xmin) xmin=temp1->Xdata;
		else if(temp1->Xdata>xmax) xmax=temp1->Xdata;

		if(temp1->Ydata<ymin) ymin=temp1->Ydata;
		else if(temp1->Ydata>ymax) ymax=temp1->Ydata;

		if(temp1->Zdata<zmin) zmin=temp1->Zdata;
		else if(temp1->Zdata>zmax) zmax=temp1->Zdata;

		if(temp1->ph<phmin) phmin=temp1->ph;
		else if(temp1->ph>phmax) phmax=temp1->ph;
		
		length++;
		temp1=temp1->next;
	}
	xbounds[0]=xmin;
	xbounds[1]=xmax;
	ybounds[0]=ymin;
	ybounds[1]=ymax;
	zbounds[0]=zmin;
	zbounds[1]=zmax;
	phbounds[0]=phmin;
	phbounds[1]=phmax;
}

//=============================================================================
//function to add data dat to curve
//=============================================================================
void curve::added(double datX,double datY,double datZ,double pht,
							double totheight){
	node	*temp;
	node 	*temp1;
	if(head==NULL){
		head=new node;
		head->Xdata=datX;
		head->Ydata=datY;
		head->Zdata=datZ;
		head->ph=pht;
		head->totht=totheight;
		head->next=NULL;
	}
	else{
		temp1=head;	
		while(temp1->next!=NULL){
			temp1=temp1->next;
		}
		temp=new node;
		temp->next=NULL;
		temp->Xdata=datX;
		temp->Ydata=datY;
		temp->Zdata=datZ;
		temp->ph=pht;
		temp->totht=totheight;
		temp1->next=temp;
	}
	updatelength();
}

//=============================================================================
//
//=============================================================================
void curve::getinfo(int &  leng){
	node* temp1;
	leng=0;
	length=0;
	temp1=head;
	if(head==NULL) cout<<endl<<"No elements in the list";
	else{
		leng=length=1;
		while(temp1->next!=NULL){
			length++;
			leng= length;
			temp1=temp1->next;
		}
		
	}
}

//=============================================================================
//
//=============================================================================
void curve::getinfo(void ){
	node* temp1;
	length=0;
	temp1=head;
	if(head==NULL) cout<<endl<<"No elements in the list";
	else{
		length=1;
		while(temp1->next!=NULL){
			length++;
			temp1=temp1->next;
		}
	}
}



//=============================================================================
//
//=============================================================================
void curve::getlength(int &  leng){
	node* temp1;
	leng=0;
	length=0;
	temp1=head;
	if(head==NULL) cout<<endl<<"No elements in the list";
	else{
		leng=length=1;
		while(temp1->next!=NULL){
			length++;
			leng= length;
			temp1=temp1->next;
		}
		
	}
}

//=============================================================================
//
//=============================================================================
void curve::writeout(void){
	FILE* out;
	out=fopen("precessed.txt","a");
	NodePtr temp1;
	temp1=head;
	int count=1;

	fprintf(out,"X Range : %lf,%lf\n",xbounds[0],xbounds[1]);
	fprintf(out,"Y Range : %lf,%lf\n",ybounds[0],ybounds[1]);
	fprintf(out,"Z Range : %lf,%lf\n",zbounds[0],zbounds[1]);
	fprintf(out,"Pile Height Range : %lf,%lf\n",phbounds[0],phbounds[1]);
	fprintf(out,"Number of locations : %d\n",length);
	while(temp1!=NULL){
		fprintf(out,"%d\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",count,
		temp1->Xdata,temp1->Ydata,temp1->Zdata,temp1->ph,temp1->totht);
		temp1=temp1->next;
		count++;
	}
	fprintf(out,"\n");
	fclose(out);
}

