#include <stdlib.h>
#include <stdio.h>
#include "lhslib.h"

int main(){
  printf("This program generates the latin hypercube sampling (lhs) sample output from Titan2D.\n");

  int Nstart;
  printf("Enter the number of initial divisions: ");
  scanf("%d",&Nstart);
  
  int Nref;
  printf("Enter the number of refinements (subdividing each interval into 2): ");
  scanf("%d",&Nref);

  int *refnum, isample, Nsample;
  double *samples, doubleswap;

  Nsample=normallhs(Nstart, Nref, &refnum, &samples);

  FILE *fp;
  char filename[128];
  double time, Vave, hmax, xcen, ycen, xspread, yspread;
  for(isample=0;isample<Nsample;isample++){
    sprintf(filename,"statout_lhs.%02d",refnum[isample]);
    fp=fopen(filename,"a");

    time=30.0+15.0*samples[isample];    if(time>60.0||time<0.0) time=-1.0;
    Vave=15.0+2.0*samples[isample];     if(Vave<0.0) Vave=0.0;
    hmax=10.0+1.0*samples[isample];     if(hmax<0.0) hmax=0.0;
    xcen=5.5e5+250.0*samples[isample];  
    ycen=1.25e6+300.0*samples[isample];
    xspread=600+100*samples[isample];   if(xspread<10.0) xspread=10.0; 
    yspread=700+100*samples[isample];   if(yspread<15.0) yspread=15.0;

    fprintf(fp,"%2d %6d %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
	    refnum[isample],isample,time,Vave,hmax,xcen,ycen,xspread,yspread);

    fclose(fp);
  }



  CDeAllocI1(refnum);
  CDeAllocD1(samples);

  return 0;
}
