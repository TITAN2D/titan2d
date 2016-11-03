#include <stdlib.h>
#include <stdio.h>
#include "lhslib.h"

int main(){
  printf("This program generates the latin hypercube sampling (lhs) sample points for the VOLUME SCALE FACTOR for Titan2D.  The resultant volume for a particular run is the product of the scale factor and the base volume entered through the Titan2D python gui\n\n");

  printf("The volume scale factor will be uniformly distributed\n");

  int Nstart;
  printf("Enter the number of initial divisions: ");
  scanf("%d",&Nstart);
  
  int Nref;
  printf("Enter the number of refinements (subdividing each interval into 2): ");
  scanf("%d",&Nref);

  int *refnum, isample, Nsample;
  double *samples, doubleswap;

  Nsample=uniformlhs(Nstart, Nref, &refnum, &samples);
  double minvol, maxvol;
  printf("Enter the min and max volume scale factor: ");
  scanf("%lf%lf",&minvol,&maxvol);
  if(minvol>maxvol){
    doubleswap=maxvol;
    maxvol=minvol;
    minvol=doubleswap;
  }
  doubleswap=maxvol-minvol;
  for(isample=0;isample<Nsample;isample++)
    samples[isample]=minvol+doubleswap*samples[isample];
  
  FILE *fp=fopen("stat_ctl.vol","w");
  for(isample=0;isample<Nsample;isample++)
    fprintf(fp,"%4d  %6d  %14.8g\n",refnum[isample],isample,samples[isample]);
  fclose(fp);

  CDeAllocI1(refnum);
  CDeAllocD1(samples);

  return 0;
}
