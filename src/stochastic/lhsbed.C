#include <stdlib.h>
#include <stdio.h>
#include "lhslib.h"

int main(){
  printf("This program generates the latin hypercube sampling (lhs) sample points for the bed friction angle for a single material to be used in Titan2D.\n");

  int distributiontype;
  printf("Please choose between a\n1) uniform distribution\n2) normal (gaussian) distribution\nenter 1 or 2: ");
  scanf("%d",&distributiontype);

  int Nstart;
  printf("Enter the number of initial divisions: ");
  scanf("%d",&Nstart);
  
  int Nref;
  printf("Enter the number of refinements (subdividing each interval into 2): ");
  scanf("%d",&Nref);

  int *refnum, isample, Nsample;
  double *samples, doubleswap;
  switch(distributiontype){
  case 1: //uniform
    Nsample=uniformlhs(Nstart, Nref, &refnum, &samples);
    double minbed, maxbed;
    printf("Enter the min and max bed friction angles in degrees: ");
    scanf("%lf%lf",&minbed,&maxbed);
    if(minbed>maxbed){
      doubleswap=maxbed;
      maxbed=minbed;
      minbed=doubleswap;
    }
    doubleswap=maxbed-minbed;
    for(isample=0;isample<Nsample;isample++)
      samples[isample]=minbed+doubleswap*samples[isample];
    break;
  case 2:
    Nsample=normallhs(Nstart, Nref, &refnum, &samples);
    double meanbed, stddevbed;
    printf("Enter the mean bed friction angle in degrees: ");
    scanf("%lf",&meanbed);

    printf("Enter the standard deviation of the bed friction angle in degrees: ");
    scanf("%lf",&stddevbed);
    for(isample=0;isample<Nsample;isample++)
      samples[isample]=meanbed+stddevbed*samples[isample];
    break;
  default:
    printf("lhsbed.C does not recognize distribution type %d\nAborting!!!",
	   distributiontype);
    exit(1);
    break;
  }
  
  double A;
  printf("Enter the number of degrees \"A\" by which the internal friction angle exceeds the bed friction angle\nintfrict = bedfrict + A\n");
  scanf("%lf",&A);

  FILE *fp=fopen("stat_ctl.bed","w");
  for(isample=0;isample<Nsample;isample++)
    fprintf(fp,"%4d  %6d  %12.6g  %12.6g\n",refnum[isample],isample,samples[isample],samples[isample]+A);
  fclose(fp);

  CDeAllocI1(refnum);
  CDeAllocD1(samples);

  return 0;
}
