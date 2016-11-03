#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lhslib.h"

int uniformlhs(int Nstart, int Nref, int **refnum, double **samples){

  int Nsample=Nstart*int(pow(2.0,Nref));
  *refnum=CAllocI1(Nsample);
  *samples=CAllocD1(Nsample);

  double *urv1=CAllocD1(Nstart), *urv2;
  int i, isample=0, iref=0;
  long idum=-1;
  ran1(&idum);
  for(i=0;i<Nstart;i++){
    urv1[i]=(i+ran1(&idum))/Nstart;
    (*refnum)[isample]=iref;
    (*samples)[isample]=urv1[i];
    isample++;
  }

  int Nbatchsize=Nstart, iurv1;
  for(iref=1;iref<=Nref;iref++){
    Nbatchsize*=2;
    urv2=CAllocD1(Nbatchsize);
    for(i=0,iurv1=0;i<Nbatchsize;i+=2,iurv1++,isample++){
      (*refnum)[isample]=iref;
      if(urv1[iurv1]<=(i+1.0)/Nbatchsize){
	urv2[i]=urv1[iurv1];
	(*samples)[isample]=urv2[i+1]=(i+1.0+ran1(&idum))/Nbatchsize;
      }
      else{
	(*samples)[isample]=urv2[i]=(i+ran1(&idum))/Nbatchsize;
	urv2[i+1]=urv1[iurv1];
      }
    }
    CDeAllocD1(urv1);
    urv1=urv2;
  }


#ifdef DEBUGUNI
  for(isample=0;isample<Nsample;isample++)
    printf("%3d  %8.6f  %8.6f\n",isample,(*samples)[isample],urv1[isample]);
#endif
  
  CDeAllocD1(urv1);
  return Nsample;
}

int normallhs(int Nstart, int Nref, int **refnum, double **samples){

  int Nsample=uniformlhs(Nstart,Nref,refnum,samples);
  double A, invsqrt2=1.0/sqrt(2.0), sqrt2pi=sqrt(8.0*atan(1.0));
  double mean=0, stddev=0, skew=0, kurtosis=0;
  int isample, inewt;
  for(isample=0;isample<Nsample;isample++){
    A=(*samples)[isample]-0.5;
    (*samples)[isample]=0;
    for(inewt=0;inewt<10;inewt++)
      (*samples)[isample]+=(A-0.5*erf((*samples)[isample]*invsqrt2))
	*sqrt2pi*exp(0.5*(*samples)[isample]*(*samples)[isample]);
    mean+=(*samples)[isample];
  }

  mean/=Nsample;
  for(isample=0;isample<Nsample;isample++){
    stddev+=pow((*samples)[isample]-mean,2);
    skew+=pow((*samples)[isample]-mean,3);
    kurtosis+=pow((*samples)[isample]-mean,4);
  }

  stddev=sqrt(stddev/Nsample);
  skew/=(Nsample*stddev*stddev*stddev);
  kurtosis/=(Nsample*stddev*stddev*stddev*stddev);
  
  printf("lhs generated standard normal statistics\nmean    =%f\nstddev  =%f\nskewness=%f\nkurtosis=%f\n",
	 mean,stddev,skew,kurtosis);


  return Nsample;
}
