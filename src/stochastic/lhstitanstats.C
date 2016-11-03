#include <stdlib.h>
#include <stdio.h>
#include "lhslib.h"
#include <math.h>

int main(){
  char filename[128];
  FILE *fpin, *fpout;
  int iref=0, refnum, isample, Nsample=0, i,j;
  double timereached, Vave, hmax, xcen, ycen, xspread, yspread;
  double probreached=0.0, statmoments[6][3];
  double stats[3];

  for(i=0;i<6; i++)
    for(j=0;j<3; j++)
      statmoments[i][j]=0.0;

  sprintf(filename,"statout_lhs.%02d",iref);
  fpin=fopen(filename,"r");
  fpout=fopen("statout.plot","w");
  fprintf(fpout,"#Nsamp prob_reached    Vave_mean Vave_std_dev    Vave_skew    hmax_mean hmax_std_dev    hmax_skew    xcen_mean xcen_std_dev    xcen_skew    ycen_mean ycen_std_dev    ycen_skew xspread_mean xspr_std_dev    xspr_skew yspread_mean yspr_std_dev    yspr_skew\n");

  while(fpin!=NULL){
    //printf("%s\n",filename);

    do{
      fscanf(fpin,"%2d%6d%lf%lf%lf%lf%lf%lf%lf\n",
	     &refnum,&isample,&timereached,&Vave,&hmax,&xcen,&ycen,&xspread,&yspread);

      //printf("isample=%d\n",isample);
      if(timereached>=0.0) probreached+=1.0;

      Nsample++;
      for(i=0; i<3; i++){
	statmoments[0][i]+=pow(Vave,i+1.0);
	statmoments[1][i]+=pow(hmax,i+1.0);
	statmoments[2][i]+=pow(xcen,i+1.0);
	statmoments[3][i]+=pow(ycen,i+1.0);
	statmoments[4][i]+=pow(xspread,i+1.0);
	statmoments[5][i]+=pow(yspread,i+1.0);
      }
      
    }while(!feof(fpin));

    /*      fscanf(fpin,"%2d%6d%lf%lf%lf%lf%lf%lf%lf\n",
	    &refnum,&isample,&timereached,&Vave,&hmax,&xcen,&ycen,&xspread,&yspread);
	    }
    */
    
    fclose(fpin);

    fprintf(fpout,"%6d %12.6g",Nsample,probreached/Nsample);
    for(i=0;i<6;i++){
      for(j=0;j<3;j++) stats[j]=statmoments[i][j]/Nsample;
      stats[2]=stats[2]-3*stats[1]*stats[0]+2*stats[0]*stats[0]*stats[0];
      stats[1]=sqrt(stats[1]-stats[0]*stats[0]);
      stats[2]/=stats[1]*stats[1]*stats[1];
      fprintf(fpout," %12.6g %12.6g %12.6g",stats[0],stats[1],stats[2]);
    }
    fprintf(fpout,"\n");

    iref++;
    sprintf(filename,"statout_lhs.%02d",iref);
    fpin=fopen(filename,"r");
  }

  fclose(fpout);



  return 0;
}
