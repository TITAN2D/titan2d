/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: useful_lib.C 131 2007-06-07 19:58:23Z dkumar $ 
 */

/*#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
*/ 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"useful_lib.h"

// start for ran1()
#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS) 
//end for ran1()


/************************************************************/
/************************************************************/
/******** multdimensional memory allocation routines ********/
/************************************************************/
/************************************************************/

/* The naming convention
   the starting "C" indicates it uses C/C++ as opposed to
   FORTRAN indices (they start from zero instead of from 
   1).  "Alloc" stands for memory allocation. alternatively
   "DeAlloc" means memory deallocation.  "I" stands for 
   integer.  Alternatively "C" means character, "F" means
   floating point, "D" means double precision floating 
   point, "U" means unsigned integer.  The trailing integer 
   stands for the dimension of the array, i.e. "1" means 
   it's 1 dimensional, "2" means it's 2 dimensional, etc. */

//should add 4th dimensional allocation routines


//unsigned's
unsigned *CAllocU1(int N1){
  return((unsigned *) malloc(N1*sizeof(unsigned)));
}

unsigned **CAllocU2(int N1, int N2){
  int i;
  unsigned **U2;
  
  U2= (unsigned **) malloc(N1*sizeof(unsigned *));
  *U2= (unsigned *) malloc(N1*N2*sizeof(unsigned));
  for(i=1; i<N1; i++) U2[i]=U2[i-1]+N2;
  
  return(U2);
}

unsigned ***CAllocU3(int N1, int N2, int N3){
  int i,j;
  unsigned ***U3;

  U3=(unsigned ***) malloc(N1*sizeof(unsigned **));
  *U3=(unsigned **) malloc(N1*N2*sizeof(unsigned *));
  **U3=(unsigned *) malloc(N1*N2*N3*sizeof(unsigned));
  
  for(i=0;i<N1;i++){
    if(i!=0){
      U3[i]=U3[i-1]+N2;          
      U3[i][0]=U3[i-1][N2-1]+N3;} 
    for(j=1;j<N2;j++)
      U3[i][j]=U3[i][j-1]+N3;}  

  return(U3);
}

//int's
int *CAllocI1(int N1){
  return((int *) malloc(N1*sizeof(int)));
}

int **CAllocI2(int N1, int N2){
  int i;
  int **I2;
  
  I2= (int **) malloc(N1*sizeof(int *));
  *I2= (int *) malloc(N1*N2*sizeof(int));
  for(i=1; i<N1; i++) I2[i]=I2[i-1]+N2;
  
  return(I2);
}

int ***CAllocI3(int N1, int N2, int N3){
  int i,j;
  int ***I3;

  I3=(int ***) malloc(N1*sizeof(int **));
  *I3=(int **) malloc(N1*N2*sizeof(int *));
  **I3=(int *) malloc(N1*N2*N3*sizeof(int));
  
  for(i=0;i<N1;i++){
    if(i!=0){
      I3[i]=I3[i-1]+N2;          
      I3[i][0]=I3[i-1][N2-1]+N3;} 
    for(j=1;j<N2;j++)
      I3[i][j]=I3[i][j-1]+N3;}  

  return(I3);
}

//float's
float *CAllocF1(int N1){
  return((float *) malloc(N1*sizeof(float)));
}

float **CAllocF2(int N1, int N2){
  int i;
  float **F2;
  
  F2= (float **) malloc(N1*sizeof(float *));
  *F2= (float *) malloc(N1*N2*sizeof(float));
  for(i=1; i<N1; i++) F2[i]=F2[i-1]+N2;
  
  return(F2);
}

float ***CAllocF3(int N1, int N2, int N3){
  int i,j;
  float ***F3;

  F3=(float ***) malloc(N1*sizeof(float **));
  *F3=(float **) malloc(N1*N2*sizeof(float *));
  **F3=(float *) malloc(N1*N2*N3*sizeof(float));
  
  for(i=0;i<N1;i++){
    if(i!=0){
      F3[i]=F3[i-1]+N2;          
      F3[i][0]=F3[i-1][N2-1]+N3;} 
    for(j=1;j<N2;j++)
      F3[i][j]=F3[i][j-1]+N3;}  

  return(F3);
}


//double's
double *CAllocD1(int N1){
  return((double *) malloc(N1*sizeof(double)));
}

double **CAllocD2(int N1, int N2){
  int i;
  double **D2;
  
  D2= (double **) malloc(N1*sizeof(double *));
  *D2= (double *) malloc(N1*N2*sizeof(double));
  for(i=1; i<N1; i++) D2[i]=D2[i-1]+N2;

  return(D2);
}

double ***CAllocD3(int N1, int N2, int N3){
  int i,j;
  double ***D3;

  D3=(double ***) malloc(N1*sizeof(double **));
  *D3=(double **) malloc(N1*N2*sizeof(double *));
  **D3=(double *) malloc(N1*N2*N3*sizeof(double));
  
  for(i=0;i<N1;i++){
    if(i!=0){
      D3[i]=D3[i-1]+N2;          
      D3[i][0]=D3[i-1][N2-1]+N3;} 
    for(j=1;j<N2;j++)
      D3[i][j]=D3[i][j-1]+N3;}  

  return(D3);
}

double ****CAllocD4(int N1, int N2, int N3, int N4){
  int i;
  double ****D4;

  D4=(double ****) malloc(N1*sizeof(double ***));
  *D4=(double ***) malloc(N1*N2*sizeof(double **));
  **D4=(double **) malloc(N1*N2*N3*sizeof(double *));
  ***D4=(double *) malloc(N1*N2*N3*N4*sizeof(double));
  
  for(i=1;i<N1*N2*N3;i++) D4[0][0][i]=D4[0][0][i-1]+N4;
  for(i=1;i<N1*N2   ;i++) D4[0][i]=D4[0][i-1]+N3;
  for(i=1;i<N1      ;i++) D4[i]=D4[i-1]+N2;

  return(D4);
}


//char's
char *CAllocC1(int N1){
  return((char *) malloc(N1*sizeof(char)));
}

char **CAllocC2(int N1, int N2){
  int i;
  char **C2;
  
  C2= (char **) malloc(N1*sizeof(char *));
  *C2= (char *) malloc(N1*N2*sizeof(char));
  for(i=1; i<N1; i++) C2[i]=C2[i-1]+N2;
  
  return(C2);
}


char ***CAllocC3(int N1, int N2, int N3){
  int i,j;
  char ***C3;

  C3=(char ***) malloc(N1*sizeof(char **));
  *C3=(char **) malloc(N1*N2*sizeof(char *));
  **C3=(char *) malloc(N1*N2*N3*sizeof(char));
  
  for(i=0;i<N1;i++){
    if(i!=0){
      C3[i]=C3[i-1]+N2;          
      C3[i][0]=C3[i-1][N2-1]+N3;} 
    for(j=1;j<N2;j++)
      C3[i][j]=C3[i][j-1]+N3;}  

  return(C3);
}


//deallocation
//unsigned's
void CDeAllocU1(unsigned *U1){
  free(U1);
  return;
}

void CDeAllocU2(unsigned **U2){
  free(*U2);
  free(U2);
  return;
}

void CDeAllocU3(unsigned ***U3){
  free(**U3);
  free(*U3);
  free(U3);
  return;
}

//int's
void CDeAllocI1(int *I1){
  free(I1);
  return;
}

void CDeAllocI2(int **I2){
  free(*I2);
  free(I2);
  return;
}

void CDeAllocI3(int ***I3){
  free(**I3);
  free(*I3);
  free(I3);
  return;
}

//float's
void CDeAllocF1(float *F1){
  free(F1);
  return;
}

void CDeAllocF2(float **F2){
  free(*F2);
  free(F2);
  return;
}

void CDeAllocF3(float ***F3){
  free(**F3);
  free(*F3);
  free(F3);
  return;
}

//double's
void CDeAllocD1(double *D1){
  free(D1);
  return;
}

void CDeAllocD2(double **D2){
  free(*D2);
  free(D2);
  return;
}

void CDeAllocD3(double ***D3){
  free(**D3);
  free(*D3);
  free(D3);
  return;
}

void CDeAllocD4(double ****D4){
  free(***D4);
  free(**D4);
  free(*D4);
  free(D4);
  return;
}

//char's
void CDeAllocC1(char *C1){
  free(C1);
  return;
}

void CDeAllocC2(char **C2){
  free(*C2);
  free(C2);
  return;
}

void CDeAllocC3(char ***C3){
  free(**C3);
  free(*C3);
  free(C3);
  return;
}

/***********************************/
/***********************************/
/**** string handling functions ****/
/***********************************/
/***********************************/

//allocate adequate space for a character array and copy a string to it 
char *allocstrcpy(const char *str){
  char *C1=CAllocC1(strlen(str)+1); //+1 is for terminating null
  return(strcpy(C1,str));
}
  
/************************************************************/
/************************************************************/
/**** functions to simplify input/output to binary files ****/
/************************************************************/
/************************************************************/

FILE *fopen_bin(char *filename, char *mode){
  FILE *fp;
  char fullmode[10]; //only expect 3 char + NULL, 10 is safety
  
  // if they didn't put a 'b' in the mode add it for them
  if(strchr(mode,(int)'b')) sprintf(fullmode,"%s",mode);
  else sprintf(fullmode,"%sb",mode);
  
  fp=fopen(filename, fullmode); //open the file
  setbuf(fp,(char *) NULL); //buffer the file 
  return(fp);
}

//read a single unsigned integer
void freadU(FILE *fp, unsigned *U){
  fread(U,sizeof(unsigned),1,fp);
  return;
}

//read a single integer
void freadI(FILE *fp, int *I){
  fread(I,sizeof(int),1,fp);
  return;
}

//read a single float
void freadF(FILE *fp, float *F){
  fread(F,sizeof(float),1,fp);
  return;
}

//read a single float to a double
void freadF2D(FILE *fp, double *D){
  float F;
  fread(&F,sizeof(float),1,fp);
  *D=F;
  return;
}

//read a single double
void freadD(FILE *fp, double *D){
  fread(D,sizeof(double),1,fp);
  return;
}

//read a single character
void freadC(FILE *fp, char *C){
  fread(C,sizeof(char),1,fp);
  return;
}

/* read a binary string: warning you must pass it an
   pointer NOT a statically allocated array */
void freadstring(FILE *fp, char **str){ 
  int num_char;
  freadI(fp,&num_char);
  *str=CAllocC1(num_char); 
  fread(*str,sizeof(char),num_char,fp);
  return;
}

//write a single unsigned integer
void fwriteU(FILE *fp, unsigned U){ 
  fwrite(&U,sizeof(unsigned),1,fp);
  return;
}

//write a single integer
void fwriteI(FILE *fp, int I){ 
  fwrite(&I,sizeof(int),1,fp);
  return;
}

//write a single float
void fwriteF(FILE *fp, float F){ 
  fwrite(&F,sizeof(float),1,fp);
  return;
}

//write a single double
void fwriteD(FILE *fp, double D){ 
  fwrite(&D,sizeof(double),1,fp);
  return;
}

//write a single character
void fwriteC(FILE *fp, char C){ 
  fwrite(&C,sizeof(char),1,fp);
  return;
}

//write a binary string 
void fwritestring(FILE *fp, char *str){ 
  int num_char=strlen(str)+1; 
  fwriteI(fp,num_char);
  fwrite(str,sizeof(char),num_char,fp);
  return;
}

/************************************************************/
/************************************************************/
/************* sorting and searching functions **************/
/************************************************************/
/************************************************************/

/* sort integer array and eliminate duplicate entries... should 
   overload this for unsigned ints, longs, floats and doubles */
void unique_sort(int *I1, int *N1){
  int i,j;
  int intswap;
  
  //sort in ascending order
  //should replace this linear sort with a call to qsort 
  for(j=0; j<*N1; j++)
    for(i=*N1-1; i>j ; i--)
      if(I1[j]>I1[i]){
	intswap=I1[j];
	I1[j]=I1[i];
	I1[i]=intswap;}
  
  //eliminate duplicate entries
  j=0;
  for(i=1; i<*N1; i++)
    if(I1[j]!=I1[i]){
      j++;
      if(i!=j) I1[j]=I1[i];}
  
  *N1=j+1;
  
  return;
}

void unique_sort_d(double *D1, int *N1){
  int i,j;
  double doubleswap;
  
  //sort in ascending order
  //should replace this linear sort with a call to qsort 
  for(j=0; j<*N1; j++)
    for(i=*N1-1; i>j ; i--)
      if(D1[j]>D1[i]){
	doubleswap=D1[j];
	D1[j]=D1[i];
	D1[i]=doubleswap;}
  
  //eliminate duplicate entries
  j=0;
  for(i=1; i<*N1; i++)
    if(D1[j]!=D1[i]){
      j++;
      if(i!=j) D1[j]=D1[i];}
  
  *N1=j+1;
  
  return;
}

/* search an integer array for a value and return it's indice
   really should use bsearch() in the stdlib instead of this 
   if you know the array is sorted in ascending order */
int searchI1(int *I1, int N1, int findme){
  int i;
  
  /* could replace this linear search with a bisection search
     to speed it up a little but then it would only work for
     sorted arrays ... see bsearch(), however I don't know 
     what bsearch() will do if the value your searching for 
     isn't there */
  for(i=0; i<N1; i++)
    if(I1[i]==findme) return(i);
  
  return(N1); /* if you can't find the entry return the
		 size of the array as an error code this
		 simulates the output of a linear search */
}


/* search an integer array for a value and return it's indice
   really should use bsearch() in the stdlib instead of this 
   if you know the array is sorted in ascending order */
int searchD1(double *D1, int N1, double findme){
  int i;
  
  /* could replace this linear search with a bisection search
     to speed it up a little but then it would only work for
     sorted arrays ... see bsearch(), however I don't know 
     what bsearch() will do if the value your searching for 
     isn't there */
  for(i=0; i<N1; i++)
    if(D1[i]==findme) return(i);
  
  return(N1); /* if you can't find the entry return the
		 size of the array as an error code this
		 simulates the output of a linear search */
}



/* taken from numerical recipes
   Minimal  random number generator of Park and Miller with Bays-Durham 
   shu e and added safeguards. Returns a uniform random deviate between 
   0.0 and 1.0 (exclusive of the endpoint values). Call with idum a 
   negative integer to initialize; thereafter, do not alter idum between
   successive deviates in a sequence. RNMX should approximate the largest
   floating value that is less than 1. 
*/

double ran1(long *idum) { 
  int j; 
  long k; 
  static long iy=0; 
  static long iv[NTAB]; 
  double temp; 
  if (*idum <= 0 || !iy) { //Initialize. 
    if (-(*idum) < 1) *idum=1; //Be sure to prevent idum = 0. 
    else *idum = -(*idum); 
    for (j=NTAB+7;j>=0;j--) { //Load the shuffle table (after 8 warm-ups). 
      k=(*idum)/IQ; 
      *idum=IA*(*idum-k*IQ)-IR*k; 
      if (*idum < 0) *idum += IM; 
      if (j < NTAB) iv[j] = *idum; 
    } 
    iy=iv[0]; 
  } 
  k=(*idum)/IQ; //Start here when not initializing. 
  *idum=IA*(*idum-k*IQ)-IR*k; //Compute idum=(IA*idum) % IM without over-  flows by Schrage s method. 
  if (*idum < 0) *idum += IM; 
  j=iy/NDIV; //Will be in the range 0..NTAB-1. 
  iy=iv[j]; //Output previously stored value and refill the shuffl e table. 
  iv[j] = *idum; 
  if ((temp=AM*iy) > RNMX) 
    return RNMX; //Because users don't expect endpoint values. 
  else 
    return temp; 
}

double key2double(unsigned key[2]) {
  EightBytes temp8;
  temp8.u[0]=key[0];
  temp8.u[1]=key[1];
  return(temp8.d);
}

unsigned* double2key(double D, unsigned key[2]) {
  EightBytes temp8;
  temp8.d=D;
  key[0]=temp8.u[0];
  key[1]=temp8.u[1];
  return(key);
}

double sign(double a) {
  if(a<0.0) return(-1.0); 
  else if(a>0.0) return(1.0);
  else return(0.0);
}
