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
 * $Id: hdfdefs.h 2 2003-08-13 19:26:11Z sorokine $ 
 */

#define DATASETNAME "ExtendibleArray1" 
#define DATASETNAME1 "ExtendibleArray2" 
#define RANK        2
#define LENGTH 1  
#define HEIGHT 10000
#define WIDTH 10
#define HEIGHT1    2500
#define aHEIGHT 10
#define aWIDTH 4
#define aHEIGHT1    2

float      data[HEIGHT][WIDTH];
float      data1[HEIGHT1][WIDTH];
float      adata[aHEIGHT][aWIDTH];
float      adata1[aHEIGHT1][aWIDTH];
hsize_t    maxdim[2] = {H5S_UNLIMITED, WIDTH}; 
hsize_t    amaxdim[2] = {H5S_UNLIMITED,aWIDTH};
hsize_t    dim[2] = {HEIGHT, WIDTH};   
hsize_t    adim[2] = {aHEIGHT, aWIDTH};   
hid_t      cparms,acparms;
hsize_t    newsize[2]={HEIGHT, WIDTH};
hsize_t    anewsize[2]={aHEIGHT, aWIDTH};
void       *tbuf = NULL;
  

int w[100],aw[100],r[100],ar[100],l[100],al[100],al1[100],l1[100];
