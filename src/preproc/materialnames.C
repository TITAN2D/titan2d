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
 * $Id: materialnames.C 17 2003-11-24 16:01:18Z kdalbey $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "GisApi.h"

/* this program receives the following as command line input when 
   called by titan_gui.py
   1) the full path of the GIS database
   2) the location
   3) the mapset
   4) the material property raster map name

   it loads the map dictated by the command line input, retrieves
   a list of material names from the map, and prints the list to
   a file named "materialnames.dat"
*/

int main(int argc, char** argv){
  if((argc != 5)) {
    printf("You entered %d arguments, materialnames.C requires 4 arguments.  They are:\n 1) the full path of the GIS database,\n 2) the location,\n 3) the mapset,\n 4) the material property raster map name.\n",argc-1);
    exit(1);}

  /*  printf("[%s]\n",argv[0]);
  printf("[%s]\n",argv[1]);
  printf("[%s]\n",argv[2]);
  printf("[%s]\n",argv[3]); */

  char *gis_matmap=(char *) malloc((strlen(argv[4])+5)*sizeof(char));
  strcpy(gis_matmap,argv[4]);
  strcat(gis_matmap+strlen(argv[4]),"_Mat");

  if(!Initialize_Raster_data(argv[1],argv[2],argv[3],gis_matmap)){
    //write the material names data file to be read by titan_gui.py
    FILE *fp=fopen("materialnames.dat","w");

    //determine the number of materials
    int nummat; //the number of materials
    Get_raster_categories(&nummat);
    fprintf(fp,"%d\n",nummat);

    //retrieve the list of material names from the map
    char materialname[200];

    /* note material id (imat) =0 is the "null material", it is not 
       used on the map "ColimaSmallSRTM_Mat" (which is the first test
       map) but could be used on other maps where the material is not 
       known, for instance ground covered by water */
    for(int imat=1;imat<=nummat;imat++){ 
      Get_raster_category_name(imat,materialname);
      fprintf(fp,"%s\n",materialname);}

    fclose(fp);

    Delete_Raster_data();

  }
  else{
    printf("Couldn't initialize the GIS Material Property information.\n");
    exit(1);}

  return(0);
}
