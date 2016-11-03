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
 * $Id: update_topo.C 134 2007-06-07 20:05:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#define WORKDIR "."

int update_topo(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		int myid, int nump, MatProps* matprops, TimeProps* timeprops, 
		MapNames* mapnames){

  //printf("update_topo.C 1\n"); fflush(stdout);

  char gis_update_map[100]="\0";
  char yada[256];
  if(myid==0){
    //printf("update_topo.C 2\n"); fflush(stdout);
     sprintf(yada,"%s/%s",WORKDIR,"updatetopo.txt");
     ifstream inp(yada,ios::in);
    //printf("update_topo.C 3\n"); fflush(stdout);
    if(!inp.fail()){ 
      inp>>gis_update_map;
      inp.close();}}
  //printf("update_topo.C 4\n"); fflush(stdout);

  if(nump>0) MPI_Bcast(gis_update_map,100,MPI_CHAR,0,MPI_COMM_WORLD);
  //printf("update_topo.C 5\n"); fflush(stdout);

  //check for a update map to open
  if(strlen(gis_update_map)>0){
    double tick= MPI_Wtime();
    //long tick=(long) time(NULL);
    Delete_GIS_data();
    Initialize_GIS_data(mapnames->gis_main  , mapnames->gis_sub, 
			mapnames->gis_mapset, gis_update_map);
    /*    Update_GIS_data(mapnames->gis_main  , mapnames->gis_sub, 
	  mapnames->gis_mapset, gis_update_map); */


    /* update the nodes */
    int num_buckets=HT_Node_Ptr->get_no_of_buckets();
    int ibucket;
    HashEntry *entryp;
    //visit every bucket
    for(ibucket=0; ibucket<num_buckets; ibucket++){
      entryp = *(HT_Node_Ptr->getbucketptr() + ibucket);
      //visit every node in this bucket
      while(entryp){
	Node* NdTemp= (Node*) entryp->value;
	assert(NdTemp);
	NdTemp->set_elevation(matprops);
	entryp = entryp->next;}
    }
    

    /* update the elements */
    num_buckets=HT_Elem_Ptr->get_no_of_buckets();

    //visit every bucket
    for(ibucket=0; ibucket<num_buckets; ibucket++){
      entryp = *(HT_Elem_Ptr->getbucketptr() + ibucket);
      
      //visit every element in this bucket
      while(entryp){	
	Element *EmTemp = (Element*)entryp->value;
	assert(EmTemp);
	
	if(EmTemp->get_adapted_flag()>0){
	  //update the topography in same order as in element2.C
	  //double eldif=EmTemp->get_elevation();
	  EmTemp->calc_topo_data(matprops);
	  //eldif=(EmTemp->get_elevation()-eldif)*matprops->LENGTH_SCALE;
	  //if(fabs(eldif)>1.0) printf("update_topo() after-before=%g\n",eldif);

	  EmTemp->calc_gravity_vector(matprops);
	  EmTemp->calc_d_gravity(HT_Elem_Ptr);
	}
	entryp = entryp->next;
      }
    }//closes: for(int ibucket=0; ibucket<num_buckets; ibucket++)
    double tock=MPI_Wtime()-tick;
    //long tock=(long) time(NULL);
    if(nump>1){
      /* long tempin[2], tempout[2];
      tempin[0]=-tick;
      tempin[1]= tock;
      MPI_Reduce(tempin, tempout, 2, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
      tick=-tempout[0];
      tock=tempout[1]; */
      MPI_Reduce(&tock, &tick, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      tock=tick;

    }

    if(myid==0){
      sprintf(yada,"%s/%s",WORKDIR,"updatetopo.done");
      FILE *fp=fopen(yada,"w");
      fprintf(fp,"%s",gis_update_map);
      fclose(fp);

      sprintf(yada,"%s/%s",WORKDIR,"updatetopo.time");
      fp=fopen(yada,"a");
      int hours, minutes; double seconds;
      timeprops->chunktime(&hours,&minutes,&seconds);
      fprintf(fp,"%s timestep=%d simtime=(%02d:%02d%g) (hrs:min:sec) timing=%g (sec)\n",gis_update_map,timeprops->iter,hours,minutes,seconds,tock);
      fclose(fp);}

    return(1); //successfull updated topography
  }//closes: if(strlen(gis_update_map>0))
  else 
    return(0); //no update gis map to open
}
	       
