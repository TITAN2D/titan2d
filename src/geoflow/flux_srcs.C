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
 * $Id: init_piles.C,v 1.4 2004/08/11 16:14:20 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

/* The cells are marked if the are with in the flux-source region,
   even if in future time. This will reduce the work later when 
   material contribution of each cell is calculated.
   
   Some of the cells marked active will go passive with refinement
   but none of the cells marked passive should go active
 */

void mark_flux_region(HashTable* ElemTable, HashTable* NodeTable,
		      MatProps *matprops, FluxProps *fluxprops, TimeProps *timeprops)
{

  int buckets=ElemTable->get_no_of_buckets();

  if(fluxprops->MaxInfluxNow(matprops,timeprops)>0.0) {
    /*
    if(timeprops->timesec()>30.0){
      printf("flux not stopped after 30 seconds\n");
      assert(0);
    }
    */
    // mdj 2007-04
    HashEntry *entryptr;
    Element *EmTemp;
#pragma omp parallel for private(entryptr,EmTemp)
    for(int ibuck=0; ibuck<buckets; ibuck++) {
      entryptr=*(ElemTable->getbucketptr()+ibuck);

      while (entryptr) {
	EmTemp = (Element*)entryptr->value;
	assert(EmTemp);
	if(EmTemp->get_adapted_flag()>0)
	  //if this element doesn't belong on this processor don't involve
	  EmTemp->calc_flux(NodeTable,fluxprops,timeprops);
	    
	entryptr=entryptr->next;
      }
    }
  }
  else{

    // mdj 2007-04
    HashEntry *entryptr;
    Element *EmTemp;
#pragma omp parallel for private(entryptr,EmTemp)
    for(int ibuck=0; ibuck<buckets; ibuck++) {
      entryptr=*(ElemTable->getbucketptr()+ibuck);

      while (entryptr) {
	EmTemp = (Element*)entryptr->value;
	assert(EmTemp);
	if(EmTemp->get_adapted_flag()>0)
	  //if this element doesn't belong on this processor don't involve
	  EmTemp->zero_influx();
	    
	entryptr=entryptr->next;
      }
    }
  }


  return;
}


void adapt_fluxsrc_region(HashTable *ElemTable, HashTable *NodeTable, 
			  MatProps *matprops, PileProps *pileprops,
			  FluxProps *fluxprops, TimeProps *timeprops, 
			  double dt, int myid, int adaptflag) {

  mark_flux_region(ElemTable,NodeTable,matprops,fluxprops,timeprops);

  if(//(adaptflag)&&
     (fluxprops->IfAnyStart(timeprops))&&
     (timeprops->iter>0) //iteration zero flux adaptation happens at same time as pile adaptation don't waste cpu work by doing it a second time
     )
    //initial_H_adapt adapts the grid for initial pile and flux sources
    initial_H_adapt(ElemTable,NodeTable,0,matprops,pileprops,fluxprops,
		    timeprops,5);

  return;
}



void Element::calc_flux(HashTable *NodeTable,
			FluxProps *fluxprops, TimeProps *timeprops) {

  int no_of_sources=fluxprops->no_of_sources;
  unsigned node[9][2];
  double temp_coef;
  double major,minor,dswap,xcoord,ycoord;
  int ikey,inode,isrc;
  double curr_time=timeprops->time;
  double start_time, end_time;
  /*
  printf("no_of_sources=%g time=%g\n",no_of_sources,timeprops->time);
  
  for(int isrc=0; isrc<no_of_sources; isrc++)
    printf("influx[%d]=%g\n",isrc,fluxprops->influx[isrc]);
  */
  Influx[0]=Influx[1]=Influx[2]=0.0;

  //get corner and edge nodes
  for(inode=0; inode<8; inode++)
    for(ikey=0; ikey<KEYLENGTH; ikey++)
      node[inode][ikey]=*(this->getNode()+inode*KEYLENGTH+ikey);
                                                                                                                                                             
  for(ikey=0; ikey<KEYLENGTH; ikey++)
  // center node key is same as element key
    node[8][ikey]=*(this->pass_key()+ikey);

  for(inode=0; inode<9; inode++) {
    double temp_coef2=temp_coef=0.0;
    double sum_flux_xmom_ymom[3]={0.0,0.0,0.0};
    Node* ndtemp=(Node*) NodeTable->lookup(&node[inode][0]);
    assert(ndtemp);
    xcoord=*(ndtemp->get_coord());
    ycoord=*(ndtemp->get_coord()+1);
      
    for(int isrc=0; isrc<no_of_sources; isrc++) {
      // normalize start and end time
      start_time=fluxprops->start_time[isrc];
      end_time =fluxprops->end_time[isrc];

      if((fluxprops->start_time[isrc] <= curr_time) && 
	 (fluxprops->end_time[isrc]   >= curr_time)) {

	/*
	if(timeprops->timesec()>30.0){
	  printf("flux not stopped after 30 seconds\n");
	  assert(0);
	}
	*/  
	  
	
	// "undo" elliptical pile rotation ... from (x,y)->(major,minor) 
        major=xcoord-fluxprops->xCen[isrc];
        minor=ycoord-fluxprops->yCen[isrc];
        dswap=( major*fluxprops->cosrot[isrc]+
                minor*fluxprops->sinrot[isrc])/fluxprops->majorrad[isrc];
        minor=(-major*fluxprops->sinrot[isrc]+
                minor*fluxprops->cosrot[isrc])/fluxprops->minorrad[isrc];
        major=dswap;
        if (major*major+minor*minor<1.0) {
          if (inode<4)
            temp_coef=1.0/16.0;
          if (inode>=4 && inode<8)
            temp_coef=1.0/8.0;
          if (inode==8) {
            temp_coef=1.0/4.0;
	    //printf("xcoord=%g ycoord=%g isrc=%d major=%g minor=%g\n",
	    //   xcoord,ycoord,isrc,major,minor);
	  }
        }
        /* 
         * distance from cell-center is used regardless of cell-center 
         * being inside or outside the source region
         */

	//hydrograph flux starts at max rate decays linearly to zero at end of duration
	double tempflux=(fluxprops->influx[isrc])*
	  (1.0-(curr_time-fluxprops->start_time[isrc])/
	   (fluxprops->end_time[isrc]-fluxprops->start_time[isrc]));
        temp_coef*=(1.0-major*major-minor*minor)*tempflux;
        //temp_coef*=(1.0-major*major-minor*minor)*(fluxprops->influx[isrc]);
	
	if(temp_coef>0.0) {
	  sum_flux_xmom_ymom[0]+=temp_coef;
	  sum_flux_xmom_ymom[1]+=temp_coef*fluxprops->xVel[isrc];
	  sum_flux_xmom_ymom[2]+=temp_coef*fluxprops->yVel[isrc];
	}

        if(temp_coef > temp_coef2) {
          temp_coef2=temp_coef;

	  //printf("xcoord=%g ycoord=%g isrc=%d major=%g minor=%g\n",
	  //	 xcoord,ycoord,isrc,major,minor);
	}

	
      }

    }
    if(temp_coef2 > 0.0) {
      Influx[0]+=temp_coef2;

      //if(sum_flux_xmom_ymom[0]>0.0) {
	Influx[1]+=sum_flux_xmom_ymom[1]*temp_coef2/sum_flux_xmom_ymom[0];
	Influx[2]+=sum_flux_xmom_ymom[2]*temp_coef2/sum_flux_xmom_ymom[0];
	//}
    }
      

  }
  
  if(!(Influx[0]>=0.0)){
    printf("error in Influx[0] assignment\n");
    assert(Influx[0]>=0.0);
  } 
  Awet=(Influx[0]>0.0)?1.0:0.0;


  return;
}

