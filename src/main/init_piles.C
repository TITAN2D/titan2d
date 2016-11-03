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
 * $Id: init_piles.C 134 2007-06-07 20:05:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#define PARABALOID
//#define CYLINDER
//#define PLANE
//#define CASITA
//#define POPO
//#define ID1
//#define ID2
int get_elem_elev(HashTable *HT_Node_Ptr, MatProps *matprops,
		  Element *EmTemp, double *elevation);

void print_grid(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, MatProps* matprops) {

  FILE *fp=fopen("gridplot00.txt","w");
  

  int num_buckets=HT_Elem_Ptr->get_no_of_buckets();
  for(int ibucket=0; ibucket<num_buckets; ibucket++){
    
    HashEntry *entryp = *(HT_Elem_Ptr->getbucketptr() + ibucket);

    //check every element in bucket
    while(entryp){	
      Element *EmTemp = (Element*)entryp->value;
      assert(EmTemp);
      
      if(EmTemp->get_adapted_flag()>0){
	double elevation;
	get_elem_elev(HT_Node_Ptr,matprops,EmTemp,&elevation);

	fprintf(fp,"%20.14g %20.14g %20.14g\n",
		(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE,
		(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE,
		elevation);
      }
      entryp = entryp->next;
    }
  }

  fclose(fp);

  exit(1);

  return;
}

void elliptical_pile_height(HashTable* HT_Node_Ptr, Element *EmTemp, 
			    MatProps* matprops_ptr, PileProps* pileprops_ptr);

void init_piles(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		int myid, int numprocs, int adaptflag,
		MatProps* matprops, TimeProps* timeprops_ptr, 
		MapNames* mapnames, PileProps* pileprops, 
		FluxProps *fluxprops, StatProps* statprops) { 

  unsigned nodes[9][KEYLENGTH], *node_key;
  int num_buckets=HT_Elem_Ptr->get_no_of_buckets();

  if(!adaptflag)
    H_adapt_to_level(HT_Elem_Ptr,HT_Node_Ptr,matprops,pileprops,
		     fluxprops,timeprops_ptr,REFINE_LEVEL);
#if defined PARABALOID || defined CYLINDER 
  if(adaptflag)
    initial_H_adapt(HT_Elem_Ptr, HT_Node_Ptr, 0, matprops, 
		    pileprops,fluxprops,timeprops_ptr,4); 
#else
  for(int ibucket=0; ibucket<num_buckets; ibucket++)
  {
    HashEntry *entryp = *(HT_Elem_Ptr->getbucketptr() + ibucket);
    //check every element in bucket
    while(entryp){	
      Element *EmTemp = (Element*)entryp->value;
      assert(EmTemp);
      if(EmTemp->get_adapted_flag()>0)
      {
	//put in the pile height right here...
	double* ndcoord = EmTemp->get_coord();
	double pile_height=0.0;
	double radius_sq;
	EmTemp->put_height(pileheight);
      } 
      entryp = entryp->next;
    }
  }
#endif //end "#if defined PARABALOID || defined CYLINDER"


  move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr,timeprops_ptr);
  slopes(HT_Elem_Ptr, HT_Node_Ptr, matprops);

  /* initial calculation of actual volume on the map */

  double realvolume=0.0, depositedvol=0.0, forcebed=0.0, meanslope=0.0;
  double epsilon[2]={matprops->epsilon, matprops->epsilon};

  HashEntryPtr* buck = HT_Elem_Ptr->getbucketptr();
  for(int ibucket=0; ibucket<HT_Elem_Ptr->get_no_of_buckets(); ibucket++)
    if(*(buck+ibucket)) {

      HashEntryPtr currentPtr = *(buck+ibucket);
      while(currentPtr) {

	Element* Curr_El=(Element*)(currentPtr->value);
        assert(Curr_El);

	if(Curr_El->get_adapted_flag()>0) { //if this is a refined element don't involve!!!
	  double *dxy=Curr_El->get_dx();
	  double dvol=dxy[0]*dxy[1]**(Curr_El->get_state_vars());
	  realvolume+=dvol;
	
	  Curr_El->put_kactxy(epsilon);
	  Curr_El->calc_stop_crit(matprops);
          if (Curr_El->get_stoppedflags()==2)
            depositedvol += dvol;
          
	  double resolution=0, xslope=0, yslope=0;
	  Get_max_resolution(&resolution);
	  Get_slope(resolution,
		    *((Curr_El->get_coord()))*matprops->LENGTH_SCALE,
		    *((Curr_El->get_coord())+1)*matprops->LENGTH_SCALE,
		    &xslope,&yslope);
	  double slope=sqrt(xslope*xslope+yslope*yslope); 

	  forcebed+=dvol*9.8/sqrt(1.0+slope*slope)*
	    tan(matprops->bedfrict[Curr_El->get_material()]);
	}
	currentPtr=currentPtr->next;      	    
      }
    }

  double tempin[3], tempout[3];
  tempin[0]=realvolume;
  tempin[1]=forcebed;
  tempin[2]=depositedvol;

  MPI_Reduce(tempin,tempout,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  statprops->realvolume=tempout[0]*(matprops->HEIGHT_SCALE)
    *(matprops->LENGTH_SCALE)*(matprops->LENGTH_SCALE);
  statprops->outflowvol=0.0;
  statprops->erodedvol=0.0;
  statprops->depositedvol=tempout[2]*(matprops->HEIGHT_SCALE)
    *(matprops->LENGTH_SCALE)*(matprops->LENGTH_SCALE);

  statprops->forceint=0.0;
  statprops->forcebed=tempout[1]/tempout[0];
  return;
}









/*******************************************************************/ 
/* assign height to point of an elliptical (in (x,y)) shaped pile, */
/* the pile can be either parabolic (in the z direction) or be     */
/* cylindrical (have uniform pile height)                          */
/*******************************************************************/
void elliptical_pile_height(HashTable* HT_Node_Ptr, Element *EmTemp, 
			    MatProps* matprops, PileProps* pileprops){
			      
  unsigned nodes[9][2]; 

  //get corner and edge nodes
  unsigned *node_key=EmTemp->getNode();
  for(int inode=0; inode<8; inode++)
    for(int ikey=0; ikey<KEYLENGTH; ikey++)
      nodes[inode][ikey]=node_key[inode*KEYLENGTH+ikey];
  
  //get center node
  node_key=EmTemp->pass_key();
  for(int ikey=0; ikey<KEYLENGTH; ikey++)
    nodes[8][ikey]=node_key[ikey];

  double node_pile_height[9];
  double sum_node_pile_height[9];
  double sum_node_xmom[9];
  double sum_node_ymom[9];
  double height;
  double vfract = 0.;

  for(int inode=0;inode<9;inode++)
  {
    //get pile height at each node...
    Node* ndtemp = (Node*) HT_Node_Ptr->lookup(&nodes[inode][0]);
    double* ndcoord = ndtemp->get_coord();
    
    // for multiple piles which may overlap, the highest value is used..
    node_pile_height[inode] = 0.0;
    sum_node_pile_height[inode] = 0.0;
    sum_node_xmom[inode] = 0.0;
    sum_node_ymom[inode] = 0.0;

    //check each pile to see which has max height at this node
    for(int ipile=0;ipile<pileprops->numpiles;ipile++)
    {
      if (pileprops->vol_fract[ipile] > vfract)
        vfract = pileprops->vol_fract[ipile];
        
      //get position relative to pile center
      double major=ndcoord[0]-pileprops->xCen[ipile];
      double minor=ndcoord[1]-pileprops->yCen[ipile];

      /* "undo" elliptical pile rotation ... from (x,y)->(major,minor) 
	 also make  nondimensional (by dividing by major and minor radius) */
      double doubleswap=(major* pileprops->cosrot[ipile]+
			 minor* pileprops->sinrot[ipile])/
	pileprops->majorrad[ipile];

      minor =(-major* pileprops->sinrot[ipile]+
	       minor* pileprops->cosrot[ipile])/
	pileprops->minorrad[ipile];
      major = doubleswap;

      /* calculate pile height based on non dimensional position relative to
	 center of pile */
#ifdef PARABALOID      
      height = pileprops->pileheight[ipile]*(1.-major*major-minor*minor);
#elif defined CYLINDER

      if (major*major+minor*minor<1.0)
	height = pileprops->pileheight[ipile];
      else height =0.0;
      
#endif
      height=(height>=0.0)?height:0.0;

      sum_node_pile_height[inode]+=height;
      sum_node_xmom[inode]+=height*(pileprops->initialVx[ipile]);
      sum_node_ymom[inode]+=height*(pileprops->initialVy[ipile]);

      if(node_pile_height[inode] < height)
	node_pile_height[inode] = height;
    }
    if(sum_node_pile_height[inode] <= GEOFLOW_TINY)
      sum_node_xmom[inode]=sum_node_ymom[inode]=0.0;
    else{
      sum_node_xmom[inode]*=height/sum_node_pile_height[inode];
      sum_node_ymom[inode]*=height/sum_node_pile_height[inode];
      //these are now the averaged momentums at each node
    }
  }


  /* The pile_height value assigned is an "area" weighted average over the 
     element's 9 nodes.  The element is divided into 4 squares, and each 
     corner of each of the 4 squares count once.  Because the center node 
     is repeated 4 times it's weight is 4 times as much as the element's 
     corner nodes which are not repeated; each edge node is repeated 
     twice */
  double pileheight=(//corner nodes
		     node_pile_height[0]+node_pile_height[1]+
		     node_pile_height[2]+node_pile_height[3]+
		     //edge nodes 
		     2.0*(node_pile_height[4]+node_pile_height[5]+
			  node_pile_height[6]+node_pile_height[7])+
		     //center node
		     4.0*node_pile_height[8]
		     )/16.0;

  double xmom=(//corner nodes
	       sum_node_xmom[0]+sum_node_xmom[1]+
	       sum_node_xmom[2]+sum_node_xmom[3]+
	       //edge nodes 
	       2.0*(sum_node_xmom[4]+sum_node_xmom[5]+
		    sum_node_xmom[6]+sum_node_xmom[7])+
	       //center node
	       4.0*sum_node_xmom[8]
	       )/16.0;

  double ymom=(//corner nodes
	       sum_node_ymom[0]+sum_node_ymom[1]+
	       sum_node_ymom[2]+sum_node_ymom[3]+
	       //edge nodes 
	       2.0*(sum_node_ymom[4]+sum_node_ymom[5]+
		    sum_node_ymom[6]+sum_node_ymom[7])+
	       //center node
	       4.0*sum_node_ymom[8]
	       )/16.0;
  EmTemp->put_height_mom(pileheight,vfract,xmom,ymom);
  return;
}
