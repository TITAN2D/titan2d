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
 * $Id: datread.C 134 2007-06-07 20:05:52Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

#ifdef SUNOS
extern "C" void initial_(int*, double*, double*);
#endif

#ifdef IBMSP
extern "C" void initial(int*, double*, double*);
#endif

#ifdef CRAY
extern "C" void INITIAL(int*, double*, double*);
#endif


void Read_data(int myid, MatProps* matprops_ptr, PileProps* pileprops_ptr, 
	       StatProps* statprops_ptr, TimeProps* timeprops_ptr, FluxProps* fluxprops,
	       int* adaptflag_ptr, int* viz_flag_ptr, int* order_flag_ptr,
	       MapNames *mapnames_ptr, DISCHARGE* discharge_ptr, 
	       OutLine* outline_ptr, int *srctype)
{ 
  /*************************************************************************/
  //regular pile info input

  //read in elliptical (in x-y coordinates) pile info
  ifstream inD2("simulation.data", ios::in);
  if(inD2.fail()) {
    printf("can't find simulation.data file\n");
    exit(0);
  }
  inD2>>*srctype; // is it a flux source or pile source?
  int numpiles = 0;  // the number of separate files to use
  int no_of_sources=0;
  int isrc;

  if ((*srctype)&0x1)
    inD2>>numpiles;
  if ((*srctype)&0x2)
    inD2>>no_of_sources;

  if ((*srctype)&0x1)
  {
    pileprops_ptr->allocpiles(numpiles);

    double rotang;
    double doubleswap1, doubleswap2;
    double maxphi=0.;
    double minphi=HUGE_VAL;
    int imat;

    for(isrc=0;isrc<numpiles;isrc++) {
      inD2>>pileprops_ptr->pileheight[isrc]; // pile height

      // solid-volume fraction
      inD2>>pileprops_ptr->vol_fract[isrc]; 
      // search for min-max phi
      if (pileprops_ptr->vol_fract[isrc] > maxphi)
        maxphi = pileprops_ptr->vol_fract[isrc];
      if (pileprops_ptr->vol_fract[isrc] < minphi)
        minphi = pileprops_ptr->vol_fract[isrc];

      inD2>>pileprops_ptr->xCen[isrc];       // pile x-center
      inD2>>pileprops_ptr->yCen[isrc];       // pile y-center
      inD2>>pileprops_ptr->majorrad[isrc];   /* pile "major axis" radius
					        x-radius if no rotation */
      inD2>>pileprops_ptr->minorrad[  isrc]; /* pile minor axis radius
					       y-radius if no rotation */
      inD2>>rotang; // counter clockwise pile rotation angle
      pileprops_ptr->cosrot[isrc]=cos(rotang*PI/180.0);
      pileprops_ptr->sinrot[isrc]=sin(rotang*PI/180.0);

      inD2>>doubleswap1; //flow speed
      inD2>>doubleswap2; //flow direction measured counter clockwise from x axis
      doubleswap2*=PI/180.0;
      pileprops_ptr->initialVx[isrc]=doubleswap1*cos(doubleswap2);
      pileprops_ptr->initialVy[isrc]=doubleswap1*sin(doubleswap2);

    }

    // cut-off extremes
    assert(minphi <= maxphi);
    if ( minphi < 0.2 )
      matprops_ptr->flow_type = FLUID_FLOW;
    else if ( maxphi > 0.9 )
      matprops_ptr->flow_type = DRY_FLOW;
    else
      matprops_ptr->flow_type = TWOPHASE;
  }

  if ((*srctype)&0x2)
  {
    double rotang=0;
    double vel,vel_angle;
    fluxprops->allocsrcs(no_of_sources);

    for (isrc=0; isrc<no_of_sources; isrc++)
    {
      inD2>>fluxprops->influx[    isrc];
      inD2>>fluxprops->start_time[isrc];
      inD2>>fluxprops->end_time[  isrc];
      inD2>>fluxprops->xCen[      isrc];
      inD2>>fluxprops->yCen[      isrc];
      inD2>>fluxprops->majorrad  [isrc];
      inD2>>fluxprops->minorrad[  isrc];
      inD2>>rotang;
      fluxprops->cosrot[isrc]=cos(rotang*PI/180);
      fluxprops->sinrot[isrc]=sin(rotang*PI/180);
      inD2>>vel;
      inD2>>vel_angle;
      vel_angle*=PI/180.;
      fluxprops->xVel[isrc]=vel*cos(vel_angle);
      fluxprops->yVel[isrc]=vel*sin(vel_angle);
      
    }

  }
  if (!srctype)
  {
    printf("ERROR: No material source was defined");
    exit(1);
  }
  

  /*************************************************************************/
  //over ride regular input with statistical sample run data

  FILE* fp;
  int lhsref, lhsid, ifstatbed=0;
  int i;
  double statbed, statint;  //friction angles
  if((fp=fopen("statin.bed","r"))!=NULL) {
    ifstatbed=1;
    fscanf(fp,"%d%d%lf%lf",&(statprops_ptr->lhs.refnum),
	   &(statprops_ptr->lhs.runid),&statbed,&statint);
    fclose(fp);
    statbed*=PI/180.0;
    statint*=PI/180.0;
  }
  else if((fp=fopen("statin.vol","r"))!=NULL) {
    double volumescale;
    fscanf(fp,"%d%d%lf%lf",&(statprops_ptr->lhs.refnum),
	   &(statprops_ptr->lhs.runid),&volumescale);
    fclose(fp);
    volumescale=pow(volumescale,1/3.0);

    for(i=0;i<numpiles;i++) {
      pileprops_ptr->pileheight[i]*=volumescale;
      pileprops_ptr->majorrad[i]*=volumescale;
      pileprops_ptr->minorrad[i]*=volumescale;
    }
  }

  //PCQ
  int isample=-1, ifpcqvolbed=0;
  double pcqbedfrict;
  fp=fopen("sample.number","r");
  if(fp!=NULL) {
    fscanf(fp,"%d",&isample);
    fclose(fp);
    statprops_ptr->lhs.runid=isample;
    double volume;
    int intswap;
    char samplefilename[256];
    sprintf(samplefilename,"dirfluxsample.%06d",isample);
    fp=fopen(samplefilename,"r");
    if(fp!=NULL) {
      assert(fluxprops->no_of_sources);
      fluxprops->no_of_sources=1;

      char stringswap[4096];
      fgets(stringswap,4096,fp);  //Nsample
      fscanf(fp,"isample=%d\n",&intswap);
      assert(intswap==isample);
      fgets(stringswap,4096,fp);  //weight
      fgets(stringswap,4096,fp);  //Nrand
      fgets(stringswap,4096,fp);  //Nvol
      fgets(stringswap,4096,fp);  //Ndir
      fgets(stringswap,4096,fp);  //randvar
      fscanf(fp,"volume=%lf\n",&volume);
      double vel=sqrt(fluxprops->xVel[0]*fluxprops->xVel[0]+
		      fluxprops->yVel[0]*fluxprops->yVel[0]);
      double vel_angle;
      fscanf(fp,"direction=%lf\n",&vel_angle);
      vel_angle*=PI/180.;
      fluxprops->xVel[0]=vel*cos(vel_angle);
      fluxprops->yVel[0]=vel*sin(vel_angle);
      fluxprops->start_time[0]=0.0;
      fscanf(fp,"flux duration=%lf\n",&(fluxprops->end_time[0]));
      fscanf(fp,"init center flux=%lf\n",&(fluxprops->influx[0]));
      printf("Vol=%g [m^3]\n",volume);
    }
  }

  statprops_ptr->runid=statprops_ptr->lhs.runid;
  /*************************************************************************/
  //scaling info

  ifstream inscale("scale.data",ios::in);
  if(inscale.fail() == 0) {  
    inscale>>matprops_ptr->LENGTH_SCALE;
    inscale>>matprops_ptr->HEIGHT_SCALE; 
    inscale>>matprops_ptr->GRAVITY_SCALE;
  }
  else{  // if this file doesn't exist, assume no scaling
    printf("Can't find scale.data for processor %d\n",myid);
    matprops_ptr->LENGTH_SCALE = 1;
    matprops_ptr->HEIGHT_SCALE = 1; 
    matprops_ptr->GRAVITY_SCALE = 1;
  }
  inscale.close();
  matprops_ptr->epsilon = matprops_ptr->HEIGHT_SCALE/matprops_ptr->LENGTH_SCALE;

  //this is used in ../geoflow/stats.C ... might want to set 
  //MAX_NEGLIGIBLE_HEIGHT to zero now that we have "good" thin 
  //layer control, need to reevaluate this, we should also 
  //reevaluate after we implement a "good" local stopping criteria

  matprops_ptr->MAX_NEGLIGIBLE_HEIGHT=matprops_ptr->HEIGHT_SCALE/10000.0;

  double TIME_SCALE=
    sqrt(matprops_ptr->LENGTH_SCALE/matprops_ptr->GRAVITY_SCALE);

  //non-dimensionalize the inputs
  double VELOCITY_SCALE=
    sqrt(matprops_ptr->LENGTH_SCALE*matprops_ptr->GRAVITY_SCALE);

  double diameter = 0.005;
  double vterm = pow(diameter,2.)*(matprops_ptr->den_solid -
                 matprops_ptr->den_fluid)*matprops_ptr->GRAVITY_SCALE/
                 (18.*matprops_ptr->viscosity);
  matprops_ptr->v_terminal = vterm;
  double smallestpileradius=HUGE_VAL;

  for(isrc=0;isrc<pileprops_ptr->numpiles;isrc++) { 
    pileprops_ptr->pileheight[isrc]/=matprops_ptr->HEIGHT_SCALE;
    pileprops_ptr->xCen[      isrc]/=matprops_ptr->LENGTH_SCALE;
    pileprops_ptr->yCen[      isrc]/=matprops_ptr->LENGTH_SCALE;
    pileprops_ptr->majorrad[  isrc]/=matprops_ptr->LENGTH_SCALE;
    pileprops_ptr->minorrad[  isrc]/=matprops_ptr->LENGTH_SCALE;
    pileprops_ptr->initialVx[ isrc]/=VELOCITY_SCALE;
    pileprops_ptr->initialVy[ isrc]/=VELOCITY_SCALE;

    if(smallestpileradius>pileprops_ptr->majorrad[isrc])
      smallestpileradius =pileprops_ptr->majorrad[isrc];
    
    if(smallestpileradius>pileprops_ptr->minorrad[isrc])
      smallestpileradius =pileprops_ptr->minorrad[isrc];
  }

  for(isrc=0;isrc<fluxprops->no_of_sources;isrc++) {
    fluxprops->influx[    isrc]*=TIME_SCALE/(matprops_ptr->HEIGHT_SCALE);
    fluxprops->start_time[isrc]/=TIME_SCALE;
    fluxprops->end_time[  isrc]/=TIME_SCALE;
    fluxprops->xCen[      isrc]/=matprops_ptr->LENGTH_SCALE;
    fluxprops->yCen[      isrc]/=matprops_ptr->LENGTH_SCALE;
    fluxprops->majorrad[  isrc]/=matprops_ptr->LENGTH_SCALE;
    fluxprops->minorrad[  isrc]/=matprops_ptr->LENGTH_SCALE;
    fluxprops->xVel[      isrc]/=VELOCITY_SCALE;
    fluxprops->yVel[      isrc]/=VELOCITY_SCALE;

    if(smallestpileradius>fluxprops->majorrad[isrc])
      smallestpileradius =fluxprops->majorrad[isrc];

    if(smallestpileradius>fluxprops->minorrad[isrc])
      smallestpileradius =fluxprops->minorrad[isrc];
  }
  //bob

  matprops_ptr->smallest_axis=2.0*smallestpileradius;
  inD2>>matprops_ptr->number_of_cells_across_axis;
  
  /*************************************************************************/
  /* the non-dimensional velocity stopping criteria is an idea that 
     didn't work for anything other than a slumping pile on a horizontal 
     surface, it's only still here because I didn't want to bother
     with removing it.  --Keith Dalbey 2005.10.28

     kappa is a to be determined constant, calculation stops when 
     v*=v_ave/v_slump<kappa (or perhaps v* < kappa/tan(intfrict)) */
  double kappa=1.0;   //should eventually move to a header file
  double gravity=9.8; //[m/s^2]
  matprops_ptr->Vslump=1.0; //kappa*sqrt(gravity*max_init_height);


  /*************************************************************************/
  //time related info

  int maxiter;
  double maxtime, timeoutput, timesave;
  inD2>>maxiter;
  inD2>>maxtime;
  inD2>>timeoutput;
  inD2>>timesave;

  timeprops_ptr->inittime(maxiter,maxtime,timeoutput,timesave,TIME_SCALE);


  /*************************************************************************/
  //flags

  inD2>>*(adaptflag_ptr);
  inD2>>*viz_flag_ptr;
  inD2>>*order_flag_ptr;

  /*************************************************************************/
  // read in GIS information

  char gis_main[300];
  char gis_sub[100];
  char gis_mapset[100];
  char gis_map[100];
  int  extramaps;

  inD2>>gis_main;
  inD2>>gis_sub;
  inD2>>gis_mapset;
  inD2>>gis_map;
  inD2>>extramaps;

  mapnames_ptr->assign(gis_main,gis_sub,gis_mapset,gis_map,extramaps);

  i = Initialize_GIS_data(gis_main, gis_sub, gis_mapset, gis_map);
  if(i!= 0) {
    printf("Problem with GIS on processor %d\n",myid);
    exit(0);
  }

  /*************************************************************************/
  //test point information

  inD2>>statprops_ptr->hxyminmax;
  if(statprops_ptr->hxyminmax==-1.0)
    statprops_ptr->hxyminmax=
      matprops_ptr->MAX_NEGLIGIBLE_HEIGHT*10.0;
  else if(statprops_ptr->hxyminmax<=0.0) {
    printf("bogus edge height=%g read in from simulation.data\n(edge height = -1 is the flag for using the default height)\nExitting!!\n",
	   statprops_ptr->hxyminmax);
    exit(0);
  }
  statprops_ptr->hxyminmax/=matprops_ptr->HEIGHT_SCALE;

  inD2>>statprops_ptr->heightifreach;
  char chargarb1[64], chargarb2[64];
  if(statprops_ptr->heightifreach!=-2) {

    //default test height is 10 time the maximum negligible height
    if(statprops_ptr->heightifreach==-1)
      statprops_ptr->heightifreach=
	matprops_ptr->MAX_NEGLIGIBLE_HEIGHT*10.0;

    statprops_ptr->heightifreach/=matprops_ptr->HEIGHT_SCALE; 

    inD2>>statprops_ptr->xyifreach[0]>>statprops_ptr->xyifreach[1];
    statprops_ptr->xyifreach[0]/=(matprops_ptr->LENGTH_SCALE);
    statprops_ptr->xyifreach[1]/=(matprops_ptr->LENGTH_SCALE);
  }
  else{
    inD2>>chargarb1>>chargarb2;
    statprops_ptr->heightifreach=
      statprops_ptr->xyifreach[0]=
      statprops_ptr->xyifreach[1]=HUGE_VAL;
  }

  //to get rid on uninitiallized memory error in saverun() (restart.C)
  statprops_ptr->forceint=statprops_ptr->forcebed=0.0;
  

  /*************************************************************************/
  //the discharge plane section starts here

  int iplane, num_planes;
  double **planes;
  inD2>>num_planes;
  //printf("num_planes=%d\n",num_planes);
  if(num_planes>0) {
    planes=CAllocD2(num_planes,4);
    
    for(iplane=0;iplane<num_planes;iplane++) {
      inD2>>planes[iplane][0]>>planes[iplane][2]>>planes[iplane][1]>>planes[iplane][3];
      //printf("plane %d: (%16.10g,%16.10g) (%16.10g,%16.10g)\n",iplane,planes[iplane][0],planes[iplane][2],planes[iplane][1],planes[iplane][3]);
      planes[iplane][0]/=matprops_ptr->LENGTH_SCALE;
      planes[iplane][1]/=matprops_ptr->LENGTH_SCALE;
      planes[iplane][2]/=matprops_ptr->LENGTH_SCALE;
      planes[iplane][3]/=matprops_ptr->LENGTH_SCALE; 
      //printf("plane %d: (%16.10g,%16.10g) (%16.10g,%16.10g)\n",iplane,planes[iplane][0],planes[iplane][2],planes[iplane][1],planes[iplane][3]);
    }

  }
  
  discharge_ptr->init(num_planes,planes);

  if(num_planes>0) CDeAllocD2(planes);
  //the discharge plane section ends here

  inD2.close(); 


  /*************************************************************************/

  //read in material properties
  fp=fopen("frict.data","r");
  fscanf(fp,"%d\n",&(matprops_ptr->material_count));
  
  matprops_ptr->matnames=(char **) malloc((matprops_ptr->material_count+1)*sizeof(char *));
  matprops_ptr->bedfrict=CAllocD1(matprops_ptr->material_count+1);  
  matprops_ptr->tanbedfrict=CAllocD1(matprops_ptr->material_count+1);
  char stringswap[512];
  int imat; 
  double doubleswap;
  for(imat=1;imat<=matprops_ptr->material_count;imat++) {
    fgets(stringswap,512,fp);
    matprops_ptr->matnames[imat]=allocstrcpy(stringswap);
    fscanf(fp,"%lf %lf\n",&doubleswap,&(matprops_ptr->bedfrict[imat]));
    matprops_ptr->bedfrict[imat]*=PI/180.0;
    matprops_ptr->tanbedfrict[imat]=tan(matprops_ptr->intfrict);

    if(imat==1) {
      matprops_ptr->intfrict=doubleswap*PI/180.0;
      matprops_ptr->tanintfrict=tan(matprops_ptr->intfrict);
    }

  }

  if(matprops_ptr->material_count>1) {
    char *gis_matmap=(char *) malloc((strlen(gis_map)+5)*sizeof(char));
    strcpy(gis_matmap,gis_map);
    strcat(gis_matmap+strlen(gis_map),"_Mat");
  
    if(Initialize_Raster_data(gis_main, gis_sub, gis_mapset, gis_matmap)) {
      printf("Problem with GIS Material on processor %d\n",myid);
      exit(0);}

    free(gis_matmap);

    int nummat; 
    Get_raster_categories(&nummat);
    if(nummat!=matprops_ptr->material_count) {
      printf("frict.data has %d materials but material map has %d, aborting\n",matprops_ptr->material_count,nummat);
      exit(0);
    }
  }

  /* physically we don't know how to specify a changing internal friction 
     angle, we don't know how the contents of an avalanche changes as it 
     picks up new material.  We don't know how to keep track of the material 
     it picks up.  So we say that the internal friction angle is a constant */


  //replace frict.data values of bed friction with ones from a
  //statistics sample run file read in above
  if(ifstatbed) {
    matprops_ptr->intfrict = statint;
    for(int imat=1;imat<=matprops_ptr->material_count;imat++)
      matprops_ptr->bedfrict[imat] = statbed;
  }

  if(ifpcqvolbed) { //printf("PCQ bedfrict yada\n");
    matprops_ptr->bedfrict[1]=pcqbedfrict*PI/180.0;
    
    double doubleswap=sqrt(tan(matprops_ptr->bedfrict[1]));
    for(int imat=2;imat<=matprops_ptr->material_count;imat++) {
      matprops_ptr->bedfrict[imat] =matprops_ptr->bedfrict[1];
    }
  }

  /*************************************************************************/
  //to read in outline parameters here when it has been added


  return;
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
//this reads in the funky grid, ignoring the material properties at the 
//end of the file, those are now read from frict.data in Read_data()


void Read_grid(int myid, int numprocs, 
	       HashTable** NodeTable, HashTable** ElemTable, 
	       MatProps* matprops_ptr, OutLine* outline_ptr)
{
  int Node_Num, Elem_Num;
  
  int NODE_TABLE_SIZE=400000;

  //char  filename[14] = "lsh4800xx.inp";
  char  filename[14] = "funkyxxxx.inp";
  //unsigned min_key[KEYLENGTH];
  //unsigned max_key[KEYLENGTH];
  double doublekeyrange[2];
  double XRange[2];
  double YRange[2];
  unsigned key[KEYLENGTH];
  double coord[DIMENSION];
  double height;
  Node* NodeP;
  int i, j, k;
  
  // read in nodal data
  sprintf(filename,"funky%04d.inp",myid);
  FILE* fp; 
  
  fp=fopen_bin(filename,"r");
  if(!fp) {
    printf("Can't open file for %d \n", myid);
    exit(0);}
  
  int version, DoublesFromFloats;
  freadI(fp,&version);
  
  switch(version) {
  case 20061109:
    DoublesFromFloats=1;
    break;
  case 20061110:
    DoublesFromFloats=0;
    break;
  default:
    printf("Read_data() does not recognize binary funkyxxxx.inp version %d\n",
	   version);
    exit(1);
    break;}
  
  freadI(fp,&Node_Num);
  
  if(DoublesFromFloats) {
    for(i=0;i<KEYLENGTH;i++) freadF2D(fp,&(doublekeyrange[i]));
    
    freadF2D(fp,&(XRange[0]));  //min x
    freadF2D(fp,&(XRange[1]));  //max x
    freadF2D(fp,&(YRange[0]));  //min y
    freadF2D(fp,&(YRange[1]));} //max y
  else{
    for(i=0;i<KEYLENGTH;i++) freadD(fp,&(doublekeyrange[i]));

    freadD(fp,&(XRange[0]));  //min x
    freadD(fp,&(XRange[1]));  //max x
    freadD(fp,&(YRange[0]));  //min y
    freadD(fp,&(YRange[1]));} //max y
  
  double xminmax[2], yminmax[2];
  for(i=0;i<2;i++) {
    XRange[i] = XRange[i]/matprops_ptr->LENGTH_SCALE;
    xminmax[i]=XRange[i];
  }
  for(i=0;i<2;i++) {
    YRange[i] = YRange[i]/matprops_ptr->LENGTH_SCALE;
    yminmax[i]=YRange[i];
  }
  
  *NodeTable  = new HashTable(doublekeyrange,NODE_TABLE_SIZE, 2017, XRange, YRange,0); 

  for(i=0;i<Node_Num;i++) {
    for(j=0;j<KEYLENGTH;j++) freadU(fp,&(key[j]));
    
    if(DoublesFromFloats)
      for(j=0;j<DIMENSION;j++) freadF2D(fp,&(coord[j]));
    else
      for(j=0;j<DIMENSION;j++) freadD(fp,&(coord[j]));

    for(j=0;j<2;j++)
      coord[j] = coord[j]/matprops_ptr->LENGTH_SCALE;
    NodeP = new Node(key, coord, matprops_ptr);
    (*NodeTable)->add(key, NodeP);}   
  
  //done reading in node data
  //start reading in element data

  int EL_TABLE_SIZE=100000;
  
  //char  filename[14] = "lsh4800xx.inp";
  int material, elm_loc[2];
  unsigned opposite_brother[2];
  
  Element* Quad9P;
  float*   value=new float[2];
  BC*      baddress[4];
  void*    p;
  int*     assocp;/*--*/
  unsigned*    keyP;

  unsigned nodes[9][2];
  unsigned neigh[4][2];
  int      neighbor_proc[4];


  int temp2;
  int interflag;

  freadI(fp,&Elem_Num);  //--number of the elements assigned to the proc

  *ElemTable  = new HashTable(doublekeyrange,EL_TABLE_SIZE, 503, XRange, YRange,0); 
  for(int ielem=0;ielem<Elem_Num;ielem++) {      

    for(j=0; j<9; j++) 
      for(k=0; k<KEYLENGTH; k++)
	freadU(fp,&(nodes[j][k]));

    interflag = 0;//---switch for interface
    for(j=0; j<4; j++) {

      freadI(fp,&(neighbor_proc[j]));//--read the neighbor info

      if(neighbor_proc[j]!=-1)//--if there is neighbor(-1 means the edge is bound)
	{
	  if(neighbor_proc[j]!=myid) //--the neighbor belongs to other proc
	    interflag = 1;//--switch is used for avoiding nominating neighbor twice
	  
	  for(k=0; k<KEYLENGTH; k++)
	    freadU(fp,&(neigh[j][k]));//--read the left parts of the key
	}
      
      else//--there is no neighbor 
	for(k=0; k<KEYLENGTH; k++) neigh[j][k]=0;
    }
    
    BC* bcptr = 0;
    int bcf = 0;

    //.....the essential boundary conditions....
    
    for(j=0; j < 4; j++) {
      freadI(fp,&temp2);

      if(temp2!=-1)//--there is bound constraint
	{
	  if(!bcf)
	    bcptr = new BC();
	  bcptr->type[j] = 1; //--intialize type

	  /* "value" is a FLOAT so DON'T use freadD when DoublesFromFloats
	     is false (and obviously don't use freadD when it's true 
	     either) */
	  for(k=0; k<2; k++)		     
	    freadF(fp,&(bcptr->value[j][0][k]));//--j: edge number

	  bcf = 1;
	}  
      }
    
    //.....the natural boundary conditions.....
    for(j=0; j < 4; j++) {

      freadI(fp,&temp2);

      if(temp2!=-1)//--there is bound constraint
	{
	  if(!bcf)
	    bcptr = new BC();
	  if(bcptr->type[j]==0)
	    bcptr->type[j] = 2; //--intialize type
	  else
	    bcptr->type[j] =3; //--intialize type
	  
	  /* "value" is a FLOAT so DON'T use freadD when DoublesFromFloats
	     is false (and obviously don't use freadD when it's true 
	     either) */
	  for(k=0; k<2; k++)		     
	    freadF(fp,&(bcptr->value[j][1][k]));//--j: edge number

	  bcf = 1;
	}
    }

    freadI(fp,&(elm_loc[0]));
    freadI(fp,&(elm_loc[1]));
    freadU(fp,&(opposite_brother[0]));
    freadU(fp,&(opposite_brother[1]));
    freadI(fp,&material);

    double pile_height = 0.0;

    if(!bcf) bcptr = NULL; //--this element is not on the bound
    Quad9P=new Element(nodes, neigh, neighbor_proc, bcptr, material, elm_loc, 
		       pile_height, myid, opposite_brother);
    (*ElemTable)->add(nodes[8], Quad9P);
    
    Quad9P->find_positive_x_side(*NodeTable);
    Quad9P->calculate_dx(*NodeTable);
  }

  /************************************************************/
  /* need to change this so that smallest cell size is chosen */
  /* based on cube root of volume with failsafe for a minimum */
  /* of 2 cells across smallest pile/flux source axis instead */
  /* of basing it on just the smallest pile/flux source axis  */
  /* which is what we're doing now, will need to print out a  */
  /* warning that the calculation may be slow when the        */
  /* failsafe is activated                                    */
  /************************************************************/

  double dx[2]={*(Quad9P->get_dx()+0),*(Quad9P->get_dx()+1)}; 
  double DX=dx[0];
  if(dx[0]<dx[1]) DX=dx[1];


  REFINE_LEVEL=Quad9P->get_gen()+
    ceil(log(DX*
	     (matprops_ptr->number_of_cells_across_axis)/
	     (matprops_ptr->smallest_axis)
	     )/
	 log(2.0));


  if(REFINE_LEVEL<0) REFINE_LEVEL=0;
  printf("REFINE_LEVEL=%d\n",REFINE_LEVEL);



  //set the nodal type information
  Element *EmTemp;
  Node *NdTemp;
  int inode;
  int num_buck=(*ElemTable)->get_no_of_buckets();
  HashEntryPtr* buck =(*ElemTable)->getbucketptr();
  for(int i=0; i<num_buck; i++)
    if(*(buck+i)){
      
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){

	EmTemp=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;     
	assert(EmTemp);

	EmTemp->put_myprocess(myid);

	NdTemp=(Node*) (*NodeTable)->lookup(EmTemp->pass_key());
	assert(NdTemp);
	NdTemp->putinfo(BUBBLE);

	for(inode=0;inode<4;inode++) {
	  NdTemp=(Node*) 
	    (*NodeTable)->lookup(EmTemp->getNode()+inode*KEYLENGTH);
	  assert(NdTemp);
	  NdTemp->putinfo(CORNER);
	}

	for(inode=4;inode<8;inode++) {
	  NdTemp=(Node*) 
	    (*NodeTable)->lookup(EmTemp->getNode()+inode*KEYLENGTH);
	  assert(NdTemp);
	  NdTemp->putinfo(SIDE);
	}
      }
    }

  //initialize the flow "outline" map (maximum pileheight in every cell throught simulation is recorded)
  //printf("dx=%g dy=%g XRange={%g,%g} YRange={%g,%g}\n",*(Quad9P->get_dx()),*(Quad9P->get_dx()+1),xminmax[0],xminmax[1],yminmax[0],yminmax[1]);
#ifdef MAX_DEPTH_MAP
  outline_ptr->init(Quad9P->get_dx(),REFINE_LEVEL-Quad9P->get_gen(),
		    xminmax,yminmax);
#endif

  delete []value;

  //assert(0);

}


