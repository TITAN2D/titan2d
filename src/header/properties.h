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
 * $Id: properties.h 129 2007-06-07 19:54:28Z dkumar $ 
 */

#ifndef MAX_DEPTH_MAP
#define MAX_DEPTH_MAP
#endif

#ifndef PROPS
#define PROPS


#include<stdlib.h>
#include<stdio.h> //useful_lib.h has a function that return pointer to FILE
#include "../useful/useful_lib.h"
#include<time.h>  //for TimeProps
#include<math.h>
#include "../gisapi/GisApi.h"

//! LHS stands for Latin Hypercube Sampling, it is a constrained sampling method whose convergence can be much faster than monte carlo 
struct LHS_Props{
  
  //! the refinement number, the number of times each of the original intervals have been subdivided into 2 intervals in each direction 
  int refnum;     

  //! the unique identifier for each sample run 
  int runid;      
 
  //! this constructor initializes refnum and runid to -1
  LHS_Props() {
    refnum=runid=-1;
  } 

};


//! the StatProps structure holds statistics about the flow
struct StatProps{
  //note all means are mass/volume averages

  //! job number for monte carlo or lhs simulations
  int    runid;

  //! x coordinate of pile centroid
  double xcen;

  //! y coordinate of pile centroid
  double ycen;

  //! variance of location of pile material in the x direction
  double xvar;

  //! variance of location of pile material in the y direction
  double yvar;

  //! mean distance from the point (0,0), this was the mean pile starting point for the ordered reduced number of runs monte carlo method validation 
  double rmean;

  //! area covered by pile of thickness greater than cutoffheight 
  double area;

  //! mean speed
  double vmean;

  //! mean x velocity
  double vxmean;     

  //! mean y velocity
  double vymean;     

  //! mean slope in the direction of velocity a negative number	indicates the flow is heading uphill 
  double slopemean;  

  //! nondimensionalized mean speed
  double vstar;
  
  //! volume of ALL material, NOT used for other statistics
  double realvolume;

  //! STAT_VOL_FRAC*realvolume, all statistics computed in ../geoflow/stats.C are based on statvolume not realvolume
  double statvolume;   

  //! volume of material that have flown off the map
  double outflowvol;

  //! volume of material that has been eroded
  double erodedvol;

  //! volume of material that is currently deposited  
  double depositedvol;

  //! pile height of contour line that encloses statvolume
  double cutoffheight;

  //! an estimate of the radius of the pile
  double piler;

  //! current spatial maximum of pile height
  double hmax;

  //! current spatial maximum of speed
  double vmax;

  //! the integrated magnitude of acceleration due to internal friction force (based on realvolume not statvolume)
  double forceint;

  //! the same thing but for bed friction
  double forcebed;

  //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
  double heightifreach; 

  //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
  double xyifreach[2];

  //! check if a flow of height heightifreach reaches the point xyifreach and record the first time it does in timereached, if the value is -1 the flow has not reached the point with the specified height.
  double timereached;

  //! xyminmax holds the minimum and maximum x and y coordinates where the pile height is greater than hxyminmax
  double xyminmax[4];

  //! xyminmax holds the minimum and maximum x and y coordinates where the pile height is greater than hxyminmax
  double hxyminmax;     

  //! the latin hypercube sampling specific stats
  LHS_Props lhs; 

  //! the constructor initializes a few statistics
  StatProps() {
    timereached=-1.0; 
    xcen=ycen=xvar=yvar=rmean=area=vmean=vxmean=vymean=slopemean=vstar=0.0;
    realvolume=statvolume=outflowvol=erodedvol=depositedvol=cutoffheight=0.0;
    piler=hmax=vmax=forceint=forcebed=0.0;
    heightifreach=xyifreach[0]=xyifreach[1]=timereached=0.0;
    xyminmax[0]=xyminmax[1]=xyminmax[2]=xyminmax[3]=hxyminmax=0.0;
    lhs.refnum=lhs.runid=-1;
  }


  

};

//! the PileProps structure holds the pile properties read in in Read_data() so the pile can be placed at the proper locations shortly thereafter in init_piles() 
struct PileProps{
  //! the number of piles
  int numpiles;

  //!array holding pile height
  double *pileheight;

  //! array of volume -fractions
  double *vol_fract;

  //!array holding x coordinate of pile center
  double *xCen;

  //!array holding y coordinate of pile center
  double *yCen;

  //!array holding the major (x before rotation) radius
  double *majorrad;

  //!array holding the minor (y before rotation) radius
  double *minorrad;

  //!array holding the cosine of the rotation angle
  double *cosrot;

  //!array holding the sine of the rotation angle;
  double *sinrot;

  //!array holding the initial x speed of the pile
  double *initialVx;

  //!array holding the initial y speed of the pile
  double *initialVy;   

  //! this constuctor initializes the number of piles to zero.
  PileProps() {
    numpiles=0;
  }

  //! function allocates space for the pile data
  void allocpiles(int numpiles_in) {
    numpiles=numpiles_in;
    pileheight=CAllocD1(numpiles);
    vol_fract= CAllocD1(numpiles);
    xCen      =CAllocD1(numpiles);
    yCen      =CAllocD1(numpiles);
    majorrad  =CAllocD1(numpiles);
    minorrad  =CAllocD1(numpiles);
    cosrot    =CAllocD1(numpiles);
    sinrot    =CAllocD1(numpiles);
    initialVx =CAllocD1(numpiles);
    initialVy =CAllocD1(numpiles);

  }

  //! this function deallocates the dynamically out array members of the PileProps structure
  ~PileProps() {
    if(numpiles>0) {
      CDeAllocD1(pileheight);
      CDeAllocD1(vol_fract);
      CDeAllocD1(xCen);
      CDeAllocD1(yCen);
      CDeAllocD1(majorrad);
      CDeAllocD1(minorrad);
      CDeAllocD1(cosrot);
      CDeAllocD1(sinrot);
      CDeAllocD1(initialVx);
      CDeAllocD1(initialVy); 
    }     
  }

};



/*************************************************************************/
/* the gis map properties                                                */
/*************************************************************************/


//! this structure holds the path to and name of the GIS map and also a flag to say if there are any extra maps, such as a material properties map, associated with the DEM
struct MapNames{

  //! gis main directory
  char *gis_main;    

  //! gis sub directory
  char *gis_sub;     

  //! gis map set
  char *gis_mapset;  

  //! the actual gis map
  char *gis_map;

  //! extra maps: 0=none  +1=2^0=<gis_map>_Mat  MATerial map  +2^(bit#-1)=<gis_map>_Xxx  not yet used
  int extramaps;     
  
  //! this destructor calls MapNames::clear() which deallocates the dynamically allocated arrays in the MapNames structure.
  ~MapNames() { clear(); return; }

  //! this function allocates space for and assigns the information about the GIS map
  void assign(char *gis_main_in, char *gis_sub_in, char *gis_mapset_in, 
      char *gis_map_in, int extramaps_in) {
    gis_main=allocstrcpy(gis_main_in);
    gis_sub=allocstrcpy(gis_sub_in);
    gis_mapset=allocstrcpy(gis_mapset_in);
    gis_map=allocstrcpy(gis_map_in);
    extramaps=extramaps_in;
    return;
  }

  //! this function deallocates the dynamically allocated arrays in the MapNames structure.
  void clear() {
    free(gis_main);
    free(gis_sub);
    free(gis_mapset);
    free(gis_map);
    return;
  }
};

/**************************************************************************/
/* the time properties                                                    */
/**************************************************************************/

//! this structure holds all the information about time and timestepping
struct TimeProps{

  //! the maximum # of iterations (a.k.a. time steps) before the simulation ends 
  int    maxiter;


  //! the current number of iterations
  int    iter;       

  //! the maximum amount of time (in seconds) before the simulation ends 
  double maxtime;    

  //! the non-dimensional maxtime
  double ndmaxtime;

  //! the amount of time (in seconds) between subsequent outputs (with one exception... when the simulation ends one final output is performed and that one will be less than "timeoutput" seconds after the previous one 
  double timeoutput;


  //! the amount of time (in seconds) between subsequent saves (with one exception... when the simulation ends one final save is performed and that one will be less than "timeoutput" seconds after the previous one 
  double timesave;

  //! count of the number of times output has been done
  int    ioutput;

  //! count of the number of times saving has been done
  int    isave;

  //! the non-dimensional time at which the next output should occur
  double ndnextoutput; 

  //! the non-dimensional time at which the next save should occur 
  double ndnextsave; 

  //! the value used to nondimensionalize time
  double TIME_SCALE;

  //! the non-dimensional current time
  double time;

  //! the non-dimensional time step
  double dtime;

  //! velocity measure/Vslump used as a STOPPING CRITERIA,which is while it's stored under TimeProps (also members of MatProps are assigned once, are permanent for the run), see MatProps struct below for Vslump, see ../main/datread.C for initialization of Vslump     
  double vstarmax;

  //! wallclock time shortly after titan starts running
  time_t starttime;  

  //! this function initializes the time properties at the beginning of the simulation
  void inittime(int maxiterin, double maxtimein, double timeoutputin, 
		double timesavein, double TIME_SCALEin) {
    maxiter=maxiterin;
    maxtime=maxtimein;
    timeoutput=timeoutputin;
    timesave=timesavein;
    TIME_SCALE=TIME_SCALEin;
    ndmaxtime=maxtime/TIME_SCALE;
    ndnextoutput=timeoutput/TIME_SCALE;
    ndnextsave=timesave/TIME_SCALE;
    iter=0;
    ioutput=0;
    isave=0;
    time=0.0;
    dtime=0.0;
    vstarmax=0.0;
   }

  //! this function increments the time step, after, if neccessary, decreasing the time step to land evenly on the next time to output save or end the simulation
  void incrtime(double *dt) {
    // first reduce dt to hit output or end time "exactly"
    if(time+*dt>ndnextoutput) *dt=ndnextoutput-time;
    if(time+*dt>ndnextsave)   *dt=ndnextsave  -time;
    if(time+*dt>ndmaxtime)    *dt=ndmaxtime   -time;
    dtime=*dt;
    // then increment time
    time+=*dt; 
    iter++;}


  int  ifstart() {return(iter==0);} //! checks if it's before first time step
  int  iffirst() {return(iter==1);} //! checks if it's at first time step

  //! checks if the simulation should end due to running out of time, running out of timesteps or meeting a legacy "stopped" criteria
  int  ifend(double vstar) {
    if(vstar>vstarmax) vstarmax=vstar;
    return((time>=ndmaxtime)||(iter>maxiter)||((vstarmax>2.0)&&!(vstar>1.0)));}

  //! checks if the simulation has passed 1/10th of the maximum time allowed
  int  ifcheckstop() {return(time>ndmaxtime/10.0);}

  //! checks if the restart file should be saved now
  int  ifsave() {
    if(time>=ndnextsave) {
      isave++; //using isave eliminates roundoff
      ndnextsave=((isave+1)*timesave)/TIME_SCALE;
      return(1);}
    else return(0);}

  //! checks if the output files should be written now
  int  ifoutput() {
    if(time>=ndnextoutput) {
      ioutput++; //using ioutput eliminates roundoff
      ndnextoutput=((ioutput+1)*timeoutput)/TIME_SCALE;
      return(1);}
    else return(0);}

  //! chunk simulated time into hours minutes and seconds
  void chunktime(int *hours, int *minutes, double *seconds) {
    double dimtime=time*TIME_SCALE;
    *hours   = ((int) dimtime)/3600;
    *minutes = (((int) dimtime)%3600)/60;
    *seconds = dimtime -(double)(*hours*3600+*minutes*60);}

  //! return the simulated time in seconds
  double timesec() {return(time*TIME_SCALE);}
    
};

/*****************************************************************************/
//! this struct holds constants for material properties as well as other constants note that the material id tag (used as the indice for material properties... matname, bedfrict) as returned by Get_raster_id() (a GIS function call) starts from 1 and not from 0 so arrays must be one element larger
/*****************************************************************************/
struct MatProps{
  //! the "maximum" number of cells across the smallest pile/flux-source minor axis
  int number_of_cells_across_axis;

  //! the smallest minor radius
  double smallest_axis;

  //! the number of different materials on the map
  int     material_count; 

  //! the names of each material
  char  **matnames;

  //! phi_{int}, the internal friction angle (must be GREATER than the bedfriction angle)
  double  intfrict;

  //! tan(phi_{int}), tangent of the internal friction angle
  double  tanintfrict;

  //! phi_{bed}, the bed friction angle, must be LESS than the internal friction angle and should be greater than about 8 degrees, this minimum angle may change once Keith's local stopping criteria is enforced
  double *bedfrict;

  //! tan(phi_{bed}), tangent of the bed friction angle
  double *tanbedfrict;

  //! v_f, legacy not used
  double  porosity;

  //! density
  double  den_solid;

  //! fluid density
  double  den_fluid;

  //! fluid viscosity
  double  viscosity;

  //! terminal velocity of a single solid particle in fluid medium
  double  v_terminal;

  //! scaling value, ratio of HEIGHT_SCALE to LENGTH_SCALE
  double  epsilon;

  //! slope limiting stuff
  double  gamma;

  //! Flow type flag
  int     flow_type;

  //! length scaling factor
  double  LENGTH_SCALE;

  //! height scaling factor
  double  HEIGHT_SCALE;

  //! gravity scaling factor
  double  GRAVITY_SCALE;

  //! cells with flow below this height are neglected for purposes of calculating statistics
  double  MAX_NEGLIGIBLE_HEIGHT; 

  /*! 
   *  Used in Bin Yu's legacy global stopping criteria, this never worked right 
   *  except for initially cylindrical piles on a horizontal plane that were released 
   *  and slumped... there Bin Yu's global stopping criteria worked great 
   */
  double  Vslump; 

  //! to get the flow to start moving you have to decrease the friction angles, this variable is used to do that
  double  frict_tiny;     

  /*! this constructor allocates initial properties unfortunately the properties 
   *  aren't known at the time this is called so dummy values are fed in instead 
   *  for many if not all of these
   */
  MatProps(int material_countin, char **matnamesin, 
	   double intfrictin, double *bedfrictin,
           double porosityin, double muin, double rhoin, double rhofin,
	   double epsilonin, double gammain, double frict_tinyin,
	   double lscale, double hscale, double gscale) {

    material_count=material_countin;

    if(material_count>0) {
      //dynamic memory allocation... see useful_lib.C for CAlloc?#()
      matnames=(char **) malloc((material_count+1)*sizeof(char *));
      bedfrict=CAllocD1(material_count+1);
      tanbedfrict=CAllocD1(material_count+1);

      for(int imat=1; imat<=material_count; imat++) {
	matnames[imat]=allocstrcpy(matnamesin[imat]);  //see useful_lib.C
	bedfrict[imat] = bedfrictin[imat];
	tanbedfrict[imat] = tan(bedfrict[imat]);}}

    intfrict=intfrictin;
    tanintfrict=tan(intfrict);
    porosity = porosityin;
    viscosity = muin;
    den_solid = rhoin;
    den_fluid = rhofin;
    epsilon = epsilonin;
    gamma = gammain;
    frict_tiny = frict_tinyin;
    LENGTH_SCALE = lscale;
    HEIGHT_SCALE = hscale;
    GRAVITY_SCALE = gscale;
  }

  //! this destructor deallocates the arrays of bed friction angles and their tangents
  ~MatProps() {
    if(material_count>0) {
      CDeAllocD1(bedfrict);
      CDeAllocD1(tanbedfrict);
      for(int imat=1; imat<=material_count; imat++) free(matnames[imat]);
      free(matnames);
    }
    return;
  }
};


//! the OutLine Structure holds the maximum throughout time flow depth at every spatial point
struct OutLine{
  // the least squares interpolated height didn't work for some reason
  // that is didn't produce a meaningful pileheight contour map so it has
  // been commented out.  that is everything to deal with pileheight2
//! number of cells in the x direction on the map
  int Nx;

  //! number of cells in the y direction on the map 
  int Ny;

  //! length of a cell in the x direction
  double dx;

  //! length of a cell in the y direction
  double dy; 

  //! min and max x coordinate on the map
  double xminmax[2];

  //! min and max y coordinate on the map
  double yminmax[2];

  //! dynamically allocated 2 dimensional array holding the maximum throughout time pileheight at every point
  double **pileheight;
  //double **pileheight2;

  //! this is the OutLine constructor it initializes the number of cells to zero
  OutLine() { Nx=Ny=0; return; }

  //! this is the OutLine it deallocates the 2 dimensional array holding maximum throughout time pileheight in every cell on the map
  ~OutLine() {
    if((Nx>0)&&(Ny>0)) CDeAllocD2(pileheight);
    return;
  }

  //! this function initializes the OutLine map/2-dimensional array 
  void init(double *dxy, int power, double *XRange, double *YRange) {
    int ix,iy;

    if(power<0) power=0;

    dx=dxy[0]/pow(2.0,power);
    dy=dxy[1]/pow(2.0,power);
    //printf("dx=%g dy=%g  XRange={%g,%g} YRange={%g,%g}\n",dx,dy,XRange[0],XRange[1],YRange[0],YRange[1]);

    xminmax[0]=XRange[0]; xminmax[1]=XRange[1];
    yminmax[0]=YRange[0]; yminmax[1]=YRange[1];

    Nx=(int)( (XRange[1]-XRange[0])/dx +0.5); //round to nearest integer
    Ny=(int)( (YRange[1]-YRange[0])/dy +0.5); //round to nearest integer

    while(Nx*Ny>1024*1024){
      dx*=2.0;
      dy*=2.0;

      Nx=(int)( (XRange[1]-XRange[0])/dx +0.5); //round to nearest integer
      Ny=(int)( (YRange[1]-YRange[0])/dy +0.5); //round to nearest integer
    }
    printf("Outline init: Nx=%d Ny=%d Nx*Ny=%d\n",Nx,Ny,Nx*Ny);


    pileheight =CAllocD2(Ny,Nx);
    //pileheight2=CAllocD2(Ny,Nx);
    for(iy=0;iy<Ny;iy++)
      for(ix=0;ix<Nx;ix++) {
	pileheight[ iy][ix]=0.0;
	//pileheight2[iy][ix]=0.0;
      }      
    return;
  }

  //! this function reinitializes the OutLine map/2-dimensional array during restart
  void init2(double *dxy, double *XRange, double *YRange) {
    int ix,iy;

    dx=dxy[0];
    dy=dxy[1];

    xminmax[0]=XRange[0]; xminmax[1]=XRange[1];
    yminmax[0]=YRange[0]; yminmax[1]=YRange[1];

    Nx=(int)( (XRange[1]-XRange[0])/dx +0.5); //round to nearest integer
    Ny=(int)( (YRange[1]-YRange[0])/dy +0.5); //round to nearest integer

    pileheight =CAllocD2(Ny,Nx);
    //pileheight2=CAllocD2(Ny,Nx);
    for(iy=0;iy<Ny;iy++)
      for(ix=0;ix<Nx;ix++) {
	pileheight[ iy][ix]=0.0;
	//pileheight2[iy][ix]=0.0;
      }      

    return;
  }


  //! this function updates the maximum throughout time pileheight in every cell covered by an arbitrary element
  void update(double xstart, double xstop, double ystart, double ystop, double height, double h2[6]) {
    int ixstart=(int)((xstart-xminmax[0])/dx+0.5);
    int ixstop=(int)((xstop-xminmax[0])/dx+0.5);
    int iystart=(int)((ystart-yminmax[0])/dy+0.5);
    int iystop=(int)((ystop-yminmax[0])/dy+0.5);

    /*    double x,y,height2;
    double xc=0.5*(xstart+xstop);
    double yc=0.5*(ystart+ystop);
    */  
    if(ixstart<0) ixstart=0;
    if(ixstop==ixstart) {
      ixstart=(int) ((xstart-xminmax[0])/dx);
      ixstop=ixstart+1;
    }
    if(ixstop>Nx) ixstop=Nx;


    if(iystart<0) iystart=0;
    if(iystop==iystart) {
      iystart=(int) ((ystart-yminmax[0])/dy);
      iystop=iystart+1;
    }
    if(iystop>Ny) iystop=Ny;

    for(int iy=iystart;iy<iystop;iy++)
      for(int ix=ixstart;ix<ixstop;ix++) {
	/*	x=(ix+0.5)*dx-xc;
	y=(iy+0.5)*dy-yc;
	height2=h2[0]+h2[1]*x+h2[2]*y+h2[3]*x*y+h2[4]*x*x+h2[5]*y*y;
	if(height2>pileheight2[iy][ix]) pileheight2[iy][ix]=height2; */
	if(height >pileheight[ iy][ix]) pileheight[ iy][ix]=height;
      }
    return;
  }

  //! this function outputs the maximum throughout time map of pileheights to the file pileheightrecord.xxxxxx
  void output(MatProps* matprops_ptr, StatProps* statprops_ptr) {
    int ix,iy;
    char filename[256];
    sprintf(filename,"pileheightrecord.%06d",statprops_ptr->runid);
    FILE *fp=fopen(filename,"w");

    //FILE *fp=fopen("outline.pileheight","w");
    fprintf(fp,"Nx=%d: X={%20.14g,%20.14g}\nNy=%d: Y={%20.14g,%20.14g}\nPileheight=\n",
                Nx,xminmax[0]*matprops_ptr->LENGTH_SCALE,
                xminmax[1]*matprops_ptr->LENGTH_SCALE,Ny,
                yminmax[0]*matprops_ptr->LENGTH_SCALE,
                yminmax[1]*matprops_ptr->LENGTH_SCALE);
    for(iy=0; iy<Ny; iy++) 
    {
      for(ix=0; ix<Nx-1; ix++) fprintf(fp,"%g ",pileheight[iy][ix]*matprops_ptr->HEIGHT_SCALE);
      fprintf(fp,"%g\n",pileheight[iy][ix]*matprops_ptr->HEIGHT_SCALE);
    }
    fclose(fp);

    sprintf(filename,"elevation.grid");
    fp=fopen(filename,"w");
    fprintf(fp,"Nx=%d: X={%20.14g,%20.14g}\nNy=%d: Y={%20.14g,%20.14g}\nPileheight=\n",Nx,xminmax[0]*matprops_ptr->LENGTH_SCALE,xminmax[1]*matprops_ptr->LENGTH_SCALE,Ny,yminmax[0]*matprops_ptr->LENGTH_SCALE,yminmax[1]*matprops_ptr->LENGTH_SCALE);
    double yy, xx, res=dx+dy, elevation;
    int ierr;
    for(iy=0; iy<Ny; iy++) {
      yy=((iy+0.5)*dy+yminmax[0])*matprops_ptr->LENGTH_SCALE;
      for(ix=0; ix<Nx-1; ix++) {
	xx=((ix+0.5)*dx+xminmax[0])*matprops_ptr->LENGTH_SCALE;
	ierr = Get_elevation(res, xx, yy, &elevation);
	fprintf(fp,"%g ",elevation);
      }
      fprintf(fp,"%g\n",pileheight[iy][ix]*matprops_ptr->HEIGHT_SCALE);
    }
    fclose(fp);


    /*    fp=fopen("outline.pileheight2","w");
    fprintf(fp,"Nx=%d: X={%20.14g,%20.14g}\nNy=%d: Y={%20.14g,%20.14g}\nPileheight=\n",Nx,xminmax[0]*matprops_ptr->LENGTH_SCALE,xminmax[1]*matprops_ptr->LENGTH_SCALE,Ny,yminmax[0]*matprops_ptr->LENGTH_SCALE,yminmax[1]*matprops_ptr->LENGTH_SCALE);
    for(iy=0; iy<Ny; iy++) {
      for(ix=0; ix<Nx-1; ix++) fprintf(fp,"%g ",pileheight2[iy][ix]*matprops_ptr->HEIGHT_SCALE);
      fprintf(fp,"%g\n",pileheight2[iy][ix]*matprops_ptr->HEIGHT_SCALE);
    }
    fclose(fp);
    */
    return;
    



  }
  
  //! this function reads in the previous map of maximum throughout time pileheight stored in the file pileheightrecord.xxxxxx during restart
  void reload(MatProps* matprops_ptr, StatProps* statprops_ptr) {
    int ix,iy;
    char filename[256];
    sprintf(filename,"pileheightrecord.%06d",statprops_ptr->runid);
    FILE *fp=fopen(filename,"r");

    if(fp==NULL)
      printf("pileheightrecord.%06d can not be found.\nRestarting from zero instead\n",statprops_ptr->runid);
    else{
      
      int Nxtemp, Nytemp;
      double xminmaxtemp[2], yminmaxtemp[2];
      fscanf(fp,"Nx=%d: X={%lf,%lf}\nNy=%d: Y={%lf,%lf}\nPileheight=\n",&Nxtemp,(xminmaxtemp+0),(xminmaxtemp+1),&Nytemp,(yminmaxtemp+0),(yminmaxtemp+1));

      if((Nxtemp==Nx)&&
	 (fabs(xminmaxtemp[0]/matprops_ptr->LENGTH_SCALE-xminmax[0])<=fabs(xminmax[0])/10000000000.0)&&
	 (fabs(xminmaxtemp[1]/matprops_ptr->LENGTH_SCALE-xminmax[1])<=fabs(xminmax[1])/10000000000.0)&&
	 (Nytemp==Ny)&&
	 (fabs(yminmaxtemp[0]/matprops_ptr->LENGTH_SCALE-yminmax[0])<=fabs(yminmax[0])/10000000000.0)&&
	 (fabs(yminmaxtemp[1]/matprops_ptr->LENGTH_SCALE-yminmax[1])<=fabs(yminmax[1])/10000000000.0)) {

	for(iy=0; iy<Ny; iy++) 
	  for(ix=0; ix<Nx; ix++) {
	    fscanf(fp,"%lf",pileheight[iy]+ix);
	    pileheight[iy][ix]/=matprops_ptr->HEIGHT_SCALE;
	  }
      }
      else {
	printf("the pileheightrecord.%06d that is present does not match the restart file.\nRestarting from zero instead\n",statprops_ptr->runid);
	printf("Nx=%d Nxtemp=%d\n",Nx,Nxtemp);
	printf("Ny=%d Nytemp=%d\n",Ny,Nytemp);
	printf("xmin=%20.14g xmintemp=%20.14g  xmax=%20.14g, xmaxtemp=%20.14g\n",
	       xminmax[0]*matprops_ptr->LENGTH_SCALE,xminmaxtemp[0],
	       xminmax[1]*matprops_ptr->LENGTH_SCALE,xminmaxtemp[1]);
	printf("ymin=%20.14g ymintemp=%20.14g  ymax=%20.14g, ymaxtemp=%20.14g\n",
	       yminmax[0]*matprops_ptr->LENGTH_SCALE,yminmaxtemp[0],
	       yminmax[1]*matprops_ptr->LENGTH_SCALE,yminmaxtemp[1]);
	exit(0);
      }
	       
      fclose(fp);
    }

    return;
  }

  //! this function deallocates the 2 dimensional array of maximum throughout time pileheights
  void dealloc() {
    CDeAllocD2(pileheight);
    //CDeAllocD2(pileheight2);
    return;
  }

};


/**********************************************************************/
//! this structure is for the calculation of volume that flows through user specified discharge planes.  The sign of the volume indicates which direction the flow went and follows the right hand rule convention.  velocity cross (point b-point a) is the sign of the flow through planes.  This means if you surround the only pile, specify the points defining the discharge planes in counter clockwise order, the flow "out of the box" will be positive.  if you specify the points in clockwise order flow "out of the box" will be negative.
struct DISCHARGE {

  //! the number of discharge planes
  int num_planes;   

  //! the discharge planes are lines (with planes normal to the surface intersecting the surface passing through the lines), this holds a lot of information associated with each planes, a lot of precomputed quantities to make updating the flux through the planes fast.
  double **planes;  

  //! this constructor initializes the number of planes to zero
  DISCHARGE() {
    num_planes=0;
    return;
  }

  //! this destructor deallocates the planes information
  ~DISCHARGE() {
    if(num_planes>0) CDeAllocD2(planes);
    return;
  }

  //reinitialized in load_run()
  //! this function initializes the planes information (to zero flux through the planes) and precomputes a number of quantities to make updating the flux through the planes fast
  void init(int num_planes_in, double **planes_in) {
    double planestemp[4];
    
    //printf("num_planes_in=%d   num_planes=%d\n",num_planes_in,num_planes);
    num_planes=num_planes_in;
    //printf("num_planes_in=%d   num_planes=%d\n",num_planes_in,num_planes);
    if(num_planes>0) {
      planes=CAllocD2(num_planes,10);    

      for(int iplane=0; iplane<num_planes; iplane++) {
	/*      if((planes_in[iplane][1]<planes_in[iplane][0])||
	 ((planes_in[iplane][1]==planes_in[iplane][0])&&
	 (planes_in[iplane][3]<planes_in[iplane][2]))) {
	 planestemp[0]=planes_in[iplane][1];
	 planestemp[1]=planes_in[iplane][0];
	 planestemp[2]=planes_in[iplane][3];
	 planestemp[3]=planes_in[iplane][2];}
	 else{ */
	planestemp[0]=planes_in[iplane][0];
	planestemp[1]=planes_in[iplane][1];
	planestemp[2]=planes_in[iplane][2];
	planestemp[3]=planes_in[iplane][3];
	//}

	planes[iplane][0]=planestemp[0]; //xa
	planes[iplane][1]=planestemp[1]; //xb
	planes[iplane][2]=planestemp[2]; //ya
	planes[iplane][3]=planestemp[3]; //yb
	planes[iplane][4]=planes[iplane][1]-planes[iplane][0]; //xb-xa
	planes[iplane][5]=planes[iplane][3]-planes[iplane][2]; //yb-ya
	planes[iplane][6]= //(xb-xa)^2+(yb-ya)^2
	  planes[iplane][4]*planes[iplane][4]+
	  planes[iplane][5]*planes[iplane][5];
	planes[iplane][7]= //ya*(xb-xa)-xa*(yb-ya)
	  planes[iplane][3]*planes[iplane][4]-
	  planes[iplane][1]*planes[iplane][5];
	planes[iplane][8]=
	  ((fabs(planes[iplane][4])+fabs(planes[iplane][5]))*
	   (fabs(planes[iplane][4])+fabs(planes[iplane][5])))/
	  planes[iplane][6];
	
	planes[iplane][9]=0.0; //discharge through the plane
	
	//printf("plane %d: (%16.10g,%16.10g) (%16.10g,%16.10g)\n",iplane,planes[iplane][0],planes[iplane][2],planes[iplane][1],planes[iplane][3]);
	
	
      }
    }
    return;
  }
 
  //! this function updates the flux through the discharge planes
  void update(double nodes[9][2], double *statevars, double dt) {
    //FILE* fp=fopen("dischargedebug","a");

    double doubleswap1, doubleswap2, doubleswap3;
    double nearestpoint[2], intersectpoint[2][2];
    int iplane, iintersect, icorner1, icorner2;
    double dxnode[4], dynode[4];
    double dist2nearest2;
    double halfsidelength=//assumes grid is not slanted
      (fabs(nodes[6][0]-nodes[4][0])+fabs(nodes[6][1]-nodes[4][1])+
       fabs(nodes[7][0]-nodes[5][0])+fabs(nodes[7][1]-nodes[5][1]))*0.25;
    double halfsidelength2=halfsidelength*halfsidelength;
    double err=halfsidelength/pow(2.0,19.0);

    for(icorner1=0;icorner1<4;icorner1++) {
      icorner2=(icorner1+1)%4;
      dxnode[icorner1]=nodes[icorner2][0]-nodes[icorner1][0];
      dynode[icorner1]=nodes[icorner2][1]-nodes[icorner1][1];
    }
    


    for(iplane=0;iplane<num_planes;iplane++) {
      /* nearest x & y along line (not necessaryly on line segment) a and b
	 and the elements center node
	x=(-(yb-ya)*(ya*(xb-xa)-xa*(yb-ya))+(xb-xa)*(y8*(yb-ya)+x8*(xb-xa)))/
	  ((xb-xa)^2+(yb-ya)^2);
	y=( (xb-xa)*(ya*(xb-xa)-xa*(yb-ya))+(yb-ya)*(y8*(yb-ya)+x8*(xb-xa)))/
	  ((xb-xa)^2+(yb-ya)^2);
      */	
      doubleswap1=nodes[8][1]*planes[iplane][5]+nodes[8][0]*planes[iplane][4];
      //         =y8         *(yb-ya)          +x8         *(xb-xa)
      nearestpoint[0]=
	(-planes[iplane][5]*planes[iplane][7]
	 +planes[iplane][4]*doubleswap1)/planes[iplane][6];
      nearestpoint[1]=
	(+planes[iplane][4]*planes[iplane][7]
	 +planes[iplane][5]*doubleswap1)/planes[iplane][6];
      
      dist2nearest2=
	(nodes[8][0]-nearestpoint[0])*(nodes[8][0]-nearestpoint[0])+
	(nodes[8][1]-nearestpoint[1])*(nodes[8][1]-nearestpoint[1]);

      //check if line interesects with cell (not counting a single corner)
      if(dist2nearest2<halfsidelength2*planes[iplane][8]) {
	/* fprintf(fp,"nodecoord=(%8f,%8f) nearestpoint=(%8f,%8f)",
		nodes[8][0],nodes[8][1],nearestpoint[0],nearestpoint[1]);
	*/

	//it does intersect... if not terminated by endpoints
	
	iintersect=0; //stop when you have 2 unique "intersections" 
	//(the end of a discharge plane line segment is considered to 
	//be an "intersection")
	for(icorner1=0;icorner1<4;icorner1++) {
	  doubleswap1=
	    nodes[icorner1][1]*dxnode[icorner1]-
	    nodes[icorner1][0]*dynode[icorner1];
	  doubleswap2=	      
	    dxnode[icorner1]*planes[iplane][5]-
	    dynode[icorner1]*planes[iplane][4];

	  if((doubleswap2<0)||(doubleswap2>0)) {

	    intersectpoint[iintersect][0]=
	      (planes[iplane][4]*doubleswap1-
	       dxnode[icorner1]*planes[iplane][7])/
	      doubleswap2;
	    intersectpoint[iintersect][1]=
	      (planes[iplane][5]*doubleswap1-
	       dynode[icorner1]*planes[iplane][7])/
	      doubleswap2;
	    
	    int ifprint=0;
	    if(((intersectpoint[iintersect][0]<planes[iplane][0]-err)&&
		(planes[iplane][0]<planes[iplane][1])) ||
	       ((intersectpoint[iintersect][0]>planes[iplane][0]+err)&&
		(planes[iplane][0]>planes[iplane][1])) ||
	       ((intersectpoint[iintersect][1]<planes[iplane][2]-err)&&
		(planes[iplane][2]<planes[iplane][3])) ||
	       ((intersectpoint[iintersect][1]>planes[iplane][2]+err)&&
		(planes[iplane][2]>planes[iplane][3]))) {
	      intersectpoint[iintersect][0]=planes[iplane][0];
	      intersectpoint[iintersect][1]=planes[iplane][2];
	    }
	    else if(((intersectpoint[iintersect][0]<planes[iplane][1]-err)&&
		     (planes[iplane][1]<planes[iplane][0])) ||
		    ((intersectpoint[iintersect][0]>planes[iplane][1]+err)&&
		     (planes[iplane][1]>planes[iplane][0])) ||
		    ((intersectpoint[iintersect][1]<planes[iplane][3]-err)&&
		     (planes[iplane][3]<planes[iplane][2])) ||
		    ((intersectpoint[iintersect][1]>planes[iplane][3]+err)&&
		     (planes[iplane][3]>planes[iplane][2]))) {
	      intersectpoint[iintersect][0]=planes[iplane][1];
	      intersectpoint[iintersect][1]=planes[iplane][3];
	    }

	    if(((((intersectpoint[iintersect][0]-intersectpoint[0][0])*
		  (intersectpoint[iintersect][0]-intersectpoint[0][0])+
		  (intersectpoint[iintersect][1]-intersectpoint[0][1])*
		  (intersectpoint[iintersect][1]-intersectpoint[0][1]))>
		 (halfsidelength/512.0))||
		(iintersect==0))&&(doubleswap2!=0.0)) iintersect++;
	    
	    if(iintersect==2) {

	      if((intersectpoint[1][0]-intersectpoint[0][0])*
		 (planes[iplane][1]-planes[iplane][0])+
		 (intersectpoint[1][1]-intersectpoint[0][1])*
		 (planes[iplane][3]-planes[iplane][2]) < 0.0) {
		doubleswap3=intersectpoint[0][0];
		intersectpoint[0][0]=intersectpoint[1][0];
		intersectpoint[1][0]=doubleswap3;
		doubleswap3=intersectpoint[0][1];
		intersectpoint[0][1]=intersectpoint[1][1];
		intersectpoint[1][1]=doubleswap3;
	      }

	      break; //we've got both intersection points
	    }
	  }
	} // for(icorner1=0;icorner1<4;icorner1++)	 
	
	if(iintersect==1) {
	  //a plane end point is in this cell
	  
	  doubleswap1= //dist from center node to 1st plane endpoint
	    (nodes[8][0]-planes[iplane][0])*(nodes[8][0]-planes[iplane][0])+
	    (nodes[8][1]-planes[iplane][2])*(nodes[8][1]-planes[iplane][2]);
	  
	  doubleswap2= //dist from center node to 2nd plane endpoint
	    (nodes[8][0]-planes[iplane][1])*(nodes[8][0]-planes[iplane][1])+
	    (nodes[8][1]-planes[iplane][3])*(nodes[8][1]-planes[iplane][3]);

	  if(doubleswap1<=doubleswap2) { //it's the 1st end point
	    intersectpoint[1][0]=planes[iplane][0];
	    intersectpoint[1][1]=planes[iplane][2];}
	  else { //it's the 2nd end point
	    intersectpoint[1][0]=planes[iplane][1];
	    intersectpoint[1][1]=planes[iplane][3];}
	  
	  iintersect=2; //we consider plane end point to be an "intersection"
	}

	if(iintersect==2) {

	  /*
	  if((intersectpoint[0][0]>intersectpoint[1][0])||
	     ((intersectpoint[0][0]==intersectpoint[1][0])&&
	      (intersectpoint[0][1]>intersectpoint[1][1]))) {
	    //need to reorder the intersection points to have correct sign 
	    //of discharge
	    
	    //swap x
	    doubleswap1=intersectpoint[0][0];
	    intersectpoint[0][0]=intersectpoint[1][0];
	    intersectpoint[1][0]=doubleswap1;
	    
	    //swap y
	    doubleswap1=intersectpoint[0][1];
	    intersectpoint[0][1]=intersectpoint[1][1];
	    intersectpoint[1][1]=doubleswap1;
	  }
	  */

	  //discharge += dt*(hVx*dy-hVy*dx)
	  planes[iplane][9]+=dt*
	    (statevars[1]*(intersectpoint[1][1]-intersectpoint[0][1])-
	     statevars[2]*(intersectpoint[1][0]-intersectpoint[0][0]));
	} // if(iintersect==2)

      } //if(halfsidelength*(absdy+absdx)>=absdx*absdx+absdy*absdy)

    } //for(iplane=0;iplane<num_planes;iplane++)

    return;
  }

  //! this function deallocates the array holding the information about the discharge planes
  void dealloc() {
    num_planes=0;
    CDeAllocD2(planes);
    return;
  }
};

//! The FluxProps Structure holds all the data about extrusion flux sources (material flowing out of the ground) they can become active and later deactivate at any time during the simulation.  There must be at least 1 initial pile or one flux source that is active at time zero, otherwise the timestep will be set to zero and the simulation will never advance. 
struct FluxProps
{
  //! number of extrusion flux sources
  int no_of_sources;

  //! array holding the influx rates for all extrusion flux sources, if calculation is not scaled this quantity has units of [m/s]
  double *influx;


  //! array holding the activation time for all extrusion flux sources
  double *start_time;

  //! array holding the deactivation time for all extrusion flux sources
  double *end_time;

  //! array holding the x coordinate of the center of each the extrusion flux source
  double *xCen;

  //! array holding the y coordinate of the center of each the extrusion flux source
  double *yCen;

  //! array holding the radius along the major axis of each the extrusion flux source, the 2D pile shape is elliptical, in 3D it's a paraboloid
  double *majorrad;

  //! array holding the radius along the minor axis of each the extrusion flux source, the 2D pile shape is elliptical, in 3D it's a paraboloid
  double *minorrad;

  //! array holding the cosine of the rotation angle of each the extrusion flux source, the 2D pile shape is elliptical, and the ellipse's major axis does not have to be aligned with the x axis.
  double *cosrot;

  //! array holding the sine of the rotation angle of each the extrusion flux source, the 2D pile shape is elliptical, and the ellipse's major axis does not have to be aligned with the x axis.
  double *sinrot;

  //! array holding the x component of velocity of material extruding from each flux source
  double *xVel;

  //! array holding the y component of velocity of material extruding from each flux source
  double *yVel;

  FluxProps(){
    no_of_sources=0;
  }
  //! this function allocates space for all the extrusion rate fluxes, Dinesh Kumar added it, Keith modified it slightly
  void allocsrcs (int nsrcs)
  {
    no_of_sources=nsrcs;
    influx  =new double[nsrcs];
    xCen    =new double[nsrcs];
    yCen    =new double[nsrcs];
    majorrad=new double[nsrcs];
    minorrad=new double[nsrcs];
    cosrot  =new double[nsrcs];
    sinrot  =new double[nsrcs];
    start_time  =new double[nsrcs];
    end_time    =new double[nsrcs];
    xVel    =new double[nsrcs];
    yVel    =new double[nsrcs];
  }

  //! this function returns 1 if any flux sources become active during the current timestep, this is used to trigger "initial adaptation" of the flux source area, Keith wrote this function
  int IfAnyStart(TimeProps *timeprops_ptr) {

    for(int isrc=0; isrc<no_of_sources; isrc++)
      if(((timeprops_ptr->time-timeprops_ptr->dtime<=start_time[isrc])&&
	  (start_time[isrc]<timeprops_ptr->time))||
	 ((timeprops_ptr->iter==0)&&
	  (start_time[isrc]==0.0))
	 )
	return(1);
 
    return(0);
  }

  //! this function returns the maximum of all currently active extrusion fluxes, Keith wrote this function
  double MaxInfluxNow(MatProps *matprops_ptr, TimeProps *timeprops_ptr) {
    double tempinflux, maxinflux=0.0;
    

    for(int isrc=0; isrc<no_of_sources; isrc++)
      if(((start_time[isrc]<=timeprops_ptr->time)&&
	  (timeprops_ptr->time<=end_time[isrc]))||
	 ((timeprops_ptr->iter==0)&&
	  (start_time[isrc]==0))
	 ){
	tempinflux=sqrt(influx[isrc]*influx[isrc]+
			(xVel[isrc]*xVel[isrc]+
			 yVel[isrc]*yVel[isrc])/
			((matprops_ptr->epsilon)*(matprops_ptr->epsilon)));
	if(tempinflux>maxinflux)
	  maxinflux=tempinflux;
      }
       
    return(maxinflux);
  }

};

#endif
