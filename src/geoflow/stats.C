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
 * $Id: stats.C 232 2012-03-27 00:33:41Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"


/* STAT_VOL_FRAC is only here temporarily until a good value for it
   is found (because changing geoflow.h requires recompiling all of
   titan), it will then be moved to ../header/geoflow.h */
#define STAT_VOL_FRAC 0.95


void calc_stats(HashTable* El_Table, HashTable* NodeTable, int myid, MatProps* matprops, 
		TimeProps* timeprops, StatProps* statprops, 
		DISCHARGE* discharge, double d_time)
{
  int i, iproc;
  double area=0.0, max_height=0.0;
  double cutoffvolume;     /* the desired volume of material to take 
			      statistics from, this is portion is the one 
			      with the largest heights, for example 
                              sampling the 95% (by volume) of the pile with 
			      the largest heights, the fraction is set by
			      STAT_VOL_FRAC define statement  */
  double cutoffheight=0.0; /* a pile height criteria equivalent to the 
			      cutoffvolume, this criteria is found at each 
			      iteration */
  double statvolume;       /* volume where pile thickness >= cutoffheight 
			      statvolume >= cutoffvolume */
  double realvolume=0.0;   /* the total volume on this processor or across 
			      processors */
  double slopevolume=0.0;  /* volume used to determine the average slope 
			      in the direction of velocity, 
			      slopevolume<=statvolume because can't count
			      cells with zero velocity */
  double testpointheight=statprops->heightifreach;
  double testpointx = statprops->xyifreach[0];
  double testpointy = statprops->xyifreach[1];
  double testpointdist2;
  double testpointmindist2;
  testpointmindist2 =pow(2.0,30.0);//HUGE_VAL;
  int testpointreach=0;

  double testvolume=0.0;
  double slope_ave=0.0;
  double v_max=0.0, v_ave=0.0, vx_ave=0.0, vy_ave=0.0, g_ave;
  double xC=0.0, yC=0.0, rC=0.0, piler2=0.0;
  double xVar=0.0, yVar=0.0;
  //assume that mean starting location is at (x,y) = (1,1)
  double xCen=1.2; //0; //1.0/(matprops->LENGTH_SCALE);
  double yCen=0.3; //0; //1.0/(matprops->LENGTH_SCALE);
  double dVol, dA;
  double min_height=matprops->MAX_NEGLIGIBLE_HEIGHT;
  double temp;
  int numproc;
  double xyminmax[4];
  for(i=0;i<4;i++) xyminmax[i]=HUGE_VAL;

  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  /* need to allocate space to store nonzero pile heights and 
     volumes, to do that we first have to count the number of 
     nonzero elements */

  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  int num_nonzero_elem=0, *all_num_nonzero_elem;
  for(i=0; i<num_buck; i++)
    if(*(buck+i)){

      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){

	Element* Curr_El=(Element*)(currentPtr->value);
	if((Curr_El->get_adapted_flag()>0)&&
	   (myid==Curr_El->get_myprocess()))
	  if(*(Curr_El->get_state_vars())>GEOFLOW_TINY)
	    num_nonzero_elem++;
	    
	currentPtr=currentPtr->next;      	    
      }
    }
  /******************************************************************/
  /*** the rewrite involves a sort which increases the cost/time  ***/
  /*** significantly and does not appear to provide a significant ***/
  /*** increase in accuracy of the statistics (it probably makes  ***/
  /*** them LESS precise/repeatable actually) so the section of   ***/
  /*** code has been disabled by the following preprocessor       ***/
  /*** directive.                                                 ***/
  /******************************************************************/
#ifdef NONONONO
  printf("NONONONO\n");
  /* now we need to tell processor 0 how many there are, so it
     can allocate space for when it combine their contributions
     to come up with a cut off height */

  MPI_Request request, *proc0request;
  MPI_Status status, *proc0status;
 
  //fprintf(fp,"calc_stats() 2\n");  fflush(fp);
  if(numproc>1){
    if(myid==0){
      proc0request=(MPI_Request *) calloc(numproc,sizeof(MPI_Request));
      proc0status=(MPI_Status *) calloc(numproc,sizeof(MPI_Status));
      all_num_nonzero_elem=CAllocI1(numproc);

      all_num_nonzero_elem[0]=num_nonzero_elem;
      for(iproc=1;iproc<numproc;iproc++)
	MPI_Irecv(&all_num_nonzero_elem[iproc],1,MPI_INT,iproc,1,MPI_COMM_WORLD,&proc0request[iproc]);
    }
    else
      MPI_Isend(&num_nonzero_elem,1,MPI_INT,0,1,MPI_COMM_WORLD,&request);
  }
  //fprintf(fp,"calc_stats() 3\n");  fflush(fp);

  /* now each processor generates a list of pile heights and volumes
     for each nonzero element sorted in descending order */


  double **myprochvol=CAllocD2(num_nonzero_elem,2);
  double height;
  int iplace;
  num_nonzero_elem=0;
  realvolume=0.0;
  for(i=0; i<num_buck; i++)
    if(*(buck+i)){

      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){

	Element* Curr_El=(Element*)(currentPtr->value);
	if((Curr_El->get_adapted_flag()>0)&&
	   (myid==Curr_El->get_myprocess())){

	  height=*(Curr_El->get_state_vars()); 
	  if(height>GEOFLOW_TINY){
	    
	    /* place this height in right (sorted in descending order) 
	       spot in height/volume arrray */
	    for(iplace=num_nonzero_elem;iplace>0;iplace--)
	      if(height>myprochvol[iplace-1][0]){
		myprochvol[iplace][0]=myprochvol[iplace-1][0];
		myprochvol[iplace][1]=myprochvol[iplace-1][1];}
	      else break;
	      
	    myprochvol[iplace][0]=height;
	    myprochvol[iplace][1]=height*
	      *(Curr_El->get_dx())**(Curr_El->get_dx()+1);
	    realvolume+=myprochvol[iplace][1];
	    num_nonzero_elem++;
	  }
	}
	currentPtr=currentPtr->next;      	    
      }
    }

  /* now processor zero will use those list(s) to determine the cut
     off height */

  if(numproc==1){ //it's for a single processor so it's very simple

    cutoffvolume=STAT_VOL_FRAC*realvolume;
    statvolume  =0.0;
    i=0;
    while((statvolume<cutoffvolume)&&(i<num_nonzero_elem)){
      cutoffheight=myprochvol[i][0];
      statvolume +=myprochvol[i++][1];
    }

    CDeAllocD2(myprochvol);
  }
  else{// to get a clear idea of what I want to do take a look
       //  at the single processor version (the else to this if)
       // and ask your self what you would need to do if instead
       // of one array of sorted heights-volume pairs, you had
       // numproc of them... that is at every step you look at 
       // the largest not yet counted height from each array and
       // add the largest ones volume to the total 


    //fprintf(fp,"calc_stats() 4\n");  fflush(fp);

    //add up the real volume across all processors
    temp=realvolume;
    MPI_Reduce(&temp,&realvolume,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);




    if(myid>0){ //send my processor's list of height and volume to processor 0
      //fprintf(fp,"calc_stats() 5\n");  fflush(fp);

      MPI_Isend(myprochvol[0],num_nonzero_elem*2,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&request);
    }
    else{ //this is processor zero, myid==0
      //fprintf(fp,"calc_stats() 6\n");  fflush(fp);


      //allocate one big 3D array to hold all the lists from across 
      //all processors 
      MPI_Waitall(numproc-1,proc0request+1,proc0status+1);
      int max_num_nonzero_elem=num_nonzero_elem;
      for(iproc=1;iproc<numproc;iproc++)
	if(max_num_nonzero_elem<all_num_nonzero_elem[iproc])
	  max_num_nonzero_elem=all_num_nonzero_elem[iproc];

      double ***hvol=CAllocD3(numproc,max_num_nonzero_elem,2);

      //fprintf(fp,"calc_stats() 7\n");  fflush(fp);

      //bring the individual processors' sorted lists of height and 
      //volume (hvol) over to processor 0 
      for(iproc=1;iproc<numproc;iproc++)
	MPI_Irecv(hvol[iproc][0],all_num_nonzero_elem[iproc]*2,MPI_DOUBLE,iproc,2,MPI_COMM_WORLD,&proc0request[iproc]);

      for(i=0;i<all_num_nonzero_elem[0];i++){
	hvol[0][i][0]=myprochvol[i][0];
	hvol[0][i][1]=myprochvol[i][1];}

      //fprintf(fp,"calc_stats() 8\n");  fflush(fp);

      MPI_Waitall(numproc-1,proc0request+1,proc0status+1);

      free(proc0request);
      free(proc0status);

      //fprintf(fp,"calc_stats() 9\n");  fflush(fp);

      /* find the cut off height that results in 
	 statvolume= STAT_VOL_FRAC*realvolume 
	 statvolume= volume represented in stats
	 realvolume= total volume */

      cutoffvolume=STAT_VOL_FRAC*realvolume;
      statvolume  =0.0;

      int *icurrent=CAllocI1(numproc);

      //fprintf(fp,"calc_stats() 10\n");  fflush(fp);

      //set the current "element" on each "processor" to the one with
      //the maximum height on that processor (the first one)
      for(iproc=0;iproc<numproc;iproc++) icurrent[iproc]=0;
    
      int imaxh=0;
      //check to see if the cut off height is small enough that we have
      //equaled or exceeded the cut off volume
      while((statvolume<cutoffvolume)&&
	    (icurrent[imaxh]<max_num_nonzero_elem)){

	//find the largest height among the "current" elements and make 
	//it the cut off height 
	for(iproc=0;iproc<numproc;iproc++)
	  if((icurrent[iproc]<all_num_nonzero_elem[iproc])&&
	     (hvol[iproc][icurrent[iproc]][0]>
	      hvol[imaxh][icurrent[imaxh]][0])) imaxh=iproc;	
	cutoffheight=hvol[imaxh][icurrent[imaxh]][0];

	//add that element's volume to the total
	statvolume +=hvol[imaxh][icurrent[imaxh]++][1];
      }

      //fprintf(fp,"calc_stats() 11\n");  fflush(fp);

      CDeAllocD3(hvol);
      CDeAllocI1(all_num_nonzero_elem);
      CDeAllocI1(icurrent);

    }
    //fprintf(fp,"calc_stats() 12\n");  fflush(fp);

    CDeAllocD2(myprochvol);

    //fprintf(fp,"calc_stats() 13\n");  fflush(fp);

    MPI_Bcast(&cutoffheight,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    //fprintf(fp,"calc_stats() 14\n");  fflush(fp);

  }
#endif

  /**************************************************/
  /****** TADA!!!!!!!!! we have finally found  ******/
  /****** this iteration's cut off height (and ******/
  /****** the global volumes too) now we can   ******/
  /****** calculate the rest of the stats in a ******/
  /****** straight forward manner              ******/
  /**************************************************/

  double VxVy[2];
  for(i=0; i<num_buck; i++)
    if(*(buck+i)){

      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){

	Element* Curr_El=(Element*)(currentPtr->value);
	assert(Curr_El!=NULL);
	if((Curr_El->get_adapted_flag()>0)&&
	   (myid==Curr_El->get_myprocess())){
	      
	  double* state_vars = Curr_El->get_state_vars();
	  
	  //calculate volume passing through "discharge planes"
	  unsigned *nodes= Curr_El->getNode();
	  double nodescoord[9][2], *coord;
	  Node* node;
	  
	  for(int inode=0;inode<8;inode++) {
	    node=(Node*) NodeTable->lookup(nodes+2*inode);
	    coord=node->get_coord();
	    if((timeprops->iter==291)&&(inode==8)) {
	      printf("coord=(%g,%g) node=%u  ",coord[0],coord[1],node);
	      fflush(stdout);
	    }
	    nodescoord[inode][0]=coord[0];
	    nodescoord[inode][1]=coord[1]; 
	    if((timeprops->iter==291)&&(inode>=8)) {
	      printf("inode=%d node=%u",inode,node); fflush(stdout);
	    }
	  }
	  nodescoord[8][0]=*(Curr_El->get_coord());
	  nodescoord[8][1]=*(Curr_El->get_coord()+1);

	  discharge->update(nodescoord,state_vars,d_time);
	    


	  // rule out non physical fast moving thin layers
	  //if(state_vars[0] >= cutoffheight){

	  if(state_vars[0] > min_height){

	    if(state_vars[0] > max_height) max_height = state_vars[0];
	    
	    double* xy = Curr_El->get_coord();

	    if(state_vars[0] >= statprops->hxyminmax){
	      if(xy[0]<xyminmax[0]) //xmin
		xyminmax[0]=xy[0];
	      if(-xy[0]<xyminmax[1]) //negative xmax
		xyminmax[1]=-xy[0];
	      if(xy[1]<xyminmax[2]) //ymin
		xyminmax[2]=xy[1];
	      if(-xy[1]<xyminmax[3]) //negative ymax
		xyminmax[3]=-xy[1];
	    }


	    //to test if pileheight of depth testpointheight
	    //has reached the testpoint
	    testpointdist2=
	      (xy[0]-testpointx)*(xy[0]-testpointx)+
	      (xy[1]-testpointy)*(xy[1]-testpointy);
	    double junktest=0;
	    if(testpointdist2>0) 
	      junktest=testpointdist2;
	    if(testpointmindist2>0) 
	      junktest=testpointmindist2;
	    if(testpointdist2<testpointmindist2) {
	      testpointmindist2=testpointdist2;
	      testpointreach=((state_vars[0]>=testpointheight)?1:0);
	    }

	    dA = *(Curr_El->get_dx())**(Curr_El->get_dx()+1);
	    area+=dA;
	    dVol=state_vars[0]*dA;
	    testvolume+=dVol;
	    xC+=xy[0]*dVol;
	    yC+=xy[1]*dVol;
	    
	    xVar+=xy[0]*xy[0]*dVol;
	    yVar+=xy[1]*xy[1]*dVol;
	    piler2+=(xy[0]*xy[0]+xy[1]*xy[1])*dVol;
	    rC+=sqrt((xy[0]-xCen)*(xy[0]-xCen)+
		     (xy[1]-yCen)*(xy[1]-yCen))*dVol;
	    
	    
	    v_ave+=sqrt(state_vars[1]*state_vars[1]+state_vars[2]*state_vars[2])*dA;
	    Curr_El->eval_velocity(0.0,0.0,VxVy);


	    if((!((v_ave<=0.0)||(0.0<=v_ave)))||
	       (!((state_vars[0]<=0.0)||(0.0<=state_vars[0])))||
	       (!((state_vars[1]<=0.0)||(0.0<=state_vars[1])))||
	       (!((state_vars[2]<=0.0)||(0.0<=state_vars[2])))) {
	      //v_ave is NaN
	      printf("calc_stats(): NaN detected in element={%10u,%10u} at iter=%d\n",
		     *(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),
		     timeprops->iter);
	      printf("prevu={%12.6g,%12.6g,%12.6g}\n",
		     *(Curr_El->get_prev_state_vars()+0),
		     *(Curr_El->get_prev_state_vars()+1),
		     *(Curr_El->get_prev_state_vars()+2));
	      printf("    u={%12.6g,%12.6g,%12.6g}\n",
		     state_vars[0],state_vars[1],state_vars[2]);
	      printf("prev {hVx/h,hVy/h}={%12.6g,%12.6g}\n",
		     *(Curr_El->get_prev_state_vars()+1)/ *(Curr_El->get_prev_state_vars()+0),
		     *(Curr_El->get_prev_state_vars()+2)/ *(Curr_El->get_prev_state_vars()+0));
	      printf("this {hVx/h,hVy/h}={%12.6g,%12.6g}\n",
		     state_vars[1]/state_vars[0],
		     state_vars[2]/state_vars[0]);
	      printf("     { Vx  , Vy  }={%12.6g,%12.6g}\n",VxVy[0],VxVy[1]);
	      ElemBackgroundCheck2(El_Table,NodeTable,Curr_El,stdout);
	      assert(0);
	    }

	    temp=sqrt(VxVy[0]*VxVy[0]+VxVy[1]*VxVy[1]);
	    if(temp>v_max) v_max=temp;
	    vx_ave+=state_vars[1]*dA;
	    vy_ave+=state_vars[2]*dA;
	    
	    //these are garbage, Bin Yu wanted them when he was trying to come up 
	    //with a global stopping criteria (to stop the calculation, not the pile)
	    //volume averaged slope in the direction of velocity
	    //a negative number means the flow is headed uphill
	    double resolution=0, xslope=0, yslope=0;
	    Get_max_resolution(&resolution);
	    Get_slope(resolution,
		      *((Curr_El->get_coord()))*matprops->LENGTH_SCALE,
		      *((Curr_El->get_coord())+1)*matprops->LENGTH_SCALE,
		      &xslope,&yslope);
	    if(temp>GEOFLOW_TINY){
	      slope_ave+=-(state_vars[1]*xslope+state_vars[2]*yslope)*dA/temp;
	      slopevolume+=dVol;}
	  }
	}
	currentPtr=currentPtr->next;      	    
      }
    }


  MPI_Reduce(xyminmax,statprops->xyminmax,4,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  if(myid==0)
    {
      statprops->xyminmax[0]*=(matprops->LENGTH_SCALE);
      statprops->xyminmax[1]*=-(matprops->LENGTH_SCALE);
      statprops->xyminmax[2]*=(matprops->LENGTH_SCALE);
      statprops->xyminmax[3]*=-(matprops->LENGTH_SCALE);
    }
  
  int inttempout;
  double tempin[13], tempout[13], temp2in[2], temp2out[2]; 

  //find the minimum distance (squared) to the test point
  MPI_Allreduce(&testpointmindist2,tempout,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  //if this processor isn't the closest to the test point it doesn't count as it's flow reaching the point
  if(tempout[0]<testpointmindist2) testpointreach=0;
  
  //did the closest point to the test point get reached by the flow?
  MPI_Reduce(&testpointreach,&inttempout,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
  testpointreach=inttempout;

  tempin[0]=xC;
  tempin[1]=yC;
  tempin[2]=rC;
  tempin[3]=area;
  tempin[4]=v_ave;
  tempin[5]=vx_ave;
  tempin[6]=vy_ave;
  tempin[7]=slope_ave;
  tempin[8]=piler2;
  tempin[9]=slopevolume;
  tempin[10]=testvolume;
  tempin[11]=xVar;
  tempin[12]=yVar;

  i = MPI_Reduce(tempin, tempout, 13, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  temp2in[0]=max_height;
  temp2in[1]=v_max;
  i = MPI_Reduce(temp2in, temp2out, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      
  if(myid == 0) 
  {
    if(testpointreach&&(statprops->timereached<0.0))
      statprops->timereached=timeprops->timesec();

    double VELOCITY_SCALE=sqrt(matprops->LENGTH_SCALE* 
			       matprops->GRAVITY_SCALE);
    //dimensionalize
    statprops->xcen=tempout[0]*(matprops->LENGTH_SCALE)/tempout[10];
    statprops->ycen=tempout[1]*(matprops->LENGTH_SCALE)/tempout[10];
    statprops->xvar=tempout[11]*(matprops->LENGTH_SCALE)*(matprops->LENGTH_SCALE)/tempout[10]-(statprops->xcen)*(statprops->xcen);
    statprops->yvar=tempout[12]*(matprops->LENGTH_SCALE)*(matprops->LENGTH_SCALE)/tempout[10]-(statprops->ycen)*(statprops->ycen);
    statprops->rmean=tempout[2]*(matprops->LENGTH_SCALE)/tempout[10];
    statprops->area=tempout[3]*(matprops->LENGTH_SCALE)* 
      (matprops->LENGTH_SCALE);
    statprops->vmean=tempout[4]*VELOCITY_SCALE/tempout[10];
    statprops->vxmean=tempout[5]*VELOCITY_SCALE/tempout[10];
    statprops->vymean=tempout[6]*VELOCITY_SCALE/tempout[10];

    statprops->slopemean=(tempout[9]>0)?tempout[7]/tempout[9]:0.0;

    statprops->realvolume=realvolume*(matprops->LENGTH_SCALE)* 
      (matprops->LENGTH_SCALE)*(matprops->HEIGHT_SCALE);

    //statvolume is really testvolume which is statvolume if it's not disabled
    statprops->statvolume=tempout[10]*(matprops->LENGTH_SCALE)* 
      (matprops->LENGTH_SCALE)*(matprops->HEIGHT_SCALE);
    
    statprops->cutoffheight=cutoffheight*(matprops->HEIGHT_SCALE);
    testvolume=tempout[10]/statvolume;

    /* the factor of 3^0.5 is a safety factor, this value was chosen because
     * it makes the "radius" of a uniformly distributed line equal to half
     * the line length 
     */
    //3 standard deviations out ~ 99.5% of the material
    statprops->piler=3.0*sqrt(statprops->xvar+statprops->yvar); 
    statprops->hmax=temp2out[0]*(matprops->HEIGHT_SCALE);
    statprops->vmax=temp2out[1]*VELOCITY_SCALE;

    /* v_star is the nondimensional global average velocity by v_slump
       once v_slump HAS BEEN CALIBRATED (not yet done see ../main/datread.C) 
       the calculation will terminate when v_star reaches 1 */
    statprops->vstar=statprops->vmean/matprops->Vslump;

    /******************/
    /* output section */
    /******************/

    /* output Center Of Mass and x and y components of mean velocity to
       assist the dynamic gis update daemon */
    FILE* fp2=fopen("com.up","w");
    fprintf(fp2,"%d %g %g %g %g %g %g\n",
	    timeprops->iter,timeprops->timesec(),
	    statprops->xcen,statprops->ycen,
	    statprops->vxmean,statprops->vymean,statprops->piler);
    fclose(fp2);


    /* standard to screen output */
    d_time*=timeprops->TIME_SCALE;
    //chunk time
    int hours, minutes; double seconds;
    timeprops->chunktime(&hours, &minutes, &seconds);

    printf("At the end of time step %d the time is %d:%02d:%g (hrs:min:sec),\n \ 
           time step length is %g [sec], volume is %g [m^3],\n \
           max height is %g [m], max velocity is %g [m/s],\n \
           ave velocity is %g [m/s], v* = %g\n\n",
	   timeprops->iter, hours, minutes, seconds, d_time,
	   statprops->statvolume, statprops->hmax, statprops->vmax,
	   statprops->vmean, statprops->vstar); //, statprops->cutoffheight,
  }

  return;
}



void out_final_stats(TimeProps* timeprops, StatProps* statprops){
  //round the run time to the nearest second and chunk it  
  int walltime =(int) (time(NULL)-timeprops->starttime+0.5);
  int wallhours   =  walltime/3600;
  int wallminutes = (walltime%3600)/60;
  int wallseconds =  walltime%60;
  int cputime = (int) (clock()/CLOCKS_PER_SEC+0.5);
  int cpuhours   =  cputime/3600;
  int cpuminutes = (cputime%3600)/60;
  int cpuseconds =  cputime%60;

  char filename[256];
  sprintf(filename,"finalstats.%06d",statprops->runid);
  //printf("runid=%d\n",statprops->runid);
  FILE* fp=fopen(filename,"w");
    //print    runid  xC yC   rC  area maxh    wall time       cpu time
  fprintf(fp,"%6d   %E %E   %E   %E   %E   %E   %d:%02d:%02d   %d:%02d:%02d\n",
	  statprops->runid,statprops->xcen,statprops->ycen,
	  statprops->rmean,statprops->area,statprops->hmax,
	  statprops->timereached,
	  wallhours,wallminutes,wallseconds,
	  cpuhours,cpuminutes,cpuseconds);
  fflush(fp);
  fclose(fp);
  
  return;
}


void InsanityCheck(HashTable* El_Table, int nump, int myid, 
		  TimeProps *timeprops_ptr){

  int num_buck=El_Table->get_no_of_buckets();
  HashEntryPtr* buck = El_Table->getbucketptr();
  for(int i=0; i<num_buck; i++)
    if(*(buck+i)){

      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){

	Element* Curr_El=(Element*)(currentPtr->value);
	currentPtr=currentPtr->next;      	    

	if((Curr_El->get_refined_flag()==0)&&
	   !((Curr_El->get_adapted_flag()==NOTRECADAPTED)||
	     (Curr_El->get_adapted_flag()==NEWFATHER)||
	     (Curr_El->get_adapted_flag()==NEWSON)||
	     (Curr_El->get_adapted_flag()==BUFFER)
	     )	   
	   ){
	  printf("FUBAR 1 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}
	
	
	if((Curr_El->get_refined_flag()!=0)&&
	   ((Curr_El->get_adapted_flag()==NOTRECADAPTED)||
	    (Curr_El->get_adapted_flag()==NEWFATHER)||
	    (Curr_El->get_adapted_flag()==NEWSON)||
	    (Curr_El->get_adapted_flag()==BUFFER)
	    )
	   ){
	  printf("FUBAR 2 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}
	
	if((Curr_El->get_refined_flag()==GHOST)&&
	   !((Curr_El->get_adapted_flag()<=-NOTRECADAPTED)&&
	     (Curr_El->get_adapted_flag()>=-BUFFER)
	     )
	   ){
	  printf("FUBAR 3 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}

	if((Curr_El->get_refined_flag()!=GHOST)&&
	   ((Curr_El->get_adapted_flag()<=-NOTRECADAPTED)&&
	    (Curr_El->get_adapted_flag()>=-BUFFER)
	     )
	   ){
	  printf("FUBAR 4 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}

	if((Curr_El->get_refined_flag()>0)&&
	   !((Curr_El->get_adapted_flag()==TOBEDELETED)||
	     (Curr_El->get_adapted_flag()==OLDFATHER)||
	     (Curr_El->get_adapted_flag()==OLDSON)
	     )
	   ){
	  printf("FUBAR 5 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}

	if(!(Curr_El->get_refined_flag()>0)&&
	   ((Curr_El->get_adapted_flag()==TOBEDELETED)||
	    (Curr_El->get_adapted_flag()==OLDFATHER)||
	    (Curr_El->get_adapted_flag()==OLDSON)
	    )
	   ){
	  printf("FUBAR 6 in InsanityCheck()\nnump=%d myid=%d iter=%d time=%g[sec]\nElement={%u,%u} myprocess=%d refined=%d adapted=%d\n",nump,myid,timeprops_ptr->iter,timeprops_ptr->timesec(),*(Curr_El->pass_key()+0),*(Curr_El->pass_key()+1),Curr_El->get_myprocess(),Curr_El->get_refined_flag(),Curr_El->get_adapted_flag());
	  assert(0);
	}
      }
    }



  return;
}
