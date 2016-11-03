#include "../header/hpfem.h"
#include "./curve.C"
#include "./contfun.C"
/*=========================================================
	Web Visualization modeule
NEW: Feb 08 2004:  changed to write output files in correcsponding
	
OUT:	in topo		web_gridxxxx.out file in ASCII
	in flow		web_outfxxxxyyy.out in ASCII

Also..python script runit.py in bin directory modifies to 
create separate directory structure for web visualization
	webviz----topo : contains topo files
	       |--data : contains flow files
               !--textures: contains textures for terrain

=========================================================== */

void web_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int myid, double time,int numprocs,	MatProps* matprops, 
		TimeProps* timeprops)
{
  int i,j, k;
  int element_counter = 0;
  Element* EmTemp;
  HashEntry* entryp;
  Element* EmTemp_nxt;
  HashEntry* nextpt;
  double xcurrent,ycurrent,zcurrent,pcurrent,totalht;
  double xposix,yposix,zposix,pposix,totalhtx;
  double xposiy,yposiy,zposiy,pposiy,totalhty;
  double max_pile,cutoff_pile;  

  double xnegix,ynegix,znegix,pnegix,totalhtxm;
  double xnegiy,ynegiy,znegiy,pnegiy,totalhtym;
  double xxy,yxy,zxy,pxy,totalhtxy;
  double phmin=1000.0,phmax=0.0;
  char flowpath[70];
  char topopath[70];
  char flowdir[]= "./webviz/data/";
  char topodir[]= "./webviz/topo/";
  char rundata[] = "web_outfxxxxxxxx.out";
  char infofile1[] = "./webviz/web_runinfo.info";
  char infofile2[] = "./webviz/webdata.info";
  char gridfile[] = "web_gridxxxx.out";
  int topolevels=6;
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();
  FILE*    fp1;
  FILE*    fp2;
  FILE*    fp3;

  printf("\nWebViz: Processing iteration: %d\n",timeprops->iter);

   // temp addition to write in time steps first line

   int hr,mn;
   double secd;
   timeprops->chunktime( &hr, &mn, &secd);
   printf("ch Time: %2d:%2d:%g", hr,mn,secd);
   printf("\nTime : %g\n",timeprops->timesec());
   sprintf(rundata,"web_outf%05d%03d.out",timeprops->ioutput,myid);
  

   // temp addition to write in time steps first line...COMPLETE

   strcpy(flowpath,flowdir);
   strcat(flowpath,rundata);
   if(timeprops->ifstart()){
      fp1 = fopen (flowpath,"w");
      fp2 = fopen (infofile1,"w");
      fp3 = fopen (infofile2,"w");
   }
   else{ 
      fp1 = fopen (flowpath,"a+");
      fp2 = fopen (infofile1,"a+");
   }

   if(fp1==NULL||fp2==NULL){
	printf("\n ERROR: could not open files : %s",flowpath);
	printf("\n    OR: could not open files : %s",infofile1);
	exit(1);
   }
   fprintf(fp1,"Chunk Time: %2d:%2d:%g\tIncrement:%lf\n", hr,mn,secd,timeprops->iter);
   double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE));
   double  xmin,xmax,ymin,ymax;

   for(i=0; i<e_buckets; i++){
      	entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      	while(entryp){	
	  	EmTemp = (Element*)entryp->value;
	  	assert(EmTemp);
	  	if(EmTemp->get_adapted_flag()>0){
	      		int* neigh_ptr = EmTemp->get_neigh_proc();
	      		if(neigh_ptr[1] != INIT && neigh_ptr[2] != INIT)
			element_counter++;
	    	}
	  	entryp = entryp->next;
	}
   }

   j=Get_window(&xmin,&xmax,&ymin,&ymax);
	
   int topolvl=1;	// topo resolution level: 1 is coarsest
//=============================================================
//   Creation of multi resolution terrain with keypoint regions
//=============================================================
  double resolution,spanx,spany,xincr,yincr;
  int top_lv_nox,top_lv_noy;
  double top_lv_resx,top_lv_resy;
  int max_kregions,x_kregions,y_kregions;
  double x_kregion_wid,y_kregion_wid;
  if(timeprops->ioutput==0 && myid==0){
      fprintf(fp3,"TexWd:\t0\n");
      fprintf(fp3,"TexHt:\t0\n");
	
      for(topolvl=1;topolvl<topolevels;topolvl++)
	{
	  double elevatlocation;

	  max_kregions=pow(pow(2.0,topolvl-1),2);
	  x_kregions=pow(2.0,topolvl-1);
	  y_kregions=pow(2.0,topolvl-1);

  	  FILE*    fp[max_kregions];
	  int region =1;

	  //creatin each key region 128 cell wide

	  top_lv_nox=128*pow(2.0,topolvl-1);
	  top_lv_noy=128*pow(2.0,topolvl-1);

	  spanx = xmax-xmin;
	  spany = ymax-ymin;

	  top_lv_resx=spanx/(double)top_lv_nox;
	  top_lv_resy=spany/(double)top_lv_noy;

	  xincr=spanx/top_lv_nox;
	  yincr=spany/top_lv_noy;

	  x_kregion_wid=spanx/x_kregions;
	  y_kregion_wid=spany/y_kregions;
	
	  double xlocations[top_lv_nox*top_lv_noy]; 
	  double ylocations[top_lv_nox*top_lv_noy];
	  double elevation_at_locations[(top_lv_nox)*(top_lv_noy)] ;

	  j = Get_max_resolution(&resolution);

	  if(resolution>top_lv_resx){
	      printf("\nWARNING:\nrequested resolution =\t%d",top_lv_resx);
	      printf("\nfinest possible resolution:\t%lf",resolution);
	      printf("\nCannot go further than this\n");
	      break;
	  }
	  else{
		fprintf(fp3,"Resolution %d:\t%lf\t%lf\n",topolvl,top_lv_resx,top_lv_resy);
		for(int i=0;i<x_kregions;i++){
			for(int j=0;j<y_kregions;j++){
			fprintf(fp3,"xmin:\t%lf\tymin:\t%lf\t",xmin+i*x_kregion_wid,ymin+j*y_kregion_wid);
			fprintf(fp3,"xmax:\t%lf\tymax:\t%lf\n",xmin+(i+1)*x_kregion_wid,ymin+(j+1)*y_kregion_wid);
			

			}
		}
	  }

	  if(j != 0){
	      printf("error in Get_elevation_grid\n");
	  }    

	  if(topolvl==1){
	      gridfile[8] = (topolvl % 10) + 48;
	      gridfile[9]  = (region % 1000)/100 + 48;
	      gridfile[10] = (region % 100)/10 + 48;
	      gridfile[11] = (region % 10) + 48;
	      region=region+1;
	
	      strcpy(topopath,topodir);
	      strcat(topopath,gridfile);

	      fp[region-1]=fopen(topopath,"w+");
	      if(fp1==NULL||fp2==NULL){
		printf("\n ERROR: could not open files : %s",topopath);
	      }

	      fprintf(fp2,"%d\n",numprocs);

	      for(i=0;i<top_lv_noy;i++){
	          for(j=0;j<top_lv_nox;j++){
		      Get_elevation(top_lv_resx,(xmin+j*xincr),(ymax-i*yincr),&elevatlocation);
		      fprintf(fp[region-1],"%lf\t%lf\t%lf\n",xmin+j*xincr,ymax-i*yincr,elevatlocation);
		  }
	      }	
	      fclose(fp[region-1]);
	  } // top level topo file written in one piece
	  else{

	      int xsteps=top_lv_nox/x_kregions;
	      int ysteps=top_lv_noy/y_kregions;

	      region =0;
	      for(int xx=0;xx<x_kregions;xx++){
		  for(int yy=0;yy<y_kregions;yy++){		
		      region=region+1;
		      gridfile[8] = (topolvl % 10) + 48;
		      gridfile[9]  = (region % 1000)/100 + 48;
		      gridfile[10] = (region % 100)/10 + 48;
		      gridfile[11] = (region % 10) + 48;

		      strcpy(topopath,topodir);
	              strcat(topopath,gridfile);
		      printf("\n opening : %s ",topopath);
		      fp[region-1]=fopen(topopath,"w+");
		      for(i=0;i<ysteps;i++){
			  for(j=0;j<xsteps;j++){
			      Get_elevation(top_lv_resx,(xmin+(xx*x_kregion_wid)+j*xincr),
			      (ymax-(yy*y_kregion_wid)-i*yincr),&elevatlocation);
			      fprintf(fp[region-1],"%lf\t%lf\t%lf\n",
			      (xmin+(xx*x_kregion_wid)+j*xincr),	
			      (ymax-(yy*y_kregion_wid)-i*yincr),elevatlocation);
			  }
		      }
	      	      fclose(fp[region-1]);
		  }
	      }
	  }//else 

	}// for loop topo levels
	fprintf(fp3,"%d",topolvl-1);
    }// if timestep=0

//=============================================================
//   Completing multi resolution terrain generation 
//=============================================================
  int hours, minutes;
  double seconds;
  timeprops->chunktime(&hours,&minutes,&seconds);

  int write_flag=0;

  for(i=0; i<e_buckets; i++)
    {

      entryp = *(HT_Elem_Ptr->getbucketptr() + i);

      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    {
	      int* neigh_ptr = EmTemp->get_neigh_proc();
	      // get the elevation at the nodes which is the elevation at the finest scale and doesn't change
	      Node* NdTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());
	      double elevation, elevationx, elevationy, elevationxm, elevationym, elevationxy;
	      elevation = NdTemp->get_elevation()*(matprops->LENGTH_SCALE);
			    
	      xcurrent= (*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
	      ycurrent= (*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
	      zcurrent= elevation;
	      pcurrent=  (*(EmTemp->get_state_vars()))*(matprops)->HEIGHT_SCALE;
	      totalht= elevation+(*(EmTemp->get_state_vars()))*(matprops)->HEIGHT_SCALE;


	      /* writing for current element location information complete
		 now starting getting information about neighboring  elements
		 till if(time_step==0) line 
	      */

	      int p_x_side = EmTemp->get_positive_x_side();
	      int xp, xm, yp, ym; //x plus, x minus, y plus, y minus
	      xp = p_x_side;
	      switch(p_x_side) {
	      case 0:
		xm = 2;
		yp = 1; 
		ym = 3;
		break;
	      case 1:
		xm = 3;
		yp = 2;
		ym = 0;
		break;
	      case 2:
		xm = 0;
		yp = 3;
		ym = 1;
		break;
	      case 3:
		xm = 1;
		yp = 0;
		ym = 2;
		break;
	      }

	      if(neigh_ptr[xp] != INIT && neigh_ptr[yp] != INIT) {
		Element *EmTempx, *EmTempy, *EmTempxm, *EmTempym;
		EmTempx = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+xp*KEYLENGTH);
		EmTempxm = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+xm*KEYLENGTH);
		EmTempy = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+yp*KEYLENGTH);
		EmTempym = (Element*) HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+ym*KEYLENGTH);

		// getting elevation at x location

		NdTemp = (Node*) HT_Node_Ptr->lookup(EmTempx->pass_key());

		if(NdTemp != NULL) // if EmTempx is a ghost cell the bubble node may not be here
		  elevationx = NdTemp->get_elevation()*(matprops->LENGTH_SCALE);
		else 
		  {  // the bubble node is not here so we get the topo data from GIS
		    j = Get_max_resolution(&resolution);
		    if(j != 0) {
		      printf("error in Get_max_resolution\n");
		      exit(1);
		    }    
		    double xcoord = *(EmTempx->get_coord())*(matprops->LENGTH_SCALE);
		    double ycoord = *(EmTempx->get_coord()+1)*(matprops->LENGTH_SCALE);
		    j = Get_elevation(resolution, xcoord, ycoord, &elevationx);
		    if(j != 0) {
		      printf("error in Get_elevation\n");
		      exit(1);
		    }   
		    else{
		      xposix=xcoord;
		      yposix=ycoord;
		      zposix=elevationx;
		    }
		  }// else endss

		// Getting elevation at Y location


		NdTemp = (Node*) HT_Node_Ptr->lookup(EmTempy->pass_key());
		if(NdTemp != NULL) // if EmTempy is a ghost cell the bubble node may not be here
		  elevationy = NdTemp->get_elevation()*(matprops->LENGTH_SCALE);
		else {  // the bubble node is not here so we get the topo data from GIS
		  double resolution;
		  j = Get_max_resolution(&resolution);
		  if(j != 0) {
		    printf("error in Get_max_resolution\n");
		    exit(1);
		  }    
		  double xcoord = *(EmTempy->get_coord())*(matprops->LENGTH_SCALE);
		  double ycoord = *(EmTempy->get_coord()+1)*(matprops->LENGTH_SCALE);
		  j = Get_elevation(resolution, xcoord, ycoord, &elevationy);
		  if(j != 0) {
		    printf("error in Get_elevation\n");
		    exit(1);
		  }   
		
		  else{
		    xposiy=xcoord;
		    yposiy=ycoord;
		    zposiy=elevationy;
		  }

		}   // else ends for y
		

		if(timeprops->ifstart()){
		    fprintf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\n",xcurrent,ycurrent,
						zcurrent,pcurrent,totalht);
		    if(pcurrent<phmin && pcurrent!=0.00)phmin=pcurrent;
		    if(pcurrent>phmax)phmax=pcurrent;

		}
		else{
		    ifstream inD2 (infofile1);
		    if(inD2.fail()){
		      printf("cant find run information file for web %s \n",infofile1);
		      exit(0);
		    }

		    max_pile = 0;  // the number of separate files to use
			
		    int line_number=0;
		    double readpilemax;
		    double readpilemin;

		    while(!inD2.eof()){
				inD2>>readpilemax;	
				inD2>>readpilemin;	
		    }
		
		    inD2.close();		
		    max_pile=readpilemax;
		    cutoff_pile=0.08*max_pile;

		    // setting skip flag : if (xnext- xcurrent) > 2% spanx , skip	 
		    bool skip_flag=0;
		    if((0.002*spany)<(xposiy-ycurrent)) skip_flag=1;
		    else skip_flag=0;
		    
		    if(!skip_flag){
			if(pcurrent>=cutoff_pile){
			    if(write_flag){
				if(pcurrent<phmin && pcurrent!=0.0)phmin=pcurrent;
				if(pcurrent>phmax)phmax=pcurrent;

				fprintf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\n",xcurrent,
				       ycurrent,zcurrent,pcurrent,totalht);
			    }
			    if(write_flag) write_flag=0;
			    else write_flag=1;
			}
		     } ///if skip flag
		    else printf("\n Watch the Skip Flag in webviz!! \n");
		}
	  
	      } // if loop  neigh ptr[xp]
	  
	    }  //if loop get adapted

	  entryp = entryp->next; 
	  
	} //while loop entryp

    } //for loop buckets

  fprintf(fp2,"%lf\t%lf\n",phmax,phmin);
  fprintf ( fp1, "\n" );
  fclose ( fp1 );
  fclose ( fp2 );
  if(timeprops->ifstart())fclose ( fp3 );

  return;
} 


/*****************************************************************************
This Program reads the flow data created by web_output module of TITAN2D and 
arranges all the files in ordered manner according to x,y locations. 
ALSO does filtering of too densely located points
Modiefied to read ultiprocessor files last three digits are for processor number

This is the latest and most versatile version, finalized before my defence.

IN:	flow files from TITAN2D named as 'web_outfxxxxxYYY.out'
OUT:	rearranged flow files not renamed 'web_outfxxxxx.pp1'
          ( saying post processing 1 done )

NOTE:	*Original flow files are not destroyed. 
*The files created by this code are read by flow projection code which
         corrects elevation values.
	
USAGE:   $executable-name  run-directory-path max-timesteps  time-step-incr
******************************************************************************/
void web_simplify(TimeProps* timeprops){
  
  int selection=9;
  char FLOW_FILE_PATH[100] ;		// example : "./cnt_colimaR5k/"
  char FLOW_FILE[] = "/webviz/data/web_outfxxxxx000.out";
  char INFO_FILE[] = "/webviz/web_runinfo.info";
  char OUT_FILE[]=   "/webviz/data/web_outfxxxxx.pp1";
  char OUT_FILE_PATH[100];
  char GRID_FILE[100];
  char FLOW_INFO[60];	// complete referance to flow info file
  int time_start ;
  int time_step_incr;
  int time_max;
  int start_index = 21;	
  double phmax;
  double phmin;
  int max_processors;
  
 
  strcpy(FLOW_FILE_PATH,"./");
  time_max= timeprops->ioutput;
  //time_step_incr= atoi(argv[3]);
  time_start = 1;

  strcpy(FLOW_INFO,FLOW_FILE_PATH);
  strcat(FLOW_INFO,INFO_FILE);	
	
  ifstream maxmin(FLOW_INFO);
  if(maxmin==NULL){
    cout<<endl<<"Unable to find flow info file...giving up!"<<endl;
    exit(1);
  }
  for(int time_step=1;time_step<=time_max;time_step++){
    char WORKING_FLOW_FILE[70];
    char lineOne[100];
    double ContLvl[5];
    if(time_step==time_start){ 
        maxmin>>max_processors;
        maxmin>>phmax;
        maxmin>>phmin;
    }
    maxmin>>phmax;
    maxmin>>phmin;

    double incr=(phmax-phmin)/6;
    for(int i=0;i<5;i++) ContLvl[i]=phmin+double(i+1)*incr;
//Creating curve to store individual contours
    curve table;
    double tempX,tempY,tempZ,ph,totht;
    double filetemp[5];

    for(int original_proc=0;original_proc<max_processors;original_proc++){

    sprintf(FLOW_FILE,"/webviz/data/web_outf%05d%03d.out",time_step,original_proc);
	
    	strcpy(WORKING_FLOW_FILE, FLOW_FILE_PATH);
    	strcat(WORKING_FLOW_FILE, FLOW_FILE);

    //=======================================================
    // File Reading section 
    //=======================================================
    	cout<<endl<<"Reading file "<<WORKING_FLOW_FILE;
        ifstream infile;
   	infile.open(WORKING_FLOW_FILE,ios::in);
    	if(infile==NULL){
      		cout<<endl<<"Error opening file..bye";
 	        exit(1);
    	}
	infile.getline(lineOne,100);
	int tempcount=0;
    	while(infile){
	      for(int i=0;i<5;i++) infile>>filetemp[i];
      	      tempX=filetemp[0];	
	      tempY=filetemp[1];	
	      tempZ=filetemp[2];	
	      ph=filetemp[3];
	      totht=filetemp[4];
      	     table.arrange(tempX,tempY,tempZ,ph,totht);
	     tempcount++;
    	}
    }

    //========================================================================
    //  Point selection and filtering...none (deleted)
    //========================================================================
    double Xminmax[2],Zminmax[2],Yminmax[2],Phminmax[2];
    table.getextents( Xminmax,Yminmax,Zminmax,Phminmax);

    //========================================================================
    //  Creation of new grid  ...   assuming 100*100
    //========================================================================

    node EvenGrid[GW][GW];

    double XfRange=Xminmax[1]-Xminmax[0];				  
    double YfRange=Yminmax[1]-Yminmax[0];				  
    double SrchRad=10;

    double XfIncr= XfRange/double(GW);
    double YfIncr= YfRange/double(GW);

    double  PhArray[GW][GW];

    for(int XC=0;XC<GW;XC++){
      for(int YC=0;YC<GW;YC++){
	EvenGrid[XC][YC].Xdata=XfIncr*double(XC)+Xminmax[0];
	EvenGrid[XC][YC].Ydata=YfIncr*double(YC)+Yminmax[0];
      }
    }

    strcpy(GRID_FILE,"./");
    strcat(GRID_FILE,"/GridFile.txt");	
    FILE* GridFile;
    GridFile=fopen(GRID_FILE,"w");

    //=======================================================================
    //	SetPile function calculates average of pile heights of points which lie
    //  withinn a proximity circle of radius 10 (hard coded) around the grid 
    //  location
    //=======================================================================
    //cout<<endl<<"Interpolating points, creating grid";
    for(int XC=0;XC<GW;XC++){
      for(int YC=0;YC<GW;YC++){
	SetPile(table,EvenGrid,XC,YC);
	PhArray[XC][YC]=EvenGrid[XC][YC].ph;
	fprintf(GridFile,"%lf\t%lf\t%lf\n",EvenGrid[XC][YC].Xdata,EvenGrid[XC][YC].Ydata,PhArray[XC][YC]);
      }
    }
    //cout<<"..done";

    //========================================================================
    // Any contouring function should come right after this because now we have
    // a 2D array names PhArray[100][100] ready with pile heights listed at 
    // each node location
    //=========================================================================
		
    //=========================================================================
    //       inserting array of pile heights(PhArray,contourheight)
    //=========================================================================
    sprintf(OUT_FILE,"/webviz/data/web_outf%05d.pp1",time_step);
    strcpy(OUT_FILE_PATH,"./");
    strcat(OUT_FILE_PATH,OUT_FILE);
    //cout<<endl<<"Opengin  file for writing : "<<OUT_FILE_PATH;
    FILE* contfile;
    contfile=fopen(OUT_FILE_PATH,"wb");	
    if(contfile==NULL){
      cout<<endl<<"Error opening file " <<OUT_FILE_PATH;
      exit(1);
    }
    //else cout << " ..successful";
	
    Contour contour_inst= Contour(100);
    contour_inst.LoadData(PhArray,EvenGrid,ContLvl,contfile);
    //cout<<endl<<"Data Loading complete"<<endl;
    contour_inst.contfun();
    //cout<<endl<<"Contouring complete"<<endl;
    fclose(contfile);
    //cout<<"done"<<endl;
	
  }//timestep loop
  

}

void web_correct(TimeProps* timeprops){

  char FLOW_FILE[] = "web_outfxxxxx.pp1";
  char outrundata[] = "lxxtxxxxx.asc"; //l(level number )t(time step).asc

  int time_step = 0;
  int time_step_incr=1;			//default:being overwritten in main
  int time_max=timeprops->ioutput;			//default:being overwritten in main

  char  FLOW_FILE_PATH[100] ;
  char  TOPO_FILE[100];
  char WORKING_FLOW_FILE[100];
  int start_index=8;
  int NUM_ROWS = 0, NUM_COLS = 0;

  /****Data file variables******/
  int line_count = 0;
  double *X, *Y, *Z;
  double *topoX, *topoY, *topoZ;
  double flowX, flowY, flowZ, newZ;
  double scaleZ=1.05;
  int topo_line_count = 0;
  int flow_line_count = 0;
  /****function declarations******
       void read_process_flow_file (void);
       void find_quad (void);
       void read_topo_file (void);
       void calculate_XY_extents (void);
       ****function declarations******/
  struct locationData{
    double X;
    double Y;
    double Z;
    double Ph;
  };
  locationData currentL;

  // from read topo files function
  FILE *topo_fp;
  double local_topoX, local_topoY, local_topoZ;
  double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
  int i,j,k;
  char string_x[14], string_y[14], string_z[14];

  int ii = 0;
  int rows_flag = -1;

  // read_process_flow _files
  FILE *flow_fp, *new_flow_fp;
  char local_flowX[15], local_flowY[15], local_flowZ[15], local_pile_ht[15],
    local_actual_flowZ[15];
  double junkZ, junk_pht;
  int iii;
	
  int thisLevel=1;
  int thisKeyRegion=1;

  strcpy(FLOW_FILE_PATH,"./");
  strcat(FLOW_FILE_PATH,"/webviz/data/");
  strcpy(TOPO_FILE,"./");
  strcat(TOPO_FILE,"/webviz/topo/web_grid1001.out");

// open the webinfo file and find out how many resolution levels are present
  char infofilename1[]="./webviz/webdata.info";
  FILE* infofile1;
  infofile1=fopen(infofilename1,"r");
  char line[100];
  fgets(line,sizeof(line),infofile1);

  // void read_topo_file (void)
  // {
  topo_fp = fopen (TOPO_FILE, "r");
  while (!feof (topo_fp)){
    fscanf (topo_fp, "%lf%lf%lf", &local_topoX, &local_topoY, &local_topoZ);
    if (local_topoX > xmax)
      xmax = local_topoX;
    if (local_topoX < xmin)
      xmin = local_topoX;

    if (local_topoY > ymax)
      ymax = local_topoY;
    if (local_topoY < ymax)
      ymax = local_topoY;

    if (local_topoZ > zmax)
      zmax = local_topoZ;
    if (local_topoZ < zmax)
      zmax = local_topoZ;

    topo_line_count++;
  }
  fclose (topo_fp);

  topoX = (double *) malloc ((topo_line_count-1) * sizeof (double));
  topoY = (double *) malloc ((topo_line_count-1) * sizeof (double));
  topoZ = (double *) malloc ((topo_line_count-1) * sizeof (double));

  topo_fp = fopen (TOPO_FILE, "r");

  for (i = 0; i < topo_line_count - 1; i++){
    fscanf(topo_fp, "%lf%lf%lf", &topoX[i],&topoY[i],&topoZ[i]);
  }
  fclose (topo_fp);
  // } 

  // void calculate_XY_extents (void)
  // {
  //calculating number of rows & columns

  if (topoX[ii] == topoX[ii + 1])
    rows_flag = 1;
  else if (topoY[ii] == topoY[ii + 1])
    rows_flag = 0;

  for (ii = 0; ii < topo_line_count; ii++){
    if (rows_flag == 1){
      if (topoX[ii] == topoX[ii + 1])
	NUM_ROWS = NUM_ROWS + 1;
      else
	break;
    }
    else if (rows_flag == 0){
      if (topoY[ii] == topoY[ii + 1])
	NUM_COLS = NUM_COLS + 1;
      else
	break;
    }
  }
  if (rows_flag == 1)
    NUM_COLS = topo_line_count / NUM_ROWS;
  if (rows_flag == 0)
    NUM_ROWS = topo_line_count / NUM_COLS;

  NUM_ROWS--;
  NUM_COLS++;
  printf ("\nWebViz:Simplify: NUM_ROWS = %d\tNUM_COLS = %d\n", NUM_ROWS, NUM_COLS);
  //}
  // worked Perfect


  //void read_process_flow_file (void)
  //{
  /* Which flow file to read */  

  time_step++;
  while (time_step <= time_max) {

    sprintf(FLOW_FILE,"web_outf%05d.pp1",time_step);
    int length=strlen(FLOW_FILE);
    char tmp_FLOW_FILE[length];
    strcpy(tmp_FLOW_FILE,FLOW_FILE);

    //sprintf(outrundata,"web_outf%05d.asc",time_step);
    sprintf(outrundata,"l%02dt%05d.asc",thisLevel,time_step);
    flow_line_count = 0;
    strcpy(WORKING_FLOW_FILE,FLOW_FILE_PATH);
    strcat(WORKING_FLOW_FILE,tmp_FLOW_FILE);
    printf ("WebViz:Simplify: Processing %s\n", WORKING_FLOW_FILE);

    flow_fp = fopen (WORKING_FLOW_FILE, "r");
    //========================================================================
    /* Read selected flow file */  
    //========================================================================
    strcpy(WORKING_FLOW_FILE,FLOW_FILE_PATH);
    strcat(WORKING_FLOW_FILE,outrundata);

    new_flow_fp = fopen (WORKING_FLOW_FILE, "w");
  
    while (!feof (flow_fp)) {
      fscanf (flow_fp, "%lf%lf%lf", &flowX, &flowY, &junk_pht);
      flow_line_count++;
	  
      //find_quad();
      // find quad
      int i, j,a ,b,c,d,n, status;
      int i2, j2, a2, b2, c2, d2, n2;
      double A, B, C, maxZ, minZ, yline;
	
      status=0;
      for(j=0; j< NUM_ROWS-1; j++){
	for(iii=0; iii<NUM_COLS-1; iii++){
	  n=NUM_COLS;
	  a=j*n+iii;
	  b=j*n+iii+n;
	  c=j*n+1+iii+n;
	  d=j*n+iii+1;

	  if ((flowX<=topoX[d]) && (flowX<=topoX[c]) && 
	      (flowX>=topoX[a]) &&(flowX>=topoX[b])      && 
	      (flowY<=topoY[a]) && (flowY<=topoY[d])     && 
	      (flowY>=topoY[b]) && (flowY>=topoY[c])) {
	    //=====================================================================
	    // Find Maximum Z of four quad points 
	    //====================================================================

	    /* Determine which triangle of quad flow point it in */
	    yline=((topoY[d]-topoY[b])/(topoX[d]-topoX[b]))
	      *(flowX-topoX[b])+topoY[b];
			
	    /* Flow point in top triangle */
	    if (yline<flowY) {
	      /* Make plane out of points a,b,d */
	      A=(topoY[b]-topoY[a])*(topoZ[d]-topoZ[a])
		-(topoY[d]-topoY[a])*(topoZ[b]-topoZ[a]);
	      B=-((topoX[b]-topoX[a])*(topoZ[d]-topoZ[a])
		  -(topoZ[b]-topoZ[a])*(topoX[d]-topoX[a]));
	      C=(topoX[b]-topoX[a])*(topoY[d]-topoY[a])
		-(topoY[b]-topoY[a])*(topoX[d]-topoX[a]);
	      newZ=(topoZ[a]-(A*(flowX-topoX[a])+
			      B*(flowY-topoY[a]))/C)*scaleZ;
	      status=1;
	      break;
	      /* Flow point in bottom triangle */
	    }// if block yline
	    else {
	      /* Make plane out of points b,c,d */
	      A=(topoY[c]-topoY[b])*(topoZ[d]-topoZ[b])
		-(topoY[d]-topoY[b])*(topoZ[c]-topoZ[b]);
	      B=-((topoX[c]-topoX[b])*(topoZ[d]-topoZ[b])
		  -(topoZ[c]-topoZ[b])*(topoX[d]-topoX[b]));
	      C=(topoX[c]-topoX[b])*(topoY[d]-topoY[b])
		-(topoY[c]-topoY[b])*(topoX[d]-topoX[b]);
	      newZ=(topoZ[b]-(A*(flowX-topoX[b])
			      +B*(flowY-topoY[b]))/C)*scaleZ;
	      status=1;
	      break;
	    }//else block yline
		

	  }//if loop 
	}//for loop columns
	if (status==1) break;
      }// for loop rows
      //}// find_quad function
      // find quad
      //======================================================================
      /* Write new flow point */
      //======================================================================

      currentL.X=flowX;
      currentL.Y=flowY;
      currentL.Z=newZ;
      currentL.Ph=junk_pht;

      fprintf(new_flow_fp, "%lf\t%lf\t%lf\t%lf\n", flowX, flowY,newZ,junk_pht);	
    }
    fclose (flow_fp);
    fclose (new_flow_fp);
   
    time_step++;
  }   
  //}
}
