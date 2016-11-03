
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
 * $Id: restart.C,v 1.4 2004/08/11 15:54:22 kdalbey Exp $ 
 */

#include "../header/hpfem.h"
//#define DEBUGSAVESPLIT
//#define DEBUGSAVEHEADER
#define NUM_CHAR_IN_SAVE_HEADER 16384 //equiv to 4096 integers of space available for save header

int loadrun(int myid, int numprocs, HashTable** NodeTable, 
	    HashTable** ElemTable, MatProps* matprops_ptr,  
	    TimeProps* timeprops_ptr, MapNames *mapnames_ptr, 
	    int *adaptflag_ptr, int *order_flag_ptr, 
	    StatProps* statprops_ptr, DISCHARGE* discharge_ptr,
	    OutLine* outline_ptr)
/* need to handle StatProps and DISCHARGE (should be able to reload the outline from the outline file, not yet programmed in final form)
 */

{
  
  char filename[64];
  sprintf(filename,"restart%04d.this",myid);
  //printf("filename=\"%s\"\n",filename);
  FILE *fp=fopen(filename,"rb");
  
  if(fp==NULL) return(0);


  char tempstring[512];

  fourbytes  temp4;
  eightbytes temp8;

  int numchar;
  fread(&numchar,sizeof(int),1,fp);
  char *header=CAllocC1(numchar);
  fread(header,sizeof(char),numchar,fp);

  //need to change number and use 4 bytes and 8bytes unions
  int Itemp=0, itemp; //first 4 bytes is left for number of bytes in rest of header takes up

  //check restart's file version
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(strncmp(header+Itemp,"one-phase:20060801",temp4.i)) {
    printf("in loadrun() processor %d opened file \"%s\" it version is \"%s\" instead of \"one-phase:20051027\" aborting\n",myid,filename,strncpy(tempstring,header+Itemp,temp4.i));
    fflush(stdout);
    assert(0);
  }
  Itemp+=temp4.i;

  //processor info... for safety check
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(temp4.i!=myid) {
    printf("in loadrun() processor %d opened file \"%s\" was written by processor %d instead of %d, aborting\n",myid,filename,temp4.i,myid);
    fflush(stdout);
    assert(0);
  }

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(temp4.i!=numprocs) {
    printf("in loadrun() processor %d opened file \"%s\" originally run on %d processors, current job is using %d processors, aborting\n",myid,filename,temp4.i,numprocs);
    fflush(stdout);
    assert(0);
  }

  //number of equations...for safety check
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(temp4.i!=EQUATIONS) {
    printf("in loadrun() processor %d opened file \"%s\" which originally had %d equations, this version of titan uses %d equations, aborting\n",myid,filename,temp4.i,EQUATIONS);
    fflush(stdout);
    assert(0);
  }

  //flags
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  *adaptflag_ptr=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  *order_flag_ptr=temp4.i;

  //TIME info
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  timeprops_ptr->iter=temp4.i;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  timeprops_ptr->time=temp8.d; //double

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  timeprops_ptr->TIME_SCALE=temp8.d; //double

  
  timeprops_ptr->ioutput=
    floor(timeprops_ptr->time*timeprops_ptr->TIME_SCALE/
	  timeprops_ptr->timeoutput);
  timeprops_ptr->ndnextoutput=
    ((timeprops_ptr->ioutput+1)*timeprops_ptr->timeoutput)/
    timeprops_ptr->TIME_SCALE;

  timeprops_ptr->isave=
    floor(timeprops_ptr->time*timeprops_ptr->TIME_SCALE/
	  timeprops_ptr->timesave);
  timeprops_ptr->ndnextsave=
    ((timeprops_ptr->isave+1)*timeprops_ptr->timesave)/
    timeprops_ptr->TIME_SCALE;


  //printf("iter=%d timesec=%g:  ioutput=%d nextoutput=%g:   isave=%d nextsave=%g\n",timeprops_ptr->iter,timeprops_ptr->time*timeprops_ptr->TIME_SCALE,timeprops_ptr->ioutput,timeprops_ptr->ndnextoutput*timeprops_ptr->TIME_SCALE,timeprops_ptr->isave,timeprops_ptr->ndnextsave*timeprops_ptr->TIME_SCALE);
  //fflush(stdout);

  //bob

  //GIS file info
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(strncmp(header+Itemp,mapnames_ptr->gis_main,temp4.i)) {
    printf("in loadrun() processor %d opened file \"%s\" old GIS main directory=\"%s\" new GIS main directory=\"%s\" proceeding anyway\n",myid,filename,strncpy(tempstring,header+Itemp,temp4.i),mapnames_ptr->gis_main);
    //fflush(stdout);
  }
  Itemp+=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(strncmp(header+Itemp,mapnames_ptr->gis_sub,temp4.i)) {
    printf("in loadrun() processor %d opened file \"%s\" old GIS sub directory=\"%s\" new GIS sub directory=\"%s\" proceeding anyway\n",myid,filename,strncpy(tempstring,header+Itemp,temp4.i),mapnames_ptr->gis_sub);
    //fflush(stdout);
  }
  Itemp+=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(strncmp(header+Itemp,mapnames_ptr->gis_mapset,temp4.i)) {
    printf("in loadrun() processor %d opened file \"%s\" old GIS mapset =\"%s\" new GIS mapset =\"%s\" proceeding anyway\n",myid,filename,strncpy(tempstring,header+Itemp,temp4.i),mapnames_ptr->gis_mapset);
    //fflush(stdout);
  }
  Itemp+=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(strncmp(header+Itemp,mapnames_ptr->gis_map,temp4.i)) {
    printf("in loadrun() processor %d opened file \"%s\" old GIS map =\"%s\" new GIS map =\"%s\" aborting\n",myid,filename,strncpy(tempstring,header+Itemp,temp4.i),mapnames_ptr->gis_map);
    fflush(stdout);
    assert(0);
  }
  Itemp+=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(temp4.i!=mapnames_ptr->extramaps) {
    printf("in loadrun() processor %d opened file \"%s\" old had extramaps =%d new has extramaps =%d, aborting\n",myid,filename,temp4.i,mapnames_ptr->extramaps);
    fflush(stdout);
    assert(0);
  }
  
  //material and scaling info
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->LENGTH_SCALE=temp8.d; //double

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->HEIGHT_SCALE=temp8.d; //double

  matprops_ptr->epsilon=
    matprops_ptr->HEIGHT_SCALE/matprops_ptr->LENGTH_SCALE;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->GRAVITY_SCALE=temp8.d; //double

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->MAX_NEGLIGIBLE_HEIGHT=temp8.d; //double

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->Vslump=temp8.d; //double

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  matprops_ptr->frict_tiny=temp8.d; //double

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  if(temp4.i!=matprops_ptr->material_count) {
    printf("in loadrun() processor %d opened file \"%s\" old had %d materials new had=%d materials, aborting\n",myid,filename,temp4.i,matprops_ptr->material_count);
    fflush(stdout);
    assert(0);
  }
  int imat;
  //skip over the material names
  for(imat=1;imat<=matprops_ptr->material_count;imat++) {
    for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
    Itemp+=temp4.i;
  }

  //bed friction
  for(imat=1;imat<=matprops_ptr->material_count;imat++) {
    for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
    if(temp8.d!=matprops_ptr->bedfrict[imat]) {
      printf("in loadrun() processor %d opened file \"%s\" for material %d new bed friction=%g [deg] old bed friction =%g [deg]\n",myid,filename,imat,temp8.d*180.0/PI,matprops_ptr->bedfrict[imat]*180.0/PI);
      fflush(stdout);
      assert(0);
    }
  }
    
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  if(temp8.d!=matprops_ptr->intfrict) {
    printf("in loadrun() processor %d opened file \"%s\" new internal friction=%g [deg] old internal friction =%g [deg]\n",myid,filename,temp8.d*180.0/PI,matprops_ptr->intfrict*180.0/PI);
    fflush(stdout);
  }

  //load the statistics
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  statprops_ptr->runid=temp4.i;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->xcen=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->ycen=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->xvar=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->yvar=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->rmean=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->area=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->vmean=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->vxmean=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->vymean=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->slopemean=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->vstar=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->realvolume=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->statvolume=temp8.d;

  //start new
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->outflowvol=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->erodedvol=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->depositedvol=temp8.d;
  //end new

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->cutoffheight=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->piler=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->hmax=temp8.d;
 
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->vmax=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->forceint=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->forcebed=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->heightifreach=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->xyifreach[0]=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->xyifreach[1]=temp8.d;

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->timereached=temp8.d;

  int i;
  for(i=0;i<4;i++) {
    for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
    statprops_ptr->xyminmax[i]=temp8.d;
  }

  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  statprops_ptr->hxyminmax=temp8.d;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  statprops_ptr->lhs.refnum=temp4.i;

  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  statprops_ptr->lhs.runid=temp4.i;

  //DISCHARGE PLANES
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  discharge_ptr->num_planes=temp4.i;

  if(discharge_ptr->num_planes>0)
    CDeAllocD2(discharge_ptr->planes); //was initialized in Read_data
  discharge_ptr->planes=CAllocD2(discharge_ptr->num_planes,10);
  
  int iplane, icol;
  for(iplane=0;iplane<discharge_ptr->num_planes;iplane++) 
    for(icol=0;icol<10;icol++) {
      for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
      discharge_ptr->planes[iplane][icol]=temp8.d;
    }

  //MAXIMUM PILE HEIGHT: A.K.A. FLOW OUTLINE
  double dxy[2], XRange[2], YRange[2];
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  dxy[0]=temp8.d;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  dxy[1]=temp8.d;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  XRange[0]=temp8.d;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  XRange[1]=temp8.d;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  YRange[0]=temp8.d;
  
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  YRange[1]=temp8.d;
  
  outline_ptr->init2(dxy, XRange, YRange);
  if(myid==0) outline_ptr->reload(matprops_ptr,statprops_ptr);

  /*  if( myid !=0 )
      MPI_Recv(&done, 1, MPI_INT, myid-1, TECTAG, MPI_COMM_WORLD, &status); */

  int ibucket;
  HashEntryPtr entryp;
  Element* EmTemp;
  Node* NdTemp;

  double doublekeyrange[2];
  double Xrange[2], Yrange[2];

  //the number of NODES in table
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  int Node_Num=temp4.i;

  //Node Hash Table doublekeyrange
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  doublekeyrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  doublekeyrange[1]=temp8.d;
  
  //Node Table Size
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  int NODE_TABLE_SIZE=temp4.i;

  //Node Hash Table Xrange 
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Xrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Xrange[1]=temp8.d;

  //Node Hash Table Yrange 
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Yrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Yrange[1]=temp8.d;

  //recreate the node hashtable
  *NodeTable  = new HashTable(doublekeyrange, NODE_TABLE_SIZE, 2017, Xrange, Yrange,1); 


  //the number of ELEMENTS in table 
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  int Elem_Num=temp4.i;

  //Elem Hash Table doublekeyrange
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  doublekeyrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  doublekeyrange[1]=temp8.d;

  //Elem Table Size
  for(itemp=0;itemp<4;itemp++) temp4.c[itemp]=header[Itemp++];
  int ELEM_TABLE_SIZE=temp4.i; //NBUCKETS

  //Elem Hash Table Xrange 
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Xrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Xrange[1]=temp8.d;

  //Elem Hash Table Yrange 
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Yrange[0]=temp8.d;
  for(itemp=0;itemp<8;itemp++) temp8.c[itemp]=header[Itemp++];
  Yrange[1]=temp8.d;

  if(Itemp!=numchar) {
    printf("DANGER!!! load_run() uses %d instead of %d bytes, load_run not compattible with save_run that wrote the file %s on processor %d\n",Itemp,numchar,filename,myid);
    fflush(stdout);
    assert(0);
  }

  //recreate the element hashtable
  *ElemTable  = new HashTable(doublekeyrange, ELEM_TABLE_SIZE, 503, Xrange, Yrange,1); 

  /*****************************************************************/

  //read in all the nodes
  Node* NodeP;
  int inode=-1;
  for(inode=0;inode<Node_Num;inode++) {
    if(inode<0)
      printf("inode=%d\n",inode);
    NodeP = new Node(fp,matprops_ptr);
    (*NodeTable)->add(NodeP->pass_key(),NodeP);
  }
  printf("inode=%d Node_Num=%d NodeP=%u\n",inode,Node_Num,NodeP);

  unsigned tempkey[2];
  tempkey[0]=3590592512; tempkey[1]=0;


  //read in all the elements
  Element* ElemP;
  int ielem=0;
  int maxgen=0;
  for(ielem=0;ielem<Elem_Num;ielem++) {
    if(ielem<0)
      printf("ielem=%d\n",ielem);
    ElemP = new Element(fp,*NodeTable,matprops_ptr,myid);
    (*ElemTable)->add(ElemP->pass_key(),ElemP);
    if(ElemP->get_gen()>maxgen) maxgen=ElemP->get_gen();
    //if((*(ElemP->pass_key()+0)==tempkey[0])&&(*(ElemP->pass_key()+1)==tempkey[1]))
    //printf("YAAAAAAAAAAAADAAAAAAAAAAAAAAAAA\n");

  }
  printf("ielem=%d Elem_Num=%d ElemP=%u\n",ielem,Elem_Num,ElemP);

  double dx=*(ElemP->get_dx()+0), dy=*(ElemP->get_dx()+1); 
  if(dx<dy) dx=dy;

  REFINE_LEVEL=ElemP->get_gen()+
    ceil(log(dx*
	     (matprops_ptr->number_of_cells_across_axis)/
	     (matprops_ptr->smallest_axis)
	     )/
	 log(2.0));

  if(maxgen>REFINE_LEVEL) REFINE_LEVEL=maxgen;

  maxgen=REFINE_LEVEL;
  MPI_Allreduce(&maxgen,&REFINE_LEVEL,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  //calc_d_gravity
  int no_of_elm_buckets = (*ElemTable)->get_no_of_buckets();
  for(ibucket=0; ibucket<no_of_elm_buckets; ibucket++) {          
    entryp = *((*ElemTable)->getbucketptr() + ibucket);
    while(entryp){	
      EmTemp = (Element*)entryp->value; assert(EmTemp);

      if(EmTemp->get_adapted_flag()>0) EmTemp->calc_d_gravity(*ElemTable);
      
      entryp = entryp->next;
    }
  }

  move_data(numprocs, myid, *ElemTable, *NodeTable,timeprops_ptr);

  return(1);
}

/*************************************************************************/
/*************************************************************************/

void saverun(HashTable** NodeTable, int myid, int numprocs, 
	     HashTable** ElemTable, MatProps* matprops_ptr,  
	     TimeProps* timeprops_ptr, MapNames *mapnames_ptr, 
	     int adaptflag, int order_flag, 
	     StatProps *statprops_ptr, DISCHARGE *discharge_ptr,
	     OutLine* outline_ptr, int *savefileflag) {

  move_data(numprocs, myid, *ElemTable, *NodeTable,timeprops_ptr);
  
  char filename[64], file2[64], file3[64], file4[64], file5[64], file6[64];
  sprintf(filename,"restart%04d.%1d",myid,(*savefileflag+1)%2);
  sprintf(file2,"header%04d.%1d",myid,(*savefileflag+1)%2);
  sprintf(file3,"node%04d.%1d",myid,(*savefileflag+1)%2);
  sprintf(file4,"elem%04d.%1d",myid,(*savefileflag+1)%2);
  sprintf(file5,"header%04d.txt.%1d",myid,(*savefileflag+1)%2); 
  FILE *fp=fopen(filename,"wb");
#ifdef DEBUGSAVESPLIT
  FILE *fp2=fopen(file2,"wb");
#endif
#ifdef DEBUGSAVEHEADER
  FILE *fp3=fopen(file5,"w");
#endif
  sprintf(file6,"elem%04d.txt.%1d",myid,(*savefileflag+1)%2);
  FILE* fptxt=fopen(file6,"w");

  fourbytes  temp4;
  eightbytes temp8;
  char header[NUM_CHAR_IN_SAVE_HEADER];  

  //need to change number and use 4 bytes and 8bytes unions
  int Itemp=4, itemp, imat, i; //first 4 bytes is left for number of bytes in rest of header takes up

  //label the restart file version... for safety check
  char versionstring[64];
  sprintf(versionstring,"one-phase:20060801");
  //save how many characters are in version string
  temp4.i=strlen(versionstring); 
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
  //save version string
  for(itemp=0;itemp<temp4.i;itemp++) header[Itemp++]=versionstring[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"versionstring: %d: \"%s\"\n",strlen(versionstring),versionstring);
#endif

  //processor info... for safety check
  temp4.i=myid;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

  temp4.i=numprocs;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

  //number of equations...for safety check
  temp4.i=EQUATIONS;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

  //flags
  temp4.i=adaptflag;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

  temp4.i=order_flag;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"myid=%d numprocs=%d\nEQUATIONS=%d\nadaptflag=%d order_flag=%d\n",myid,numprocs,EQUATIONS,adaptflag,order_flag);
#endif

  //TIME info
  temp4.i=timeprops_ptr->iter;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

  temp8.d=timeprops_ptr->time; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=timeprops_ptr->TIME_SCALE; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"iter=%d time=%g TIME_SCALE=%g\n",timeprops_ptr->iter,timeprops_ptr->time,timeprops_ptr->TIME_SCALE);
#endif

  //GIS file info
  temp4.i=strlen(mapnames_ptr->gis_main);
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
  for(itemp=0;itemp<temp4.i;itemp++)
    header[Itemp++]=mapnames_ptr->gis_main[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"gis_main: %d: \"%s\"\n",temp4.i,mapnames_ptr->gis_main);
#endif

  temp4.i=strlen(mapnames_ptr->gis_sub);
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
  for(itemp=0;itemp<temp4.i;itemp++)
    header[Itemp++]=mapnames_ptr->gis_sub[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"gis_sub: %d: \"%s\"\n",temp4.i,mapnames_ptr->gis_sub);
#endif

  temp4.i=strlen(mapnames_ptr->gis_mapset);
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
  for(itemp=0;itemp<temp4.i;itemp++)
    header[Itemp++]=mapnames_ptr->gis_mapset[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"gis_mapset: %d: \"%s\"\n",temp4.i,mapnames_ptr->gis_mapset);
#endif

  temp4.i=strlen(mapnames_ptr->gis_map);
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
  for(itemp=0;itemp<temp4.i;itemp++) 
    header[Itemp++]=mapnames_ptr->gis_map[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"gis_map: %d: \"%s\"\n",temp4.i,mapnames_ptr->gis_map);
#endif

  temp4.i=mapnames_ptr->extramaps;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"extramaps=%d:\n",temp4.i);
#endif

  //material and scaling info
  temp8.d=matprops_ptr->LENGTH_SCALE; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp8.d=matprops_ptr->HEIGHT_SCALE; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp8.d=matprops_ptr->GRAVITY_SCALE; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp8.d=matprops_ptr->MAX_NEGLIGIBLE_HEIGHT; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp8.d=matprops_ptr->Vslump; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp8.d=matprops_ptr->frict_tiny; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  

  temp4.i=matprops_ptr->material_count;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"LENGTH_SCALE=%g HEIGHT_SCALE=%g GRAVITY_SCALE=%g\n",
	  matprops_ptr->LENGTH_SCALE,matprops_ptr->HEIGHT_SCALE,
	  matprops_ptr->GRAVITY_SCALE);
  fprintf(fp3,"max_neg_height=%g Vslump=%g frict_tiny=%g\nmat_count=%d\n",
	  matprops_ptr->MAX_NEGLIGIBLE_HEIGHT,matprops_ptr->Vslump,
	  matprops_ptr->frict_tiny,matprops_ptr->material_count);
#endif

  //save the material names
  for(imat=1;imat<=matprops_ptr->material_count;imat++) {
    //count and save the NUMBER of characters in this name
    temp4.i=strlen(matprops_ptr->matnames[imat]);
    for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];

    //save the characters in this name
    for(itemp=0;itemp<temp4.i;itemp++) 
      header[Itemp++]=matprops_ptr->matnames[imat][itemp];

#ifdef DEBUGSAVEHEADER
    fprintf(fp3,"material %d: %d: \"%s\"\n",imat,temp4.i,matprops_ptr->matnames[imat]);
#endif
  }

  //bed friction
  for(imat=1;imat<=matprops_ptr->material_count;imat++) {
    temp8.d=matprops_ptr->bedfrict[imat]; //double
    for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
    fprintf(fp3,"material %d: bedfrict=%g\n",imat,temp8.d);
#endif
  }

  temp8.d=matprops_ptr->intfrict; //double
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];  
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"intfrict=%g\n",temp8.d);
#endif

  //save the statistics
  temp4.i=statprops_ptr->runid;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"runid=%d\n",statprops_ptr->runid);
#endif

  temp8.d=statprops_ptr->xcen;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->ycen;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"xcen=%g ycen=%g\n",statprops_ptr->xcen,statprops_ptr->ycen);
#endif

  temp8.d=statprops_ptr->xvar;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->yvar;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"xvar=%g yvar=%g\n",statprops_ptr->xvar,statprops_ptr->yvar);
#endif

  temp8.d=statprops_ptr->rmean;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->area;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"rmean=%g area=%g\n",statprops_ptr->rmean,statprops_ptr->area);
#endif

  temp8.d=statprops_ptr->vmean;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->vxmean;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->vymean;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"vmean=%g vxmean=%g vymean=%g\n",statprops_ptr->vmean,statprops_ptr->vxmean,statprops_ptr->vymean);
#endif

  temp8.d=statprops_ptr->slopemean;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->vstar;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"slopemean=%g vstar=%g\n",statprops_ptr->slopemean,statprops_ptr->vstar);
#endif

  temp8.d=statprops_ptr->realvolume;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->statvolume;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->outflowvol;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->erodedvol;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->depositedvol;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->cutoffheight;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"realvolume=%g statvolume=%g outflowvol=%g erodedvol=%g depositedvol=%g cutoffheight=%g\n",
	  statprops_ptr->realvolume,statprops_ptr->statvolume,
	  statprops_ptr->outflowvol,statprops_ptr->erodedvol,
	  statprops_ptr->depositedvol,statprops_ptr->cutoffheight);
#endif

  temp8.d=statprops_ptr->piler;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->hmax;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
 
  temp8.d=statprops_ptr->vmax;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"piler=%g hmax=%g vmax=%g\n",statprops_ptr->piler,statprops_ptr->hmax,statprops_ptr->vmax);
#endif

  temp8.d=statprops_ptr->forceint;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->forcebed;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"forceint=%g forcebed=%g\n",statprops_ptr->forceint,statprops_ptr->forcebed);
#endif

  temp8.d=statprops_ptr->heightifreach;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->xyifreach[0];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->xyifreach[1];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];

  temp8.d=statprops_ptr->timereached;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"to reach: height=%g x=%g y=%g time=%g\nminmax: xy={",
	  statprops_ptr->heightifreach,statprops_ptr->xyifreach[0],
	  statprops_ptr->xyifreach[1],statprops_ptr->timereached);
#endif
  for(i=0;i<4;i++) {
    temp8.d=statprops_ptr->xyminmax[i];
    for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
    fprintf(fp3," %g",temp8.d);
#endif
  }

  temp8.d=statprops_ptr->hxyminmax;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"}  h=%g\n",temp8.d);
#endif

  temp4.i=statprops_ptr->lhs.refnum;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"lhs: refnum=%d",temp4.i);
#endif

  temp4.i=statprops_ptr->lhs.runid;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3," runid=%d\n",temp4.i);
#endif  


  //DISCHARGE PLANES
  temp4.i=discharge_ptr->num_planes;
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"discharge: numplanes=%d",temp4.i);
#endif

  int iplane, icol;
  for(iplane=0;iplane<discharge_ptr->num_planes;iplane++) {
#ifdef DEBUGSAVEHEADER
    fprintf(fp3,"\nplane[%d]=",iplane);
#endif
    for(icol=0;icol<10;icol++) {
      temp8.d=discharge_ptr->planes[iplane][icol];
#ifdef DEBUGSAVEHEADER
      fprintf(fp3," %g",temp8.d);
#endif
      for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
    }
  }

  //MAXIMUM PILE HEIGHT: A.K.A. FLOW OUTLINE
  temp8.d=outline_ptr->dx;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
    
  temp8.d=outline_ptr->dy;
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
    
  temp8.d=outline_ptr->xminmax[0];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
  
  temp8.d=outline_ptr->xminmax[1];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
  
  temp8.d=outline_ptr->yminmax[0];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
    
  temp8.d=outline_ptr->yminmax[1];
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"\nMaxPileHeight/Outline: dx=%g dy=%g  (%g<=x<=%g) (%g<=y<=%g)\n",
	  outline_ptr->dx,outline_ptr->dy,
	  outline_ptr->xminmax[0],outline_ptr->xminmax[1],
	  outline_ptr->yminmax[0],outline_ptr->yminmax[1]);
#endif

  //node and element HASHTABLES

  int ibucket;
  int no_of_elm_buckets = (*ElemTable)->get_no_of_buckets();
  int no_of_node_buckets = (*NodeTable)->get_no_of_buckets();
  HashEntryPtr entryp;
  Element* EmTemp;
  Node* NdTemp;




  //count the number of NODES to save
  temp4.i=0;
  for(ibucket=0;ibucket<no_of_node_buckets;ibucket++) {
    entryp = *((*NodeTable)->getbucketptr() + ibucket);
    while(entryp) {
      temp4.i++;
      entryp = entryp->next;
    }
  }
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"\nnumnodes=%d\n",temp4.i);
#endif

  //Node Hash Table doublekeyrange
  temp8.d=*((*NodeTable)->get_doublekeyrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"node: doublekeyrange=[%20g,",temp8.d);
#endif

  temp8.d=*((*NodeTable)->get_doublekeyrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%20g] ",temp8.d);
#endif

  //get node table size;
  temp4.i=(*NodeTable)->get_nbuckets();
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"table_size=%d\n",temp4.i);
#endif

  //Node Hash Table Xrange 
  temp8.d=*((*NodeTable)->get_Xrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"node: xrange=[%g ",temp8.d);
#endif

  temp8.d=*((*NodeTable)->get_Xrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%g] ",temp8.d);
#endif

  //Node Hash Table Yrange 
  temp8.d=*((*NodeTable)->get_Yrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"yrange=[%g ",temp8.d);
#endif

  temp8.d=*((*NodeTable)->get_Yrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%g]\n",temp8.d);
#endif




  //count the number of ELEMENTS to save
  temp4.i=0;
  for(ibucket=0; ibucket<no_of_elm_buckets; ibucket++) {          
    entryp = *((*ElemTable)->getbucketptr() + ibucket);
    while(entryp){	
      EmTemp = (Element*)entryp->value;
      assert(EmTemp);
      if(EmTemp->get_myprocess()==myid) {
	if(EmTemp->get_adapted_flag()>0) {
	  temp4.i++;
	}
      }
      entryp = entryp->next;
    }
  }
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"numelem=%d\n",temp4.i);
#endif

 //Element Hash Table MinKeys
  //Node Hash Table doublekeyrange
  temp8.d=*((*ElemTable)->get_doublekeyrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"elem: doublekeyrange=[%20g,",temp8.d);
#endif

  temp8.d=*((*ElemTable)->get_doublekeyrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%20g] ",temp8.d);
#endif

  //get elem table size;
  temp4.i=(*ElemTable)->get_nbuckets();
  for(itemp=0;itemp<4;itemp++) header[Itemp++]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"table_size=%d\n",temp4.i);
#endif

  //Element Hash Table Xrange 
  temp8.d=*((*ElemTable)->get_Xrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"elem: xrange=[%g ",temp8.d);
#endif

  temp8.d=*((*ElemTable)->get_Xrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%g] ",temp8.d);
#endif

  //Element Hash Table Yrange 
  temp8.d=*((*ElemTable)->get_Yrange()+0);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"yrange=[%g ",temp8.d);
#endif

  temp8.d=*((*ElemTable)->get_Yrange()+1);
  for(itemp=0;itemp<8;itemp++) header[Itemp++]=temp8.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"%g]\n",temp8.d);
#endif
  header[Itemp]=0; //this should not do anything, should not be saved, just making sure the "character string" terminates here but it probably already has so it is not possible to print the WHOLE thing with a %s

  if(Itemp>=NUM_CHAR_IN_SAVE_HEADER) {
    printf("In restart.C while saving, ran out of room in header;\nItemp=%d>=NUM_CHAR_IN_SAVE_HEADER=%d\nABORTING!!!!\n",Itemp,NUM_CHAR_IN_SAVE_HEADER);
    fflush(stdout);
    exit(1);
  }

  //tell the file loader how big the header is
  temp4.i=Itemp-4; //-4 because they will read the first 4 bytes to find out how big the rest is
  for(itemp=0;itemp<4;itemp++) header[itemp]=temp4.c[itemp];
#ifdef DEBUGSAVEHEADER
  fprintf(fp3,"numbytes=%d",temp4.i);
#endif
  //actually write the header
  fwrite(header,sizeof(char),Itemp,fp);
#ifdef DEBUGSAVESPLIT
  fwrite(header,sizeof(char),Itemp,fp2);
  fflush(fp2); fclose(fp2); fp2=fopen(file3,"wb");
#endif
#ifdef DEBUGSAVEHEADER
  fflush(fp3); //text header is written
  fclose(fp3);
#endif  


  /*****************************************************************/

  //what remains to be done is loop through the nodes and then the 
  //elements calling the "save" member functions for each one that
  //should be written

  for(ibucket=0;ibucket<no_of_node_buckets;ibucket++) {
    entryp = *((*NodeTable)->getbucketptr() + ibucket);
    while(entryp) {
      NdTemp = (Node*) (entryp->value);

      NdTemp->save_node(fp);
#ifdef DEBUGSAVESPLIT
      NdTemp->save_node(fp2);
#endif
      entryp = entryp->next;
    }
  }
#ifdef DEBUGSAVESPLIT
  fflush(fp2); fclose(fp2); fp2=fopen(file4,"wb");
#endif

  for(ibucket=0; ibucket<no_of_elm_buckets; ibucket++) {          
    entryp = *((*ElemTable)->getbucketptr() + ibucket);
    while(entryp) {	
      EmTemp = (Element*)entryp->value; assert(EmTemp);

      if(EmTemp->get_myprocess()==myid) {
	if(EmTemp->get_adapted_flag()>0) {
	  EmTemp->save_elem(fp,fptxt);
#ifdef DEBUGSAVESPLIT
	  EmTemp->save_elem(fp2,fptxt);
#endif
	}
      }
      entryp = entryp->next;
    }
  }

#ifdef DEBUGSAVESPLIT
  fflush(fp2); fclose(fp2); 
#endif
  fflush(fp);  fclose(fp);  fclose(fptxt);
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  //flip savefileflag to indicate that THIS is the most recent completed file
  *savefileflag=(*savefileflag+1)%2;

}
