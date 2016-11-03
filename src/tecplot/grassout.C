#include "../header/hpfem.h"

void grass_sites_header_output(TimeProps* timeprops){
  char filename[24];  

  sprintf(filename,"grass_sites%08d.hdr",timeprops->iter);
  FILE *fp=fopen(filename,"w");

  fprintf(fp,"name|titangrassout.site\n");
  fprintf(fp,"desc|time=%E iter=%d\n",timeprops->timesec(),timeprops->iter);
  fprintf(fp,"labels|east north elevation pile_height vx vy x_mom y_mom\n");
  fprintf(fp,"form||||\%\%\%\%\%\n");

  fclose(fp);

  return;
}

void grass_sites_proc_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			     int myid, MatProps* matprops, 
			     TimeProps* timeprops)
{
  int ielembucket;
  int numelembucket=HT_Elem_Ptr->get_no_of_buckets();

  double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE));
  double momentum_scale = matprops->HEIGHT_SCALE * velocity_scale; // scaling factor for the momentums


  char filename[24];

  sprintf(filename,"grass_sites%08d.%03d",timeprops->iter,myid);
  FILE *fp=fopen(filename,"w");

  /***************************************************************/
  /* print out the variables at the bubble node of every element */
  /* on this processor                                           */
  /***************************************************************/

  //check every bucket for elements
  for(ielembucket=0; ielembucket<numelembucket; ielembucket++){
    HashEntry *entryp = *(HT_Elem_Ptr->getbucketptr() + ielembucket);

    //check every element in a bucket
    while(entryp){

      Element *EmTemp = (Element*)entryp->value;
      assert(EmTemp);

      //if the current element is an active one on this processor
      if(EmTemp->get_adapted_flag()>0) {

	unsigned *nodes = EmTemp->getNode();
	double *state_vars=EmTemp->get_state_vars();

	Node *NodeTemp = (Node*) HT_Node_Ptr->lookup(nodes+8*KEYLENGTH);

	double pile_height=state_vars[0]*(matprops->HEIGHT_SCALE);
	double x_mom=state_vars[1]*momentum_scale;
	double y_mom=state_vars[2]*momentum_scale;
	double VxVy[2];
	EmTemp->eval_velocity(0.0,0.0,VxVy);
	//double vx=(pile_height>GEOFLOW_TINY)?x_mom/pile_height:0.0;
	//double vy=(pile_height>GEOFLOW_TINY)?y_mom/pile_height:0.0;
	double vx=(pile_height>GEOFLOW_TINY)?VxVy[0]*velocity_scale:0.0;
	double vy=(pile_height>GEOFLOW_TINY)?VxVy[1]*velocity_scale:0.0;
	
	//print x,y,z,h,Vx,Vy,h*Vx,h*Vy
	fprintf(fp,"%g|%g|%g|%%%g %%%g %%%g %%%g %%%g\n",
		(*(NodeTemp->get_coord()))*(matprops->LENGTH_SCALE), //x
		(*(NodeTemp->get_coord()+1))*(matprops->LENGTH_SCALE), //y
		NodeTemp->get_elevation()*(matprops->LENGTH_SCALE), //elevation
		pile_height,vx,vy,x_mom,y_mom);}

      entryp = entryp->next;}}

  fclose(fp);

  return;
}


