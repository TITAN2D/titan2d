#include "../header/hpfem.h"

void grass_sites_header_output(TimeProps* timeprops, const char * output_prefix)
{
    char filename[24];
    
    sprintf(filename, "%s/grass_sites%08d.hdr", output_prefix, timeprops->iter);
    FILE *fp = fopen(filename, "w");
    
    fprintf(fp, "name|titangrassout.site\n");
    fprintf(fp, "desc|time=%E iter=%d\n", timeprops->timesec(), timeprops->iter);
    fprintf(fp, "labels|east north elevation pile_height vx vy x_mom y_mom\n");
    fprintf(fp, "form||||\%\%\%\%\%\n");
    
    fclose(fp);
    
    return;
}

void grass_sites_proc_output(ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, int myid, MatProps* matprops,
                             TimeProps* timeprops, const char * output_prefix)
{
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;
    
    double velocity_scale = sqrt(matprops->scale.length * (matprops->scale.gravity));
    double momentum_scale = matprops->scale.height * velocity_scale; // scaling factor for the momentums
            
    char filename[24];
    
    sprintf(filename, "%s/grass_sites%08d.%03d", output_prefix, timeprops->iter, myid);
    FILE *fp = fopen(filename, "w");
    
    /***************************************************************/
    /* print out the variables at the bubble node of every element */
    /* on this processor                                           */
    /***************************************************************/

    //check every bucket for elements
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element *EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            
            //if the current element is an active one on this processor
            if(EmTemp->adapted_flag() > 0)
            {
	        // HT_Node_Ptr->lookup(EmTemp->node_key(8)) returns a null pointer.
	        // Per header/hashtab.h, the n9th (8 out of 0->8) node is the bubble node
	        // it's key is not stored separately since it has the same key as the element
                //Node *NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(8));
	        Node *NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->key());
                
                double pile_height = EmTemp->state_vars(0) * (matprops->scale.height);
                double x_mom = EmTemp->state_vars(1) * momentum_scale;
                double y_mom = EmTemp->state_vars(2) * momentum_scale;
                double VxVy[2];
                EmTemp->eval_velocity(0.0, 0.0, VxVy);
                //double vx=(pile_height>GEOFLOW_TINY)?x_mom/pile_height:0.0;
                //double vy=(pile_height>GEOFLOW_TINY)?y_mom/pile_height:0.0;
                double vx = (pile_height > GEOFLOW_TINY) ? VxVy[0] * velocity_scale : 0.0;
                double vy = (pile_height > GEOFLOW_TINY) ? VxVy[1] * velocity_scale : 0.0;
                
                //print x,y,z,h,Vx,Vy,h*Vx,h*Vy
                fprintf(fp, "%g|%g|%g|%%%g %%%g %%%g %%%g %%%g\n",
                        (NodeTemp->coord(0)) * (matprops->scale.length), //x
                        (NodeTemp->coord(1)) * (matprops->scale.length), //y
                        NodeTemp->elevation() * (matprops->scale.length), //elevation
                        pile_height, vx, vy, x_mom, y_mom);
            }
        }
    }
    
    fclose(fp);
    
    return;
}

