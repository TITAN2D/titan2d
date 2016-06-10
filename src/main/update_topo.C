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
 * $Id: update_topo.C 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"
#define WORKDIR "."

int update_topo(ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, int myid, int nump, MatProps* matprops,
                TimeProps* timeprops, MapNames* mapnames)
{
    
    char gis_update_map[100] = "\0";
    char yada[256];
    
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;
    
    if(myid == 0)
    {
        sprintf(yada, "%s/%s", WORKDIR, "updatetopo.txt");
        ifstream inp(yada, ios::in);
        if(!inp.fail())
        {
            inp >> gis_update_map;
            inp.close();
        }
    }

#ifdef USE_MPI
    if(nump > 0)
        MPI_Bcast(gis_update_map, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif //USE_MPI
    //printf("update_topo.C 5\n"); fflush(stdout);
    
    //check for a update map to open
    if(strlen(gis_update_map) > 0)
    {
        double tick = MPI_Wtime();
        Delete_GIS_data();
        Initialize_GIS_data(mapnames->gis_main.c_str(), mapnames->gis_sub.c_str(), mapnames->gis_mapset.c_str(), gis_update_map);
        
        /* update the nodes */
        int num_buckets = HT_Node_Ptr->get_no_of_buckets();
        //visit every bucket
        //@NodesSingleLoop
        for(int i = 0; i < HT_Node_Ptr->elenode_.size(); i++)
        {
            if(HT_Node_Ptr->status_[i]>=0)
                HT_Node_Ptr->elenode_[i].elevation(matprops);
        }
        
        /* update the elements */
        //@ElementsBucketDoubleLoop
        for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
        {
            for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
            {
                Element *EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
                
                if(EmTemp->adapted_flag() > 0)
                {
                    //update the topography in same order as in element2.C
                    //double eldif=EmTemp->get_elevation();
                    EmTemp->calc_topo_data(matprops);
                    //eldif=(EmTemp->get_elevation()-eldif)*matprops->LENGTH_SCALE;
                    //if(fabs(eldif)>1.0) printf("update_topo() after-before=%g\n",eldif);
                    
                    EmTemp->calc_gravity_vector(matprops);
                    EmTemp->calc_d_gravity(HT_Elem_Ptr);
                }
            }
        }	  //closes: for(int ibucket=0; ibucket<num_buckets; ibucket++)
        double tock = MPI_Wtime() - tick;
        //long tock=(long) time(NULL);
#ifdef USE_MPI
        if(nump > 1)
        {
            /* long tempin[2], tempout[2];
             tempin[0]=-tick;
             tempin[1]= tock;
             MPI_Reduce(tempin, tempout, 2, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
             tick=-tempout[0];
             tock=tempout[1]; */
            MPI_Reduce(&tock, &tick, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            tock = tick;
            
        }
#endif //USE_MPI
        
        if(myid == 0)
        {
            sprintf(yada, "%s/%s", WORKDIR, "updatetopo.done");
            FILE *fp = fopen(yada, "w");
            fprintf(fp, "%s", gis_update_map);
            fclose(fp);
            
            sprintf(yada, "%s/%s", WORKDIR, "updatetopo.time");
            fp = fopen(yada, "a");
            int hours, minutes;
            double seconds;
            timeprops->chunktime(&hours, &minutes, &seconds);
            fprintf(fp, "%s timestep=%d simtime=(%02d:%02d%g) (hrs:min:sec) timing=%g (sec)\n", gis_update_map,
                    timeprops->iter, hours, minutes, seconds, tock);
            fclose(fp);
        }
        
        return (1); //successfull updated topography
    } //closes: if(strlen(gis_update_map>0))
    else
        return (0); //no update gis map to open
}

