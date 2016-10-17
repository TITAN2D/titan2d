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
 * $Id: init_piles.C 229 2012-02-26 22:11:07Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"

#include "../header/titan_simulation.h"

int get_elem_elev(NodeHashTable *HT_Node_Ptr, MatProps *matprops, Element *EmTemp, double &elevation);

void print_grid(ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, MatProps* matprops)
{
    
    FILE *fp = fopen("gridplot00.txt", "w");
    
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;

    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element *EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            
            if(EmTemp->adapted_flag() > 0)
            {
                double elevation;
                get_elem_elev(HT_Node_Ptr, matprops, EmTemp, elevation);
                
                fprintf(fp, "%20.14g %20.14g %20.14g\n", (EmTemp->coord(0)) * (matprops)->scale.length,
                        (EmTemp->coord(1)) * (matprops)->scale.length, elevation);
            }
        }
    }
    
    fclose(fp);
    
    assert(0);
    
    return;
}


void cxxTitanSimulation::init_piles()
{

    MatProps* matprops_ptr = get_matprops();
    FluxProps* fluxprops_ptr = get_fluxprops();
    TimeProps* timeprops_ptr = get_timeprops();
    StatProps* statprops_ptr = get_statprops();

    ElementsHashTable* HT_Elem_Ptr=get_HT_Elem();
    NodeHashTable* HT_Node_Ptr=get_HT_Node();

    PileProps* pileprops_ptr=get_pileprops();

    unsigned nodes[9][KEYLENGTH], *node_key;
    
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;
    

    PileProps::PileType pile_type= pileprops_ptr->get_default_piletype();

    int i;
    bool allPilesAreElliptical=true;
    for(i=0;i<pileprops_ptr->numpiles;i++)
    {
        if(!(pileprops_ptr->pile_type[i] == PileProps::PARABALOID || pileprops_ptr->pile_type[i] == PileProps::CYLINDER))
            allPilesAreElliptical=false;
    }
    if(pileprops_ptr->numpiles>0)pile_type= pileprops_ptr->pile_type[0];



    if(!adapt)
    {
        H_adapt_to_level(HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr, pileprops_ptr, fluxprops_ptr, timeprops_ptr, REFINE_LEVEL);
		HT_Elem_Ptr->updateLocalElements();
		HT_Elem_Ptr->updateNeighboursIndexes();
    }

    if(allPilesAreElliptical)
    {
        if(adapt)
            initial_H_adapt(HT_Elem_Ptr, HT_Node_Ptr, 0, matprops_ptr, pileprops_ptr, fluxprops_ptr, timeprops_ptr, 4);
    }
    else
    {
        printf("It seems this type of piles have hardcoded coordinates\n");
        assert(0);
        //@ElementsBucketDoubleLoop
        for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
        {
            for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
            {
                Element *EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
                
                if(EmTemp->adapted_flag() > 0)
                {
                    //put in the pile height right here...
                    double pile_height = 0.0;
                    double radius_sq;
                    switch (pile_type)
                    {

                        case PileProps::PLANE:
                            radius_sq = pow(EmTemp->coord(0) - 76., 2) + pow(EmTemp->coord(1) - 80., 2);
                            if(radius_sq < 35.)
                                pile_height = 10 * (1. - radius_sq / 35.);
                            break;
                        case PileProps::CASITA:
                            radius_sq = pow(EmTemp->coord(0) - 504600., 2) + pow(EmTemp->coord(1) - 1402320., 2);
                            if(radius_sq < 30000.)
                                pile_height = 15 * (1. - radius_sq / 30000.);
                            break;
                        case PileProps::POPO: //popo topo
                            radius_sq = pow(EmTemp->coord(0) - 537758. / matprops_ptr->scale.length, 2)
                                    + pow(EmTemp->coord(1) - 2100910. / matprops_ptr->scale.length, 2);
                            if(radius_sq < (10000. / matprops_ptr->scale.length))
                                pile_height = 1. - radius_sq / (10000. / matprops_ptr->scale.length);
                            break;
                        case PileProps::ID1: // iverson and denlinger experiments I -- as pictured

                            if(EmTemp->coord(0) < 53.345 && EmTemp->coord(0) > 45.265 && EmTemp->coord(1) > -10. && EmTemp->coord(1) < 300.)
                            {
                                if(EmTemp->coord(0) < 51.148)
                                    pile_height = 3.5912 * (1.0 - (51.148 - EmTemp->coord(0)) / 5.8832);
                                else
                                    pile_height = 3.59 * (53.345 - EmTemp->coord(0)) / 2.1967;
                                if(pile_height < 0)
                                    pile_height = 0;
                            }
                            break;
                        case PileProps::ID2: //iverson and denlinger experiments II -- 90 angle with plane
                            if(EmTemp->coord(0) < 53.345 / matprops_ptr->scale.length && EmTemp->coord(0)
                                    > 46.45 / matprops_ptr->scale.length)
                                pile_height = 4.207255
                                        * (1.0 - (53.345 / matprops_ptr->scale.length - EmTemp->coord(0)) / 6.895
                                                * matprops_ptr->scale.length);
                            break;
                        default:
                            printf("Danger no recognized pile type defined in init_piles.C\n");
                            assert(0);
                    }
                    EmTemp->put_height(pile_height);

                }
            }
        }
    } //end "#if defined PARABALOID || defined CYLINDER"
    move_data(numprocs, myid, HT_Elem_Ptr, HT_Node_Ptr, timeprops_ptr);
    
    //update temporary arrays of elements/nodes pointers
    HT_Node_Ptr->flushNodeTable();
    HT_Elem_Ptr->flushElemTable();
    HT_Elem_Ptr->updateLocalElements();
    HT_Elem_Ptr->updateNeighboursIndexes();
    
    slopes(HT_Elem_Ptr, HT_Node_Ptr, matprops_ptr);
    
    /* initial calculation of actual volume on the map */

    double realvolume = 0.0, depositedvol = 0.0, forcebed = 0.0, meanslope = 0.0;
    double epsilon[DIMENSION];
    for(i=0;i<DIMENSION;i++)
        epsilon[i]=matprops_ptr->scale.epsilon;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            Element* Curr_El  = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(Curr_El->adapted_flag() > 0)
            { //if this is a refined element don't involve!!!
                double dvol = Curr_El->dx(0) * Curr_El->dx(1) * Curr_El->state_vars(0);
                realvolume += dvol;
                Curr_El->set_kactxy(epsilon);

                Curr_El->calc_stop_crit(matprops_ptr,integrator);
                if(Curr_El->stoppedflags() == 2)
                    depositedvol += dvol;

                double resolution = 0, xslope = 0, yslope = 0;
                Get_max_resolution(&resolution);
                Get_slope(resolution, Curr_El->coord(0) * matprops_ptr->scale.length,
                          Curr_El->coord(1) * matprops_ptr->scale.length, xslope, yslope);
                double slope = sqrt(xslope * xslope + yslope * yslope);

                forcebed += dvol * 9.8 / sqrt(1.0 + slope * slope)
                        * tan(matprops_ptr->bedfrict[Curr_El->material()]);

            }
        }
    }
    
    double tempin[3], tempout[3];
    tempin[0] = realvolume;
    tempin[1] = forcebed;
    tempin[2] = depositedvol;
#ifdef USE_MPI
    MPI_Reduce(tempin, tempout, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else //USE_MPI
    for(int i9=0;i9<3;++i9)tempout[i9]=tempin[i9];
#endif //USE_MPI
    
    statprops_ptr->realvolume = tempout[0] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length) * (matprops_ptr->scale.length);
    statprops_ptr->outflowvol = 0.0;
    statprops_ptr->erodedvol = 0.0;
    statprops_ptr->depositedvol = tempout[2] * (matprops_ptr->scale.height) * (matprops_ptr->scale.length)
                              * (matprops_ptr->scale.length);
    
    statprops_ptr->forceint = 0.0;
    statprops_ptr->forcebed = tempout[1] / tempout[0];
    
    return;
}

