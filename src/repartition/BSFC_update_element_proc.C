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
 * $Id: BSFC_update_element_proc.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "./repartition_BSFC.h"

// this routine updates the processor that an element is assigned to

void BSFC_update_element_proc(int myid, int numprocs, ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr,
                              BSFC_VERTEX_PTR sfc_vert_ptr)
{
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;
    
    int i, j, k;
    Element* EmTemp;
    
    j = 0;
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
            {
                if(EmTemp->key() == (unsigned) 74470341)
                    k = i;
                if(EmTemp->new_old() == 1)
                {  //only this element
                    if(EmTemp->myprocess() != sfc_vert_ptr[j].destination_proc)
                    { // this element will get moved to a new processor
                        EmTemp->set_myprocess(sfc_vert_ptr[j].destination_proc);
                    }
                    j++;
                }
                else if(EmTemp->new_old() > 1)
                {  //multiple elements connected by a constrained node
                   //check for constrained nodes on the vertex nodes
                    k = 4;
                    while (k < 8)
                    {
                        int ll = 1;
                        Node* ndtemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(k));
                        if(ndtemp->info() == S_S_CON)
                        {
                            BSFC_combine_elements(k - 4, EmTemp, HT_Elem_Ptr, HT_Node_Ptr,
                                                  sfc_vert_ptr[j].destination_proc);
                            k = 8;  //exit out of the loop because we found a constrained node...
                        }
                        k++;
                    }
                    j++;
                }
            }
        }
    }
    
    return;
}

