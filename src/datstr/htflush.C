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
 * $Id: htflush.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"

void htflush(ElementsHashTable* ht_elem_ptr, NodeHashTable* ht_node_ptr, int option)
{
    
    int i, j, k;
    Element* EmTemp;
    Node* NdTemp;
    unsigned KeyTemp[KEYLENGTH];
    unsigned* keyP;
    void* p;
    int* dofP;
    double* sol;
    int e_buckets = ht_elem_ptr->get_no_of_buckets();
    //int n_buckets = ht_node_ptr->get_no_of_buckets();
    
    switch (option)
    {
        case 1:
            for(i = 0; i < ht_elem_ptr->elenode_.size(); i++)
            {
                if(ht_elem_ptr->status_[i]>=0)
                {
                    ht_elem_ptr->elenode_[i].set_new_old(OLD);
                }
            }
            break;
        /*case 2:
            for(i = 0; i < n_buckets; i++)
            {
                entryp = *(ht_node_ptr->getbucketptr() + i);
                while (entryp)
                {
                    NdTemp = (Node*) (entryp->value);
                    //NdTemp->putinfo(INIT);
                                        
                    entryp = entryp->next;
                }
            }*/
            /*for(i = 0; i < e_buckets; i++)
            {
                entryp = *(ht_elem_ptr->getbucketptr() + i);
                while (entryp)
                {
                    EmTemp = (Element*) (entryp->value);
                    for(j = 0; j < 8; j++)
                    {
                        EmTemp->put_recv_flag(j, 0);
                        EmTemp->put_send_flag(j, 0);
                    }
                    entryp = entryp->next;
                }
            }*/
    }
}
