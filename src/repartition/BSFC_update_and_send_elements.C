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
 * $Id: BSFC_update_and_send_elements.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/exvar.h"
#include "./repartition_BSFC.h"

void create_element(ElemPack* elem2, ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, int myid);

void BSFC_update_and_send_elements(int myid, int numprocs, ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr,
                                   int time_step)
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
            if(EmTemp->refined_flag())
            {
                EmTemp->set_new_old(-1);
            }
        }
    }
    // now figure out what processors need to have neighbor info updated 
    // and what processors to send elements to
    int* send_info = new int[numprocs * 2]; // number of {neigh_info, elements} for each proc
    for(i = 0; i < 2 * numprocs; i++)
        send_info[i] = 0;
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
            {
                if(myid != EmTemp->myprocess())
                { // this element will get moved to a new processor
                  // neigh info
                    for(j = 0; j < 8; j++)
                        if(EmTemp->neigh_proc(j) != myid && EmTemp->neigh_proc(j) >= 0)
                            send_info[EmTemp->neigh_proc(j) * 2] += 1;
                    
                    // element info
                    send_info[EmTemp->myprocess() * 2 + 1] += 1;
                    
                }
            }
        }
    }
    
    int* temp_info = new int[numprocs * 2];
    for(i = 0; i < 2 * numprocs; i++)
        temp_info[i] = send_info[i];
    int* recv_info = new int[numprocs * 2];
    i = MPI_Alltoall(temp_info, 2, MPI_INT, recv_info, 2, MPI_INT, MPI_COMM_WORLD);
    delete[] temp_info;
    //debug stuff
    /*  for(i=0;i<numprocs;i++)
     if(myid != i) 
     printf("proc %d send %d neighs and %d elms and receive %d neighs and %d elms from proc %d\n",
     myid, send_info[i*2], send_info[i*2+1], recv_info[i*2], recv_info[i*2+1], i); */
    //end debug stuff
    //first post all of the receives...
    MPI_Request* recv_request = new MPI_Request[2 * numprocs];
    int recv_count[2] =
    { 0, 0 };
    for(i = 0; i < numprocs; i++)
    {
        recv_count[0] += recv_info[2 * i];
        recv_count[1] += recv_info[2 * i + 1];
    }
    int counter_recv[2] =
    { 0, 0 };
    unsigned* recv_neigh_array = new unsigned[recv_count[0] * (2 * KEYLENGTH + 1)]; // 2 keys + 1 proc
    ElemPack* recv_elm_array = new ElemPack[recv_count[1]];
    int neigh_tag = 99565, elm_tag = 88476;  // random tag numbers
            
    for(i = 0; i < numprocs; i++)
    {
        if(recv_info[2 * i] != 0)
        {  // receive neighbor info
            j = MPI_Irecv((recv_neigh_array + counter_recv[0]), recv_info[2 * i] * (2 * KEYLENGTH + 1), MPI_UNSIGNED, i,
                          neigh_tag, MPI_COMM_WORLD, (recv_request + i * 2));
            counter_recv[0] += recv_info[2 * i] * (2 * KEYLENGTH + 1);
        }
        if(recv_info[2 * i + 1] != 0)
        { // receive elements
            j = MPI_Irecv((recv_elm_array + counter_recv[1]), recv_info[2 * i + 1], ELEMTYPE, i, elm_tag,
                          MPI_COMM_WORLD, (recv_request + i * 2 + 1));
            counter_recv[1] += recv_info[2 * i + 1];
        }
    }
    // done posting the receives
    
    int send_count = 0;
    for(i = 0; i < numprocs; i++)
        send_count += send_info[2 * i];
    
    unsigned* send_neigh_array = new unsigned[send_count * (2 * KEYLENGTH + 1)];
    for(i = 0; i < send_count * (2 * KEYLENGTH + 1); i++)
        send_neigh_array[i] = 5;
    int* counter_send_proc = new int[numprocs];
    counter_send_proc[0] = 0;
    for(i = 1; i < numprocs; i++)
        counter_send_proc[i] = counter_send_proc[i - 1] + send_info[2 * (i - 1)];
    // pack neighbor information to send out
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
            {
                if(myid != EmTemp->myprocess())
                { // this element will get moved to a new processor
                    for(j = 0; j < 8; j++)
                        if(EmTemp->neigh_proc(j) != myid && EmTemp->neigh_proc(j) > -1)
                        {
                            //printf("update neighbor is %u %u on proc %d to proc %d \n", *(EmTemp->get_neighbors()+j*KEYLENGTH+0), *(EmTemp->get_neighbors()+j*KEYLENGTH+1), myid, neigh_proc(j));

                            SET_OLDKEY((&(send_neigh_array[counter_send_proc[EmTemp->neigh_proc(j)] * (2 * KEYLENGTH + 1)])),
                                    EmTemp->neighbor(j));
                            SET_OLDKEY((&(send_neigh_array[counter_send_proc[EmTemp->neigh_proc(j)] * (2 * KEYLENGTH + 1) + 2 ])),
                                    EmTemp->key());

                            send_neigh_array[counter_send_proc[EmTemp->neigh_proc(j)] * (2 * KEYLENGTH + 1) + 4] =
                                    (unsigned) EmTemp->myprocess();
                            
                            counter_send_proc[EmTemp->neigh_proc(j)] += 1;
                        }
                }
            }
        }
    }
    // send out neighbor information
    int counter = 0;
    MPI_Request* send_request = new MPI_Request[numprocs];
    for(i = 0; i < numprocs; i++)
        if(send_info[2 * i] != 0)
        {
            j = MPI_Isend((send_neigh_array + counter * (2 * KEYLENGTH + 1)), send_info[2 * i] * (2 * KEYLENGTH + 1),
                          MPI_UNSIGNED, i, neigh_tag, MPI_COMM_WORLD, (send_request + i));
            counter += send_info[2 * i];
        }
    // update neighbor information on this processor
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
                if(myid != EmTemp->myprocess())
                { // this element will get moved to a new processor
                  // neigh info
                    for(j = 0; j < 8; j++)
                    {
                        if(EmTemp->neigh_proc(j) == myid)
                        {
                            Element* EmNeigh = (Element*) HT_Elem_Ptr->lookup(EmTemp->neighbor(j));
                            k = EmNeigh->which_neighbor(EmTemp->key());
                            EmNeigh->set_neigh_proc(k, EmTemp->myprocess());
                            EmTemp->set_neigh_proc(j, EmNeigh->myprocess());
                        }
                    }
                }
        }
    }
    
// wait to receive new neighbor info and then update neighbor info
    // from neighbor information received from other processors
    MPI_Status status;
    counter = 0;
    for(i = 0; i < numprocs; i++)
        if(recv_info[2 * i] != 0)
        {
            j = MPI_Wait((recv_request + 2 * i), &status);
            for(j = counter; j < counter + recv_info[2 * i]; j++)
            {
                EmTemp = (Element*) HT_Elem_Ptr->lookup(sfc_key_from_oldkey(recv_neigh_array + j * (2 * KEYLENGTH + 1)));
                k = EmTemp->which_neighbor(sfc_key_from_oldkey(recv_neigh_array + j * (2 * KEYLENGTH + 1) + KEYLENGTH));
                EmTemp->set_neigh_proc(k, recv_neigh_array[j * (2 * KEYLENGTH + 1) + 2 * KEYLENGTH]);
            }
            counter += recv_info[2 * i];
        }
    delete[] recv_neigh_array;
    // wait and then delete sent info that was already received
    for(i = 0; i < numprocs; i++)
        if(send_info[2 * i] != 0)
            j = MPI_Wait((send_request + i), &status);
    
    delete[] send_neigh_array;
    
    ///////////////////////////////////////////////////////
    // can now migrate any elements which need to be moved
    ///////////////////////////////////////////////////////
    
    // first pack elements in array
    send_count = 0;
    for(i = 0; i < numprocs; i++)
        send_count += send_info[2 * i + 1];
    ElemPack* send_elm_array = new ElemPack[send_count];
    
    counter_send_proc[0] = 0;
    for(i = 1; i < numprocs; i++)
        counter_send_proc[i] = counter_send_proc[i - 1] + send_info[2 * i - 1];
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(!EmTemp->refined_flag())
            {
                if(myid != EmTemp->myprocess())
                { // this element will get moved to a new processor
                    int myprocess = EmTemp->myprocess();
                    assert(myprocess >= 0 && myprocess < numprocs);
                    Pack_element(EmTemp, (send_elm_array + counter_send_proc[myprocess]), HT_Node_Ptr, myprocess);
                    counter_send_proc[myprocess] += 1;
                    HT_Elem_Ptr->removeElement(EmTemp);
                    ielm--;
                }
            }
        }
    }
    delete[] counter_send_proc;
    
    // now send out packed elements
    counter = 0;
    for(i = 0; i < numprocs; i++)
        if(send_info[2 * i + 1] != 0)
        {
            j = MPI_Isend((send_elm_array + counter), send_info[2 * i + 1], ELEMTYPE, i, elm_tag, MPI_COMM_WORLD,
                          (send_request + i));
            counter += send_info[2 * i + 1];
        }
    
    counter = 0;
    for(i = 0; i < numprocs; i++)
        if(recv_info[2 * i + 1] != 0)
        {
            j = MPI_Wait((recv_request + i * 2 + 1), &status);
            for(j = 0; j < recv_info[2 * i + 1]; j++)
            {
                create_element((recv_elm_array + counter), HT_Elem_Ptr, HT_Node_Ptr, myid);
                counter++;
            }
        }
    delete[] recv_elm_array;
    delete[] recv_request;
    delete[] recv_info;
    for(i = 0; i < numprocs; i++)
        if(send_info[2 * i + 1] != 0)
            j = MPI_Wait((send_request + i), &status);
    delete[] send_request;
    delete[] send_info;
    delete[] send_elm_array;
    
    return;
}

void delete_unused_elements_nodes(ElementsHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, int myid)
{
    int i, j;
    Element* EmTemp, *son;
    int no_of_node_buckets = HT_Node_Ptr->get_no_of_buckets();
    Node* NdTemp;
    
    int no_of_buckets = HT_Elem_Ptr->get_no_of_buckets();
    vector<HashEntryLine> &bucket=HT_Elem_Ptr->bucket;
    tivector<Element> &elenode_=HT_Elem_Ptr->elenode_;
    
    // initialize the node flags
    //@NodesSingleLoop
    for(i = 0; i < HT_Node_Ptr->elenode_.size(); i++)
    {
        if(HT_Node_Ptr->status_[i]>=0)
            HT_Node_Ptr->elenode_[i].id(0);
    }
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->refined_flag() != 0)
            {  // not an active element
                EmTemp->void_bcptr();  // don't remove bc's
                HT_Elem_Ptr->removeElement(EmTemp);
                ielm--;
            }
            else
            {  //active element on this processor -- flag all nodes as being used
                for(j = 0; j < 8; j++)
                {
                    NdTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(j));
                    NdTemp->id(1);
                }
                NdTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->key());
                NdTemp->id(1);
            }
        }
    }
    
    // delete the nodes that aren't used
    for(int i = 0; i < HT_Node_Ptr->elenode_.size(); i++)
    {
        if(HT_Node_Ptr->status_[i]>=0 && HT_Node_Ptr->elenode_[i].id() == 0)
        {
            HT_Node_Ptr->remove(HT_Node_Ptr->elenode_[i].key());
        }
    }
    return;
}
