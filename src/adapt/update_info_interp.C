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
 * $Id: update_info_interp.C 127 2007-06-07 19:48:25Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include "../header/hpfem.h"
#include "../header/refined_neighbor_info.h"

void update_neighbor_interprocessor(NodeHashTable* HT_Elem_Ptr, NodeHashTable* HT_Node_Ptr, refined_neighbor* refined_start,
                                    int myid, int numprocs)
{
#ifdef USE_MPI
    int SIDE_SONS[4][2] =
    {
    { 0, 1 },
      { 1, 2 },
      { 2, 3 },
      { 3, 0 } };
    int i;
    Element* EmTemp;
    Element* Son;
    Node* NodeTemp;
    unsigned* son_key;
    SFC_Key change_neighbor[8];
    int which_neighbor;
    int neighbor_gen;
    int my_gen;
    int neighbor_proc;
    int mate_proc;
    refined_neighbor* refined_current = refined_start->next;
    refined_neighbor_pack** send_buffer = new refined_neighbor_pack*[numprocs];
    refined_neighbor_pack** recv_buffer = new refined_neighbor_pack*[numprocs];
    for(i = 0; i < numprocs; i++)
    {
        send_buffer[i] = NULL;
        recv_buffer[i] = NULL;
    }
    int* send_counter = new int[numprocs];
    int* recv_counter = new int[numprocs];
    int refinedtag = 616; //Andrew's b.day 06.16.
    int atag = 1130; //Andrew's name day 11.30.
    
    MPI_Status status;
    MPI_Request* recv_request = new MPI_Request[numprocs];
    MPI_Request* request = new MPI_Request[numprocs];
    MPI_Request* request2 = new MPI_Request[numprocs];
    extern MPI_Datatype REFINED_INFO;
    
    for(i = 0; i < numprocs; i++)
        send_counter[i] = 0;
    
    //counting how many packs will be sent to each proc
    refined_current = refined_start->next;
    while (refined_current)
    {
        send_counter[refined_current->target_proc]++;
        refined_current = refined_current->next;
    }
    
    //creating the space for the packs
    
    int index1 = 0;
    int* index2 = new int[numprocs];
    for(i = 0; i < numprocs; i++)
    {
        send_buffer[i] = new refined_neighbor_pack[send_counter[i]];
        index2[i] = 0;
    }
    
    //filling up the packs
    refined_current = refined_start->next;
    
    while (refined_current)
    {
        index1 = refined_current->target_proc;
        
        send_buffer[index1][index2[index1]].orig_gen = refined_current->orig_gen;
        for(int j = 0; j < KEYLENGTH; j++)
        {
            send_buffer[index1][index2[index1]].target_element[j] = refined_current->target_element[j];
            send_buffer[index1][index2[index1]].old_neighbor[j] = refined_current->old_neighbor[j];
            
            for(int k = 0; k < 2; k++)
                send_buffer[index1][index2[index1]].sons[k][j] = refined_current->sons[k][j];
        }
        
        refined_current = refined_current->next;
        index2[index1]++;
        
    }
    
    for(i = 1; i < numprocs; i++)
    {
        
        mate_proc = i ^ myid; //bitwise xor
        assert(mate_proc < numprocs);
        
        MPI_Isend(&send_counter[mate_proc], 1, MPI_INT, mate_proc, atag, MPI_COMM_WORLD, (request + mate_proc));
        MPI_Isend(send_buffer[mate_proc], send_counter[mate_proc], REFINED_INFO, mate_proc, refinedtag, MPI_COMM_WORLD,
                  (request2 + mate_proc));
        
        MPI_Recv(&recv_counter[mate_proc], 1, MPI_INT, mate_proc, atag, MPI_COMM_WORLD, &status);
        recv_buffer[mate_proc] = new refined_neighbor_pack[recv_counter[mate_proc]];
        
        MPI_Irecv(recv_buffer[mate_proc], recv_counter[mate_proc], REFINED_INFO, mate_proc, refinedtag, MPI_COMM_WORLD,
                  &recv_request[mate_proc]);
        
    } //end of communication
    
    for(i = 1; i < numprocs; i++)
    {
        mate_proc = i ^ myid;
        MPI_Wait(&recv_request[mate_proc], &status);
        
        for(int k = 0; k < recv_counter[mate_proc]; k++)
        {
            EmTemp = (Element*) HT_Elem_Ptr->lookup(sfc_key_from_oldkey(recv_buffer[mate_proc][k].target_element));
            assert(EmTemp);
            which_neighbor = EmTemp->which_neighbor(sfc_key_from_oldkey(recv_buffer[mate_proc][k].old_neighbor));
            assert(which_neighbor < 4);
            neighbor_gen = recv_buffer[mate_proc][k].orig_gen;
            assert(neighbor_gen == recv_buffer[mate_proc][k].orig_gen);
            my_gen = EmTemp->generation();
            neighbor_proc = EmTemp->neigh_proc(which_neighbor);
            
            //dealing with the different cases then
            if(!EmTemp->refined_flag())
            {
                assert(my_gen >= neighbor_gen);
                if(my_gen == neighbor_gen)
                {
                    SET_NEWKEY(change_neighbor[0],recv_buffer[mate_proc][k].sons[0]);
                    SET_NEWKEY(change_neighbor[1],recv_buffer[mate_proc][k].sons[1]);
                    EmTemp->change_neighbor(change_neighbor, which_neighbor, neighbor_proc, 10);
                    NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(which_neighbor + 4));
                    assert(NodeTemp);
                    NodeTemp->info(S_C_CON); //theoretically the order doesn't have to be changed
                }
                else
                {
                    SET_NEWKEY(change_neighbor[0],recv_buffer[mate_proc][k].sons[0]);
                    SET_NEWKEY(change_neighbor[1],recv_buffer[mate_proc][k].sons[1]);
                    EmTemp->change_neighbor(change_neighbor, which_neighbor, neighbor_proc, 11);
                    NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(which_neighbor + 4));
                    assert(NodeTemp);
                    if(!(NodeTemp->info() == S_S_CON))
                    {
                        printf("which_neighbor=%d NodeTemp->info==%d\n", which_neighbor, NodeTemp->info());
                    }
                    
                    assert(NodeTemp->info() == S_S_CON);
                    NodeTemp->info(SIDE);
                    
                    NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(which_neighbor));
                    if(NodeTemp->info() == S_C_CON)
                        NodeTemp->info(CORNER);
                    else
                    {
                        NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(1 + which_neighbor));
                        if(NodeTemp->info() == S_C_CON)
                            NodeTemp->info(CORNER);
                        
                    }
                    
                }
                
            }
            else //this element was also refined
            {
                assert(my_gen == neighbor_gen);
                NodeTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->node_key(which_neighbor + 4));
                assert(NodeTemp->info() == S_C_CON);
                NodeTemp->info(CORNER);
                
                for(int l = 0; l < 2; l++)
                {
                    Son = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(SIDE_SONS[which_neighbor][l]));
                    SET_NEWKEY(change_neighbor[0],recv_buffer[mate_proc][k].sons[l]);
                    //SET_NEWKEY(change_neighbor[1],recv_buffer[mate_proc][k].sons[l+1]);
                    Son->change_neighbor(change_neighbor, which_neighbor, neighbor_proc, 11);
                    NodeTemp = (Node*) HT_Node_Ptr->lookup(Son->node_key(which_neighbor + 4));
                    assert(NodeTemp);
                    NodeTemp->info(SIDE);
                }
            }
            
        }
    }
    
    for(i = 0; i < numprocs; i++)
    {
        if(send_buffer[i])
            delete[] (send_buffer[i]);
        if(recv_buffer[i])
            delete[] (recv_buffer[i]);
    }
    delete[] send_buffer;
    delete[] recv_buffer;
    delete[] send_counter;
    delete[] recv_counter;
    delete[] recv_request;
    delete[] index2;
    for(i = 1; i < numprocs; i++)
    {
        mate_proc = i ^ myid; //bitwise xor
        MPI_Wait((request + mate_proc), &status);
        MPI_Wait((request2 + mate_proc), &status);
    }
    delete[] request;
    delete[] request2;
    return;
#endif //USE_MPI
}
