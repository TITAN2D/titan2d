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
 * $Id: data_update.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

void write_node_info(HashTable* HT_Node_Ptr, Element* Em1, Element* Em2, int side, int case_flag);

void write_node_info_ext(HashTable* HT_Node_Ptr, Element* Em, int start, int mid);

void data_update(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, Recv* RecvHead, int myid, int numprocs, int h_count)
{
    int i, j, k;
    
    Recv* recv;
    Recv* recv_new;
    
    Element* EmTemp;
    Element* SonTemp1;
    Element* SonTemp2;
    
    
    void* p;
    
    /*--scan the Recv list--*/
    recv = RecvHead->next;
    while (recv)
    {
        EmTemp = recv->targetP;
        if(recv->sender_gen == EmTemp->generation())
        {
            assert(recv->side < 4);
            if(EmTemp->refined_flag()) //-- if the receptor was refined
            {
                /*-- case 1:(refer to my note book to know the definition of each case)
                 sender and receptor were both refined and they are on the same generation-----*/
                /*-- case 3:
                 sender was not refined but receptor was and they are on the same generation---*/;
                
                switch (recv->side)
                {
                    
                    case 0: //-- filling the neigbor info for son 0 and 1
                    
                        SonTemp1 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(0)); //--son 0, side 0
                        SonTemp2 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(1)); //--son 1, side 0
                                
                        if(recv->sender_refined) //--sender was refined, case 1
                        {
                            
                            SonTemp1->set_neighbor(0, sfc_key_from_oldkey(recv->sender_son2));
                            if(recv->sender_order[1] > SonTemp1->order(0))
                                SonTemp1->set_order(0, recv->sender_order[1]);
                            SonTemp1->set_neigh_proc(0, recv->sender_id);
                            SonTemp1->get_neigh_gen(0, recv->sender_gen + 1); //--999
                                    
                            SonTemp2->set_neighbor(0, sfc_key_from_oldkey(recv->sender_son1));
                            if(recv->sender_order[0] > SonTemp2->order(0))
                                SonTemp2->set_order(0, recv->sender_order[0]);
                            SonTemp2->set_neigh_proc(0, recv->sender_id);
                            SonTemp2->get_neigh_gen(0, recv->sender_gen + 1); //--999
                                    
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 0, 1);
                            
                        }
                        else //--sender was not refined, case 3
                        {
                            
                            SonTemp1->set_neighbor(0, EmTemp->neighbor(0));
                            SonTemp1->set_neigh_proc(0, recv->sender_id);
                            SonTemp1->get_neigh_gen(0, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp1->order(0))
                                SonTemp1->set_order(0, recv->sender_order[0]);
                            
                            SonTemp2->set_neigh_proc(0, recv->sender_id);
                            SonTemp2->set_neighbor(0, sfc_key_from_oldkey(recv->sender));
                            SonTemp2->get_neigh_gen(0, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(0))
                                SonTemp2->set_order(0, recv->sender_order[0]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 0, 3);
                            
                        }
                        break;
                        
                    case 1: //-- son 1 and 2
                    
                        SonTemp1 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(1)); //--son 1, side 1
                        SonTemp2 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(2)); //--son 2, side 1
                                
                        if(recv->sender_refined)
                        {
                            
                            SonTemp1->set_neighbor(1, sfc_key_from_oldkey(recv->sender_son2));
                            SonTemp1->set_neigh_proc(1, recv->sender_id);
                            SonTemp1->get_neigh_gen(1, recv->sender_gen + 1); //--999
                            if(recv->sender_order[1] > SonTemp1->order(1))
                                SonTemp1->set_order(1, recv->sender_order[1]);
                            
                            SonTemp2->set_neighbor(1, sfc_key_from_oldkey(recv->sender_son1));
                            SonTemp2->set_neigh_proc(1, recv->sender_id);
                            SonTemp2->get_neigh_gen(1, recv->sender_gen + 1); //--999
                            if(recv->sender_order[0] > SonTemp2->order(1))
                                SonTemp2->set_order(1, recv->sender_order[0]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 1, 1);
                            
                        }
                        else
                        {
                            
                            SonTemp1->set_neighbor(1, sfc_key_from_oldkey(recv->sender));
                            SonTemp1->set_neigh_proc(1, recv->sender_id);
                            SonTemp1->get_neigh_gen(1, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(1))
                                SonTemp1->set_order(1, recv->sender_order[0]);
                            
                            SonTemp2->set_neighbor(1, sfc_key_from_oldkey(recv->sender));
                            SonTemp2->set_neigh_proc(1, recv->sender_id);
                            SonTemp2->get_neigh_gen(1, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(1))
                                SonTemp1->set_order(1, recv->sender_order[0]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 1, 3);
                            
                        }
                        break;
                        
                    case 2: //-- son 2 and 3
                    
                        SonTemp1 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(2)); //--son 2, side 2
                        SonTemp2 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(3)); //--son 3, side 2
                                
                        if(recv->sender_refined)
                        {
                            
                            SonTemp1->set_neighbor(2, sfc_key_from_oldkey(recv->sender_son2));
                            SonTemp1->set_neigh_proc(2, recv->sender_id);
                            SonTemp1->get_neigh_gen(2, recv->sender_gen + 1); //--999
                            if(recv->sender_order[1] > SonTemp1->order(2))
                                SonTemp1->set_order(2, recv->sender_order[1]);
                            
                            SonTemp2->set_neighbor(2, sfc_key_from_oldkey(recv->sender_son1));
                            SonTemp2->set_neigh_proc(2, recv->sender_id);
                            SonTemp2->get_neigh_gen(2, recv->sender_gen + 1); //--999
                            if(recv->sender_order[0] > SonTemp2->order(2))
                                SonTemp2->set_order(2, recv->sender_order[1]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 2, 1);
                            
                        }
                        else
                        {
                            
                            SonTemp1->set_neighbor(2, sfc_key_from_oldkey(recv->sender));
                            SonTemp1->set_neigh_proc(2, recv->sender_id);
                            SonTemp1->get_neigh_gen(2, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(2))
                                SonTemp1->set_order(2, recv->sender_order[0]);
                            
                            SonTemp2->set_neighbor(2, sfc_key_from_oldkey(recv->sender));
                            SonTemp2->set_neigh_proc(2, recv->sender_id);
                            SonTemp2->get_neigh_gen(2, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(2))
                                SonTemp1->set_order(2, recv->sender_order[0]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 2, 3);
                            
                        }
                        break;
                        
                    case 3: //-- son 3 and 0
                    
                        SonTemp1 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(3)); //--son 3, side 3
                        SonTemp2 = (Element*) HT_Elem_Ptr->lookup(EmTemp->son(0)); //--son 0, side 3
                                
                        if(recv->sender_refined)
                        {
                            
                            SonTemp1->set_neighbor(3, sfc_key_from_oldkey(recv->sender_son2));
                            SonTemp1->set_neigh_proc(3, recv->sender_id);
                            SonTemp1->get_neigh_gen(3, recv->sender_gen + 1); //--999
                            if(recv->sender_order[1] > SonTemp1->order(3))
                                SonTemp1->set_order(3, recv->sender_order[1]);
                            
                            SonTemp2->set_neighbor(3, sfc_key_from_oldkey(recv->sender_son1));
                            SonTemp2->set_neigh_proc(3, recv->sender_id);
                            SonTemp2->get_neigh_gen(3, recv->sender_gen + 1); //--999
                            if(recv->sender_order[0] > SonTemp2->order(3))
                                SonTemp2->set_order(3, recv->sender_order[1]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 3, 1);
                            
                        }
                        else
                        {
                            
                            SonTemp1->set_neighbor(3, sfc_key_from_oldkey(recv->sender));
                            SonTemp1->set_neigh_proc(3, recv->sender_id);
                            SonTemp1->get_neigh_gen(3, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(3))
                                SonTemp1->set_order(3, recv->sender_order[0]);
                            
                            SonTemp2->set_neighbor(3, sfc_key_from_oldkey(recv->sender));
                            SonTemp2->set_neigh_proc(3, recv->sender_id);
                            SonTemp2->get_neigh_gen(3, recv->sender_gen); //--999
                            if(recv->sender_order[0] > SonTemp2->order(3))
                                SonTemp1->set_order(3, recv->sender_order[0]);
                            
                            write_node_info(HT_Node_Ptr, SonTemp1, SonTemp2, 3, 3);
                            
                        }
                        break;
                        
                }
                
            }
            else
            {
                
                /*--case 2:
                 sender was refined but receptor was not and they are on the same generation--*/

                if(recv->sender_refined)
                {
                    
                    switch (recv->side % 4)
                    //--refer to the definition to understand the following
                    {
                        
                        case 0:

                            EmTemp->set_neigh_proc(0, recv->sender_id);
                            EmTemp->set_neigh_proc(4, recv->sender_id);
                            EmTemp->set_neighbor(0, sfc_key_from_oldkey(recv->sender_son2));
                            EmTemp->set_neighbor(4, sfc_key_from_oldkey(recv->sender_son1));
                            EmTemp->get_neigh_gen(0, recv->sender_gen + 1); //--999
                            EmTemp->get_neigh_gen(4, recv->sender_gen + 1); //--999
                                    
                            if(recv->sender_order[0] > EmTemp->order(0))
                                EmTemp->set_order(0, recv->sender_order[0]);
                            if(recv->sender_order[1] > EmTemp->order(0))
                                EmTemp->set_order(0, recv->sender_order[1]);
                            write_node_info(HT_Node_Ptr, EmTemp, EmTemp, 0, 2);
                            break;
                            
                        case 1:
                            EmTemp->set_neigh_proc(1, recv->sender_id);
                            EmTemp->set_neigh_proc(5, recv->sender_id);
                            EmTemp->set_neighbor(1, sfc_key_from_oldkey(recv->sender_son2));
                            EmTemp->set_neighbor(5, sfc_key_from_oldkey(recv->sender_son1));
                            EmTemp->get_neigh_gen(1, recv->sender_gen + 1); //--999
                            EmTemp->get_neigh_gen(5, recv->sender_gen + 1); //--999
                                    
                            if(recv->sender_order[0] > EmTemp->order(0))
                                EmTemp->set_order(1, recv->sender_order[0]);
                            if(recv->sender_order[1] > EmTemp->order(0))
                                EmTemp->set_order(1, recv->sender_order[1]);
                            write_node_info(HT_Node_Ptr, EmTemp, EmTemp, 1, 2);
                            break;
                            
                        case 2:
                            EmTemp->set_neigh_proc(2, recv->sender_id);
                            EmTemp->set_neigh_proc(6, recv->sender_id);
                            EmTemp->set_neighbor(2, sfc_key_from_oldkey(recv->sender_son2));
                            EmTemp->set_neighbor(6, sfc_key_from_oldkey(recv->sender_son1));
                            EmTemp->get_neigh_gen(2, recv->sender_gen + 1); //--999
                            EmTemp->get_neigh_gen(6, recv->sender_gen + 1); //--999
                                    
                            if(recv->sender_order[0] > EmTemp->order(0))
                                EmTemp->set_order(2, recv->sender_order[0]);
                            if(recv->sender_order[1] > EmTemp->order(0))
                                EmTemp->set_order(2, recv->sender_order[1]);
                            write_node_info(HT_Node_Ptr, EmTemp, EmTemp, 2, 2);
                            break;
                            
                        case 3:
                            EmTemp->set_neigh_proc(3, recv->sender_id);
                            EmTemp->set_neigh_proc(7, recv->sender_id);
                            EmTemp->set_neighbor(3, sfc_key_from_oldkey(recv->sender_son2));
                            EmTemp->set_neighbor(7, sfc_key_from_oldkey(recv->sender_son1));
                            EmTemp->get_neigh_gen(3, recv->sender_gen + 1); //--999
                            EmTemp->get_neigh_gen(7, recv->sender_gen + 1); //--999
                                    
                            if(recv->sender_order[0] > EmTemp->order(0))
                                EmTemp->set_order(3, recv->sender_order[0]);
                            if(recv->sender_order[1] > EmTemp->order(0))
                                EmTemp->set_order(3, recv->sender_order[1]);
                            write_node_info(HT_Node_Ptr, EmTemp, EmTemp, 3, 2);
                            break;
                            
                    }
                }
                
                else /*--both sender and receptor are not refined and they are on the same generation--*/
                {
                    if(recv->sender_order[0] > EmTemp->order(recv->side))
                        EmTemp->set_order(recv->side, recv->sender_order[0]);
                    p = HT_Node_Ptr->lookup(EmTemp->node_key(recv->side + 4)); //-- cc
                    assert(p);
                    Node* NdTemp = (Node*) p;
                    NdTemp->order(EmTemp->order(recv->side));
                }
            }
            
        }
        
        else //-- sender's generation is diffrent from receptor's
        {
            if(recv->sender_gen > EmTemp->generation()) //--sender is smaller than receptor
            {
                assert(!recv->sender_refined); //-- debug window
                if(EmTemp->refined_flag()) //-- happened only in h-refinement
                {
                    /*-- search, which son of EmTemp? a is the result --*/
                    int a = recv->side;
                    if(a == 7)
                        a = 0;
                    else if(a > 3)
                        a = a - 3;
                    
                    SonTemp1 = (Element*) (HT_Elem_Ptr->lookup(EmTemp->son(a)));
                    
                    int start;
                    int end;
                    int mid;
                    
                    /*-- find out which side of this son be modified--*/
                    if(a == recv->side)
                    {
                        start = recv->side;
                        end = start + 1;
                        if(end == 4)
                            end = 0;
                        mid = start + 4;
                        write_node_info_ext(HT_Node_Ptr, SonTemp1, end, mid);
                    }
                    else
                    {
                        start = recv->side - 4;
                        mid = start + 4;
                        write_node_info_ext(HT_Node_Ptr, SonTemp1, start, mid);
                    }
                    
                    SonTemp1->set_neighbor(start, sfc_key_from_oldkey(recv->sender));
                    SonTemp1->get_neigh_gen(start, recv->sender_gen);
                    SonTemp1->set_neigh_proc(start, recv->sender_id); //--maybe surplus
                                             
                }
                else //-- p enrichment
                {
                    assert(recv->sender_gen == EmTemp->generation() + 1); //-- debug window
                            
                    int start = recv->side % 4;
                    
                    Node* NdTemp = (Node*) (HT_Node_Ptr->lookup(EmTemp->node_key(start + 4)));
                    
                    if(recv->sender_order[0] > EmTemp->order(start))
                        EmTemp->set_order(start, recv->sender_order[0]);
                    
                    int final_order = EmTemp->order(start);
                    
                    NdTemp->order(final_order);
                }
                
            }
            
            if(recv->sender_gen < EmTemp->generation()) //--sender is bigger than receptor
            {
                if(recv->sender_refined)
                {
                    assert(recv->side < 4);
                    int start = recv->side;
                    int end = start + 1;
                    if(end == 4)
                        end = 0;
                    
                    Node* NdTemp = (Node*) (HT_Node_Ptr->lookup(EmTemp->node_key(start)));
                    
                    if(NdTemp->getinfo() == S_S_CON)
                    {
                        EmTemp->set_neighbor(start, sfc_key_from_oldkey(recv->sender_son1));
                        write_node_info_ext(HT_Node_Ptr, EmTemp, start, start + 4);
                        
                    }
                    else
                    {
                        int mid = -1;
                        NdTemp = (Node*) (HT_Node_Ptr->lookup(EmTemp->node_key(end)));
                        EmTemp->set_neighbor(start, sfc_key_from_oldkey(recv->sender_son2));
                        if(end == 0)
                            mid = 7;
                        else
                            mid = end + 3;
                        write_node_info_ext(HT_Node_Ptr, EmTemp, end, mid);
                        
                    }
                    
                    EmTemp->get_neigh_gen(start, recv->sender_gen + 1);
                    EmTemp->set_neigh_proc(start, recv->sender_id);
                    
                }
                else
                {
                    assert(recv->sender_gen == EmTemp->generation() - 1); //-- debug window
                            
                    int start = recv->side;
                    int end = start + 1;
                    if(end == 4)
                        end = 0;
                    
                    Node* NdTemp = (Node*) (HT_Node_Ptr->lookup(EmTemp->node_key(start)));
                    if(NdTemp->getinfo() != S_C_CON)
                    {
                        NdTemp = (Node*) (HT_Node_Ptr->lookup(EmTemp->node_key(end)));
                        //assert(NdTemp->getinfo()==S_C_CON);
                    }
                    
                    if(recv->sender_order[0] > EmTemp->order(start))
                        EmTemp->set_order(start, recv->sender_order[0]);
                    
                    int final_order = EmTemp->order(start);
                    
                    NdTemp->order(final_order);
                }
                
            }
            
        }
        
        recv = recv->next;
    }
}

void write_node_info(HashTable* HT_Node_Ptr, Element* Em1, Element* Em2, int side, int case_flag)
{
    void* p;
    Node* NdTemp1;
    Node* NdTemp2;
    Node* NdTemp3;
    int order;
    
    if(Em1 == Em2) //--case 2
    {
        p = HT_Node_Ptr->lookup(Em1->node_key(side + 4));
        NdTemp1 = (Node*) p;
        NdTemp1->putinfo(S_C_CON);
    }
    else
    {
        switch (side)
        {
            case 0: //--side 0 of son 0, 1
                p = HT_Node_Ptr->lookup(Em1->node_key(0));
                NdTemp1 = (Node*) p; //-- node 0 of son 0
                p = HT_Node_Ptr->lookup(Em1->node_key(1));
                NdTemp2 = (Node*) p; //-- node 1 of son 0
                p = HT_Node_Ptr->lookup(Em1->node_key(4));
                NdTemp3 = (Node*) p; //-- node 4 of son 0
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER); //-- node 0 of son 0 in case 1 or case 3 is a corner
                    NdTemp2->putinfo(CORNER); //-- node 1 of son 0 in case 1 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 4 of son 0 in case 1 is a side
                    NdTemp3->order(Em1->order(0));
                }
                else if(case_flag == 3)
                {
                    NdTemp1->putinfo(CORNER);
                    NdTemp2->putinfo(S_C_CON); //-- node 1 of son 0 in case 3 is a constrained node, has side dof
                    order = Em1->order(0) > Em2->order(0) ? Em1->order(0) : Em2->order(0);
                    NdTemp2->order(order);
                    NdTemp3->putinfo(S_S_CON);  //-- node 4 of son 0 in case 3 has no dof
                }
                
                p = HT_Node_Ptr->lookup(Em2->node_key(0));
                NdTemp1 = (Node*) p;  //-- node 0 of son 1
                p = HT_Node_Ptr->lookup(Em2->node_key(1));
                NdTemp2 = (Node*) p;  //-- node 1 of son 1
                p = HT_Node_Ptr->lookup(Em2->node_key(4));
                NdTemp3 = (Node*) p;  //-- node 4 of son 1
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 0 of son 1 in case 1 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 1 of son 1 in case 1 or case 3 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 4 of son 1 in case 1 is a side
                    NdTemp3->order(Em2->order(0));
                }
                else if(case_flag == 3)
                {
                    NdTemp2->putinfo(CORNER); //-- node 0 of son 1 in case 1 is a corner
                    //NdTemp1->putinfo(S_C_CON);//-- node 0 of son 1 in case 3 is a constrained node, has side dof, same node as node 1 of son 0
                    NdTemp3->putinfo(S_S_CON);	      //-- node 4 of son 0 in case 3 has no dof
                }
                
                break;
                
            case 1:	      //--side 1 of son 1, 2
                p = HT_Node_Ptr->lookup(Em1->node_key(1));
                NdTemp1 = (Node*) p;	      //-- node 1 of son 1
                p = HT_Node_Ptr->lookup(Em1->node_key(2));
                NdTemp2 = (Node*) p;	      //-- node 2 of son 1
                p = HT_Node_Ptr->lookup(Em1->node_key(5));
                NdTemp3 = (Node*) p;	      //-- node 5 of son 1
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);	      //-- node 1 of son 1 in case 1 or case 3 is a corner
                    NdTemp2->putinfo(CORNER);	      //-- node 2 of son 1 in case 1 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 5 of son 1 in case 1 is a side node
                    NdTemp3->order(Em1->order(1));
                }
                else if(case_flag == 3)
                {
                    NdTemp1->putinfo(CORNER);
                    NdTemp2->putinfo(S_C_CON); //-- node 2 of son 1 in case 3 is a constrained node, has side dof
                    order = Em1->order(1) > Em2->order(1) ?
                            Em1->order(1) : Em2->order(1);
                    NdTemp2->order(order);
                    NdTemp3->putinfo(S_S_CON);  //-- node 5 of son 1 in case 3 has no dof
                }
                
                p = HT_Node_Ptr->lookup(Em2->node_key(1));
                NdTemp1 = (Node*) p;  //-- node 1 of son 2
                p = HT_Node_Ptr->lookup(Em2->node_key(2));
                NdTemp2 = (Node*) p;  //-- node 2 of son 2
                p = HT_Node_Ptr->lookup(Em2->node_key(5));
                NdTemp3 = (Node*) p;  //-- node 5 of son 2
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 1 of son 2 in case 1 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 2 of son 2 in case 1 or case 3 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 5 of son 2 in case 1 is a side
                    NdTemp3->order(Em2->order(1));
                }
                else if(case_flag == 3)
                {
                    NdTemp2->putinfo(CORNER);
                    //NdTemp1->putinfo(S_C_CON);//-- node 1 of son 2 in case 3 is a constrained node, has side dof
                    NdTemp3->putinfo(S_S_CON);  //-- node 5 of son 2 in case 3 has no dof
                }
                
                break;
                
            case 2:  //--side 2 of son 2, 3
                p = HT_Node_Ptr->lookup(Em1->node_key(2));
                NdTemp1 = (Node*) p;  //-- node 2 of son 2
                p = HT_Node_Ptr->lookup(Em1->node_key(3));
                NdTemp2 = (Node*) p;  //-- node 3 of son 2
                p = HT_Node_Ptr->lookup(Em1->node_key(6));
                NdTemp3 = (Node*) p;  //-- node 6 of son 2
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 2 of son 2 in case 1 or case 3 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 3 of son 2 in case 1 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 6 of son 2 in case 1 is a side node
                    NdTemp3->order(Em1->order(2));
                }
                else if(case_flag == 3)
                {
                    NdTemp1->putinfo(CORNER);
                    NdTemp2->putinfo(S_C_CON); //-- node 3 of son 2 in case 3 is a constrained node, has side dof
                    order = Em1->order(2) > Em2->order(2) ?
                            Em1->order(2) : Em2->order(2);
                    NdTemp2->order(order);
                    NdTemp3->putinfo(S_S_CON);  //-- node 6 of son 2 in case 3 has no dof
                }
                
                p = HT_Node_Ptr->lookup(Em2->node_key(2));
                NdTemp1 = (Node*) p;  //-- node 2 of son 3
                p = HT_Node_Ptr->lookup(Em2->node_key(3));
                NdTemp2 = (Node*) p;  //-- node 3 of son 3
                p = HT_Node_Ptr->lookup(Em2->node_key(6));
                NdTemp3 = (Node*) p;  //-- node 6 of son 3
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 2 of son 3 in case 1 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 3 of son 3 in case 1 or case 3 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 6 of son 3 in case 1 is a side
                    NdTemp3->order(Em2->order(2));
                }
                else if(case_flag == 3)
                {
                    NdTemp2->putinfo(CORNER);
                    //NdTemp1->putinfo(S_C_CON);//-- node 2 of son 3 in case 3 is a constrained node, has side dof
                    NdTemp3->putinfo(S_S_CON);  //-- node 6 of son 3 in case 3 has no dof
                }
                
                break;
                
            case 3:
                p = HT_Node_Ptr->lookup(Em1->node_key(3));
                NdTemp1 = (Node*) p;  //-- node 3 of son 3
                p = HT_Node_Ptr->lookup(Em1->node_key(0));
                NdTemp2 = (Node*) p;  //-- node 0 of son 3
                p = HT_Node_Ptr->lookup(Em1->node_key(7));
                NdTemp3 = (Node*) p;  //-- node 7 of son 3
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 3 of son 3 in case 1 or case 3 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 0 of son 3 in case 1 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 7 of son 3 in case 1 is a side node
                    NdTemp3->order(Em1->order(3));
                }
                else if(case_flag == 3)
                {
                    NdTemp1->putinfo(CORNER);
                    NdTemp2->putinfo(S_C_CON); //-- node 0 of son 3 in case 3 is a constrained node, has side dof
                    order = Em1->order(3) > Em2->order(3) ?
                            Em1->order(3) : Em2->order(3);
                    NdTemp2->order(order);
                    NdTemp3->putinfo(S_S_CON);  //-- node 7 of son 3 in case 3 has no dof
                }
                
                p = HT_Node_Ptr->lookup(Em2->node_key(3));
                NdTemp1 = (Node*) p;  //-- node 3 of son 0
                p = HT_Node_Ptr->lookup(Em2->node_key(0));
                NdTemp2 = (Node*) p;  //-- node 0 of son 0
                p = HT_Node_Ptr->lookup(Em2->node_key(7));
                NdTemp3 = (Node*) p;  //-- node 7 of son 0
                if(case_flag == 1)
                {
                    NdTemp1->putinfo(CORNER);  //-- node 3 of son 0 in case 1 is a corner
                    NdTemp2->putinfo(CORNER);  //-- node 0 of son 0 in case 1 or case 3 is a corner
                    NdTemp3->putinfo(SIDE);  //-- node 7 of son 0 in case 1 is a side
                    NdTemp3->order(Em2->order(3));
                }
                else if(case_flag == 3)
                {
                    NdTemp2->putinfo(CORNER);
                    //NdTemp1->putinfo(S_C_CON);//-- node 3 of son 0 in case 3 is a constrained node, has side dof
                    NdTemp3->putinfo(S_S_CON);  //-- node 7 of son 0 in case 3 has no dof
                }
                
                break;
        }
    }
}

void write_node_info_ext(HashTable* HT_Node_Ptr, Element* Em, int start, int mid)
{
    
    Node* NdTemp;
    
    int order = Em->order(mid - 4);
    
    NdTemp = (Node*) (HT_Node_Ptr->lookup(Em->node_key(start)));
    assert(NdTemp);
    NdTemp->putinfo(CORNER);
    NdTemp->order(1);
    
    NdTemp = (Node*) (HT_Node_Ptr->lookup(Em->node_key(mid)));
    assert(NdTemp);
    NdTemp->putinfo(SIDE);
    NdTemp->order(order);
    
}
