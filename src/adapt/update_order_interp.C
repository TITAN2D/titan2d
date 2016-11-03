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
 * $Id: update_order_interp.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#include "../header/enriched_neighbor_info.h"

void update_order_interp(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			 enriched_neighbor* enriched_start, int myid, int numprocs)
{
  int i;
  Element* EmTemp;
  Element* EmTemp2;
  Node*    NodeTemp;
  int which_neighbor;
  int neighbor_gen;
  int my_gen;
  int my_order;
  int final_order;
  int neigh_order;
  int neighbor_proc;
  int mate_proc;
  enriched_neighbor* enriched_current=enriched_start->next;
  enriched_neighbor_pack** send_buffer=new enriched_neighbor_pack*[numprocs];
  enriched_neighbor_pack** recv_buffer=new enriched_neighbor_pack*[numprocs];
  for(i=0;i<numprocs;i++) {
    send_buffer[i] = NULL;
    recv_buffer[i] = NULL;
  }  int* send_counter=new int[numprocs];
  int* recv_counter=new int[numprocs];
  int enrichedtag=75;//Andrew's b.year
  int atag=24;//Andrew's age

  MPI_Status status;
  MPI_Request* recv_request=new MPI_Request[numprocs];
  MPI_Request* request = new MPI_Request[numprocs];
  MPI_Request* request2 = new MPI_Request[numprocs];
  extern MPI_Datatype ENRICHED_INFO;

  for(i=0; i<numprocs; i++)
    send_counter[i]=0;

  //counting how many packs will be sent to each proc
  enriched_current=enriched_start->next;
  while(enriched_current)
    {
      send_counter[enriched_current->target_proc]++;
      enriched_current=enriched_current->next;
    }
  
  //creating the space for the packs

  int index1=0;
  int* index2=new int[numprocs];
  for(i=0; i<numprocs; i++)
    {
      send_buffer[i]=new enriched_neighbor_pack[send_counter[i]];
      index2[i]=0;
    }
  
  //filling up the packs
  enriched_current=enriched_start->next;

  while(enriched_current)
    {
      index1=enriched_current->target_proc;
      
      send_buffer[index1][index2[index1]].order=enriched_current->order;
      for(int j=0; j<KEYLENGTH; j++)
	{
	  send_buffer[index1][index2[index1]].target_element[j]=enriched_current->target_element[j];
	  send_buffer[index1][index2[index1]].neighbor[j]=enriched_current->neighbor[j];

	}

      enriched_current=enriched_current->next;
      index2[index1]++;

    }

  for(i=1; i<numprocs; i++)
    {

      mate_proc=i^myid;//bitwise xor
      assert(mate_proc < numprocs);    

      MPI_Isend(&send_counter[mate_proc], 1, MPI_INT, mate_proc, atag, MPI_COMM_WORLD, (request+mate_proc));
      MPI_Isend(send_buffer[mate_proc], send_counter[mate_proc], ENRICHED_INFO, mate_proc, enrichedtag, MPI_COMM_WORLD, (request2+mate_proc));
      
      MPI_Recv(&recv_counter[mate_proc], 1, MPI_INT, mate_proc,atag, MPI_COMM_WORLD, &status);
      recv_buffer[mate_proc]=new enriched_neighbor_pack[recv_counter[mate_proc]];
      
      MPI_Irecv(recv_buffer[mate_proc], recv_counter[mate_proc] , ENRICHED_INFO, mate_proc, enrichedtag, MPI_COMM_WORLD, &recv_request[mate_proc]);

	 	 	    
    }//end of communication

  for(i=1; i<numprocs; i++)
    {
      mate_proc=i^myid;
      MPI_Wait(&recv_request[mate_proc], &status);
      
      for(int k=0; k<recv_counter[mate_proc]; k++)
	{

	  EmTemp=(Element*)HT_Elem_Ptr->lookup(recv_buffer[mate_proc][k].target_element);
	  assert(EmTemp);
	  my_gen=EmTemp->get_gen();
	  which_neighbor=EmTemp->which_neighbor(recv_buffer[mate_proc][k].neighbor);
	  neighbor_gen=*(EmTemp->get_neigh_gen()+which_neighbor);
	  neigh_order=recv_buffer[mate_proc][k].order;
	  my_order=*(EmTemp->get_order()+which_neighbor%4);

	  final_order = (my_order>=neigh_order) ? my_order:neigh_order;
	  if(final_order > my_order)
	    {
	      EmTemp->put_order(which_neighbor%4, final_order);
	      NodeTemp=(Node*)HT_Node_Ptr->lookup(EmTemp->getNode()+(which_neighbor%4+4)*KEYLENGTH);
	      assert(NodeTemp->get_order()==my_order);
	      NodeTemp->put_order(final_order);

	      //added acbauer
	      if(NodeTemp->getinfo() == S_S_CON)
		{
		  Node* TargetNode=(Node*)
		    (HT_Node_Ptr->lookup(EmTemp->getNode()+(which_neighbor%4)*KEYLENGTH));
		  
		  if(TargetNode->getinfo()!=S_C_CON)
		    
		    if(which_neighbor%4 !=3) TargetNode=(Node*)(HT_Node_Ptr->lookup(EmTemp->getNode()+(which_neighbor%4-1)*KEYLENGTH));
		    else TargetNode=(Node*)(HT_Node_Ptr->lookup(EmTemp->getNode()));		

		  TargetNode->put_order(final_order);
		}  //done acbauer
	    

	      if(neighbor_gen < my_gen)//if there is another little `twin'
		{
		  int next, pre;
		  next=(which_neighbor%4)+1;
		  if(next==4) next=0;
		  EmTemp2=(Element*)HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+next*KEYLENGTH);
		  if(!EmTemp2 || EmTemp->get_gen()!=EmTemp->get_gen()) {
		    int* elm_loc = EmTemp->get_elm_loc();
		    int* elm_loc2 = EmTemp2->get_elm_loc();
		    if( elm_loc[0]/2 != elm_loc2[0]/2 || elm_loc[1]/2 != elm_loc2[1]/2)  //elements don't have the same parents
		      {
			pre=(which_neighbor%4)-1;
			if(pre==-1)pre=3;
			EmTemp2=(Element*)HT_Elem_Ptr->lookup(EmTemp->get_neighbors()+pre*KEYLENGTH);
		      }
		  }
		  assert(EmTemp2);
		  
		  //assert(*(EmTemp2->get_neighbors()+(which_neighbor%4)*KEYLENGTH)==recv_buffer[mate_proc][k].neighbor[0]);
		  //assert(*(EmTemp2->get_neighbors()+(which_neighbor%4)*KEYLENGTH+1)==recv_buffer[mate_proc][k].neighbor[1]);
		  EmTemp2->put_order(which_neighbor%4, final_order);
		  
		  NodeTemp=(Node*)HT_Node_Ptr->lookup(EmTemp2->getNode()+(which_neighbor%4+4)*KEYLENGTH);
		  assert(NodeTemp->get_order()==my_order);
		  NodeTemp->put_order(final_order);
		  
		}
	    }  

	}
      

    }
  for(i=0;i<numprocs;i++) {
    if(send_buffer[i])
      delete [](send_buffer[i]);
    if(recv_buffer[i])
      delete [](recv_buffer[i]);
  }
  delete []send_buffer;
  delete []recv_buffer;
  delete []send_counter;
  delete []recv_counter;
  delete []recv_request;
  delete []index2;
  for(i=1; i<numprocs; i++)    {
    mate_proc=i^myid;//bitwise xor
    MPI_Wait((request+mate_proc), &status);
    MPI_Wait((request2+mate_proc), &status);
  }
    delete []request;
    delete []request2;

  return;
}
