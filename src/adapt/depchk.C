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
 * $Id: depchk.C 127 2007-06-07 19:48:25Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

#define NumTriggerRef 256

void depchk(Element* EmTemp, HashTable* El_Table, HashTable* NodeTable, 
	    int* ifg, ElemPtrList* RefinedList)

  /*---
    refined[] stores the address of ready-for-refinement element of the sub-domain
    refined_temp[] stores the address of ready-for-refinement element triggered by one element refinement
    count is counting the number of refinement of the subdomain
    j is counting the number of refinement triggered by one element refinement
    ---------------*/
{

  int i, j, k;
  Element* element;
  Element* Neigh;
  ElemPtrList TempList(384);
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  TempList.add(EmTemp); 

  j = 0; k = 0;
  element = EmTemp;//-- EmTemp is the trigger of this round of refinement

  while(element&&(j<NumTriggerRef)) {//--element is temporary varible
   
    for(i=0;i<4;i++) {//-- checking the four neighbors to identify which must be refined
	
      int neigh_proc = *(element->get_neigh_proc()+i);
      
      if((neigh_proc!=-1)&&(neigh_proc!=-2)) {//-- if there is a neighbor
	   
	      
	Neigh = (Element*)(El_Table->lookup(element->get_neighbors()+i*KEYLENGTH));
	/*
	if(!Neigh){
	  printf("this element is {%u,%u} and is missing his %dth neighbor\n ",
		 *(element->pass_key()+0),*(element->pass_key()+1),i);
	  ElemBackgroundCheck(El_Table,NodeTable,element->pass_key(),stdout);
	  fflush(stdout);
	  printf("that is... the missing element is {%u,%u}\n",
		 *(element->get_neighbors()+i*KEYLENGTH+0),
		 *(element->get_neighbors()+i*KEYLENGTH+1));
	  fflush(stdout);
	  ElemBackgroundCheck(El_Table,NodeTable,element->get_neighbors()+i*KEYLENGTH,stdout);
	  fflush(stdout);
	}
	*/
	      
	//assert(Neigh);
	if(Neigh != NULL && neigh_proc == myid) {//-- if this neighbor is in the same proc as element is
		
	  if(element->get_gen()>Neigh->get_gen()) {
	    //-- if the neighbor is bigger, then it must be refined
	    
	    if((Neigh->get_adapted_flag()==NOTRECADAPTED)||
	       (Neigh->get_adapted_flag()==NEWFATHER)
	       ) {
	      int flag = 1;
	      for(int m=0;m<TempList.get_num_elem();m++)
		if(compare_key(TempList.get_key(m), 
			       Neigh->pass_key())) {
		    flag = 0; 
		    break;
		}
	      
	      if(flag) {//-- if this neighbor has not yet been marked
		
		j++;
		TempList.add(Neigh);
	      }
	    }
	    else if(Neigh->get_adapted_flag()!=OLDFATHER){
	      *ifg=0;
	      TempList.trashlist();
	      break;
	    }
	  }
	  
	}
	else {//-- need neighbor's generation infomation
	  
	  if(element->get_gen()>*(element->get_neigh_gen()+i)) {//--stop this round of refinement
	    
	    *ifg = 0;
	    TempList.trashlist();
	    break;
	  }
	}
      }			    
      
    }
    if(!*ifg) break;
    k++;
    element = TempList.get(k);//--check next
  }
  

  //copy TempList to RefinedList
  if(*ifg) {
    
    if(j<NumTriggerRef)//-- NumTriggerRef is the maximum tolerence of related refinement
      for(int m=0;m<TempList.get_num_elem();m++) {
	int sur = 0;
	for(int mi=0;mi<RefinedList->get_num_elem();mi++)
	  if(sur=compare_key(RefinedList->get_key(mi), 
			     TempList.get_key(m)))
	    break;
	if(!sur)
	  RefinedList->add(TempList.get(m));
      }
    else {
	
      *ifg = 0;//-- refuse to do the refinement
    }
  }

  for(int m=0;m<TempList.get_num_elem()-1;m++)
    for(int mi=m+1;mi<TempList.get_num_elem();mi++)
      assert(!compare_key(TempList.get_key(m),TempList.get_key(mi)));
 
  return;
}



#ifdef DISABLED
void depchk(Element* EmTemp, HashTable* El_Table, int* ifg, Element* refined[], int* count)

  /*---
    refined[] stores the address of ready-for-refinement element of the sub-domain
    refined_temp[] stores the address of ready-for-refinement element triggered by one element refinement
    count is counting the number of refinement of the subdomain
    j is counting the number of refinement triggered by one element refinement
    ---------------*/
{

  int i, j, k;
  Element* element;
  Element* Neigh;
  Element* refined_temp[128];
  void* p;
  int myid, numprocs;
  unsigned send_buf[4*KEYLENGTH];
  unsigned recv_buf[4*KEYLENGTH];

  MPI_Status     status;
  MPI_Request    request;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  for(j=1;j<128;j++) refined_temp[j] = NULL;
  refined_temp[0] = EmTemp;
 
  j = 0; k = 0;
  element = EmTemp;//-- EmTemp is the trigger of this round of refinement

  while(element&&(j<10))//--element is temporary varible
    {
      for(i=0;i<4;i++)//-- checking the four neighbors to identify which must be refined
	{
	  int neigh_proc = *(element->get_neigh_proc()+i);
	  
	  if((neigh_proc!=-1)&&(neigh_proc!=-2))//-- if there is a neighbor
	    {
	      
	      Neigh = (Element*)(El_Table->lookup(element->get_neighbors()+i*KEYLENGTH));
	      
	      assert(Neigh);
	      if(Neigh != NULL && neigh_proc == myid) //-- if this neighbor is in the same proc as element is
		{
		  if((!Neigh->get_refined_flag())&&element->get_gen()>Neigh->get_gen())//-- if the neighbor is bigger, then it must be refined
		    {
		      int flag = 1; int m = 0;
		      while(refined_temp[m]) 
			{ 
			  if(compare_key(refined_temp[m]->pass_key(), Neigh->pass_key()))
			    {flag = 0; break;}
			  else m++;
			  assert(m<128);
			}
		      
		      if(flag)//-- if this neighbor has not yet been marked
			{
			  j++;
			  assert(j<128);
			  refined_temp[j] = Neigh;
			}
		    }
		}
	      else//-- need neighbor's generation infomation
		{
		  if(element->get_gen()>*(element->get_neigh_gen()+i))//--stop this round of refinement
		    {
		      *ifg = 0;
		      for(int m=0;m<128;m++) refined_temp[m] = NULL;
		      break;
		    }
		}
	    }			    
			    
	}
      if(!*ifg) break;
      k++;
      element = refined_temp[k];//--check next
    }
  
  if(*ifg)
    {
      if(j<10)//-- 10 is the maximum tolerence of related refinement
	{
	  int m = 0;
	  while(refined_temp[m])
	    {
	      int sur = 0; int mi = 0;
	      while(refined[mi])
		{
		  if(compare_key(refined[mi]->pass_key(), refined_temp[m]->pass_key())) //-- KEYLENGTH should be considered
		    {
		      sur = 1;
		      break;
		    }
		  mi++;
		}
	      if(!sur)
		{
		  refined[*count] = refined_temp[m];
		  *count = *count+1;
		}
	      m++;
	    }
	}
      else
	{
	  *ifg = 0;//-- refuse to do the refinement
	}
    }
}

#endif

