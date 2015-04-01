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
# include <config.h>
#endif
 
#include "../header/hpfem.h"

void htflush(HashTable* ht_elem_ptr, HashTable*  ht_node_ptr, int option)
{

  int i, j, k;
  Element*       EmTemp;
  Node*          NdTemp;
  unsigned       KeyTemp[KEYLENGTH];
  HashEntryPtr   entryp;
  unsigned*      keyP;
  void*          p;
  int*           dofP;
  double*        sol;
  int e_buckets=ht_elem_ptr->get_no_of_buckets();
  int n_buckets=ht_node_ptr->get_no_of_buckets();

  switch(option)
    {
    case 1:
      for(i=0;i<e_buckets;i++)
	{
	  entryp = *(ht_elem_ptr->getbucketptr() + i);
	  while(entryp)
	    { 
	      EmTemp = (Element*)(entryp->value);
	      EmTemp->put_new_old(OLD);
	      entryp = entryp->next;
	    }
	}
      break;
    case 2:
      for(i=0;i<n_buckets;i++)
	{
	  entryp = *(ht_node_ptr->getbucketptr() + i);
	  while(entryp)
	    { 
	      NdTemp = (Node*)(entryp->value);
	      //NdTemp->putinfo(INIT);
	      NdTemp->putdof(INIT,INIT );
	      NdTemp->putglnum(INIT);
	      NdTemp->put_reconstructed(0);


	      sol = NdTemp->getsol();
	      /*if(sol)
		delete sol;*/

	      if((!NdTemp->get_sol_deleted())&&(sol != NULL))
		{		  
		  sol = NdTemp->getsol();
		  delete sol;
		  NdTemp->put_sol_deleted(1);
		}

	      entryp = entryp->next;
	    }
	}
      for(i=0;i<e_buckets;i++)
	{
	  entryp = *(ht_elem_ptr->getbucketptr() + i);
	  while(entryp)
	    { 
	      EmTemp = (Element*)(entryp->value);
	      for(j=0;j<8;j++) {
		EmTemp->put_recv_flag(j, 0);
		EmTemp->put_send_flag(j, 0);
	      }
	      entryp = entryp->next;
	    }
	}
    }
}
