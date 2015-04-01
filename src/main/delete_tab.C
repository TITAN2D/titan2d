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
 * $Id: delete_tab.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 

#include "../header/hpfem.h"

void Delete_Table(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr)
{

  int i, j, k;
  HashEntryPtr entryp;

  int elements = HT_Elem_Ptr->get_no_of_buckets();
  int nodes    = HT_Node_Ptr->get_no_of_buckets();   
  for(i=0;i<elements;i++){ 
    entryp = *(HT_Elem_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	Element* EmTemp = (Element*)(entryp->value);
	delete EmTemp;
	entryp = entryp->next;
      }
  }

  for(i=0;i<nodes;i++){ 
    entryp = *(HT_Node_Ptr->getbucketptr() + i); 
    while(entryp)
      { 
	Node* NdTemp = (Node*)(entryp->value);
	delete NdTemp;
	entryp = entryp->next;
      }
  }

  delete HT_Elem_Ptr;
  delete HT_Node_Ptr;

}
