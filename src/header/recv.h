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
 * $Id: recv.h 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifndef RECV_H
#define RECV_H

class Recv
{
 public:
  Element*   targetP;
  unsigned   sender[KEYLENGTH];
  unsigned   target[KEYLENGTH];
  unsigned   sender_son1[KEYLENGTH];
  unsigned   sender_son2[KEYLENGTH];
  int        sender_order[2];
  int        sender_gen;
  int        sender_refined;
  int        side;
  int        sender_id;
  Recv*      pre;
  Recv*      next;
  int        myid;
  void       identify ();

  Recv(HashTable* ht_elem_ptr, unsigned* recv_buf, int assoc, int iam)
    {
      int i, j;
      side = INIT;
      myid = iam;
      for(i=0;i<KEYLENGTH;i++) sender[i] = *(recv_buf+i);
      for(i=KEYLENGTH;i<2*KEYLENGTH;i++) sender_son1[i-KEYLENGTH] = *(recv_buf+i);
      for(i=2*KEYLENGTH;i<3*KEYLENGTH;i++) sender_son2[i-2*KEYLENGTH] = *(recv_buf+i);
      sender_order[0] = *(recv_buf+i);
      sender_order[1] = *(recv_buf+i+1);
      sender_gen = *(recv_buf+i+2);
      sender_refined = *(recv_buf+i+3);
      for(i=0;i<KEYLENGTH;i++) target[i] = recv_buf[5*KEYLENGTH+i];
      targetP = (Element*)ht_elem_ptr->lookup(target);
      sender_id = assoc;
      identify();
      pre =NULL;
      next = NULL;
    }
  Recv(){
    next   = NULL;
    pre    = NULL;
  }
  ~Recv(){
     if(next)
       next->pre = pre;
     if(pre)
       pre->next = next;
  }

};

#endif
