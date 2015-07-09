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
 * $Id: hpfem.C 235 2012-03-29 07:04:35Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#if HAVE_HDF5_H
#include "../header/hd5calls.h"
#endif

//#define DEBUG
//#define LOAD_BAL_DEBUG  //turns on a whole mess of mpi barriers so it makes run time more sensitive to load imbalances i.e. more sensitive to the load balance weights, it just makes it easier to adjust the constants.
//#define PERFTEST
#define TARGETPROC  -1


#ifdef TWO_PHASES
void checkelemnode2(HashTable *El_Table, HashTable *NodeTable, int myid, FILE *fpdebug, double loc)
{
    unsigned elemdebugkey2a[2] =
    { 695876266, 2863311530 };
    unsigned nodedebugkey2a[2] =
    { 695852110, 3303820997 };
    
    if(fpdebug == NULL)
        fpdebug = stdout;
    
    if(myid == TARGETPROC)
    {
        fprintf(fpdebug, "**************************\nmyid=%d location=%g\n", myid, loc);
        ElemBackgroundCheck(El_Table, NodeTable, elemdebugkey2a, fpdebug);
        //NodeBackgroundCheck(El_Table,NodeTable,nodedebugkey2a,fpdebug);
        
        fflush(fpdebug);
    }
    
    return;
}
#endif

