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

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

#include <Python.h>

#include "../header/hpfem.h"

#if HAVE_HDF5_H
#include "../header/hd5calls.h"
#endif

#include "../header/titan2d.h"

#include "../header/titan2d_utils.h"

#include "../header/titan_simulation.h"

extern "C" void init_cxxtitan();

int main(int argc, char *argv[])
{
    int myid, master, numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);
    
    //set program full name
    int i;
    char executable[1028];
    char buffer[1028];
    
    i = readlink("/proc/self/exe", executable, 1028);
    executable[i] = '\0';
    Py_SetProgramName(executable);
    
    Py_Initialize();
    
    if(myid == 0)
    {
        printf("Titan2d %d.%d.%d\n", TITAN2D_VERSION_MAJOR, TITAN2D_VERSION_MINOR, TITAN2D_VERSION_REVISION);
        printf("Executable location: %s\n", executable);
        printf("-------------------------------------------------------------------------------\n");
        printf("\n");
    }
    
    PyRun_SimpleString("import os");
    PyRun_SimpleString("import sys");
    
    sprintf(buffer,
            "sys.path.insert(1,os.path.abspath(os.path.join('%s','..','..','lib','python2.7','site-packages','titan')))",
            executable);
    PyRun_SimpleString(buffer);
    
    init_cxxtitan();
    PyRun_SimpleString("from titan import *");
    
    //init some titan staff
    sfc_key_null=0;
    sfc_key_zero=0;
    
    if(argc > 0)
    {
        Py_Main(argc, argv);
    }
    else if(myid == 0)
    {
        printf("Usage:\n");
        printf("\t%s <run_script>\n", argv[0]);
    }
    
    Py_Finalize();
    if(myid == 0)
    {
        printf("Done\n");
    }
    MPI_Finalize();
    return 0;
}


