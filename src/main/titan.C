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

#define __STDC_FORMAT_MACROS 
#include <inttypes.h>

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "../header/ticore/omp_mpi.hpp"

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
    int myid=0, master, numprocs=1;
    int namelen;
    
    IF_MPI(MPI_Init(&argc, &argv);)
    IF_MPI(MPI_Comm_size(MPI_COMM_WORLD, &numprocs);)
    IF_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &myid);)
    
    //set program full name
    int i;
    char executable[1028];
    char python_home[1028];
    char buffer[1028];
    
#ifdef __linux__
    i = readlink("/proc/self/exe", executable, 1028);
    executable[i] = '\0';
#endif
#ifdef __APPLE__
    uint32_t namesize=1024;
    _NSGetExecutablePath(executable, &namesize);
    executable[namesize] = '\0';
#endif
#ifdef _WIN32
    somthing with GetModuleFileName()
#endif
    Py_SetProgramName(executable);

//#define MAKE_PORTABLE
#ifdef MAKE_PORTABLE
    strncpy(python_home,executable,1028);
    printf("%s\n",python_home);
    int count=0;
    for(i=strlen(python_home);i>=0;--i)
    {
        if(python_home[i]=='/')
        {
            python_home[i]='\0';
            ++count;
            if(count==2)break;
        }
    }

    printf("%s\n",python_home);
    sprintf(python_home+strlen(python_home),"/lib/titan2d_dep");
    printf("%s\n",python_home);
    Py_SetPythonHome(python_home);
#else

#endif
    printf("Executable location: %s\n", executable);
    printf("python home: %s\n", Py_GetPythonHome());
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
    
    //init some titan staff should be in function titan_init
    sfc_key_null=0;
    sfc_key_zero=0;
    elementsHashTable=nullptr;
    nodeHashTable=nullptr;
    
    init_TiH5();

    bool printUsage=false;

    //variables for restart command line
    bool restart=false;
    string restartFilename="";
    string new_gis_main="";
    int addIterations=-1;
    double addTime=-1.0;

    if(argc > 1)
    {
    	int argc4py=0;
    	char **argv4py=new char*[argc];

    	int num_threads=1;

    	for(int i=0;i<argc;i++){
    		string arg(argv[i]);
    		//strip off arguments which should not make to python
            if(arg=="-h")
            {
                printUsage=true;//print usage and exit
                break;
            }
    		if(arg=="-nt")
    		{

    			if((i+1>=argc) || (!isdigit(argv[i+1][0])))
    			{
    				argc4py=0;
    				printUsage=true;//print usage and exit
    				printf("Error: Incorrect command line format\n");
    				break;
    			}
    			++i;

    			num_threads=atoi(argv[i]);

    	       //set openmp threads
#ifdef _OPENMP
    			printf("Set threads number to %d\n\n",num_threads);
    			omp_set_num_threads(num_threads);
#else
    			if(num_threads>1)
    				printf("Warning: This version was compiled without openmp support!\n");

    			num_threads=1;
    			printf("Set threads number to %d\n\n",num_threads);
#endif
    			continue;
    		}
            if(arg=="-restart")
            {
                for(int i=0;i<argc;i++)
                {
                    arg=argv[i];
                    if(arg=="-add-iter")
                    {
                        ++i;
                        arg=argv[i];
                        addIterations=stoi(arg);
                    }
                    else if(arg=="-add-time")
                    {
                        ++i;
                        arg=argv[i];
                        addTime=stoi(arg);
                    }
                    else if(arg=="-gis-main")
                    {
                        ++i;
                        arg=argv[i];
                        new_gis_main=arg;
                    }
                    else
                    {
                        restartFilename=argv[i];
                    }
                }
                if(restartFilename=="")
                {
                    printUsage=true;//print usage and exit
                    printf("Error: Incorrect command line format\n");
                    break;
                }

                restart=true;
                continue;
            }
    		argv4py[argc4py]=argv[i];
    		++argc4py;
    	}

#ifdef _OPENMP
    	threads_number=omp_get_max_threads();
    	printf("Max threads number is %d\n\n",threads_number);
#else
        threads_number=1;
#endif

        if(printUsage)
        {
            //print usage and exit
        }
        else if(restart)
        {
            //restart
            printf("restart %s %d %f\n",restartFilename.c_str(),addIterations,addTime);
            cxxTitanSimulation sim;
            sim.load_restart(restartFilename,new_gis_main);

            if(addIterations>0)
            {
                sim.timeprops.maxiter+=addIterations;
            }
            if(addTime>0)
            {
                sim.timeprops.maxtime+=addTime/sim.timeprops.TIME_SCALE;
            }
            sim.run(true);

        }
        else if(argc4py > 1)
        {
            //execute user script
    		Py_Main(argc4py, argv4py);
    	}
    	else
    	{
    	    //print usage and exit
    		printUsage=true;
    	}

        delete [] argv4py;
    }
    else
    {
    	printUsage=true;//print usage and exit
    }
    if(printUsage)
    {
		if(myid == 0)
		{
			printf("Usage:\n");
			printf("\t%s [-nt <number of threads>] [-restart [-add-iter <iter>] [-add-time <time>] [-gis-main <new position of gis>] <restart file> | <run_script>]\n", argv[0]);
		}
    }
    
    Py_Finalize();
    if(myid == 0)
    {
        printf("Done\n");
    }
    IF_MPI(MPI_Finalize());
    return 0;
}


