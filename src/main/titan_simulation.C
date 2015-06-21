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
 */

#include <stdio.h>
#include <mpi.h>
#include "../header/titan_simulation.h"

MaterialMap::MaterialMap()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}
MaterialMap::~MaterialMap()
{

}
void MaterialMap::print0()
{
    if(myid!=0){
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    int i;
    printf("Material map:\n");
    for(i=0;i<name.size();i++){
        printf("%d %s %f %f\n",i,name[i].c_str(),intfrict[i],bedfrict[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
int MaterialMap::get_material_count(){
    return name.size();
}
cxxTitanSimulation::cxxTitanSimulation()
{
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gis_format=-1;

    topomain="";
    toposub="";
    topomapset="";
    topomap="";

    region_limits_set=false;

    min_location_x=0.0;
    max_location_x=0.0;
    min_location_y=0.0;
    max_location_y=0.0;

    MPI_Barrier(MPI_COMM_WORLD);
}
cxxTitanSimulation::~cxxTitanSimulation()
{

}
int hpfem();
void cxxTitanSimulation::run()
{
	printf("cxxTitanSimulation::run\n");
	MPI_Barrier(MPI_COMM_WORLD);
	hpfem();
	MPI_Barrier(MPI_COMM_WORLD);

}

void cxxTitanSimulation::input_summary()
{
    if(myid!=0){
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    printf("GIS:\n");
    printf("\tgis_format: %d\n",gis_format);

    printf("\ttopomain: %s\n",topomain.c_str());
    printf("\ttoposub: %s\n",toposub.c_str());
    printf("\ttopomapset: %s\n",topomapset.c_str());
    printf("\ttopomap: %s\n",topomap.c_str());

    printf("\tregion_limits_set %d\n",(int)region_limits_set);

    MPI_Barrier(MPI_COMM_WORLD);
}
