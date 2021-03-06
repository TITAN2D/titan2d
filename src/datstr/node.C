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
 * $Id: node.C 233 2012-03-27 18:30:40Z dkumar $ 
 */
//#define DEBUG_SAVE_NODE
#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#ifndef NODE_C
#define NODE_C
#include "../header/elenode.hpp"
#include "../gisapi/GisApi.h"
#include "../header/properties.h"

#include "../header/ticore/omp_mpi.hpp"
#include <assert.h>

void Node::init()
{
    id(0);
    connection_id(-1);
    info(INIT);
    
    for(int i = 0; i < DIMENSION; i++)
        coord(i,0.0);
    
    //set_key(sfc_key_null);
    zero_flux();
    elevation(0.0);
}

void Node::init(const SFC_Key& keyi, double* coordi, MatProps* matprops_ptr)
{
    int i;
    id(0); /* keith added this so save_node wouldn't write an uninitialized
     variable and valgrind wouldn't flag an error.  id is used in
     ../repartition/BSFC_update_and_send_elements.C */
    
    connection_id(-1);
    
    info(INIT);
    
    for(i = 0; i < DIMENSION; i++)
        coord(i,coordi[i]);
    
    set_key(keyi);
    zero_flux();
    // find the max resolution of the GIS info and then get the elevation at this node
    double resolution = 0;
    double xcoord = coord(0) * (matprops_ptr->scale.length);
    double ycoord = coord(1) * (matprops_ptr->scale.length);
    i = Get_max_resolution(&resolution);
    if(i != 0)
    {
        printf("error in Get_max_resolution\n");
        assert(0);
    }
    i = Get_elevation(resolution, xcoord, ycoord, elevation_ref());
    if(i != 0)
    {
        printf("error in Get_elevation(%d) r=%g (x,y)=(%g,%g) e=%g\n", i, resolution, xcoord, ycoord, elevation());
        assert(0);
    }
    elevation(elevation() / matprops_ptr->scale.length);
    /*  if((unsigned) 1548032885 == key[0])
     printf("created the missing node...\n"); */

    /*
     if((coord[0]==0)||(coord[1]==0)){
     printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
     assert(coord[0]*coord[1]);
     }
     */
}

void Node::init(const SFC_Key& keyi, const double *coordi, const int inf, const int ord, const MatProps *matprops_ptr)  //for refined
{
    int i;
    id(0); /* keith added this so save_node wouldn't write an uninitialized
     variable and valgrind wouldn't flag an error.  id is used in
     ../repartition/BSFC_update_and_send_elements.C */
    
    connection_id(-1);
    info(INIT);
    
    for(i = 0; i < DIMENSION; i++)
        coord(i, coordi[i]);
    
    set_key(keyi);
    
    info(inf);
    //geoflow stuff
    zero_flux();
    // find the max resolution of the GIS info and then get the elevation at this node
    double resolution = 0;
    i = Get_max_resolution(&resolution);
    if(i != 0)
    {
        printf("error in Get_max_resolution\n");
        assert(0);
    }
    double xcoord = coord(0) * (matprops_ptr->scale.length);
    double ycoord = coord(1) * (matprops_ptr->scale.length);
    i = Get_elevation(resolution, xcoord, ycoord, elevation_ref());
    if(i != 0)
    {
        printf("error in Get_elevation\n");
        assert(0);
    }
    elevation(elevation() / matprops_ptr->scale.length);
    /*  if((unsigned) 1548032885 == key[0])
     printf("created the missing node111111...\n");*/

    /*
     if((coord[0]==0)||(coord[1]==0)){
     printf("node={%u,%u} coord=(%20g,%20g)\n",key[0],key[1],coord[0],coord[1]);
     assert(coord[0]*coord[1]);
     }
     */
    return;
}

void Node::init(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada)
{
    int i;
    id(0); /* keith added this so save_node wouldn't write an uninitialized
     variable and valgrind wouldn't flag an error.  id is used in
     ../repartition/BSFC_update_and_send_elements.C */
    
    connection_id(-1);
    
    info(INIT);
    
    for(i = 0; i < DIMENSION; i++)
        coord(i, coordi[i]);
    
    set_key(keyi);
    
    info(inf);
    //geoflow stuff
    zero_flux();
    elevation(elev);
    /*
     if((coord[0]==0)||(coord[1]==0)){
     printf("inode=%d node={%u,%u} coord=(%20g,%20g)\n",yada,key[0],key[1],coord[0],coord[1]);
     assert(coord[0]*coord[1]);
     }
     */
    return;
}

Node::~Node()
{
}

void Node::elevation(MatProps* matprops_ptr)
{
    double resolution = 0;
    double xcoord = coord(0) * (matprops_ptr->scale.length);
    double ycoord = coord(1) * (matprops_ptr->scale.length);
    int i = Get_max_resolution(&resolution);
    if(i != 0)
    {
        printf("error in Get_max_resolution\n");
        assert(0);
    }
    i = Get_elevation(resolution, xcoord, ycoord, elevation_ref());
    if(i != 0)
    {
        printf("error in Get_elevation\n");
        assert(0);
    }
    elevation(elevation() / matprops_ptr->scale.length);
    
}

void Node::save_node(FILE* fp)
{
    
    FourBytes temp4;
    EightBytes temp8;
    unsigned writespace[13];
    
    int Itemp = 0, itemp;
    sfc_key_write_to_space(key(),writespace,Itemp);
    assert(Itemp == 2);
/*#ifdef DEBUG_SAVE_NODE
    FILE *fpdb=fopen("save_node.debug","w");
    fprintf(fpdb,"key=%u %u\n",key[0],key[1]);
#endif*/
    
    for(itemp = 0; itemp < DIMENSION; itemp++)
    {
        temp8.d = coord(itemp);
        writespace[Itemp++] = temp8.u[0];
        writespace[Itemp++] = temp8.u[1];
    }
#ifdef DEBUG_SAVE_NODE
    fprintf(fpdb,"coord=%g %g\n",coord(0),coord(1));
#endif
    assert(Itemp == 6);
    
    temp4.i = id();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 7);
#ifdef DEBUG_SAVE_NODE
    fprintf(fpdb,"id=%d\n",id());
#endif
    
    temp4.i = info();
    writespace[Itemp++] = temp4.u;
    assert(Itemp == 8);
#ifdef DEBUG_SAVE_NODE
    fprintf(fpdb,"info=%d\n",info());
#endif
    
    /* these are Legacy and are not used
     temp4.i=order;
     writespace[Itemp++]=temp4.u;
     assert(Itemp==9);
     #ifdef DEBUG_SAVE_NODE
     fprintf(fpdb,"order=%d\n",order); 
     #endif

     temp4.i=dof[0];
     writespace[Itemp++]=temp4.u;
     temp4.i=dof[1];
     writespace[Itemp++]=temp4.u;
     assert(Itemp==11);
     #ifdef DEBUG_SAVE_NODE
     fprintf(fpdb,"dof=%d %d\n",dof[0],dof[1]); 
     #endif

     temp4.i=glnum;
     writespace[Itemp++]=temp4.u;
     assert(Itemp==12);
     #ifdef DEBUG_SAVE_NODE
     fprintf(fpdb,"glnum=%d\n",glnum); 
     #endif

     temp4.i=reconstructed;
     writespace[Itemp++]=temp4.u;
     assert(Itemp==13);
     #ifdef DEBUG_SAVE_NODE
     fprintf(fpdb,"reconstructed=%d\n",reconstructed); 
     #endif
     */

#ifdef DEBUG_SAVE_NODE
    fclose(fpdb);
#endif
    fwrite(writespace, sizeof(unsigned), Itemp, fp);
    
    return;
}

void Node::init(FILE* fp, MatProps* matprops_ptr)
{      
    FourBytes temp4;
    EightBytes temp8;
    //unsigned readspace[13];
    unsigned readspace[8];
    int Itemp = 0, itemp;
    
    //fread(readspace,sizeof(unsigned),13,fp);
    fread(readspace, sizeof(unsigned), 8, fp);
    
    //KEYLENGTH should be 2 but put it in a loop to make it generic.
    set_key(sfc_key_read_from_space(readspace,Itemp));
    assert(Itemp == 2);
    
    //DIMENSION should be 2 but put it in a loop to make it generic.
    for(itemp = 0; itemp < DIMENSION; itemp++)
    {
        temp8.u[0] = readspace[Itemp++];
        temp8.u[1] = readspace[Itemp++];
        coord(itemp, temp8.d);
    }
    assert(Itemp == 6);
    
    temp4.u = readspace[Itemp++];
    id(temp4.i);
    assert(Itemp == 7);
    
    temp4.u = readspace[Itemp++];
    info(temp4.i);
    assert(Itemp == 8);
    
    /* these are legacy and are not used
     temp4.u=readspace[Itemp++];
     order=temp4.i;
     assert(Itemp==9);

     temp4.u=readspace[Itemp++];
     dof[0]=temp4.i;
     temp4.u=readspace[Itemp++];
     dof[1]=temp4.i;
     assert(Itemp==11);

     temp4.u=readspace[Itemp++];
     glnum=temp4.i;
     assert(Itemp==12);

     temp4.u=readspace[Itemp++];
     reconstructed=temp4.i;
     assert(Itemp==13);
     */
    // find the max resolution of the GIS info and then get the elevation at this node
    double resolution = 0;
    double xcoord = coord(0) * (matprops_ptr->scale.length);
    double ycoord = coord(1) * (matprops_ptr->scale.length);
    int i = Get_max_resolution(&resolution);
    if(i != 0)
    {
        printf("error in Get_max_resolution\n");
        assert(0);
    }
    i = Get_elevation(resolution, xcoord, ycoord, elevation_ref());
    if(i != 0)
    {
        printf("error in Get_elevation\n");
        assert(0);
    }
    elevation(elevation() / matprops_ptr->scale.length);
    
    return;
}

#endif

