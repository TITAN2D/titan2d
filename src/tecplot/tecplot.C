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
 * $Id: tecplot.C 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#define TECPLOTASCII
#define TARGETPROCA -1
#define TARGETPROCB -1

#ifndef TECPLOTASCII
//generates a corrupt binary tecplot file... it doesn't work
#define num_var 8
#endif

void get_elem_orient(Element * EmTemp, int *xm, int *xp, int *ym, int *yp);
int get_elem_elev(ElementsHashTable * NodeTable, MatProps * matprops, Element * EmTemp, double *elevation);

int get_ll_polygon(ElementsHashTable * El_Table, NodeHashTable * NodeTable, int myid, Element * EmArray[4]);
int get_ur_tri(ElementsHashTable * El_Table, NodeHashTable * NodeTable, int myid, Element * EmArray[4]);

int print_bubble_node(ElementType elementType, FILE * fp, NodeHashTable * NodeTable, MatProps * matprops, Element * EmTemp);

void DumpString(FILE * fp, char *str);

void testkey(ElementsHashTable * El_Table, Element * EmTemp)
{
    Element *EmTemp2 = (Element *) El_Table->lookup(EmTemp->key());
    if(EmTemp != EmTemp2)
    {
        printf("Element can't look up self with key");
        printf("\n");
    }
    /*if((key[0] == 486561007) && (key[1] == 1010586963))
    {
        printf("danger key present");
        printf("\n");
    }*/
    
    return;
}

void testkey2(ElementsHashTable * El_Table)
{
    unsigned key[2];
    key[0] = 486561007;
    key[1] = 1010586963;
    /* if(((Element*) El_Table->lookup(key))==0) {
     printf("danger Element went bad");
     printf("\n");
     } */

    return;
}

void tecplotter(ElementType elementType,ElementsHashTable * El_Table, NodeHashTable * NodeTable, MatProps * matprops, TimeProps * timeprops,
                MapNames * mapnames, double v_star)
{
    int i_buck, i_neigh;          //indices
    int xp, xm, yp, ym;           //x plus, x minus, y plus, y minus directions
    int gen;          //level of grid refinement
    int numprocs;                 //number of processes
    int myid;        //process id's
            
    int done = 1;
    unsigned key[2];
    int TECTAG = 123;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    Element *EmTemp, *EmArray[4];
    
    /*There are 2 ways to draw triangles
     1) as a separate ZONE of triangles 
     2) as quads with duplicate points
     I chose the 2nd option.  Why?  Because I want all the elements
     from a given processor to be in a single zone.  Also duplicate 
     points for every element that uses them have been eliminated.
     That makes the "boundary" option in tecplot work right, but 
     more importantly the nodes the triangles use aren't repeated 
     (you get the triangles for free apart from the element 
     conectivity data) and nodes for quads are reduced by about a
     factor of 4.  This reduces file the size by about a factor of 4 
     when grid adaptation is used, and by about a factor of 3 when 
     it's not.  

     MODIFICATION: 20070115: Now Each And Every BUBBLE Node of Active 
     and Ghost Elements are printed once in the order they appear in 
     the Element HashTable.  This eliminates the need to sort (quick
     or otherwise) the elements speeding up tecplot, at the cost of 
     slightly increasing the file size on multiprocessor runs 
     (adding BUBBLE nodes of GHOST elements along the right and upper
     edges of the collection of elements this processor "owns")
     */
    //move_data(numprocs, myid, El_Table, NodeTable,timeprops);
    
    if(myid == TARGETPROCA)
    {
        printf("at tecplotter 1.0\n");
        fflush(stdout);
    }
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //count the number of BUBBLE nodes, and tecplot quads and tecplot 
    //triangles I need to draw
    int num_tec_node = 0, num_tec_tri = 0, num_tec_quad = 0;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            EmArray[0] = EmTemp;
            if((EmTemp->adapted_flag() != TOBEDELETED) && (abs(EmTemp->adapted_flag()) <= BUFFER))
            {
                num_tec_node++;
                EmTemp->set_ithelem(num_tec_node);
                
                if(EmTemp->adapted_flag() > TOBEDELETED)
                {
                    switch (get_ll_polygon(El_Table, NodeTable, myid, EmArray))
                    {
                        case 4:    //we should draw a lower left quad
                            num_tec_quad++;
                            break;
                        case 3:    //we should draw a lower left triangle
                            num_tec_tri++;
                            break;
                        case 0:    //we should NOT draw a lower left polygon
                            break;
                        default:   //something wierd happened report this to the "FBI"
                            assert(0);
                    }           //switch(get_ll_polygon(El_Table,NodeTable,myid,EmArray))
                    
                    //match the lower left triangle at the lower left corner of 
                    //a different subdomain with an upper right triangle in this 
                    //subdomain
                    if(get_ur_tri(El_Table, NodeTable, myid, EmArray))
                        num_tec_tri++;
                    
                    gen = EmTemp->generation();
                    //one triangle is added to an element for each side whose
                    //neighboring element is more refined than current element 
                    //provided the triangle wouldn't fall on a domain boundary
                    for(i_neigh = 0; i_neigh < 4; i_neigh++)
                        if((EmTemp->neigh_gen(i_neigh) > gen) && (EmTemp->neigh_proc(i_neigh) != INIT))
                            num_tec_tri++;
                }               //if(EmTemp->get_adapted_flag()>TOBEDELETED)
            }                   //if((EmTemp->get_adapted_flag()!=TOBEDELETED)&&
        }                       //while(entryp)
    }                           //for(i_buck=0;i_buck<e_buckets;i_buck++)
    
    int num_tec_elem = num_tec_quad + num_tec_tri;
    int num_tec_elem2 = 0;
    
    char filename[20];
    sprintf(filename, "tecpl%02d%08d.tec", myid, timeprops->iter);
    FILE *fp = fopen(filename, "w");
    
    //print the tecplot header
    int hrs, mins;
    double secs;
    timeprops->chunktime(&hrs, &mins, &secs);
    
    fprintf(fp, "TITLE= \" %s: time %d:%02d:%g (hrs:min:sec), V*=%g\"\n", mapnames->gis_map.c_str(), hrs, mins, secs, v_star);
    if(elementType == ElementType::TwoPhases)
    {
        fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"PILE_HEIGHT\","
                "\"VOL_FRACT\", \"SOLID_X_MOMENTUM\", \"SOLID_Y_MOMENTUM\","
                "\"FLUID_X_MOMENTUM\", \"FLUID_Y_MOMENTUM\","
                "\"ELEVATION\", \"SOLID_SPEED\", \"FLUID_SPEED\" \n");
    }
    if(elementType == ElementType::SinglePhase)
    {
        fprintf(fp,
                "VARIABLES = \"X\", \"Y\", \"Z\", \"PILE_HEIGHT\", \"X_MOMENTUM\", \"Y_MOMENTUM\", \"ELEVATION\", \"CORRECTSPEED\" \n"); //}
    }
    fprintf(fp, "\n");
    fprintf(fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", num_tec_node, num_tec_elem);
    
    //print the element/BUBBLE-node information
    
    int num_missing_bubble_node = 0;      //for legacy debugging, I (Keith)
    //believe the missing BUBBLE node problem has been solved, it 
    //doesn't seem to be a problem anymore
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if((EmTemp->adapted_flag() != TOBEDELETED) && (abs(EmTemp->adapted_flag()) <= BUFFER))
                num_missing_bubble_node += print_bubble_node(elementType, fp, NodeTable, matprops, EmTemp);
        }  //while(entryp)
    }  //for(i_buck=0;i_buck<e_buckets;i_buck++)
    
    //print the tecplot element connectivity data
    int tec_nodes_in_tec_quad[4];
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            EmArray[0] = EmTemp;
            if((EmTemp->adapted_flag() > TOBEDELETED) && (EmTemp->adapted_flag() <= BUFFER))
            {
                gen = EmTemp->generation();
                
                switch (get_ll_polygon(El_Table, NodeTable, myid, EmArray))
                {
                    case 4:        //we should draw a lower left quad
                        tec_nodes_in_tec_quad[0] = EmArray[2]->ithelem();
                        tec_nodes_in_tec_quad[1] = EmArray[0]->ithelem();
                        tec_nodes_in_tec_quad[2] = EmArray[1]->ithelem();
                        tec_nodes_in_tec_quad[3] = EmArray[3]->ithelem();
                        
                        fprintf(fp, "%d %d %d %d\n", tec_nodes_in_tec_quad[0], tec_nodes_in_tec_quad[1],
                                tec_nodes_in_tec_quad[2], tec_nodes_in_tec_quad[3]);
                        
                        num_tec_elem2++;      //counter of number of drawing elements
                        break;
                    case 3:        //we should draw a lower left triangle
                        tec_nodes_in_tec_quad[0] = EmArray[2]->ithelem();
                        tec_nodes_in_tec_quad[1] = EmArray[0]->ithelem();
                        tec_nodes_in_tec_quad[2] = tec_nodes_in_tec_quad[3] = EmArray[1]->ithelem();
                        
                        fprintf(fp, "%d %d %d %d\n", tec_nodes_in_tec_quad[0], tec_nodes_in_tec_quad[1],
                                tec_nodes_in_tec_quad[2], tec_nodes_in_tec_quad[3]);
                        
                        num_tec_elem2++;      //counter of number of drawing elements
                        break;
                    case 0:        //we should NOT draw a lower left polygon
                        break;
                    default:       //something weird happened report this to the "FBI"
                        assert(0);
                }               //switch(get_ll_polygon(El_Table, myid, EmArray))
                
                /**********************************************************/
                /* match the lower left triangle at the lower left corner */
                /* of a different subdomain with an upper right one       */
                /**********************************************************/

                if(get_ur_tri(El_Table, NodeTable, myid, EmArray))
                {
                    tec_nodes_in_tec_quad[0] = EmArray[1]->ithelem();
                    tec_nodes_in_tec_quad[1] = EmArray[2]->ithelem();
                    tec_nodes_in_tec_quad[2] = tec_nodes_in_tec_quad[3] = EmArray[0]->ithelem();
                    
                    fprintf(fp, "%d %d %d %d\n", tec_nodes_in_tec_quad[0], tec_nodes_in_tec_quad[1],
                            tec_nodes_in_tec_quad[2], tec_nodes_in_tec_quad[3]);
                    
                    num_tec_elem2++;      //counter of number of drawing elements
                }               //if(get_ur_tri(El_Table,myid,EmArray))
                
                /**********************************************************/
                /* add triangles for elements with more refined neighbors */
                /* (aka split neighbors) this fills in the holes          */
                /**********************************************************/
                tec_nodes_in_tec_quad[0] = EmTemp->ithelem();
                for(i_neigh = 0; i_neigh < 4; i_neigh++)
                    if((EmTemp->neigh_gen(i_neigh) > gen) && (EmTemp->neigh_proc(i_neigh) != INIT))
                    {
                        /* draw triangles as quads with the last corner of
                         the triangle repeated.  The three corners are
                         1) this element
                         2) the primary neighbor (one of sides 0->3) 
                         3) the secondary neighbor (the primary +4 side) */

                        EmArray[1] = EmArray[2] = EmArray[3] = NULL;
                        
                        //The Primary Neighbor
                        EmArray[1] = (Element *) El_Table->lookup(EmTemp->neighbor(i_neigh));
                        assert(EmArray[1]);
                        tec_nodes_in_tec_quad[1] = EmArray[1]->ithelem();
                        
                        //The Secondary Neighbor
                        EmArray[2] = (Element *) El_Table->lookup(EmTemp->neighbor(i_neigh + 4));
                        assert(EmArray[2]);
                        tec_nodes_in_tec_quad[2] = tec_nodes_in_tec_quad[3] = EmArray[2]->ithelem();
                        
                        fprintf(fp, "%d %d %d %d\n", tec_nodes_in_tec_quad[0], tec_nodes_in_tec_quad[1],
                                tec_nodes_in_tec_quad[2], tec_nodes_in_tec_quad[3]);
                        
                        num_tec_elem2++;    //counter of number of drawing elements
                    }  //if((neigh_gen[i_neigh]>gen)&&(neigh_proc[i_neigh]!=INIT))
            }  //if((EmTemp->get_adapted_flag()>TOBEDELETED)&& ...
        }  //while(entryp)
    }  //for(i_buck=0;i_buck<e_buckets;i_buck++)
    fflush(fp);
    fclose(fp);
    assert(num_tec_elem2 == num_tec_elem);
    
    MPI_Barrier (MPI_COMM_WORLD);
    
    return;
}

/*
 ***********************************************************
 ***********************************************************
 *********** get orientation of an element *****************
 ***********************************************************
 ***********************************************************
 */
void get_elem_orient(Element * EmTemp, int *xm, int *xp, int *ym, int *yp)
{
    //this does the same thing as Andy's switch statement and is shorter
    *xp = EmTemp->positive_x_side();
    *xm = (2 + *xp) % 4;
    *yp = (1 + *xp) % 4;
    *ym = (3 + *xp) % 4;
    return;
}

/*
 ***********************************************************
 ***********************************************************
 ************* get elevation for an element's **************
 ********************** bubble node ************************
 ***********************************************************
 ***********************************************************
 */

int get_elem_elev(NodeHashTable * NodeTable, MatProps * matprops, Element * EmTemp, double &elevation)
{
    int j;
    double resolution, xcoord, ycoord;
    Node *NdTemp;
    
    NdTemp = (Node *) NodeTable->lookup(EmTemp->key());
    
    if(NdTemp == NULL)
    {
        j = Get_max_resolution(&resolution);
        if(j != 0)
        {
            printf("error in Get_max_resolution\n");
            exit(1);
        }
        
        xcoord = EmTemp->coord(0) * (matprops->LENGTH_SCALE);
        ycoord = EmTemp->coord(1) * (matprops->LENGTH_SCALE);
        j = Get_elevation(resolution, xcoord, ycoord, elevation);
        if(j != 0)
        {
            printf("error in Get_elevation\n");
            exit(1);
        }
        return (1);
    }
    else
    {
        elevation = NdTemp->elevation() * (matprops->LENGTH_SCALE);
        return (0);
    }
}

/* Keith Wrote this 20061120
 ************************************************************************
 * get lower, left, and lower-left neighbors: store them in EmArray[1..3]
 * returns 4 if should draw a lower left quad
 * returns 3 if should draw a lower left triangle
 * returns 0 if shouldn't draw a lower left polygon
 ************************************************************************
 */
int get_ll_polygon(ElementsHashTable * El_Table, NodeHashTable * NodeTable, int myid, Element * EmArray[4])
{
    
    int xm, xp, ym, yp;           //this element's  left right down up 
    int zmxm, zmxp, zmym, zmyp;   //left/down neighbor's left right down up
    int gen = EmArray[0]->generation();
    
    
    int yada = 0;
    
    for(int ielem = 1; ielem < 4; ielem++)
        EmArray[ielem] = NULL;
    
    //yada++;
    /*
     if((*(EmArray[0]->pass_key()+0)==905969664)&&
     (*(EmArray[0]->pass_key()+1)==0))
     printf("stop me\n");
     */
    //yada++;
    //determine which ways are left (xm) and down (ym)
    get_elem_orient(EmArray[0], &xm, &xp, &ym, &yp);
    
    int neigh_proc_xm = EmArray[0]->neigh_proc(xm);
    int neigh_proc_xm_p_4 = EmArray[0]->neigh_proc(xm+4);
    int neigh_proc_ym = EmArray[0]->neigh_proc(ym);
    
    int neigh_gen_xm=EmArray[0]->neigh_gen(xm);
    int neigh_gen_ym=EmArray[0]->neigh_gen(ym);
    
    
    if(((neigh_proc_xm == INIT) || (neigh_proc_ym == INIT)) || (((neigh_gen_xm < gen) || (neigh_gen_ym < gen))
            && (EmArray[0]->which_son() != 0)))
        return 0;
    
    //yada++;
    
    EmArray[1] = (Element *) El_Table->lookup(EmArray[0]->neighbor(xm + 4));
    EmArray[2] = (Element *) El_Table->lookup(EmArray[0]->neighbor(ym));
    
    if(!(((neigh_proc_xm_p_4 == myid) || ((neigh_proc_xm_p_4 == -2) && (neigh_proc_xm == myid))) || (neigh_proc_ym
            == myid)))
        //can't get to lower left neighbor by taking either route
        //because lower and left neighbors are both ghost cells
        //so draw a lower left triangle instead of a lower left quad
        return 3;
    
    if((neigh_proc_xm_p_4 == myid) || ((neigh_proc_xm_p_4 == -2) && (neigh_proc_xm == myid)))
    {
        //can get to lower left neighbor by going left then down
        assert(EmArray[1]);
        get_elem_orient(EmArray[1], &zmxm, &zmxp, &zmym, &zmyp);
        if(EmArray[1]->neigh_proc(zmym) == -1)
        {
            EmArray[1] = EmArray[2] = EmArray[3] = NULL;
            return 0;
        }
        EmArray[3] = (Element *) El_Table->lookup(EmArray[1]->neighbor(zmym + 4));
        //printf("left then down\n");
        //yada++;
        
    }
    else
    {                           //neigh_proc_ym==myid
        //printf("down then left\n");
        //can get to lower left neighbor by going down then left
        assert(EmArray[2]);
        get_elem_orient(EmArray[2], &zmxm, &zmxp, &zmym, &zmyp);
        if(EmArray[2]->neigh_proc(zmxm) == -1)
        {
            EmArray[1] = EmArray[2] = EmArray[3] = NULL;
            return 0;
        }
        EmArray[3] = (Element *) El_Table->lookup(EmArray[2]->neighbor(zmxm));
    }
    
    if(!(EmArray[3]))
    {
        printf("myid=%d xm=%d xp=%d ym=%d yp=%d zmxm=%d zmxp=%d zmym=%d zmyp=%d EmArray[0]->key={\n", myid, xm,
               xp, ym, yp, zmxm, zmxp, zmym, zmyp);
        cout<<EmArray[0]->key()<<"}\n";
        scanf("%d", &yada);
    }
    assert(EmArray[3]);          //make this full proof
    return 4;
}

/* Keith Wrote this 20061120
 ************************************************************************
 * get right and upper neighbors: store them in EmArray[1..2]
 * returns 3 if should draw an upper right triangle
 * returns 0 otherwise
 ************************************************************************
 */

int get_ur_tri(ElementsHashTable * El_Table, NodeHashTable * NodeTable, int myid, Element * EmArray[4])
{
    int xm, xp, ym, yp;           //this  element's  left right down up 
    int xpxm, xpxp, xpym, xpyp;   //right neighbor's left right down up
    int ypxm, ypxp, ypym, ypyp;   //up    neighbor's left right down up
    int gen;
    
    assert(EmArray[0]);
    EmArray[1] = EmArray[2] = EmArray[3] = NULL;
    
    //determine which ways are left (xm) and down (ym)
    get_elem_orient(EmArray[0], &xm, &xp, &ym, &yp);
    gen = EmArray[0]->generation();
    int xpfirst = xp;
    if(EmArray[0]->neigh_gen(xp) > gen)
        xp += 4;                    //account for more refined neighbors
    
    
    //make sure it's not a domain boundary
    if((EmArray[0]->neigh_proc(xpfirst) != INIT) && (EmArray[0]->neigh_proc(yp) != INIT))
    {
        
        EmArray[2] = (Element *) El_Table->lookup(EmArray[0]->neighbor(yp));
        if(!EmArray[2])
        {
            ElemBackgroundCheck(El_Table, NodeTable, EmArray[0]->key(),stdout);
            ElemBackgroundCheck(El_Table, NodeTable, EmArray[0]->neighbor(yp),stdout);
            
        }
        assert(EmArray[2]);
        
        get_elem_orient(EmArray[2], &ypxm, &ypxp, &ypym, &ypyp);
        
        EmArray[1] = (Element *) El_Table->lookup(EmArray[0]->neighbor(xp));
        if(!EmArray[1])
        {
            ElemBackgroundCheck(El_Table, NodeTable, EmArray[0]->key(),stdout);
            ElemBackgroundCheck(El_Table, NodeTable, EmArray[0]->neighbor(xp),stdout);
        }
        assert(EmArray[1]);
        get_elem_orient(EmArray[1], &xpxm, &xpxp, &xpym, &xpyp);
        //account for more refined neighbor's neighbors
        if(EmArray[1]->neigh_gen(xpyp) > EmArray[0]->neigh_gen(xp))
            xpyp += 4;
        
        //check if it's a lower left corner of another subdomain
        if((EmArray[2]->neigh_proc(ypxp) != INIT) && (EmArray[2]->neigh_proc(ypxp) != EmArray[0]->neigh_proc(xp))
           && (EmArray[2]->neigh_proc(ypxp) != EmArray[0]->neigh_proc(yp)) && (EmArray[2]->neigh_proc(ypxp) == EmArray[1]->neigh_proc(xpyp)))
            return (3);
    }
    
    EmArray[1] = EmArray[2] = EmArray[3] = NULL;
    
    return (0);
}

/*
 ***********************************************************
 ***********************************************************
 ************ print the element's bubble node **************
 ***********************************************************
 ***********************************************************
 */
//for use when writing ascii tecplot files
int print_bubble_node(ElementType elementType,FILE * fp, NodeHashTable * NodeTable, MatProps * matprops, Element * EmTemp)
{
    int num_missing_bubble_node;
    double velocity_scale, momentum_scale, elevation;
    
    velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities
    momentum_scale = matprops->HEIGHT_SCALE * velocity_scale;     // scaling factor for the momentums
            
    num_missing_bubble_node = get_elem_elev(NodeTable, matprops, EmTemp, elevation);

    if(elementType == ElementType::TwoPhases)
    {
        double Vel[4];
        double volf;
        if(EmTemp->state_vars(0) > GEOFLOW_TINY)
            volf = (EmTemp->state_vars(1) / (EmTemp->state_vars(0)));
        else
            volf = 0;
        EmTemp->eval_velocity(0.0, 0.0, Vel);
        fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e\n", (EmTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                (EmTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                elevation + (EmTemp->state_vars(0)) * (matprops)->HEIGHT_SCALE,
                EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, volf,
                EmTemp->state_vars(2) * momentum_scale, EmTemp->state_vars(3) * momentum_scale,
                EmTemp->state_vars(4) * momentum_scale, EmTemp->state_vars(5) * momentum_scale,
                elevation, sqrt(Vel[0] * Vel[0] + Vel[1] * Vel[1]) * velocity_scale,
                sqrt(Vel[2] * Vel[2] + Vel[3] * Vel[3]) * velocity_scale);
    }
    if(elementType == ElementType::SinglePhase)
    {
        double VxVy[2];
        EmTemp->eval_velocity(0.0, 0.0, VxVy);
        fprintf(fp, "%e %e %e %e %e %e %e %e\n", (EmTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                (EmTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                (elevation + (EmTemp->state_vars(0)) * (matprops)->HEIGHT_SCALE),
                EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, EmTemp->state_vars(1) * momentum_scale,
                EmTemp->state_vars(2) * momentum_scale, elevation,
                sqrt(VxVy[0] * VxVy[0] + VxVy[1] * VxVy[1]) * velocity_scale);
    }
    return (num_missing_bubble_node);
}

/*
 ***********************************************************
 ***********************************************************
 ******** write a binary string the way tecplot does *******
 ***********************************************************
 ***********************************************************
 */

/* the name DumpString is what tecplot uses in 
 tecplot9/util/preplot/preplot.c
 see useful_lib.C for fwriteI() */
void DumpString(FILE * fp, char *str)
{
    int i = 0;
    while (str[i] != '\0')
        fwriteI(fp, (int) str[i++]);
    fwriteI(fp, 0);
    return;
}

/*
 ***********************************************************
 ***********************************************************
 *********** geoflow visualization output ******************
 ***********************************************************
 ***********************************************************
 */

void viz_output(ElementType elementType,ElementsHashTable * El_Table, NodeHashTable * NodeTable, int myid, int numprocs, MatProps * matprops,
                TimeProps * timeprops, MapNames * mapnames)
{
    int i, k;
    int element_counter = 0;
    Element *EmTemp;
    char filename[18] = "viz_outputxxx.out";
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    FILE *fp;
    
    // global descriptor -- list of file names
    if(myid == 0 && timeprops->ifstart())
    {
        
        char gl_filename[19] = "viz_filenames.out";
        fp = fopen(gl_filename, "w");
        fprintf(fp, "%d    !number of files\n", numprocs);
        
        for(i = 0; i < numprocs; i++)
        {
            filename[10] = (i % 1000) / 100 + 48;
            filename[11] = (i % 100) / 10 + 48;
            filename[12] = (i % 10) + 48;
            
            fprintf(fp, "%s\n", filename);
        }
        fprintf(fp, "\n");
        fclose(fp);
    }
    
    filename[10] = (myid % 1000) / 100 + 48;
    filename[11] = (myid % 100) / 10 + 48;
    filename[12] = (myid % 10) + 48;
    double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities
            
    // actual data output
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
                element_counter++;
        }
    }
    //information below will probably be changed later
    int points = element_counter;
    int variables = 7;
    int inputs = 12;
    int parameters = 2;
    
    if(timeprops->ifstart())
        fp = fopen(filename, "w");
    else
        fp = fopen(filename, "a+");
    fprintf(fp, "%f %d %d %d %d\n", timeprops->cur_time, parameters, points, variables, inputs);
    fprintf(fp,
            "\"key unsigned[%d]\", \"coordinate double[%d]\", \"terrain_slope double[%d]\", \"pile_height double[1]\", \"pile_velocity double[%d]\", \"processor int[1]\", \"FOM double[1]\" \n",
            2, 3, 2, 2);
    // output parameters
    fprintf(fp, "time_step int[1] %d \n", timeprops->iter);
    fprintf(fp, "number_of_processors int[1] %d \n", numprocs);
    
    //output field data
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
            {
                double height = EmTemp->state_vars(0);
                Node *NdTemp = (Node *) NodeTable->lookup(EmTemp->key());
                if(height < GEOFLOW_TINY)
                {
                    double zero = 0;
                    fprintf_sfc_key(fp, EmTemp->key());
                    fprintf(fp, "%f %f %f %f %f %f %f %f %d %f \n",
                            (EmTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                            (EmTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                            (NdTemp->elevation()) * (matprops)->LENGTH_SCALE, EmTemp->zeta(0),
                            EmTemp->zeta(1), zero, zero, zero, myid, zero);
                }
                else
                {
                    if(elementType == ElementType::TwoPhases)
                    {
                        double Vel[4];
                        EmTemp->eval_velocity(0.0, 0.0, Vel);
                        fprintf_sfc_key(fp, EmTemp->key());
                        fprintf(fp, " %f %f %f %f %f %f %f %f %d %f \n", (EmTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (EmTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                (NdTemp->elevation()) * (matprops)->LENGTH_SCALE, EmTemp->zeta(0),
                                EmTemp->zeta(1), height * (matprops)->HEIGHT_SCALE,
                                velocity_scale * Vel[0], //*EmTemp->state_vars(1)/EmTemp->state_vars(0),
                                velocity_scale * Vel[1], //EmTemp->state_vars(2)/EmTemp->state_vars(0),
                                myid, (double) (height / 20.));  //EmTemp->state_vars(0)/20 is just a filler item
                    }
                    if(elementType == ElementType::SinglePhase)
                    {
                        double VxVy[2];
                        EmTemp->eval_velocity(0.0, 0.0, VxVy);
                        fprintf_sfc_key(fp, EmTemp->key());
                        fprintf(fp, "%f %f %f %f %f %f %f %f %d %f \n",
                                (EmTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (EmTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                (NdTemp->elevation()) * (matprops)->LENGTH_SCALE, EmTemp->zeta(0),
                                EmTemp->zeta(1), height * (matprops)->HEIGHT_SCALE,
                                velocity_scale * VxVy[0],     //*EmTemp->state_vars(1)/EmTemp->state_vars(0),
                                velocity_scale * VxVy[1],    //EmTemp->state_vars(2)/EmTemp->state_vars(0),
                                myid, (double) (height / 20.));       //EmTemp->state_vars(0)/20 is just a filler item
                    }
                }
            }
        }
    }
    
    fprintf(fp, "\n");
    fclose(fp);
    
    return;
}

/*************************************
 **************************************
 **************************************
 *  output mesh only
 **************************************
 **************************************
 *************************************/
void meshplotter(ElementsHashTable * El_Table, NodeHashTable * NodeTable, MatProps * matprops, TimeProps * timeprops,
                 MapNames * mapnames, double v_star)
{
    int myid, i;
    int numprocs;
    int material;
    int done = 1;
    int TECTAG = 123;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int element_counter = 0;
    Element *EmTemp;
    Node *NodeTemp;
    char filename[256];
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 1.0\n");
        fflush(stdout);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    sprintf(filename, "mshpl%02d%08d.tec", myid, timeprops->iter);
    
    int order;
    
    double momentum_scale = matprops->HEIGHT_SCALE * sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums
            
    FILE *fp;
    
    //  if ( myid == 0 )
    // {
    fp = fopen(filename, "w");
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 2.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    int hours, minutes;
    double seconds;
    timeprops->chunktime(&hours, &minutes, &seconds);
    
    fprintf(fp, "TITLE= \" %s (MESH OUTPUT) time %d:%02d:%g (hrs:min:sec), V*=%g\"\n", mapnames->gis_map.c_str(), hours,
            minutes, seconds, v_star);
    
    fprintf(fp,
            "VARIABLES = \"X\" \"Y\" \"NODAL_ELEVATION\" \"PROC\" \"PILE_HEIGHT\" \"XMOMENTUM\" \"YMOMENTUM\" \"KEY0\" \"KEY1\" \"GENERATION\" \"SON\" \"ELM_ELEVATION\" \"XSLOPE\" \"YSLOPE\" \"XCURVATURE\" \"YCURVATURE\" \"ELMLOC1\" \"ELMLOC2\"");
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 3.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
                element_counter++;
        }
    }
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 4.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    fprintf(fp, "\n");
    fprintf(fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", element_counter * 4, element_counter);
    
    int elements = El_Table->get_no_of_buckets();
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 4.1\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
            {
                order = 1;
                for(int k = 0; k < 5; k++)
                {
                    int help_order = EmTemp->order(k);
                    if(help_order > order)
                        order = help_order;
                }
                
                double err = sqrt(EmTemp->el_error(0));
                for(int j = 0; j < 4; j++)
                {
                    NodeTemp = (Node *) NodeTable->lookup(EmTemp->node_key(j));
                    assert(NodeTemp);
                    //int* dof = NodeTemp->getdof();
                    int jj = j;
                    if(1)
                    {           //NodeTemp->getinfo() != S_C_CON) {
                        unsigned tmpkey[KEYLENGTH];
                        SET_OLDKEY(tmpkey,EmTemp->key());
#ifdef TWO_PHASES 
                        fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n",
                                (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                NodeTemp->elevation() * (matprops->LENGTH_SCALE), myid,
                                EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, EmTemp->state_vars(2) * momentum_scale,
                                EmTemp->state_vars(3) * momentum_scale, tmpkey[0],tmpkey[1],
                                EmTemp->generation(), EmTemp->which_son(),
                                EmTemp->elevation() * (matprops->LENGTH_SCALE), EmTemp->zeta(0),
                                EmTemp->zeta(1), EmTemp->curvature(0) / (matprops->LENGTH_SCALE),
                                EmTemp->curvature(1) / (matprops->LENGTH_SCALE), EmTemp->elm_loc(0),
                                EmTemp->elm_loc(1));
#else
                        fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n",
                                (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                NodeTemp->elevation() * (matprops->LENGTH_SCALE), myid,
                                EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, EmTemp->state_vars(1) * momentum_scale,
                                EmTemp->state_vars(2) * momentum_scale, tmpkey[0],tmpkey[1],
                                EmTemp->generation(), EmTemp->which_son(),
                                EmTemp->elevation() * (matprops->LENGTH_SCALE), EmTemp->zeta(0),
                                EmTemp->zeta(1), EmTemp->curvature(0) / (matprops->LENGTH_SCALE),
                                EmTemp->curvature(1) / (matprops->LENGTH_SCALE), EmTemp->elm_loc(0),
                                EmTemp->elm_loc(1));
#endif
                    }
                    else
                    { // S_C_CON will have a discontinuity in the elevation so fix that by interpolation
                        double elev;
                        int neighside, mynode;
                        if(j == EmTemp->which_son() + 1)
                        {
                            mynode = j - 1;
                            neighside = j;
                        }
                        else if(j == EmTemp->which_son() - 1)
                        {
                            mynode = j + 1;
                            neighside = j - 1;
                            if(neighside < 0)
                                neighside = 3;
                        }
                        else if(EmTemp->which_son() == 0)
                        {
                            mynode = 0;
                            neighside = 2;
                        }
                        else if(EmTemp->which_son() == 3)
                        {
                            mynode = 3;
                            neighside = 0;
                        }
                        else
                        {
                            mynode = 1;
                            neighside = 1;
                            //assert(0);
                        }
                        Node *NodeTemp2 = (Node *) NodeTable->lookup(EmTemp->node_key(mynode));
                        assert(NodeTemp2);
                        
                        elev = .5 * NodeTemp2->elevation();
                        Element *EmTemp2 = (Element *) El_Table->lookup(
                                (EmTemp->neighbor(neighside)));
                        assert(EmTemp2);
                        
                        NodeTemp2 = (Node *) NodeTable->lookup(EmTemp2->node_key(j));
                        elev += .5 * NodeTemp2->elevation();

                        unsigned tmpkey[KEYLENGTH];
                        SET_OLDKEY(tmpkey,EmTemp->key());
#ifdef TWO_PHASES 
                        fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n",
                                (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                elev * (matprops->LENGTH_SCALE), myid, EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE,
                                EmTemp->state_vars(2) * momentum_scale, EmTemp->state_vars(3) * momentum_scale, tmpkey[0],tmpkey[1], EmTemp->generation(), EmTemp->which_son(),
                                EmTemp->elevation() * (matprops->LENGTH_SCALE), EmTemp->zeta(0),
                                EmTemp->zeta(1), EmTemp->curvature(0) / (matprops->LENGTH_SCALE),
                                EmTemp->curvature(1) / (matprops->LENGTH_SCALE), EmTemp->elm_loc(0),
                                EmTemp->elm_loc(1));
#else
                        fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d %e %e %e %e %e %d %d\n",
                                (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                                (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                                elev * (matprops->LENGTH_SCALE), myid, EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE,
                                EmTemp->state_vars(1) * momentum_scale, EmTemp->state_vars(2) * momentum_scale, tmpkey[0],tmpkey[1], EmTemp->generation(), EmTemp->which_son(),
                                EmTemp->elevation() * (matprops->LENGTH_SCALE), EmTemp->zeta(0),
                                EmTemp->zeta(1), EmTemp->curvature(0) / (matprops->LENGTH_SCALE),
                                EmTemp->curvature(1) / (matprops->LENGTH_SCALE), EmTemp->elm_loc(0),
                                EmTemp->elm_loc(1));
#endif
                    }
                    
                }
            }
        }
    }
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 5.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    fprintf(fp, "\n");
    
    for(i = 0; i < element_counter; i++)
    {
        for(int j = 0; j < 4; j++)
            fprintf(fp, "%d ", i * 4 + j + 1);
        fprintf(fp, "\n");
    }
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 6.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    fclose(fp);
    
    if(myid == TARGETPROCB)
    {
        printf("at meshplotter 7.0\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}

/*************************************
 **************************************
 **************************************
 *  output mesh for geoflow viz people
 **************************************
 **************************************
 *************************************/
void vizplotter(ElementsHashTable * El_Table, NodeHashTable * NodeTable, MatProps * matprops, TimeProps * timeprops)
{
    int myid, i;
    int numprocs;
    int material;
    int done = 1;
    int TECTAG = 123;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int element_counter = 0;
    Element *EmTemp;
    Node *NodeTemp;
    SFC_Key *nodes;
    char filename[20];            //="vizplotxxxxxxxx.plt";
    sprintf(filename, "vizplot%08d.plt", timeprops->iter);
    
    /*filename[14] = (which % 10) + 48;
     filename[13] = (which % 100)/10 + 48;
     filename[12] = (which % 1000)/100 + 48;
     filename[11] = (which % 10000)/1000 + 48;
     filename[10] = (which % 100000)/10000 + 48;
     filename[9] = (which % 1000000)/100000 + 48;
     filename[8] = (which % 10000000)/1000000 + 48;
     filename[7] = (which % 100000000)/10000000 + 48; */
    int order;
    double momentum_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums
            
    FILE *fp;
    
    if(myid == 0)
    {
        fp = fopen(filename, "w");
        fprintf(fp, "TITLE= \"MESH OUTPUT\"\n");
        fprintf(fp,
                "VARIABLES = \"X\", \"Y\", \"Z\", \"PROC\" \"PILE_HEIGHT\" \"XMOMENTUM\" \"YMOMENTUM\" \"KEY0\" \"KEY1\" \"GENERATION\" \"SON\"");
    }
    else
    {
        MPI_Recv(&done, 1, MPI_INT, myid - 1, TECTAG, MPI_COMM_WORLD, &status);
        //outData.open(filename, ios::app);
        fp = fopen(filename, "a+");
    }
    
    int no_of_buckets = El_Table->get_no_of_buckets();
    vector<HashEntryLine> &bucket=El_Table->bucket;
    tivector<Element> &elenode_=El_Table->elenode_;
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
                element_counter++;
        }
    }
    
    //outData<<'\n';
    fprintf(fp, "\n");
    
    //outData<<"ZONE N="<<element_counter*4<<", E="<<element_counter<<", F=FEPOINT, ET=QUADRILATERAL"<<'\n';
    fprintf(fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n", element_counter * 4, element_counter);
    
    //@ElementsBucketDoubleLoop
    for(int ibuck = 0; ibuck < no_of_buckets; ibuck++)
    {
        for(int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++)
        {
            EmTemp = &(elenode_[bucket[ibuck].ndx[ielm]]);
            if(EmTemp->adapted_flag() > 0)
            {
                order = 1;
                for(int k = 0; k < 5; k++)
                {
                    int help_order = EmTemp->order(k);
                    if(help_order > order)
                        order = help_order;
                }
                
                double err = sqrt(EmTemp->el_error(0));
                for(int j = 0; j < 4; j++)
                {
                    NodeTemp = (Node *) NodeTable->lookup(EmTemp->node_key(j));
                    unsigned tmpkey[KEYLENGTH];
                    SET_OLDKEY(tmpkey,EmTemp->key());
                    //int* dof = NodeTemp->getdof();
#ifdef TWO_PHASES 
                    fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d\n",
                            (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                            (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                            (NodeTemp->elevation()) * (matprops)->LENGTH_SCALE, myid,
                            EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, EmTemp->state_vars(2) * momentum_scale,
                            EmTemp->state_vars(3) * momentum_scale, tmpkey[0],tmpkey[1],
                            EmTemp->generation(), EmTemp->which_son());
#else
                    fprintf(fp, "%e %e %e %d %e %e %e %u %u %d %d\n",
                            (NodeTemp->coord(0)) * (matprops)->LENGTH_SCALE,
                            (NodeTemp->coord(1)) * (matprops)->LENGTH_SCALE,
                            (NodeTemp->elevation()) * (matprops)->LENGTH_SCALE, myid,
                            EmTemp->state_vars(0) * (matprops)->HEIGHT_SCALE, EmTemp->state_vars(1) * momentum_scale,
                            EmTemp->state_vars(2) * momentum_scale, tmpkey[0],tmpkey[1],
                            EmTemp->generation(), EmTemp->which_son());
#endif
                    
                }
            }
        }
    }
    
    //outData<<'\n'; 
    fprintf(fp, "\n");
    
    for(i = 0; i < element_counter; i++)
    {
        for(int j = 0; j < 4; j++)
            //outData<<i*4+j+1<<' '; 
            fprintf(fp, "%d ", i * 4 + j + 1);
        
        //outData<<'\n'; 
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    
    if(myid != numprocs - 1)
        MPI_Send(&done, 1, MPI_INT, myid + 1, TECTAG, MPI_COMM_WORLD);
    
}
