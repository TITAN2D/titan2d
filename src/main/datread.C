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
 * $Id: datread.C 233 2012-03-27 18:30:40Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

#ifdef SUNOS
extern "C" void initial_(int*, double*, double*);
#endif

#ifdef IBMSP
extern "C" void initial(int*, double*, double*);
#endif

#ifdef CRAY
extern "C" void INITIAL(int*, double*, double*);
#endif


//this reads in the funky grid, ignoring the material properties at the 
//end of the file, those are now read from frict.data in Read_data()
void Read_grid(int myid, int numprocs, HashTable** NodeTable, ElementsHashTable** ElemTable, MatProps* matprops_ptr,
               OutLine* outline_ptr)
{
    int Node_Num, Elem_Num;
    
    int NODE_TABLE_SIZE = 400000;
    
    //char  filename[14] = "lsh4800xx.inp";
    char filename[14] = "funkyxxxx.inp";
    //unsigned min_key[KEYLENGTH];
    //unsigned max_key[KEYLENGTH];
    double doublekeyrange[2];
    double XRange[2];
    double YRange[2];
    unsigned key[KEYLENGTH];
    double coord[DIMENSION];
    double height;
    Node* NodeP;
    int i, j, k;
    
    // read in nodal data
    sprintf(filename, "funky%04d.inp", myid);
    FILE* fp;
    
    fp = fopen_bin(filename, "r");
    if(!fp)
    {
        printf("Can't open file for %d \n", myid);
        exit(0);
    }
    
    int version, DoublesFromFloats;
    freadI(fp, &version);
    
    switch (version)
    {
        case 20061109:
            DoublesFromFloats = 1;
            break;
        case 20061110:
            DoublesFromFloats = 0;
            break;
        default:
            printf("Read_data() does not recognize binary funkyxxxx.inp version %d\n", version);
            exit(1);
            break;
    }
    
    freadI(fp, &Node_Num);
    
    if(DoublesFromFloats)
    {
        for(i = 0; i < KEYLENGTH; i++)
            freadF2D(fp, &(doublekeyrange[i]));
        
        freadF2D(fp, &(XRange[0]));  //min x
        freadF2D(fp, &(XRange[1]));  //max x
        freadF2D(fp, &(YRange[0]));  //min y
        freadF2D(fp, &(YRange[1]));
    } //max y
    else
    {
        for(i = 0; i < KEYLENGTH; i++)
            freadD(fp, &(doublekeyrange[i]));
        
        freadD(fp, &(XRange[0]));  //min x
        freadD(fp, &(XRange[1]));  //max x
        freadD(fp, &(YRange[0]));  //min y
        freadD(fp, &(YRange[1]));
    } //max y
    
    double xminmax[2], yminmax[2];
    for(i = 0; i < 2; i++)
    {
        XRange[i] = XRange[i] / matprops_ptr->LENGTH_SCALE;
        xminmax[i] = XRange[i];
    }
    for(i = 0; i < 2; i++)
    {
        YRange[i] = YRange[i] / matprops_ptr->LENGTH_SCALE;
        yminmax[i] = YRange[i];
    }
    
    *NodeTable = new HashTable(doublekeyrange, NODE_TABLE_SIZE, XRange, YRange);
    
    for(i = 0; i < Node_Num; i++)
    {
        for(j = 0; j < KEYLENGTH; j++)
            freadU(fp, &(key[j]));
        
        if(DoublesFromFloats)
            for(j = 0; j < DIMENSION; j++)
                freadF2D(fp, &(coord[j]));
        else
            for(j = 0; j < DIMENSION; j++)
                freadD(fp, &(coord[j]));
        
        for(j = 0; j < 2; j++)
            coord[j] = coord[j] / matprops_ptr->LENGTH_SCALE;
        NodeP = new Node(sfc_key_from_oldkey(key), coord, matprops_ptr);
        (*NodeTable)->add(sfc_key_from_oldkey(key), NodeP);
    }
    (*NodeTable)->print0();
    //done reading in node data
    //start reading in element data
    
    int EL_TABLE_SIZE = 100000;
    
    //char  filename[14] = "lsh4800xx.inp";
    int material, elm_loc[2];
    unsigned opposite_brother[2];
    
    Element* Quad9P;
    float* value = new float[2];
    BC* baddress[4];
    void* p;
    int* assocp;/*--*/
    unsigned* keyP;
    
    unsigned nodes_old[9][2];
    unsigned neigh_old[4][2];
    SFC_Key nodes[9];
    SFC_Key neigh[9];
    int neighbor_proc[4];
    
    int temp2;
    int interflag;
    
    freadI(fp, &Elem_Num);  //--number of the elements assigned to the proc
           
    *ElemTable = new ElementsHashTable(doublekeyrange, EL_TABLE_SIZE, XRange, YRange, *NodeTable);
    for(int ielem = 0; ielem < Elem_Num; ielem++)
    {
        
        for(j = 0; j < 9; j++)
        {
            for(k = 0; k < KEYLENGTH; k++)
                freadU(fp, &(nodes_old[j][k]));
            SET_NEWKEY(nodes[j],nodes_old[j]);
        }
        
        interflag = 0;  //---switch for interface
        for(j = 0; j < 4; j++)
        {
            
            freadI(fp, &(neighbor_proc[j]));  //--read the neighbor info
                   
            if(neighbor_proc[j] != -1)  //--if there is neighbor(-1 means the edge is bound)
            {
                if(neighbor_proc[j] != myid) //--the neighbor belongs to other proc
                    interflag = 1; //--switch is used for avoiding nominating neighbor twice
                            
                for(k = 0; k < KEYLENGTH; k++)
                    freadU(fp, &(neigh_old[j][k])); //--read the left parts of the key
            }
            
            else
                //--there is no neighbor 
                for(k = 0; k < KEYLENGTH; k++)
                    neigh_old[j][k] = 0;

            SET_NEWKEY(neigh[j],neigh_old[j]);
        }
        
        BC* bcptr = 0;
        int bcf = 0;
        
        //.....the essential boundary conditions....
        
        for(j = 0; j < 4; j++)
        {
            freadI(fp, &temp2);
            
            if(temp2 != -1) //--there is bound constraint
            {
                if(!bcf)
                    bcptr = new BC();
                bcptr->type[j] = 1; //--intialize type
                        
                /* "value" is a FLOAT so DON'T use freadD when DoublesFromFloats
                 is false (and obviously don't use freadD when it's true 
                 either) */
                for(k = 0; k < 2; k++)
                    freadF(fp, &(bcptr->value[j][0][k])); //--j: edge number
                           
                bcf = 1;
            }
        }
        
        //.....the natural boundary conditions.....
        for(j = 0; j < 4; j++)
        {
            
            freadI(fp, &temp2);
            
            if(temp2 != -1) //--there is bound constraint
            {
                if(!bcf)
                    bcptr = new BC();
                if(bcptr->type[j] == 0)
                    bcptr->type[j] = 2; //--intialize type
                else
                    bcptr->type[j] = 3; //--intialize type
                            
                /* "value" is a FLOAT so DON'T use freadD when DoublesFromFloats
                 is false (and obviously don't use freadD when it's true 
                 either) */
                for(k = 0; k < 2; k++)
                    freadF(fp, &(bcptr->value[j][1][k])); //--j: edge number
                           
                bcf = 1;
            }
        }
        
        freadI(fp, &(elm_loc[0]));
        freadI(fp, &(elm_loc[1]));
        freadU(fp, &(opposite_brother[0]));
        freadU(fp, &(opposite_brother[1]));
        freadI(fp, &material);
        double pile_height = 0.0;
        
        if(!bcf)
            bcptr = NULL; //--this element is not on the bound
        Quad9P = (*ElemTable)->generateElement(nodes, neigh, neighbor_proc, bcptr, material, elm_loc, pile_height, myid,
                                               sfc_key_from_oldkey(opposite_brother));
        (*ElemTable)->add(nodes[8], Quad9P);
        Quad9P->find_positive_x_side(*NodeTable);
        Quad9P->calculate_dx(*NodeTable);
    }
    (*ElemTable)->print0();
    /************************************************************/
    /* need to change this so that smallest cell size is chosen */
    /* based on cube root of volume with failsafe for a minimum */
    /* of 2 cells across smallest pile/flux source axis instead */
    /* of basing it on just the smallest pile/flux source axis  */
    /* which is what we're doing now, will need to print out a  */
    /* warning that the calculation may be slow when the        */
    /* failsafe is activated                                    */
    /************************************************************/

    double dx[2] =
    { *(Quad9P->get_dx() + 0), *(Quad9P->get_dx() + 1) };
    double DX = dx[0];
    if(dx[0] < dx[1])
        DX = dx[1];
    
    REFINE_LEVEL = Quad9P->get_gen()
            + ceil(log(DX * (matprops_ptr->number_of_cells_across_axis) / (matprops_ptr->smallest_axis)) / log(2.0));
    
    //(mdj)2007-04-12 if(REFINE_LEVEL<3) REFINE_LEVEL=3;
    if(REFINE_LEVEL < 0)
        REFINE_LEVEL = 0;
    
    //set the nodal type information
    Element *EmTemp;
    Node *NdTemp;
    int inode;
    int num_buck = (*ElemTable)->get_no_of_buckets();
    HashEntryPtr* buck = (*ElemTable)->getbucketptr();
    for(int i = 0; i < num_buck; i++)
        if(*(buck + i))
        {
            
            HashEntryPtr currentPtr = *(buck + i);
            while (currentPtr)
            {
                
                EmTemp = (Element*) (currentPtr->value);
                currentPtr = currentPtr->next;
                assert(EmTemp);
                
                EmTemp->put_myprocess(myid);
                
                NdTemp = (Node*) (*NodeTable)->lookup(EmTemp->key());
                assert(NdTemp);
                NdTemp->putinfo(BUBBLE);
                
                for(inode = 0; inode < 4; inode++)
                {
                    NdTemp = (Node*) (*NodeTable)->lookup(EmTemp->getNode()[inode]);
                    assert(NdTemp);
                    NdTemp->putinfo(CORNER);
                }
                
                for(inode = 4; inode < 8; inode++)
                {
                    NdTemp = (Node*) (*NodeTable)->lookup(EmTemp->getNode()[inode]);
                    assert(NdTemp);
                    NdTemp->putinfo(SIDE);
                }
            }
        }
    
#ifdef MAX_DEPTH_MAP
    outline_ptr->init(Quad9P->get_dx(), REFINE_LEVEL - Quad9P->get_gen(), xminmax, yminmax);
#endif
    
    delete[] value;
}

