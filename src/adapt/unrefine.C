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
 * $Id: unrefine.C 150 2007-06-27 20:28:42Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"
#include "../header/hadapt.h"
#include "../header/hadapt_inline.hpp"


//#define MIN_GENERATION -1
#define TARGET_PROC -1
/*
int IfMissingElem(HashTable* El_Table, int myid, int iter, int isearch)
{
    
    return (0);
    
    if(iter < 4)
        return (0);
    
    int i, yada = 0;
    HashEntryPtr* buck = El_Table->getbucketptr();
    for(i = 0; i < El_Table->get_no_of_buckets(); i++)
        if(*(buck + i))
        {
            HashEntryPtr currentPtr = *(buck + i);
            while (currentPtr)
            {
                Element* Curr_El = (Element*) currentPtr->value;
                currentPtr = currentPtr->next;
                if((*(Curr_El->pass_key() + 0) == 2110300160) && (*(Curr_El->pass_key() + 1) == 0))
                {
                    printf("search %d found it, adapted=%d\n", isearch, Curr_El->get_adapted_flag());
                    //printf("search %d found it, adapted=%d, enter an int\n",isearch,Curr_El->get_adapted_flag());
                    //scanf("%d",&yada);
                    return (0);
                }
            }
        }
    
    return (1);
}*/
HAdaptUnrefine::HAdaptUnrefine(ElementsHashTable* _ElemTable, NodeHashTable* _NodeTable,TimeProps* _timeprops, MatProps* _matprops):
   EleNodeRef(_ElemTable,_NodeTable),
   ElemProp(_ElemTable->ElemProp),
   matprops_ptr(_matprops),
   timeprops_ptr(_timeprops)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    nodesToDelete.resize(threads_number);
    elementsToDelete.resize(threads_number);
    for(int ithread=0;ithread<threads_number;++ithread)
    {
        nodesToDelete[ithread].resize(0);
        elementsToDelete[ithread].resize(0);
    }

    brothers_to_unrefine_ndx.resize(threads_number);
}

void HAdaptUnrefine::unrefine(const double target)
{
    int time_step = timeprops_ptr->iter;
    PROFILING3_DEFINE(pt_start);
    PROFILING3_START(pt_start);
    
    //convinience references
    tivector<ContentStatus> &status=ElemTable->status_;
    tivector<int> &adapted=ElemTable->adapted_;
    tivector<int> &generation=ElemTable->generation_;
    tivector<double> &el_error=ElemTable->el_error_[0];
    tivector<int> &refined=ElemTable->refined_;
    tivector<int> &which_son=ElemTable->which_son_;

    Element* Curr_El;
    for(int ithread=0;ithread<threads_number;++ithread)
    {
        brothers_to_unrefine_ndx[ithread].resize(0);
    }

    NewFatherList.resize(0);
    OtherProcUpdate.resize(0);
    
    ASSERT2(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine PRE",false)==0);

    //-------------------go through all the elements of the subdomain------------------------
    int no_of_buckets = ElemTable->get_no_of_buckets();
    vector<HashEntryLine> &bucket=ElemTable->bucket;
    tivector<Element> &elenode_=ElemTable->elenode_;

    //@ElementsSingleLoopNoStatusCheck
    for(ti_ndx_t ndx=0;ndx<ElemTable->size();++ndx)
    {
        //don't need to check if element schedule for deletion
        if(adapted[ndx] == NEWFATHER)adapted[ndx] = NOTRECADAPTED;
    }
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_init,pt_start);
    
    // start unrefinement
    find_brothers_to_unrefine__create_new_fathers(target);
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_find,pt_start);
    //updateBrothersIndexes for new fathers, as brothers which are also new fathers was not handled during initiation
    ElemTable->updateBrothersIndexes(true);
    ASSERT3(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine After new fathers creation",false)==0);
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_updateBrothersIndexes,pt_start);

    int iproc;
    time_t tic, toc;
    
    //assert(!IfMissingElem(El_Table, myid, time_step, 0));
    unrefine_neigh_update();
    ASSERT3(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine After unrefine_neigh_update",false)==0);
    //assert(!IfMissingElem(El_Table, myid, time_step, 1));
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_unrefine_neigh_update,pt_start);
    
    unrefine_interp_neigh_update();
    if(numprocs>1)ElemTable->update_neighbours_ndx_on_ghosts(true);
    ASSERT3(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine After unrefine_interp_neigh_update",false)==0);
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_unrefine_interp_neigh_update,pt_start);
    
    delete_oldsons();
    if(numprocs>1)ElemTable->update_neighbours_ndx_on_ghosts(true);
    ASSERT3(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine After delete_oldsons",false)==0);
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_delete_oldsons,pt_start);

    move_data(numprocs, myid, ElemTable, NodeTable, timeprops_ptr);
    if(numprocs>1)ElemTable->update_neighbours_ndx_on_ghosts(true);
    ASSERT2(ElemTable->checkPointersToNeighbours("HAdaptUnrefine::unrefine After move_data",false)==0);
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_update_neighbours_ndx_on_ghosts,pt_start);
    
    ElemProp->calc_wet_dry_orient2();
    PROFILING3_STOPADD_RESTART(HAdaptUnrefine_unrefine_calc_wet_dry_orient,pt_start);
    return;
}

void HAdaptUnrefine::delete_oldsons()
{
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();
        nodesToDelete[ithread].resize(0);
        elementsToDelete[ithread].resize(0);
        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(int iupdate = 0; iupdate < NewFatherList.size(); iupdate++)
        {
            int ison, isonneigh, ineigh, inode;
            Element *EmSon, *EmNeigh;
            Node* NdTemp;
            Element *EmFather = &(ElemTable->elenode_[NewFatherList[iupdate]]);
            ASSERT3(EmFather);
        
            EmFather->set_refined_flag(0);

    #ifdef DEB3
            for(inode = 4; inode < 8; inode++)
            {
                NdTemp = (Node *) NodeTable->lookup(EmFather->node_key(inode));
                ASSERT3(NdTemp);
            }
    #endif
            inode = 8;
            ASSERT3(EmFather->node_bubble_ndx()==NodeTable->lookup_ndx(EmFather->key()));
            NdTemp = &(NodeTable->elenode_[EmFather->node_bubble_ndx()]);
            
            for(ison = 0; ison < 4; ison++)
            {
                EmSon = &(ElemTable->elenode_[EmFather->son_ndx(ison)]);
                ASSERT3(EmFather->son_ndx(ison)==ElemTable->lookup_ndx(EmFather->son(ison)));
                ASSERT3(EmSon->adapted_flag()==OLDSON);
                EmSon->set_adapted_flag(TOBEDELETED);

                //delete son's bubble nodes
                ASSERT3(EmSon->node_bubble_ndx()==NodeTable->lookup_ndx(EmFather->son(ison)));
                nodesToDelete[ithread].push_back(EmSon->node_bubble_ndx());

                //delete son to son edge nodes
                inode = (ison + 1) % 4 + 4;
                ASSERT3(EmSon->node_key_ndx(inode)==NodeTable->lookup_ndx(EmSon->node_key(inode)));
                nodesToDelete[ithread].push_back(EmSon->node_key_ndx(inode));

                //check 2 other edge nodes per son and delete if necessary
                //the first
                isonneigh = ison;
                inode = isonneigh + 4;
                ineigh = ison;

                //EmNeigh=(Element *) El_Table->lookup(EmFather->neighbor[ineigh]);
                if((EmFather->neigh_gen(ineigh) == EmFather->generation()) || (EmFather->neigh_proc(ineigh) == -1))
                {
                    if(ti_ndx_not_negative(EmSon->node_key_ndx(inode)) && NodeTable->status_[EmSon->node_key_ndx(inode)]>=0)
                    {
                        //delete if not deleted previously
                        ASSERT3(EmSon->node_key_ndx(inode)==NodeTable->lookup_ndx(EmSon->node_key(inode)));
                        nodesToDelete[ithread].push_back(EmSon->node_key_ndx(inode));
                    }
                }
                else if(EmFather->neigh_proc(ineigh) == myid)
                {
                    ASSERT3(EmFather->neigh_gen(ineigh) == EmFather->generation() + 1);
                    ASSERT3(EmSon->node_key_ndx(inode)==NodeTable->lookup_ndx(EmSon->node_key(inode)));
                    NodeTable->info_[EmSon->node_key_ndx(inode)]=S_S_CON;
                }

                //the second
                isonneigh = (ison + 3) % 4;
                inode = isonneigh + 4;
                ineigh = inode;

                //EmNeigh=(Element *) El_Table->lookup(EmFather->neighbor[ineigh]);
                if((EmFather->neigh_gen(ineigh) == EmFather->generation()) || (EmFather->neigh_proc(ineigh % 4) == -1))
                {
                    if(ti_ndx_not_negative(EmSon->node_key_ndx(inode)) && NodeTable->status_[EmSon->node_key_ndx(inode)]>=0)
                    {
                        //delete if not deleted previously
                        ASSERT3(EmSon->node_key_ndx(inode)==NodeTable->lookup_ndx(EmSon->node_key(inode)));
                        nodesToDelete[ithread].push_back(EmSon->node_key_ndx(inode));
                    }
                    ASSERT3(EmFather->node_key_ndx(inode)==NodeTable->lookup_ndx(EmFather->node_key(inode)));
                    NodeTable->info_[EmFather->node_key_ndx(inode)]=SIDE;
                }
                else if(EmFather->neigh_proc(ineigh) == myid)
                {
                    ASSERT3(EmSon->node_key_ndx(inode)==NodeTable->lookup_ndx(EmSon->node_key(inode)));
                    NodeTable->info_[EmSon->node_key_ndx(inode)]=S_S_CON;

                    ASSERT3(EmFather->node_key_ndx(inode)==NodeTable->lookup_ndx(EmFather->node_key(inode)));
                    NodeTable->info_[EmFather->node_key_ndx(inode)]=S_C_CON;
                }

                //Now delete this oldson Element
                elementsToDelete[ithread].push_back(EmSon->ndx());
                EmFather->son_ndx(ison,ti_ndx_doesnt_exist);
            }
            inode = 8;
            
            ASSERT3(EmFather->node_bubble_ndx()==NodeTable->lookup_ndx(EmFather->key()));
            NodeTable->info_[EmFather->node_bubble_ndx()]=BUBBLE;
        }
    }
    merge_vectors_from_threads_to0_omp(nodesToDelete);
    merge_vectors_from_threads_to0_omp(elementsToDelete);

    NodeTable->removeNodes(&(nodesToDelete[0][0]),nodesToDelete[0].size());
    ElemTable->removeElements(&(elementsToDelete[0][0]),elementsToDelete[0].size());

    return;
}
void HAdaptUnrefine::unrefine_neigh_update()
{
    int iupdate, ineigh, isonA, isonB, ineighme, ikey;
    Element *EmNeigh, *EmFather;
    ti_ndx_t EmNeighNdx, EmFatherNdx;
    
    //loop through the NEWFATHER elements
    for(iupdate = 0; iupdate < NewFatherList.size(); iupdate++)
    {
        
        //I'm a NEWFATHER I'm going to update my neighbors with
        //my information and if he's a NEWFATHER too I'm going
        //to and update my information about him
        EmFatherNdx=NewFatherList[iupdate];
        EmFather = &(ElemTable->elenode_[EmFatherNdx]);
        ASSERT2(EmFather); //Help I've been abducted call the FBI!!!

        for(ineigh = 0; ineigh < 8; ineigh++)
            if(EmFather->neigh_proc(ineigh) == myid)
            {
                //only update the information of on processor neighbors in
                //this function.
                EmNeighNdx = EmFather->neighbor_ndx(ineigh);//ElemTable->lookup_ndx(EmFather->neighbor(ineigh));
                EmNeigh = &(ElemTable->elenode_[EmNeighNdx]);
                ASSERT2(EmNeighNdx == ElemTable->lookup_ndx(EmFather->neighbor(ineigh))); //Somebody has abducted my neighbor call the FBI!!!
                
                if(EmNeigh->adapted_flag() != NEWFATHER)
                {
                    //If I knew my neighbor was a NEWFATHER that means
                    //my neighbor already updated his and my neighbor
                    //information about each other
                    
                    if(EmNeigh->adapted_flag() == OLDSON)
                    {
                        //I am introduced to a NEWFATHER neighbor by his OLDSON
                        ASSERT2(EmNeigh->father_ndx() == ElemTable->lookup_ndx(EmNeigh->father_by_ref())); //Somebody has abducted my neighbor call the FBI!!!*/
                        EmNeighNdx = EmNeigh->father_ndx();//ElemTable->lookup_ndx(EmNeigh->father_by_ref());
                        EmNeigh = &(ElemTable->elenode_[EmNeighNdx]);
                    }
                    
                    //One of my OLDSONs will introduce my neighbor to me
                    isonA = ineigh % 4;
                    isonB = (isonA + 1) % 4;
                    
                    for(ineighme = 0; ineighme < 4; ineighme++)
                        if((EmNeigh->neighbor(ineighme)==EmFather->son(isonA)) || (
                                EmNeigh->neighbor(ineighme)==EmFather->son(isonB)))
                            break;
                    if(!(ineighme < 4))
                        printf("DANGER\n");
                    assert(ineighme < 4);
                    
                    //Give my neighbor my contact information
                    EmNeigh->get_neigh_gen(ineighme, EmFather->generation());
                    EmNeigh->get_neigh_gen(ineighme + 4, EmFather->generation());
                    
                    EmNeigh->set_neigh_proc(ineighme + 4, -2);
                    
                    if(EmNeigh->adapted_flag() == NEWFATHER)
                    {
                        //if my neighbor is a NEWFATHER, I need to update my
                        //information about him too
                        EmFather->get_neigh_gen(ineigh, EmNeigh->generation());
                        EmFather->get_neigh_gen(ineigh + 4, EmNeigh->generation());
                        
                        EmFather->set_neigh_proc(ineigh + 4, -2);
                        
                        //EmNeigh is inside this if() so only one loop through
                        //KEYLENGTH is made, it saves on overhead
                            EmFather->set_neighbor(ineigh, EmNeigh->key());
                            EmFather->set_neighbor(ineigh + 4, EmNeigh->key());
                            EmNeigh->set_neighbor(ineighme, EmFather->key());
                            EmNeigh->set_neighbor(ineighme + 4, EmFather->key());

                            EmFather->neighbor_ndx(ineigh, EmNeigh->ndx());
                            EmFather->neighbor_ndx(ineigh + 4, EmNeigh->ndx());
                            EmNeigh->neighbor_ndx(ineighme, EmFather->ndx());
                            EmNeigh->neighbor_ndx(ineighme + 4, EmFather->ndx());
                    }
                    else
                    {
                        //EmNeigh is inside this if() so only one loop through
                        //KEYLENGTH is made, it saves on overhead
                        EmNeigh->set_neighbor(ineighme, EmFather->key());
                        EmNeigh->set_neighbor(ineighme + 4, EmFather->key());

                        EmNeigh->neighbor_ndx(ineighme, EmFather->ndx());
                        EmNeigh->neighbor_ndx(ineighme + 4, EmFather->ndx());
                    }
                }
            }
    }
    
    return;
}
void HAdaptUnrefine::find_brothers_to_unrefine__create_new_fathers(double target)
{
    #pragma omp parallel
    {
        int ithread=omp_get_thread_num();
        nodesToDelete[ithread].resize(0);
        elementsToDelete[ithread].resize(0);
        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t ndx = 0; ndx < ElemTable->size(); ++ndx)
        {
            if(status_[ndx]>=0 && adapted_[ndx]==NOTRECADAPTED)
            {
                //if this is a refined element don't involve!!!

                // if this if the original element, don't unrefine.  only son 0 checks for unrefinement!
                if((which_son_[ndx] == 0) && (generation_[ndx] > MIN_GENERATION))
                {
                    find_brothers_to_unrefine(ndx, target, ithread);
                }
            }
        }
        //convinience reference to 0th element
        vector< array<ti_ndx_t,4> > &m_brothers_to_unrefine_ndx=brothers_to_unrefine_ndx[0];
        //!allocate new farthers
        if(ithread==0)
        {
            merge_vectors_from_threads_to0(brothers_to_unrefine_ndx);

            ti_ndx_t N=m_brothers_to_unrefine_ndx.size();
            for(ti_ndx_t i=0;i<N;++i)
            {
                // we want to unrefine this element...
                //first we create the father element
                ti_ndx_t new_father = ElemTable->generateAddElement_ndx(node_key_[0][m_brothers_to_unrefine_ndx[i][2]]);
                ASSERT2(ti_ndx_not_negative(new_father));// a copy of the parent should always be on the same process as the sons
                NewFatherList.push_back(new_father);
            }
        }
        #pragma omp barrier

        //!init new farthers
        ti_ndx_t N=m_brothers_to_unrefine_ndx.size();
        #pragma omp for schedule(dynamic,TITAN2D_DINAMIC_CHUNK)
        for(ti_ndx_t i=0;i<N;++i)
        {
            // we want to unrefine this element...
            //first we create the father element
            ti_ndx_t new_father=NewFatherList[i];
            ElemTable->elenode_[new_father].init(m_brothers_to_unrefine_ndx[i].data(), NodeTable, ElemTable, matprops_ptr);
        }
        #pragma omp barrier
        if(ithread==0)
        {
            for(ti_ndx_t i=0;i<N;++i)
            {
                ti_ndx_t new_father=NewFatherList[i];
                for(int ineigh = 0; ineigh < 8; ineigh++)
                {
                    if((ElemTable->neigh_proc_[ineigh][new_father] >= 0) && (ElemTable->neigh_proc_[ineigh][new_father] != myid))
                    {
                        OtherProcUpdate.push_back(new_father);
                        break;
                    }
                }
            }
        }
    }
}

void HAdaptUnrefine::find_brothers_to_unrefine(ti_ndx_t ndx,double target, int ithread)
{
    ASSERT2(which_son_[ndx]==0);

    int j;
    int unrefine_flag = 1;
    ti_ndx_t new_father;
    std::array<ti_ndx_t,4> bros_ndx;
    /*if(opposite_brother_flag_[ndx] == 0)
    {
        elements_[ndx].find_opposite_brother(ElemTable);
        if(opposite_brother_flag_[ndx] == 0)
            return 0;
    }*/
    int i = 0;
    while (i < 4 && unrefine_flag == 1)
    {
        ti_ndx_t bro_ndx;
#ifdef DEB3
        //extra checking
        bro_ndx=ElemTable->lookup_ndx(brothers_[i][ndx]);
        if(i==0)
        {
            assert(bro_ndx==ndx);
        }
        else if(i==1)
        {
            if(neighbor_ndx_[1][ndx]==neighbor_ndx_[5][ndx])
                assert(bro_ndx==neighbor_ndx_[1][ndx]);
        }
        else if(i==2)
        {
            if(neighbor_ndx_[2][neighbor_ndx_[1][ndx]]==neighbor_ndx_[6][neighbor_ndx_[1][ndx]])
                assert(bro_ndx==neighbor_ndx_[2][neighbor_ndx_[1][ndx]]);
        }
        else if(i==3)
        {
            if(neighbor_ndx_[2][ndx]==neighbor_ndx_[6][ndx])
                assert(bro_ndx==neighbor_ndx_[2][ndx]);
        }
#endif
        if(i==0)
        {
            //0th brother is itself
            bro_ndx=ndx;
        }
        else if(i==1)
        {
            //i.e. neighboring element is same generation
            if(neighbor_ndx_[1][ndx]==neighbor_ndx_[5][ndx])
                bro_ndx=neighbor_ndx_[1][ndx];
            else
                return;
        }
        else if(i==2)
        {
            //i.e. neighboring element is same generation
            if(neighbor_ndx_[2][neighbor_ndx_[1][ndx]]==neighbor_ndx_[6][neighbor_ndx_[1][ndx]])
                bro_ndx=neighbor_ndx_[2][neighbor_ndx_[1][ndx]];
            else
                return;
        }
        else if(i==3)
        {
            //i.e. neighboring element is same generation
            if(neighbor_ndx_[2][ndx]==neighbor_ndx_[6][ndx])
                bro_ndx=neighbor_ndx_[2][ndx];
            else
                return;
        }

        if(ti_ndx_negative(bro_ndx)) //|| EmTemp->refined != 0)
            return;

        if(adapted_[bro_ndx] != NOTRECADAPTED) //this should be sufficient
            return;
        bros_ndx[i] = bro_ndx;
        if(myprocess_[bro_ndx] != myid)
            return; //should not be necessary because of "adapted" check

        unrefine_flag = check_unrefinement(bro_ndx, target);
        i++;
    }

    if(unrefine_flag == 1)
    {
        brothers_to_unrefine_ndx[ithread].push_back(bros_ndx);
    }
    return;
}

//make this a Node and Element friend fucntion
void HAdaptUnrefine::unrefine_interp_neigh_update()
{
#ifdef USE_MPI
    if(numprocs < 2)
        return;
    
    int *num_send = CAllocI1(numprocs), *num_recv = CAllocI1(numprocs), *isend = CAllocI1(numprocs);
    int ierr, iopu, iproc, ineigh, ison, ikey, neigh_proc; //iopu stands for i other processor update
    int send_tag = 061116; //2006 November 16, the day I coded this function.
    MPI_Request* request = new MPI_Request[2 * numprocs];
    Element* EmFather;
    Node *NdTemp;
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 1.0\n", myid);
    fflush(stdout);
    
    for(iproc = 0; iproc < numprocs; iproc++)
        num_send[iproc] = num_recv[iproc] = isend[iproc] = 0;
    
    for(iopu = 0; iopu < OtherProcUpdate.size(); iopu++)
    {
        EmFather = &(ElemTable->elenode_[OtherProcUpdate[iopu]]);
        assert(EmFather);
        assert(EmFather->myprocess() == myid);
        for(ineigh = 0; ineigh < 8; ineigh++)
        {
            neigh_proc = EmFather->neigh_proc(ineigh);
            if((neigh_proc >= 0) && (neigh_proc != myid))
                num_send[neigh_proc]++;
        }
    }
    num_send[myid] = 0; //better paranoid then dead;
            
    //send to each processor the number of his elements he has to update because of me
    //receive from each processor the number of elements I have to update because of him
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, MPI_COMM_WORLD);
    
    if(myid == TARGET_PROC)
        for(iproc = 0; iproc < numprocs; iproc++)
        {
            printf("myid=%2d to/from proc %2d num_send=%6d num_recv=%6d\n", myid, iproc, num_send[iproc],
                   num_recv[iproc]);
            fflush(stdout);
        }
    MPI_Barrier (MPI_COMM_WORLD);
    if((myid != TARGET_PROC) && (TARGET_PROC < numprocs) && (TARGET_PROC > 0))
    {
        printf("myid=%2d to/from proc %2d num_send=%6d num_recv=%6d\n", myid, TARGET_PROC, num_send[TARGET_PROC],
               num_recv[TARGET_PROC]);
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*
     for(iproc=0;iproc<nump;iproc++)
     printf("myid=%d neighproc=%d numsend=%d numrecv=%d\n",
     myid,iproc,num_send[iproc],num_recv[iproc]);
     exit(0);
     */
    num_recv[myid] = 0; //better paranoid then dead;
            
    int max_num_send = 0, max_num_recv = 0;
    
    for(iproc = 0; iproc < numprocs; iproc++)
    {
        if(num_send[iproc] > max_num_send)
            max_num_send = num_send[iproc];
        if(num_recv[iproc] > max_num_recv)
            max_num_recv = num_recv[iproc];
    }
    
    unsigned **send, **recv;
    
    if(max_num_send > 0)
        send = CAllocU2(numprocs, 4 * KEYLENGTH * max_num_send);
    if(max_num_recv > 0)
        recv = CAllocU2(numprocs, 4 * KEYLENGTH * max_num_recv);
    
    for(iproc = 0; iproc < numprocs; iproc++)
        if((iproc != myid) && (num_recv[iproc] > 0))
            ierr = MPI_Irecv((void *) recv[iproc], 4 * KEYLENGTH * num_recv[iproc], MPI_UNSIGNED, iproc,
                             send_tag + iproc, MPI_COMM_WORLD, (request + numprocs + iproc));
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 2.0\n", myid);
    fflush(stdout);
    
    for(iopu = 0; iopu < OtherProcUpdate.size(); iopu++)
    {
        EmFather = &(ElemTable->elenode_[OtherProcUpdate[iopu]]);
        assert(EmFather);
        for(ineigh = 0; ineigh < 8; ineigh++)
        {
            neigh_proc = EmFather->neigh_proc(ineigh);
            if((neigh_proc >= 0) && (neigh_proc != myid))
            {
                
                switch (ineigh)
                {
                    case 0:
                    case 7:
                        ison = 0;
                        break;
                    case 1:
                    case 4:
                        ison = 1;
                        break;
                    case 2:
                    case 5:
                        ison = 2;
                        break;
                    case 3:
                    case 6:
                        ison = 3;
                        break;
                    default:
                        assert(0);
                }
                
                NdTemp = (Node*) NodeTable->lookup(EmFather->node_key(ineigh % 4 + 4));
                assert(NdTemp);
                
                if(EmFather->neigh_gen(ineigh) - 1 == EmFather->generation())
                    NdTemp->info(S_C_CON);
                
                //The element I want my neighbor to update
                SET_OLDKEY((&(send[neigh_proc][(4 * isend[neigh_proc] + 0) * KEYLENGTH])), EmFather->neighbor(ineigh));
                //the OLDSON who will introduce his NEWFATHER to his
                //neighbor on another processor
                SET_OLDKEY((&(send[neigh_proc][(4 * isend[neigh_proc] + 1) * KEYLENGTH])), EmFather->son(ison));
                //the NEWFATHER element
                SET_OLDKEY((&(send[neigh_proc][(4 * isend[neigh_proc] + 2) * KEYLENGTH])), EmFather->key());
                SET_OLDKEY((&(send[neigh_proc][(4 * isend[neigh_proc] + 3) * KEYLENGTH])), NdTemp->key());
                
                isend[neigh_proc]++;
            }
        }
    }
    for(iproc = 0; iproc < numprocs; iproc++)
        assert(isend[iproc] == num_send[iproc]);
    
    CDeAllocI1(isend);
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 3.0\n", myid);
    fflush(stdout);
    
    for(iproc = 0; iproc < numprocs; iproc++)
        if((iproc != myid) && (num_send[iproc] > 0))
            ierr = MPI_Isend((void *) send[iproc], 4 * KEYLENGTH * num_send[iproc], MPI_UNSIGNED, iproc,
                             send_tag + myid, MPI_COMM_WORLD, (request + iproc));
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 4.0\n", myid);
    fflush(stdout);
    
    int NumProcsNotRecvd, ifrecvd;
    MPI_Status status;
    Element *EmTemp;
    if(max_num_recv > 0)
    {
        do
        {
            
            for(iproc = 0; iproc < numprocs; iproc++)
                if((iproc != myid) && (num_recv[iproc] > 0))
                {
                    if(myid == TARGET_PROC)
                        printf("myid=%d iproc=%2d unref_interp_neigh_update 5.0\n", myid, iproc);
                    fflush(stdout);
                    
                    //only check processors I haven't already handled
                    ifrecvd = 0;
                    //printf("myid=%d before Test",myid); fflush(stdout);
                    MPI_Test(request + numprocs + iproc, &ifrecvd, &status);
                    //printf("myid=%d after Test",myid); fflush(stdout);
                    if(ifrecvd)
                    {
                        if(myid == TARGET_PROC)
                            printf("myid=%d iproc=%2d unref_interp_neigh_update 6.0\n", myid, iproc);
                        fflush(stdout);
                        
                        //I have just received new data from a neighboring processor
                        //I need to update some elements on my interprocessor
                        //boundary

                        for(iopu = 0; iopu < num_recv[iproc]; iopu++)
                        {
                            //one by one check my Element's that my neighbor processor
                            //says I need to update
                            
                            //Hi I'm EmTemp
                            EmTemp = (Element*) ElemTable->lookup(sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 0) * KEYLENGTH])));
                            
                            //my old neighbor will introduce his NEWFATHER to me
                            for(ineigh = 0; ineigh < 8; ineigh++)
                                if(EmTemp->neighbor(ineigh)==sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 1) * KEYLENGTH])))
                                    break;
                            assert(ineigh < 8); //I don't know this Element pretending to be
                            //my neighbor, I'm calling the FBI to report a suspicious
                            //person... err.. suspicious Element
                            int ineighmod4 = ineigh % 4;
                            
                            if(EmTemp->adapted_flag() == OLDSON)
                            {
                                //I'm moving out too so I'll introduce my NEWFATHER to my
                                //neighbor's NEWFATHER who my neighbor just introduced to me
                                
                                //we can use my father's key directly instead of having to
                                //use get_father() because we know that my father's key
                                //was assigned to me in unrefine_elements()
                                EmFather = (Element*) ElemTable->lookup(EmTemp->father_by_ref());
                                assert(EmFather);
                                
                                //which means my father will only have 1 neighbor on that side
                                EmFather->set_neigh_proc(ineighmod4 + 4, -2);
                                
                                EmFather->set_neighbor(ineighmod4, sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 2) * KEYLENGTH])));
                                EmFather->set_neighbor(ineighmod4 + 4, sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 2) * KEYLENGTH])));
                                
                                //EmFather->neighbor_ndx(ineighmod4, ElemTable->lookup_ndx(EmFather->neighbor(ineighmod4)));
                                //EmFather->neighbor_ndx(ineighmod4 + 4, ElemTable->lookup_ndx(EmFather->neighbor(ineighmod4 + 4)));

                                //I know my neighbor on the other processor is the same
                                //generation as me because we were both unrefined and only
                                //one of us would have been able to unrefine if we were of
                                //different generations, that means our NEWFATHERs are the
                                //same generation as each other too
                                EmFather->get_neigh_gen(ineighmod4, EmFather->generation());
                                EmFather->get_neigh_gen(ineighmod4 + 4, EmFather->generation());
                                
                                //all the other neighbor information remains the same but
                                //must delete 2 nodes I'll do 1 now and the other one was
                                //just done previous or will be deleted next
                                NdTemp = (Node *) NodeTable->lookup(EmTemp->node_key(ineighmod4 + 4));
                                if(NdTemp)
                                {
                                    NodeTable->removeNode(NdTemp);
                                }
                                NdTemp = (Node *) NodeTable->lookup(EmFather->node_key(ineighmod4 + 4));
                                NdTemp->info(SIDE);
                                
                            }
                            else if(EmTemp->adapted_flag() >= NOTRECADAPTED)
                            {
                                //my neighbor on the other processor was of either my
                                //generation or one generation younger/higher than me,
                                //that makes
                                //his father either one generation lower than me or of
                                //my generation. Either way his FATHER will be the only
                                //neighbor I have on that side
                                
                                EmTemp->set_neigh_proc(ineighmod4 + 4, -2);
                                
                                EmTemp->set_neighbor(ineighmod4 + 4, sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 2) * KEYLENGTH])));
                                EmTemp->set_neighbor(ineighmod4, sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 2) * KEYLENGTH])));
                                
                                //EmFather->neighbor_ndx(ineighmod4 + 4, ElemTable->lookup_ndx(EmFather->neighbor(ineighmod4 + 4)));
                                //EmFather->neighbor_ndx(ineighmod4, ElemTable->lookup_ndx(EmFather->neighbor(ineighmod4)));

                                EmTemp->get_neigh_gen(ineighmod4 + 4, EmTemp->neigh_gen(ineighmod4) - 1);
                                EmTemp->get_neigh_gen(ineighmod4, EmTemp->neigh_gen(ineighmod4) - 1);
                                //(recv[iproc][(4*iopu+3)*KEYLENGTH+0]?-1:1)*
                                //(recv[iproc][(4*iopu+3)*KEYLENGTH+1]);
                                
                                //all the other neighbor information remains the same but may need
                                //to delete some nodes and update side node
                                //titan doesn't use node order but for a continuous galerkin code
                                //you may need to reset the node order here, which would be extra
                                //information to be communicated accross.  If you enforce the
                                //minimum generation to be zero then you only need 1 unsigned
                                //number to hold the generation or you could use a union of
                                //int and unsigned an the first bit of the unsigned would be the
                                //sign of the int so you wouldn't need 2 unsigneds to communicate
                                //the generation.
                                NdTemp = (Node *) NodeTable->lookup(EmTemp->node_key(ineighmod4 + 4));
                                assert(NdTemp);
                                if(EmTemp->neigh_gen(ineighmod4) == EmTemp->generation())
                                    NdTemp->info(SIDE);
                                else
                                {	      //my neighbor was my generation but his father moved in
                                    NdTemp->info(S_S_CON);
                                    
                                    NdTemp = (Node*) NodeTable->lookup(sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 3) * KEYLENGTH])));
                                    if(!NdTemp)
                                    {
                                        ElemBackgroundCheck(ElemTable, NodeTable,
                                                            sfc_key_from_oldkey(&(recv[iproc][(4 * iopu + 3) * KEYLENGTH])), stdout);
                                        assert(NdTemp);
                                    }
                                    NdTemp->info(S_C_CON);
                                }
                            }
                            else
                                //EmTemp is missing... he may have been abducted
                                //Call the FBI to report a missing person! err... make that
                                //Call the FBI to report a missing Element!
                                assert(0);
                            
                        }
                        num_recv[iproc] = 0;
                        if(myid == TARGET_PROC)
                            printf("myid=%d iproc=%2d unref_interp_neigh_update 7.0\n", myid, iproc);
                        fflush(stdout);
                        
                    }
                }
            if(myid == TARGET_PROC)
                printf("myid=%d unref_interp_neigh_update 8.0\n", myid);
            fflush(stdout);
            
            NumProcsNotRecvd = 0;
            for(iproc = 0; iproc < numprocs; iproc++)
                if((iproc != myid) && (num_recv[iproc] > 0))
                    NumProcsNotRecvd++;
            
        }
        while (NumProcsNotRecvd > 0);
    }
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 9.0\n", myid);
    fflush(stdout);
    
    if(max_num_recv > 0)
        CDeAllocU2(recv);
    CDeAllocI1(num_recv);
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 10.0\n", myid);
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(max_num_send > 0)
        CDeAllocU2(send);
    CDeAllocI1(num_send);
    delete[] request;
    
    if(myid == TARGET_PROC)
        printf("myid=%d unref_interp_neigh_update 11.0\n", myid);
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
#endif //USE_MPI
}
