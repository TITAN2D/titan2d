/*
 * =====================================================================================
 *
 *       Filename:  xdfm_write.cc
 *
 *    Description:  XDMF writer generates paraview readable output
 *                  1. writes titan output to HDF5 file 
 *                  2. writes metadata to XML file
 *                  TODO 1: Need to re-refine UNREFINED cells 
 *
 *        Created:  05/15/2007 12:17:29 PM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <hdf5.h>
#include <mpi.h>

#include <hashtab.h>
#include <element2.h>
#include <node.h>
#include <GisApi.h>
#include <hd5calls.h>


//just an innocent declaration. defined at end of the file
double interpolate_elv(HashTable *,HashTable *,Element *,int ,unsigned *);
void xdmf_fopen (char *);
void xdmf_fclose (ofstream &);

int write_xdmf(HashTable *El_Table, HashTable *NodeTable,TimeProps *timeprops_ptr,
      MatProps *matprops_ptr, MapNames *mapnames, const int mode)
{
   int i,j,k;
   double *state_vars,*coord,elevation;
   
   // scaling factor for the momentums
   double momentum_scale=matprops_ptr->HEIGHT_SCALE*
      sqrt(matprops_ptr->LENGTH_SCALE*(matprops_ptr->GRAVITY_SCALE)); 
   unsigned *nodes;
   static double time_prev;

   /* generate XML file if required */
   int myid; 
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   ofstream xmlf;
   char filename[20];
   sprintf(filename,"xdmf%03d.xmf",myid);
   if (mode == XDMF_NEW)
      xdmf_fopen (filename);
   xmlf.open (filename, ios::app);

   /* save time-value to compare at end */
   if (mode != XDMF_CLOSE) 
     time_prev = timeprops_ptr->timesec();

   /* don't rewrite the last time-step */
   /* just close and return */
   if ((mode == XDMF_CLOSE) && 
       (timeprops_ptr->timesec() - time_prev < 0.0001))
   {
      xdmf_fclose (xmlf);
      xmlf.close ();
      return 0;
   }

   Element *EmTemp=NULL;
   Node *NodeTemp=NULL;
   HashEntry *entryptr;
   // reset connection-ids
   int buckets = NodeTable->get_no_of_buckets();
   for (i = 0; i < buckets; i++)
   {
     entryptr = *(NodeTable->getbucketptr () + i);
     while (entryptr)
     {
        NodeTemp = (Node *) entryptr->value;
        NodeTemp->put_con_id (-1);
        entryptr = entryptr->next;
      }
   }

   /* Add connection information to each node */
   int num_elem = 0;
   int num_node = 0;
   int id=0;
   buckets = El_Table->get_no_of_buckets();
   for (i = 0; i < buckets; i++)
   {
      entryptr = *(El_Table->getbucketptr()+i);
      while (entryptr)
      {
         EmTemp = (Element *) entryptr->value;
         if (!EmTemp->get_refined_flag ())
         {
            num_elem++;
            nodes = EmTemp->getNode();
            for (j = 0; j < 4; j++)
            {
               NodeTemp = (Node *) NodeTable->lookup(nodes + j * KEYLENGTH);
               if ( NodeTemp->get_con_id () < 0 )
               {
                  NodeTemp->put_con_id(id++);
                  num_node++;
               }
            }
         }
         entryptr = entryptr->next;
      }
   }

   /* allocate memory for connections and element verticies */
   double * xyz = new double [num_node * 3];
   int * connections = new int [num_elem * 4];

   /* and other output variables */
   double * pheight = new double [num_elem];
   double * xmom = new double [num_elem];
   double * ymom = new double [num_elem];

   /* scan HashTable and store coordinates and variables in vectors */
   int ielm = 0;
   for (i=0; i<buckets; i++)
   {
      entryptr=*(El_Table->getbucketptr()+i);
      while (entryptr)
      {
	 EmTemp= (Element *) entryptr->value;
	 if (!(EmTemp->get_refined_flag()))
	 {
	    state_vars=EmTemp->get_state_vars();
	    pheight[ielm] = state_vars[0]*matprops_ptr->HEIGHT_SCALE;
	    xmom[ielm] = state_vars[1] * momentum_scale;
	    ymom[ielm] = state_vars[2] * momentum_scale;
	    nodes = EmTemp->getNode ();
	    for (j = 0; j < 4; j++)
	    {
	       NodeTemp=(Node *) NodeTable->lookup(nodes+j*KEYLENGTH);
               int inode = NodeTemp->get_con_id();
	       coord = NodeTemp->get_coord();
	       xyz[3*inode] = coord[0] * matprops_ptr->LENGTH_SCALE;
	       xyz[3*inode + 1] = coord[1] * matprops_ptr->LENGTH_SCALE;
	       connections[4*ielm + j] = inode;
		  // Side-Corner node will have discontinuity in elevation
	       if (NodeTemp->getinfo() != S_C_CON)
	       { 
	           elevation = NodeTemp->get_elevation()*matprops_ptr->LENGTH_SCALE;
	           xyz[3*inode + 2] = elevation;
	       }
	       // hence <get a fake one> interpolate from neighbor
	       else
	       {
	         elevation = interpolate_elv(El_Table, NodeTable, EmTemp, j, nodes) *
		               matprops_ptr->LENGTH_SCALE;
		 xyz[3*inode + 2] = elevation;
	       }
	    }
            ielm++;
	 }
	 entryptr=entryptr->next;
      }
   }
  
   /* 
    * Write HDF5 File 
    */
   char hdf5file[20];
   sprintf(hdf5file,"xdmf%02d%08d.h5",myid,timeprops_ptr->iter);
   hid_t h5fid=GH5_fopen(hdf5file,'n');
   
   // Write Mesh data
   GH5_write_mesh_data (h5fid, num_elem, num_node, connections, xyz);


   /*  FREE ALLOCATED MEMORY */
   // xyz memory
   delete []xyz;
   // conns memory
   delete []connections;

   /* Write Variable Data */
   //Pile Height
   GH5_write_state_vars (h5fid, num_elem, pheight, "PILE_HEIGHT");

   //X-Momentum
   GH5_write_state_vars (h5fid, num_elem, xmom, "XMOMENTUM");

   //Y-Momentum
   GH5_write_state_vars (h5fid, num_elem, ymom, "YMOMENTUM");

   GH5_fclose (h5fid);

   // free-memory
   delete [] pheight;
   delete [] xmom;
   delete [] ymom;

   /*xmlf <<  add comments about when the data was generated */
   xmlf << "<Grid Name=\""<<mapnames->gis_map<<"\">"<<endl;
   xmlf << "<Time Value=\""<<timeprops_ptr->timesec()<<"\" />"<<endl;
   // connectivity data
   xmlf << "<Topology Type=\"QUADRILATERAL\" Dimensions=\""<<
     num_elem<<"\" Order=\"0 3 2 1\">"<<endl;
   xmlf << "<DataItem Name=\"Connections\" DataType=\"Int\" Precision=\"4\""<<endl;
   xmlf << "Dimensions=\""<<num_elem<<" 4\" Format=\"HDF\">"<<endl;
   xmlf << "\t\t"<<hdf5file<<":/Mesh/Connections"<<endl;
   xmlf << "</DataItem>"<<endl;
   xmlf << "</Topology>"<<endl;
   // grid-points
   xmlf << "<Geometry Type=\"XYZ\">"<<endl;
   xmlf << "<DataItem Name=\"Coordinates\" DataType=\"Float\" Precision=\"8\""<<endl;
   xmlf << "Dimensions=\""<<num_node<<" 3\" Format=\"HDF\">"<<endl;
   xmlf << "\t\t"<<hdf5file<<":/Mesh/Points"<<endl;
   xmlf << "</DataItem>"<<endl;
   xmlf << "</Geometry>"<<endl;
   // pile-height
   xmlf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Pile Height\">"<<endl;
   xmlf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
   xmlf << "Dimensions=\""<<num_elem<<" 1\" Format=\"HDF\">"<<endl;
   xmlf << "\t\t"<<hdf5file<<":/Properties/PILE_HEIGHT"<<endl;
   xmlf << "</DataItem>"<<endl;
   xmlf << "</Attribute>"<<endl;
   // x-momentum
   xmlf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"X_Momentum\">"<<endl;
   xmlf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
   xmlf << "Dimensions=\""<<num_elem<<" 1\" Format=\"HDF\">"<<endl;
   xmlf << "\t\t"<<hdf5file<<":/Properties/XMOMENTUM"<<endl;
   xmlf << "</DataItem>"<<endl;
   xmlf << "</Attribute>"<<endl;
   // y-momentum
   xmlf << "<Attribute Type=\"Scalar\" Center=\"Cell\" Name=\"Y_Momentum\">"<<endl;
   xmlf << "<DataItem DataType=\"Float\" Precision=\"8\" ";
   xmlf << "Dimensions=\""<<num_elem<<" 1\" Format=\"HDF\">"<<endl;
   xmlf << "\t\t"<<hdf5file<<":/Properties/YMOMENTUM"<<endl;
   xmlf << "</DataItem>"<<endl;
   xmlf << "</Attribute>"<<endl;
   xmlf << "</Grid>"<<endl;
   
   if ( mode == XDMF_CLOSE )
     xdmf_fclose (xmlf);

   xmlf.close();  // close till next output is written

   return 0;
}

/* open an XML file for Metadata and write header information */
void xdmf_fopen (char *filename)
{
  
   ofstream xmf(filename, ios::out);
   if ( !xmf.is_open() )
   {
      cout<<"ERROR: unable to open file for writing,"<<endl<<
	    "check disc space or check permissions"<<endl;
      exit(1);
   }
   xmf <<"<?xml version=\"1.0\" ?>"<< endl;
   xmf <<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]> "<<endl;
   xmf <<"<Xdmf Version=\"2.0\">"<< endl;
   xmf <<"<Domain>"<< endl;
   xmf <<"<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">"<< endl;
   xmf.close ();
   return;
}

/* close XML header file */
void xdmf_fclose (ofstream & xmlf)
{
   xmlf <<"</Grid>"<<endl;
   xmlf <<"</Domain>"<<endl;
   xmlf <<"</Xdmf>"<<endl;
}

/*
 * interpolate_elv is exactly the same code that is
 * being used in meshplotter(...), in tecplot.C
 */
double interpolate_elv(HashTable *El_Table,HashTable *NodeTable,Element *EmTemp,int j, unsigned *nodes)
{
   double elev;
   int neighside, mynode;
   if(j==EmTemp->get_which_son()+1) 
   {
      mynode = j - 1;
      neighside = j;
   }
   else if(j==EmTemp->get_which_son()-1) 
   {
      mynode = j+1;
      neighside = j-1;
      if(neighside < 0)
	 neighside = 3;
   }
   else if(EmTemp->get_which_son()==0) 
   {
     mynode = 0;
     neighside = 2;
   }
   else if(EmTemp->get_which_son()==3) 
   {
      mynode = 3;
      neighside = 0;
   }
   else 
   {
      mynode = 1;
      neighside = 1;
   }
   Node* NodeTemp2 = (Node*) NodeTable->lookup(nodes+mynode*KEYLENGTH);
   elev = 0.5*NodeTemp2->get_elevation();
   Element* EmTemp2 = (Element*) El_Table->lookup((EmTemp->get_neighbors()+KEYLENGTH*neighside));
   NodeTemp2 = (Node*) NodeTable->lookup(EmTemp2->getNode()+j*KEYLENGTH);
   elev += 0.5*NodeTemp2->get_elevation();
   return elev;
}
