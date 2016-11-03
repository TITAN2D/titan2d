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
 * $Id: hdftable.C 136 2007-06-07 20:18:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#if HAVE_LIBHDF5 == 1 && HAVE_LIBHDF5_HL == 1

/* Adding support for simple hdf table formatted output file with incremental data -- Abani*/
#include "../header/hpfem.h"
#include "H5TB.h"

void hdfviz_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int time_step, double time, int myid, int numprocs,
		MatProps* matprops)
{
  int i, k;
  int element_counter = 0;
  Element* EmTemp;
  HashEntry* entryp;
  char filename[17] = "viz_outputxxx.h5";
  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();
  double velocity_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the velocities
  double momentum_scale = sqrt(matprops->LENGTH_SCALE * (matprops->GRAVITY_SCALE)); // scaling factor for the momentums
  double time_scale = sqrt(matprops->GRAVITY_SCALE/matprops->LENGTH_SCALE );
  
#define NFIELDS 8
#define NRECORDS 1
  typedef struct gridpoint 
  {
    double xcoord;
    double ycoord;
    double elevation;
    double pile;
    int time_st;
    float tim;
    double xmomentum;
    double ymomentum;
  } gridpoint;

  gridpoint dst_buf;
  size_t dst_size =  sizeof( gridpoint );
  size_t dst_offset[NFIELDS] = { HOFFSET( gridpoint, xcoord ),
				 HOFFSET( gridpoint, ycoord ),
				 HOFFSET( gridpoint, elevation ),
				 HOFFSET( gridpoint, pile ),
				 HOFFSET( gridpoint, time_st ),
				 HOFFSET( gridpoint, tim ),
				 HOFFSET( gridpoint, xmomentum),
				 HOFFSET( gridpoint, ymomentum) };
  
  size_t dst_sizes[NFIELDS] = { sizeof( dst_buf.xcoord),
				sizeof( dst_buf.xcoord),
				sizeof( dst_buf.elevation),
				sizeof( dst_buf.pile),
				sizeof( dst_buf.time_st),
				sizeof( dst_buf.tim),
				sizeof( dst_buf.xmomentum),
				sizeof( dst_buf.ymomentum),
  };
  
  //gridpoint firstpoint = { 0., 0., 0., 0., 0,0.f,0.,0.};
  gridpoint * firstpoint = NULL;
  gridpoint nextpoint;
  
  /* Define field information */
  const char *field_names[NFIELDS]  = { "XCoord","Ycoord", "Elevation", "Pile", "Time step","Time", "XMomentum","Ymomentum"};
  hid_t      field_type[NFIELDS];
  hid_t      string_type;
  hid_t      file_id;
  hsize_t    chunk_size = 10;
  int        *fill_data = NULL;
  int        compress  = 1;
  herr_t     status; 
  /* Initialize the field field_type */
  
  field_type[0] =  H5T_NATIVE_DOUBLE;
  field_type[1] =  H5T_NATIVE_DOUBLE;
  field_type[2] =  H5T_NATIVE_DOUBLE;
  field_type[3] =  H5T_NATIVE_DOUBLE;
  field_type[4] =  H5T_NATIVE_INT;
  field_type[5] =  H5T_NATIVE_FLOAT;
  field_type[6] =  H5T_NATIVE_DOUBLE;
  field_type[7] =  H5T_NATIVE_DOUBLE;
  
  filename[10] = (myid % 1000)/100 + 48;
  filename[11] = (myid % 100)/10 + 48;
  filename[12] = (myid % 10) + 48;
  
  if (time_step == 0)
    {  /* Create a new file using default properties. */
      
      
      file_id = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
      
      status=H5TBmake_table( " TITAN OUT DATA", file_id, "Table1",(hsize_t) NFIELDS, (hsize_t)0, dst_size,   field_names, dst_offset, field_type, chunk_size, fill_data, compress, firstpoint  );
      
      /* Close the file. */
      H5Fclose( file_id );
    }
  
  
  /* open the file for incremental writes */
  
  file_id = H5Fopen(filename, H5F_ACC_RDWR,H5P_DEFAULT);
  
  // fill up data in on record and write 
  
  //output field data
  for(i=0; i<e_buckets; i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0)
	    {
	      double* state_vars = EmTemp->get_state_vars();
	      Node* NdTemp = (Node*) HT_Node_Ptr->lookup(EmTemp->pass_key());
	      if(state_vars[0] < GEOFLOW_TINY) {
		double zero = 0;
	
		nextpoint.xcoord=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
		nextpoint.ycoord=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE; 
		nextpoint.elevation=(NdTemp->get_elevation())*(matprops)->LENGTH_SCALE;
		nextpoint.pile=	state_vars[0]*(matprops)->HEIGHT_SCALE;
		nextpoint.time_st=time_step;
		nextpoint.tim = time*time_scale;
		nextpoint.xmomentum=	zero;
		nextpoint.ymomentum=	zero;

   /* append one record */
		if (time_step ==0)
		  status=H5TBappend_records( file_id, "Table1",(hsize_t)1, dst_size, dst_offset, &nextpoint );
		
	      }
	      else
		nextpoint.xcoord=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
	      nextpoint.ycoord=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
	        nextpoint.elevation=(NdTemp->get_elevation())*(matprops)->LENGTH_SCALE;
		nextpoint.pile=state_vars[0]*(matprops)->HEIGHT_SCALE, 
		nextpoint.time_st=time_step;
		nextpoint.tim = time*time_scale;
		nextpoint.xmomentum=momentum_scale*state_vars[1];
		nextpoint.ymomentum=momentum_scale*state_vars[2];

		/* append one record */
		status=H5TBappend_records( file_id, "Table1",(hsize_t)1, dst_size, dst_offset, &nextpoint );
	    }
	  
	  entryp = entryp->next; 
	  
	} 
    } 
  
  
  H5Fclose(file_id);
}
#endif
