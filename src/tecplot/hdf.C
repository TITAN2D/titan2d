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
 * $Id: hdf.C 136 2007-06-07 20:18:23Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#if HAVE_LIBHDF5 == 1 && HAVE_LIBHDF5_HL == 1
#include "../header/hpfem.h"

// hdf library  header file
#include "hdf5.h"
#include "hdfdefs.h"

/*********************************************************************************/
//  HDF OUTPUT  
/********************************************************************************/

void hdf_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int time_step, double time, int myid, int numprocs,MatProps* matprops)
{
  int i, k;
  int element_counter = 0;
  Element* EmTemp;
  HashEntry* entryp;
  char datavalue[] = "1-X CORDINATE,2-Y CORDINATE,3-Z CORDINATE, 4- PILE HEIGHT, 5-X VELOCITY,6- Y VELOCITY,7-SLOPE(x dir),8-SLOPE(y-dir),9-Key unsigned1,10-Key unsigned2";
  char datavalue1[] = "1-Time Step, 2-Time ,3-No of points, 4-No. of rows in dataset";
  char filename1[17] = "hdf_outputxxx.h5";
  char afilename[18] = "hdf_aoutputxxx.h5";

  int e_buckets=HT_Elem_Ptr->get_no_of_buckets();
  hid_t      file,aid3,aid4,atype;
  hid_t       afile;

  hid_t      dataset, memspace, space, filespace,xfer,attribute,attribute1;
  hid_t       adataset, amemspace, aspace, afilespace;
  herr_t     status;
  hsize_t    chkdim[3]={1,HEIGHT, WIDTH };
  hsize_t    achkdim[3]={1,aHEIGHT, aWIDTH };
  hssize_t offset[3]={0, 0, 0};
  hssize_t aoffset[3]={0, 0, 0};
  hsize_t count1[3]={LENGTH, HEIGHT1, WIDTH};
  hsize_t acount1[3]={LENGTH, aHEIGHT1, aWIDTH};

  
  // global descriptor -- list of file names

  filename1[10] = (myid % 1000)/100 + 48;
  filename1[11] = (myid % 100)/10 + 48;
  filename1[12] = (myid % 10) + 48;

  afilename[11] = (myid % 1000)/100 + 48;
  afilename[12] = (myid % 100)/10 + 48;
  afilename[13] = (myid % 10) + 48;

  // actual data output
  for(i=0; i<e_buckets; i++)
    {
      entryp = *(HT_Elem_Ptr->getbucketptr() + i);
      while(entryp)
	{	
	  EmTemp = (Element*)entryp->value;
	  assert(EmTemp);
	  if(EmTemp->get_adapted_flag()>0) 
	    element_counter++;
	  
	  entryp = entryp->next;
	}
    }
  //information below will probably be changed later

  int points = element_counter;
  int variables = 7;
  int inputs = 12;
  int parameters = 2;
 
  if(time_step == 0){



    w[myid]=0;
    aw[myid]=0;
    r[myid]=0;
    ar[myid]=0;
    l[myid]=0;
    al[myid]=0;
    l1[myid]=0;
    al1[myid]=0;


  }
 
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
             

	      if(state_vars[0] < GEOFLOW_TINY) 
		{
		  double zero = 0;
		  
		  
		  if(l[myid]<HEIGHT)
		    {
			    
		      //data[0][l[myid]][0]=EmTemp->coord[0];
		      data[0][l[myid]][0]=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
		      //data[0][l[myid]][1]= EmTemp->coord[1];
		      data[0][l[myid]][1]=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
		      data[0][l[myid]][2]=EmTemp->get_elevation();
			
		      data[0][l[myid]][3]=	zero;
		      data[0][l[myid]][4]=  zero;
		      data[0][l[myid]][5]=	zero;
		      data[0][l[myid]][6]=*(EmTemp->get_zeta());
		      data[0][l[myid]][7]= *(EmTemp->get_zeta()+1);
		      data[0][l[myid]][8]=*(EmTemp->pass_key());
		      data[0][l[myid]][9]=  *(EmTemp->pass_key()+1);
		    }
		  else
		    {
			 
		      l1[myid] = l[myid]-HEIGHT1*w[myid]-HEIGHT*r[myid];

		      //data1[0][l1[myid]][0]=EmTemp->coord[0];
		      //data1[0][l1[myid]][1]= EmTemp->coord[1];

		      data1[0][l1[myid]][0]=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
			data1[0][l1[myid]][1]=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
		      data1[0][l1[myid]][2]=EmTemp->get_elevation();
		      data1[0][l1[myid]][3]= zero;
		      data1[0][l1[myid]][4]= zero;
		      data1[0][l1[myid]][5]= zero ;
		      data1[0][l1[myid]][6]=*(EmTemp->get_zeta());
		      data1[0][l1[myid]][7]= *(EmTemp->get_zeta()+1);
		      data1[0][l1[myid]][8]=*(EmTemp->pass_key());
		      data1[0][l1[myid]][9]=  *(EmTemp->pass_key()+1);

		    }
		   
			  
		}
		
	      else
		{
		  
	
		  if(l[myid]<HEIGHT)
		    {
			
		      //data[0][l[myid]][0]=EmTemp->coord[0];
		      //data[0][l[myid]][1]= EmTemp->coord[1];
		      data1[0][l1[myid]][0]=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
			data1[0][l1[myid]][1]=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
		      data[0][l[myid]][2]=EmTemp->get_elevation();
		      data[0][l[myid]][3]=	state_vars[0];
		      data[0][l[myid]][4]=  (state_vars[1]/state_vars[0]);
		      data[0][l[myid]][5]=	(state_vars[2]/state_vars[0]) ;
		      data[0][l[myid]][6]=*(EmTemp->get_zeta());
		      data[0][l[myid]][7]= *(EmTemp->get_zeta()+1);
		      data[0][l[myid]][8]=*(EmTemp->pass_key());
		      data[0][l[myid]][9]=  *(EmTemp->pass_key()+1);

		    }   
		  else
		    {
		      l1[myid] = l[myid]-HEIGHT1*w[myid]-HEIGHT*r[myid];

			
		      //data1[0][l1[myid]][0]=EmTemp->coord[0];
		      //data1[0][l1[myid]][1]= EmTemp->coord[1];
		      data1[0][l1[myid]][0]=(*(EmTemp->get_coord()))*(matprops)->LENGTH_SCALE;
			data1[0][l1[myid]][1]=(*(EmTemp->get_coord()+1))*(matprops)->LENGTH_SCALE;
		      data1[0][l1[myid]][2]=EmTemp->get_elevation();
		      data1[0][l1[myid]][3]=	state_vars[0];
		      data1[0][l1[myid]][4]= (state_vars[1]/state_vars[0]);
		      data1[0][l1[myid]][5]=	(state_vars[2]/state_vars[0]) ;
		      data1[0][l1[myid]][6]=*(EmTemp->get_zeta());
		      data1[0][l1[myid]][7]= *(EmTemp->get_zeta()+1);
		      data1[0][l1[myid]][8]=*(EmTemp->pass_key());
		      data1[0][l1[myid]][9]=  *(EmTemp->pass_key()+1);
		    }  
		    
               
		}
	      
	      
	      l[myid]++;
	      // writing of datasets 

	      if(l[myid]==HEIGHT)
		{
		  space = H5Screate_simple (RANK, dim, maxdim);
                  file = H5Fcreate(filename1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		  printf ("H5Fcreate: %i\n", file);
		  cparms = H5Pcreate (H5P_DATASET_CREATE);
		  printf ("H5Pcreate: %i\n", cparms);
		  status = H5Pset_chunk ( cparms, RANK, chkdim);
		  printf ("H5Pset_chunk: %i\n", status);
		  tbuf = (float*) malloc (500000);
		  xfer = H5Pcreate (H5P_DATASET_XFER);
		  status= H5Pset_buffer (xfer, (hsize_t)500000, tbuf, NULL); 
		  dataset = H5Dcreate(file, DATASETNAME, H5T_NATIVE_FLOAT, space, cparms);
		  printf ("H5Dcreate: %i\n", dataset);
		  status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, xfer, data);
		  printf ("H5Dwrite: %i\n", status);

		  aid4 = H5Screate(H5S_SCALAR);
		  atype = H5Tcopy(H5T_C_S1);
		  H5Tset_size(atype, 150);
		  attribute1 = H5Acreate(dataset, "Units",atype , aid4, H5P_DEFAULT);

		  status = H5Awrite(attribute1,atype, datavalue);
                  status=H5Sclose(aid4);
		  status=H5Aclose(attribute1);
		  status= H5Pclose (xfer);
		  printf ("H5Pclose: %i\n", status);
		  status= H5Pclose (cparms);
		  printf ("H5Pclose: %i\n", status);
		  status=H5Sclose(space);
		  printf ("H5Sclose: %i\n", status);
 
		  status=H5Dclose(dataset);
		  printf ("H5Dclose: %i\n", status);
		  status=H5Fclose(file);
		  printf ("H5Fclose: %i\n", status);

		  free (tbuf);  
		  r[myid]++;
		}
		
	   
	      if((l[myid]>HEIGHT)&&(l[myid]%HEIGHT1)==0)
		{
		  file = H5Fopen (filename1, H5F_ACC_RDWR, H5P_DEFAULT);
		  dataset = H5Dopen (file, DATASETNAME);
		  offset[1] = newsize[1]; 
		  newsize[1]=newsize[1]+ HEIGHT1;
		  status = H5Dextend (dataset, newsize);
		  filespace = H5Dget_space (dataset);
		  status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
						count1, NULL);  
		  memspace= H5Screate_simple (RANK, count1, NULL); 
		  status = H5Dwrite (dataset, H5T_NATIVE_FLOAT, memspace, filespace,
				     H5P_DEFAULT, data1);  
		  status = H5Sclose (memspace);
		  status = H5Dclose (dataset);
		  status = H5Sclose (filespace);
		  status = H5Fclose (file);
		  printf ("H5Fclose: %i\n", status);
		  w[myid]++;
		}  
		
	    }     	
	  entryp = entryp->next; 
	} 
    
    }
  
  // writing attributes of the datasets 
           
  if(al[myid]<aHEIGHT)
    {
      adata[0][al[myid]][0]=time_step;
      adata[0][al[myid]][1]=  time;
      adata[0][al[myid]][2]=points;
      adata[0][al[myid]][3]=HEIGHT*r[myid]+HEIGHT1*w[myid];
    }  
  else
    {
      al1[myid]=al[myid]-aHEIGHT1*aw[myid]-ar[myid]*aHEIGHT;

      adata1[0][al1[myid]][0]=time_step;
      adata1[0][al1[myid]][1]=  time;
      adata1[0][al1[myid]][2]=points;
      adata1[0][al1[myid]][3]= HEIGHT*r[myid]+HEIGHT1*w[myid];
    } 
   
         
 
  al[myid]++;

  if(al[myid]==aHEIGHT)
    {
      aspace = H5Screate_simple (RANK, adim, amaxdim);
      /* Create the file  */
      afile = H5Fcreate(afilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      printf ("H5Fcreate: %i\n", afile);
      acparms = H5Pcreate (H5P_DATASET_CREATE);
      printf ("H5Pcreate: %i\n", acparms);
      status = H5Pset_chunk ( acparms, RANK, achkdim);
      printf ("H5Pset_chunk: %i\n", status);
      adataset = H5Dcreate(afile, DATASETNAME, H5T_NATIVE_FLOAT, aspace,acparms );
      printf ("H5Dcreate: %i\n", adataset);
      status = H5Dwrite(adataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT , adata);
      aid3 = H5Screate(H5S_SCALAR);
      atype = H5Tcopy(H5T_C_S1);
      H5Tset_size(atype, 150);
      attribute = H5Acreate(adataset, "Units",atype , aid3, H5P_DEFAULT);

      status = H5Awrite(attribute,atype, datavalue1);
      status=H5Sclose(aid3);
      status=H5Aclose(attribute);
      status=H5Sclose(aspace);
      printf ("H5Sclose: %i\n", status);
      status=H5Dclose(adataset);
      printf ("H5Dclose: %i\n", status);
      status=H5Fclose(afile);
      printf ("H5Fclose: %i\n", status);
      ar[myid]++;
    }
    

 

  if((al[myid]>aHEIGHT)&&(al[myid]%aHEIGHT1)==0)
    {
      afile = H5Fopen (afilename, H5F_ACC_RDWR, H5P_DEFAULT);
      adataset = H5Dopen (afile, DATASETNAME);
      aoffset[1] = anewsize[1]; 
      anewsize[1]=anewsize[1]+ aHEIGHT1;
      status = H5Dextend (adataset, anewsize);
      afilespace = H5Dget_space (adataset);
      status = H5Sselect_hyperslab (afilespace, H5S_SELECT_SET, aoffset, NULL,
				    acount1, NULL); 
      amemspace= H5Screate_simple (RANK, acount1, NULL); 
      status = H5Dwrite (adataset, H5T_NATIVE_FLOAT, amemspace, afilespace,
			 H5P_DEFAULT, adata1);     
      status = H5Sclose (amemspace);
      status = H5Dclose (adataset);
      status = H5Sclose (afilespace);
      status = H5Fclose (afile);
      printf ("H5Fclose: %i\n", status);
      aw[myid]++;
    
    }  
    

 
  return;
    
} 
#endif












  
  
  
  
      
