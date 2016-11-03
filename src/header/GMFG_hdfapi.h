/*
 * =====================================================================================
 *
 *       Filename:  GMFG_hdf5api.h
 *
 *    Description:  TITAN specific API Calls to hdf5
 *
 *        Created:  05/17/2007 10:26:15 AM EDT
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

#ifndef __GMFG_HDF5API__
#define __GMFG_HDF5API__

const int XDMF_NEW=0;
const int XDMF_OLD=1;
const int XDMF_CLOSE=2;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_HDF5_H
#include <hdf5.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! GH5_openfile opens hdf file for reading or writing.
 * mode can have "n" or "o" values
 * mode "n" opens a new file for writing data, if the file
 * of same name already exists, it fails and stops TITAN
 * mode "o" open existing file in read-only mode
 * If it fails to open the file, it stops TITAN
 */
hid_t GH5_openfile(const char *filename, char mode);

/*! Closed hdf file. Every file must be closed
 * after data writing is finished, or system will
 * run out of file descriptors
 */
inline void GH5_closefile(hid_t fp)
{
   herr_t status=H5Fclose(fp);
}

/*! Writes the mesh to hdf file. 
 * GH5_openfile must be called before this fucntion
 * fp is File id to which data is being written
 * conns specifies no. of connections, must be same as no of elements
 * points specifies no. of XYZ points, should be same as no. of nodes
 * conndata is 2D array of INTEGERS ("Num_Elems x 4") 
 * ptsdata is 2D array of doubles ("Num_Nodes x 3")
 */
void GH5_write_mesh_data(hid_t fp,int conns,int points,int *conndata,double *ptsdata);

/*! Writes state variables to hdf file.
 * GH5_openfile must be called before this fucntion
 * fp is File id to which data is being written
 * num_elms is total number if elements being written
 * state_vars is array of doubles of size (size=Num_Elems)
 * var_name is a strings, it contains names of the variables
 */
void GH5_write_state_vars(hid_t fp,int num_elms,double *state_var,const char *var_names);


#ifdef __cplusplus
}
#endif
/*************************************************************
 *
 *      FOLLOWING CALLS ARE NOT EXPECTED TO BE MADE
 *      FROM TITAN. THEY ARE TO BE USED IN API ITSELF
 *
 * ***********************************************************/

/*! Creates a group within file. If first argument is
 * a file id, name should contain absolute path.
 * If first argument is a group id, relative path
 * can be used
 */
hid_t GH5_open_group(hid_t fp, const char *name);


/*! Creates a dataset to which actual data is to be written 
 *  gid is intended to be group id, to which data is to be written.
 *  however it can also be a file id, in that case dsetname must be absolute path.
 *  types are defined in GMFG_hdfapi.h, they represents type of the data in dataset
 *  rank is order of of the array, i.e. 1D,2D etc
 *  dims contains dimensions, i.e. rows, cols etc
 */
hid_t GH5_createdataset(hid_t gid, hid_t spcid, const char *dsetname, unsigned type);

#endif
