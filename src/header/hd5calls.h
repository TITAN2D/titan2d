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

const int XDMF_NEW = 0;
const int XDMF_OLD = 1;
const int XDMF_CLOSE = 2;

#include <hdf5.h>



#ifdef __cplusplus
#include <string>
extern "C"
{
#endif

/*! GH5_openfile opens hdf file for reading or writing.
 * mode can have "n" or "o" values
 * mode "n" opens a new file for writing data, if the file
 * of same name already exists, it fails and stops TITAN
 * mode "o" open existing file in read-only mode
 * If it fails to open the file, it stops TITAN
 */
hid_t GH5_fopen(const char *filename, char mode);

/*! Closed hdf file. Every file must be closed
 * after data writing is finished, or system will
 * run out of file descriptors
 */
inline void GH5_fclose(hid_t fp)
{
    herr_t status = H5Fclose(fp);
}

/*! Writes the mesh to hdf file. 
 * GH5_openfile must be called before this fucntion
 * fp is File id to which data is being written
 * conns specifies no. of connections, must be same as no of elements
 * points specifies no. of XYZ points, should be same as no. of nodes
 * conndata is 2D array of INTEGERS ("Num_Elems x 4") 
 * ptsdata is 2D array of doubles ("Num_Nodes x 3")
 */
void GH5_write_mesh_data(hid_t fp, int conns, int points, int *conndata, double *ptsdata);

/*! Writes state variables to hdf file.
 * GH5_openfile must be called before this fucntion
 * fp is File id to which data is being written
 * num_elms is total number if elements being written
 * state_vars is array of doubles of size (size=Num_Elems)
 * var_name is a strings, it contains names of the variables
 */
void GH5_write_state_vars(hid_t fp, int num_elms, double *state_var, const char *var_names);

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



#ifdef __cplusplus
#include <H5Cpp.h>

extern H5::EnumType datatypeElementType;

/**
 * initialize all hdf5 enum datatypes
 */
void init_TiH5();

#define TiH5_writeIntAttribute(group,value) TiH5_writeScalarAttribute(group, &value, #value,H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT)
#define TiH5_readIntAttribute(group,value) TiH5_readScalarAttribute(group, &value, #value,H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT)
#define TiH5_writeDoubleAttribute(group,value) TiH5_writeScalarAttribute(group, &value, #value,H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE)
#define TiH5_readDoubleAttribute(group,value) TiH5_readScalarAttribute(group, &value, #value,H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE)

#define TiH5_writeBoolAttribute(group,value) TiH5_writeBoolAttribute__(group, value, #value)
#define TiH5_readBoolAttribute(group,value) TiH5_readBoolAttribute__(group, value, #value)
#define TiH5_writeScalarDataTypeAttribute(group,value,type) TiH5_writeScalarDataTypeAttribute__(group, &value, #value,type)
#define TiH5_readScalarDataTypeAttribute(group,value,type) TiH5_readScalarDataTypeAttribute__(group, &value, #value,type)

#define TiH5_writeIntArrayAttribute(group,value,size) TiH5_writeArrayAttribute(group, size, value, #value, H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);
#define TiH5_readIntArrayAttribute(group,value,size) TiH5_readArrayAttribute(group, size, value, #value, H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);
#define TiH5_writeDoubleArrayAttribute(group,value,size) TiH5_writeArrayAttribute(group, size, value, #value, H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE);
#define TiH5_readDoubleArrayAttribute(group,value,size) TiH5_readArrayAttribute(group, size, value, #value, H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE);

#define TiH5_writeStringAttribute(group,value) TiH5_writeStringAttribute__(group, &value, #value)
#define TiH5_readStringAttribute(group,value) TiH5_readStringAttribute__(group, &value, #value)

#define TiH5_writeStringAttribute(group,value,length) TiH5_writeStringAttribute__(group, &value, #value,length)
#define TiH5_readStringAttribute(group,value,length) TiH5_readStringAttribute__(group, &value, #value,length)


inline void TiH5_writeScalarAttribute(H5::Group &group, const void *value, const char *name, const H5::DataType& typeRecord,const H5::DataType& typeNative)
{
    hsize_t dims = 1;
    // Create the data space for the attribute.
    H5::DataSpace attr_dataspace = H5::DataSpace (1, &dims );

    // Create a dataset attribute.
    H5::Attribute attribute = group.createAttribute(name, typeRecord, attr_dataspace);

    // Write the attribute data.
    attribute.write(typeNative, value);
}
inline void TiH5_readScalarAttribute(const H5::Group &group, void *value, const char *name, const H5::DataType& typeRecord,const H5::DataType& typeNative)
{
    // Create a dataset attribute.
    H5::Attribute attribute = group.openAttribute(name);

    // Read the attribute data.
    attribute.read(typeNative, value);
}
inline void TiH5_writeBoolAttribute__(H5::Group &group, const int value, const char *name)
{
    TiH5_writeScalarAttribute(group,&value,name,H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);
}
inline void TiH5_readBoolAttribute__(const H5::Group &group, bool &value, const char *name)
{
    int rvalue;
    TiH5_readScalarAttribute(group,&rvalue,name,H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);
    if(rvalue==0)value=false;
    else value=true;
}

inline void TiH5_writeScalarDataTypeAttribute__(H5::Group &group, const void *value, const char *name, const H5::DataType& type)
{
    hsize_t dims = 1;
    H5::DataSpace attr_dataspace = H5::DataSpace (1, &dims );

    // Create a dataset attribute.
    H5::Attribute attribute = group.createAttribute(name, type, attr_dataspace);

    // Write the attribute data.
    attribute.write(type, value);
}
inline void TiH5_readScalarDataTypeAttribute__(const H5::Group &group, void *value, const char *name, const H5::DataType& type)
{
    H5::Attribute attribute = group.openAttribute(name);
    attribute.read(type, value);
}

inline void TiH5_writeArrayAttribute(H5::Group &group, const hsize_t dims, const void *value, const char *name, const H5::DataType& typeRecord,const H5::DataType& typeNative)
{
    // Create the data space for the attribute.
    H5::DataSpace attr_dataspace = H5::DataSpace (1, &dims );

    // Create a dataset attribute.
    H5::Attribute attribute = group.createAttribute(name, typeRecord, attr_dataspace);

    // Write the attribute data.
    attribute.write(typeNative, value);
}
inline void TiH5_readArrayAttribute(const H5::Group &group, const hsize_t dims, void *value, const char *name, const H5::DataType& typeRecord,const H5::DataType& typeNative)
{
    H5::Attribute attribute = group.openAttribute(name);
    attribute.read(typeNative, value);
}
inline void TiH5_writeStringAttribute__(H5::Group &group, const std::string &value, const char *name, const int length=256)
{
    if(group.attrExists(name))
        group.removeAttr(name);
    H5::DataSpace attr_dataspace = H5::DataSpace (H5S_SCALAR);
    H5::StrType type(H5::PredType::C_S1, length);
    const H5std_string strwritebuf(value);
    H5::Attribute attribute = group.createAttribute(name, type, attr_dataspace);
    attribute.write(type, strwritebuf);
}
inline void TiH5_readStringAttribute__(const H5::Group &group, std::string &value, const char *name, const int length=256)
{
    H5::Attribute attribute = group.openAttribute(name);
    H5::StrType type(H5::PredType::C_S1, length);
    H5std_string strreadbuf ("");
    attribute.read(type, strreadbuf);
    value=strreadbuf;
}
#endif
#endif
