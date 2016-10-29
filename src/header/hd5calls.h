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
const int XDMF_ONLYCLOSE = 3;

#define __STDC_FORMAT_MACROS
#include <hdf5.h>




#include <string>


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



#include <H5Cpp.h>

#include <iostream>
#include <string>
#include <vector>
#include <array>

#include <string.h>
#include <assert.h>
extern H5::EnumType datatypeElementType;
extern H5::EnumType datatypePileType;
extern H5::EnumType datatypeOutLineInitSize;

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

#define TiH5_writeDataTypeArrayAttribute(group,value,type,size) TiH5_writeArrayAttribute(group, size, value, #value, type, type);
#define TiH5_readDataTypeArrayAttribute(group,value,type,size) TiH5_readArrayAttribute(group, size, value, #value, type, type);

#define TiH5_writeIntVectorAttribute(group,value,size) TiH5_writeArrayAttribute(group, size, &(value[0]), #value, H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);
#define TiH5_readIntVectorAttribute(group,value,size) {value.resize(size);TiH5_readArrayAttribute(group, size, &(value[0]), #value, H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT);}

#define TiH5_writeDoubleVectorAttribute(group,value,size) TiH5_writeArrayAttribute(group, size, &(value[0]), #value, H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE);
#define TiH5_readDoubleVectorAttribute(group,value,size) {value.resize(size);TiH5_readArrayAttribute(group, size, &(value[0]), #value, H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE);}

#define TiH5_writeDataTypeVectorAttribute(group,value,type,size) TiH5_writeArrayAttribute(group, size, &(value[0]), #value, type, type);
#define TiH5_readDataTypeVectorAttribute(group,value,type,size) {value.resize(size);TiH5_readArrayAttribute(group, size, &(value[0]), #value, type, type);}


#define TiH5_writeStringAttribute(group,value) TiH5_writeStringAttribute__(group, value, #value)
#define TiH5_readStringAttribute(group,value) TiH5_readStringAttribute__(group, value, #value)


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
    if(dims>0)
    {
        // Create the data space for the attribute.
        H5::DataSpace attr_dataspace = H5::DataSpace (1, &dims );

        // Create a dataset attribute.
        H5::Attribute attribute = group.createAttribute(name, typeRecord, attr_dataspace);

        // Write the attribute data.
        attribute.write(typeNative, value);
    }
}
inline void TiH5_readArrayAttribute(const H5::Group &group, const hsize_t dims, void *value, const char *name, const H5::DataType& typeRecord,const H5::DataType& typeNative)
{
    if(dims>0)
    {
        H5::Attribute attribute = group.openAttribute(name);
        attribute.read(typeNative, value);
    }
}
#define TiH5_writeVectorArrayAttribute(group,value) TiH5_writeVectorArrayAttribute__(group, value, #value);
#define TiH5_readVectorArrayAttribute(group,value)  TiH5_readVectorArrayAttribute__(group, value, #value);

template<std::size_t N>
inline void TiH5_writeVectorArrayAttribute__(H5::Group &group, const std::vector<std::array<double,N> > &value, const char *name)
{
    const hsize_t dims[2]={value.size(),N};
    double *v=nullptr;
    if(dims[0]>0)
    {
        if(dims[0]>1)
        {
            v=new double[dims[0]*dims[1]];
            for(std::size_t i=0;i<dims[0];++i)
                for(std::size_t j=0;j<dims[1];++j)
                    v[i*dims[1]+j]=value[i][j];
        }
        // Create the data space for the attribute.
        H5::DataSpace attr_dataspace = H5::DataSpace (2, dims );

        // Create a dataset attribute.
        H5::Attribute attribute = group.createAttribute(name, H5::PredType::IEEE_F64LE, attr_dataspace);

        // Write the attribute data.
        if(dims[0]>1)
            attribute.write(H5::PredType::NATIVE_DOUBLE, v);
        else
            attribute.write(H5::PredType::NATIVE_DOUBLE, &(value[0][0]));

        if(dims[0]>1)delete [] v;
    }
}
template<std::size_t N>
inline void TiH5_readVectorArrayAttribute__(const H5::Group &group, std::vector<std::array<double,N> > &value, const char *name)
{
    if(!group.attrExists(name))
    {
        value.resize(0);
        return;
    }

    H5::Attribute attribute = group.openAttribute(name);
    H5::DataSpace attr_dataspace=attribute.getSpace();


    hsize_t dims[2];
    attr_dataspace.getSimpleExtentDims(dims);
    if(dims[1]!=N)
    {
        std::cout<<"ERROR:dimentions in hdf do not match the program\n";
        assert(0);
    }

    value.resize(dims[0]);
    if(dims[0]>0)
    {
        if(dims[0]>1)
        {
            double *v=new double[dims[0]*dims[1]];
            attribute.read(H5::PredType::NATIVE_DOUBLE, v);
            for(std::size_t i=0;i<dims[0];++i)
                for(std::size_t j=0;j<dims[1];++j)
                    value[i][j]=v[i*dims[1]+j];
            delete [] v;
        }
        else
        {
            attribute.read(H5::PredType::NATIVE_DOUBLE, &(value[0][0]));
        }
    }
}
#define TiH5_writeVectorStringAttribute(group,value) TiH5_writeVectorStringAttribute__(group, value, #value)
inline void TiH5_writeVectorStringAttribute__(H5::Group &group, const std::vector<std::string> &value, const char *name)
{
    if(value.size()==0)
    {
        return;
    }

    std::size_t max_string_size=0;
    for(std::size_t i=0;i<value.size();++i)
    {
        if(value[i].size()>max_string_size)max_string_size=value[i].size();
    }
    ++max_string_size;

    char *v=new char[value.size()*max_string_size];

    for(std::size_t i=0;i<value.size();++i)
    {
        memcpy(v+i*max_string_size,value[i].c_str(),value[i].size());
        for(std::size_t j=value[i].size();j<max_string_size;++j)
            v[i*max_string_size+j]='\0';
    }
    H5::StrType strtype=H5::StrType(H5::PredType::C_S1,max_string_size);
    H5::CompType datatype=H5::CompType(max_string_size);
    datatype.insertMember("string",0,strtype);
    hsize_t dims[1]={value.size()};
    H5::DataSpace dataspace = H5::DataSpace (1, dims );
    H5::Attribute attribute = group.createAttribute(name, datatype, dataspace);
    attribute.write(datatype, v);

    delete [] v;
}
#define TiH5_readVectorStringAttribute(group,value) TiH5_readVectorStringAttribute__(&group, value, #value)
inline void TiH5_readVectorStringAttribute__(const H5::Group *group, std::vector<std::string> &value, const char *name)
{
    if(!group->attrExists(name))
    {
        value.resize(0);
        return;
    }

    H5::Attribute attribute = group->openAttribute(name);
    H5::DataSpace dataspace=attribute.getSpace();
    H5::DataType datatype=attribute.getDataType();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);

    value.resize(dims[0]);
    std::size_t max_string_size=datatype.getSize();

    char *v=new char[value.size()*max_string_size];
    attribute.read(datatype, v);

    for(std::size_t i=0;i<value.size();++i)
    {
        value[i]=v+i*max_string_size;
    }

    delete [] v;
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

#define TiH5_writeArrayDataSet(group,value,size) TiH5_writeArray2DDataSet__(group, value, #value,size)
inline void TiH5_writeArrayDataSet__(H5::Group &group, const double *value, const char *name, const hsize_t dims)
{
    hsize_t chunk_dims=10240;
    if(chunk_dims>dims)chunk_dims=dims;

    // Create the data space for the dataset
    H5::DataSpace dataspace(1, &dims);

    // Modify dataset creation property to enable chunking
    H5::DSetCreatPropList plist;
    plist.setChunk(1, &chunk_dims);
    // Set ZLIB (DEFLATE) Compression using level 6.
    plist.setDeflate(6);

    // Create the dataset.
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::IEEE_F64LE, dataspace,plist);
    // Write the attribute data.
    dataset.write(value, H5::PredType::NATIVE_DOUBLE);
}
#define TiH5_writeArray2DDataSet(group,value,size1,size2) TiH5_writeArray2DDataSet__(group, value, #value,size1,size2)
inline void TiH5_writeArray2DDataSet__(H5::Group &group, const double *value, const char *name, const hsize_t dims1, const hsize_t dims2)
{
    hsize_t chunk_dims[2]={64,64};
    const hsize_t dims[2]={dims1,dims2};

    if(chunk_dims[0]>dims[0])chunk_dims[0]=dims[0];
    if(chunk_dims[1]>dims[1])chunk_dims[1]=dims[1];

    // Create the data space for the dataset
    H5::DataSpace dataspace(2, dims);

    // Modify dataset creation property to enable chunking
    H5::DSetCreatPropList plist;
    plist.setChunk(2, chunk_dims);
    // Set ZLIB (DEFLATE) Compression using level 6.
    plist.setDeflate(6);

    // Create the dataset.
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::IEEE_F64LE, dataspace, plist);
    // Write the attribute data.
    dataset.write(value, H5::PredType::NATIVE_DOUBLE);
}
#define TiH5_readArray2DDataSet(group,value,size1,size2) TiH5_readArray2DDataSet__(group, value, #value,size1,size2)
inline void TiH5_readArray2DDataSet__(const H5::Group &group, double *value, const char *name, const hsize_t dims1, const hsize_t dims2)
{
    if(dims1*dims2>0)
    {
        hsize_t dims[2]={dims1,dims2};
        H5::DataSet dataset = group.openDataSet(name);
        H5::DataSpace dataspace=dataset.getSpace();

        dataspace.getSimpleExtentDims(dims);

        if(dims[0]!=dims1||dims[1]!=dims2)
        {
            std::cout<<"Error: dimension(s) in hdf file is different when in program\n";
            assert(0);
        }
        dataset.read(value, H5::PredType::NATIVE_DOUBLE);
    }
}






#endif
