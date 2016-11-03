#! /bin/bash

## check for MPI installation
MPCC=
MPCC=$(which mpicc 2> /dev/null)  
MPCPP=
MPCPP=`which mpiCC 2> /dev/null` || MPCPP=`which mpicxx 2> /dev/null` 
mpi_home=
if [ "x$MPCPP" == x ]; then
    echo '-----------------------------------------------------------'
    echo ' MPI installtion not found. MPI is required to build titan'
    echo '-----------------------------------------------------------'
    echo 'Enter absolute path to MPI installtion directory if installed in non-standard place.'
    echo '[ctrl-c to quit]'
    read mpi_home
fi

## check for HDF5
hdf5_dir=
hcc=`which h5cc`
tempdir=`$hcc -showconfig | fgrep 'Installation point' | awk '{print $3}'`
echo $tempdir
header='include/hdf5.h'
if [ -f $tempdir/$header ]; then
    hdf5_dir=$tempdir
#look for hdf5 in some usual places
elif [ -f /usr/$header ]; then

    hdf5_dir=/usr
elif [ -f /usr/local/$header ]; then
    hdf5_dir=/usr/local
# ${HOME}
elif [ -f ${HOME}/$header ]; then
    hdf5_dir=${HOME}
# /opt
elif [ -f /opt/$header ]; then
    hdf5_dir=/opt
elif

if [[ "x$hdf5_dir" == x ]]; then
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+  Unable to locate HDF5 installtion. It might not be installed  +"
    echo "+  or installed in unusual path. PARAVIEW support needs hdf5     +"
    echo "+  Please provide absolute path to hdf5 installtion.             +"
    echo "+  HDF5 is freely downloadle from http://hdf.ncsa.uiuc.edu/HDF5/ +"
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "Continue without HDF5? (yes/no) [no]"
    read response
    if [[ "x$response" != xyes ]]; then
        echo "If HDF5 is installed in non-standard path, "
        echo "  provide /full/path/to/hdf5 ...(ctrl-c to QUIT)"
        read hdf5_dir
        if [ ! -f $hdf5_dir/$header ]; then
            echo "ERROR: Could not find hdf5.h in $hdf5_dir"
            exit 1;
        fi
    hdf5_flags='-D H5_USE_16_API'
    else
        hdf5_dir=no
    fi
fi

# run configure
MPI=
if test -n "$MPCPP"; then
   MPI="CC=$MPCC CXX=$MPCPP"
else
   MPI="--with-mpi=$mpi_home"
fi
   
./configure $MPI --with-hdf5=$hdf5_dir CPPFLAGS=$hdf5_flags
make && make install
