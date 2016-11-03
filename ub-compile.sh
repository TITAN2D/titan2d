#!/bin/bash
#
# Automatically sets compilation parameters for computers ta UB
# by A. Sorokine <sorokine@buffalo.edu>, based on Makefiles
#
# $Id: ub-compile.sh 138 2007-06-07 20:54:54Z dkumar $
#

host=`hostname`
make_jobs=''

case "$host" in 
#    alpha*)
#        mpi_home=${HOME}/mpich-1.2.6
#    ;;
        
    ganesha*)
	mpi_home=/eng/mae/ganesha4/mpich-1.2.1
    ;;

    alpha*|rainier*)
	mpi_home=${HOME}/mpich2-1.0.4
#        mpi_home=/usr/local/mpich-1.2.6
    ;;

    popo*|colima*|elchichon*|sthelens*)
	mpi_home=/usr/local/mpich-1.2.4/ch_p4/
	hdf_home=/volcano/elchichon1/HDF5/gnu/
	hdfhl_home=/volcano/elchichon1/HDF5/hl-gnu/
    ;;

    casita*)
	 mpi_home=${HOME}/mpich-1.2.6
    ;;
#       mpi_home=/usr/local/mpich-1.2.5/
#	hdf_home=/volcano/elchichon1/HDF5/gnu/
#	hdfhl_home=/volcano/elchichon1/HDF5/hl-gnu/

    nash*)
	mpi_home=/util/mpich-gm/gnu/ch_gm
	hdf_home=/util/hdf/
	hdfhl_home=/util/hdf_hl/
    ;;

    joplin*)
	mpi_home=/util/mpich-gm/gnu/ch_gm
	hdf_home=/util/hdf/
	hdfhl_home=/util/hdf_hl/
    ;;

    crosby*)
	hdf_home=/util/HDF/HDF5/n32pp
	hdfhl_home=/util/HDF/HDF5_hl/n32
    ;;

    lennon.*)
	export LDFLAGS=-lmpi
	make_jobs=-j5
    ;;

    bono*|u2*|pe2950-dev*)
#	mpi_home=/util/mpich/1.2.6/gcc-3.4.3/ch_p4
#	mpi_home=/util/mpich/1.2.7p1-fixed/gcc-3.4.6/ch_gm
#	mpi_home=/util/openmpi/1.2/gnu-3.4.6
	export CXXFLAGS="-openmp -openmp_report2"
#	mpi_home=/util/openmpi/1.2/intel-9
	mpi_home=/util/mpich/1.2.7p1-fixed/intel-9.1/ch_gm
    ;;

    localhost*)
	mpi_home=/usr/local/mpich-1.2.5
	hdf_home=/usr/local/5-1.4.5-linux/
	hdfhl_home=/usr/local/hdf5_hl-ia32/
    ;;

    *)
	echo No host-specific settings for $host, using defaults
    ;;

esac

./configure --prefix=`pwd` \
    `test -n "$mpi_home" && echo "--with-mpi=$mpi_home"` \
    `test -z "$mpi_home" && echo "--disable-mpi-compilers"` \
    --without-grass --without-hdf5 && make $make_jobs && make install
#    `test -n "$hdf_home" && echo "--with-hdf5=$hdf_home"` \
#    `test -n "$hdfhl_home" && echo "--with-hdf5_hl=$hdfhl_home"` \
