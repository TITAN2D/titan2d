# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(titan, 2.0.2)
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_SRCDIR([src/main/hpfem.C])
AC_CONFIG_HEADER([config.h])
AC_PREFIX_DEFAULT([${PWD}])
AC_LANG_CPLUSPLUS

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_F77_LIBRARY_LDFLAGS
AC_SYS_LARGEFILE

## checks for libraries
#======================================================================================
# checks for mpich installation
AC_CHECK_LIB([pthread],[pthread_create])
AC_ARG_WITH([mpi],[AS_HELP_STRING([--with-mpi],[MPI is required to build titan])])

if test -n "$with_mpi"; then
    CFLAGS="$CFLAGS -I$with_mpi/include"
    CXXFLAGS="$CXXFLAGS -I$with_mpi/include"
    CPPFLAGS="$CPPFLAGS -I$with_mpi/include"
    LDFLAGS="$LDFLAGS -L$with_mpi/lib"
    AC_CHECK_LIB([mpich],[MPI_Init],, \
                   AC_CHECK_LIB([mpi],[MPI_Init],, \
                                  AC_MSG_FAILURE([could not find MPI libraries \
                                                 in specified location])))
else
    AC_CHECK_LIB([mpich],[MPI_Init],, \
                   AC_CHECK_LIB([mpi],[MPI_Init],, \
                                  AC_MSG_FAILURE([could not find MPI libraries \
                                                 run configure --with-mpi=/path/to/mpi dir])))
fi

## checks for HDF5 intalltion
#------------------------------------------------------------------------------------
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5],[HDF5 is needed for paraview support. \
                            Provide absolute HDF5 location --with-hdf5=/path/to/hdf5])])
if test "x$with_hdf5" != xno; then
    if test -n "x$with_hdf5"; then
        CFLAGS="$CFLAGS -I$with_hdf5/include"
        CXXFLAGS="$CXXFLAGS  -I$with_hdf5/include"
        CPPFLAGS="$CPPFLAGS  -I$with_hdf5/include -DH5_USE_16_API"
        LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
        AC_CHECK_LIB([hdf5],[H5open])
    fi
    AC_CHECK_LIB([hdf5],[H5open],AC_DEFINE([HAVE_HDF5],[1],[Build XDMF output]),
                  AC_MSG_FAILURE([HDF5 is strongly recommended for paraview support.]
                  [Re-run configure --with-hdf5=/path/to/hdf5.]
                  [Re-run configure --without-hdf5 to get rid of this message]))
fi

## checks for GRASS
#------------------------------------------------------------------------------------
AC_ARG_WITH([grass],
            [AS_HELP_STRING([--with-grass],[grass is needed for grass sites output.
                            Provide absolute grass location --with-grass=/path/to/grass])])

if test -n "$with_grass"; then
    CFLAGS="$CFLAGS -I$with_grass/include"
    CXXFLAGS="$CXXFLAGS  -I$with_grass/include"
    CPPFLAGS="$CPPFLAGS  -I$with_grass/include"
    LDFLAGS="$LDFLAGS -L$with_grass/lib"
fi
AC_CHECK_LIB([gis],[G_gisinit])
              
#=============================================================
# enable debug
AC_ARG_ENABLE(debug,
              AS_HELP_STRING([--enable-debug],
                             [enable debugging (disabled by default)]),
              [HAVE_DEBUG=$enableval], [HAVE_DEBUG=yes])

if test "x$enable_debug" = "xyes" ; then

# Try to strip out optimization.
  CXXFLAGS=`echo "$CXXFLAGS " | sed "s/-O2 //g"`
  CFLAGS=`echo "$CFLAGS " | sed "s/-O2 //g"`
  CXXFLAGS=`echo "$CXXFLAGS" | sed "s/-O //g"`
  CFLAGS=`echo "$CFLAGS" | sed "s/-O //g"`

# remove existing -g to avoid duplicates.
  CXXFLAGS=`echo "$CXXFLAGS" | sed "s/-g //"`
  CFLAGS=`echo "$CFLAGS" | sed "s/-g //"`

#  -g and -DDEBUG
  CXXFLAGS="$CXXFLAGS -O0 -g"
  CFLAGS="$CFLAGS -O0 -g"

  AC_MSG_RESULT([yes, modified CFLAGS="$CFLAGS"])
else
  AC_MSG_RESULT([no, CFLAGS="$CFLAGS"])
fi

#-------------------------------------------------------------------------------------
AC_CHECK_LIB(m,matherr)
AC_CHECK_LIB(z,gzdopen)
AC_CHECK_LIB(nsl, inet_ntoa)

## Checks for header files
#=====================================================================================
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h mpi.h hdf5.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

## Checks for library functions.
#======================================================================================
AC_HEADER_MAJOR
AC_CHECK_FUNCS([floor pow sqrt strchr strdup])
#======================================================================================
AC_CONFIG_FILES([Makefile
                 scripts/Makefile
                 src/Makefile
                 src/adapt/Makefile
                 src/datstr/Makefile
                 src/geoflow/Makefile
                 src/gisapi/Makefile
                 src/header/Makefile
                 src/main/Makefile
                 src/preproc/Makefile
                 src/repartition/Makefile
                 src/stochastic/Makefile
                 src/tecplot/Makefile
                 src/useful/Makefile
                 src/vectordatapreproc/Makefile])
AC_OUTPUT
