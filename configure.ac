# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(titan, 3.0.0)
AM_INIT_AUTOMAKE

## Host specific stuff
#===================================================================================
AC_CANONICAL_HOST
case $host in
  *-*-irix*)
    CFLAGS="$CFLAGS -woff all"
    CXXFLAGS="$CXXFLAGS -LANG:std -woff all"
    LDFLAGS="$LDFLAGS -LANG:std -w -lc -lC -lscs -lmpi++ "
    CXX="CC"
    F77="f77"
    ;;
  *-*-linux*)
    dnl better remove -w
    test_CFLAGS=${CFLAGS+set}
    if test "$test_CFLAGS" != set; then
    	CFLAGS="-g -w -O2"
    fi
    
    test_CXXFLAGS=${CXXFLAGS+set}
    if test "$test_CXXFLAGS" != set; then
    	CXXFLAGS="-g -w -O2 -std=c++11"
    else
    	CXXFLAGS="-std=c++11 $CXXFLAGS"
    fi
    
    CPPFLAGS="$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"
    
    test_LDFLAGS=${LDFLAGS+set}
    if test "$test_LDFLAGS" != set; then
    	LDFLAGS="-g -w -O2"
    fi
    ;;
esac

AC_CONFIG_SRCDIR([src/main/titan.C])
AC_CONFIG_HEADER([config.h])
AC_PREFIX_DEFAULT([${PWD}])
AC_LANG_CPLUSPLUS

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_F77_LIBRARY_LDFLAGS
AC_SYS_LARGEFILE

AM_PATH_PYTHON([2.7],, [:])
PYTHON_BIN_DIR=`AS_DIRNAME($PYTHON)`
PYTHON_ROOT=`AS_DIRNAME($PYTHON_BIN_DIR)`
AC_SUBST([PYTHON_ROOT],[$PYTHON_ROOT])




## checks for libraries
#======================================================================================
# checks for mpich installation
AC_CHECK_LIB([pthread],[pthread_create])

# in all likelyhood pkg-module is not going to find anything,
# this check has been included to build .rpm and .deb
#PKG_CHECK_MODULES([MPICH2], [mpich2-ch3])
CPPFLAGS="$CPPFLAGS $MPICH2_CFLAGS"
LDFLAGS="$LDFLAGS $MPICH2_LIBS"

# if MPI is installed in non-standard location
AC_ARG_WITH([mpi-includes],[AS_HELP_STRING([--with-mpi-includes],[path to mpi headers])])
if  test -n "$with_mpi_includes"; then
    CPPFLAGS="$CPPFLAGS -I$with_mpi_includes"
fi

AC_ARG_WITH([mpi-libs],[AS_HELP_STRING([--with-mpi-libs],[path to mpi libraries])])
if  test -n "$with_mpi_libs"; then
    LDFLAGS="$LDFLAGS -L$with_mpi_libs"
fi

#AC_CHECK_LIB([mpich],[MPI_Init],, \
#                   AC_CHECK_LIB([mpi],[MPI_Init],, \
#                                  AC_MSG_FAILURE([could not find MPI libraries.])))
## checks for HDF5 intalltion
#------------------------------------------------------------------------------------
AC_CHECK_LIB([z], [gzread])
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5],[HDF5 is needed for paraview support. \
                            Provide absolute HDF5 location --with-hdf5=/path/to/hdf5])])
if test "x$with_hdf5" != "xno"; then
    if test -n "$with_hdf5"; then
        CPPFLAGS="$CPPFLAGS  -I$with_hdf5/include -DH5_USE_16_API"
        LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
    fi
    AC_CHECK_LIB([hdf5],[H5open],,
                  AC_MSG_FAILURE([HDF5 is strongly recommended for paraview support.]
                  [Re-run configure --with-hdf5=/path/to/hdf5.]
                  [Re-run configure --without-hdf5 to get rid of this message]))
fi
AM_CONDITIONAL([BUILD_HDF5], test "x$with_hdf5" != "xno")

## checks for GRASS
#------------------------------------------------------------------------------------
AC_ARG_WITH([grass],
            [AS_HELP_STRING([--with-grass],[grass is needed for grass sites output.
                            Provide absolute grass location --with-grass=/path/to/grass])])

if test -n "$with_grass"; then
    CPPFLAGS="$CPPFLAGS  -I$with_grass/include"
    LDFLAGS="$LDFLAGS -L$with_grass/lib"
fi
AC_CHECK_LIB([gis],[G_gisinit])

#-------------------------------------------------------------------------------------
AC_ARG_WITH(gdal,
            AS_HELP_STRING([--with-gdal],[for non-grass DEM input formats (default=yes)]),,
            [with_gdal="yes"])

AC_ARG_WITH(gdal-config,
            AS_HELP_STRING([--with-gdal-config],
                           [--with-gdal-config=</full/path/to>/gdal-config if not in path]))

# check for gdal-config
AC_CHECK_PROG([GDAL_CONFIG],[gdal-config],[gdal-config],[no])

# user's choice should override default
if test  "x$with_gdal_config" != "x"; then
   GDAL_CONFIG=$with_gdal_config
fi
case "$with_gdal" in
    yes)
        # for building RPMS and DEBIAN packages
        GDAL_CFLAGS=`$GDAL_CONFIG --cflags`
        GDAL_LIBS=`$GDAL_CONFIG --libs`
        CPPFLAGS="$CPPFLAGS $GDAL_CFLAGS" 
        LDFLAGS="$LDFLAGS $GDAL_LIBS"
        AC_CHECK_LIB([gdal], [GDALOpen],[gdal_lib="yes"])
        AC_CHECK_HEADERS([gdal_priv.h, ogr_spatialref.h],[gdal_inc="yes"])
        ;;
    no)
        gdal_lib=
        gdal_inc=
        ;;
    *)
        CPPFLAGS="$CPPFLAGS -I$with_gdal/include"
        LDFLAGS="$LDFLAGS -L$with_gdal/lib"
        AC_CHECK_LIB([gdal], [GDALOpen],[gdal_lib="yes"],
                     [AC_MSG_FAILURE([couldn't find GDAL libraries])])
        AC_CHECK_HEADERS([gdal_priv.h, ogr_spatialref.h],[gdal_inc="yes"],
                     [AC_MSG_FAILURE([couldn't find GDAL headers])])
        ;;
esac
AC_MSG_CHECKING([for GDAL availability])
if test "x$gdal_inc" = "xyes" -a  "x$gdal_lib" = "xyes"; then
    AC_DEFINE([HAVE_GDAL], [1], [If Gdal is available])
    AC_MSG_RESULT([yes])
else
    AC_DEFINE([HAVE_GDAL], [0], [If Gdal is available])
    AC_MSG_RESULT([no])
fi
AM_CONDITIONAL([WITH_GDAL], test "x$with_gdal" != "xno")


#-------------------------------------------------------------------------------------
AC_ARG_VAR([SWIG],[path to swig, python api generator])
AC_PATH_PROG([SWIG],[swig],[])
if test "x$SWIG" = "x"; then
	AC_MSG_ERROR([swig is not found, set SWIG variable!])
fi
SWIG_VERSION=`$SWIG -version|grep ersion|cut -d ' ' -f 3`
SWIG_VERSION_MAJOR=`echo $SWIG_VERSION|cut -d . -f1`
if test "x$SWIG_VERSION_MAJOR" != "x3"; then
	AC_MSG_ERROR([swig is too old should be at least 3.0, set SWIG variable for custom location!])
fi
AC_MSG_NOTICE([swig version: $SWIG_VERSION])

#-------------------------------------------------------------------------------------
AC_CHECK_LIB(m,matherr)
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
                 doc/Makefile
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
                 src/vectordatapreproc/Makefile
                 src/main/titanvars.sh
                 tests/Makefile
                 tests/ticore/Makefile])


AC_OUTPUT
