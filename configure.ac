# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(titan2d, 3.1.0)
AM_INIT_AUTOMAKE

## Host specific stuff
#===================================================================================
AC_CANONICAL_HOST
case $host in
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
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_SYS_LARGEFILE

AM_PATH_PYTHON([2.7],, [:])
PYTHON_BIN_DIR=`AS_DIRNAME($PYTHON)`
PYTHON_ROOT=`AS_DIRNAME($PYTHON_BIN_DIR)`
AC_SUBST([PYTHON_ROOT],[$PYTHON_ROOT])

## debug option
#======================================================================================

## checks for libraries
#======================================================================================

# checks for omp
#------------------------------------------------------------------------------------
AC_ARG_ENABLE(openmp,
	[AS_HELP_STRING([--enable-openmp],[Compile with OpenMP support (default:disabled)])],
	[],[enable_openmp=no]
)
AS_IF([test "x$enable_openmp" != xno],
	[
		AC_OPENMP
		CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
		AC_CHECK_FUNCS([omp_get_num_threads],,AC_MSG_FAILURE([Can not compile simple OpenMP program]))
	]
)
#AC_CHECK_LIB([pthread],[pthread_create])
# checks for mpich installation
#------------------------------------------------------------------------------------
AC_ARG_ENABLE(mpi,
	[AS_HELP_STRING([--enable-mpi],[Compile with MPI support, CXX variable should point to MPI enabled compiler (default:disabled)])],
	[],[enable_mpi=no]
)
AS_IF([test "x$enable_mpi" != xno],
	[
		CPPFLAGS="$CPPFLAGS -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK"
		AC_CHECK_FUNCS([MPI_Init],,AC_MSG_FAILURE([Can not compile simple MPI program]))
	]
)  

# checks for HDF5 installation
#------------------------------------------------------------------------------------
AC_CHECK_LIB([z], [gzread],,AC_MSG_FAILURE([Can not find libz]),)
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5],[HDF5 is needed for paraview support. \
                            Provide absolute HDF5 location --with-hdf5=/path/to/hdf5])],)

if test -n "$with_hdf5"; then
    CPPFLAGS="$CPPFLAGS  -I$with_hdf5/include"
    LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
fi

AC_CHECK_LIB([hdf5],[H5open],,
              AC_MSG_FAILURE([HDF5 is strongly recommended for paraview support.]
              [Re-run configure --with-hdf5=/path/to/hdf5.]
              [Re-run configure --without-hdf5 to get rid of this message]))
LIBS="-lhdf5_cpp $LIBS"


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
        GDAL_PATH=`$GDAL_CONFIG --prefix`
        CPPFLAGS="$CPPFLAGS $GDAL_CFLAGS" 
        LIBS="$GDAL_LIBS $LIBS"
        AC_CHECK_LIB([gdal], [GDALOpen],[gdal_lib="yes"])
        AC_CHECK_HEADERS([gdal_priv.h, ogr_spatialref.h],[gdal_inc="yes"])
        ;;
    no)
        gdal_lib=
        gdal_inc=
        GDAL_PATH=
        ;;
    *)
        CPPFLAGS="$CPPFLAGS -I$with_gdal/include"
        LIBS="$GDAL_LIBS $LIBS"
        GDAL_PATH="$with_gdal"
        AC_CHECK_LIB([gdal], [GDALOpen],[gdal_lib="yes"],
                     [AC_MSG_FAILURE([couldn't find GDAL libraries])])
        AC_CHECK_HEADERS([gdal_priv.h, ogr_spatialref.h],[gdal_inc="yes"],
                     [AC_MSG_FAILURE([couldn't find GDAL headers])])
        ;;
esac
AC_SUBST([GDAL_PATH],[$GDAL_PATH])
AC_MSG_CHECKING([for GDAL availability])
if test "x$gdal_inc" = "xyes" -a  "x$gdal_lib" = "xyes"; then
    AC_DEFINE([HAVE_GDAL], [1], [If Gdal is available])
    AC_MSG_RESULT([yes])
else
    AC_DEFINE([HAVE_GDAL], [0], [If Gdal is available])
    AC_MSG_RESULT([no])
fi
AM_CONDITIONAL([WITH_GDAL], test "x$with_gdal" != "xno")

#SWIG
#-------------------------------------------------------------------------------------
AC_ARG_WITH(swig,
            AS_HELP_STRING([--with-swig=<path to swig executable>],
            [Specify path to SWIG executable. ]
            [If set will (re)generate python API, otherwise will try to use previously generated sources. ]
            [Note that github repository do not track them.]),)
#only developers/advanced user should call swig
#the public source code releases should come with swig generated files
if test -n "$with_swig"; then
    SWIG=$with_swig
    AC_PATH_PROG([SWIG],[swig],[])
	if test "x$SWIG" = "x"; then
		AC_MSG_ERROR([swig is not found, set SWIG variable!])
	fi
	SWIG_VERSION=`$SWIG -version|grep ersion|cut -d ' ' -f 3`
	SWIG_VERSION_MAJOR=`echo $SWIG_VERSION|cut -d . -f1`
	if test "x$SWIG_VERSION_MAJOR" != "x3"; then
		AC_MSG_ERROR([swig is too old should be at least 3.0, set SWIG variable for custom location!])
	fi
	AC_MSG_NOTICE([swig: $SWIG])
	AC_MSG_NOTICE([swig version: $SWIG_VERSION])
else
    #check for swig generated files presence
    AC_CHECK_FILE($srcdir/src/main/cxxtitan_wrap.cxx,
    	AC_MSG_NOTICE([found cxxtitan_wrap.cxx]),
    	AC_MSG_FAILURE([Couldn't find SWIG generated python API sources.]
    	[Set --with-swig to SWIG to regenerate interface.])
    )
    AC_CHECK_FILE($srcdir/src/main/cxxtitan.py,
    	AC_MSG_NOTICE([found cxxtitan_wrap.cxx]),
    	AC_MSG_FAILURE([Couldn't find SWIG generated python API sources.]
    	[Set --with-swig to SWIG to regenerate interface.])
    )
fi           
AM_CONDITIONAL([WITH_SWIG], test "x$with_swig" != "xno")

## checks for java
#------------------------------------------------------------------------------------
AC_ARG_VAR(JAVA,[java interpretator (needed for java GUI)])
AC_ARG_VAR(JAVAC,[java compiler (needed for java GUI)])
AC_ARG_VAR(JAR,[java packager (needed for java GUI)])
AC_ARG_ENABLE(java-gui,
	[AS_HELP_STRING([--disable-java-gui],[Do not compile java GUI (default compile)])],
	[],
	[enable_java_gui=yes]
)
AS_IF([test "x$enable_java_gui" != xno],
	[
		AS_IF([test "x$JAVA" == x], [AC_PATH_PROG([JAVA],[java])])
		AS_IF([test "x$JAVAC" == x], [AC_PATH_PROG([JAVAC],[javac])])
		AS_IF([test "x$JAR" == x], [AC_PATH_PROG([JAR],[jar])])
		
		JAVA_VERSION1=`$JAVA -version 2>&1 |grep 'java version'|cut -d '"' -f 2|cut -d '.' -f 1`
		JAVA_VERSION2=`$JAVA -version 2>&1 |grep 'java version'|cut -d '"' -f 2|cut -d '.' -f 2`
		
		AS_IF([test "x$JAVA_VERSION1" == x1], [], [AC_MSG_ERROR([Java should be 1.7 or 1.8!])])
		AS_IF([test "x$JAVA_VERSION2" == x7 || test "x$JAVA_VERSION2" == x8], [], [AC_MSG_ERROR([Java should be 1.7 or 1.8!])])
		
		AC_MSG_RESULT("HELLO JAVA")
		
		AC_MSG_NOTICE([java $JAVA])
		AC_MSG_NOTICE([javac $JAVAC])
		AC_MSG_NOTICE([jar $JAR])
		
	]
)

#-------------------------------------------------------------------------------------
AC_CHECK_LIB(m,matherr)

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
                 gui/Makefile
                 doc/Makefile
                 doc/help/Makefile
                 doc/help/html/Makefile
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
                 tests/ticore/Makefile
                 share/Makefile])


AC_OUTPUT

## Configuration Summary
#======================================================================================
AC_MSG_RESULT()
AC_MSG_RESULT(Configuration Summary:)
AC_MSG_RESULT([AS_HELP_STRING([CXX],[$CXX])])
AC_MSG_RESULT([AS_HELP_STRING([CXXFLAGS],[$CXXFLAGS])])
AC_MSG_RESULT([AS_HELP_STRING([LDFLAGS],[$LDFLAGS])])
AS_IF([test "x$enable_openmp" != xno],
	[AC_MSG_RESULT([AS_HELP_STRING([OpenMP],[enabled])])],
	[AC_MSG_RESULT([AS_HELP_STRING([OpenMP],[disabled])])]
)
AS_IF([test "x$enable_mpi" != xno],
	[AC_MSG_RESULT([AS_HELP_STRING([MPI],[enabled])])],
	[AC_MSG_RESULT([AS_HELP_STRING([MPI],[disabled])])]
)
AS_IF([test "x$enable_java_gui" != xno],
	[
		AC_MSG_RESULT([AS_HELP_STRING([java GUI],[enabled])])
		AC_MSG_RESULT([AS_HELP_STRING([JAVA],[$JAVA])])
		AC_MSG_RESULT([AS_HELP_STRING([JAVAC],[$JAVAC])])
		AC_MSG_RESULT([AS_HELP_STRING([JAR],[$JAR])])
	],
	[
		AC_MSG_RESULT([AS_HELP_STRING([java GUI],[disabled])])
	]
)



