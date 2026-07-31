# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(titan2d, 4.2.0)
AM_INIT_AUTOMAKE

## Check is portable package is made
#===================================================================================
AC_ARG_ENABLE(portable,
	[AS_HELP_STRING([--enable-portable],[Make universal binaries (default: disabled)])],
	[],[enable_portable=no]
)
AS_IF([test "x$enable_portable" != xno],
	[
	    TITAN2D_PORTABLE=yes
		AC_SUBST(MAKE_PORTABLE)
		LDFLAGS="$LDFLAGS -Wl,-rpath,'\$\$ORIGIN/../lib:\$\$ORIGIN/../lib/titan2d_dep/lib'"
	],
	[TITAN2D_PORTABLE=no]
)
AC_SUBST([TITAN2D_PORTABLE],[$TITAN2D_PORTABLE])

## Check if titan dependencies are provided
#===================================================================================
AC_ARG_WITH([titan2d-dep],
            [AS_HELP_STRING([--with-titan2d-dep],[Specify path to Titan2d dependencies package. \
                            Provide absolute HDF5 location --with-hdf5=/path/to/titan2d_ep])],
            [],[with_titan2d_dep=no]
)
AS_IF([test "x$with_titan2d_dep" != xno],
	[
		case $with_titan2d_dep in
		/*)
		    ;;
		*)
			AC_MSG_ERROR([with_titan2d_dep should be full path!])
			;;
		esac
		AC_MSG_RESULT([Titan2D dependencies package is at $with_titan2d_dep])
		AS_IF([test "x$PYTHON" == x],
		[
			PYTHON=$with_titan2d_dep/bin/python
		],[])
		TITAN2D_DEP=$with_titan2d_dep
		AS_IF([test "x$enable_portable" == xno],
			[
				LDFLAGS="$LDFLAGS -Wl,-rpath,'\$\$ORIGIN/../lib:\$(TITAN2D_DEP)/lib'"
			]
		)
		AC_CHECK_FILE([$with_titan2d_dep/bin/python],[],
		[AC_MSG_ERROR([Titan2D dependency package ($with_titan2d_dep) do not contain python. Check path.])]
		)
	],
	[TITAN2D_DEP=no]
) 
AC_SUBST([TITAN2D_DEP],[$TITAN2D_DEP])
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
      if test "x$enable_portable" != xno; then
        CXXFLAGS="-g -w -O3 -mtune=generic -std=c++11"
      else
        CXXFLAGS="-g -w -O2 -std=c++11"
      fi
    else
      CXXFLAGS="-std=c++11 $CXXFLAGS"
    fi

    test_LDFLAGS=${LDFLAGS+set}
    if test "$test_LDFLAGS" != set; then
      LDFLAGS="-g -w -O2"
    fi
    ;;
  *-*-darwin*)
    test_CXXFLAGS=${CXXFLAGS+set}
    if test "$test_CXXFLAGS" != set; then
      CXXFLAGS="-g -O2 -std=c++11"
    else
      CXXFLAGS="-std=c++11 $CXXFLAGS"
    fi
    ;;
esac

AS_IF([test "x$enable_portable" != xno],
	[CPPFLAGS=-DMAKE_PORTABLE $CPPFLAGS]
	[]
)  

AC_CONFIG_SRCDIR([src/main/titan.C])
AC_CONFIG_HEADER([src/header/titan_config.h])
AC_PREFIX_DEFAULT([${PWD}])
AC_LANG_CPLUSPLUS

# Checks for programs.
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

AM_PATH_PYTHON([3.8],, [:])
PYTHON_BIN_DIR=`AS_DIRNAME($PYTHON)`
PYTHON_ROOT=`AS_DIRNAME($PYTHON_BIN_DIR)`
AC_PATH_PROGS([PYTHON_CONFIG],
	[python$PYTHON_VERSION-config python3-config python-config],
	[no],
	[$PYTHON_BIN_DIR:$PATH]
)
AS_IF([test "x$PYTHON_CONFIG" == xno],
	[AC_MSG_FAILURE([Could not find a matching python-config for $PYTHON])]
)

PYTHON_CPPFLAGS=`$PYTHON_CONFIG --includes`
PYTHON_EMBED_LDFLAGS=`$PYTHON_CONFIG --ldflags --embed 2>/dev/null`
AS_IF([test "x$PYTHON_EMBED_LDFLAGS" == x],
	[PYTHON_EMBED_LDFLAGS=`$PYTHON_CONFIG --ldflags`]
)

AC_SUBST([PYTHON_ROOT],[$PYTHON_ROOT])
AC_SUBST([PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
AC_SUBST([PYTHON_EMBED_LDFLAGS],[$PYTHON_EMBED_LDFLAGS])

AH_TEMPLATE([TITAN_PYTHONHOME],[eventually PYTHONHOME])
AC_DEFINE_UNQUOTED([TITAN_PYTHONHOME],["$PYTHON_ROOT"])
## debug option
#======================================================================================

## if titan dependencies are provided set common flags
#===================================================================================
AS_IF([test "x$with_titan2d_dep" != xno],
	[
		CPPFLAGS="$CPPFLAGS  -I$with_titan2d_dep/include"
	    LDFLAGS="$LDFLAGS -L$with_titan2d_dep/lib"
	],
)

## checks for libraries
#======================================================================================

# checks for OpenMP
#------------------------------------------------------------------------------------
AC_ARG_ENABLE(openmp,
	[AS_HELP_STRING([--enable-openmp],[Compile with OpenMP support (default: disabled)])],
	[],[enable_openmp=no]
)
AS_IF([test "x$enable_openmp" != xno],
	[
		AC_OPENMP
		AS_IF([test "x$ac_cv_prog_cxx_openmp" = xunsupported],
			[
				for titan_libomp_prefix in /opt/homebrew/opt/libomp /usr/local/opt/libomp; do
					AS_IF([test -f "$titan_libomp_prefix/include/omp.h" && test -f "$titan_libomp_prefix/lib/libomp.dylib"],
						[
							AC_MSG_NOTICE([using Homebrew libomp at $titan_libomp_prefix])
							OPENMP_CXXFLAGS="-Xpreprocessor -fopenmp"
							CPPFLAGS="$CPPFLAGS -I$titan_libomp_prefix/include"
							LDFLAGS="$LDFLAGS -L$titan_libomp_prefix/lib -Wl,-rpath,$titan_libomp_prefix/lib"
							LIBS="$LIBS -lomp"
							ac_cv_prog_cxx_openmp="$OPENMP_CXXFLAGS"
							break
						]
					)
				done
			]
		)
		AS_IF([test "x$ac_cv_prog_cxx_openmp" = xunsupported],
			[AC_MSG_FAILURE([Can not find OpenMP support for $CXX. On macOS install Homebrew libomp or use a compiler/runtime pair that provides OpenMP.])]
		)
		CPPFLAGS="$CPPFLAGS -DTITAN2D_OPENMP_RUNTIME=1"
		titan_save_CXXFLAGS="$CXXFLAGS"
		CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
		AC_CHECK_FUNCS([omp_get_num_threads],,AC_MSG_FAILURE([Can not compile simple OpenMP program]))
		CXXFLAGS="$titan_save_CXXFLAGS"
	]
)

# checks for MPI installation
#------------------------------------------------------------------------------------
AC_ARG_ENABLE(mpi,
	[AS_HELP_STRING([--enable-mpi],[Compile with MPI support, CXX variable should point to MPI enabled compiler (default: disabled)])],
	[],[enable_mpi=no]
)
AS_IF([test "x$enable_mpi" != xno],
	[
		CPPFLAGS="$CPPFLAGS -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK"
		AC_CHECK_FUNCS([MPI_Init],,AC_MSG_FAILURE([Can not compile simple MPI program]))
	]
)  

# checks for libz
#------------------------------------------------------------------------------------
AC_CHECK_LIB([z], [gzread],,AC_MSG_FAILURE([Can not find libz]),)

# checks for libm
#-------------------------------------------------------------------------------------
AC_CHECK_LIB(m,matherr,)

# checks for HDF5 installation
#------------------------------------------------------------------------------------
AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5],[Specify path to hdf5 library. \
                    Provide absolute HDF5 location --with-hdf5=/path/to/hdf5])],
    [],[with_hdf5=no]
)

AS_IF([test "x$with_hdf5" != xno],
	[
	    CPPFLAGS="$CPPFLAGS  -I$with_hdf5/include"
	    LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
	]
)
AC_CHECK_LIB([hdf5],[H5open],,
              AC_MSG_FAILURE([can not compile against HDF5 library. If HDF5 not in standard paths re-run configure --with-hdf5=/path/to/hdf5.])
)
AC_CHECK_LIB([hdf5_cpp],[H5open],,
              AC_MSG_FAILURE([HDF5 should be compiled with c++ support.Re-run configure --with-hdf5=/path/to/hdf5.])
)
AC_MSG_CHECKING([whether HDF5 C++ API links with $CXX])
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM(
		[[#include <H5Cpp.h>
#include <string>]],
		[[std::string name("Titan2DConfigureCheck");
H5::Exception::dontPrint();
H5::H5File file(name, H5F_ACC_TRUNC);
file.createGroup(std::string("group"));]]
	)],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_FAILURE([HDF5 C++ API does not link with $CXX. Use an HDF5 C++ library built with the same C++ ABI/compiler, or use Apple clang with Homebrew HDF5 on macOS.])]
)
AC_SUBST([HDF5_PATH],[$with_hdf5])

## checks for GDAL
#-------------------------------------------------------------------------------------
AC_ARG_WITH(gdal,
    [AS_HELP_STRING([--with-gdal],[GDAL is used for non-grass DEM input formats. \
        Specify path to GDAL library. \
        Provide absolute GDAL location --with-gdal=/path/to/gdal (default: try to find in standard paths)])],
    [],
    [with_gdal="try_to_find"]
)

AS_IF([test "x$with_gdal" != xno],
	[titan_gdal_requested=yes],
	[titan_gdal_requested=no]
)

GDAL_CONFIG=no
GDAL_CFLAGS=
GDAL_LIBS=
GDAL_PATH=
titan_gdal_found=no
titan_gdal_required=no

AS_IF([test "x$titan_gdal_requested" = xyes],
	[
		titan_gdal_search_path="$PATH"
		AS_IF([test "x$with_gdal" != xtry_to_find],
			[
				titan_gdal_required=yes
				AS_IF([test ! -d "$with_gdal"],
					[AC_MSG_FAILURE([Requested GDAL prefix does not exist: $with_gdal])]
				)
				titan_gdal_search_path="$with_gdal/bin:$titan_gdal_search_path"
			]
		)
		AS_IF([test "x$with_titan2d_dep" != xno],
			[titan_gdal_search_path="$with_titan2d_dep/bin:$titan_gdal_search_path"]
		)

		# Prefer gdal-config when it is available.
		AC_PATH_PROG([GDAL_CONFIG], [gdal-config], [no], [$titan_gdal_search_path])
		AS_IF([test "x$GDAL_CONFIG" != xno],
			[
				GDAL_CFLAGS=`$GDAL_CONFIG --cflags`
				GDAL_LIBS=`$GDAL_CONFIG --libs`
				titan_gdal_bindir=`dirname "$GDAL_CONFIG"`
				GDAL_PATH=`cd "$titan_gdal_bindir/.." && pwd -P`
				titan_gdal_found=yes
			],
			[
				# Some modern packages only publish pkg-config metadata.
				AC_PATH_PROG([PKG_CONFIG], [pkg-config], [no])
				AS_IF([test "x$PKG_CONFIG" != xno],
					[
						titan_gdal_pkg_path="$PKG_CONFIG_PATH"
						AS_IF([test "x$with_gdal" != xtry_to_find],
							[titan_gdal_pkg_path="$with_gdal/lib/pkgconfig:$with_gdal/share/pkgconfig:$titan_gdal_pkg_path"]
						)
						AS_IF([PKG_CONFIG_PATH="$titan_gdal_pkg_path" $PKG_CONFIG --exists gdal],
							[
								GDAL_CFLAGS=`PKG_CONFIG_PATH="$titan_gdal_pkg_path" $PKG_CONFIG --cflags gdal`
								GDAL_LIBS=`PKG_CONFIG_PATH="$titan_gdal_pkg_path" $PKG_CONFIG --libs gdal`
								GDAL_PATH=`PKG_CONFIG_PATH="$titan_gdal_pkg_path" $PKG_CONFIG --variable=prefix gdal`
								titan_gdal_found=yes
							]
						)
					]
				)
			]
		)
	]
)

AC_MSG_CHECKING([for GDAL availability])
AS_IF([test "x$titan_gdal_found" = xyes],
	[
		with_gdal="$GDAL_PATH"
		CPPFLAGS="$CPPFLAGS $GDAL_CFLAGS"
		LIBS="$GDAL_LIBS $LIBS"
		AC_DEFINE([HAVE_GDAL], [1], [If Gdal is available])
		AC_MSG_RESULT([yes])
	],
	[
		GDAL_CFLAGS=
		GDAL_LIBS=
		GDAL_PATH=
		with_gdal=no
		AC_DEFINE([HAVE_GDAL], [0], [If Gdal is available])
		AC_MSG_RESULT([no])
		AS_IF([test "x$titan_gdal_required" = xyes],
			[AC_MSG_FAILURE([GDAL was requested but neither gdal-config nor gdal.pc was found])]
		)
	]
)
AC_SUBST([GDAL_PATH],[$GDAL_PATH])
AC_SUBST([GDAL_CFLAGS],[$GDAL_CFLAGS])
AC_SUBST([GDAL_LIBS],[$GDAL_LIBS])
AM_CONDITIONAL([WITH_GDAL], test "x$titan_gdal_found" = "xyes")


#SWIG
#-------------------------------------------------------------------------------------
AC_ARG_WITH(swig,
            AS_HELP_STRING([--with-swig=<path to swig home directory>],
            [Specify path to SWIG package. ]
            [If set will (re)generate python API, otherwise will try to use previously generated sources. ]
            [Note that github repository do not track them.]
            [Provide absolute path to swig top directory, --with-swig=/path/to/swig/home (default: without)]),
            [],[with_swig=$with_titan2d_dep])
#only developers/advanced user should call swig
#the public source code releases should come with swig generated files
AS_IF([test "x$with_swig" != xno],
[
    AC_PATH_PROGS([SWIG],[swig swig3 swig3.0],[no],[$with_swig/bin:$PATH])
	if test "x$SWIG" = "xno"; then
		AC_MSG_ERROR([swig is not found, set --with-swig variable!])
	fi
	SWIG_VERSION=`$SWIG -version|grep ersion|cut -d ' ' -f 3`
	SWIG_VERSION_MAJOR=`echo $SWIG_VERSION|cut -d . -f1`
	if test "$SWIG_VERSION_MAJOR" -lt 3; then
		AC_MSG_ERROR([swig is too old should be at least 3.0, set SWIG variable for custom location!])
	fi
	
	AS_IF([test "x$with_titan2d_dep" != xno],
	[
		SWIG="SWIG_LIB=$with_titan2d_dep/share/swig/3.0.5 $SWIG"
	]
	) 
	
	AC_MSG_NOTICE([swig: $SWIG])
	AC_MSG_NOTICE([swig version: $SWIG_VERSION])
],
[
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
])          
AM_CONDITIONAL([WITH_SWIG], test "x$with_swig" != "xno")

## set rpath if without_titan2d_dep
#------------------------------------------------------------------------------------
AS_IF([test "x$with_titan2d_dep" == xno],
	[
		PYTHON_LIBDIR=`$PYTHON -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR') or '')"`
		RPATH=$PYTHON_LIBDIR
		AS_IF([test "x$with_hdf5" != xno],
			[RPATH=$with_hdf5/lib:$RPATH]
		)
		AS_IF([test "x$with_gdal" != xno],
			[RPATH=$RPATH:$with_gdal/lib]
		)
		LDFLAGS="$LDFLAGS -Wl,-rpath,'\$\$ORIGIN/../lib:$RPATH'"
	]
	
)

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
		AS_IF([test "x$JAVA" == x], [AC_PATH_PROG([JAVA],[java],[no],[$with_titan2d_dep/lib/java-1.7/bin:$PATH])])
		AS_IF([test "x$JAVAC" == x], [AC_PATH_PROG([JAVAC],[javac],[no],[$with_titan2d_dep/lib/java-1.7/bin:$PATH])])
		AS_IF([test "x$JAR" == x], [AC_PATH_PROG([JAR],[jar],[no],[$with_titan2d_dep/lib/java-1.7/bin:$PATH])])
		
		AS_IF([test "x$JAVA" == xno], [AC_MSG_FAILURE([Can not find java, set JAVA variable or --disable-java-gui])])
		AS_IF([test "x$JAVAC" == xno], [AC_MSG_FAILURE([Can not find javac, set JAVAC variable or --disable-java-gui])])
		AS_IF([test "x$JAR" == xno], [AC_MSG_FAILURE([Can not find jar, set JAR variable or --disable-java-gui])])
		
		JAVA_VERSION1=`$JAVA -version 2>&1 |grep 'openjdk\|java version'|cut -d '"' -f 2|cut -d '.' -f 1`
		JAVA_VERSION2=`$JAVA -version 2>&1 |grep 'openjdk\|java version'|cut -d '"' -f 2|cut -d '.' -f 2`
		
		AS_IF([test "x$JAVA_VERSION1" == x1], [], [AC_MSG_ERROR([Java should be 1.7 or 1.8!])])
		AS_IF([test "x$JAVA_VERSION2" == x7 || test "x$JAVA_VERSION2" == x8], [], [AC_MSG_ERROR([Java should be 1.7 or 1.8!])])
		
		AC_MSG_NOTICE([java $JAVA])
		AC_MSG_NOTICE([javac $JAVAC])
		AC_MSG_NOTICE([jar $JAR])
		
	]
)
AM_CONDITIONAL([ENABLE_JAVA_GUI], [test "$enable_java_gui" = "yes"])

# compile for vhub
#------------------------------------------------------------------------------------
AC_ARG_ENABLE(vhub,
	[AS_HELP_STRING([--enable-vhub],[Compile for VHub (default: disabled)])],
	[],[enable_vhub=no]
)
AS_IF([test "x$enable_vhub" != xno],
	[TITAN2D_ON_VHUB=yes],
	[TITAN2D_ON_VHUB=no]
)
AC_SUBST([TITAN2D_ON_VHUB],[$TITAN2D_ON_VHUB])

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

## Clean up some stuff
#======================================================================================
LDFLAGS=`echo $LDFLAGS |sed 's|-L/usr/lib64 ||g;s|-L/usr/lib ||g'`

#======================================================================================
AC_CONFIG_FILES([Makefile
                 gui/Makefile
                 gui/titan_gui.sh
                 doc/Makefile
                 doc/help/Makefile
                 doc/help/html/Makefile
                 scripts/Makefile
                 scripts/driverKML.sh
                 scripts/driverMultiProcessor.sh
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
AC_MSG_RESULT([AS_HELP_STRING([CPPFLAGS],[$CPPFLAGS])])
AC_MSG_RESULT([AS_HELP_STRING([LDFLAGS],[$LDFLAGS])])
AC_MSG_RESULT([AS_HELP_STRING([LIBS],[$LIBS])])
AS_IF([test "x$enable_portable" != xno],
	[
	AC_MSG_RESULT([AS_HELP_STRING([portable],[yes])])
	]
)
AS_IF([test "x$enable_openmp" != xno],
	[AC_MSG_RESULT([AS_HELP_STRING([OpenMP],[enabled])])],
	[AC_MSG_RESULT([AS_HELP_STRING([OpenMP],[disabled])])]
)
AS_IF([test "x$enable_mpi" != xno],
	[AC_MSG_RESULT([AS_HELP_STRING([MPI],[enabled])])],
	[AC_MSG_RESULT([AS_HELP_STRING([MPI],[disabled])])]
)
AS_IF([test "x$with_gdal" != xno],
	[
		AC_MSG_RESULT([AS_HELP_STRING([GDAL],[enabled])])
		AC_MSG_RESULT([AS_HELP_STRING([GDAL_CONFIG],[$GDAL_CONFIG])])
		AC_MSG_RESULT([AS_HELP_STRING([GDAL_LIBS],[$GDAL_LIBS])])
	],
	[AC_MSG_RESULT([AS_HELP_STRING([GDAL],[disabled])])]
)
AC_MSG_RESULT([AS_HELP_STRING([SWIG],[$SWIG])])

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
AC_MSG_RESULT([AS_HELP_STRING([TITAN2D_DEP],[$TITAN2D_DEP])])
AC_MSG_RESULT([AS_HELP_STRING([TITAN2D_PORTABLE],[$TITAN2D_PORTABLE])])
AC_MSG_RESULT([AS_HELP_STRING([TITAN2D_ON_VHUB],[$TITAN2D_ON_VHUB])])
