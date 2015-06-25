#
# automake template for src/main directory
#
# $Id: Makefile.am 134 2007-06-07 20:05:52Z dkumar $
#

bin_PROGRAMS = titan
bin_SCRIPTS = titanvars.sh

titan_SOURCES = titan.C cxxtitan_wrap.cxx titan_simulation.C compare_key.C delete_tab.C hpfem.C datread.C hilbert.C init_piles.C update_topo.C restart.C

titan_CXXFLAGS = -fpic -I$(srcdir)/../gisapi -I@PYTHON_ROOT@/include/python@PYTHON_VERSION@

titan_LDADD = ../preproc/libtitan_preprocess.a ../gisapi/libgisapi.a ../adapt/libadapt.a @FLIBS@ \
  ../datstr/libdatstr.a ../geoflow/libgeoflow.a ../useful/libuseful.a \
  ../repartition/librepartition.a ../tecplot/libtecplot.a \
  ../vectordatapreproc/libvectordatpreproc.a \
  -L@PYTHON_ROOT@/lib -lpython@PYTHON_VERSION@ -lutil

#titan_simulation_CXXFLAGS = -fpic

cxxtitan_wrap.cxx:cxxtitan.i ../header/titan_simulation.h ../preproc/preproc.h
	@SWIG@ -c++ -v -python -o $(srcdir)/cxxtitan_wrap.cxx -oh $(srcdir)/cxxtitan.py $(srcdir)/cxxtitan.i



pyfilesdir=$(pythondir)/titan
pyfiles_PYTHON=cxxtitan.py titan.py

