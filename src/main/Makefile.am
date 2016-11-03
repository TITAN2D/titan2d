#
# automake template for src/main directory
#
# $Id: Makefile.am 134 2007-06-07 20:05:52Z dkumar $
#

bin_PROGRAMS = titan

titan_SOURCES = compare_key.C delete_tab.C hpfem.C datread.C hilbert.C init_piles.C update_topo.C restart.C

titan_CXXFLAGS = -I../gisapi/

titan_LDADD = ../gisapi/libgisapi.a ../adapt/libadapt.a @FLIBS@ \
  ../datstr/libdatstr.a ../geoflow/libgeoflow.a ../useful/libuseful.a \
  ../repartition/librepartition.a ../tecplot/libtecplot.a


