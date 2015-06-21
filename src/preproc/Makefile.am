#
# automake template for preproc/ directory
#
# $Id: Makefile.am 12 2003-11-07 17:58:49Z kdalbey $
#

bin_PROGRAMS = titan_materialnames
#titan_preprocess

titan_materialnames_SOURCES = materialnames.C

titan_materialnames_CXXFLAGS = -I$(srcdir)/../gisapi -I$(srcdir)/../useful

titan_materialnames_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 

noinst_LIBRARIES = libtitan_preprocess.a

libtitan_preprocess_a_SOURCES =  = preprocess.C boundary_preproc.C createfunky.C element_preproc.C \
	hilbert.C node_preproc.C \
	boundary_preproc.h element_preproc.h node_preproc.h \
	../header/FileFormat.h

libtitan_preprocess_a_CXXFLAGS = -I$(srcdir)/../gisapi -I$(srcdir)/../useful

#titan_preprocess_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 


