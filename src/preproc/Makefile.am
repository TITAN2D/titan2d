#
# automake template for preproc/ directory
#
# $Id: Makefile.am 12 2003-11-07 17:58:49Z kdalbey $
#

bin_PROGRAMS = titan_materialnames titan_preprocess

titan_materialnames_SOURCES = materialnames.C

titan_materialnames_CXXFLAGS = -I../gisapi/ -I../useful/

titan_materialnames_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 

titan_preprocess_SOURCES = preprocess.C boundary.C createfunky.C element2.C \
	hilbert.C node.C \
	boundary.h element.h node.h \
	../header/FileFormat.h

titan_preprocess_CXXFLAGS = -I../gisapi/ -I../useful/

titan_preprocess_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 


