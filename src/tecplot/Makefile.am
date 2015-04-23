#
# automake template for src/tecplot directory
#
# $Id: Makefile.am 233 2012-03-27 18:30:40Z dkumar $
#

noinst_LIBRARIES = libtecplot.a

libtecplot_a_SOURCES = outsum.C tecplot.C grassout.C

INCLUDES=-I../header/ -I../gisapi

if BUILD_HDF5
    libtecplot_a_SOURCES += xdfm_write.C hd5calls.c
endif
