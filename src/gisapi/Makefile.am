#
# automake template for src/gisapi directory
#
# $Id: Makefile.am 233 2012-03-27 18:30:40Z dkumar $
#

noinst_LIBRARIES = libgisapi.a

libgisapi_a_SOURCES = GisApi.C GisAscFile.C GisBinFile.C GisCats.C \
	GisRasterHdr.C GisColors.C GisLines.C GisLabels.C GisSPRFile.C \
	GisGrid.C GisTriFile.C \
	GisApi.h GisAscFile.h GisBinFile.h GisCats.h GisRasterHdr.h \
	GisColors.h GisLines.h GisLabels.h GisSPRFile.h GisGrid.h GisTriFile.h

if WITH_GDAL
  libgisapi_a_SOURCES += GdalApi.cc gmfg_GdalApi.h
endif
