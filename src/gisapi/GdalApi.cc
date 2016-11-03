#include <iostream>
#include <cstdlib>
using namespace std;

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "GisApi.h"
#include "gmfg_GdalApi.h"

extern Gis_Head gishead;
extern Gis_Grid gis_grid;

int
Initialize_GDAL_data(const char *fullGispath)
{
   GDALDataset *dataset;
   int nx, ny;
   double adfGeoTransform[6];
   Gis_Head gishead;

   // Initialize GDAL
   GDALAllRegister();

   // open GIS dataset
   dataset = (GDALDataset *) GDALOpen (fullGispath, GA_ReadOnly);
   if ( dataset == NULL )
   {
      cout << "GIS ERROR: Unable to open GIS DATA" << endl;
      return (-4);
   }

   // check if the data is in UTM coordinates
   OGRSpatialReference *oSRS = new OGRSpatialReference(NULL);
   char *pszProjection = (char *) dataset->GetProjectionRef();
   if ( oSRS->importFromWkt(&pszProjection) == CE_None )
   {
      int *pbxmax;
      int zone = oSRS->GetUTMZone(pbxmax);
      if ( zone == 0 )
      {
	 cout <<"FATAL ERROR:"<<endl;
	 cout <<"The mapset seems to be in coordinates, other than UTM."<<endl;
         cout <<"TITAN needs the mapset to be in UTM coordinates." <<endl;
	 cout <<"Consider using \"gdalwarp\" to transform."<<endl;
	 return (-4);
      }
   }
   else
   {
      cout <<"FATAL ERROR:"<<endl;
      cout <<"Can't read the projection information from mapset" << endl;
      return (-4);
   }
   nx = dataset->GetRasterXSize();
   ny = dataset->GetRasterYSize();
   if ( nx < 1 || ny < 1)
   {
      cout <<"FATAL ERROR:"<<endl;
      cout <<"GDAL could not read, metadata from the mapset."<<endl;
      return (-4);
   }
   gishead.ncols = nx;
   gishead.nrows = ny;

   // read the window information
   if ( dataset->GetGeoTransform( adfGeoTransform ) == CE_None )
   {
      if ( adfGeoTransform[2] != 0 || adfGeoTransform[4] != 0 )
      {
         cout <<"FATAL ERROR:"<<endl; 
	 cout <<"Raster Map is not North-up."<<endl;
         cout <<"Consider using \"gdalwarp\" to transform"<<endl;
         return (-4);
      }
      gishead.xmin = adfGeoTransform[0];
      gishead.wresolution = fabs(adfGeoTransform[1]);
      gishead.xmax = gishead.xmin + nx * gishead.wresolution;
      gishead.ymax = adfGeoTransform[3];
      gishead.resolution = fabs(adfGeoTransform[5]);
      gishead.ymin = gishead.ymax - gishead.resolution * ny;
      gishead.compressed = 1; // doesn't mean anything if GDAL is bing used
      gis_grid.ghead = gishead;
      gis_grid.ghead.datafile = strdup(fullGispath);
      return 0;
   }
   return -4;
}

int
load_GDAL_data()
{
   GDALDataset *dataset;
   GDALRasterBand  *poBand;
   GDALAllRegister();
   int nrows, ncols;
   const char *fullGispath = gis_grid.ghead.datafile;

   // open GIS dataset
   dataset = (GDALDataset *) GDALOpen (fullGispath, GA_ReadOnly);

   // allocate memory
   nrows = gis_grid.ghead.nrows;
   ncols = gis_grid.ghead.ncols;
   if ( !(gis_grid.elev = alloc_float_matrix (nrows,ncols)))
      return -3; // memory error

   poBand = dataset->GetRasterBand( 1 );
   if( poBand->GetColorTable() != NULL )
   {
       cout<<"Band has a color table with "
             <<poBand->GetColorTable()->GetColorEntryCount()<<" entries."<<endl;
       cout<<"Can't process RGB information, yet."<<endl;
       return -4;
   }

   // check if the raster band size matches with the size in header
   int nXSize = poBand->GetXSize();
   int nYSize = poBand->GetYSize();
   if ( nXSize != ncols || nYSize != nrows );
   {
      cout<<"Size mismatch between, header and raster band data."<<endl;
      return -4;
   }
   // read the data, 
   for (int i=0; i<nrows; i++)
      if (poBand->RasterIO( GF_Read,0, i, ncols, 1, gis_grid.elev[i], 
                        ncols, 1, GDT_Float32, 0, 0 ) != CE_None )
      {
         cout <<"FATAL ERROR:"<<endl;
         cout <<"Failed to read data from mapset."<<endl;
         return -4;
      }
   return 0;
}
