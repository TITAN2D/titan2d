/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:
 * Description:
 *
 *******************************************************************
 * $Id: GisApi.C 214 2009-07-14 22:28:13Z dkumar $
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "GisApi.h"
#include "GisBinFile.h"
#include "GisRasterHdr.h"
#include "GisCats.h"
#include "GisColors.h"
#include "GisLabels.h"
#include "GisLines.h"

Gis_Grid gis_grid, gis_grid2;
Gis_Raster gis_rast;
Gis_Image gis_image;
Gis_Vector gis_vector;
double vectorDataScale = 50000.;

int igiscall=0;

int load_GIS_data();
int set_from_header(GisRasterHdr& gisHeader, Gis_Head& aHeadStruct);
int clear_gis_grid();
int clear_gis_rast();
int clear_gis_image();
char **alloc_char_matrix( int nrows, int ncols );
float **alloc_float_matrix( int rows, int cols );
int free_char_matrix( char **m );
int free_float_matrix( float **m );
void get_grid(double resolution, double xmin, double xmax, double ymin, double ymax, float** ingrid, double* outgrid);
void get_int_grid(double resolution, double xmin, double xmax, double ymin, double ymax, char** ingrid, int* outgrid);
void get_rgb_grid(double resolution, double xmin, double xmax, double ymin, double ymax, unsigned char** ingrid, unsigned char* r, unsigned char* g, unsigned char* b);
char **set_cats(GisCats& g_cats);

int print_grid();
int calculate_slope();
int calculate_curvature();
int find_min_max();

/***************************************************************/
/* SELECTION OF DATA */
/***************************************************************/
int Initialize_Vector_data( char* GISDbase, char* location, char* mapset, char* vector_file)
{
  char gisPath[200];
  char gisFullPath[250];
#if defined WIN32
  char* gisSlash = "\\";
#else
  char* gisSlash = "/";
#endif

  gis_vector.glabels = 0;
  gis_vector.glines  = 0;
  
  if ( GISDbase && location && mapset && vector_file )
    {
      strcpy (gisPath, GISDbase);
      sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
      sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%svector%s%s",gisFullPath,gisSlash,vector_file);
      
      GisLabels* gisLabels = new GisLabels (gisFullPath);
      GisLines*  gisLines  = new GisLines  (gisFullPath);
      
      if ( gisLabels->numberOfLabels() > 0 )
	gisLines->setIndices(*gisLabels, vectorDataScale);
      
      gis_vector.glabels = gisLabels;
      gis_vector.glines  = gisLines;
    }
  return 0;
}

int Delete_Vector_data()
{
  if ( gis_vector.glabels )
    delete gis_vector.glabels;
  if ( gis_vector.glines )
    delete gis_vector.glines;
  
  return 0;
}

int Get_vector_n_lines(int* n_lines)
{
  if ( gis_vector.glines )
    {
      *n_lines = gis_vector.glines->numberOfLines();
      return 0;
    }
  return -4;
}

int Get_vector_line_type(int line_index, int* line_type)
{
  if ( gis_vector.glines )
    return gis_vector.glines->getIndex(line_index, line_type);
  
  return -4;
}

int Get_vector_line_label(int line_index, string *line_str)
{
  if ( gis_vector.glines ){
    gis_vector.glines->getLabel(line_index, line_str);
    //    cout<<"in Get_vector_line_label .."<<*line_str<<endl;
    return 0;
  }
  
  return -4;
}

int Get_vector_line_size(int line_index, int* line_size)
{
  if ( gis_vector.glines )
    return gis_vector.glines->getLineSize(line_index, line_size);
  
  return -4;
}

int Get_vector_line(int line_index, double* line_x, double* line_y)
{
  if ( gis_vector.glines )
    return gis_vector.glines->getLine(line_index, line_x, line_y);
  
  return -4;
}

void Set_vector_scale(double scale)
{ vectorDataScale = scale; }

int Initialize_GIS_data (char* GISDbase, char* location, char* mapset, char* raster_file)
{
  int nrows, ncols;
  
  char gisPath[200];
  char gisFullPath[250];
#if defined WIN32
  char* gisSlash = "\\";
#else
  char* gisSlash = "/";
#endif
  
  clear_gis_grid();
  
  if ( GISDbase && location && mapset && raster_file )
    {
      strcpy (gisPath, GISDbase);
      sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
      sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scellhd%s%s",gisFullPath,gisSlash,raster_file);
      
      GisRasterHdr gisHeader (gisFullPath);
      Gis_Head gHeadStruct;
      if ( ( ! gisHeader.good() ) ||
	   ( set_from_header(gisHeader,gHeadStruct) != 0 ) )
	return -4;
      gis_grid.ghead = gHeadStruct;
      
      nrows = gisHeader.Rows();
      ncols = gisHeader.Cols();
      
      if ( nrows < 1 || ncols < 1 )
	return -4;
      
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%sfcell%s%s",gisFullPath,gisSlash,raster_file);
      gis_grid.ghead.datafile = strdup(gisFullPath);
      
      return 0;
    }
  return -4;
}

//Keith Added
int Update_GIS_data (char* GISDbase, char* location, char* mapset, char* raster_file){

  char gisPath[200];
  char gisFullPath[250];
#if defined WIN32
  char* gisSlash = "\\";
#else
  char* gisSlash = "/";
#endif

  if ( GISDbase && location && mapset && raster_file ){
    
    strcpy (gisPath, GISDbase);
    sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
    sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
    strcpy(gisFullPath,gisPath);
    sprintf(gisFullPath,"%scellhd%s%s",gisFullPath,gisSlash,raster_file);
    //printf("\n%s\n",gisFullPath);
    /* this is a nasty hack need to get only whether is compressed or not
       don't need the rest of the header... assumed to be identicle to
       original should really allow it to be different than original but
       that would take more work do it later if needed*/
    GisRasterHdr gisHeader (gisFullPath);

    int nrows=gis_grid.ghead.nrows;
    int ncols=gis_grid.ghead.ncols;

    if ( ! ( gis_grid2.elev = alloc_float_matrix ( nrows, ncols ) ) )
      return -3;	// memory error


    strcpy(gisFullPath,gisPath);
    sprintf(gisFullPath,"%sfcell%s%s",gisFullPath,gisSlash,raster_file);

    GisBinFile binFile (gisFullPath);

    if ( binFile.good() ){
      binFile.setEndian ("big");
      binFile.setDataSize (4);
      binFile.setIsInteger (false);

      binFile.isCompressed (gisHeader.isCompressed());
      binFile.nRows(nrows);
      binFile.nCols(ncols);

      for (int row = 0; row < nrows; row++)
	if ( ! binFile.readRow (row, gis_grid2.elev[row]) )
	  return -4;

      for(int row = 0; row < nrows; row++)
	for(int col = 0; row < ncols; col++)
	  if(fabs(gis_grid2.elev[row][col]-gis_grid.elev[row][col])>1.0)
	    printf("row=%d, col=%d, elevation old=%g new=%g\n",
		   row,col,gis_grid.elev[row][col],gis_grid2.elev[row][col]);

      return 0;
    }
  }
  return -4;
}

int load_GIS_data ()
{
  int row;
  int nrows, ncols;
  char* fullGISDataFilePath = gis_grid.ghead.datafile;
  
  GisBinFile binFile (fullGISDataFilePath,"r");
  if ( binFile.good() )
    {
      binFile.setEndian ("big");
      binFile.setDataSize (4);
      binFile.setIsInteger (false);
      binFile.isCompressed (gis_grid.ghead.compressed==1);
      // 0 = uncompressed, 1 = compressed
      
      nrows = gis_grid.ghead.nrows;
      ncols = gis_grid.ghead.ncols;
      binFile.nRows(nrows);
      binFile.nCols(ncols);
      
      if ( ! ( gis_grid.elev = alloc_float_matrix ( nrows, ncols ) ) )
	return -3;	// memory error
      
      for (row = 0; row < nrows; row++)
	if ( ! binFile.readRow (row, gis_grid.elev[row]) )
	  return -4;
      
      gis_grid.ghead.zmin = G_API_BIGFLOAT;
      gis_grid.ghead.zmax = -G_API_BIGFLOAT;
      return 0;
    }
  return -4;
}

int Initialize_Raster_data (char* GISDbase, char* location, char* mapset, char* raster_file)
{
  int nrows, ncols;
  int row;
  
  char gisPath[200];
  char gisFullPath[250];
#if defined WIN32
  char* gisSlash = "\\";
#else
  char* gisSlash = "/";
#endif
  clear_gis_rast();

  if ( GISDbase && location && mapset && raster_file )
    {
      strcpy (gisPath, GISDbase);
      sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
      sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scellhd%s%s",gisFullPath,gisSlash,raster_file);

      GisRasterHdr gisHeader (gisFullPath);
      Gis_Head gHeadStruct;
      if ( ( ! gisHeader.good() ) ||
	   ( set_from_header(gisHeader,gHeadStruct) != 0 ) )
	return -4;
      gis_rast.ghead = gHeadStruct;
      
      nrows = gisHeader.Rows();
      ncols = gisHeader.Cols();
      
      if ( nrows < 1 || ncols < 1 )
	return -4;
      
      if ( ! ( gis_rast.cvals = alloc_char_matrix ( nrows, ncols ) ) )
	return -3;	// memory error
      
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scats%s%s",gisFullPath,gisSlash,raster_file);
      
      GisCats g_cats (gisFullPath);
      if ( ! g_cats.good() )
	return -4;
      else
	{
	  gis_rast.ncats = g_cats.mumberOfCats();
	  gis_rast.cnames = set_cats(g_cats);
	}
      
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scell%s%s",gisFullPath,gisSlash,raster_file);
      
      GisBinFile binFile (gisFullPath);
      if ( binFile.good() )
	{
	  binFile.setEndian ("big");
	  binFile.setDataSize (gisHeader.cellFormat()+1);
	  binFile.setIsInteger (true);
	  
	  binFile.isCompressed (gisHeader.isCompressed());
	  binFile.nRows(gisHeader.Rows());
	  binFile.nCols(gisHeader.Cols());
	  
	  for (row = 0; row < nrows; row++)
	    if ( ! binFile.readRow (row, gis_rast.cvals[row]) )
	      return -4;
	  return 0;
	}
    }
  return -4;
  
}

int Initialize_Image_data (char* GISDbase, char* location, char* mapset, char* raster_file)
{
  int nrows, ncols;
  int row;
  
  char gisPath[200];
  char gisFullPath[250];
#if defined WIN32
  char* gisSlash = "\\";
#else
  char* gisSlash = "/";
#endif
  
  clear_gis_image();
  
  if ( GISDbase && location && mapset && raster_file )
    {
      strcpy (gisPath, GISDbase);
      sprintf(gisPath,"%s%s%s%s",gisPath,gisSlash,location,gisSlash);
      sprintf(gisPath,"%s%s%s",gisPath,mapset,gisSlash);
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scellhd%s%s",gisFullPath,gisSlash,raster_file);
      
      GisRasterHdr gisHeader (gisFullPath);
      Gis_Head gHeadStruct;
      if ( ( ! gisHeader.good() ) ||
	   ( set_from_header(gisHeader,gHeadStruct) != 0 ) )
	return -4;
      gis_image.ghead = gHeadStruct;
      
      nrows = gisHeader.Rows();
      ncols = gisHeader.Cols();
      
      if ( nrows < 1 || ncols < 1 )
	return -4;

      if ( ! ( gis_image.ivals = (unsigned char**) alloc_char_matrix ( nrows, ncols ) ) )
	return -3;	// memory error

      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scolr%s%s",gisFullPath,gisSlash,raster_file);
      
      GisColors g_colors (gisFullPath);
      if ( ! g_colors.good() )
	return -4;
      else
	{
	  int i;
	  for (i = 0; i < 256; i++)
	    g_colors.getColor(i, gis_image.rlut[i], gis_image.glut[i], gis_image.blut[i] );
	}
      
      strcpy(gisFullPath,gisPath);
      sprintf(gisFullPath,"%scell%s%s",gisFullPath,gisSlash,raster_file);
      
      GisBinFile binFile (gisFullPath);
      if ( binFile.good() )
	{
	  binFile.setEndian ("big");
	  binFile.setDataSize (gisHeader.cellFormat()+1);
	  binFile.setIsInteger (true);
	  
	  binFile.isCompressed (gisHeader.isCompressed());
	  binFile.nRows(gisHeader.Rows());
	  binFile.nCols(gisHeader.Cols());
	  
	  for (row = 0; row < nrows; row++)
	    if ( ! binFile.readRow (row, (char *) gis_image.ivals[row]) )
	      return -4;
	  
	  return 0;
	}
    }
  return -4;
}

int Delete_GIS_data()
{
  clear_gis_grid();
  return 0;
}

int Delete_Raster_data()
{
  clear_gis_rast();
  return 0;
}

int Delete_Image_data()
{
  clear_gis_image();
  return 0;
}

/***************************************************************/
/*BASIC INFORMATION RECOVERY */
/***************************************************************/
int Get_image_xmax(double resolution, double* xmax)
{
  if ( gis_image.ghead.resolution > resolution )
    return -1;
  *xmax = gis_image.ghead.xmax;
  return 0;
}

int Get_image_xmin(double resolution, double* xmin)
{
  if ( gis_image.ghead.resolution > resolution )
    return -1;
  *xmin = gis_image.ghead.xmin;
  return 0;
}

int Get_image_ymax(double resolution, double* ymax)
{
  if ( gis_image.ghead.resolution > resolution )
    return -1;
  *ymax = gis_image.ghead.ymax;
  return 0;
}

int Get_image_ymin(double resolution, double* ymin)
{
  if ( gis_image.ghead.resolution > resolution )
    return -1;
  *ymin = gis_image.ghead.ymin;
  return 0;
}

int Get_raster_categories(int* n_categories)
{
  *n_categories = gis_rast.ncats;
  return 0;
}

int Get_raster_category_name(int category_id, char* category_name)
{
  if (category_id <= gis_rast.ncats || category_id > 0)
    strcpy(category_name,gis_rast.cnames[category_id-1]);
  else
    category_name[0] = 0;
  return 0;
}

int Get_raster_xmax(double resolution, double* xmax)
{
  if ( gis_rast.ghead.resolution > resolution )
    return -1;
  *xmax = gis_rast.ghead.xmax;
  return 0;
}

int Get_raster_xmin(double resolution, double* xmin)
{
  if ( gis_rast.ghead.resolution > resolution )
    return -1;
  *xmin = gis_rast.ghead.xmin;
  return 0;
}

int Get_raster_ymax(double resolution, double* ymax)
{
  if ( gis_rast.ghead.resolution > resolution )
    return -1;
  *ymax = gis_rast.ghead.ymax;
  return 0;
}

int Get_raster_ymin(double resolution, double* ymin)
{
  if ( gis_rast.ghead.resolution > resolution )
    return -1;
  *ymin = gis_rast.ghead.ymin;
  return 0;
}

int Get_xmax(double resolution, double* xmax)
{
  if ( gis_grid.ghead.resolution > resolution )
    return -1;
  *xmax = gis_grid.ghead.xmax;
  return 0;
}

int Get_xmin(double resolution, double* xmin)
{
  if ( gis_grid.ghead.resolution > resolution )
    return -1;
  *xmin = gis_grid.ghead.xmin;
  return 0;
}

int Get_ymax(double resolution, double* ymax)
{
  if ( gis_grid.ghead.resolution > resolution )
    return -1;
  *ymax = gis_grid.ghead.ymax;
  return 0;
}

int Get_ymin(double resolution, double* ymin)
{
  if ( gis_grid.ghead.resolution > resolution )
    return -1;
  *ymin = gis_grid.ghead.ymin;
  return 0;
}

int Get_elev_min(double resolution, double* elevmin)
{
  int status;
  if ( gis_grid.ghead.zmin > gis_grid.ghead.zmax )
    {
      status = find_min_max();
      if (status != 0)
	return status;
    }
  
  *elevmin = gis_grid.ghead.zmin;
  return 0;
}

int Get_elev_max(double resolution, double* elevmax)
{
  int status;
  
  if ( gis_grid.ghead.zmin > gis_grid.ghead.zmax )
    {
      status = find_min_max();
      if (status != 0)
	return status;
    }
  
  *elevmax = gis_grid.ghead.zmax;
  return 0;
}

int Get_image_resolution(double* resolution)
{
  *resolution = gis_image.ghead.resolution;
  return 0;
}

int Get_image_nrows(int *rows)
{
  *rows = gis_image.ghead.nrows;
  return 0;
}

int Get_image_ncols(int *cols)
{
  *cols = gis_image.ghead.ncols;
  return 0;
}

int Get_raster_window(double *xmin, double *xmax, double *ymin, double *ymax)
{
  if ( gis_rast.ghead.xmin < gis_rast.ghead.xmax &&
       gis_rast.ghead.ymin < gis_rast.ghead.ymax )
    {
      *xmin = gis_rast.ghead.xmin;
      *xmax = gis_rast.ghead.xmax;
      *ymin = gis_rast.ghead.ymin;
      *ymax = gis_rast.ghead.ymax;
      return 0;
    }
  return -4;
}

int Get_image_window(double *xmin, double *xmax, double *ymin, double *ymax)
{
  if ( gis_image.ghead.xmin < gis_image.ghead.xmax &&
       gis_image.ghead.ymin < gis_image.ghead.ymax )
    {
      *xmin = gis_image.ghead.xmin;
      *xmax = gis_image.ghead.xmax;
      *ymin = gis_image.ghead.ymin;
      *ymax = gis_image.ghead.ymax;
      return 0;
    }
  return -4;
}

int Get_window(double* xmin, double* xmax, double* ymin, double* ymax)
{
  if ( gis_grid.ghead.xmin < gis_grid.ghead.xmax &&
       gis_grid.ghead.ymin < gis_grid.ghead.ymax )
    {
      *xmin = gis_grid.ghead.xmin;
      *xmax = gis_grid.ghead.xmax;
      *ymin = gis_grid.ghead.ymin;
      *ymax = gis_grid.ghead.ymax;
      return 0;
    }
  return -4;
}

int Get_raster_resolution(double* resolution)
{
  *resolution = gis_rast.ghead.resolution;
  return 0;
}

int Get_raster_nrows(int *rows)
{
  *rows = gis_rast.ghead.nrows;
  return 0;
}

int Get_raster_ncols(int *cols)
{
  *cols = gis_rast.ghead.ncols;
  return 0;
}

int Get_max_resolution(double* resolution)
{
  *resolution = gis_grid.ghead.resolution;
  return 0;
}

int Get_number_of_rows(int *rows)
{
  *rows = gis_grid.ghead.nrows;
  return 0;
}

int Get_number_of_columns(int *cols)
{
  *cols = gis_grid.ghead.ncols;
  return 0;
}

/***************************************************************/
/*GETTING VALUES FOR SINGLE POINTS*/
/***************************************************************/
double interpolate_bilinear_at ( double resolution, double x, double y, float** ingrid )
{
  double dx, dy, dx1, dy1, p1, p2;
  int row, col;

  col = (int)x;
  row = (int)y;

  dx  = x - (double)col;
  dy  = y - (double)row;
  dx1 = 1.0 - dx;
  dy1 = 1.0 - dy;
  p1  = ingrid[row  ][col  ] * dy1 +
	ingrid[row+1][col  ] * dy;
  p2  = ingrid[row  ][col+1] * dy1 +
        ingrid[row+1][col+1] * dy;
  return ( p1 * dx1 + p2 * dx );
}

int Get_image(double resolution, double x, double y, unsigned char* r, unsigned char* g, unsigned char* b)
{
  int row, col;
  unsigned char i;
  
  if ( x >= gis_image.ghead.xmin && x <= gis_image.ghead.xmax &&
       y >= gis_image.ghead.ymin && y <= gis_image.ghead.ymax )
    {
      x = ( x - gis_image.ghead.xmin ) / gis_image.ghead.resolution;
      y = ( gis_image.ghead.ymax - y ) / gis_image.ghead.resolution;
      col = (int)x;
      col = col >= gis_image.ghead.ncols ? gis_image.ghead.ncols-1 : col;
      col = col < 0 ? 0 : col;
      row = (int)y;
      row = row >= gis_image.ghead.nrows ? gis_image.ghead.nrows-1 : row;
      row = row < 0 ? 0 : row;
      i = gis_image.ivals[row][col];
      *r = gis_image.rlut[i];
      *g = gis_image.glut[i];
      *b = gis_image.blut[i];
    }
  return 0;
}

int Get_raster_id(double resolution, double x, double y, int* category_id)
{
  int row, col;
  
  if ( x >= gis_rast.ghead.xmin && x <= gis_rast.ghead.xmax &&
       y >= gis_rast.ghead.ymin && y <= gis_rast.ghead.ymax )
    {
      x = ( x - gis_rast.ghead.xmin ) / gis_rast.ghead.resolution;
      y = ( gis_rast.ghead.ymax - y ) / gis_rast.ghead.resolution;
      col = (int)x;
      col = col >= gis_rast.ghead.ncols ? gis_rast.ghead.ncols-1 : col;
      col = col < 0 ? 0 : col;
      row = (int)y;
      row = row >= gis_rast.ghead.nrows ? gis_rast.ghead.nrows-1 : row;
      row = row < 0 ? 0 : row;
      *category_id = gis_rast.cvals[row][col];
    }
  return 0;
}

int Get_elevation(double resolution, double x, double y, double* elev)
{
  int status;
  
  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  if ( x >= gis_grid.ghead.xmin && x <= gis_grid.ghead.xmax &&
       y >= gis_grid.ghead.ymin && y <= gis_grid.ghead.ymax )
    {
      x = ( x - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution;
      y = ( gis_grid.ghead.ymax - y ) / gis_grid.ghead.resolution;
      if ( (int)y >= gis_grid.ghead.nrows-1 ||
	   (int)x >= gis_grid.ghead.ncols-1 )
	{
	  if((int)y >= gis_grid.ghead.nrows-1) y=gis_grid.ghead.nrows-1;
	  if((int)x >= gis_grid.ghead.ncols-1) x=gis_grid.ghead.ncols-1;
	  *elev = gis_grid.elev[(int)y][(int)x];
	}
      else
	*elev = interpolate_bilinear_at ( gis_grid.ghead.resolution, x, y, gis_grid.elev );
      //		*elev = gis_grid.elev[row][col];
    }
  return 0;
}

int Get_slope(double resolution, double x, double y, double* xslope, double* yslope)
{
  int status;
  
  if ( gis_grid.xslope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  
  if ( x >= gis_grid.ghead.xmin && x <= gis_grid.ghead.xmax &&
       y >= gis_grid.ghead.ymin && y <= gis_grid.ghead.ymax )
    {
      //		row = ( gis_grid.ghead.ymax - y ) / gis_grid.ghead.resolution;
      //		col = ( x - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution;
      //		*xslope = gis_grid.xslope[row][col];
      //		*yslope = gis_grid.yslope[row][col];
      x = ( x - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution;
      y = ( gis_grid.ghead.ymax - y ) / gis_grid.ghead.resolution;
      if ( (int)y >= gis_grid.ghead.nrows-1 ||
	   (int)x >= gis_grid.ghead.ncols-1 )
	{
	  if((int)y >= gis_grid.ghead.nrows-1) y=gis_grid.ghead.nrows-1;
	  if((int)x >= gis_grid.ghead.ncols-1) x=gis_grid.ghead.ncols-1;
	  *xslope = gis_grid.xslope[(int)y-1][(int)x-1];
	  *yslope = gis_grid.yslope[(int)y-1][(int)x-1];
	}
      else
	{
	  *xslope = interpolate_bilinear_at ( gis_grid.ghead.resolution, x, y, gis_grid.xslope );
	  *yslope = interpolate_bilinear_at ( gis_grid.ghead.resolution, x, y, gis_grid.yslope );
	}
      
	}
  return 0;
}

int Get_curvature(double resolution, double x, double y, double* xcurv, double* ycurv)
{
  int row, col, i, j;
  int status;
  float** xxcurv;
  float** yycurv;
  double x1, y1;
  
  xxcurv = alloc_float_matrix(2,2);
  yycurv = alloc_float_matrix(2,2);
  if ( !xxcurv || !yycurv)
    return -3;
  
  if ( gis_grid.xslope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
	}

  if ( x >= gis_grid.ghead.xmin && x <= gis_grid.ghead.xmax &&
       y >= gis_grid.ghead.ymin && y <= gis_grid.ghead.ymax )
    {
      row = (int)(( gis_grid.ghead.ymax - y ) / gis_grid.ghead.resolution);
      if ( row == 0 )
	row++;
      if ( row == gis_grid.ghead.nrows - 1 )
	row--;
      col = (int)(( x - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution);
      if ( col == 0 )
	col++;
      if ( col == gis_grid.ghead.ncols - 1 )
	col--;
      if ( col >= (gis_grid.ghead.ncols - 3) || row >= gis_grid.ghead.nrows - 3 )
	{
	  *xcurv = ( ( gis_grid.xslope[row-1][col+1] - gis_grid.xslope[row-1][col-1] ) +
		     2 * ( gis_grid.xslope[row  ][col+1] - gis_grid.xslope[row  ][col-1] ) +
		     ( gis_grid.xslope[row+1][col+1] - gis_grid.xslope[row+1][col-1] ) ) /
	    ( 8 * gis_grid.ghead.resolution );
	  
	  *ycurv = ( ( gis_grid.yslope[row-1][col-1] - gis_grid.yslope[row+1][col-1] ) +
		     2 * ( gis_grid.yslope[row-1][col  ] - gis_grid.yslope[row+1][col  ] ) +
		     ( gis_grid.yslope[row-1][col+1] - gis_grid.yslope[row+1][col+1] ) ) /
	    ( 8 * gis_grid.ghead.resolution );
	}
      else
	{
	  for (i = 0; i < 2; i++)
	    {
	      for (j = 0; j < 2; j++)
		{
		  xxcurv[i][j] = ( ( gis_grid.xslope[row+i-1][col+j+1] - gis_grid.xslope[row+i-1][col+j-1] ) +
				   2 * ( gis_grid.xslope[row+i][col+j+1] - gis_grid.xslope[row+i][col+j-1] ) +
				   ( gis_grid.xslope[row+i+1][col+j+1] - gis_grid.xslope[row+i+1][col+j-1] ) ) /
		    ( 8 * (float)gis_grid.ghead.resolution );
		  yycurv[i][j] = ( ( gis_grid.yslope[row+i-1][col+j-1] - gis_grid.yslope[row+i+1][col+j-1] ) +
				   2 * ( gis_grid.yslope[row+i-1][col+j] - gis_grid.yslope[row+i+1][col+j] ) +
				   ( gis_grid.yslope[row+i-1][col+j+1] - gis_grid.yslope[row+i+1][col+j+1] ) ) /
		    ( 8 * (float)gis_grid.ghead.resolution );
		}
	    }
	  x1 = x - gis_grid.ghead.xmin - gis_grid.ghead.resolution*(double)col;
	  y1 = gis_grid.ghead.ymax - gis_grid.ghead.resolution*(double)row - y;
	  x1 = x1/gis_grid.ghead.resolution;
	  y1 = y1/gis_grid.ghead.resolution;
	  *xcurv = interpolate_bilinear_at ( gis_grid.ghead.resolution, x1, y1, xxcurv );
	  *ycurv = interpolate_bilinear_at ( gis_grid.ghead.resolution, x1, y1, yycurv );
	  
	}
    }
  free_float_matrix (xxcurv);
  free_float_matrix (yycurv);
  return 0;
}

/***************************************************************/
/* GETTING VALUES FOR MULTIPLE POINTS */
/***************************************************************/
int Get_image_array(double* resolution, double* x, double* y, unsigned char* r, unsigned char* g, unsigned char* b, int number_of_locations)
{
  double x1, y1;
  int row, col;
  int i, j;
  for ( i = 0; i < number_of_locations; i++ )
    {
      if ( x[i] >= gis_image.ghead.xmin && x[i] <= gis_image.ghead.xmax &&
	   y[i] >= gis_image.ghead.ymin && y[i] <= gis_image.ghead.ymax )
	{
	  x1 = ( x[i] - gis_image.ghead.xmin ) / gis_image.ghead.resolution;
	  y1 = ( gis_image.ghead.ymax - y[i] ) / gis_image.ghead.resolution;
	  
	  col = (int)x1;
	  col = col >= gis_image.ghead.ncols ? gis_image.ghead.ncols-1 : col;
	  col = col < 0 ? 0 : col;
	  row = (int)y1;
	  row = row >= gis_image.ghead.nrows ? gis_image.ghead.nrows-1 : row;
	  row = row < 0 ? 0 : row;
	  j = gis_image.ivals[row][col];
	  r[i] = gis_image.rlut[j];
	  g[i] = gis_image.glut[j];
	  b[i] = gis_image.blut[j];
	}
    }
  return 0;
}

int Get_raster_id_array(double* resolution, double* x, double* y, int* category_id, int number_of_locations)
{
  double x1, y1;
  int row, col;
  int i;
  for ( i = 0; i < number_of_locations; i++ )
    {
      if ( x[i] >= gis_rast.ghead.xmin && x[i] <= gis_rast.ghead.xmax &&
	   y[i] >= gis_rast.ghead.ymin && y[i] <= gis_rast.ghead.ymax )
	{
	  x1 = ( x[i] - gis_rast.ghead.xmin ) / gis_rast.ghead.resolution;
	  y1 = ( gis_rast.ghead.ymax - y[i] ) / gis_rast.ghead.resolution;
	  
	  col = (int)x1;
	  col = col >= gis_rast.ghead.ncols ? gis_rast.ghead.ncols-1 : col;
	  col = col < 0 ? 0 : col;
	  row = (int)y1;
	  row = row >= gis_rast.ghead.nrows ? gis_rast.ghead.nrows-1 : row;
	  row = row < 0 ? 0 : row;
	  category_id[i] = gis_rast.cvals[row][col];
	}
    }
  return 0;
}

int Get_elevation_array(double* resolution, double* x, double* y, double* elev, int number_of_locations)
{
  double x1, y1;
  int i;
  int status;
  
  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  for ( i = 0; i < number_of_locations; i++ )
    {
      if ( x[i] >= gis_grid.ghead.xmin && x[i] <= gis_grid.ghead.xmax &&
	   y[i] >= gis_grid.ghead.ymin && y[i] <= gis_grid.ghead.ymax )
	{
	  //			row = ( gis_grid.ghead.ymax - y[i] ) / gis_grid.ghead.resolution;
	  //			col = ( x[i] - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution;
	  //			elev[i] = gis_grid.elev[row][col];
	  x1 = ( x[i] - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution;
	  y1 = ( gis_grid.ghead.ymax - y[i] ) / gis_grid.ghead.resolution;
	  if ( (int)x1 >= gis_grid.ghead.nrows-1 ||
	       (int)y1 >= gis_grid.ghead.ncols-1 )
	    elev[i] = gis_grid.elev[(int)y1-1][(int)x1-1];
	  else
	    elev[i] = interpolate_bilinear_at ( gis_grid.ghead.resolution, x1, y1, gis_grid.elev);
	}
    }
  return 0;
}

int Get_slope_array(double* resolution, double* x, double* y, double* xslope, double* yslope, int number_of_locations)
{
  int row, col;
  int i, status;
  
  if ( gis_grid.xslope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  
  for ( i = 0; i < number_of_locations; i++ )
    {
      if ( x[i] >= gis_grid.ghead.xmin && x[i] <= gis_grid.ghead.xmax &&
	   y[i] >= gis_grid.ghead.ymin && y[i] <= gis_grid.ghead.ymax )
	{
	  row = (int)(( gis_grid.ghead.ymax - y[i] ) / gis_grid.ghead.resolution);
	  col = (int)(( x[i] - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution);
	  xslope[i] = gis_grid.xslope[row][col];
	  yslope[i] = gis_grid.yslope[row][col];
	}
    }
  return 0;
}

int Get_curvature_array(double* resolution, double* x, double* y, double* xcurv, double* ycurv, int number_of_locations)
{
  int row, col;
  int i, status;
  
  if ( gis_grid.slope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  if ( gis_grid.xcurv == 0 )
    {
      status = calculate_curvature();
      if ( status != 0 )
	return status;
    }
  for ( i = 0; i < number_of_locations; i++ )
    {
      if ( x[i] >= gis_grid.ghead.xmin && x[i] <= gis_grid.ghead.xmax &&
	   y[i] >= gis_grid.ghead.ymin && y[i] <= gis_grid.ghead.ymax )
	{
	  row = (int)(( gis_grid.ghead.ymax - y[i] ) / gis_grid.ghead.resolution);
	  col = (int)(( x[i] - gis_grid.ghead.xmin ) / gis_grid.ghead.resolution);
	  xcurv[i] = gis_grid.xcurv[row][col];
	  ycurv[i] = gis_grid.ycurv[row][col];
	}
    }
  return 0;
}

int Get_elevation_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* elev)
{
  int status;
  
  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.elev, elev);
  return 0;
}

int Get_raster_id_grid(double resolution, double xmin, double xmax, double ymin, double ymax, int* rIds)
{
  get_int_grid( resolution, xmin, xmax, ymin, ymax, gis_rast.cvals, rIds);
  return 0;
}

int Get_image_grid(double resolution, double xmin, double xmax, double ymin, double ymax, unsigned char* r, unsigned char* g, unsigned char* b)
{
  get_rgb_grid( resolution, xmin, xmax, ymin, ymax, gis_image.ivals, r, g, b);
  return 0;
}

int Get_slope_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* slope)
{
  int status;
  
  if ( gis_grid.slope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  
  get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.slope, slope);
  return 0;
}

int Get_curvature_grid(double resolution, double xmin, double xmax, double ymin, double ymax, double* xcurv, double* ycurv)
{
  int status;
  
  if ( gis_grid.slope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  
  if ( gis_grid.xcurv == 0 )
    {
      status = calculate_curvature();
      if ( status != 0 )
	return status;
    }
  
  get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.xcurv, xcurv);
  get_grid( resolution, xmin, xmax, ymin, ymax, gis_grid.ycurv, ycurv);
  
  return 0;
}

/***************************************************************/
/* INTERNAL USE FUNCTIONS */
/***************************************************************/

void get_grid(double resolution, double xmin, double xmax, double ymin, double ymax, float** ingrid, double* outgrid)
{
  int row, col, i, j, li;
  int irow, icol;			/* initial row and column */
  int frow, fcol;			/* final row and column */
  double icold, irowd;
  double resr;
  
  icold = (xmin - gis_grid.ghead.xmin) / gis_grid.ghead.resolution;
  icol  = (int) icold;
  
  irowd = (gis_grid.ghead.ymax - ymax) / gis_grid.ghead.resolution;
  irow  = (int) irowd;
  
  fcol  = (int) ( (xmax - gis_grid.ghead.xmin) / gis_grid.ghead.resolution);
  frow  = (int) ( (gis_grid.ghead.ymax - ymin) / gis_grid.ghead.resolution);
  
  resr = resolution/gis_grid.ghead.resolution;
  
  li = 0;
  j = 0;
  row = irow;
  if (fabs(resr - 1.) < 0.0000001)
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      outgrid[li++] = ingrid[row][col++];
	    }
	  row++;
	}
    }
  else
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      outgrid[li++] = ingrid[row][col];
	      col = (int) (icold + (double)(++i) * resr);
	    }
	  row = (int) (irowd  + (double)(++j) * resr);
	}
    }
}

void get_int_grid(double resolution, double xmin, double xmax, double ymin, double ymax, char** ingrid, int* outgrid)
{
  int row, col, i, j, li;
  int irow, icol;			/* initial row and column */
  int frow, fcol;			/* final row and column */
  double icold, irowd;
  double resr;
  
  icold = (xmin - gis_rast.ghead.xmin) / gis_rast.ghead.resolution;
  icol  = (int) icold;
  
  irowd = (gis_rast.ghead.ymax - ymax) / gis_rast.ghead.resolution;
  irow  = (int) irowd;
  
  fcol  = (int) ( (xmax - gis_rast.ghead.xmin) / gis_rast.ghead.resolution);
  frow  = (int) ( (gis_rast.ghead.ymax - ymin) / gis_rast.ghead.resolution);

  resr = resolution/gis_rast.ghead.resolution;
  
  li = 0;
  j = 0;
  row = irow;
  if (fabs(resr - 1.) < 0.0000001)
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      outgrid[li++] = ingrid[row][col++];
	    }
	  row++;
	}
    }
  else
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      outgrid[li++] = ingrid[row][col];
	      col = (int) (icold + (double)(++i) * resr);
	    }
	  row = (int) (irowd  + (double)(++j) * resr);
	}
    }
}

void get_rgb_grid(double resolution, double xmin, double xmax, double ymin, double ymax, unsigned char** ingrid, unsigned char* r, unsigned char* g, unsigned char* b)
{
  int row, col, i, j, li;
  int irow, icol;			/* initial row and column */
  int frow, fcol;			/* final row and column */
  double icold, irowd;
  double resr;
  
  icold = (xmin - gis_image.ghead.xmin) / gis_image.ghead.resolution;
  icol  = (int) icold;
  
  irowd = (gis_image.ghead.ymax - ymax) / gis_image.ghead.resolution;
  irow  = (int) irowd;
  
  fcol  = (int) ( (xmax - gis_image.ghead.xmin) / gis_image.ghead.resolution);
  frow  = (int) ( (gis_image.ghead.ymax - ymin) / gis_image.ghead.resolution);
  
  resr = resolution/gis_image.ghead.resolution;
  
  li = 0;
  j = 0;
  row = irow;
  if (fabs(resr - 1.) < 0.0000001)
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      r[li] = gis_image.rlut[ ingrid[row][col] ];
	      g[li] = gis_image.glut[ ingrid[row][col] ];
	      b[li++] = gis_image.blut[ ingrid[row][col++] ];
	    }
	  row++;
	}
    }
  else
    {
      while ( row < frow )
	{
	  col = icol;
	  i = 0;
	  while ( col < fcol )
	    {
	      r[li] = gis_image.rlut[ ingrid[row][col] ];
	      g[li] = gis_image.glut[ ingrid[row][col] ];
	      b[li++] = gis_image.blut[ ingrid[row][col] ];
	      col = (int) (icold + (double)(++i) * resr);
	    }
	  row = (int) (irowd  + (double)(++j) * resr);
	}
    }
}

int clear_gis_grid()
{
  free(gis_grid.ghead.datafile);
  
  gis_grid.ghead.zmin = 0;
  gis_grid.ghead.zmax = 0;
  gis_grid.ghead.xmin = 0;
  gis_grid.ghead.xmax = 0;
  gis_grid.ghead.ymin = 0;
  gis_grid.ghead.ymax = 0;
  gis_grid.ghead.resolution = 0;
  gis_grid.ghead.wxmin = 0;
  gis_grid.ghead.wxmax = 0;
  gis_grid.ghead.wymin = 0;
  gis_grid.ghead.wymax = 0;
  gis_grid.ghead.wresolution = 0;
  
  free_float_matrix ( gis_grid.elev );
  free_float_matrix ( gis_grid.xslope );
  free_float_matrix ( gis_grid.yslope );
  free_float_matrix ( gis_grid.slope );
  free_float_matrix ( gis_grid.xcurv );
  free_float_matrix ( gis_grid.ycurv );
  
  gis_grid.elev = 0;
  gis_grid.xslope = 0;
  gis_grid.yslope = 0;
  gis_grid.slope = 0;
  gis_grid.xcurv = 0;
  gis_grid.ycurv = 0;

    return 0;
}

int clear_gis_rast()
{
  gis_rast.ghead.zmin = 0;
  gis_rast.ghead.zmax = 0;
  gis_rast.ghead.xmin = 0;
  gis_rast.ghead.xmax = 0;
  gis_rast.ghead.ymin = 0;
  gis_rast.ghead.ymax = 0;
  gis_rast.ghead.resolution = 0;
  gis_rast.ghead.wxmin = 0;
  gis_rast.ghead.wxmax = 0;
  gis_rast.ghead.wymin = 0;
  gis_rast.ghead.wymax = 0;
  gis_rast.ghead.wresolution = 0;
  
  gis_rast.ncats = 0;
  
  free_char_matrix ( gis_rast.cvals );
  free_char_matrix ( gis_rast.cnames );
  
  gis_rast.cvals = 0;
  gis_rast.cnames = 0;
  
  return 0;
}

int clear_gis_image()
{
  gis_image.ghead.zmin = 0;
  gis_image.ghead.zmax = 0;
  gis_image.ghead.xmin = 0;
  gis_image.ghead.xmax = 0;
  gis_image.ghead.ymin = 0;
  gis_image.ghead.ymax = 0;
  gis_image.ghead.resolution = 0;
  gis_image.ghead.wxmin = 0;
  gis_image.ghead.wxmax = 0;
  gis_image.ghead.wymin = 0;
  gis_image.ghead.wymax = 0;
  gis_image.ghead.wresolution = 0;
  
  free_char_matrix ( (char **) gis_image.ivals );
  
  gis_image.ivals =0;
  
  return 0;
}

float **alloc_float_matrix( int nrows, int ncols )
{
  float **m=0;
  int i;
  
  if ( m = (float **) calloc (nrows, sizeof(float *)) )
    {
      if ( m[0] = (float *) calloc (nrows*ncols, sizeof(float)) )
	{
	  for (i = 1; i < nrows; i++)
	    m[i] = m[i-1] + ncols;
	}
      else
	return 0;
      return m;
    }
  return 0;
}

char **alloc_char_matrix( int nrows, int ncols )
{
  char **m=0;
  int i;
  
  if ( m = (char **) calloc (nrows, sizeof(char *)) )
    {
      if ( m[0] = (char *) calloc (nrows*ncols, sizeof(char)) )
	{
	  for (i = 1; i < nrows; i++)
	    m[i] = m[i-1] + ncols;
	}
      else
	return 0;
      return m;
    }
  return 0;
}

int free_float_matrix( float **m )
{
  if ( m )
    {
      if ( m[0] )
	free ( m[0] );
      free (m);
    }
  return 0;
}

int free_char_matrix( char **m )
{
  if ( m )
    {
      if ( m[0] )
	free ( m[0] );
      free (m);
    }
  return 0;
}

int calculate_slope()
{
  int row, col;
  int status;
  
  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  
  free_float_matrix ( gis_grid.xslope );
  free_float_matrix ( gis_grid.yslope );
  free_float_matrix ( gis_grid.slope );
  
  free_float_matrix ( gis_grid.xcurv );
  free_float_matrix ( gis_grid.ycurv );
  
  if ( ! ( gis_grid.xslope = alloc_float_matrix ( gis_grid.ghead.nrows, gis_grid.ghead.ncols ) ) )
    return -3;	/*memory error*/
  
  if ( ! ( gis_grid.yslope = alloc_float_matrix ( gis_grid.ghead.nrows, gis_grid.ghead.ncols ) ) )
    return -3;	/*memory error*/
  
  if ( ! ( gis_grid.slope = alloc_float_matrix ( gis_grid.ghead.nrows, gis_grid.ghead.ncols ) ) )
    return -3;	/*memory error*/
  
  for (row = 1; row < gis_grid.ghead.nrows - 1; row++)
    {
      for (col = 1; col < gis_grid.ghead.ncols - 1; col++)
	{
	  gis_grid.xslope[row][col] = ( ( gis_grid.elev[row-1][col+1] - gis_grid.elev[row-1][col-1] ) +
					2 * ( gis_grid.elev[row  ][col+1] - gis_grid.elev[row  ][col-1] ) +
					( gis_grid.elev[row+1][col+1] - gis_grid.elev[row+1][col-1] ) ) /
	    ( 8 * (float)gis_grid.ghead.resolution );
	  gis_grid.yslope[row][col] = ( ( gis_grid.elev[row-1][col-1] - gis_grid.elev[row+1][col-1] ) +
					2 * ( gis_grid.elev[row-1][col  ] - gis_grid.elev[row+1][col  ] ) +
					( gis_grid.elev[row-1][col+1] - gis_grid.elev[row+1][col+1] ) ) /
	    ( 8 * (float)gis_grid.ghead.resolution );
	  gis_grid.slope[row][col] = (float)sqrt ( gis_grid.xslope[row][col]*gis_grid.xslope[row][col] +
						   gis_grid.yslope[row][col]*gis_grid.yslope[row][col] );
	}
    }
  for (col = 1; col < gis_grid.ghead.ncols - 1; col++)
    {
      gis_grid.xslope[0][col] = gis_grid.xslope[1][col];
      gis_grid.xslope[gis_grid.ghead.nrows-1][col] = gis_grid.xslope[gis_grid.ghead.nrows-2][col];
      gis_grid.yslope[0][col] = gis_grid.yslope[1][col];
      gis_grid.yslope[gis_grid.ghead.nrows-1][col] = gis_grid.yslope[gis_grid.ghead.nrows-2][col];
      gis_grid.slope[0][col] = gis_grid.slope[1][col];
      gis_grid.slope[gis_grid.ghead.nrows-1][col] = gis_grid.slope[gis_grid.ghead.nrows-2][col];
    }
  for (row = 0; row < gis_grid.ghead.nrows; row++)
    {
      gis_grid.xslope[row][0] = gis_grid.xslope[row][1];
      gis_grid.xslope[row][gis_grid.ghead.ncols-1] = gis_grid.xslope[row][gis_grid.ghead.ncols-2];
      gis_grid.yslope[row][0] = gis_grid.yslope[row][1];
      gis_grid.yslope[row][gis_grid.ghead.ncols-1] = gis_grid.yslope[row][gis_grid.ghead.ncols-2];
      gis_grid.slope[row][0] = gis_grid.slope[row][1];
      gis_grid.slope[row][gis_grid.ghead.ncols-1] = gis_grid.slope[row][gis_grid.ghead.ncols-2];
    }
  return 0;
}

int calculate_curvature()
{
  int row, col;
  int status;
  
  if ( gis_grid.xslope == 0 )
    {
      status = calculate_slope();
      if ( status != 0 )
	return status;
    }
  
  free_float_matrix ( gis_grid.xcurv );
  free_float_matrix ( gis_grid.ycurv );
  
  if ( ! ( gis_grid.xcurv = alloc_float_matrix ( gis_grid.ghead.nrows, gis_grid.ghead.ncols ) ) )
    return -3;	/*memory error*/

  if ( ! ( gis_grid.ycurv = alloc_float_matrix ( gis_grid.ghead.nrows, gis_grid.ghead.ncols ) ) )
    return -3;	/*memory error*/

  for (row = 1; row < gis_grid.ghead.nrows - 1; row++)
    {
      for (col = 1; col < gis_grid.ghead.ncols - 1; col++)
	{
	  gis_grid.xcurv[row][col] = ( ( gis_grid.slope[row-1][col+1] - gis_grid.slope[row-1][col-1] ) +
				       2 * ( gis_grid.slope[row  ][col+1] - gis_grid.slope[row  ][col-1] ) +
				       ( gis_grid.slope[row+1][col+1] - gis_grid.slope[row+1][col-1] ) ) /
	    ( 8 * (float)gis_grid.ghead.resolution );
	  gis_grid.ycurv[row][col] = ( ( gis_grid.slope[row-1][col-1] - gis_grid.slope[row+1][col-1] ) +
				       2 * ( gis_grid.slope[row-1][col  ] - gis_grid.slope[row+1][col  ] ) +
				       ( gis_grid.slope[row-1][col+1] - gis_grid.slope[row+1][col+1] ) ) /
	    ( 8 * (float)gis_grid.ghead.resolution );
	}
    }
  for (col = 1; col < gis_grid.ghead.ncols - 1; col++)
    {
      gis_grid.xcurv[0][col] = gis_grid.xcurv[1][col];
      gis_grid.xcurv[gis_grid.ghead.nrows-1][col] = gis_grid.xcurv[gis_grid.ghead.nrows-2][col];
      gis_grid.ycurv[0][col] = gis_grid.ycurv[1][col];
      gis_grid.ycurv[gis_grid.ghead.nrows-1][col] = gis_grid.ycurv[gis_grid.ghead.nrows-2][col];
    }
  for (row = 0; row < gis_grid.ghead.nrows; row++)
    {
      gis_grid.xcurv[row][0] = gis_grid.xcurv[row][1];
      gis_grid.xcurv[row][gis_grid.ghead.ncols-1] = gis_grid.xcurv[row][gis_grid.ghead.ncols-2];
      gis_grid.ycurv[row][0] = gis_grid.ycurv[row][1];
      gis_grid.ycurv[row][gis_grid.ghead.ncols-1] = gis_grid.ycurv[row][gis_grid.ghead.ncols-2];
    }
  return 0;
}

int find_min_max()
{
  int nrows, ncols;
  int row, col;
  int status;
  
  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  nrows = gis_grid.ghead.nrows;
  ncols = gis_grid.ghead.ncols;
  
  for (row = 0; row < nrows; row++)
    {
      for (col = 0; col < ncols; col++)
	{
	  if ( gis_grid.elev[row][col] > gis_grid.ghead.zmax )
	    gis_grid.ghead.zmax = gis_grid.elev[row][col];
	  if ( gis_grid.elev[row][col] < gis_grid.ghead.zmin )
	    gis_grid.ghead.zmin = gis_grid.elev[row][col];
	}
    }
  return 0;
}

int print_grid()
{
  int row, col;
  int status;

  if ( gis_grid.elev == 0 )
    {
      status = load_GIS_data();
      if ( status != 0 )
	return status;
    }
  
  if ( ! gis_grid.elev )
    return -4;
  
  for (row = 0; row < gis_grid.ghead.nrows; row++)
    {
      fprintf (stdout, "row %d\n", row );
      for (col = 0; col < gis_grid.ghead.ncols; col++)
	fprintf (stdout, "%12.8f ", gis_grid.elev[row][col] );
      fprintf (stdout, "\n" );
    }
  fprintf (stdout, "\n" );
  fprintf (stdout, "x=%12.8f y=%12.8f \n", gis_grid.ghead.wxmin, gis_grid.ghead.wymin);
  
  return 0;
}

int set_from_header(GisRasterHdr& gisHeader, Gis_Head& aHeadStruct)
{
  //	Copy the resolutions
  double xres = gisHeader.XRes();
  double yres = gisHeader.YRes();
  if ( fabs(xres- yres)>=xres/10.0 )
    return -4;	// resolution must be the same
  aHeadStruct.wresolution = aHeadStruct.resolution = xres;
  
  //	Copy the edges of the region
  //	Set rows and cols
  aHeadStruct.nrows = aHeadStruct.wnrows = gisHeader.Rows();
  aHeadStruct.ncols = aHeadStruct.wncols = gisHeader.Cols();
  
  aHeadStruct.wxmin = aHeadStruct.xmin = gisHeader.West();
  //	Use this instead of
  //	gis_grid.ghead.wxmax = gis_grid.ghead.xmax = gisHeader.east
  //	to ensure correct ncols
  aHeadStruct.wxmax = aHeadStruct.xmax = gisHeader.West() + gisHeader.Cols() * xres;
  aHeadStruct.wymin = aHeadStruct.ymin = gisHeader.South();
  //	Use this instead of
  //	gis_grid.ghead.wymax = gis_grid.ghead.ymax = gisHeader.north
  //	to ensure correct nrows
  aHeadStruct.wymax = aHeadStruct.ymax = gisHeader.South() + gisHeader.Rows() * yres;
  aHeadStruct.compressed = gisHeader.isCompressed();
  return 0;
}

char **set_cats(GisCats& g_cats)
{
  char **c=0;
  int i;
  
  if ( c = (char **) calloc (g_cats.mumberOfCats(), sizeof(char *)) )
    {
      for (i=1; i<=g_cats.mumberOfCats(); i++)
	{
	  c[i-1] = strdup(g_cats.category(i));
	  //cout << i << ":" << string(c[i-1]) << endl;
	}
      
      return c;
    }
  return 0;
}
