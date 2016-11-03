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
 * $Id: GisApi.h 214 2009-07-14 22:28:13Z dkumar $
 */

#ifndef GIS_API_H
#define GIS_API_H

#include <string>
using namespace std;

class GisLabels;
class GisLines;

//! structure holding the GIS header
typedef struct
{
  //! minimum z value
  double zmin;

  //! maximum z value
  double zmax;

  //! minimum x value
  double xmin;

  //! maximum x value
  double xmax;

  //! minimum y value
  double ymin;

  //! maximum y value
  double ymax;

  //! GIS cell edge size
  double resolution;

  //! number of rows of GIS cells or points
  int nrows;

  //! number of columns of GIS cells or points
  int ncols;
  double wxmin;
  double wxmax;
  double wymin;
  double wymax;
  double wresolution;
  int wnrows;
  int wncols;

  //! I believe this is a string holding a GIS filename, maybe a complete path.
  char *datafile;

  //! flag saying wether the file is compressed or not
  int compressed;
} Gis_Head;

//! structure holding the GIS terrain elevation data
typedef struct
{
  //! GIS header
  Gis_Head ghead;

  //! array holding elevation data at GIS grid points
  float **elev;

  //! array holding x direction slope data at GIS grid points
  float **xslope;

  //! array holding y direction slope data at GIS grid points
  float **yslope;

  //! array holding slope (actually I think this is a mixed directions second order partial derivative first order x, first order y)) data at GIS grid points
  float **slope;

  //! array holding curvature (actually I think this is second derivative rather than curvature) in the x direction data at GIS grid points
  float **xcurv;

  //! array holding (actually I think this is second derivative rather than curvature) in the y direction data at GIS grid points
  float **ycurv;
} Gis_Grid;

//! structure holding the GIS material map
typedef struct
{
  //! GIS header
  Gis_Head ghead;

  //! number of categories (materials possibly plus 1 for 0 state)
  int ncats;

  //! array holding the category values (which material at each GIS grid point)
  char **cvals;

  //! array of pointers to strings holding the category (material) names
  char **cnames;
} Gis_Raster;

//! structure holding information about a matching GIS image, for the gmfg viewer 
typedef struct
{
  //! GIS header
  Gis_Head ghead;

  //! array of integer (color) values
  unsigned char **ivals;

  //! red ??lower upper triangular??
  unsigned char rlut[256];

  //! green ??lower upper triangular??
  unsigned char glut[256];

  //! blue ??lower upper triangular??
  unsigned char blut[256];
} Gis_Image;

//! structure holding information about GIS vector data (putting roads, bridges, buildings etc on the image map), for the gmfg viewer 
typedef struct
{
  GisLabels *glabels;
  GisLines *glines;
} Gis_Vector;

#define G_API_MAXFLOAT 1.0e31
#define G_API_BIGFLOAT 1.0e32

/***************************************************************/
/* Error codes:
	 0 -> everything okay
	-1 -> resolution finer than available information
	-2 -> requested (x,y) location is not in the data set
	-3 -> memory problem
	-4 -> something else wrong*/
/***************************************************************/

/***************************************************************/
/* SELECTION OF DATA */
/***************************************************************/
//! initialize the GIS digital elevation map
int Initialize_GIS_data (char *GISDbase, char *location, char *mapset,
                         char *raster_file);

//! initialize the GIS material map
int Initialize_Raster_data (char *GISDbase, char *location, char *mapset,
                            char *raster_file);

//! initialize the matching visual image map to cover the terrain (texture wrapping) for the gmfg viewer
int Initialize_Image_data (char *GISDbase, char *location, char *mapset,
                           char *raster_file);

//! initialize the matching vector data (where do you put symbols representing roads, bridges, buildings etc on the map) for the gmfg viewer
int Initialize_Vector_data (char *GISDbase, char *location, char *mapset,
                            char *vector_file);

/* Select mapset at location in database
Input:
	GISDbase - Full path of database
	location - Location
	mapset - Mapset
	raster_file - Raster map name
Output:
Return:
	0 if error, 1 otherwise */

//! deletes the GIS elevation map data
int Delete_GIS_data ();

//! deletes the GIS material map data
int Delete_Raster_data ();

//! deletes the image map (for texture wrapping the terrain in the gmfg viewer) data
int Delete_Image_data ();

//! delets the vector data (where to put symbols on the map representing roads, bridges, buildings, etc in the gmfg viewer)
int Delete_Vector_data ();

/* Clears all GIS_data
Input:
Output:
Return:
	0 if error, 1 otherwise */

/*BASIC INFORMATION RECOVERY */
//! Return number of categories in a raster map Input: none, Output: n_categories - number of categories, Return: 0 if OK, see table otherwise
int Get_raster_categories (int *n_categories);



//! Return number of lines in a vector map, Input: none, Output: n_lines - number of lines, Return: 0 if OK, see table otherwise 
int Get_vector_n_lines (int *n_lines);


//! Return line type for line with index line_index, Input: line_index - index of line, Output: line_type - line type, Return:  0 if OK, see table otherwise 
int Get_vector_line_type (int line_index, int *line_type);


//! Return line label for line with index line_index, Input: line_index - index of line, Output: line_str - line label, Return: 0 if OK, see table otherwise
int Get_vector_line_label (int line_index, string * line_str);

//! Return number of elements for line with index line_index, Input: line_index - index of line, Output: line_size - number of line elements, Return: 0 if OK, see table otherwise 
int Get_vector_line_size (int line_index, int *line_size);

//! Return number of lines in a vector map, Input: line_index - index of line, Output: line_x - line elements x coordinate, line_y - line elements y coordinate, Return: 0 if OK, see table otherwise
int Get_vector_line (int line_index, double *line_x, double *line_y);

//! Return category name the given category number, Input: category_id - category number, Output: category_name - category name, Return: 0 if OK, see table otherwise 
int Get_raster_category_name (int category_id, char *category_name);

//! Return extents of original elevation data grid, Input: resolution - resolution, Output: xmax - maximum X coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_xmax (double resolution, double *xmax);

//! Return extents of original raster (material map) grid, Input: resolution - resolution, Output: xmax - maximum X coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_raster_xmax (double resolution, double *xmax);

//! Return extents of original image (for texture wrapping) grid, Input: resolution - resolution, Output: xmax - maximum X coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_image_xmax (double resolution, double *xmax);

//! Return extents of original elevation map grid, Input: resolution - resolution, Output:      xmin - minimum X coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_xmin (double resolution, double *xmin);

//! Return extents of original raster (material map) grid, Input: resolution - resolution, Output:      xmin - minimum X coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_raster_xmin (double resolution, double *xmin);

//! Return extents of original image (for texture wrapping) grid, Input: resolution - resolution, Output:       xmin - minimum X coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_image_xmin (double resolution, double *xmin);

//! Return extents of original elevatio map grid, Input: resolution - resolution, Output: ymax - maximum Y coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_ymax (double resolution, double *ymax);

//! Return extents of original raster (material map) grid, Input: resolution - resolution, Output: ymax - maximum Y coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_raster_ymax (double resolution, double *ymax);

//! Return extents of original image (for texture wrapping) grid, Input: resolution - resolution, Output: ymax - maximum Y coordinate of original grid, Return: 0 if OK, see table otherwise 
int Get_image_ymax (double resolution, double *ymax);

//! Return extents of original elevation map grid, Input: resolution - resolution, Output: ymin - minimum Y coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_ymin (double resolution, double *ymin);

//! Return extents of original raster (material map) grid, Input: resolution - resolution, Output: ymin - minimum Y coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_raster_ymin (double resolution, double *ymin);

//! Return extents of original image (for texture wrapping) grid, Input: resolution - resolution, Output: ymin - minimum Y coordinate of original grid, Return: 0 if OK, see table otherwise
int Get_image_ymin (double resolution, double *ymin);

//! Return minimum elevation of original grid, Input: resolution - resolution, Output: elevmin - minimum elevation of original grid, Return: 0 if OK, see table otherwise 
int Get_elev_min (double resolution, double *elevmin);

//! Return maximum elevation of original grid, Input: resolution - resolution, Output: elevmax - maximum elevation of original grid, Return: 0 if OK, see table otherwise
int Get_elev_max (double resolution, double *elevmax);


//! Return extents of active window of the elevation map, Input: none, Output: xmin - minimum X coordinate of active region, xmax - maximum X coordinate of active region, ymin - minimum Y coordinate of active region, ymax - maximum Y coordinate of active region, Return: 0 if error, 1 otherwise 
int Get_window (double *xmin, double *xmax, double *ymin, double *ymax);

//! Return extents of active window of the raster (material) map, Input: none, Output: xmin - minimum X coordinate of active region, xmax - maximum X coordinate of active region, ymin - minimum Y coordinate of active region, ymax - maximum Y coordinate of active region, Return: 0 if error, 1 otherwise 
int Get_raster_window (double *xmin, double *xmax, double *ymin,
                       double *ymax);

//! Return extents of active window of the image (for texture wrapping), Input: none, Output: xmin - minimum X coordinate of active region, xmax - maximum X coordinate of active region, ymin - minimum Y coordinate of active region, ymax - maximum Y coordinate of active region, Return: 0 if error, 1 otherwise 
int Get_image_window (double *xmin, double *xmax, double *ymin, double *ymax);

//! Return original resolution of elevation map grid, Input: none, Output: resolution - resolution, Return: 0 if OK, see table otherwise 
int Get_max_resolution (double *resolution);

//! Return original resolution of raster (material) map grid, Input: none, Output: resolution - resolution, Return: 0 if OK, see table otherwise 
int Get_raster_resolution (double *resolution);

//! Return original resolution of image (for texture wrapping) grid, Input: none, Output: resolution - resolution, Return: 0 if OK, see table otherwise 
int Get_image_resolution (double *resolution);


//! Return number of rows of original elevation map grid, Input: none, Output: rows - number of rows of original grid, Return: 0 if error, 1 otherwise
int Get_number_of_rows (int *rows);

//! Return number of rows of original raster (material) map grid, Input: none, Output: rows - number of rows of original grid, Return: 0 if error, 1 otherwise
int Get_raster_nrows (int *rows);

//! Return number of rows of original image (for texture wrapping) grid, Input: none, Output: rows - number of rows of original grid, Return: 0 if error, 1 otherwise
int Get_image_nrows (int *rows);


//! Return number of columns of original elevation map grid, Input: none, Output: cols - number of columns of original grid, Return: 0 if error, 1 otherwise 
int Get_number_of_columns (int *cols);

//! Return number of columns of original raster (material) map grid, Input: none, Output: cols - number of columns of original grid, Return: 0 if error, 1 otherwise 
int Get_raster_ncols (int *cols);

//! Return number of columns of original image (for texture wrapping) grid, Input: none, Output: cols - number of columns of original grid, Return: 0 if error, 1 otherwise 
int Get_image_ncols (int *cols);


/***************************************************************/

/***************************************************************/
/*GETTING VALUES FOR SINGLE POINTS*/
//! Return elevation at point XY of original grid, Input: resolution - resolution, x - point X coordinate, y - Point Y coordinate, Output: elev - elevation at point XY of original grid, Return: 0 if OK, see table otherwise
int Get_elevation (double resolution, double x, double y, double *elev);

//! Return slope at point XY of original grid, Input: resolution - resolution, x - point X coordinate, y - Point Y coordinate, Output: xslope - slope at point XY of original grid in X direction, yslope - slope at point XY of original grid in Y direction, Return: 0 if OK, see table otherwise 
int Get_slope (double resolution, double x, double y, double *xslope,
               double *yslope);

//! Return curvature at point XY of original grid, Input: resolution - resolution, x - point X coordinate, y - Point Y coordinate, Output: xcurv - curvature at point XY of original grid in X direction, ycurv - curvature at point XY of original grid in Y direction, Return: 0 if OK, see table otherwise 
int Get_curvature (double resolution, double x, double y, double *xcurv,
                   double *ycurv);

//! Return category number at point XY of raster map, Input: resolution - resolution, x - point X coordinate, y - Point Y coordinate, Output: category_id - category number at point XY of raster map, Return: 0 if OK, see table otherwise
int Get_raster_id (double resolution, double x, double y, int *category_id);

//! Return RGB at point XY of image, Input: resolution - resolution, x - point X coordinate, y - Point Y coordinate, Output: r - R component at point XY of image, g - G component at point XY of image, b - B component at point XY of image, Return: 0 if OK, see table otherwise
int Get_image (double resolution, double x, double y, unsigned char *r,
               unsigned char *g, unsigned char *b);

/***************************************************************/
/* GETTING VALUES FOR MULTIPLE POINTS */
/***************************************************************/
//! Return elevation at points XY of original grid, Input: resolution - resolution, x - array with points X coordinate, y - array with points Y coordinate, number_of_locations - number of points, Output: elev - elevation at points XY of original grid, Return: 0 if OK, see table otherwise
int Get_elevation_array (double *resolution, double *x, double *y,
                         double *elev, int number_of_locations);

//! Return slope at points XY of original grid, Input: resolution - resolution, x - array with points X coordinate, y - array with points Y coordinate, number_of_locations - number of points, Output: xslope - slopes at points XY of original grid in X direction, yslope - slopes at points XY of original grid in Y direction, Return: 0 if OK, see table otherwise
int Get_slope_array (double *resolution, double *x, double *y, double *xslope,
                     double *yslope, int number_of_locations);

//! Return curvature at points XY of original grid, Input: resolution - resolution, x - array with points X coordinate, y - array with points Y coordinate, number_of_locations - number of points, Output: xcurv - curvature at point XY of original grid in X direction, ycurv - curvature at point XY of original grid in Y direction, Return: 0 if OK, see table otherwise
int Get_curvature_array (double *resolution, double *x, double *y,
                         double *xcurv, double *ycurv,
                         int number_of_locations);

//! Return category number at points XY of raster map, Input: resolution - resolution, x - array with points X coordinate, y - array with points Y coordinate, number_of_locations - number of points, Output: category_id - category number at points XY of raster map, Return: 0 if OK, see table otherwise
int Get_raster_id_array (double *resolution, double *x, double *y,
                         int *category_id, int number_of_locations);

//! Return RGB at point XY of image, Input: resolution - resolution, x - array with points X coordinate, y - array with points Y coordinate, number_of_locations - number of points, Output: r - R component at points XY of image, g - G component at points XY of image, b - B component at points XY of image, Return: 0 if OK, see table otherwise
int Get_image_array (double *resolution, double *x, double *y,
                     unsigned char *r, unsigned char *g, unsigned char *b,
                     int number_of_locations);

//! Return elevation at points XY of original grid, Input: resolution - resolution, xmin - minimum X coordinate of points window, xmax - maximum X coordinate of points window, ymin - minimum Y coordinate of points window, ymax - maximum Y coordinate of points window, Output: elev - elevation at points XY of selected window, Return: 0 if OK, see table otherwise
int Get_elevation_grid (double resolution, double xmin, double xmax,
                        double ymin, double ymax, double *elev);

//! Return elevation at points XY of original grid, Input: resolution - resolution, xmin - minimum X coordinate of points window, xmax - maximum X coordinate of points window, ymin - minimum Y coordinate of points window, ymax - maximum Y coordinate of points window, Output: elev - elevation at points XY of selected window, Return: 0 if OK, see table otherwise
int Get_slope_grid (double resolution, double xmin, double xmax, double ymin,
                    double ymax, double *slope);

//! Return elevation at points XY of original grid, Input: resolution - resolution, xmin - minimum X coordinate of points window, xmax - maximum X coordinate of points window, ymin - minimum Y coordinate of points window, ymax - maximum Y coordinate of points window, Output: elev - elevation at points XY of selected window, Return: 0 if OK, see table otherwise
int Get_curvature_grid (double resolution, double xmin, double xmax,
                        double ymin, double ymax, double *xcurv,
                        double *ycurv);

//! Return category number at points XY of original grid, Input: resolution - resolution, xmin - minimum X coordinate of points window, xmax - maximum X coordinate of points window, ymin - minimum Y coordinate of points window, ymax - maximum Y coordinate of points window, Output: category_id - category number at points XY of selected window, Return: 0 if OK, see table otherwise
int Get_raster_id_grid (double resolution, double xmin, double xmax,
                        double ymin, double ymax, int *category_id);

//! Return RGB at points XY of original grid, Input: resolution - resolution, xmin - minimum X coordinate of points window, xmax - maximum X coordinate of points window, ymin - minimum Y coordinate of points window, ymax - maximum Y coordinate of points window, Output: r - R component at points XY of selected window, g - G component at points XY of selected window, b - B component at points XY of selected window, Return: 0 if OK, see table otherwise
int Get_image_grid (double resolution, double xmin, double xmax, double ymin,
                    double ymax, unsigned char *r, unsigned char *g,
                    unsigned char *b);

/***************************************************************/
/* SETTING GLOBAL PARAMETERS */
/***************************************************************/

//! Set vector scale, used to calculate precision, Input: scale - Scale of the vector data, ex.: 10000, API pre-defines scale at 50000 
void Set_vector_scale (double scale);

#endif
