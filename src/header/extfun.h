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
 * $Id: extfun.h 211 2009-06-16 20:02:10Z dkumar $ 
 */

#ifndef EXTFUN_H
#define EXTFUN_H
#define DEBUG_HEADER
const int QUADNODES = 9;

//! this function checks for any and all possible mesh errors, i.e. it checks if the mesh is legal, it says nothing about the quality of a legal mesh, you must have ghost information present before performing this check, WARNING THIS CHECK TAKES A LOT OF TIME, ONLY USE IT TO DEBUG.
void AssertMeshErrorFree(HashTable *El_Table, HashTable* NodeTable,
				  int numprocs, int myid,double loc);

//! investigate an Element, question his "friends and family" about him.
void ElemBackgroundCheck(HashTable* El_Table, HashTable* NodeTable,
			 unsigned *debugkey, FILE *fp);

void ElemBackgroundCheck2(HashTable* El_Table, HashTable* NodeTable,
			  void *EmDebug, FILE *fp);


//! investigate a Node question his "friends and family" about him.
void NodeBackgroundCheck(HashTable* El_Table, HashTable* NodeTable,
			 unsigned *debugkey, FILE *fp);

//! this function loops through all the elements on this processor and (by calling other functions) checks which elements satisfy criteria for being okay to unrefine, if they can be it unrefines them.
void unrefine(HashTable* El_Table, HashTable* NodeTable, double target, int myid, int nump, TimeProps* timeprops_ptr, MatProps* matprops_ptr);

void delete_oldsons(HashTable* El_Table, HashTable* NodeTable,
		    int myid, void *EmFather);

void refine_neigh_update(HashTable* El_Table, HashTable* NodeTable, 
			 int numprocs, int myid, void* RefinedList,
			 TimeProps* timeprops_ptr);

void unrefine_neigh_update(HashTable* El_Table, HashTable* NodeTable, int myid, void* NewFatherList);


//void unrefine_neigh_update(HashTable* El_Table, int myid, int NumNewFathers, Element** NewFatherList);

void  unrefine_interp_neigh_update(HashTable* El_Table, HashTable* NodeTable,
				   int nump, int myid, void* OtherProcUpdate);
//void  unrefine_interp_neigh_update(HashTable* El_Table, HashTable* NodeTable,				   int nump, int myid, int NumOtherProcUpdate, 				   Element **OtherProcUpdate);

//! only used in debugging
int IfMissingElem(HashTable* El_Table, int myid, int iter, int isearch);

//! only used in debugging
void InsanityCheck(HashTable* El_Table, int nump, int myid, 
		   TimeProps *timeprops_ptr);

//! this function implements 1 time step which consists of (by calling other functions) computing spatial derivatives of state variables, computing k active/passive and wave speeds and therefore timestep size, does a finite difference predictor step, followed by a finite volume corrector step, and lastly computing statistics from the current timestep's data.
void step(HashTable* El_Table, HashTable* NodeTable, int myid, int nump,
	  MatProps* matprops_ptr, TimeProps* timeprops_ptr, 
	  PileProps *pileprops_ptr, FluxProps *fluxprops, 
	  StatProps* statprops_ptr, int* order_flag, 
	  OutLine* outline_ptr, DISCHARGE* discharge, int adaptflag);

//! this function deletes unused elements and nodes, this is called durring grid adaptation in hadpt.C
void delete_unused_elements_nodes(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
				  int myid);
//! this function is used for dynamic replacement of a GIS digital elevation map, for example, when you're running a large simulation on multiple processors about a real life scenario you're expecting to occur soon and a channel collapses.  If in the simulation the flow has not yet a gotten to the channel this could allow you to replace the DEM with another one in which the channel has collapsed and continue running the simulation without having to restart from the beginning.  while this code updates the map, the logic and external programs used to in real time decide if the map should be updated have not been fully developed, this is the damd (data manager daemon, area of research in which Dr. Matt Jones is Participating).
int update_topo(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		int myid, int nump, MatProps* matprops, TimeProps* timeprops, 
		MapNames* mapnames);

/* this before the addition of the restart capability... 
   it has now been split into Read_data() and Read_grid() below
extern void Read_data(int* Node_Num, HashTable** NodeTable, int myid, 
		      int numprocs, int* Elem_Num, HashTable** ElemTable, 
		      MatProps* matprops_ptr, TimeProps* timeprops_ptr, 
		      MapNames* mapnames_ptr, int* adaptflag_ptr, 
		      int* viz_flag_ptr, int* order_flag_ptr, 
		      StatProps* statprops_ptr, PileProps* pileprops_ptr,
		      OutLine* outline_ptr, DISCHARGE* discharge_ptr);
*/

//! this function reads in the input data (excluding the "funky" grid) at the start of a run, whether or not run is a restart.
extern void Read_data(int imat, MatProps* matprops_ptr, PileProps* pileprops_ptr, 
		      StatProps* statprops_ptr, TimeProps* timeprops_ptr, FluxProps *fluxprops, 
		      int* adaptflag_ptr, int* viz_flag_ptr, 
		      int* order_flag_ptr, MapNames* mapnames_ptr, 
		      DISCHARGE* discharge_ptr, OutLine* outline_ptr,
                      int *srctype);

//! this function reads in the "funky" grid at the start of an original run but not during restart.  This used to be part of Read_data() before Keith seperated them when adding the restart capability.  It is my (Keith's) opinion that this should be torn out and along with the preprocessor rewritten into a new format that is a lot more like what happens during the restart, this would significantly reduce the startup time for large runs.
extern void Read_grid(int myid, int numprocs, 
		      HashTable** NodeTable, HashTable** ElemTable, 
		      MatProps* matprops_ptr, OutLine* outline_ptr);

//! this function loads the restart file, recreates the hashtables and restores the saved nodes and elements.  Only one readstatement per Node is performed and one or two per Element depending upon the Element's boundary conditions so it is very fast.  Keith, who wrote this, believes a slightly cleaner solution is to add/move functionality to useful_lib.h and useful_lib.C to pack/unpack variables into an unsigned array, which is what should be done if Read_grid is ever rewritten.
extern int loadrun(int myid, int numprocs, HashTable** NodeTable, 
		   HashTable** ElemTable, MatProps* matprops_ptr,  
		   TimeProps* timeprops_ptr, MapNames *mapnames_ptr, 
		   int *adaptflag_ptr, int *order_flag_ptr, 
		   StatProps* statprops_ptr, DISCHARGE* discharge_ptr,
		   OutLine* outline_ptr);

//! this function writes a restart file, all the non Node, non Element data, thfor example material properties, statistics, and hastable information.  A loop through the hastables call member functions Node::save_node() and Element::save_elem() which each save 1 Node or Element in a single write statement to the restart file so this is VERY fast.  However it could be rewritten in a slightly cleaner fashion by adding/moving functionality to useful_lib.h and useful_lib.C to pack/unpack variables into an unsigned array, which is what should be done if Read_grid is ever rewritten.
extern void saverun(HashTable** NodeTable, int myid, int numprocs, 
		    HashTable** ElemTable, MatProps* matprops_ptr,  
		    TimeProps* timeprops_ptr, MapNames *mapnames_ptr, 
		    int adaptflag, int order_flag, 
		    StatProps *statprops_ptr, DISCHARGE *discharge_ptr,
		    OutLine* outline_ptr, int *savefileflag);

//! this function intializes the piles, by commenting/uncommenting define statements you can switch from parabaloid to elliptical cylinder shaped piles, or even a hard coded pileshapes written to match particular experiments.  Adaptive remeshing and pile reinitialization helps detect small piles and refine around pile edges to obtain a more accurate initial solution and speed up the first few timesteps before adaptive refinement and unrefinement would otherwise occur.  
extern void init_piles(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		       int myid, int numprocs, int adaptflag,
		       MatProps* matprops, TimeProps* timeprops, 
		       MapNames* mapnames, PileProps* pileprops, 
		       FluxProps* fluxprops, StatProps* statprops);


//! this function performs adaptive refinement at timestep zero for refining initial piles and whenever a flux source is activated.
void  initial_H_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int h_count, 
		      MatProps* matprops_ptr, PileProps *pileprops_ptr,
		      FluxProps *fluxprops_ptr, TimeProps* timeprops_ptr,
		      int num_buffer_layer);

//! this function refines all elements whos generation is less than refinelevel, until they are of generation refinelevel and then places the flux sources and, if it is at timestep zero, initial piles.
void  H_adapt_to_level(HashTable* El_Table, HashTable* NodeTable, 
		       MatProps* matprops_ptr, PileProps* pileprops_ptr,
		       FluxProps *fluxprops_ptr, TimeProps* timeprops_ptr, 
		       int refinelevel);

//! this is the normal grid adaptive refinement function it also refreshes the flux sources
void  H_adapt(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int h_count, 
	      double target, MatProps* matprops_ptr, FluxProps *fluxprops_ptr,
	      TimeProps* timeprops_ptr, int num_buffer_layer);


//! this function flushes the hashtables, it is called during grid adaptation in hadpt.C
extern void  htflush(HashTable*, HashTable*, int);

//! Pack_element() is a friend function of the Element and Node classes that packs relevant information from an element "sendel" into a smaller data structure (ElemPack) to be sent by an mpi call to another processor, this is used when exchanging ghost cell information or repartitioning.
extern void Pack_element(void *sendel, ElemPack *elem,
			 HashTable* HT_Node_Ptr, int destination_proc);
//extern void Pack_element(Element* sendel, ElemPack** elemptr, HashTable* HT_Node_Ptr, int s_f);

//! Create new MPI datatype: ElemPack type definition in struct.h so structures of ElemPack and NeighborPack can be sent and received 
extern void MPI_New_Datatype();

//! this function repartitions (redistributes) the number of elements on each processor so they all have approximately the same ammount of work to do.  it is called in hpfem.C and init_piles.C 
extern void repartition(HashTable*, HashTable*, int);
//extern void repartition(HashTable*, HashTable*);

//! the replacement for repartition(), this function repartitions (redistributes) the number of elements on each processor so they all have approximately the same ammount of work to do
extern void repartition2(HashTable* El_Table, HashTable* NodeTable,
			 TimeProps* timeprops_ptr);

//! this function creates the elements listed in recv_array and adds them to the Element HashTable, it will fail an assertion if you tell it to create an Element that already exists, it is called by repartion2(), which requires the deletion of ghost elements first.
extern void IncorporateNewElements(HashTable* El_Table, HashTable* NodeTable,
				   int myid, int num_recv, 
				   ElemPack *recv_array, 
				   TimeProps* timeprops_ptr);

//! quicksort into ascending order, according to matching double precision numbers, the array of pointers to data
extern void q_sort_data(double *numbers, void **data, int left, int right);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void smooth(HashTable*, HashTable*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void smooth_II(HashTable*, HashTable*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void Delete_Table(HashTable*, HashTable*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void all_check(HashTable* eltab, HashTable* ndtab, int myid, int m);

//! at every output interval this function writes a line of revelant statistics about the flow to a file named output_summary.######, among the information contained in the file is which times of output match with which timestep numbers (which are used in the file names of the tecplot output files).  Keith wrote this.
extern void output_summary(TimeProps* timeprops, StatProps* statprops, 
			   int savefileflag);

//! Adam Stinton wanted some particular outputs, such as coordinates, speed, an pile height at locations of pile centroid, maximum pileheight and maximum flow speed, so I (Keith) wrote this for him.  The function writes to a file named "flow_dynamics.stats" (Adam's choice of file name).  The information is useful enough that I included it in the mainline version of the code but we are getting too many output files, it's about time to combine a bunch of them into one file (prefereably output_summary.######)
void OUTPUT_ADAM_STATS(HashTable* El_Table, MatProps* matprops_ptr,
		       TimeProps* timeprops_ptr, StatProps* statprops_ptr);

//! titan now has the capability to caclulate the flux through discharge planes, the sign of the flux follows a righthand rule convention, net v cross ds is non negative for a series os segments "ds" drawn counter clockwise surrounding the pile(s) (out of the box is positive into the box in negative), Keith wrote this
extern void output_discharge(MatProps* matprops, TimeProps* timeprops, 
		      DISCHARGE* discharge, int myid);

//! this function outputs statistics for with a collection of probabilistic runs, this output should probably be combined with the finalstats.###### file, Keith wrote this
extern void output_stoch_stats(MatProps* matprops, StatProps* statprops);

//! this function writes text tecplot output files in the tecplxxxxxxxx.plt format.  Keith rewrote this function to eliminate a lot of bugs shortly after he started the project, since then files were 1/3 the size they were previously.
extern void tecplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		       MatProps* matprops, TimeProps* timeprops, 
		       MapNames* mapnames, double v_star);

//! this function writes text tecplot output files in the mshplxxxxxxxx.plt format.  This is largely untouched since before I (Keith) joined the GMFG, just minor changes.  This is the preferred (by Professor Patra) format of output for debugging purposes even though tecplxxxxxxxx.plt create nicer images.
extern void meshplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
			MatProps* matprops, TimeProps* timeprops, 
			MapNames* mapnames, double v_star);

//! one of "pady's" output functions, since Keith never met "pady" this is probably long out of date
extern void vizplotter(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, 
		       MatProps* matprops, TimeProps* timeprops);

//! another of "pady's" output functions, since Keith never met "pady" this is probably long out of date
void viz_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		       int myid, int numprocs, MatProps* matprops, 
		       TimeProps* timeprops, MapNames* mapnames);
#if HAVE_LIBHDF5 == 1 && HAVE_LIBHDF5_HL == 1
//! the hdf output never worked right (never got debugged)
void hdfviz_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int time_step, double time, int myid, int numprocs, MatProps*);
//! the hdf output never worked right (never got debugged)
void hdf_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int time_step, double time, int myid, int numprocs, MatProps*);
#endif

//! eXtensible Data Model and Format (http://www.arl.hpc.mil/ice/) is a Paraview readable data format 
int write_xdmf(HashTable *El_Table, HashTable *NodeTable,TimeProps *timeprops_ptr,
      MatProps *matprops_ptr, MapNames *mapnames, const int mode);

//! this is Amrita's output function, Keith wrote it to her specifications, it is for use with the gmfg viewer, which Daniel rewrote during the summer of 2006 to remove a lot dependencies and use basically only opengl calls.  This makes the viewer compatible with almost every linux machine and a lot easier (read as possible) to install
void incr_tri_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		     int myid, int numprocs, MatProps* matprops, 
		     TimeProps* timeprops, double v_star);

//! this is one of the small output file size functions, this output format is meant for use with a web browser based viewer, this is now out of date.
void web_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
		int myid, double time ,int numprocs,	MatProps* matprops, 
		TimeProps* timeprops);
//! this is one of the small output file size functions, this output format is meant for use with a web browser based viewer, this is now out of date.
void web_simplify(TimeProps* timeprops);
//! this is one of the small output file size functions, this output format is meant for use with a web browser based viewer, this is now out of date.
void web_correct(TimeProps* timeprops);

//! this function writes the header for grass sites style output, the grass sites output format is correct and works, but importing the data into GIS packages such as ARCGIS is non trivial until you know the trick to it.  Keith does not know the trick, I just wrote this in the format Alex Sorokine specified.
void grass_sites_header_output(TimeProps* timeprops);

//! this function writes one processors grass sites style output, the grass sites output format is correct and works, but importing the data into GIS packages such as ARCGIS is non trivial until you know the trick to it.  Keith does not know the trick, I (Keith) just wrote this in the format Alex Sorokine specified.
void grass_sites_proc_output(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr,
			     int myid, MatProps* matprops, 
			     TimeProps* timeprops);


//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void check_p_order(HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern int howmanyelements(HashTable*);	

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern double random_error(HashTable*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void search_object ( HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, unsigned* key , int myid );

//! this function compares 2 hashtable keys and says if they're the same or different.
extern int compare_key(unsigned*, unsigned*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void list_elements(HashTable*, HashTable*, int, int);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void get_max_key(HashTable*, HashTable*);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void objects_per_slot(HashTable*, HashTable*, int, int);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void fact_matrix(double*, double*, int, int);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern void Assemble_bubble(int, int, int, int, double*, int*, 
		     HashTable*, HashTable*, int);

//! this function is legacy, it is not defined in the finite difference/volume version of titan
extern int make_block(int, int**, int**, int**, int**, int*,int*,
		HashTable*, HashTable*, int, int, int, int);

#endif
#ifdef DEBUG_HEADER  //debug stuff is in SGI/SUNOS fortran interface 
extern void view_elm(HashTable* El_Table, HashTable* NodeTable, int myid);
extern void Mat_write(int, double*, int, int, int);
extern void Mat_write_f(int, double*, int, int, int);
extern "C" void mat_view_(double*, int*, int*);
extern void mat_restore(int, int, int*, int*, int*, int*, double*);
#endif
