/* define some constants that are used in multiple files */
#define BSFC_NO_CUT 0
#define BSFC_CUT 1
#define BSFC_NOT_BALANCED 1
#define BSFC_BALANCED 0
#define BSFC_COARSE_LEVEL_FLAG 2

// new_old flags...
#define BSFC_NEW 1 //either not checked yet or no constrained node (1 is used for 1 element)
#define BSFC_OLD 0 //has a constrained node which is accounted for by another element


struct sfc_vertex {         /* BSFC vertex  */
  int destination_proc;
  int cut_bin_flag;  /* =BSFC_NO_CUT (0) if this object does not belong to
			a bin with a cut in it
			=BSFC_CUT (1) if this object belongs to a bin with
			a cut in it */
  int next_sfc_vert_index; /* used for creating a linked list of vertices 
			      (linklist is fortran style) */
  unsigned sfc_key[KEYLENGTH];   /* space-filling curve key */
  unsigned my_bin;
  float lb_weight;
  //int count; // used to put the destination proc back into the element hashtable
  
};
typedef struct sfc_vertex BSFC_VERTEX;
typedef struct sfc_vertex * BSFC_VERTEX_PTR; 

struct unstructured_communication {
  int* send_procs_ptr;
  int* recv_procs_ptr;
  int  send_count;
  int  recv_count;
  sfc_vertex* recv_sfc_vert;
  int used_flag;
};

/* declare functions in the sfc routines */
void BSFC_update_element_proc(int, Element*, HashTable*, HashTable*, int);

//! this function figures out how to "bunch" together elements that cannot be put on different processors because of a constrained node, B stands for bunch SFC stands for space filling curve 
void BSFC_combine_elements(int side, Element *EmTemp, 
			   HashTable *HT_Elem_Ptr, HashTable *HT_Node_Ptr, 
			   int destination_proc);

//void BSFC_combine_elements(int, Element*, HashTable*, HashTable*);

int BSFC_find_imbalance(float* work_percent_array, 
			float cumulative_work, 
			float total_work,
			int which_proc, 
			int numprocs);

void BSFC_single_wgt_calc_partition(int wgt_dim, float work_allocated,
				    float* total_weight_array, 
				    int* bin_proc_array, 
				    float* binned_weight_array, 
				    float* work_percent_array, 
				    float* actual_work_allocated,
				    int number_of_bins, int* number_of_cuts,
				    int current_proc, int level_flag, int*);

int BSFC_get_array_location(int number_of_bins, int number_of_bits, 
			    int prev_used_bits, BSFC_VERTEX_PTR sfc_vert_ptr);

void BSFC_get_normed_coords(double min_bounding_box[], 
			    double max_bounding_box[], 
			    double normed_coords[],
			    int num_dims, double my_coords[]);

void BSFC_create_info(double min_bounding_box[], 
		      double max_bounding_box[], int num_dims,
		      int num_local_objects, int wgt_dim, 
		      BSFC_VERTEX_PTR sfc_vert_ptr, 
		      double* coords);

void BSFC_refine_partition(int* local_balanced_flag, 
			   int *amount_of_used_bits, int num_vert_in_cut,
			   BSFC_VERTEX_PTR vert_in_cut_ptr, 
			   float* work_percent_array, float total_weight,
			   float* global_actual_work_allocated,
			   int number_of_cuts, 
			   int* ll_bins_head, float* work_prev_allocated,
			   int subbins_per_bin, int* local_balanced_flag_array,
			   int myid, int numprocs);

int BSFC_create_compare_key(unsigned sfc_key[], unsigned compare_key[], 
			    unsigned AND_operator_array[], int prev_used_bits);

int BSFC_check_refine(unsigned* sfc_key, unsigned* compare_key,
		      unsigned* AND_operator_array);

void BSFC_update_and_send_elements(int myid, int numprocs, 
				   HashTable* HT_Elem_Ptr, HashTable* HT_Node_Ptr, int);

int BSFC_pow(int intbase, int intexp);
