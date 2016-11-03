/* right now this file is used by the preprocessor (initial grid generator)
   and by ../main/datread.C  Currently the hpfem must be recompiled to 
   switch between the different types (text, binary, and coming soon... HDF)
   of initial grid files.  It is desirable to change the preprocessor so
   that the different types of initial grid files have different extensions.
   This would allow datread to DETECT the file type, and (with if statements)
   load it automatically (without recompiling) */

/* choose to pass funky directly or write an intermediate binary or text file
   note that you can choose to write the intermediate files in addition passing
   funky directly... this is to aid in debugging */
#define PASSFUNKY
//#define WRITEFUNKYBIN
//#define WRITEFUNKYDAT

/* automatic choice of which (if any) intermediate files to read depends on
   which intermediate files you chose to write... if you chose not to pass
   funky directly and not to write any, it tries to read an OLD text file
   named funky.dat */
#ifndef PASSFUNKY
#ifdef WRITEFUNKYBIN
#define READFUNKYBIN
#else
#define READFUNKYDAT
#endif
#endif

/* choose what format to write/read hpfem input files funkyxxxx.inp in, your
   choices are: 
     HDF (Hierarchical Data Format... this is the smallest) not yet implemnted
     BIN (binary... in case HDF library is not installed)
     TXT (text... this is the largest but it's useful in debugging) 
   Note that in all three cases the data file is named funkyxxxx.inp (where
   xxxx is the process id).  If more than one is defined precedence is given 
   to HDF then BIN then TXT.  
   Note that when neither HDFINPUT nor BININPUT is defined the default is to
   use text input (it is not necessary to define TXTINPUT)
*/
//#define HDFINPUT //not yet implemented
#define BININPUT
//#define TXTINPUT 

// options specific to BIN (binary) files
//#define WRITEDOUBLEASFLOAT //this can affect the solution but reduces file size
//#define DEBUGFUNKYBIN
//#define DEBUGBIN2
