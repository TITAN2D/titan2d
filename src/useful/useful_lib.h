//#ifdef __cplusplus 
//extern "C" 
//{
//#endif

/* multi dimensional memmory allocation routines
   naming conventions
   CAlloc__() the leading C stands for C style indices (they start from 0 
      rather than 1... at one point I also had a list of FAlloc__() funtions
      whose indices started from 1), Alloc stands _ALLOC_ation, the __ is a 
      letter number pair that indicate the type and number of dimensions of 
      the array, specifically...
      U# = #th dimensional _U_nsigned array
      I# = #th dimensional _I_nteger array
      F# = #th dimensional _F_loat array
      D# = #th dimensional _D_ouble array
      C# = #th dimensional _C_har array
   CDeAlloc__() follows the same conventions as CAlloc__(), where DeAlloc
      stands for _DEALLOC_ation */
unsigned   *CAllocU1(int N1);
unsigned  **CAllocU2(int N1, int N2);
unsigned ***CAllocU3(int N1, int N2, int N3);
int        *CAllocI1(int N1);
int       **CAllocI2(int N1, int N2);
int      ***CAllocI3(int N1, int N2, int N3);
float      *CAllocF1(int N1);
float     **CAllocF2(int N1, int N2);
float    ***CAllocF3(int N1, int N2, int N3);
double     *CAllocD1(int N1);
double    **CAllocD2(int N1, int N2);
double   ***CAllocD3(int N1, int N2, int N3);
double  ****CAllocD4(int N1, int N2, int N3, int N4);
char       *CAllocC1(int N1);
char      **CAllocC2(int N1, int N2);
char     ***CAllocC3(int N1, int N2, int N3);

void CDeAllocU1(unsigned   *U1);
void CDeAllocU2(unsigned  **U2);
void CDeAllocU3(unsigned ***U3);
void CDeAllocI1(int        *I1);
void CDeAllocI2(int       **I2);
void CDeAllocI3(int      ***I3);
void CDeAllocF1(float      *F1);
void CDeAllocF2(float     **F2);
void CDeAllocF3(float    ***F3);
void CDeAllocD1(double     *D1);
void CDeAllocD2(double    **D2);
void CDeAllocD3(double   ***D3);
void CDeAllocD4(double  ****D4);
void CDeAllocC1(char       *C1);
void CDeAllocC2(char      **C2);
void CDeAllocC3(char     ***C3);

//string handeling functions

char *allocstrcpy(const char *str);


//input/output to binary binary files
//fopen_bin also buffers the binary file
FILE *fopen_bin(char *filename, char *mode);

void freadU(FILE *fp, unsigned *U);
void freadI(FILE *fp, int    *I);
void freadF(FILE *fp, float  *F);
void freadF2D(FILE *fp, double *D);
void freadD(FILE *fp, double *D);
void freadC(FILE *fp, char   *C);
void freadstring(FILE *fp, char **str);

void fwriteU(FILE *fp, unsigned U);
void fwriteI(FILE *fp, int    I);
void fwriteF(FILE *fp, float  F);
void fwriteD(FILE *fp, double D);
void fwriteC(FILE *fp, char   C);
void fwritestring(FILE *fp, char *str);

//sorting and searching functions
void unique_sort(int *a, int *N);
void unique_sort_d(double *a, int *N);
int searchI1(int *I1, int N1, int findme);
int searchD1(double *D1, int N1, double findme);

//0 to 1 random number generator 
double ran1(long *idum);

//some unions for changing variable types while preserving the original bits


#ifndef UNIONBYTES
#define UNIONBYTES

typedef union fourbytes {
  char c[4];
  unsigned u;
  int i;
  float f;
} FourBytes; 

typedef union eightbytes {
  char c[8];  
  unsigned u[2];
  int i[2];
  float f[2];
  double d;
} EightBytes; 

double key2double(unsigned key[2]);
unsigned* double2key(double D, unsigned key[2]);

//miscellaneous
double sign(double a);

#endif

//miscellaneous
/* how isnan() works, when nan is compared to anything 0 is returned
   any non nan number will either be greater than or equal to or less 
   than zero */
//static inline int isnan(double x){return(!((x>=0.0)||(x<0.0)));};
//static inline int isfinite(double x){return(fabs(x)<HUGE_VAL);};

//#ifdef __cplusplus
//}
//#endif
