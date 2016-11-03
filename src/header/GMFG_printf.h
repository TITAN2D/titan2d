#ifndef __GMFG_PRINT_H__
#define __GMFG_PRINT_H__

#include <stdarg.h>
#ifdef PRINT_DEBUG_FILES
#  define gmfg_printf fprintf
#else 
#  define gmfg_printf noprint
#endif
inline void noprint(FILE *fname,...) { return; };

#endif
