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
 * $Id: compare_key.C 2 2003-08-13 19:26:11Z sorokine $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"

// compare key compares 2 keys and returns 0 if they are not the same and 1 if they are
int compare_key(unsigned* key1, unsigned* key2)
{
  int i;
 
  for ( i = 0; i < KEYLENGTH; i++ )
    if ( *(key1 + i) != *(key2 + i) )
      return 0;

  
  return 1;
}
