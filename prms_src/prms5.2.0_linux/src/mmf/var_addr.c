/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : var_addr
 * COMMENT  : returns a pointer to a PUBVAR struct which contains the given key
 *            returns NULL if key not found
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define VAR_ADDR_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : var_addr
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
PUBVAR *var_addr (char *key) { 
  PUBVAR **vars;
  int lowcomp, midcomp, highcomp;
  long low, mid, high;

  if (Mnvars == 0) return NULL; /* no variables to locate */

  /*
   * get vars from Mvarbase, the global pointer
   */

  vars = Mvarbase;

  /*
   * search between 0 and Mnvars-1
   */

  low = 0;
  high = Mnvars-1;

   lowcomp = strcmp(vars[low]->key, key);
   if (!lowcomp) {
      return vars[low];
   }

   if (lowcomp > 0) return NULL; /* key out of limits */

   highcomp = strcmp(vars[high]->key, key);

   if (!highcomp) {
      return vars[high];
  }

  if (highcomp < 0) return NULL; /* key out of limits */

  /*
   * the basic search uses bisection
   */

  while (low != high) {

    mid = (low + high) / 2;
    midcomp = strcmp(vars[mid]->key, key);

    if (!midcomp) {
       return vars[mid];
    } else {
      if ((mid == low) || (mid == high)) {   /* the search has closed to */
	return NULL;                         /* width 1 without success  */
      } else {
	if (midcomp < 0) {   /* reset low or high as appropriate */
	  low = mid;
	} else {
	  high = mid;
	}
      }
    }

  }

  /* if no match found, return null */

  return NULL;

}

