/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : sort_vars
 * COMMENT  : sorts the pubvar array so that the key for each
 *            structure is in increasing alphabetical order
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#ifdef MALLOC_FUNC_CHECK
#include <malloc_dbg.h>
#endif

#define SORT_VARS_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : sort_vars
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void sort_vars (void) {

  PUBVAR **vars;
  PUBVAR *tmpvar;
  int i, j;


  /*
   * get vars from varbase, the global pointer
   */

  vars =  Mvarbase;

  for (i = Mnvars-2; i >= 0; i--) {

    for (j =  0; j <= i; j++) {

      if(strcmp(vars[j]->key,vars[j+1]->key) > 0) {

	tmpvar = vars[j];
	vars[j] = vars[j+1];
	vars[j+1] = tmpvar;

      }

    }

  }
/*
  printf("sort_vars\n");
  for (i = 0; i < Mnvars; i++) {
          printf("I: %ld %s\n",i,vars[i]->key);
      }
*/

}

