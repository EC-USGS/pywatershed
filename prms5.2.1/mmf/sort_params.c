/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : sort_params
 * COMMENT  : sorts the param array so that the key for each structure
 *            is in increasing alphabetical order
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define SORT_PARAMS_C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : sort_params
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void sort_params (void) {

	PARAM **params;
	PARAM *tmpparam;
	int i, j;

	params =  Mparambase;

/*
**	Make a array of the unsorted parameter order.
*/
	if (!unsort_params) {
		unsort_params = (PARAM **)malloc (Mnparams * sizeof (PARAM *));
		for (i = 0; i < Mnparams; i++)
			unsort_params[i] = params[i];
	}

/*
**	Sort the parameter data base
*/
	for (i = Mnparams-2; i >= 0; i--) {
		for (j =  0; j <= i; j++) {
			if(strcmp(params[j]->key,params[j+1]->key) > 0) {
				tmpparam = params[j];
				params[j] = params[j+1];
				params[j+1] = tmpparam;
			}
		}
	}
}
