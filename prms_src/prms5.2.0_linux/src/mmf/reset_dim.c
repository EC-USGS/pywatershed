/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : reset_dim
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define RESET_DIM_C 0
#define VALUE_CASE 0
#define MIN_CASE 1
#define MAX_CASE 2
#define NCASES 3

#include <stdio.h>
#include <stdlib.h>
#include "mms.h"

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static void resize_param (PARAM *, long, long, long, long);

/*--------------------------------------------------------------------*\
 | FUNCTION     : reset_dim
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void reset_dim (DIMEN *dim, long nnew) {
	int dimen_used;
	long iparam, ivar, idimen, i, j;
	long size_new, nold;
	long dimen_num;

	long *lptr_max, *lptr_value;

	PARAM *param;
	PUBVAR *var;

	if (nnew == dim->value) return;

/*
**	reset the dimension
*/
	nold = dim->value;
	dim->value = nnew;

/*
**	if existing entry had any index names, free the excess ones, if new size is
**	smaller than previous size. Otherwise, fill the new ones with null strings.
*/

	if (nnew > nold) {
		if (dim->names) {
			dim->names = (char **)realloc ((char **)dim->names,
				nnew * (sizeof(char *)));
			for (i = nold ; i < nnew; i++)
				dim->names[i] = NULL;
		}

		if (dim->notes) {
			dim->notes = (char **)realloc ((char **)dim->notes,
				nnew * (sizeof(char *)));
			for (i = nold ; i < nnew; i++)
				dim->notes[i] = NULL;
		}	

	} else {
		for (i = nnew + 1; i < nold; i++) {
			if (dim->names && dim->names[i]) {
				dim->names[i] = NULL;
			}

			if (dim->notes && dim->notes[i]) {
				dim->notes[i] = NULL;
			}
		}

		if (nnew) {
			if (dim->names)
				dim->names = (char **)realloc ((char **)dim->names,
					nnew * sizeof (char *));
			if (dim->notes)
				dim->notes = (char **)realloc ((char **)dim->notes,
					nnew * sizeof (char *));
		} else {
			dim->names = NULL;
			dim->notes = NULL;
		}
	}

/*
* search through params for parameters which use this dimension
* and resize
*/

	for (iparam = 0; iparam < Mnparams; iparam++) {
		param = Mparambase[iparam];
		dimen_used = FALSE;
		size_new = 1;
		dimen_num = 1;

		for (idimen = 0; idimen < param->ndimen; idimen++) {
			size_new *= param->dimen[idimen]->value;
			if (dim == param->dimen[idimen]) {
				dimen_num = idimen;
				dimen_used = TRUE;
			}
		}
/*
* if this dimension is used by this parameter, resize the parameter
* array
*/
		if (dimen_used) {

/*
* if size_new is zero, set size_new to 1 so that there is at least one
* entry in the parameters data base. This is necesary so that the
* default, maximum and minimum values will be retained for use when
* the size is set to a non-zero value
*/
			if (size_new == 0)
				size_new = 1;

			resize_param (param, dimen_num, nold, nnew, size_new);
			param->size = size_new;
		}
	}

/*
* if a param is bounded by this dimension,
* reset the maximum values accordingly, and the set the current
* values to the maximum if they exceed it
*/
	for (iparam = 0; iparam < Mnparams; iparam++) {
		param = Mparambase[iparam];
		if((param->bound_status == M_BOUNDED) &&
		 					(param->bound_dimen == dim)) {
/*
 (void)fprintf (stderr,"check bound max for %s\n", param->name);
 (void)fprintf (stderr,"   dim = %s;   bound dim = %s\n", dim->name, param->bound_dimen->name);
*/
			lptr_value = (long *) param->value;
			lptr_max = (long *) param->max;

			for (j = 0; j < param->size; j++) {
				lptr_max[j] = dim->value;
/*
 (void)fprintf (stderr,"   j = %d\n", j);
*/
				if (lptr_value[j] > lptr_max[j])
					lptr_value[j] = lptr_max[j];
			}
		}
	}

/*
* search through vars for variables which use this dimension
* and reset size
*/
   for (ivar = 0; ivar < Mnvars; ivar++) {
      var = Mvarbase[ivar];
      if (!(var->private)) {
         dimen_used = FALSE;
         size_new = 1;

         for (idimen = 0; idimen < var->ndimen; idimen++) {
            size_new *= var->dimen[idimen]->value;
            if (dim == var->dimen[idimen]) {
               dimen_num = idimen;
               dimen_used = TRUE;
            }
         }
/*
* if this dimension is used by this variable, resize the variable
*/
         if (dimen_used) {
            var->size = size_new;
         }
      }
   }
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : resize_param
 | COMMENT		: resizes and repacks param array to take account of
 |                  a change in the value of a dimension
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void resize_param (PARAM *param, long dimen_num, long nold, long nnew, long size_new) {

	char *aptr_prev = NULL, *aptr_new = NULL;

	long i, j, iframe, inew, iold, icase;
	long nframes;
	long blocksize;
	long new_framesize, old_framesize;
	long new_index, old_index;
	long *lptr_prev = NULL, *lptr_new = NULL;

	float *fptr_prev = NULL, *fptr_new = NULL;

	double *dptr_prev = NULL, *dptr_new = NULL;

/*
* compute the number of frames
*/

	nframes = 1;

	for (i = dimen_num + 1; i < param->ndimen; i++)
		nframes *= param->dimen[i]->value;

/*
* compute the block size
*/

	blocksize = 1;

	for (i = 0; i < dimen_num; i++)
		blocksize *= param->dimen[i]->value;

/*
* compute the old and new frame sizes
*/

	old_framesize = blocksize * nold;
	new_framesize = blocksize * nnew;

/*
**	resize the value_desc
*/
//	if (size_new)
//		param->value_desc = (char **) realloc (param->value_desc,
//			size_new * sizeof (char *));
//
//	for (i = param->size; i < size_new; i++)
//		param->value_desc[i] = NULL;

/*
* copy the data
*/
	for (icase = 0; icase < NCASES; icase++) {
		switch (icase) {
			case VALUE_CASE:
				aptr_prev = param->value;
				break;

			case MIN_CASE:
				aptr_prev = param->min;
				break;

			case MAX_CASE:
				aptr_prev = param->max;
				break;
		}

		switch (param->type) {
			case M_LONG:
				lptr_prev = (long *) aptr_prev;
				aptr_new = (char *) umalloc (size_new * sizeof(long));
				lptr_new = (long *) aptr_new;
				break;

			case M_FLOAT:
				fptr_prev = (float *) aptr_prev;
				aptr_new = (char *) umalloc (size_new * sizeof(float));
				fptr_new = (float *) aptr_new;
				break;

			case M_DOUBLE:
				dptr_prev = (double *) aptr_prev;
				aptr_new = (char *) umalloc (size_new * sizeof(double));
				dptr_new = (double *) aptr_new;
				break;

		} /* switch (param->type) */

		for (iframe = 0; iframe < nframes; iframe++) {
			for (inew = 0; inew < nnew; inew++) {
				if (inew < nold)
					iold = inew;
				else
					iold = nold - 1;

				for (j = 0; j < blocksize; j++) {
					new_index = j + inew * blocksize + iframe * new_framesize;
					old_index = j + iold * blocksize + iframe * old_framesize;

					switch (param->type) {
						case M_LONG:
							lptr_new[new_index] = lptr_prev[old_index];
							break;

						case M_FLOAT:
							fptr_new[new_index] = fptr_prev[old_index];
							break;

						case M_DOUBLE:
							dptr_new[new_index] = dptr_prev[old_index];
							break;

					} /* switch (param->type) */
				} /* j */
			} /* inew */
		} /* iframe */

		switch (icase) {
			case VALUE_CASE:
				param->value = aptr_new;
				break;

			case MIN_CASE:
				param->min = aptr_new;
				break;

			case MAX_CASE:
				param->max = aptr_new;
				break;

		} /* switch (icase) */
	} /* icase */
}
