/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : load_param
 * COMMENT  : Stores the parameter value, minima and maxima at the
 *            required address.  Uses str_to_vals to decode the strings and
 *            store the values. This routine mainly handles the error conditions.
 *            Examples of legal strings for this routine are given in str_to_vals.c
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define LOAD_PARAM_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : load_param
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : Returns 0 if successful, 1 otherwise.
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long load_param (PARAM *param) {

	long i;
	double *dval, *dmin, *dmax, *ddef;
	float *fval, *fmin, *fmax, *fdef;
	int *lval, *lmin, *lmax, *ldef;
	char **sval, **sdef;    // 2016-01-13 PAN: added string pointers

	if (param->type == M_DOUBLE) {
		param->value = (char *)umalloc (param->size * sizeof (double));
		param->def = (char *)umalloc (param->size * sizeof (double));
		param->min = (char *)umalloc (param->size * sizeof (double));
		param->max = (char *)umalloc (param->size * sizeof (double));
	} else if (param->type == M_FLOAT) {
		param->value = (char *)umalloc (param->size * sizeof (float));
		param->def = (char *)umalloc (param->size * sizeof (float));
		param->min = (char *)umalloc (param->size * sizeof (float));
		param->max = (char *)umalloc (param->size * sizeof (float));
	} else if (param->type == M_LONG) {
		param->value = (char *)umalloc (param->size * sizeof (int));
		param->def = (char *)umalloc (param->size * sizeof (int));
		param->min = (char *)umalloc (param->size * sizeof (int));
		param->max = (char *)umalloc (param->size * sizeof (int));
	} else if (param->type == M_STRING) {
		param->value = (char *)umalloc (param->size * sizeof (char *));
		param->def = (char *)umalloc (param->size * sizeof (char *));
		param->min = (char *)umalloc (param->size * sizeof (char *));
		param->max = (char *)umalloc (param->size * sizeof (char *));
	}

/*
* decode minima
*/
	if (param->bound_status == M_BOUNDED) {
		lmin = (int *)(param->min);	
		for (i = 0; i < param->size; i++)
			*lmin++ = 0;
	} else {
		if (str_to_vals (param->min_string, param->size,
									param->type, param->min)) {
			(void)fprintf (stderr, "Parameter is '%s'\n", param->key);
			(void)fprintf (stderr, "Decoding minimum values.\n");
			(void)fprintf (stderr, "Encoded string is:\n'%s'\n", param->min_string);
			return (1);
		}
	}

/*
* decode maxima
*/
	if (param->bound_status == M_BOUNDED) {
		lmax = (int *)(param->max);	
		for (i = 0; i < param->size; i++)
			*lmax++ = (long)(param->bound_dimen->value);
	} else {
		if (str_to_vals (param->max_string, param->size,
									param->type, param->max)) {
			(void)fprintf (stderr,"Parameter is '%s'\n", param->key);
			(void)fprintf (stderr,"Decoding maximum values.\n");
			(void)fprintf (stderr,"Encoded string is:\n'%s'\n",param->max_string);
			return (1);
		}
	}

/*
* decode default values
*/
	if (str_to_vals (param->value_string, param->size, param->type,
								param->def)) {
		(void)fprintf(stderr,"Parameter is '%s'\n", param->key);
		(void)fprintf(stderr,"Decoding default values.\n");
		(void)fprintf(stderr,"Encoded string is:\n'%s'\n",param->value_string);
		return(1);
	}

	switch (param->type) {
		case M_DOUBLE:
			dval = (double *)param->value;
			ddef = (double *)param->def;
			for (i = 0; i < param->size; i++)
				*dval++ = *ddef++;
			break;

		case M_FLOAT:
			fval = (float *)param->value;
			fdef = (float *)param->def;
			for (i = 0; i < param->size; i++)
				*fval++ = *fdef++;
			break;

		case M_LONG:
			lval = (int *)param->value;
			ldef = (int *)param->def;
			for (i = 0; i < param->size; i++)
				*lval++ = *ldef++;
			break;

        // 2016-01-13 PAN: Added case for string parameters
		case M_STRING:
			sval = (char **) param->value;
			sdef = (char **) param->def;
			for (i = 0; i < param->size; i++) {
                *(sval++) = strdup(*(sdef++));
            }
			break;
	}

/*
* check that the defaults lie within the min and max range
*/
	switch (param->type) {

	case M_DOUBLE:

		dval = (double *) param->value;
		dmin = (double *) param->min;
		dmax = (double *) param->max;

		for (i = 0; i < param->size; i++) {

			if (dmin[i] > dmax[i]) {
				(void)fprintf(stderr,
					"ERROR: minimum value exceeds maximum value.\n");
				(void)fprintf(stderr, "Parameter is: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "Default minimum and maximum values are:\nMin: '%s'\nMax: '%s'\n",
				    param->min_string, param->max_string);
//				(void)fprintf(stderr, "The problem is with posn no %ld.\n", i+1);
				(void)fprintf(stderr,
				    "Assigned minimum value = %lf, Specified maximum value = %lf\n", dmin[i], dmax[i]);
				return(1);
			}

			if (dval[i] < dmin[i] || dval[i] > dmax[i]) {
				(void)fprintf(stderr,
					"\nERROR: assigned value is out of range for Parameter: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "       Default: '%s'\n       Minimum: '%s'\n       Maximum: '%s'\n",
				    param->value_string, param->min_string, param->max_string);
				(void)fprintf(stderr, "       Assigned values are:\n");
				(void)fprintf(stderr,
				    "       Value = %lf, Minimum = %lf, Maximum = %lf\n",
				    dval[i], dmin[i], dmax[i]);
				return(1);
			}

		}

		break;

	case M_FLOAT:

		fval = (float *) param->value;
		fmin = (float *) param->min;
		fmax = (float *) param->max;

		for (i = 0; i < param->size; i++) {

			if (fmin[i] > fmax[i]) {
				(void)fprintf(stderr,
					"ERROR: minimum value exceeds maximum value.\n");
				(void)fprintf(stderr, "Parameter is: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "Default minimum and maximum values are:\nMin: '%s'\nMax: '%s'\n",
				    param->min_string, param->max_string);
//				(void)fprintf(stderr, "The problem is with posn no %ld.\n", i+1);
				(void)fprintf(stderr,
				    "Assigned minimum = %f, maximum = %f\n", fmin[i], fmax[i]);
				return(1);
			}

			if (fval[i] < fmin[i] || fval[i] > fmax[i]) {
				(void)fprintf(stderr,
					"\nERROR: assigned value is out of range for Parameter: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "       Default: '%s'\n       Minimum: '%s'\n       Maximum: '%s'\n",
				    param->value_string, param->min_string, param->max_string);
				(void)fprintf(stderr, "       Assigned values are:\n");
				(void)fprintf(stderr,
				    "       Value = %f, Minimum = %f, Maximum = %f\n",
				    fval[i], fmin[i], fmax[i]);
				return(1);
			}

		}

		break;

	case M_LONG:

		lval = (int *) param->value;
		lmin = (int *) param->min;
		lmax = (int *) param->max;

		for (i = 0; i < param->size; i++) {

			if (lmin[i] > lmax[i]) {
				(void)fprintf(stderr,
					"ERROR: minimum value exceeds maximum value.\n");
				(void)fprintf(stderr, "Parameter is: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "Default Minimum and maximum values are:\nMin: '%s'\nMax: '%s'\n",
				    param->min_string, param->max_string);
//				(void)fprintf(stderr, "The problem is with posn no %ld.\n", i+1);
				(void)fprintf(stderr,
				    "Assigned minimum = %d, maximum = %d\n", lmin[i], lmax[i]);
				return(1);
			}

			if (lval[i] < lmin[i] || lval[i] > lmax[i]) {
				(void)fprintf(stderr,
					"\nERROR: assigned value is out of range for Parameter: '%s'\n", param->key);
				(void)fprintf(stderr,
				    "       Default: '%s'\n       Minimum: '%s'\n       Maximum: '%s'\n",
				    param->value_string, param->min_string, param->max_string);
				(void)fprintf(stderr, "       Assigned values are:\n");
				(void)fprintf(stderr,
				    "       Value = %d, Minimum = %d, Maximum = %d\n",
				    lval[i], lmin[i], lmax[i]);
				return(1);
			}
		}
		break;

	case M_STRING:
// Nothing to check
		break;

	}
	return(0);
}
