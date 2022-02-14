/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : control_array - generic, returns (char *) as a generic pointer
 *            control_larray - returns long *
 *            control_farray - returns float *
 *            control_darray - returns double *
 *            control_sarray - returns char ** - string
 *            These return pointers to particular elements in a control array.
 * COMMENT  : control_array routines
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define CONTROL_ARRAY_C
#include <stdlib.h>
#include "mms.h"

/**************************************************************************
 * control_array.c: 
 *
 * returns a pointer to a particular entry in a CONTROL struct
 *
 * index is 0-based, max size-1
 *
 **************************************************************************/

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_array
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *control_array (char *key, long ind) {
 
  CONTROL *control;

  if ((control = control_addr(key)) == NULL) {
    (void)fprintf(stderr, 
	    "ERROR - control_array - key '%s' not found.\n", key);
    exit(1);
  }

  if (ind >= control->size) {
    (void)fprintf(stderr, 
	    "ERROR - control_array - ind %ld too high for %s.\n", ind, key);
    (void)fprintf(stderr, 
	    "Max ind is %ld.\n", control->size-1);
    exit(1);
  }

	switch (control->type) {
		case M_DOUBLE:
			return (char *) ((double *)(control->start_ptr) + ind * sizeof(double));

		case M_FLOAT:
			return (char *) ((float *)(control->start_ptr) + ind * sizeof(float));

		case M_LONG:
			return (char *) ((long *)(control->start_ptr) + ind * sizeof(long));

		case M_STRING:
			printf ("control_array: key = %s ind = %ld val = %s\n", key, ind, *((char **)control->start_ptr + ind));
//			return (char *) (((char **)(control->start_ptr)) + (ind * sizeof(char *)));
			return *((char **)control->start_ptr + ind);
	}

	return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_larray
 | COMMENT		: returns a pointer to a long entry 
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long *control_larray (char *key, long ind) {
  return ((long *) control_array(key, ind));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_farray
 | COMMENT		: returns a pointer to a float entry in control array
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
float *control_farray (char *key, long ind) {
  return ((float *) control_array(key, ind));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_darray
 | COMMENT		: returns a pointer to a double entry in control array
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double *control_darray (char *key, long ind) {
  return ((double *) control_array(key, ind));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_sarray
 | COMMENT		: returns a pointer to a string entry in control array
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *control_sarray (char *key, long ind) {
  return control_array(key, ind);
}
