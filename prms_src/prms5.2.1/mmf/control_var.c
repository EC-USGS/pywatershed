/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : control_var - generic, returns (char *) as a generic pointer
 *            control_lvar - returns long *
 *            control_fvar - returns float *
 *            control_dvar - returns double *
 *            control_svar - returns char ** - string
 *            returns pointers to various control array entries
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define CONTROL_VAR_C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_var
 | COMMENT		: returns a pointer to the start of the variable
 |			( or first element in * the array) in a CONTROL struct
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *control_var (char *key) {
 
  CONTROL *control;

  if ((control = control_addr(key)) == NULL) {
    (void)fprintf(stderr, 
	    "ERROR - control_var - key '%s' not found.\n", key);
    exit(1);
  }
  return (char *) control->start_ptr;

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_lvar
 | COMMENT		: returns a pointer to a long variable
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long *control_lvar (char *key) {
  return ((long *) control_var(key));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_fvar
 | COMMENT		: returns a pointer to a float variable
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
float *control_fvar (char *key) {
  return ((float *) control_var(key));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_dvar
 | COMMENT		: returns a pointer to a double variable
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double *control_dvar (char *key) {
  return ((double *) control_var(key));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_svar
 | COMMENT		: returns a pointer to a string variable
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char **control_svar (char *key) {
  return ((char **) control_var(key));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_string_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_string_ (char *retval, char *tag, ftnlen len, ftnlen tlen) {
	char *foo;

	foo = (char *) umalloc(tlen + 1);
	strncpy(foo, tag, tlen);
	foo[tlen] = '\0';

	memset (retval, ' ', len);
	strncpy (retval, *control_svar(foo), len);
	return 0;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_string_array_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_string_array_ (char *retval, char *tag, int *index, ftnlen len, ftnlen tlen) {
	char *foo;
    char **strings;
    int i;

	foo = (char *) umalloc(tlen + 1);
	strncpy(foo, tag, tlen);
	foo[tlen] = '\0';

    strings = (char **) control_var(foo);
    i = *index - 1;
	strncpy (retval, *(strings+i), len);
	return 0;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_integer_
 | COMMENT		: returns a long variable value
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_integer_ (int *retval, char *key, ftnlen len) {
	char *foo;
	long *longs, intVal;

	foo = (char *) umalloc(len + 1);
	strncpy(foo, key, len);
	foo[len] = '\0';

	longs = (long *) control_var(foo);
	intVal = *(longs);
	*retval = (int)intVal;
	return 0;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_real_
 | COMMENT		: returns a float variable value
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_real_(float* retval, char* key, ftnlen len) {
	char* foo;

	foo = (char*)umalloc(len + 1);
	strncpy(foo, key, len);
	foo[len] = '\0';

	*retval = *control_fvar(foo);
	return 0;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_integer_array_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_integer_array_ (int *retval, int *index, char *key, ftnlen tlen) {
	char *foo;
	long intVal;
    long *longs;
    int i;

	foo = (char *) umalloc(tlen + 1);
	strncpy(foo, key, tlen);
	foo[tlen] = '\0';

    longs = (long *) control_var(foo);
    i = *index - 1;
	intVal = *(longs+i);
	*retval = (int)intVal;
	return 0;
}

/*--------------------------------------------------------------------*\
| FUNCTION     : control_file_name_
| COMMENT		: called from fortran
| PARAMETERS   :
| RETURN VALUE :
| RESTRICTIONS :
\*--------------------------------------------------------------------*/
long control_file_name_(char *retval, ftnlen tlen) {
	char *foo;
	foo = (char *)umalloc(tlen + 1);
	strncpy(foo, MAltContFile, tlen);
	foo[tlen] = '\0';
	memset(retval, ' ', tlen);
	strncpy(retval, foo, tlen);
	return 0;
}
