/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : decl_control
 * COMMENT  : initializes a module variable entry in the memory database
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define DECL_CONTROL_C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**************************************************************************/

/*--------------------------------------------------------------------*\
 | FUNCTION     : add_control
 | COMMENT		: This allocates a control structure and adds it to the
 |                control DB.  It also allocates the space for the variables.
 | PARAMETERS   :
 | RETURN VALUE : None
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
CONTROL *add_control (char *key, long type, long size) {
   CONTROL *cp;

/*
**	check that key does not already exist
*/

   if (control_addr (key)) {
      (void)fprintf (stderr,
         "ERROR - add_control - key '%s' already exists.\n", key);
      exit(1);
   }
// printf ("adding control parameter - key: %s type: %ld size: %ld\n", key, type, size);

/*
**  allocate space for a structure, and store pointer in controls
**  allocate space, and store control variable properties
*/
   cp = (CONTROL *) umalloc (sizeof(CONTROL));
   ADD_to_list (cont_db, (void *)cp);

   cp->key = strdup (key);
   cp->size = size;
   cp->type = type;
   cp->set_in_file = 0;

   if (type == M_STRING) {
      cp->start_ptr = (char *)umalloc (sizeof (char *) * size);
   
   } else if (type == M_LONG) {
      cp->start_ptr = (char *)umalloc (sizeof (long) * size);

   } else if (type == M_FLOAT) {
	   cp->start_ptr = (char *)umalloc (sizeof (float) * size);

   } else if (type == M_DOUBLE) {
	   cp->start_ptr = (char *)umalloc (sizeof (double) * size);

   } else {
      (void)fprintf (stderr,
         "ERROR - add_control - key '%s' don't know what type code %ld is.\n", key, type);
      exit(1);
   }

   return cp;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : decl_control
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : None
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void decl_control (char *key, long type, long size, void *valstr) {
   CONTROL *cp;

/*
**	check that key does not already exist
*/

   if (control_addr (key)) {
      (void)fprintf (stderr,
         "ERROR - decl_control - key '%s' already exists.\n", key);
      exit(1);
   }

/*
**  allocate space for a structure, and store pointer in controls
**  allocate space, and store control variable properties
*/
   cp = (CONTROL *) umalloc (sizeof(CONTROL));
   ADD_to_list (cont_db, (void *)cp);

   cp->key = key;
   cp->size = size;
   cp->type = type;
   cp->start_ptr = (char *)valstr;
   cp->set_in_file = 0;

}

void decl_control_string (char *key, char *valstr) {
   char **cp;
   cp = (char **)umalloc (sizeof (char *) * 1);
   *cp = strdup (valstr);
   decl_control (strdup (key), M_STRING, 1, cp);
}

void decl_control_string_array (char *key, long size, char *valstr) {
   char **cp;
   int i;

   cp = (char **)umalloc (sizeof (char *) * size);
   for (i = 0; i < size; i++) {
      cp[i] = strdup (valstr);
   }

   decl_control (strdup (key), M_STRING, size, cp);
}

void decl_control_int_array (char *key, long size, long *valstr) {
   long *lp;
   int i;

   lp = (long *)umalloc (sizeof (long) * size);
   for (i = 0; i < size; i++) {
      lp[i] = valstr[i];
   }

   decl_control (strdup (key), M_LONG, size, lp);
}

void decl_control_float_array (char *key, long size, float *valstr) {
   float *fp;
   int i;

   fp = (float *)umalloc (sizeof (float) * size);
   for (i = 0; i < size; i++) {
      fp[i] = (float)(valstr[i]);
   }

   decl_control (strdup (key), M_FLOAT, size, fp);
}

void decl_control_double_array (char *key, long size, double *valstr) {
   double *fp;
   int i;

   fp = (double *)umalloc (sizeof (double) * size);
   for (i = 0; i < size; i++) {
      fp[i] = (double)(valstr[i]);
   }

   decl_control (strdup (key), M_DOUBLE, size, fp);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : decl_control_
 | COMMENT		: decl_control_() is called from Fortran, sorts out args
 |                 and calls decl_control()
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void decl_control_ (char *ckey, ftnint *ctype, ftnint *csize, void *value, ftnlen klen) {
	char *key;
	long type, size;

  /*
   * copy ctype and csize to local long int
   */
	type = *ctype;
	size = *csize;

  /*
   * copy args to new strings, and terminate correctly
   */
	key = (char *) umalloc((unsigned int)(klen + 1));
	strncpy(key, ckey, (int)klen);
	key[klen] = '\0';

	decl_control(key, type, size, value);
	return;
}
