/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : readvar() to be called from C
 *            readvar_() to be called from Fortran
 *            returns 0 if success, 1 if failure
 * COMMENT  : reads the values associated with a key from an input file,
 *            and stores it in the data base
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READVAR_C
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "mms.h"

/**2************************* LOCAL MACROS ****************************/
#define MISSING_VAR -999

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : readvar_
 | COMMENT		: called from Fortran, sorts out args and calls readvar()
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long readvar_ (char *mname, char *vname, ftnlen mnamelen, ftnlen vnamelen) {
	char module[80], name[80];
	long retval;

/*
* copy args to new strings, and terminate
*/
	strncpy (module, mname, mnamelen);
	*(module + mnamelen) = '\0';

	strncpy (name, vname, vnamelen);
	*(name + vnamelen) = '\0';

/*
* call C version of readvar()
*/
	retval = readvar (module, name);

	return (retval);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : readvar
 | COMMENT		: called from C
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long readvar (char *module, char *name) {

	PUBVAR *var;
	long i, found;
	char vkey[80];
	long *long_var;
	float *float_var;
	double *double_var;

/*
* compute the key
*/

/*
  vkey = (char *) umalloc(strlen(module) + strlen(name) + 2);
  (void)strncpy(vkey, module, 80);
  strncat(strncat(vkey, ".", 80), name, 80);
*/
	strncpy (vkey, name, 80);

/*
* get pointer to variable with key
*/
	if (!(var = var_addr (vkey))) {
		(void)fprintf(stderr, "ERROR - readvar - variable not found.\n");
		(void)fprintf(stderr, "Key:   %s.\n", vkey);
		return(1);
	}

/*
* check that this is the correct variable, and that the size is
* set to the expected read count
*/
	found = -1;
	for (i = 0; i < Mnreads; i++) {
		if (var == Mcheckbase[i]->var) {
			found = i;
			break;
		}
	}

/*
	if (!(var->size))
		return (0);
*/

	if (found == -1) {
		(void)fprintf(stderr, "\nERROR: Attempting to read variable %s, which is not in Data File\n", vkey);
		return (1);
	}

/*
* data is present in file
*/

  if(var->size != Mcheckbase[found]->count) {
    (void)fprintf(stderr, "ERROR - readvar\n");
    (void)fprintf(stderr, "Reading var '%s'\n", vkey);
    (void)fprintf(stderr, "Attempting to read %ld items\n", var->size);
    (void)fprintf(stderr, "Data file has %ld items for this variable.\n",
	    Mcheckbase[found]->count);
    return(1);

  }
    
  /*
   * copy the variable from the input line into the data base,
   * according to the type, if size > 0
   */

  if (var->size > 0) {
  
    switch (var->type) {
      
    case M_LONG:
      long_var = (long *) var->value;
      for (i = 0; i < var->size; i++) {
	long_var[i] = Mcheckbase[found]->Types.valuel[i];
      }
      break;
      
    case M_FLOAT:
      float_var = (float *) var->value;
      for (i = 0; i < var->size; i++) {
	float_var[i] = Mcheckbase[found]->Types.valuef[i];
      }
      break;
      
    case M_DOUBLE:
      double_var = (double *) var->value;
      for (i = 0; i < var->size; i++) {
	double_var[i] = Mcheckbase[found]->Types.valued[i];
      }
      break;
      
    }
    

  }  /* if(var->size > 0) */


  return(0);
}
