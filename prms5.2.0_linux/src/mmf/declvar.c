/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : declvar() to be called from C
 *            declvar_() to be called from Fortran
 *            Returns 0 if successful, 1 otherwise.
 * COMMENT  : initializes a module variable entry in the memory database
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define DECLVAR_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

#define LONG 1
#define FLOAT 2
#define DOUBLE 3

/*--------------------------------------------------------------------*\
 | FUNCTION     : declvar_
 | COMMENT		: called from Fortran, sorts out args and calls declvar()
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long declvar_ (char *mname, char *vname, char *vdimen, ftnint *maxsizeptr,
	char *vtype, char *hstr, char *ustr, char *value, ftnlen mnamelen,
	ftnlen vnamelen, ftnlen vdimenlen, ftnlen vtypelen, ftnlen hlen, ftnlen ulen) {

  char *module, *name, *dimen, *type, *help, *units;
  long maxsize, retval;

  /*
   * copy maxsize to local long int
   */

  maxsize = *maxsizeptr;

  /*
   * copy args to new strings, and terminate correctly
   */

  module = (char *) umalloc((unsigned int)(mnamelen + 1));
  strncpy(module, mname, (int)mnamelen);
  module[mnamelen] = '\0';

  name = (char *) umalloc((unsigned int)(vnamelen + 1));
  strncpy(name, vname, (int)vnamelen);
  name[vnamelen] = '\0';

  dimen = (char *) umalloc((unsigned int)(vdimenlen + 1));
  strncpy(dimen, vdimen, (int)vdimenlen);
  dimen[vdimenlen] = '\0';

  type = (char *) umalloc((unsigned int)(vtypelen + 1));
  strncpy(type, vtype, (int)vtypelen);
  type[vtypelen] = '\0';

  help = (char *) umalloc((unsigned int)(hlen + 1));
  strncpy(help, hstr, (int)hlen);
  help[hlen] = '\0';

  units = (char *) umalloc((unsigned int)(ulen + 1));
  strncpy(units, ustr, (int)ulen);
  units[ulen] = '\0';

  /*
   * call C version of declvar()
   */

  retval = declvar(module, name, dimen, maxsize, type, help, units, value);

  return(retval);

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : declvar()
 | COMMENT		: is called from C
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long declvar (char *module, char *name, char *dimen, long maxsize, char *type,
	char *help, char *units, char *value) {
  int var_type;

  char *vkey;
  char *token;
  char *tmpdimen;
  long i, size;

  PUBVAR **vars, *var;

  /*
   * realloc if too large
   */

  if(Mnvars >= max_vars -1) {
	max_vars += 100;
  	Mvarbase = (PUBVAR **)urealloc ((char *)Mvarbase,
		max_vars * sizeof(PUBVAR *));
  }

  /*
   * compute the key
   */

  vkey = strdup (name);

  if (var_addr(vkey) != NULL) {
	  if (print_mode) {
	      return(0);
	  } else {
              fprintf(stderr,
	      "ERROR - declvar - key '%s' already exists.\n", vkey);
              return(1); }
  }

  /*
   * convert fortran types to C equivalents
   */

  var_type = M_LONG;
  if (!strcmp(type, "real") || !strcmp(type, "float"))
    var_type = M_FLOAT;
  else if (!strcmp(type, "double precision") || !strcmp(type, "double"))
    var_type = M_DOUBLE;

  /*
   * check that type is possible
   */

	if((var_type != M_LONG) && (var_type != M_FLOAT) && (var_type != M_DOUBLE))
		{
    	(void)fprintf(stderr,
	    	"ERROR - declvar - type '%s' is illegal.\n", type);
    	(void)fprintf(stderr, "Key is '%s'.\n", vkey);
    	(void)fprintf(stderr, "Type is '%s'.\n", type);
    	return(1);
  		}

  /*  
   * get vars from Mvarbase, the global pointer
   */


  if (Mdebuglevel >= M_FULLDEBUG) {
    (void)fprintf(stderr, "Declaring variable '%s'\n", vkey);
  }

  /*
   * allocate space for a structure, and store pointer in vars
   */
  Mnvars += 1;

  vars = Mvarbase;

  var = (PUBVAR *) umalloc (sizeof(PUBVAR));
  vars[Mnvars-1] = var; /* copy address into vars array */

  /*
   * determine dimensions
   */

  tmpdimen = strdup (dimen);

  var->ndimen = 0;

  token = strtok (tmpdimen, ",");

  while (token != (char *) NULL) {
    var->ndimen++;
    token = strtok((char *) NULL, ",");
  }

  if (var->ndimen > MAX_NDIMEN) {

    (void)fprintf(stderr, "ERROR - declvar\n");
    (void)fprintf(stderr, "Attempt to use %ld dimensions - this is too many.\n",
	    var->ndimen);
    (void)fprintf(stderr, "Max number of dimensions allowed : %d.\n", MAX_NDIMEN);
    (void)fprintf(stderr, "Key is '%s'.\n", vkey);
    return(1);

  }

  var->dimen = (DIMEN **)umalloc (var->ndimen * sizeof(DIMEN *));

  (void)strncpy (tmpdimen, dimen, strlen(dimen)+1);

  i = 0;
  token = strtok(tmpdimen, ",");

  while (token != (char *) NULL) {
    if (!(var->dimen[i] = dim_addr (token))) {
      (void)fprintf(stderr, "ERROR - declvar\n");
      (void)fprintf(stderr, "Variable '%s'\n", vkey);
      (void)fprintf(stderr, "Dimension '%s' is not declared.\n", token);
      return(1);
    }
    token = strtok ((char *) NULL, ",");
    i++;
  }

  /*
   * get the size of the variable
   */
  
  size = 1;

  for (i = 0; i < var->ndimen; i++) {

    size *= var->dimen[i]->value;

  }

  var->size = size;

  if (size > maxsize) {

    (void)fprintf(stderr,
	    "ERROR - declvar - dimension exceeds space available.\n");
    (void)fprintf(stderr, "Key is '%s'.\n", vkey);
    (void)fprintf(stderr, "Size is %ld.\n", size);
    for (i = 0; i < var->ndimen; i++) {
      (void)fprintf (stderr, "Dimension '%s' is %ld.\n",
	      var->dimen[i]->name, var->dimen[i]->value);
    }
    (void)fprintf(stderr, "Space available is %ld.\n", maxsize);
    return(1);

  }

  /*
   * allocate space, and store variable properties
   */

  var->key = vkey;
  var->module = strdup (module);
  var->name = strdup (name);

  if(var_type == M_DOUBLE)
    var->type = M_DOUBLE;
	else 
		if (var_type == M_FLOAT)
    		var->type = M_FLOAT;
  			else 
				if (var_type == M_LONG) 
    				var->type = M_LONG;
  var->value = value;
  var->help = strdup (help);
  var->units = strdup (units);
  var->private = FALSE;

  sort_vars();

  return(0);
  
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : declpri_
 | COMMENT		: called from Fortran, sorts out args and calls declpri()
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long declpri_ (char *vname, ftnint *maxsizeptr,
	char *vtype, char *value,
	ftnlen vnamelen, ftnlen vtypelen) {

  char *name, *type;
  long maxsize, retval;

  /*
   * copy maxsize to local long int
   */

  maxsize = *maxsizeptr;

  /*
   * copy args to new strings, and terminate correctly
   */

  name = (char *) umalloc((unsigned int)(vnamelen + 1));
  strncpy(name, vname, (int)vnamelen);
  name[vnamelen] = '\0';

  type = (char *) umalloc((unsigned int)(vtypelen + 1));
  strncpy(type, vtype, (int)vtypelen);
  type[vtypelen] = '\0';

  /*
   * call C version of declpri()
   */

  retval = declpri(name, maxsize, type, value);

  /*
   * free up allocated strings
   */

  return(retval);

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : declpri()
 | COMMENT		: is called from C
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long declpri (char *name, long size, char *type, char *value) {
  int var_type;

  char *vkey;

  PUBVAR **vars, *var;

  /*
   * realloc if too large
   */

  if(Mnvars >= max_vars -1) {
	max_vars += 100;
  	Mvarbase = (PUBVAR **)urealloc ((char *)Mvarbase,
		max_vars * sizeof(PUBVAR *));
  }

  /*
   * compute the key
   */

  vkey = strdup (name);

  if (var_addr(vkey) != NULL) {
    (void)fprintf(stderr,
	    "ERROR - declvar - key '%s' already exists.\n", vkey);
    return(1);
  }

  /*
   * convert fortran types to C equivalents
   */

  var_type = M_LONG;
  if (!strcmp(type, "real") || !strcmp(type, "float"))
    var_type = M_FLOAT;
  else if (!strcmp(type, "double precision") || !strcmp(type, "double"))
    var_type = M_DOUBLE;

  /*
   * check that type is possible
   */

	if((var_type != M_LONG) && (var_type != M_FLOAT) && (var_type != M_DOUBLE))
		{
    	(void)fprintf(stderr,
	    	"ERROR - declvar - type '%s' is illegal.\n", type);
    	(void)fprintf(stderr, "Key is '%s'.\n", vkey);
    	(void)fprintf(stderr, "Type is '%s'.\n", type);
    	return(1);
  		}

  /*  
   * get vars from Mvarbase, the global pointer
   */


  if (Mdebuglevel >= M_FULLDEBUG) {
    (void)fprintf(stderr, "Declaring private variable '%s'\n", vkey);
  }

  /*
   * allocate space for a structure, and store pointer in vars
   */
  Mnvars += 1;

  vars = Mvarbase;

  var = (PUBVAR *) umalloc (sizeof(PUBVAR));
  vars[Mnvars-1] = var; /* copy address into vars array */

  /*
   * get the size of the variable
   */
  
  var->size = size;

  /*
   * allocate space, and store variable properties
   */

  var->key = vkey;
  var->module = NULL;
  var->name = strdup (name);

   if(var_type == M_DOUBLE) var->type = M_DOUBLE;
   else if (var_type == M_FLOAT) var->type = M_FLOAT;
   else if (var_type == M_LONG) var->type = M_LONG;

   var->value = value;
   var->help = NULL;
   var->units = NULL;
   var->private = TRUE;

   sort_vars();

   return(0);
}

