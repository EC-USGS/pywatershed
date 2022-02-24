/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : getdimname
 * COMMENT  : The following are two routines to obtain the "ith" index name 
 *            of a dimension variable from either Fortran or C modules.
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define GETDIMNAME_C
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : getdimname_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void getdimname_ (char *name, ftnint *i, char *idxname, ftnlen namelen, ftnlen idxlen) {
  /*
   * local copies
   */

  char * lname;

  lname = (char *)malloc(namelen+1);
  strncpy(lname, name, namelen);
  lname[namelen] = '\0';

  /*
   * call c version
   */
 getdimname(lname, (*i) - 1, idxname, idxlen);
  
  idxlen = strlen(idxname);

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : getdimdesc_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void getdimdesc_ (char *name, ftnint *i, char *desc, ftnlen namelen, ftnlen desclen) {
  /*
   * local copies
   */

  char * lname;

  lname = (char *)malloc(namelen+1);
  strncpy(lname, name, namelen);
  lname[namelen] = '\0';

/*
**	call c version
*/
	getdimdesc (lname, (*i) - 1, desc, namelen+1);
	desclen = strlen (desc);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : getdimnameint_
 | COMMENT		: called from fortran
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void getdimnameint_ (char *name, ftnint *i, ftnint *idx, ftnlen namelen) {
  /*
   * local copies
   */

  char * lname;
  char idxname[80];

  lname = (char *)malloc(namelen+1);
  strncpy(lname, name, namelen);
  lname[namelen] = '\0';

  /*
   * call c version
   */
 getdimname(lname, (*i) - 1, idxname, 80);
  
  *idx = atoi(idxname);

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : getdimname
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void getdimname (char *name, long i, char *idxname, int idxlen) {
  DIMEN *dim;

  dim = dim_addr(name);
  if (!dim) {
      (void)fprintf(stderr, "ERROR - getdimname, Can't find dimension named %s\n",name);
      return;
	}
  
  if (!dim->names) {
      (void)fprintf(stderr, "ERROR - getdimname. Dimension %s has no named indices\n",name);
      return;
    }
  (void)strncpy(idxname, dim->names[i], idxlen);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : getdimdesc
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void getdimdesc (char *name, long i, char *descname, int deslen) {
	DIMEN *dim;

	dim = dim_addr(name);
	if (!dim) {
		(void)fprintf (stderr,
			"ERROR - getdimname, Can't find dimension named %s\n", name);
		(void)strncpy (descname, "", deslen);
		return;
	}
  
	if (!dim->notes) {
		(void)strncpy (descname, "", deslen);
		return;
	}

	if (dim->notes[i])
		(void)strncpy (descname, dim->notes[i], deslen);
	else
		(void)strncpy (descname, "", deslen);

}
