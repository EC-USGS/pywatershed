/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : control_addr
 * returns a pointer to a CONTROL struct which contains the given key
 * returns NULL if key not found
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define CONTROL_ADDR_C
#include <stdio.h>
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_addr
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
CONTROL *control_addr (char *key) { 
	CONTROL	*cp;
	int		i;

	for (i = 0; i < cont_db->count; i++) {
		cp = (CONTROL *)(cont_db->itm[i]);
		//printf ("control_addr: i = %d comparing %s to %s \n", i, key, cp->key);
		if (!strcmp (cp->key, key))
			return (cp);
	}

	return (NULL);
}

