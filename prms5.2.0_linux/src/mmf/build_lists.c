/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : build_lists
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define BUILD_LISTS_C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**2************************* LOCAL MACROS ****************************/

#define INPUT  1
#define OUTPUT  2

/*--------------------------------------------------------------------*\
 | FUNCTION     : ALLOC_list
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
LIST *ALLOC_list (char *name, int type, int size) {

	LIST	*list;

	list = (LIST *)malloc (sizeof (LIST));

	if (name)
		list->name = strdup (name);
	else
		list->name = NULL;

	list->size = size;
	list->count = 0;
	list->type = type;
	list->out = FALSE;
	list->user_data = NULL;
	if (size)
		list->itm = (void **)malloc (size * sizeof (void *));
	else 
		list->itm = NULL;

	return (list);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : DELETE_list
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void DELETE_list (LIST *list) {

	int		i;

	if (list->name)
		free (list->name);

	for (i = 0; i < list->count; i++)
		if (list->itm[i])
			free (list->itm[i]);

	free (list->itm);
	free (list);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : RESIZE_list
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void RESIZE_list (LIST *list, int new_size) {
	list->size = new_size;
	list->itm = (void **)realloc (list->itm, new_size * sizeof (void *));
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : ADD_to_list
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void ADD_to_list (LIST *list, void *itm) {
	if (list->count >= list->size)
		RESIZE_list (list, list->size + 100);

	list->itm[list->count++] = itm;
}

