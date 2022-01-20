/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : check_vars
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define CHECK_VARS_C
#include <string.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : CHECK_stat_vars
 | COMMENT		: Makes sure that the selected statistic variables
 |                  are valid.
 | PARAMETERS   :
 | RETURN VALUE : 0 - no bad ones are found
 |              : 1 - at least one bad one found
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *CHECK_stat_vars (void) {

	static char err_message[256];

	int		i, status = 0;
	char	**names, **elements, buf[80], *ptr;

	names = (char **) control_var("statVar_names");
	elements = (char **) control_var("statVar_element");

	for (i = 0; i < *((long *)control_var ("nstatVars")); i++) {
		(void)strncpy (buf, names[i], 80);
		ptr = strchr (buf, '.');
		if (ptr) *ptr = '\0';

		if (CheckIndices (buf, elements[i], M_VARIABLE)) {
			(void)fprintf (stderr, "ERROR - CHECK_stat_vars: %s[%s] is not a valid stat variable.\n", names[i], elements[i]);
			(void)snprintf (err_message, 256, "Set stat variables: %s[%s] is not a valid stat variable.\n", names[i], elements[i]);
			status = 1;
		}
	}

	if (status)
		return (err_message);
	else
		return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : CHECK_disp_vars
 | COMMENT		: Makes sure that the selected display variables
 |                  are valid.
 | PARAMETERS   :
 | RETURN VALUE : error message
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char * CHECK_disp_vars (void) {
       static char	err_message[256];
	int		status = 0;
	int		i, j;
	char	buf[MAXVARLEN], buf0[MAXVARLEN], buf1[MAXVARLEN], buf2[MAXVARLEN];
	char	*dv_name, *dv_index, *ptr;

	for (i = 0; i < *(control_lvar ("ndispGraphs")); i++) {
		(void)snprintf (buf0, MAXVARLEN, "ndispVars%d", i);
		(void)snprintf (buf1, MAXVARLEN, "dispVar_names%d", i);
		(void)snprintf (buf2, MAXVARLEN, "dispVar_element%d", i);

		for (j = 0; j < *(control_lvar (buf0)); j++) {
			dv_name = *((char **)(control_sarray (buf1, j)));
			(void)strncpy (buf, dv_name, MAXVARLEN);
			if ( (ptr = strchr (buf, '.')) )
				*ptr = '\0';
			dv_index = *((char **)(control_sarray (buf2, j)));
			if ((CheckIndices (buf, dv_index, M_VARIABLE))) {
				(void)fprintf (stderr, "ERROR - CHECK_disp_vars: %s[%s] from graph %d is not a valid display variable.\n", buf, dv_index, i+1);
				(void)snprintf (err_message, 256, "Set display variables: %s[%s] from graph %d is not a valid display variable.\n", buf, dv_index, i+1);
				status = 1;
			}
		}
	}

	if (status)
		return (err_message);
	else
		return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : CHECK_ani_vars
 | COMMENT      : Makes sure that the selected ani variables
 |                  are valid.
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *CHECK_ani_vars (void) {
    int     i, status = 0;
    char    **names, buf[80], *ptr;
    PUBVAR  *vaddr;
    static char err_message[256];

    names = (char **) control_var ("aniOutVar_names");

    for (i = 0; i < *((long *)control_var ("naniOutVars")); i++) {
        (void)strncpy (buf, names[i], 80);
        ptr = strchr (buf, '.');
        if (ptr) *ptr = '\0';

        if (!(vaddr = var_addr (buf))) {
           fprintf (stderr, "ERROR - CHECK_ani_vars: %s is not a valid animation output variable.\n", names[i]);
           snprintf (err_message, 256, "Set animation variables: %s is not a valid stat variable.\n", names[i]);
           status = 1;
        }
    }

    if (status)
        return (err_message);
    else
        return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : CHECK_map_vars
 | COMMENT      : Makes sure that the selected map variables
 |                  are valid.
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *CHECK_map_vars (void) {
    int     i, status = 0;
    char    **names, buf[80], *ptr;
    PUBVAR  *vaddr;
    static char err_message[256];

    names = (char **) control_var ("mapOutVar_names");

    for (i = 0; i < *((long *)control_var ("nmapOutVars")); i++) {
        (void)strncpy (buf, names[i], 80);
        ptr = strchr (buf, '.');
        if (ptr) *ptr = '\0';

        if (!(vaddr = var_addr (buf))) {
           fprintf (stderr, "ERROR - CHECK_map_vars: %s is not a valid map output variable.\n", names[i]);
           snprintf (err_message, 256, "Set map variables: %s is not a valid stat variable.\n", names[i]);
           status = 1;
        }
    }

    if (status)
        return (err_message);
    else
        return (NULL);
}

