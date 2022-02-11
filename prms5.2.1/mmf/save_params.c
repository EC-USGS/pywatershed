/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : save_params
 * COMMENT  : saves the param data base to a file. File name is passed in.
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define SAVE_PARAMS_C
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mms.h"

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static void write_parameters (FILE *, int);
static void write_dimensions (FILE *);
static void write_header (FILE *, char *);

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : save_params
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int save_params (char *param_file_name) {
	FILE *param_file;
	//PARAM *param;
	//DIMEN *dim;
	//char *ptr;
	//long i,j;
	//double	*dvalptr;
	//float	*fvalptr;
	//long	*lvalptr;

	if ((param_file = fopen (param_file_name, "w")) == NULL) {
		(void)fprintf(stderr, "ERROR - save_params - creating file '%s'\n", param_file_name);
		return(1);
	}

	write_header (param_file, "Default Parameter File generated based on active modules and any specified Parameter File(s)\n");
	write_dimensions (param_file);
	write_parameters (param_file, TRUE);
	
	fclose(param_file);
	return(0);
}

int write_preprocess_params () {
	FILE *param_file;
	char param_file_name[512];
	char   **fname;
	/*char *extension, *ptr, *ptr1;*/
	char *ptr, *ptr1;

	fname =   control_svar ("param_file");
	strncpy (param_file_name, fname[0], 512);

// Isolate the file name from the path
	ptr1 = strrchr (param_file_name, '/');

// Find the last "." in the file name
	if (!ptr1) {
		ptr = NULL;
	} else {
		ptr = strrchr (ptr1, '.');
	}

	if (!ptr) {
		ptr = param_file_name + strlen(param_file_name);
	}
	strncpy (ptr, "_preprocess.params", 512-strlen(param_file_name));


	printf ("NOTICE: preprocessed parameters are being written to file: %s\n", param_file_name);

	if ((param_file = fopen (param_file_name, "w")) == NULL) {
		(void)fprintf(stderr, "ERROR - save_params - creating file '%s'\n", param_file_name);
		return(1);
	}

	write_parameters (param_file, FALSE);
	return(0);
}

static void write_header (FILE *param_file, char *desc) {
    (void)fprintf (param_file, "%s\n", desc);
	(void)fprintf (param_file, "PRMS version 4\n");
}

static void write_dimensions (FILE *param_file) {
	DIMEN *dim;
	long i,j;
	(void)fprintf(param_file, "** Dimensions **\n");

	for (i = 0; i < dim_db->count; i++) {

		dim = (DIMEN *)(dim_db->itm[i]);

		(void)fprintf(param_file, "####\n");
		(void)fprintf(param_file, "%s\n", dim->name);
		(void)fprintf(param_file, "%ld\n", dim->value);
		for (j = 0; j < dim->value; j++) {
			if (dim->names && dim->names[j])
				(void)fprintf (param_file, "%s\n", dim->names[j]);
			if (dim->notes && dim->notes[j])
				(void)fprintf (param_file, "@%s\n", dim->notes[j]);
		}
	}
}


static void write_parameters (FILE *param_file, int writeAllParams) {
	PARAM *param;
//	char *ptr;
	long i,j;
	double	*dvalptr;
	float	*fvalptr;
//	long	*lvalptr;
	int	*lvalptr;

    // 2016-01-13 PAN: cvalptr declaration removed
    // char *cvalptr;

/*
* Write out parameter values and description if any.
*/
	if (writeAllParams) {
		(void)fprintf(param_file, "** Parameters **\n");
	}

	for (i = 0; i < Mnparams; i++) {
		param = Mparambase[i];

		if (writeAllParams || param->preprocess ) {
            if (Mdebuglevel >= M_FULLDEBUG)
                (void)fprintf (stderr, "Writing string param '%s'\n", param->key);

			(void)fprintf(param_file, "####\n");
			(void)fprintf(param_file, "%s %ld", param->key, param->column_width);
			if (param->format)
				(void)fprintf(param_file, " %s\n", param->format);
			else
				(void)fprintf (param_file, "\n");
			(void)fprintf (param_file, "%ld\n", param->ndimen);
			for (j = 0; j < param->ndimen; j++)
				(void)fprintf(param_file, "%s\n", param->dimen[j]->name);

			(void)fprintf(param_file, "%ld\n", param->size);
			(void)fprintf(param_file, "%ld\n", param->type);

			switch (param->type) {
				case M_DOUBLE:
					if (writeAllParams) {
						dvalptr = (double *) param->value;
					} else {
						dvalptr = (double *) (param->references[0]);
					}

					for (j = 0; j < param->size; j++) {
						(void)fprintf(param_file, "%.20le\n", *dvalptr);
						dvalptr++;
						//if (param->value_desc[j]) {
						 // while ((ptr = strchr (param->value_desc[j], '\n'))) {
							//*ptr = '\0';
							//(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
							//param->value_desc[j] = ptr + 1;
						 // }
						 // if (param->value_desc[j] && strlen (param->value_desc[j]))
							//(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
						//}
					}
					break;

				case M_FLOAT:
					if (writeAllParams) {
						fvalptr = (float *) param->value;
					} else {
						fvalptr = (float *) (param->references[0]);
					}

					for (j = 0; j < param->size; j++) {
						(void)fprintf(param_file, "%.12e\n", *fvalptr);
						fvalptr++;
						//if (param->value_desc[j]) {
						//  while ((ptr = strchr (param->value_desc[j], '\n'))) {
						//	*ptr = '\0';
						//	(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
						//	param->value_desc[j] = ptr + 1;
						//  }
						//  if (param->value_desc[j] && strlen (param->value_desc[j]))
						//	(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
						//}
					}
					break;

				case M_LONG:
					if (writeAllParams) {
//						lvalptr = (long *) param->value;
						lvalptr = (int *) param->value;
					} else {
//						lvalptr = (long *) (param->references[0]);
						lvalptr = (int *) (param->references[0]);
					}

					for (j = 0; j < param->size; j++) {
//						(void)fprintf(param_file, "%ld\n", *lvalptr);
						(void)fprintf(param_file, "%d\n", *lvalptr);
						lvalptr++;
						//if (param->value_desc[j]) {
						//  while ((ptr = strchr (param->value_desc[j], '\n'))) {
						//	*ptr = '\0';
						//	(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
						//	param->value_desc[j] = ptr + 1;
						//  }
						//  if (param->value_desc[j] && strlen (param->value_desc[j]))
						//	(void)fprintf (param_file, "@%s\n", param->value_desc[j]);
						//}
					}
					break;

                // 2015-12-17 PAN: Added following block to write out
                //                 parameters of type M_STRING. This
                //                 code does not handle writing out a
                //                 single parameter value.
				case M_STRING:
                    if (writeAllParams) {
//						cvalptr = (char *)param->value;
//					}
//					else {
//						cvalptr = (char *)(param->references[0]);
//					}
                        for (j = 0; j < param->size; j++) {
                            if (*((char **) param->value + j) == NULL || *((char **) param->value + j)[0] == '\0') {
                                (void)fprintf(param_file, "\n");
                            } else {
                                (void)fprintf(param_file, "%s\n", *((char **) param->value + j));
                            }
//							(void)fprintf(param_file, "%s\n", *cvalptr);
//                          cvalptr++;
			}
		}
					break;
//                    					if (writeAllParams) {
//						lvalptr = (long *) param->value;
//						lvalptr = (int *) param->value;
//					} else {
//						lvalptr = (long *) (param->references[0]);
//						lvalptr = (int *) (param->references[0]);
//					}
	}
}
	}
}

