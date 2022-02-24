/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : read_params
 * COMMENT  : reads the params data base from a file
 *            File name is passed in as an argument
 *
 * $Id$
 * Last modified 03/27/2018 RSR modified for multiple values per line separated by spaces
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_PARAMS_C
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <stdlib.h>
#include "mms.h"

static char *READ_param_head (PARAM **, int);
static char *READ_param_values (long, long, char *, char *, char[], char *, char *, int);
static char *rp (int, int);
static int checkForValidDimensions (PARAM *);
static int isDimensionIncompatable (char *, char *);
static void oneToAnySizedArray(PARAM *, char *);
static int getParamFileParamSize (PARAM *);
static char *getMapParamName(char *);
static void subbasinTo1DArray (PARAM *, PARAM *, char *);

static char *open_parameter_file();
static void close_parameter_file ();
static char *get_next_line ();
static char *error_string (char *);
static char *warning_string (char *);
static void bad_param_value_l (long, int, char *, long, long);
static void bad_param_value (double, int, char *, double, double);

static char* dimNames[] = {"nhru", "nsegment", "nrain", "ntemp", "nobs", "ngw", "nssr"};

//static char* mapParamNames[] = {"hru_subbasin", "segment_subbasin",
//	"rain_subbasin", "temp_subbasin", "obs_subbasin", "gw_subbasin",
//	"ssr_subbasin"
static char* mapParamNames[] = { "bad", "bad1",	"bad2", "bad3", "bad4", "bad5",	"bad6"
};

int nComments;
char **Comments;
static int lineNumber;
static FILE *param_file;
static char *line = NULL;
static char *key = NULL;
static char file_name[256];

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_params
 | COMMENT		: This is called from within a loop for each parameter file
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_params (char *param_file_name, int index, int mapping_flag) {
  	char *cptr;

// Get the static variables ready.
	strncpy (file_name, param_file_name, 256);
	if (key == NULL) {
		key = (char *) umalloc(max_data_ln_len * sizeof(char));
	}

	cptr = rp (index, mapping_flag);
	return (cptr);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_dims
 | COMMENT	: This is called once from MMF for the first parameter file
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_dims (char *param_file_name) {
	DIMEN *dim;
	int dim_size, i, j; 
	static char buf[256], *bp, *line_p;
	char *endptr;
	char *nch;
	int done;
    int max_comment;
  
/*
* get param name, open file
*/
	strncpy (file_name, param_file_name, 256);
    param_file = NULL;
    bp = open_parameter_file (param_file_name);
	if (bp != NULL) {
		return bp;
	}
/*
* read in run info string
*/
	if (!(line_p = get_next_line ())) {
		close_parameter_file();
		return error_string ("problems reading info");
	}
	Mparaminfo = strdup (line_p);

/*
**	See if version number is set
*/
	if (!(line_p = get_next_line ())) {
		close_parameter_file();
		return error_string ("problems reading version number");
	}

	if (!(line_p = get_next_line ())) {
		close_parameter_file();
		return error_string ("problems reading dimension label");
	}

/*
 *  Read in comments -- everything between version line and
 *  "** Dimensions **" line is a comment
 */
    max_comment = 1000;
	Comments = (char **)malloc (max_comment * sizeof (char *));
	nComments = 0;

	while (strncmp (line, "** Dimensions **", 16)) {
		if (!(line_p = get_next_line ())) {
			close_parameter_file();
			return error_string ("problems skipping comments");
		}

		if (strncmp (line, "** Dimensions **", 16)) {
            if (nComments < max_comment) {
//			     printf ("Comment line = %s\n", line);
			    Comments[nComments++] = strdup (line);
            } else {
                printf("read_dims: more than %d dimension comments.\n", max_comment);
            }
		}
	}

/*
**	Check dimension label
*/
	if (strncmp (line, "** Dimensions **", 16)) {
		close_parameter_file();
		return error_string ("** Dimensions ** label not found");
	}
  
	if (!(line_p = get_next_line ())) {
		close_parameter_file();
		return error_string ("unexpected end of file");
	}

/*
* read in dimensions
*/
	while (strncmp (line, "** Parameters **", 16)) {
		if (strncmp (line, "####", 4)) {
			close_parameter_file();
			return error_string ("expecting '####'");
		}

/*
**	Read dimension name from parameter file.
*/
		if (!(line_p = get_next_line ())) {
			close_parameter_file();
			return error_string ("trying to read dimension name");
		}

		line[strlen(line)-1] = '\0';
		dim = dim_addr (line);
		if (dim) {
/*
**	Read dimension size from parameter file.
*/
			if (!(line_p = get_next_line ())) {
				close_parameter_file();
				return error_string ("can't read dimension size");
			}

			errno = 0;
			dim_size = strtol(line, &endptr, 10);
			if (errno != 0) {
				close_parameter_file();
				return error_string ("size problem");
			}

/*
**	If necessary, reset dimension to value read from file.
*/
			if (dim->value != dim_size) {
				//reset_dim (dim, dim_size);
				dim->value = dim_size;
			}

/*
* check if there are index names below
*/
			if ((line_p = get_next_line ())) {
				if (strncmp (line, "** Parameters **", 16)) {
					if (dim->names) {
				//        free (dim->names);
						dim->names = NULL;
					}

					if (dim->notes) {
					//        free (dim->notes);
						dim->notes = NULL;
					}

					if (strncmp (line, "####", 4)) {
						dim->names = (char **)calloc (dim_size, sizeof (char *));
						dim->notes = (char **)calloc (dim_size, sizeof (char *));

						done = FALSE;
						i = 0;
						while (!done) {
							if (!strncmp (line, "####", 4)) {
								for (j = i; j < dim_size; j++) {
									dim->names[j] = NULL;
									dim->notes[j] = NULL;
								}
								done = TRUE;

							} else if (line[0] == '@') {
/* DANGER something wrong when reading dimension index.
** when i comes in as 0, i gets decremented and notes[i] is
** indexed to -1
*/
								i--;
								nch = (char *)strchr (line, '\n');
								if (nch) {
									*nch = '\0';
								}
								dim->notes[i] = strdup (&(line[1]));
								line_p = get_next_line ();
								i++;

							} else {
								nch = (char *)strchr (line, '\n');
								if (nch) {
									*nch = '\0';
								}
								dim->names[i] = strdup (line);
								line_p = get_next_line ();
								i++;
							}

							if ((i > dim_size) || ((i == dim_size) && (line[0] != '@'))) {
								done = TRUE;
							}
						}
					} else {
						dim->names = NULL;
						dim->files = NULL;
						dim->notes = NULL;
					}
				}
			} else {
				close_parameter_file();
				return NULL;
			}
		} else {
			key = (char *)strtok(line_p, " ");
			key[strlen(key)] = '\0';
			snprintf(buf, 256, "dimension '%s' is not required", key);
			fprintf (stderr,"\n%s\n", warning_string (buf));
			line_p = get_next_line ();
			line_p = get_next_line ();
		}
	}

	close_parameter_file();
	return NULL;
} // read_dims

/*--------------------------------------------------------------------*\
 | FUNCTION     : rp
 | COMMENT	:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *rp (int index, int map_flag) {
	PARAM *param;
	static char buf[256], *buf_ptr;
	char *line_p;
	char *pf_value, *mapParamName;
	int i, j, k;
	PARAM *mapping_param;

/*
* get param name, open file
*/
    buf_ptr = open_parameter_file ();
	if (buf_ptr != NULL) {
		return (buf_ptr);
	}

	get_next_line ();
	if (index == 0) {  // if index equals zero, than this parameter file has dimension stuff and we need to skip over it.
		while (strncmp (line, "** Parameters **", 16)) {
			if (!(line_p = get_next_line ())) {
				close_parameter_file();
			}
		}
		get_next_line ();
	}

/*
**	Read in parameters.
*/
	while (!feof (param_file)) {
		buf_ptr = READ_param_head (&param, map_flag);
		if (buf_ptr) {
			if (buf_ptr == (char *)-1) {
				close_parameter_file();
			  return NULL;
			} else {
				close_parameter_file();
			   return buf_ptr;
			}
		}

		if (param != NULL) {

			// If param->size (array size as defined by module) is the
			// same as param->pf_size (array size as defined by parameter
			// file, then read the values for this parameter from the
			// parameter file directly into param->value because the 
			// dimensions match.  If the sizes do not match, read the
			// values into a temporary array that is used to remap the
			// values into the correct size and shape for param->value.
			if (param->pf_size == param->size) {
				buf_ptr = READ_param_values (param->size, param->type, param->name, param->value, line, param->min_string, param->max_string, param->bound_status);

			} else {

//  Make sure that this resizing code is in sync with the ParamToolExpandor code in the oui4 code base.

				if (param->type == M_DOUBLE) {
					pf_value = (char *)umalloc (param->pf_size * sizeof (double));
				} else if (param->type == M_FLOAT) {
					pf_value = (char *)umalloc (param->pf_size * sizeof (float));
				} else if (param->type == M_LONG) {
					pf_value = (char *)umalloc (param->pf_size * sizeof (int));
				} else if (param->type == M_STRING) {
					pf_value = (char *)umalloc (param->pf_size * sizeof (char *));
				} else {
					pf_value = NULL;
				}

				buf_ptr = READ_param_values (param->pf_size, param->type, param->name, pf_value, line, param->min_string, param->max_string, param->bound_status);

				// The values read from the parameter file need to be resized to fit into the size
				// of the module array for this parameter.

				// It's easy when the size is 1.  Tested this works for floats
				if (param->pf_size == 1) {
					oneToAnySizedArray(param, pf_value);

				} else if (param->pf_ndimen == 1 && !strncmp(param->pf_dimNames[0], "nsub", 4) &&
						(!strncmp(param->dimen[0]->name, "nhru", 4) || !strncmp(param->dimen[0]->name, "nsegment", 8) ||
						!strncmp(param->dimen[0]->name, "nrain", 5) || !strncmp(param->dimen[0]->name, "ntemp", 5) ||
						!strncmp(param->dimen[0]->name, "nobs", 4) || !strncmp(param->dimen[0]->name, "ngw", 3) ||
						!strncmp(param->dimen[0]->name, "nssr", 4))) {  // subbasin to one mappable dimension

					mapParamName = getMapParamName(param->dimen[0]->name);

					mapping_param = param_addr (mapParamName);

					if (!mapping_param || !(mapping_param->read_in)) {
						snprintf (buf, 256, "\nERROR: mapping parameter %s must be set in parameter file before parameter %s\n",
							mapParamName, param->name);
						return (buf);
					}

					subbasinTo1DArray (param, mapping_param, pf_value);

				} else if (param->pf_ndimen == 1 && param->ndimen == 2) { // 1D in parameter file to 2D in module

   					// convert "nmonths" to "nhru,nmonths"
					if (!strncmp(param->pf_dimNames[0], "nmonths", 7)
                                && !strncmp(param->dimen[0]->name, "nhru", 4) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[j + (i*param->dimen[0]->value)] = ((double *)pf_value)[i];
								}
							}

						} else if (param->type == M_FLOAT) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[j + (i*param->dimen[0]->value)] = ((float *)pf_value)[i];
								}
							}

						} else if (param->type == M_LONG) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[j + (i*param->dimen[0]->value)] = ((int *)pf_value)[i];
								}
							}

						} else if (param->type == M_STRING) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (i*param->dimen[0]->value)) = strdup (pf_value + i);
								}
							}
						}

            // convert "nhru" to "nhru,nmonths"
					} else if (!strncmp(param->pf_dimNames[0], "nhru", 4)
                                && !strncmp(param->dimen[0]->name, "nhru", 4) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[k++] = ((double *)pf_value)[j];
								}
							}

						} else if (param->type == M_FLOAT) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[k++] = ((float *)pf_value)[j];
								}
							}

						} else if (param->type == M_LONG) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[k++] = ((int *)pf_value)[j];
								}
							}

						} else if (param->type == M_STRING) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (k++)) = strdup (pf_value + j);
								}
							}
						}

            // convert "nrain" to "nrain,nmonths"
					} else if (!strncmp(param->pf_dimNames[0], "nrain", 5)
                                && !strncmp(param->dimen[0]->name, "nrain", 5) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[k++] = ((double *)pf_value)[j];
								}
							}

						} else if (param->type == M_FLOAT) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[k++] = ((float *)pf_value)[j];
								}
							}

						} else if (param->type == M_LONG) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[k++] = ((int *)pf_value)[j];
								}
							}

						} else if (param->type == M_STRING) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (k++)) = strdup (pf_value + j);
								}
							}
						}

            // convert "ntemp" to "ntemp,nmonths"
					} else if (!strncmp(param->pf_dimNames[0], "ntemp", 5)
                                && !strncmp(param->dimen[0]->name, "ntemp", 5) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[k++] = ((double *)pf_value)[j];
								}
							}

						} else if (param->type == M_FLOAT) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[k++] = ((float *)pf_value)[j];
								}
							}

						} else if (param->type == M_LONG) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[k++] = ((int *)pf_value)[j];
								}
							}

						} else if (param->type == M_STRING) {
							k = 0;
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (k++)) = strdup (pf_value + j);
								}
							}
						}

            // convert "nmonths" to "nrain,nmonths"
					} else if (!strncmp(param->pf_dimNames[0], "nmonths", 7)
                                && !strncmp(param->dimen[0]->name, "nrain", 5) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[j + (i*param->dimen[0]->value)] = ((double *)pf_value)[i];
								}
							}

						} else if (param->type == M_FLOAT) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[j + (i*param->dimen[0]->value)] = ((float *)pf_value)[i];
								}
							}

						} else if (param->type == M_LONG) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[j + (i*param->dimen[0]->value)] = ((int *)pf_value)[i];
								}
							}

						} else if (param->type == M_STRING) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (i*param->dimen[0]->value)) = strdup (pf_value + i);
								}
							}
						}

            // convert "nmonths" to "ntemp,nmonths"
					} else if (!strncmp(param->pf_dimNames[0], "nmonths", 7)
                                && !strncmp(param->dimen[0]->name, "ntemp", 5) 
                                && !strncmp(param->dimen[1]->name, "nmonths", 7)) {

						if (param->type == M_DOUBLE) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((double *)(param->value))[j + (i*param->dimen[0]->value)] = ((double *)pf_value)[i];
								}
							}

						} else if (param->type == M_FLOAT) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((float *)(param->value))[j + (i*param->dimen[0]->value)] = ((float *)pf_value)[i];
								}
							}

						} else if (param->type == M_LONG) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									((int *)(param->value))[j + (i*param->dimen[0]->value)] = ((int *)pf_value)[i];
								}
							}

						} else if (param->type == M_STRING) {
							for (i = 0; i < param->dimen[1]->value; i++) {
								for (j = 0; j < param->dimen[0]->value; j++) {
									*((char **)param->value + (i*param->dimen[0]->value)) = strdup (pf_value + i);
								}
							}
						}
                    }  // end of 1D to 2D conversion code
				}
			}

			if (buf_ptr) {
				close_parameter_file();
				return (buf_ptr);
			}

			// This function copies the parameter values from the param structure
			// to the arrays in the modules.
			updateparam (param->name);
		}
	}

	close_parameter_file();
	return (NULL);
} // rp

/*--------------------------------------------------------------------*\
 | FUNCTION     : READ_param_head
 | COMMENT		: Read the preliminary stuff for the parameter.  This is
 |                 the stuff between the ####s and where the data actually
 |                 starts.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *READ_param_head (PARAM **param_ptr, int map_flag) {
	char dimen[MAXDATALNLEN];
	static char buf[256];
	char *temp, *npos, *tempfmt;
	int tempwidth, i, param_size, type;
	static long silent_flag;
	silent_flag = *control_lvar("print_debug");
	int badFlag;
	char *line_p;
  
/*
* space fwd to #### header
*/
	while (strncmp (line, "####", 4)) {
		if (!(line_p = get_next_line ())) {
			close_parameter_file();
			return (char *)-1;
		}
	}
  
/*
* get key, column width and format
*/
	if (!(line_p = get_next_line ())) {
		return error_string("Early end of file");
	}

/*
**	get the key
*/
	temp = (char *)strtok(line," ");
	npos = strchr(temp,'\n');
	if (npos) *npos = '\0';

	strncpy(key, temp, max_data_ln_len);
	key[strlen(temp)] = '\0';
	
/*
**	get the column width
*/
	temp = (char *)strtok(NULL," ");
	if (temp) {
		tempwidth = atoi(temp);
	} else {
		tempwidth = 0;
	}

/*
**	get the format
*/
	tempfmt = (char *)strtok(NULL," ");

/*
** markstro -- this check is added so that if there is just a space
**             after the width the parameter will not have a blank
**             format.
*/
	if (tempfmt && (strlen (tempfmt) < 2)) {
		tempfmt = NULL;
	}

/*
**  param is allocated by calls from the modules to declparam.
*/
	if (!map_flag) {
		*param_ptr = param_addr(key);

	} else {
		if (!strncmp (key, mapParamNames[0], strlen(key)) ||
			!strncmp (key, mapParamNames[1], strlen(key)) ||
			!strncmp (key, mapParamNames[2], strlen(key)) ||
			!strncmp (key, mapParamNames[3], strlen(key)) ||
			!strncmp (key, mapParamNames[4], strlen(key)) ||
			!strncmp (key, mapParamNames[5], strlen(key)) ||
			!strncmp (key, mapParamNames[6], strlen(key))) {
			*param_ptr = param_addr (key);

			if (*param_ptr == NULL) {  // Didn't find this mapping parameter in the parameter DB so declare one
				declparam ("read_params", key, NULL, "integer", NULL, NULL, NULL, NULL, NULL, NULL);
				*param_ptr = param_addr(key);
			}

		} else {
			*param_ptr = NULL;
		}
	}

	if (*param_ptr) {
	  /*
	  **  Set the read_in flag to true
	  */
		(*param_ptr)->read_in = 1;
/*
* save format and column width
*/
		(*param_ptr)->column_width = tempwidth;
		if (tempfmt) {
			tempfmt[strlen(tempfmt)-1] = '\0';
			if(!(*param_ptr)->format) {
				(*param_ptr)->format = (char *)(malloc(strlen(tempfmt)+1));
			} else {
				(*param_ptr)->format = (char *)(realloc((*param_ptr)->format, strlen(tempfmt) + 1));
			}   
			(void)strncpy((*param_ptr)->format, tempfmt, strlen(tempfmt)+1);
		} else {
			(*param_ptr)->format = NULL;
		}
/*
* get number of dimensions
*/
		if (!(line_p = get_next_line ())) {
			return error_string("number of dimensions");
		}

		if (isdigit(*line)) {
			(*param_ptr)->pf_ndimen = atol(line);

			if((*param_ptr)->pf_ndimen == 0) {
				return error_string("number of dimensions is 0");
			}
/*
* get number of dimensions if file format supports 2D arrays. Otherwise
* get dimension name.
*/
			(*param_ptr)->pf_dimNames = (char **)malloc ((*param_ptr)->pf_ndimen * sizeof (char *));

			for (i = 0; i < (*param_ptr)->pf_ndimen; i++) {
				if (!(line_p = get_next_line ())) {
					return error_string("number of dimensions is wrong");
				}
                strncpy (dimen, line, MAXDATALNLEN);
				dimen[strlen(dimen) - 1] = '\0';
				(*param_ptr)->pf_dimNames[i] = strdup(dimen);
			}

			if (map_flag) { // Need to set some values in the param structure for mapping parameter
				(*param_ptr)->ndimen = 1;
				(*param_ptr)->dimen = (DIMEN **)umalloc ((*param_ptr)->ndimen * sizeof (DIMEN *));
				(*param_ptr)->dimen[0] = dim_addr((*param_ptr)->pf_dimNames[0]);
			}

			badFlag = checkForValidDimensions (*param_ptr);  // 0 = good;  1 = bad

			if (badFlag) {
				return error_string("dimensions are incompatable with declaration in module");
			}

			(*param_ptr)->pf_size = getParamFileParamSize(*param_ptr);

			if (map_flag) { // Need to set some values in the param structure for mapping parameter
				(*param_ptr)->size = (*param_ptr)->pf_size;
				(*param_ptr)->value = (char *)umalloc ((*param_ptr)->size * sizeof (int)); // Mapping parameters are always integers
			}
/*
* get param size
*/
			if (!(line_p = get_next_line ())) {
				return error_string("incorrect parameter size");
			}

			if((param_size = atol(line)) == 0) {
				return error_string("incorrect parameter size");
			}

			if (param_size != (*param_ptr)->pf_size) {
				return error_string("incorrect parameter size");
			}

		} else {  //  number of dimensions not a digit
			(*param_ptr)->ndimen = 1;
			strncpy(dimen, line, strlen(line));
			dimen[strlen(line)-1] = '\0';

			if (strcmp(dimen, (*param_ptr)->dimen[0]->name)) {
				return error_string("incorrect dimension specified");
			}

			(*param_ptr)->size = getdim(dimen);
			param_size = (*param_ptr)->size;
		}
/*
* get type
*/
		if (!(line_p = get_next_line ())) {
			return error_string("incorrect data type");
		}

		if ((type = atol(line)) == 0) {
			return error_string("incorrect data type");
		}

		if (type != (*param_ptr)->type) {
			return error_string("incorrect data type");
		}
  
	} else {
		if (!map_flag) {
			if (silent_flag > -2) {
				snprintf (buf, 256, "parameter '%s' is not required", key);
				fprintf (stderr,"\n%s\n", warning_string (buf));
			}
		}
	}

	return (NULL);
} // READ_param_head

/*--------------------------------------------------------------------*\
 | FUNCTION     : READ_param_values
 | COMMENT		: Read the values and comments for the parameter.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *READ_param_values (long size, long type, char *name,
				char *value, char *line, char *min_string, char *max_string,
				int bound_status) {
	int i, j;
    int  done;
	//int	desc_count = 0;
	int repeat_count;
	char delims[] = " ";
	char *result = NULL;
	char *comp_ptr = NULL;
	static char *crap = NULL;
	static char *crap2 = NULL;
	//float foo;
	double d, min_d = 0.0, max_d = 1.0;
	char *endp;
	long l, min_l = 0, max_l = 1;
	char *line_p;
	
	if (crap == NULL) {
		crap = (char *) umalloc(max_data_ln_len * sizeof(char));
	}

	if (crap2 == NULL) {
		crap2 = (char *) umalloc(max_data_ln_len * sizeof(char));
	}

/*
** Decode the min and max string based on data type
*/
	switch (type) {
		case M_STRING:
			break;

		case M_DOUBLE:
		case M_FLOAT:
			min_d = strtod(min_string, &endp);
			max_d = strtod(max_string, &endp);
			break;

		case M_LONG:
			min_l = strtol(min_string, &endp, 0);
			max_l = strtol(max_string, &endp, 0);
			break;
	} // switch

/*
**  Space for the values and value_desc are allocated in declparam
*/
	done = FALSE;
	i = 0;
	while (!done) {
		if (!(line_p = get_next_line ())) {
			done = TRUE;

		} else if (!strncmp (line, "####", 4)) {
			done = TRUE;

		} else {
			//desc_count = 0;
			result = NULL;
			strncpy (crap, line, max_data_ln_len);

			result = strtok (crap, delims);
			while (result != NULL && !done) {
				strncpy (crap2, result, max_data_ln_len);
				comp_ptr = strchr (crap2, '*');
				if (comp_ptr == NULL){
					repeat_count = 1;
					comp_ptr = crap2;
				} else {
					*comp_ptr = '\0';
					repeat_count = atol(crap2);
					comp_ptr++;
					//foo = (float) atof(comp_ptr);
				}

				for (j = 0; j < repeat_count && !done; j++) {
					if (i < size) {
						switch (type) {
							case M_STRING:
                                comp_ptr[strlen(comp_ptr)-1] = '\0';
								*((char **)value + i) = strdup (comp_ptr);
                                i++;
								break;

							case M_DOUBLE:
								d = strtod(comp_ptr, &endp);

								if (d < min_d || d > max_d) {
									bad_param_value (d, i, name, min_d, max_d);
								}

								if (comp_ptr != endp && (*endp == '\n' || *endp == '\0')) {
									((double *)value)[i++] = d;
								} else {
									return error_string("parameter format error");
								}
								break;

							case M_FLOAT:
								d = strtod(comp_ptr, &endp);

								if (d < min_d || d > max_d) {
									bad_param_value (d, i, name, min_d, max_d);
								}

								if (comp_ptr != endp && (*endp == '\n' || *endp == '\0')) {
									((float *)value)[i++] = (float)d;
								} else {
									return error_string("parameter format error");
								}
								break;

							case M_LONG:
								l = strtol(comp_ptr, &endp, 0);

// Does not check parameter ranges if bounded
								if ((l < min_l || l > max_l) &&
										bound_status == M_UNBOUNDED) {
									bad_param_value_l (l, i, name, min_l, max_l);
								}

								if (comp_ptr != endp && (*endp == '\n' || *endp == '\0')) {
									((int *)value)[i++] = (int)l;
								} else {
									return error_string("parameter format error");
								}
								break;
						} // switch
				 
					} else { // if (i < size)
						done = TRUE;
						i++;
					} // if (i < size)
				}
				result = strtok(NULL, delims);
			} // while
		}
	}

	if (i < size) {
		return error_string("too few values read for paramter");
	} else if (i > size && !done) {
		return error_string("too many values read for paramter");
	}
	return NULL;
} // READ_param_values

/*--------------------------------------------------------------------*\
 | FUNCTION     : checkForValidDimensions
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 0 = good;  1 = bad
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static int checkForValidDimensions (PARAM *param_ptr) {
	int i, badFlag;

//	printf ("checkForValidDimensions name = %s\n", param_ptr->name);
//	printf ("   pf_ndimen = %d; module_ndimen = %d\n", (int)(param_ptr->pf_ndimen), (int)(param_ptr->ndimen));

	if (param_ptr->pf_ndimen > param_ptr->ndimen ) { // more dimensions in the parameter file is always invalid
		return 1;

	} else if (param_ptr->pf_ndimen == param_ptr->ndimen ) {
		badFlag = 1;
		for (i = 0; i < param_ptr->pf_ndimen; i++) {  // check each dimension for compatiblilty
//printf ("   1 comparing %s to %s\n", param_ptr->pf_dimNames[i], param_ptr->dimen[i]->name);
			badFlag = isDimensionIncompatable (param_ptr->pf_dimNames[i], param_ptr->dimen[i]->name); // 0 = good;  1 = bad
		}
		if (badFlag == 1) {
			return 1;
		}

	} else { // less dimensions in the parameter file than declared in the module.
		badFlag = 1;
//printf ("   2 parameter file has %d dimensions\n", param_ptr->pf_ndimen);
		for (i = 0; i < param_ptr->ndimen; i++) {  // check each dimension for compatiblilty; only need to find one that is compatable
//printf ("   2 comparing %s to %s\n", param_ptr->pf_dimNames[0], param_ptr->dimen[i]->name);
			if (badFlag == 1) {
				badFlag = isDimensionIncompatable (param_ptr->pf_dimNames[0], param_ptr->dimen[i]->name); // 0 = good;  1 = bad
			}
		}
		if (badFlag == 1) {
			return 1;
		}
	}

	//param_ptr->ndimen
	return 0;
} // checkForValidDimensions

/*--------------------------------------------------------------------*\
 | FUNCTION     : isDimensionIncompatable
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 0 = good;  1 = bad
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static int isDimensionIncompatable (char *pfDimName, char *modDimName) {
//	char *dimNames[] ={"one",
//		"ncascade",
//		"ncascdgw",
//		"nsegment",
//		"npoigages",
//		"nsub",
//		"nhrucell",
//		"ngw",
//		"nhru",
//		"nssr",
//		"nsfres",
//		"nlake",
//		"nrain",
//		"nsol",
//		"ntemp",
//		"nratetbl",
//		"nwateruse",
//		"ndepl",
//		"ndeplval",
//		"ndays",
//		"nmonths",
//		"nlapse",
//		"nobs",
//		"nsnow",
//		"nform",
//		"nevap",
//		"nsfelev",
//		"nlakeelev",
//		"nwind",
//		"nhumid",
//		"ngate",
//		"nstage",
//		"ngate2",
//		"nstage2",
//		"ngate3",
//		"nstage3",
//		"ngate4",
//		"nstage4",
//		"mxnsos",
//	};

	if (!strncmp (pfDimName, modDimName, 10)) {  // a dimension is compatable with itself
		return 0; 
	}

	if (!strncmp (pfDimName, "one", 3)) {  // "one" in the parameter file is compatable with everything
		return 0; 
	}

	// Subbasin (nsub) can be mapped to these dimensions with mapping parameter
	// "nhru" "hru_subbasin";
    // "nsegment" "segment_subbasin";
    // "nrain" "rain_subbasin";
    // "ntemp" "temp_subbasin";
    // "nobs"  "obs_subbasin";
    // "ngw" "gw_subbasin";
    // "nssr" "ssr_subbasin";
	if (!strncmp (pfDimName, "nsub", 4)) {
		if (!strncmp (modDimName, "nhru", 4)) {
			return 0;

		} else if (!strncmp (modDimName, "nsegment", 8)) {
			return 0;

		} else if (!strncmp (modDimName, "nrain", 5)) {
			return 0;

		} else if (!strncmp (modDimName, "ntemp", 5)) {
			return 0;

		} else if (!strncmp (modDimName, "nobs", 4)) {
			return 0;

		} else if (!strncmp (modDimName, "ngw", 3)) {
			return 0;

		} else if (!strncmp (modDimName, "nssr", 4)) {
			return 0;
		}
	}
	return 1;
} // isDimensionIncompatable

/*--------------------------------------------------------------------*\
 | FUNCTION     : getParamFileParamSize
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static int getParamFileParamSize (PARAM *param) {
	int i, size;

	size = 1;
	for (i = 0; i < param->pf_ndimen; i++) {  // check each dimension for size
		size = size * getdim(param->pf_dimNames[i]);
	}
	return size;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : oneToAnySizedArray
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void oneToAnySizedArray(PARAM *param, char *pf_value) {
	int i;

	if (param->type == M_DOUBLE) {
		for (i = 0; i < param->size; i++) {
			((double *)(param->value))[i] = *((double *)pf_value);
		}
	} else if (param->type == M_FLOAT) {
		for (i = 0; i < param->size; i++) {
			((float *)(param->value))[i] = *((float *)pf_value);
		}
	} else if (param->type == M_LONG) {
		for (i = 0; i < param->size; i++) {
			((int *)(param->value))[i] = *((int *)pf_value);
		}
	} else if (param->type == M_STRING) {
		for (i = 0; i < param->size; i++) {
			*((char **)param->value + i) = strdup (pf_value);
		}
	}
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : getMapParamName
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *getMapParamName(char *name) {
	char *mapParamName;
	int i;

	mapParamName = NULL;
	for (i = 0; i < (sizeof (dimNames) / sizeof (dimNames[0])); i++) {
		if (!strncmp (name, dimNames[i], strlen(dimNames[i]))) {
			mapParamName = mapParamNames[i];
		} 
	}

    return mapParamName;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : subbasinTo1DArray
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void subbasinTo1DArray (PARAM *param, PARAM *mapping_param, char *pf_value) {
	int i, map;

	if (param->type == M_DOUBLE) {
		for (i = 0; i < param->size; i++) {
			map = ((int *)(mapping_param->value))[i];
			((double *)(param->value))[i] = ((double *)pf_value)[map - 1];
		}

	} else if (param->type == M_FLOAT) {
		for (i = 0; i < param->size; i++) {
			map = ((int *)(mapping_param->value))[i];
			((float *)(param->value))[i] = ((float *)pf_value)[map - 1];
		}

	} else if (param->type == M_LONG) {
		for (i = 0; i < param->size; i++) {
			map = ((int *)(mapping_param->value))[i];
			((int *)(param->value))[i] = ((int *)pf_value)[map - 1];
		}

	} else if (param->type == M_STRING) {
		for (i = 0; i < param->size; i++) {
			map = ((int *)(mapping_param->value))[i];
			*((char **)param->value + i) = strdup (pf_value + map - 1);
		}
	}
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : open_parameter_file
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : error message otherwise NULL
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *open_parameter_file () {
    param_file = NULL;
	if ((param_file = fopen (file_name, "r")) == NULL) {
		return error_string("cannot open parameter file");
	}
	lineNumber = 0;
    return NULL;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : close_parameter_file
 | COMMENT		: local
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void close_parameter_file () {
	if (param_file != NULL) {
		fclose (param_file);
		param_file = NULL;
		lineNumber = -1;
	}
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : get_next_line
 | COMMENT		: Use this only to read line from parameter file
 | PARAMETERS   :
 | RETURN VALUE : Pointer to "line" string if read was successful
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *get_next_line () {
	char *line_p;

/* DANGER -- could have overflowed if max_data_ln_len was smaller than MAXDATALNLEN
*/
	if (line == NULL) {
        if (max_data_ln_len < MAXDATALNLEN) {
		   line = (char *) umalloc(MAXDATALNLEN * sizeof(char));
        } else {
		   line = (char *) umalloc(max_data_ln_len * sizeof(char));
        }
	}

	if ((line_p = fgets (line, MAXDATALNLEN, param_file))) {
		lineNumber++;
	}
	return line_p;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : error_string
 | COMMENT		: Generates an error message
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *error_string (char *message) {
	static char buf[256];
	snprintf (buf, 256, "ERROR: %s; file is %s; line number %d", message, file_name, lineNumber);
	return buf;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : warning_string
 | COMMENT		: Generates a warning message.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *warning_string (char *message) {
	static char buf[256];
	snprintf (buf, 256, "WARNING: %s; file is %s; line number %d", message, file_name, lineNumber);
	return buf;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : bad_param_value_l
 | COMMENT		: Generates the warning message when a long value is out of range.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void bad_param_value_l (long l, int i, char *name, long min_l, long max_l) {
	static char buf[256];

	if (*control_lvar("parameter_check_flag") > 0) {
		snprintf (buf, 256, "%s[%d] = %ld is out of range (%ld to %ld)", name, i, (long)l, (long)min_l, (long)max_l);
		fprintf (stderr, "%s\n", warning_string(buf));
	}
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : bad_param_value
 | COMMENT		: Generates the warning message when a double or float value is out of range.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void bad_param_value (double d, int i, char *name, double min_d, double max_d) {
	static char buf[256];

	if (*control_lvar("parameter_check_flag") > 0) {
		snprintf (buf, 256, "%s[%d] = %f is out of range (%f to %f)", name, i, d, min_d, max_d);
		fprintf (stderr, "%s\n", warning_string(buf));
	}
}
