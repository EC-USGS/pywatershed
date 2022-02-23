/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : read_vars
 * COMMENT  : reads the vars data base from a file.
 *            File name is passed in as an argument
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_VARS_C
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "mms.h"

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static int read_var_line (char *, char *, FILE *, char *);

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : read_vars
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int read_vars (char *var_file_name) {

	FILE *var_file;
	PUBVAR *var;
	DIMEN *dim;
	long dim_size, var_size, type, i;
	double *dvalptr;
	float *fvalptr;
	long *lvalptr;
	char line[MAXVARLEN], key[MAXVARLEN];
	char dimen[MAXVARLEN];
	char *pathname;
	char *endptr;

/*
* get var name, open file
*/
      pathname = strdup (var_file_name);

   if ((var_file = fopen (pathname, "r")) == NULL) {
      (void)fprintf(stderr, "WARNING - read_vars - cannot open file '%s'\n",
                    pathname);
      return(0);
   }

/*
* read in run info string
*/
   if (fgets(line, MAXVARLEN, var_file) == NULL) {
      fclose(var_file);
      return(0);
   }
   Mparaminfo = strdup (line);

/*
* read in last nstep
*/
   if (fgets(line, MAXVARLEN, var_file) == NULL) {
      fclose(var_file);
      return(0);
   }

   Mnsteps = strtol(&(line[11]), &endptr, 10);

/*
* read in last time step
*/
   if (fgets(line, MAXVARLEN, var_file) == NULL) {
      fclose(var_file);
      return(0);
   }

/*
   (void)fprintf (stderr,"read_vars Mnowtime->jt stirng = %s\n", &(line[17]));
   Mnowtime->jt = strtod(&(line[17]), &endptr);
   (void)fprintf (stderr,"read_vars Mnowtime->jt = %d\n", Mnowtime->jt);
*/

/*
* read in last delta time
*/
   if (fgets(line, MAXVARLEN, var_file) == NULL) {
      fclose(var_file);
      return(0);
   }

   Mdeltat = strtod(&(line[16]), &endptr);
   Mdeltanext = strtod(&(line[16]), &endptr);

/*
* read in dimensions
*/
   while (!feof(var_file)) {

/*
* space fwd to #### header
*/
   (void)strncpy(line, " ", MAXVARLEN);
   while (strncmp(line, "####", 4)) {
      if (fgets(line, MAXVARLEN, var_file) == NULL) {
         fclose(var_file);
         return(0);
      }

/*
* break if variable list starts
*/
      if(!strncmp(line, "** Variables **", strlen("** Variables **")))
         goto variables;
      }

/*
* get dimen name
*/

      if(fgets(key, MAXVARLEN, var_file) == NULL) {
         (void)fprintf(stderr, "ERROR - read_var, reading dimen name.\n");
         (void)fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
         return(1);
      }
      key[strlen(key)-1] = '\0';

      if ((dim = dim_addr(key)) == NULL) {
         (void)fprintf(stderr, "WARNING - read_vars.\n");
         (void)fprintf(stderr, "Using var file '%s'\n", pathname);
         (void)fprintf(stderr, "Dimension '%s' not declared.\n", key);
         (void)fprintf(stderr, "Variables not read from file.\n");
         fclose(var_file);
         return(0);
      } else {

/*
* get dimen size
*/
         if(fgets(line, MAXVARLEN, var_file) == NULL) {
            (void)fprintf(stderr, "ERROR - read_var, reading dimen size.\n");
            fprintf(stderr,"Early end-of-file, file '%s'\n",var_file_name);
            return(1);
         }

         errno = 0;
         dim_size = strtol(line, &endptr, 10);
         if(errno != 0) {
            (void)fprintf(stderr,
                        "ERROR - read_var, decoding size from '%s'.\n", line);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            perror(" ");
            return(1);
         }
/*
* check dimension size
*/
         if (dim->value != dim_size) {
            (void)fprintf(stderr, "WARNING - read_vars.\n");
            (void)fprintf(stderr, "Using var file '%s'\n", pathname);
            (void)fprintf(stderr, "Dimension '%s' has size %ld.\n", key, dim->value);
            (void)fprintf(stderr, "Size in var file is %ld.\n", dim_size);
            (void)fprintf(stderr, "Variables not read from file.\n");
            fclose(var_file);
            return(0);
         }
      }
   } /* while */

/*
* read in variables
*/

variables:
   while (!feof(var_file)) {

/*
* space fwd to #### header
*/
      (void)strncpy(line, " ", MAXVARLEN);
      while (strncmp(line, "####", 4)) {
         if (fgets(line, MAXVARLEN, var_file) == NULL) {
            fclose(var_file);
            return(0);
         }
      }

/*
* get key
*/
      if(fgets(key, MAXVARLEN, var_file) == NULL) {
         (void)fprintf(stderr, "ERROR - read_var, reading var key.\n");
         (void)fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
         return(1);
      }
      key[strlen(key) - 1] = '\0';

      if ((var = var_addr(key)) != NULL) {
/*
* get number of dimensions
*/
         if(fgets(line, MAXVARLEN, var_file) == NULL) {
            (void)fprintf(stderr, "ERROR - read_var, reading var ndimen.\n");
            fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
            return(1);
         }

         if((var->ndimen = atol(line)) == 0) {
            (void)fprintf(stderr,
                  "ERROR - read_var, decoding var ndimen from '%s'.\n", line);
            (void)fprintf(stderr, "Key is '%s'\n", key);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            return(1);
         }

/*
* get dimens
*/

         for (i = 0; i < var->ndimen; i++) {
            if(fgets(dimen, MAXVARLEN, var_file) == NULL) {
               (void)fprintf(stderr, "ERROR - read_var, reading var dimen.\n");
               (void)fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
               return(1);
            }
            dimen[strlen(dimen) - 1] = '\0';

            if (strcmp(dimen, "PRIVATE")) {
               if (strcmp(dimen, var->dimen[i]->name)) {
                  (void)fprintf(stderr, "ERROR - read_var, reading var dimen.\n");
                  (void)fprintf(stderr, "Expecting dimension '%s'\n", var->dimen[i]->name);
                  (void)fprintf(stderr, "Read dimension '%s'\n", dimen);
                  (void)fprintf(stderr, "Key is '%s'\n", key);
                  (void)fprintf(stderr, "File '%s'\n", var_file_name);
                  return(1);
               }
            }
         } /* i */

/*
* get var size
*/

         if(fgets(line, MAXVARLEN, var_file) == NULL) {
            (void)fprintf(stderr, "ERROR - read_var, reading var size.\n");
            (void)fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
            return(1);
         }

         errno = 0;
         var_size = strtol(line, &endptr, 10);
         if(errno != 0) {
            (void)fprintf(stderr,
                     "ERROR - read_var, decoding var size from '%s'.\n", line);
            (void)fprintf(stderr, "Key is '%s'\n", key);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            return(1);
         }
         if(var_size != var->size) {
            (void)fprintf(stderr, "ERROR - read_var, size incorrect.\n");
            (void)fprintf(stderr, "Key is '%s'\n", key);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            return(1);
         }

/*
* get type
*/
         if(fgets(line, MAXVARLEN, var_file) == NULL) {
            (void)fprintf(stderr, "ERROR - read_var, reading var type.\n");
            (void)fprintf(stderr, "Early end-of-file, file '%s'\n", var_file_name);
            return(1);
         }
         if((type = atol(line)) == 0) {
            (void)fprintf(stderr,
                  "ERROR - read_var, decoding var type from '%s'.\n", line);
            (void)fprintf(stderr, "Key is '%s'\n", key);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            return(1);
         }
         if(type != var->type) {
            (void)fprintf(stderr, "ERROR - read_var, type incorrect.\n");
            (void)fprintf(stderr, "Key is '%s'\n", key);
            (void)fprintf(stderr, "Var file '%s'.\n", var_file_name);
            return(1);
         }

/*
* read in and store the file data
*/

         switch (type) {
            case M_DOUBLE:
               dvalptr = (double *) var->value;
               for (i = 0; i < var_size; i++) {
                  if(read_var_line(key, line, var_file, var_file_name))
                     return(1);
                  dvalptr[i] = atof(line);
               }
               break;

            case M_FLOAT:
               fvalptr = (float *) var->value;
               for (i = 0; i < var_size; i++) {
                  if(read_var_line(key, line, var_file, var_file_name))
                     return(1);
                  fvalptr[i] = (float) atof(line);
               }
               break;

         case M_LONG:
            lvalptr = (long *) var->value;
            for (i = 0; i < var_size; i++) {
               if(read_var_line(key, line, var_file, var_file_name))
                  return(1);
               lvalptr[i] =  atol(line);
            }
            break;
         }
      } /* if (var ... */

   } /* while */

   fclose(var_file);

   return(0);
}

/**7****************** LOCAL FUNCTION DEFINITIONS *********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : read_var_line
 | COMMENT		: gets a line from the variable file
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static int read_var_line (char *key, char *line, FILE *var_file, char *var_file_name) {

	if (fgets(line, MAXVARLEN, var_file) == NULL) {
		(void)fprintf(stderr,
		    "ERROR - read_var, reading data.\n");
		(void)fprintf(stderr,
		    "Early end-of-file, file '%s'\n", var_file_name);
		(void)fprintf(stderr, "Key is '%s'\n", key);
		return(1);
	}

	return(0);

}
