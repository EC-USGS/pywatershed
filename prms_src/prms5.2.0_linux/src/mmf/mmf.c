/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : xmms
 * COMMENT  : main driver for xmms
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define MAIN
#define MMF_C

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include "mms.h"


/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
extern int call_modules(char *);
extern int call_setdims(void);

/*--------------------------------------------------------------------*\
  | FUNCTION     : main
  | COMMENT		: Main source for xmms
  | PARAMETERS   :
  |     int argc -- Argument count of command line parameters.
  |     char *argv[] -- Command line parameters.
  | RETURN VALUE : None
  | RESTRICTIONS : None
  \*--------------------------------------------------------------------*/
int main (int argc, char *argv[]) {
   int set_count;
   int i, j;
   char **set_name, **set_value;
   CONTROL *cp;
   char *cptr;
   double *dptr;
   float *fptr;
   long *lptr;
   char **cpt;
   struct stat stbuf;
   char	*err;
   static int      num_param_files = 0;
   char   **fname;
   char pathname[MAXPATHLEN];
   int set_size;


    /*
	**  Maximum buffer size for reading lines from files.
	**  This used to be set as a C precompiler directive.
	**  That is still the default, but now users are give.
	**  the option to set this on the command line, otherwise
	**  size still comes from the defs.h file.
	*/
    max_data_ln_len = MAXDATALNLEN;

	/*
	**  List of modules that are used by the model. This is
	**  determined by calls to declmodule
	*/
	module_db = ALLOC_list ("Module Data Base", 0, 100);

  /*
  **	parse the command-line arguments
  */
   set_count = 0;
   set_size = 100;
   set_name = (char **)malloc (set_size * sizeof (char *));
   set_value = (char **)malloc (set_size * sizeof (char *));
   parse_args (argc, argv, &set_count, set_name, set_value, set_size);

   if (MAltContFile == NULL) {
      (void)fprintf (stderr,"Usage: Set the full path to the control file using the '-C' option.\n\n");
      exit(0);
   }

   alloc_space ();

   setup_cont ();
   	err = read_control (MAltContFile);
	if (err) {
       (void)fprintf (stderr,"%s\n", err);
        exit (1);           
   }

	fname =   control_svar ("param_file");
    num_param_files = control_var_size ("param_file");

   for (i = 0; i < set_count; i++) {
      cp = control_addr (*(set_name + i));
      if (cp) {

         (void)fprintf (stderr,"\nControl variable %s set to %s.\n\n",
                 *(set_name + i), *(set_value + i));

         cptr = (char *)strtok (*(set_value + i), ",");

         j = 0;
         while (cptr) {
            if (cp->type == M_DOUBLE) {
               dptr = (double *) cp->start_ptr;
               dptr[j] = atof(cptr);
            } else if (cp->type == M_FLOAT) {
               fptr = (float *) cp->start_ptr;
               fptr[j] = (float) atof(cptr);
            } else if (cp->type == M_LONG) {
               lptr = (long *) cp->start_ptr;
               lptr[j] =  atol(cptr);
            } else if (cp->type == M_STRING) {
               cpt = (char **) cp->start_ptr;
               cpt[j] = strdup (cptr);
            }

            cptr = (char *)strtok (NULL, ",");
            j++;
         }

      } else {
         (void)fprintf (stderr,"\nControl variable %s not found -- ignored.\n\n",
                 *(set_name + i));
      }
   }

	fname =   control_svar ("param_file");
    num_param_files = control_var_size ("param_file");

    if (call_setdims()) {
	  (void)fprintf(stderr, "\nERROR: Calling function 'call_setdims'\n");
      exit (1);
    }

    /*
    **	read dimension info from parameter file
    */
    if (stat (*control_svar("param_file"), &stbuf) != -1) {
       if (stbuf.st_size) {
      } else {
	     (void)fprintf (stderr, "Parameter File: %s is empty.",
		               *control_svar("param_file"));
        exit (1);
	   }
    }
    
    err = read_dims (*control_svar("param_file"));
    if (err) {
//		(void)fprintf (stderr,"\nERROR: reading dimensions from Parameter File\n");
		fprintf (stderr,"\n%s\n", err);
        exit (1);
	}

	fname =   control_svar ("param_file");
    num_param_files = control_var_size ("param_file");

    if (call_modules("declare")) {
		(void)fprintf(stderr, "\nERROR: in declare procedure, in function 'call_modules'\n");
        exit (1);
	}
    
    /*
    **	read in parameter values from parameter file
    */
	fname =   control_svar ("param_file");
    num_param_files = control_var_size ("param_file");

	/*
	**  Look for, declare and read in mapping parameters before any of the "module" parameters
	*/
	for (i = 0; i < num_param_files; i++) {
		if (stat (fname[i], &stbuf) != -1) {
		   if (stbuf.st_size) {
		  } else {
			  (void)fprintf (stderr, "ERROR: Parameter file: %s is empty.",
						   fname[i]);
			  exit (1);
		   }
		}
	    
		err = read_params (fname[i], i, 1);
		if (err) {
			(void)fprintf (stderr,"\n%s\n", err);
			exit (1);
		}
	}

	/*
	**  Read in the parameters declared by the modules.
	*/

	for (i = 0; i < num_param_files; i++) {
		if (stat (fname[i], &stbuf) != -1) {
		   if (stbuf.st_size) {
		  } else {
			  (void)fprintf (stderr, "ERROR: Parameter file: %s is empty.",
						   fname[i]);
			  exit (1);
		   }
		}
	    
		err = read_params (fname[i], i, 0);
		if (err) {
			(void)fprintf (stderr,"\n%s\n", err);
			exit (1);
		}
	}
    
    /*
    **  get data info string into the global
    */
    err = READ_data_info ();
    if (err) (void)fprintf (stderr,"\nMMS - Warning: %s", err);

    /*
    **	get start and end time
    */
    get_times ();
    
    if (print_mode) {
      print_params();
      print_vars();
      print_model_info();
	  (void)snprintf (pathname, MAXPATHLEN, "%s.param", MAltContFile);
	  save_params (pathname);

    } else {

      BATCH_run ();
      ;
    }

    exit (0);
}
