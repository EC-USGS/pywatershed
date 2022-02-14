/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : batch_run_functions
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include "mms.h"

/**5*********************** LOCAL VARIABLES ***************************/
  static FILE *statvar_file;
  static FILE **ani_out_files;
  static long nstatVars, naniVars;
  static char **statVar_names, **statVar_element;
  static char **aniVar_names;
  static char statvar_path[MAXPATHLEN];
  static char ani_path[MAXPATHLEN];
  static char buf[256];
  static long i, j, stats_flag, ani_out_flag;
  static char  *err_message, *c;
  static char   err[256];
  static int       started;
  static PUBVAR    **ani_out_vars, *var;
  static DIMEN **ani_out_dims, *dim;
  static DIMEN *nhrudim, *ngwdim, *nssrdim, *foobar;
  static FILE **ani_var_files;
  static int num_ani_dims, found, k;
  static DATETIME start_of_data, end_of_data;

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_pre_init
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_pre_init () {

  stats_flag = *control_lvar ("statsON_OFF");
  if (stats_flag == 1)
     nstatVars = *control_lvar("nstatVars");
  else
     nstatVars = 0;

  if (!nstatVars)
    stats_flag = 0;

  ani_out_flag = *control_lvar ("aniOutON_OFF");
  if (ani_out_flag == 1)
     naniVars = *control_lvar("naniOutVars");
  else
     naniVars = 0;

  if (!naniVars)
    ani_out_flag = 0;

/*
**  Make sure that all of the selected stat and display variables are OK
**  before a run is attempted.
*/
  if (stats_flag) {
    err_message = CHECK_stat_vars();
    if (err_message) return (err_message);
  }

  if (ani_out_flag) {
    err_message = CHECK_ani_vars ();
    if (err_message) return (err_message);
  }

/*
* create stats vars linked list
*/
  if (stats_flag)
    create_vstats();

/*
* open data file
* ensure datainfo is up to date
*/
  err_message = DATA_read_init ();
  if (err_message) {
     (void)printf ("%s\n", err_message);
     return (err_message);
  }

/*
 * **  Reset run period to period of record if -por flag set
 */
  if (run_period_of_record) {
    DATA_find_end (&start_of_data, &end_of_data);
    julday (&start_of_data);
    julday (&end_of_data);

    printf ("resetting model start time to period of record %ld %ld %ld %ld %ld %ld\n",
            start_of_data.year, start_of_data.month, start_of_data.day,
            start_of_data.hour, start_of_data.min, start_of_data.sec);

    printf ("resetting model end time to period of record %ld %ld %ld %ld %ld %ld\n",
            end_of_data.year, end_of_data.month, end_of_data.day,
            end_of_data.hour, end_of_data.min, end_of_data.sec);


    Mendtime = &end_of_data;
    Mstrttime = &start_of_data;
  }

  err_message = DATA_check_start ();
  if (err_message) {
     (void)printf ("%s\n", err_message);
     return (err_message);
  }

/*
* Open statvar file, and store number of variables and variable names
*/
  if (stats_flag) {
    (void)snprintf(statvar_path, MAXPATHLEN, "%s", *((char **) control_var("stat_var_file")));

    if ((statvar_file = fopen(statvar_path, "w")) == NULL) {
      (void)snprintf (err, 256, "ERROR - single_run: Could not open statvar file '%s'\n",
		     statvar_path);
      return (err);
    }

    statVar_names = (char **) control_var("statVar_names");
    statVar_element = (char **) control_var("statVar_element");

/*
* write number of variables and statVars names to stats data file.
*/
    (void)fprintf(statvar_file,"%ld\n",nstatVars);

    for (i = 0; i < nstatVars; i++)
      (void)fprintf(statvar_file,"%s %s\n", statVar_names[i], statVar_element[i]);
  }

/*
* Open ani output files.
*/
  if (ani_out_flag) {
    aniVar_names = (char **) control_var("aniOutVar_names");
    (void)snprintf(ani_path, MAXPATHLEN, "%s", *((char **) control_var("ani_output_file")));

    ani_out_dims = (DIMEN **)malloc (naniVars * sizeof (DIMEN *));
    ani_var_files = (FILE **)malloc (naniVars * sizeof (FILE *));
    ani_out_vars = (PUBVAR **)malloc (naniVars * sizeof (PUBVAR *));

/*
**  Get the pubvars.
*/
    for (i = 0; i < naniVars; i++) {

       c = strchr (aniVar_names[i], '.');
         if (c)
           *c = '\0';

       snprintf (buf, 256, "%s.%s", ani_path, aniVar_names[i]);
       c = strchr (buf, ' ');
         if (c)
           *c = '\0';

       ani_out_vars[i] = var_addr (aniVar_names[i]);
    }

/*
**  List of unique ANIMATION dimensions.
*/
    num_ani_dims = 0;
    nhrudim = dim_addr("nhru");
    ngwdim = dim_addr("ngw");
    nssrdim = dim_addr("nssr");
    for (i = 0; i < naniVars; i++) {
       found = FALSE;
       
/* Reset ngw and nssr dimensions to nhru for the purposes of the animation file.*/
       foobar = ani_out_vars[i]->dimen[0];
       if (foobar == ngwdim || foobar == nssrdim) {
          foobar = nhrudim;
       }

       for (j = 0; j < num_ani_dims; j++) {
          if (foobar == ani_out_dims[j]) {
             found = TRUE;
          }
       }

       if (!found) {
          ani_out_dims[j] = foobar;
          num_ani_dims++;
       }
    }

/*
**  Open a file for each dimension.
*/
    ani_out_files = (FILE **)malloc (num_ani_dims * sizeof (FILE *));

    for (i = 0; i < num_ani_dims; i++) {
       snprintf (buf, 256, "%s.%s", ani_path, ani_out_dims[i]->name);
       if ((ani_out_files[i] = fopen(buf, "w")) == NULL) {
          (void)snprintf (err, 256, "ERROR - single_run: Could not open ani file '%s'\n", buf);
          return (err);
       }

       fprintf (ani_out_files[i], "#\n# Begin DBF\n");
       fprintf (ani_out_files[i], "# timestamp,#FIELD_ISODATE,10,0\n");
       fprintf (ani_out_files[i], "# %s,#FIELD_DECIMAL,10,2\n", ani_out_dims[i]->name);
    }

/*
**  Map each variable to a file.
*/
    for (i = 0; i < naniVars; i++) {
/* Reset ngw and nssr dimensions to nhru for the purposes of the animation file.*/
       foobar = ani_out_vars[i]->dimen[0];
       if (foobar == ngwdim || foobar == nssrdim) {
          foobar = nhrudim;
       }

       for (j = 0; j < num_ani_dims; j++) {
          if (foobar == ani_out_dims[j]) {
             ani_var_files[i] = ani_out_files[j];
          }
       }
    }

/*
**  Finish writing the headers.
*/
    for (i = 0; i < naniVars; i++)
       fprintf (ani_var_files[i], "# %s,#FIELD_DECIMAL,10,2\n", ani_out_vars[i]->name);

    for (i = 0; i < num_ani_dims; i++) {
       fprintf (ani_out_files[i], "# End DBF\n#\n");
    }

/*
**  Write variable name line
*/
    for (i = 0; i < num_ani_dims; i++) {
       fprintf (ani_out_files[i], "timestamp	%s", ani_out_dims[i]->name);
    }

    for (i = 0; i < naniVars; i++) {
       fprintf (ani_var_files[i], "	%s", ani_out_vars[i]->name);
    }

    for (i = 0; i < num_ani_dims; i++)
       fprintf (ani_out_files[i], "\n");

/*
**  Write variable size line
*/
    for (i = 0; i < num_ani_dims; i++) {
       fprintf (ani_out_files[i], "10d	10n");
    }

    for (i = 0; i < naniVars; i++) {
       fprintf (ani_var_files[i], "	10n");
    }

    for (i = 0; i < num_ani_dims; i++)
       fprintf (ani_out_files[i], "\n");

  }

/*
**  set initial values of nsteps global variable if they
**  don't come from the var init file
*/
    Mnsteps = 0;

/*
* initialize modules
*/

  errno = 0;
  return(NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_post_init
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_post_init () {

  initializeRuntimeGraphs();

/*
**  set initial values of nsteps global variable if they
**  don't come from the var init file
*/
    Mnsteps = 0;

  started = FALSE;

   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_pre_run
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_pre_run () {

	started = TRUE;

	errno = 0;
	return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_post_run
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_post_run () {

      if (stats_flag)
         write_vstats (statvar_file);

/*
**  Write the ANIMATION output variables to their respective files.
*/
      if (ani_out_flag) {

/*
**  Each dimension has it's own file.
*/
         for (i = 0; i < num_ani_dims; i++) {
            dim = ani_out_dims[i];
            for (j = 0; j < dim->value; j++) {

/*
**  Write the current time stamp to the dimension file.
*/
            fprintf (ani_out_files[i], "%4ld-%02ld-%02ld %10ld",
                     Mnowtime->year,
                     Mnowtime->month, Mnowtime->day, j + 1);

/*
**  Write the variable values to the dimension file.
*/
               for (k = 0; k < naniVars; k++) {
                  var = ani_out_vars[k];
/*                  if (var->dimen[0] == dim) {
*/
                     switch (var->type) {
                        case M_DOUBLE:
                           fprintf (ani_var_files[k], " %10.3e",
                                    *((double *) var->value + j));
                           break;

                        case M_FLOAT:
/*                           fprintf (ani_var_files[k], "\t%14.6e",
*/
                           fprintf (ani_var_files[k], " %10.3e",
                                    *((float *) var->value + j));
                           break;

                        case M_LONG:
// markstro test
//                           fprintf (ani_var_files[i], "\t%10ld",
//                                    *((long *) var->value + j));
                           fprintf (ani_var_files[i], "\t%10d",
                                    *((int *) var->value + j));
                           break;
                     }
/*
                  }
*/
               }
               fprintf (ani_out_files[i], "\n");
            }
         }
      }

      plotRuntimeGraphValue();

   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_pre_cleanup
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_pre_cleanup () {

/*
  if (GISProcessRunning)
    GISEndAnimation();

  if (ani_init) {
    RESET_animation_control ();
  }
*/

   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : single_run_post_cleanup
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *single_run_post_cleanup () {

/*
* free up stats vars linked list
*/
  if (stats_flag)
    free_vstats();

/*
* close files and tidy up
*/
  DATA_close ();
  //closeUserFiles();

  if (stats_flag)
    fclose(statvar_file);

/*
**  Close the ANIMATION output files
*/
   if (ani_out_flag) {
      for (i = 0; i < num_ani_dims; i++)
         fclose (ani_out_files[i]);

   }

  closeRuntimeGraphs();

  if (!started) {
    return ("Run period outside of data in file.");
  }

/*
* if required, save vars to file
*/
  if (preprocess_on) {
    write_preprocess_params ();
  }

   return (NULL);
}
