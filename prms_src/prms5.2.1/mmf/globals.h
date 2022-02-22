/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION :
 * COMMENT  :
 *
 *  $Id$
 *
-*/

#ifndef MMS_GLOBAL_H
#define MMS_GLOBAL_H

#ifdef MAIN

/*
 * for MAIN only
 */

int max_data_ln_len;                /* now possible to set this on command line */
char *MAltContFile = NULL;          /* Alt. name of control file */
long Mdebuglevel = 0;               /* the current debug level */
char *model_name = NULL;
char *executable_model = NULL;
int batch_run_mode = FALSE;         /* flag for running in batch mode  */
int run_period_of_record = FALSE;   /* flag for running entire period of
                                            record in batch mode  */
int print_mode = FALSE;
int runtime_graph_on = FALSE;
int preprocess_on = FALSE;         /* flag for running in preprocess mode */
LIST *cont_db;
LIST *dim_db;
LIST *module_db;
MODULE_DATA *current_module;
PUBVAR **Mvarbase = NULL;           /* pointer to public variables data base */
long Mnvars = 0;                    /* no of public variables in data base */
PARAM **Mparambase = NULL;          /* pointer to parameter data base */
long Mnparams = 0;                  /* no of params in data base */
READCHECK **Mcheckbase = NULL;      /* pointer to read check data base */
long Mnreads = 0;                   /* max no of calls to be made by readvar */
DATETIME *Mstrttime = NULL;         /* pointer to start time structure */
DATETIME *Mendtime = NULL;          /* pointer to end time structure */
DATETIME *Mnowtime = NULL;          /* pointer to current data time structure */
DATETIME *Mnexttime = NULL;         /* pointer to next data time structure */
char *Mparaminfo = NULL;            /* pointer to param information string */
char *Mdatainfo = NULL;             /* pointer to data information string */
PARAM **unsort_params = NULL;       /* pointer to unsorted parameters */
FILE_DATA   **fd = NULL;
long Mnsteps = 0;                   /* the number of steps so far */
double Mprevjt = -1.0;              /* the latest previous Julian time  */
double Mdeltat = 0.0;               /* the latest time step in hours */
char *Minpptr = NULL;               /* pointer to current posn in data input line*/
double Mdeltanext = 0.0;            /* the latest next time step in hours */
int  M_stop_run = 0;                /* Run switch 0 -or 1 */
STAT_LIST_TYPE *Mfirst_stat_list = NULL;     /* pointer to first entry
						in stats link list */
char *Mtypes[] = {"", "long", "float", "double", "string", "", "","", "", ""};
long ParamBaseIsDirty = FALSE;
int max_vars;
int max_params;
int max_read_vars;
int max_dims;
int max_controls;

#else

/*
 * for all functions except main
 */

extern int max_data_ln_len;
extern char *MAltContFile;
extern long Mdebuglevel;
extern char *model_name;
extern char *executable_model;
extern int batch_run_mode;
extern int run_period_of_record;
extern int print_mode;
extern int runtime_graph_on;
extern int preprocess_on;
extern LIST *cont_db;
extern LIST *dim_db;
extern LIST *module_db;
extern MODULE_DATA *current_module;
extern PUBVAR **Mvarbase;
extern long Mnvars;
extern PARAM **Mparambase;
extern long Mnparams;
extern READCHECK **Mcheckbase;
extern long Mnreads;
extern DATETIME *Mstrttime;
extern DATETIME *Mendtime;
extern DATETIME *Mnowtime;
extern DATETIME *Mnexttime;
extern char *Mparaminfo;
extern char *Mdatainfo;
extern PARAM **unsort_params;
extern FILE_DATA   **fd;
extern long Mnsteps;
extern double Mprevjt;
extern double Mdeltat;
extern char *Minpptr;
extern double Mdeltanext;
extern int  M_stop_run;
extern STAT_LIST_TYPE *Mfirst_stat_list;
extern char *Mtypes[];
extern long ParamBaseIsDirty;
extern int max_vars;
extern int max_params;
extern int max_read_vars;
extern int max_dims;
extern int max_controls;

#endif /* MAIN */

#endif /* MMS_GLOBAL_H */
