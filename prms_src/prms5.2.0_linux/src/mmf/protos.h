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

/**1************* SWITCH FOR DEFINITION AND DECLARATION ***************/
#ifndef MSYS_PROTO_H
#define MSYS_PROTO_H

/**5**************** DECLARATION EXPORTED FUNCTIONS *******************/

/***  mmf.c  **************************************************/
extern long setdims_ (void);
extern long call_modules_ (char *, int);

/***  alloc_space.c  **************************************************/
#undef EXTERN
#ifdef ALLOC_SPACE_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void alloc_space (void);

/***  check_vars.c  **************************************************/
#undef EXTERN
#ifdef CHECK_VARS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *CHECK_stat_vars (void);
EXTERN char *CHECK_disp_vars (void);
EXTERN char *CHECK_ani_vars (void);
EXTERN char *CHECK_map_vars (void);

/***  create_vstats.c  **************************************************/
#undef EXTERN
#ifdef CREATE_VSTATS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void create_vstats (void);

///***  dprint.c  **************************************************/
//#undef EXTERN
//#ifdef DPRINT_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN void dpstr_ (char *, ftnint *, ftnlen);
//EXTERN void dpstr (char *, long);
//EXTERN void dpint4_ (char *, ftnint *, ftnint *, ftnint *, ftnlen);
//EXTERN void dplong (char *, long *, long, long);
//EXTERN void dpreal_ (char *, float *, ftnint *, ftnint *, ftnlen);
//EXTERN void dpfloat (char *, float *, long, long);
//EXTERN void dpdble_ (char *, double *, ftnint *, ftnint *, ftnlen);
//EXTERN void dpdble (char *, double *, long, long);

/***  free_vstats.c  **************************************************/
#undef EXTERN
#ifdef FREE_VSTATS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void free_vstats (void);

/***  get_times.c  **************************************************/
#undef EXTERN
#ifdef GET_TIMES_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void get_times (void);

/***  get_elem_add.c  **************************************************/
#undef EXTERN
#ifdef GET_ELEM_ADD_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int CheckIndices (char *, char *, int);
EXTERN char *GetElemAddress (char *, char *, int);

/***  oprint.c  **************************************************/
#undef EXTERN
#ifdef OPRINT_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void opstr_ (char *, ftnlen);
EXTERN void opstr (char *);
//EXTERN void opint4_ (char *, ftnint *, ftnint *, ftnlen);
//EXTERN void oplong (char *, long *, long);
//EXTERN void opreal_ (char *, float *, ftnint *, ftnlen);
//EXTERN void opfloat (char *, float *, long);
//EXTERN void opdble_ (char *, double *, ftnint *, ftnlen);
//EXTERN void opdble (char *, double *, long);

///***  rosopt.c  **************************************************/
//#undef EXTERN
//#ifdef ROSOPT_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN char *rosopt (ROSEN_DATA *, float[], float[]);

///***  opinit.c  **************************************************/
//#undef EXTERN
//#ifdef OPINIT_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN char *opinit (float *, float *, int *, ROSEN_DATA *);

/***  bdry.c  **************************************************/
#undef EXTERN
#ifdef BDRY_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int bdry (float *, int *, float *, float *, int *, int *, float *);

///***  param.c  **************************************************/
//#undef EXTERN
//#ifdef PARAM_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN int param (int *, ROSEN_DATA *);

/***  coropt.c  **************************************************/
#undef EXTERN
#ifdef COROPT_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int coropt (int *, float *, float *, int *);

///***  sub1.c  **************************************************/
//#undef EXTERN
//#ifdef SUB1_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN int sub1 (int *, float *, int *, float *, float *, float *,
//      float *, float *, int *, int *, float *,
//      int *, float *, float *, int *, int *, ROSEN_DATA *);

/***  tcale.c  **************************************************/
#undef EXTERN
#ifdef TCALE_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int tcale (int *, float *, float *, float *, float *, int *);

/***  unscal.c  **************************************************/
#undef EXTERN
#ifdef UNSCAL_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int unscal (int *, float *, float *, float *, float *, int *);

///***  sub3.c  **************************************************/
//#undef EXTERN
//#ifdef SUB3_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN int sub3 (int *, float *, float *, float *, int *, float *, ROSEN_DATA *);

/***  parse_args.c  **************************************************/
#undef EXTERN
#ifdef PARSE_ARGS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void parse_args (int, char **, int *, char **, char **, int);

/***  read_datainfo.c  **************************************************/
#undef EXTERN
#ifdef READ_DATAINFO_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *read_datainfo (FILE_DATA *);

/***  read_line.c  **************************************************/
#undef EXTERN
#ifdef READ_LINE_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long read_line (void);
EXTERN char *DATA_read_init (void);
EXTERN char *DATA_check_start (void);
EXTERN void DATA_close (void);
EXTERN int control_var_size (char *);
EXTERN FILE_DATA *FILE_with_next_ts (void);
EXTERN char *EXTRACT_time (FILE_DATA *);
EXTERN int CHECK_data (int, FILE_DATA *);
EXTERN void DATA_find_end (DATETIME *, DATETIME *);
EXTERN char *READ_data_info (void);

/***  setup_cont.c  **************************************************/
#undef EXTERN
#ifdef SETUP_CONT_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void setup_cont (void);

/***  batch_run.c  **************************************************/
#undef EXTERN
#ifdef BATCH_RUN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int BATCH_run (void);

/***  graph_single_run.c  **************************************************/
#undef EXTERN
#ifdef GRAPH_SINGLE_RUN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int initializeRuntimeGraphs (void);
EXTERN int plotRuntimeGraphValue (void);
EXTERN int closeRuntimeGraphs (void);

/***  stats.c  **************************************************/
#undef EXTERN
#ifdef STATS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int stats (void);

///***  uprint.c  **************************************************/
//#undef EXTERN
//#ifdef UPRINT_C
//#define EXTERN
//#else
//#define EXTERN extern
//#endif
//
//EXTERN FILE *GetUserFile (char *, long);
//EXTERN void closeUserFiles (void);
//EXTERN void upstr_ (char *, ftnint *, char *, ftnlen, ftnlen);
//EXTERN void upstr (char *, long, char *);
//EXTERN void upint4_ (char *, ftnint *, char *, ftnint *, ftnint *, ftnlen, ftnlen);
//EXTERN void uplong (char *, long, char *, long *, long);
//EXTERN void upreal_ (char *, ftnint *, char *, float *, ftnint *, ftnlen, ftnlen);
//EXTERN void upfloat (char *, long, char *, float *, long);
//EXTERN void updble_ (char *, ftnint *, char *, double *, ftnint *, ftnlen, ftnlen);
//EXTERN void updble (char *, long, char *, double *, long);

/***  write_vstats.c  **************************************************/
#undef EXTERN
#ifdef WRITE_VSTATS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void write_vstats (FILE *);

/***  julconvert.c  **************************************************/
#undef EXTERN
#ifdef JULCONVERT_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN double getjulday(int, int, int, int, int, double);
EXTERN int dayofweek(double);
EXTERN long isleap_ (ftnint *);
EXTERN int isleap (int);

/***  build_lists.c  **************************************************/
#undef EXTERN
#ifdef BUILD_LISTS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN LIST *ALLOC_list (char *, int, int);
EXTERN void RESIZE_list (LIST *, int);
EXTERN void DELETE_list (LIST *);
EXTERN void ADD_to_list (LIST *, void *);

/***  sensitivity.c  **************************************************/
#undef EXTERN
#ifdef SENSITIVITY_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *PRMS_sens (void);
EXTERN int IN_obj_period (int, int, int);

/***  control_addr.c  **************************************************/
#undef EXTERN
#ifdef CONTROL_ADDR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN CONTROL *control_addr (char *);

/***  control_array.c  **************************************************/
#undef EXTERN
#ifdef CONTROL_ARRAY_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *control_array (char *, long);
EXTERN long *control_larray (char *, long);
EXTERN float *control_farray (char *, long);
EXTERN double *control_darray (char *, long);
EXTERN char *control_sarray (char *, long);

/***  control_var.c  **************************************************/
#undef EXTERN
#ifdef CONTROL_VAR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *control_var (char *);
EXTERN long *control_lvar (char *);
EXTERN float *control_fvar (char *);
EXTERN double *control_dvar (char *);
EXTERN char **control_svar (char *);

/***  decldim.c  **************************************************/
#undef EXTERN
#ifdef DECLDIM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long decldim_ (char *, ftnint *, ftnint *, char *, ftnlen, ftnlen);
EXTERN long decldim (char *, long, long, char *);
EXTERN long declfix (char *, long, long, char *);
EXTERN long declfix_ (char *, ftnint *, ftnint *, char *, ftnlen, ftnlen);

/***  declmodule.c    **************************************************/
EXTERN long declmodule (char *, char *, char *);

/***  decl_control.c  **************************************************/
#undef EXTERN
#ifdef DECL_CONTROL_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void decl_control (char *, long, long, void *);
EXTERN CONTROL *add_control (char *, long, long);

/***  declparam.c  **************************************************/
#undef EXTERN
#ifdef DECLPARAM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long declparam_ (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
EXTERN long declparam (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *);

EXTERN long declparam_u_ (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *, char *, long *,
       	ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
EXTERN long declparam_u (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *, char *, long *);

EXTERN long declparam_p_ (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *, char *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
EXTERN long declparam_p (char *, char *, char *, char *, char *,
	char *, char *, char *, char *, char *, char *);

/***  declvar.c  **************************************************/
#undef EXTERN
#ifdef DECLVAR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long declvar_ (char *, char *, char *, ftnint *, char *,
	char *, char *, char *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
EXTERN long declvar (char *, char *, char *, long, char *,
	char *, char *, char *);
EXTERN long declpri_ (char *, ftnint *, char *, char *, ftnlen, ftnlen);
EXTERN long declpri (char *, long, char *, char *);

/***  dim_addr.c  **************************************************/
#undef EXTERN
#ifdef DIM_ADDR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN DIMEN *dim_addr (char *);
EXTERN char *dim_notes (char *);

/***  getdim.c  **************************************************/
#undef EXTERN
#ifdef GETDIM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long getdim_ (char *, ftnlen);
EXTERN long getdim (char *);

/***  getdimname.c  **************************************************/
#undef EXTERN
#ifdef GETDIMNAME_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void getdimname_ (char *, ftnint *, char *, ftnlen, ftnlen);
EXTERN void getdimname (char *, long, char *, int);
EXTERN void getdimdesc_ (char *, ftnint *, char *, ftnlen, ftnlen);
EXTERN void getdimdesc (char *, long, char *, int);

/***  getparam.c  **************************************************/
#undef EXTERN
#ifdef GETPARAM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long updateparam (char *);
EXTERN long getparam_ (char *, char *, ftnint *, char *, double *, ftnlen, ftnlen, ftnlen);
EXTERN long getparam (char *, char *, int, char *, double *);
EXTERN long getdatainfo_ (char *, ftnlen);
EXTERN long getdatainfo (char *, ftnlen);
EXTERN long getoutname_ (char *, char *, ftnlen, ftnlen);
EXTERN long getoutname (char *, int, char *);
EXTERN long getdataname_ (char *, char *, ftnlen, ftnlen);
EXTERN long getdataname (char *, int, char *);
EXTERN long getoutdirfile_ (char *, char *, ftnlen, ftnlen);
EXTERN long getoutdirfile (char *, int, char *);
EXTERN long getuserdirfile_ (char *, char *, ftnlen, ftnlen);
EXTERN long getuserdirfile (char *, int, char *);

/***  getvar.c  **************************************************/
#undef EXTERN
#ifdef GETVAR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long getvar_ (char *, char *, ftnint *, char *, double *,
	ftnlen, ftnlen, ftnlen);
EXTERN long getvar (char *, char *, long, char *, double *);

/***  load_param.c  **************************************************/
#undef EXTERN
#ifdef LOAD_PARAM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long load_param (PARAM *);

/***  param_addr.c  **************************************************/
#undef EXTERN
#ifdef PARAM_ADDR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN PARAM *param_addr (char *);

/***  print_params.c  **************************************************/
#undef EXTERN
#ifdef PRINT_PARAMS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int print_params (void);

/***  print_vars.c  **************************************************/
#undef EXTERN
#ifdef PRINT_VARS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int print_vars (void);

/***  putvar.c  **************************************************/
#undef EXTERN
#ifdef PUTVAR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long putvar_ (char *, char *, ftnint *, char *, double *value,
	     ftnlen, ftnlen, ftnlen);
EXTERN long putvar (char *, char *, long, char *, double *);

/***  read_control.c  **************************************************/
#undef EXTERN
#ifdef READ_CONTROL_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *read_control (char *);

/***  read_params.c  **************************************************/
#undef EXTERN
#ifdef READ_PARAMS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *read_params (char *, int, int);
EXTERN char *read_dims (char *);

/***  read_vars.c  **************************************************/
#undef EXTERN
#ifdef READ_VARS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int read_vars (char *);

/***  readvar.c  **************************************************/
#undef EXTERN
#ifdef READVAR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long readvar_ (char *, char *, ftnlen, ftnlen);
EXTERN long readvar (char *, char *);

/***  reset_dim.c  **************************************************/
#undef EXTERN
#ifdef RESET_DIM_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void reset_dim (DIMEN *, long);

/***  save_control.c  **************************************************/
#undef EXTERN
#ifdef SAVE_CONTROL_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void save_control (char *);

/***  save_params.c  **************************************************/
#undef EXTERN
#ifdef SAVE_PARAMS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int save_params (char *);
EXTERN int write_preprocess_params (void);

/***  save_vars.c  **************************************************/
#undef EXTERN
#ifdef SAVE_VARS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int save_vars (char *);

/***  sort_dims.c  **************************************************/
#undef EXTERN
#ifdef SORT_DIMS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void sort_dims (void);

/***  sort_params.c  **************************************************/
#undef EXTERN
#ifdef SORT_PARAMS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void sort_params (void);

/***  sort_vars.c  **************************************************/
#undef EXTERN
#ifdef SORT_VARS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void sort_vars (void);

/***  str_to_vals.c  **************************************************/
#undef EXTERN
#ifdef STR_TO_VALS_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN long str_to_vals (char *, long, long, char *);

/***  timing.c  **************************************************/
#undef EXTERN
#ifdef TIMING_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void dattim_ (char *, ftnint *, ftnlen);
EXTERN void dattim (char *, long *);
EXTERN long julian_ (char *, char *, ftnlen, ftnlen);
EXTERN long julian (char *, char *);
EXTERN double deltim (void);
EXTERN double delnex (void);
EXTERN long getstep_ (void);
EXTERN long getstep (void);
EXTERN double djulian_ (char *, char *, ftnlen, ftnlen);
EXTERN double djulian (char *, char *);

/***  var_addr.c  **************************************************/
#undef EXTERN
#ifdef VAR_ADDR_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN PUBVAR *var_addr (char *);

/***  esp_batch_run.c  **************************************************/
#undef EXTERN
#ifdef ESP_BATCH_RUN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *ESP_batch_run (void);

/***  rosenbrock_batch_run.c  **************************************************/
#undef EXTERN
#ifdef ROSENBROCK_BATCH_RUN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void ROSENBROCK_batch_run (void);

/***  umalloc_etc.c  **************************************************/
#undef EXTERN
#ifdef UMALLOC_ETC_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN char *umalloc (unsigned);
EXTERN char *urealloc (char *, unsigned);
EXTERN char *ucalloc (unsigned, unsigned);
EXTERN void ufree (char *);

/***  julday.c  **************************************************/
#undef EXTERN
#ifdef JULDAY_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int julday (DATETIME *);

/***  matinv.c  **************************************************/
#undef EXTERN
#ifdef MATINV_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int matinv (float *, float *, int *, int *);

/***  matind.c  **************************************************/
#undef EXTERN
#ifdef MATIND_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int matind (double *, double *, int *);

/***  snort.c  **************************************************/
#undef EXTERN
#ifdef SNORT_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int snort (float *, int *, int *, int *);

/***  print_model_info.c  **************************************************/
#undef EXTERN
#ifdef PRINT_MODEL_INFO_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int print_model_info (void);

#endif
