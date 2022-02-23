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

#ifndef _STRUCTS_H
#define _STRUCTS_H

#include <stdio.h>

typedef struct dimen_t{
  char *name;
  long value;
  long max;
  char *descr;
  char **names;
  char **notes;
  FILE **files;
  long column_width;
  char *format;
  int fixed;
  int got;
} DIMEN;             /* dimension pointer structure */

typedef struct {
  char *key;
  char *module;
  char *name;
  long ndimen;
  long pf_ndimen;
  struct dimen_t **dimen;
  char **pf_dimNames;
  long size;
  long pf_size;
  long type;
  long bound_status;
  struct dimen_t *bound_dimen;
  char *value;
  char *min;
  char *max;
  char *def;
  char *descr;
  char *help;
  char *units;
  char *format;
  long column_width;
  //char **value_desc;
  char *value_string;
  char *min_string;
  char *max_string;
  char *def_string;
  long read_in;
  void **references;
  long num_references;
  long size_references;
  long preprocess;
} PARAM;                 /* parameter pointer structure */

typedef struct {
  char *key;
  long size;
  long type;
  char *start_ptr;
  long set_in_file;
} CONTROL;                 /* control variable pointer structure */

typedef struct {
  long year, month, day, hour, min, sec, jd;
  double jt;
} DATETIME;                 /* date and time structure */

typedef struct list_t {
    char        *name;
    int         size;
    int         count;
    int         type;
    int         out;
    void        *user_data;
    void        **itm;
} LIST;

typedef struct {
  char *key;
  char *module;
  char *name;
  long ndimen;
  struct dimen_t **dimen;
  long size;
  long type;
  char *help;
  char *units;
  char *value;
  int private;
} PUBVAR;                 /* public variable pointer structure */

typedef struct {
  PUBVAR *var;
  long count;
  union {
    long   * valuel;
    float  * valuef;
    double * valued;
  }Types;
} READCHECK; /* for checking the readvar function calls */

typedef struct file_data_t {
	FILE    *fp;
	char    *name;
//	char    line[MAXDATALNLEN];
	char    *line;
	char    *start_of_data;
	float   delta_t;
//	char    info[MAXDATALNLEN];
	char    *info;
	DATETIME    time;
} FILE_DATA;

typedef struct STAT_LIST_TYPE {
//  char key[MAXDATALNLEN];
	char *key;
	char *element;
    long type;
    char *value;
    struct STAT_LIST_TYPE *next;
} STAT_LIST_TYPE;   /* linked list element of stat vars */

typedef struct module_data_t {
	char    *name;
	char    *version;
	LIST    *params;
	LIST    *vars;
} MODULE_DATA;

#endif
