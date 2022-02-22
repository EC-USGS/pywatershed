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

#ifndef MMS_DEF_H
#define MMS_DEF_H

#define TRUE 1
#define FALSE 0

#define ftnint int
#define ftnlen int

#define M_LONG 1
#define M_FLOAT 2
#define M_DOUBLE 3
#define M_STRING 4

#define M_PARAMETER 0
#define M_VARIABLE  1

// Markstrom 05/24/2010 found some unused defined constants and removed them.
//#define MAX_OPT_ARRAY_SIZE 200

#define M_BOUNDED 1
#define M_UNBOUNDED 2

#define ERROR_TIME 100
//#define M_NODEBUG 0
//#define M_PARTDEBUG 1
#define M_FULLDEBUG 2

#ifndef MIN
#define MIN(a,b) (a < b) ? a : b
#endif

#ifndef MAX
#define MAX(a,b) (a > b) ? a : b
#endif

#define MAXDATALNLEN 12000 /* max no. of chars in input file line */

#define ENDOFFILE 2L
#define ENDOFDATA 1L
#define NOTENDOFDATA 0l

// Markstrom 05/24/2010 found some reads that were differnt sizes that the character arrays that
//                      were being read into. I fixed this by getting rid of all of these other
//                      string lengths and changing all string lengths to "MAXDATALNLEN" since this
//                      is the length that we need to read data lines in the data file.
//#define MAXENVLEN 256    /* max env file line length */
//#define MAXTOKLEN 128    /* max no of chars in token */
#define MAXCTRLLINELEN 256 /* max length of a line in a Control File */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256   /* max no of chars in a path */
#endif

#define MAX_NDIMEN 3       /* max no. of dimensions for a var or param */
#define MAXVARLEN 32       /* max no. of chars in variable string */
//#define MAXKEYLEN 50       /* max no. of chars in key string */
//#define MAXDIMLEN 50       /* max no. of chars in dimen string */
//#define MAXSTATVARS 200     /* max no. of statistic variables */
//#define MAXDISPVARS 200    /* max no. of display variables */
//#define MAXINFOLEN 80      /* max no. of chars in run info string */

//#define MAX_SAVE_MAP 5

#endif /* MMS_DEF_H */

