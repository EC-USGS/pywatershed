/*+
** DANGER looks for spaces before comment and sets them to null so
** data part of the line does not have any part of the comment
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : read_control
 * COMMENT  : reads the control data base from a file
 *            File name is obtained from the environment variable "mms_control_file"
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_CONTROL_C
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static char *rc (char *);
char *fgets_rc (char *, int , FILE *);

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_control
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_control (char *control_name) {
   static char *foo = NULL;
   char old[256], *cptr;

   if (foo) {
      strncpy (old, foo, 256);
      free (foo);
      foo = strdup (control_name);
   } else {
      strncpy (old, control_name, 256);
      foo = strdup (control_name);
   }

   cptr = rc (control_name);

   if (cptr) {
      rc (old);

      free (foo);
      foo = strdup (old);
   }

   return (cptr);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : rc
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static char *rc (char *control_name) {
   FILE   *control_file;
   CONTROL *cp;
   long   type, size, i;
   double   *dptr;
   float   *fptr;
   long   *lptr;
   char   line[MAXCTRLLINELEN], *key;
   static char      buf[256];

/*
* compute control path, open file
*/
   if ((control_file = fopen (control_name, "r")) == NULL) {
      (void)snprintf (buf, 256, "read_control: Couldn't open %s", control_name);
      return (buf);
   }

   if (!fgets_rc(line, MAXCTRLLINELEN, control_file)) {
      fclose (control_file);
      (void)snprintf (buf, 256, "read_control: Problems reading %s", control_name);
      return (buf);
   }

/*
**   Read in control variables.
*/
   while (!feof (control_file)) {

/*
**   Space fwd to #### header.
*/
      while (strncmp(line, "####", 4)) {
         if (fgets_rc(line, MAXCTRLLINELEN, control_file) == NULL) {
            fclose(control_file);
            return(NULL);
         }
      }
/*
**   get key
*/
      if (!fgets_rc (line, MAXCTRLLINELEN, control_file)) {
         (void)snprintf (buf, 256, "read_control: reading key; Early end-of-file");
         printf ("read_control: reading key; Early end-of-file\n");
         return (buf);
      }

	  /* Replace the end of line with a null */
      *(line + strlen (line) - 1) = '\0';
	  key = strdup (line);


/*
**   get size
*/
      if (!fgets_rc (line, MAXCTRLLINELEN, control_file)) {
         (void)snprintf (buf, 256, "read_control: reading size; key = %s", key);
         return (buf);
      }

      if ((size = atol(line)) < 0) {
         (void)snprintf (buf, 256, "read_control: negative size; key = %s, line = %s", key, line);
         return (buf);
      }

/* DANGER check that size is within some range to prevent overflow.
** This is a hack, but 1000 should be a good upper limit on the number of indexes for a control variable.
*/
      if (size > 999) {
         (void)snprintf (buf, 256, "read_control: too many control indexes; size = %ld, key = %s, line = %s", size, key, line);
         return (buf);
      }

/*
**   get type
*/
      if (!fgets_rc (line, MAXCTRLLINELEN, control_file)) {
         (void)snprintf (buf, 256, "WARNING: reading type; key = %s", key);
         return (buf);
      }

      if (!(type = atol(line))) {
         (void)snprintf (buf, 256, "WARNING: invalid type; key = %s, line = %s", key, line);
         return (buf);
      }

      cp = control_addr (key);
      if (!cp) {
         cp = add_control (key, type, size); // This is if the control variable was not set in the setupcont function
//	  printf ("   read_control E %s NOT FOUND in SETUPCONT\n", key);
     }
     
	  if (cp->set_in_file > 0) {
		   printf ("\n\nWARNING: %s is duplicated in the control file %s.\n\n", key, control_name);
	  }

//  Set the values to what was just read from the file
      cp->key = strdup(key);
      cp->type = type;
      cp->size = size;
      cp->set_in_file = 1;

      switch (type) {
         case M_DOUBLE:
			dptr = (double *)umalloc (sizeof (double) * size);
            cp->start_ptr = (void *)dptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXCTRLLINELEN, control_file) == NULL) {
                  (void)snprintf (buf, 256, "read_control: key is %s.\n, file: %s", key, control_name);
                  printf ("read_control CRASH reading control file: key is %s.\n, file: %s\n", key, control_name);
                  return (buf);
               }
               dptr[i] = atof(line);
            }
            break;

         case M_FLOAT:
			fptr = (float *)umalloc (sizeof (float) * size);
            cp->start_ptr = (void *)fptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXCTRLLINELEN, control_file) == NULL) {
                  (void)snprintf (buf, 256, "read_control: key is %s.\n, file: %s", key, control_name);
                  printf ("read_control CRASH reading control file: key is %s.\n, file: %s\n", key, control_name);
                  return (buf);
               }
               fptr[i] = (float) atof(line);
            }
            break;

         case M_LONG:
			lptr = (long *)umalloc (sizeof (long) * size);
            cp->start_ptr = (void *)lptr;
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXCTRLLINELEN, control_file) == NULL) {
                  (void)snprintf (buf, 256, "read_control: key is %s.\n, file: %s", key, control_name);
                  printf ("read_control CRASH reading control file: key is %s.\n, file: %s\n", key, control_name);
                  return (buf);
               }
               lptr[i] =  atol(line);
            }
            break;

         case M_STRING:
			cp->start_ptr = umalloc (sizeof (char *) * size);
            for (i = 0; i < size; i++) {
               if (fgets_rc(line, MAXCTRLLINELEN, control_file) == NULL) {
                  (void)snprintf (buf, 256, "read_control: key is %s.\n, file: %s", key, control_name);
                  printf ("read_control CRASH reading control file: key is %s.\n, file: %s\n", key, control_name);
                  return (buf);
               }
               line[strlen(line)-1] = '\0';
               *((char **)cp->start_ptr + i) = strdup (line);

/*			   printf ("read_control just put in string value %s\n", *((char **)cp->start_ptr + i));*/
            }
            break;
      }
   }
   fclose (control_file);
   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : fgets_rc
 | COMMENT      : replacement in read_control functions for fgets which
 |                stripps off comments.
 | PARAMETERS   :
 | RETURN VALUE :
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *fgets_rc (char *str, int num, FILE *stream) {
   char *ptr, *ptr2;
   // Four situations: (1) Line is blank, (2) line starts with //,
   //   (3) line has info and contains //, (4) line has info and
   //   does not contain //

   ptr = fgets(str, num, stream);
   if (!ptr) return ptr;

/*
**  A line that starts with "//" is a comment (as far as MMS is concerned).
*/
      if ((str[0] == '/') && (str[1] == '/')) {
         return fgets_rc (str, num, stream);

/*
**  Ignore blank lines
*/
      } else if (strlen (str) <= 1) {
         return fgets_rc (str, num, stream);
/*
** Assume anything else is a data line
*/

      } else if (strstr (str, "//")) {
/*
** comment in data line
** DANGER looks for spaces before comment and sets them to null so
** data part of the line does not have any part of the comment
*/
         ptr2 = strstr (str, "//");

/* New -- terminate data part of line where the comment starts
*/
         *(ptr2) = '\0';
         
/* Now look for spaces before the comment deliminator. If there
** are some, null them out.
*/

         ptr2--;
         while (*ptr2 != ' ') ptr2--;
         *(ptr2 + 1) = '\0';
         return ptr;

      } else {
         return ptr;
      }
}
