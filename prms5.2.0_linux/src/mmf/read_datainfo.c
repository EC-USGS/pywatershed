/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : read_datainfo
 * COMMENT  : reads the data file and updates the
 *            datainfo string and the data variable names and sizes
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_DATAINFO_C
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/*--------------------------------------------------------------------*\
 | FUNCTION     : read_datainfo
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *read_datainfo (FILE_DATA *fd) {

   static char   err_buf[512];
   long   count, nline = 0;
   static char   *line = NULL;
   static char *linecopy = NULL;
   PUBVAR   *var;
   char   *key, *countstr, *endptr;

   Mnreads = 0;

   if (line == NULL) {
	   line = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

   if (linecopy == NULL) {
	   linecopy = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

   if (!(fgets (fd->info, max_data_ln_len, fd->fp))) {
      (void)snprintf (err_buf, 512, "Can't read data file info string\n%s", fd->name);
      return (err_buf);
   }

   fd->info[strlen (fd->info) - 1] = '\0';
   nline++;

/*
**  Read the header of the data file.  The end of the header occures
**  when a line starts with at least 4 "#"s.
*/
   (void)strncpy (line, "", max_data_ln_len);
   while (strncmp (line, "####", 4)) {
      if ((fgets (line, max_data_ln_len, fd->fp)) == NULL) {
         (void)snprintf (err_buf, 512, "#### delimiter not found in data file\n%s", fd->name);
         return (err_buf);
      }

      nline++;
      if (strncmp(line, "####", 4)) {

/*
**  A line that starts with "//" is a comment (as far as MMS is concerned).
*/
         if ((line[0] == '/') && (line[1] == '/')) {
/*
            printf ("Comment: %s\n", line);
*/
/*
**  Ignore blank lines
*/
         } else if (strlen (line) <= 1) {

/*
** Assume anything else is a data line
*/
         } else {

/* 
**   Increase size of array pointers
*/
            if (Mnreads >= max_read_vars) {
               max_read_vars += 50;
               Mcheckbase = (READCHECK **)realloc(Mcheckbase, max_read_vars * 
                  sizeof(READCHECK *));
            }
      
/*
**    get key, check the var has been declared
*/

            (void)strncpy(linecopy, line, max_data_ln_len);
            key = strtok(linecopy, " \t");

            if (key == NULL) {
               (void)snprintf (err_buf, 512, "Check format at line number %ld in\n%s\n%s", nline,
                  fd->name, line);
               return (err_buf);
            }

/*
**   get size of var in data file
*/
            countstr = strtok(NULL, " \t");

            if (countstr == NULL) {
               (void)snprintf (err_buf, 512, "Check format at line number %ld in\n%s\n%s", nline,
                  fd->name, line);
               return (err_buf);
            }

            errno = 0;
            count = strtol (countstr, &endptr, 10);

				if (count > 0) {  // Old style Data files have variables with size 0. PRMS should now skip over these, as if they are not there.

					if ((var = var_addr(key)) == NULL) {
						(void)snprintf (err_buf,  512,
							"Variable %s not declared at line number %ld in\n%s\n%s",
							key, nline, fd->name, line);
						return (err_buf);
					}

	/*
	**   make space for data base entry, load pointer to var
	*/

					Mcheckbase[Mnreads] = (READCHECK *) umalloc(sizeof(READCHECK));
					Mcheckbase[Mnreads]->var = var;



					if (errno || (count < 0)) {
						(void)snprintf (err_buf, 512, "Decoding %s at line number %ld in\n%s\n%s",
							countstr, nline, fd->name, line);
						return (err_buf);
					}

					Mcheckbase[Mnreads]->count = count;

	/* 
	**   allocate enough room to read variables in, depending on variable type
	*/

					if (Mcheckbase[Mnreads]->var) {
						switch (Mcheckbase[Mnreads]->var->type) {
							case M_LONG :
								Mcheckbase[Mnreads]->Types.valuel = (long *)umalloc(count * 
									sizeof(long));
								break;
	     
							case M_FLOAT :
								Mcheckbase[Mnreads]->Types.valuef = (float *)umalloc(count *
									sizeof(float));
								break;
	     
							case M_DOUBLE :
								Mcheckbase[Mnreads]->Types.valued=(double *)umalloc(count *
									sizeof(double));
								break;

						}
					}           
					Mnreads++;
				}
			}
      }
   }

/*
**   Read first line of data
*/
   if (!(fgets (fd->line, max_data_ln_len, fd->fp))) {
      (void)snprintf (err_buf, 512,
         "read_datainfo: Data for first timestep not found in file %s\n",
         fd->name);
      return (err_buf);
   }

   return (NULL);
}
