/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : read_line
 *            Returns 1l if end of data, 0l if more data to be read,
 *            and 2l if end of file
 * COMMENT  : Reads one line from data file into a string,
 *            decodes date and time, and returns if the time
 *            is within start and end limits.
 *            Otherwise reads lines until within limits or the
 *            end of file is encountered.
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define READ_LINE_C

#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**4***************** DECLARATION LOCAL FUNCTIONS *********************/
static void INSERT_time (char *, DATETIME *);

/**5*********************** LOCAL VARIABLES ***************************/
/*
static double   prevjt = -1.0;
*/

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : read_line
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : void
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long read_line (void) {

   /*static char err[80];*/

   char   *start_point, *end_point;
   //float   initial_deltat;
   long   i,j;
   static int   start_of_data;
   static long	data_eof_flag;
   DATETIME   prevtime;
   FILE_DATA   *cur_fd;
   char   *err_ptr;
   static char *line = NULL;
   static long silent_flag;
   silent_flag = *control_lvar("print_debug");

   if (line == NULL) {
	   line = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

/*
**   get initial delta-t from control data base
*/
 //  initial_deltat = *control_fvar("initial_deltat");
   data_eof_flag = *control_lvar ("ignore_data_file_end");

   if (Mnsteps == 0) {
      start_of_data = TRUE;
      Mprevjt = 1.0;
/*
      prevjt = -1.0;
*/
   }

   prevtime = *Mnowtime;

   while (TRUE) {

/*
**  Load cur_fd with the data for the next time step.
*/
      cur_fd = FILE_with_next_ts ();

/*
**  9999 in the year field is the code for EOF. 
*/
      if (cur_fd->time.year == 9999) {
		  if (silent_flag > -2) {
			  (void)fprintf(stderr, "\nWARNING, date of end_time reached the last date in the Data File \n");
			  (void)fprintf(stderr, "         simulation stopped on: %ld %ld %ld \n", Mnowtime->year, Mnowtime->month, Mnowtime->day);
		  }
	      return ENDOFFILE;
	  }

/*
**   DANGER -- This "if" is a hack to get delta time back after storm mode
**            This is the situation :

cur_fd->time = {year = 1956, month = 2, day = 19, hour = 0, min = 0, sec = 0, 
  jd = 2435523, jt = 2435523}

*Mnowtime = {year = 1956, month = 2, day = 18, hour = 24, min = 0, sec = 0, 
  jd = 2435522, jt = 2435523}

*/

      if ((Mprevjt < 0.0) && (cur_fd->time.jt - Mprevjt <= 0.000001)) {
         Mprevjt = (double)(Mnowtime->jd);
      }

/*
**  End of DANGER
*/

/*
**   Copy time from current file into global time structure
*/
      Mnowtime->year = cur_fd->time.year;
      Mnowtime->month = cur_fd->time.month;
      Mnowtime->day = cur_fd->time.day;
      Mnowtime->hour = cur_fd->time.hour;
      Mnowtime->min = cur_fd->time.min;
      Mnowtime->sec = cur_fd->time.sec;
      Mnowtime->jd = cur_fd->time.jd;
      Mnowtime->jt = cur_fd->time.jt;

/*
**   check if data time is within limits
*/
      if (Mnowtime->jt > Mendtime->jt) {
		  /*
		  (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->year);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->year);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->month);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->month);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->day);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->day);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->hour);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->hour);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->min);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->min);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->sec);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->sec);
		 (void)fprintf (stderr,"\n\n nowtime = %ld\n", Mnowtime->jd);
		 (void)fprintf (stderr,"\n\n endtime = %ld\n", Mendtime->jd);
  		 (void)fprintf (stderr,"\n\n nowtime = %f\n", Mnowtime->jt);
		 (void)fprintf (stderr,"\n\n endtime = %f\n", Mendtime->jt); 
		 */
		  *Mnowtime = prevtime;
		 //(void)fprintf (stderr,"nowtime=endtime\n");
		 if (data_eof_flag == 1) {
			Mendtime = Mnowtime;
			return (0); }
		 //(void)fprintf (stderr,"nowtime=endtime 2\n");
         return ENDOFDATA;
      }

      if (Mnowtime->jt >= Mstrttime->jt) {
         if (start_of_data) {
            start_of_data = 0;
            //Mprevjt = Mnowtime->jt - (double)(initial_deltat / 24.0);
			Mprevjt = Mnowtime->jt - (double)(1.0);
         }

         (void)strncpy (line, cur_fd->start_of_data, max_data_ln_len);
         Mnsteps++;
/*
**  DANGER -- Mprevjt must be hacked if starting from var init file.
**            It is computed based on current time and Mdeltat (which
**            is read from var init file.
*/

/*
         (void)fprintf (stderr,"\n\n read_line Mprevjt = %f Mnowtime->jt = %f\n", Mprevjt, Mnowtime->jt);
         (void)fprintf (stderr,"     read_line dt = %f\n", Mnowtime->jt - Mprevjt );
         (void)fprintf (stderr,"          year, mon, day = %d %d %d %d %d\n", Mnowtime->year, Mnowtime->month, Mnowtime->day, Mnowtime->hour, Mnowtime->min);
*/

		 if (Mnowtime->jt < Mprevjt) {
			(void)fprintf (stderr,"\n\n read_line Mprevjt = %f Mnowtime->jt = %f\n", Mprevjt, Mnowtime->jt);
			(void)fprintf (stderr,"Current time step is before previous time step.\n");
			(void)fprintf (stderr,"The data file(s) are running backwards.\n");
			return (ERROR_TIME);
		 }

         if ((Mnowtime->jt - Mprevjt) < 0.0000115) {
/*
** DANGER This hack is to come out of the storm
*/
            (void)fprintf (stderr,"read_line:  coming out of storm. dt = 1 day\n");
            Mdeltat = 1.0;
            Mprevjt = Mnowtime->jt - Mdeltat;

         } else {
            if (Mprevjt < 0.0) {
               Mprevjt = Mnowtime->jt - Mdeltat;
            } else {
               Mdeltat = Mnowtime->jt - Mprevjt;
            }
         }
/*
         (void)fprintf (stderr,"read_line Mdeltat = %f\n", Mdeltat );
*/

/*
**  End DANGER
*/

         Minpptr = end_point = line;

/* 
**   Read variables from the line into their respective buffers
*/
         for (i = 0; i < Mnreads; i++) {
            for (j = 0; j < Mcheckbase[i]->count; j++) {
               start_point = NULL;
               if (Mcheckbase[i]->var) {
                  start_point = end_point;
                  errno = 0;
                  switch (Mcheckbase[i]->var->type) {
                     case M_LONG:
                        Mcheckbase[i]->Types.valuel[j] =
                           strtol (start_point, &end_point, 10);
                        break;
      
                     case M_FLOAT:
                        Mcheckbase[i]->Types.valuef[j] =
                           (float)strtod (start_point, &end_point);
                        break;
      
                     case M_DOUBLE:
                        Mcheckbase[i]->Types.valued[j] =
                           strtod (start_point, &end_point);
                        break;
                     }
				  	if (CHECK_data (errno, cur_fd))	 (void)fprintf (stderr,"nowtime=endtime 2\n");
                  if (CHECK_data (errno, cur_fd)) return (ENDOFDATA);

               } else {
                  (void)strtod (start_point, &end_point);
               }
            }
         }

         Mprevjt = Mnowtime->jt;

         if (!(fgets (cur_fd->line, max_data_ln_len, cur_fd->fp))) {
            fclose (cur_fd->fp);
            cur_fd->fp = NULL;
            cur_fd->time.year = 9999;
         } else if (cur_fd->line[0] == '\n') {
            fclose (cur_fd->fp);
            cur_fd->fp = NULL;
            cur_fd->time.year = 9999;
         } else {
            err_ptr = EXTRACT_time (cur_fd);
            if (err_ptr) {
               (void)fprintf (stderr,"%s\n", err_ptr);
               return (ERROR_TIME);
            }
         }

/*
**   Copy time from current file into global next time structure
*/
         if (cur_fd && cur_fd->fp) {
            Mnexttime->year = cur_fd->time.year;
            Mnexttime->month = cur_fd->time.month;
            Mnexttime->day = cur_fd->time.day;
            Mnexttime->hour = cur_fd->time.hour;
            Mnexttime->min = cur_fd->time.min;
            Mnexttime->sec = cur_fd->time.sec;
            Mnexttime->jd = cur_fd->time.jd;
            Mnexttime->jt = cur_fd->time.jt;
            Mdeltanext = Mnexttime->jt - Mnowtime->jt;
         }
         return (0);
      } else {
/*
**   Read throgh the data before the start time.
*/
         Mprevjt = Mnowtime->jt;

         if (!(fgets (cur_fd->line, max_data_ln_len, cur_fd->fp))) {
            fclose (cur_fd->fp);
            cur_fd->fp = NULL;
            cur_fd->time.year = 9999;
         }

         err_ptr = EXTRACT_time (cur_fd);
         if (err_ptr) {
               (void)fprintf (stderr,"%s\n", err_ptr);
            return (ERROR_TIME);
         }
      }
   }
/*
   return (0);
*/
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : DATA_read_init
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *DATA_read_init (void) {

   static char err[256];

   int      i;
   static int      num_data_files = 0;
   char   **fname, *err_ptr;
   static char      buf[256];
   static char   *line = NULL;

   if (line == NULL) {
	   line = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

/*
**   Clean up the old files
*/
   if (fd) {
      for (i = 0; i < num_data_files; i++)
         if ((fd[i])->fp) {
               fclose ((fd[i])->fp);
            fd[i]->fp = NULL;
            }
   }

   fname =   control_svar ("data_file");
   num_data_files = control_var_size ("data_file");

   fd = (FILE_DATA **)malloc (num_data_files * sizeof (FILE_DATA *));
    for (i = 0; i < num_data_files; i++) {
      fd[i] = (FILE_DATA *)malloc (sizeof (FILE_DATA));
	  fd[i]->line = (char *) umalloc(max_data_ln_len * sizeof(char));
	  fd[i]->info = (char *) umalloc(max_data_ln_len * sizeof(char));
    }

/*
**   Open the files.
*/
   for (i = 0; i < num_data_files; i++) {
      (fd[i])->name = strdup (fname[i]);
      if (!((fd[i])->fp = fopen (fname[i], "r"))) {
         (void)snprintf (err, 256, "DATA_read_init: can't open data file %s\n",
            fname[i]);
         return (err);
      }

      fgets (line, max_data_ln_len, (fd[i])->fp);
      while (strncmp (line, "####", 4)) {
         if (!(fgets (line, max_data_ln_len, (fd[i])->fp))) {
               (void)snprintf (buf, 256, "DATA_read_init - Spacing fwd to data - Check format of file %s.", fname[i]);
               return (buf);
            }
      }

    /*
     * initialize year
     * PJR 7/10/95
     */
       (fd[i])->time.year = 0;
      fgets ((fd[i])->line, max_data_ln_len, (fd[i])->fp);
/* DANGER
      if (err_ptr = EXTRACT_time (&(fd[i])))
*/
      err_ptr = EXTRACT_time (fd[i]);
      if (err_ptr) return (err_ptr);
    
   }
   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : READ_data_info
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : char *
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *READ_data_info (void) {

   static char err[256];

   int      i, num_data_files;
   char   **fname, *err_ptr;
   FILE_DATA  lfd;
   struct stat stbuf;

   lfd.info = (char *) umalloc(max_data_ln_len * sizeof(char));
   lfd.line = (char *) umalloc(max_data_ln_len * sizeof(char));

   fname =   (char **) control_var ("data_file");
   num_data_files = control_var_size ("data_file");

/*
**   Check the files
*/
   for (i = 0; i < num_data_files; i++) {
      if (stat (fname[i], &stbuf) == -1) {
         (void)snprintf (err, 256, "Reading Data Info: Can't open data file %s\n",
            fname[i]);
         return (err);
      } else if ((stbuf.st_mode & S_IFMT) == S_IFDIR) {
         (void)snprintf (err, 256, "Reading Data Info: Can't open data file %s\n",
            fname[i]);
         return (err);
      }
   }
/*
**   Open the files.
*/
   for (i = 0; i < num_data_files; i++) {
      lfd.name = strdup (fname[i]);
      if (!(lfd.fp = fopen (fname[i], "r"))) {
         (void)snprintf (err, 256, "DATA_read_init: can't open data file %s\n",
            fname[i]);
         return (err);
      }
      err_ptr = read_datainfo (&lfd);
      if (err_ptr) return (err_ptr);

      fclose (lfd.fp);
   }
   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : DATA_check_start
 | COMMENT      : Check if start time of model is more than a day before
 |                 the start time of the data.
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char *DATA_check_start (void) {

   int num_data_files, found, i;

   num_data_files = control_var_size ("data_file");

   found = FALSE;
   for (i = 0; (i < num_data_files) && !found; i++)
      if ((fd[i])->time.jd <= Mstrttime->jd)
         found = TRUE;

   if (!found)
      return ("Start time is before first data.");

   return (NULL);
}
/*--------------------------------------------------------------------*\
 | FUNCTION     : DATA_close
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : void 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void DATA_close (void) {
   int      i, num_data_files;

   num_data_files = control_var_size ("data_file");

   for (i = 0; i < num_data_files; i++) {
        if (((fd[i])->fp) != NULL) {
         fclose ((fd[i])->fp);
         (fd[i])->fp = NULL;
        }
    }

   fd = NULL;
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : control_var_size
 | COMMENT      : returns the size of the array
 | PARAMETERS   :
 | RETURN VALUE : int - returns the size of the array
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int control_var_size (char *key) {
   CONTROL *control;

   if (!(control = control_addr(key))) {
      (void)fprintf (stderr, 
         "control_var_size - key '%s' not found.\n", key);
      return (1);
   }
   return (control->size);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : FILE_with_next_ts
 | COMMENT      : Determine the file with the next time step.
 | PARAMETERS   : None
 | RETURN VALUE : Pointer to file data structure
 | RESTRICTIONS : None
\*--------------------------------------------------------------------*/
FILE_DATA * FILE_with_next_ts (void) {

   FILE_DATA *cur_fd, *fd_ptr;
   int      num_data_files, i;
   char   *err_ptr;

   num_data_files = control_var_size ("data_file");
   cur_fd = fd[0];

   for (i = 1; i < num_data_files; i++) {
        fd_ptr = fd[i];  
      if (fd_ptr->time.year != 9999) {
/*
**   If two files have the same julian day, assume one is a storm file
**   and one is a daily file.  Throw out the daily value.
*/
         if (fd_ptr->time.jd == cur_fd->time.jd) {
            if (fd_ptr->time.jt == cur_fd->time.jt) {
               (void)fprintf (stderr,
                  "FILE_with_next_ts: The files %s and %s both seem to contain the same storm on %ld - %ld - %ld.\n",
                  fd_ptr->name, cur_fd->name, fd_ptr->time.year,
                  fd_ptr->time.month, fd_ptr->time.day);

            } else if (fd_ptr->time.jt < cur_fd->time.jt) {

               Mprevjt = fd_ptr->time.jt;

               if (!(fgets (fd_ptr->line, max_data_ln_len, fd_ptr->fp))) {
                  fclose (fd_ptr->fp);
                      fd_ptr->fp = NULL;
                  fd_ptr->time.year = 9999;
               } else if (fd_ptr->line[0] == '\n') {
                  fclose (fd_ptr->fp);
                      fd_ptr->fp = NULL;
                  fd_ptr->time.year = 9999;
               } else {
                  err_ptr = EXTRACT_time (fd_ptr);
                  if (err_ptr) (void)fprintf (stderr,"%s\n", err_ptr);
               }


            } else {

               Mprevjt = cur_fd->time.jt;

               if (!(fgets (cur_fd->line, max_data_ln_len, cur_fd->fp))) {
                  fclose (cur_fd->fp);
                      cur_fd->fp = NULL;
                  cur_fd->time.year = 9999;
               } else if (cur_fd->line[0] == '\n') {
                  fclose (cur_fd->fp);
                      cur_fd->fp = NULL;
                  cur_fd->time.year = 9999;
               } else {
                  err_ptr = EXTRACT_time (cur_fd);
                  if (err_ptr) (void)fprintf (stderr,"%s\n", err_ptr);
               }

               cur_fd = fd_ptr;
            }

         } else if (fd_ptr->time.jd < cur_fd->time.jd) {
            cur_fd = fd_ptr;
         }
      }
   }

   return (cur_fd);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : EXTRACT_time
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
char * EXTRACT_time (FILE_DATA *data) {
   char   *start_point, *end_point;
   static char   *line = NULL;
   static char   *err_buf = NULL;

   if (line == NULL) {
	   line = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

   if (err_buf == NULL) {
	   err_buf = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

   if (data->time.year == 9999) {
	   (void)fprintf (stderr,"9990\n");
      return (NULL);
   }

   (void)strncpy (line, data->line, max_data_ln_len);

   start_point = line;

   errno = 0;
   data->time.year = strtol (start_point, &end_point, 10);

   if (data->time.year < 1800 || data->time.year > 2200) {
      (void)snprintf (err_buf, max_data_ln_len, "EXTRACT_time - year %ld out of range.\nline:%s",
         data->time.year, data->line);
      return (err_buf);
   }

   if (errno == EDOM || errno == ERANGE) {
      (void)snprintf(err_buf, max_data_ln_len, "EXTRACT_time - Decoding year line %s", start_point);
      return (err_buf);
   }

   start_point = end_point;
   errno = 0;
   data->time.month = strtol (start_point, &end_point, 10);

   if (data->time.month < 1 || data->time.month > 12) {
      (void)snprintf (err_buf, max_data_ln_len, "EXTRACT_time - month %ld out of range.\nline:%s",
         data->time.month, data->line);
      return (err_buf);
   }

   if (errno == EDOM || errno == ERANGE) {
      (void)snprintf(err_buf, max_data_ln_len, "EXTRACT_time - Decoding month line %s",start_point);
      return (err_buf);
   }

   start_point = end_point;
   errno = 0;
   data->time.day = strtol (start_point, &end_point, 10);

   if (data->time.day < 1 || data->time.day > 31) {
      (void)snprintf (err_buf, max_data_ln_len, "EXTRACT_time - day %ld out of range.\nline:%s",data->time.day, data->line);
      return (err_buf);
   }

   if (errno == EDOM || errno == ERANGE) {
      (void)snprintf (err_buf, max_data_ln_len, "Decoding day line %s", start_point);
      return (err_buf);
   }

   start_point = end_point;
   errno = 0;
   data->time.hour = strtol(start_point, &end_point, 10);

   if (data->time.hour < 0 || data->time.hour > 24) {
      (void)snprintf(err_buf, max_data_ln_len, "EXTRACT_time - hour %ld out of range.\nline:%s",data->time.hour, data->line);
      return  (err_buf);
   }

   if (errno == EDOM || errno == ERANGE)   {
      (void)snprintf (err_buf, max_data_ln_len, "Decoding hour line %s", start_point);
      return (err_buf);
   }

   start_point = end_point;
   errno = 0;
   data->time.min = strtol(start_point, &end_point, 10);


    if (data->time.min < 0 || data->time.min > 59) {
      (void)snprintf (err_buf, max_data_ln_len, "EXTRACT_time - minute %ld out of range.\nline:%s",
         data->time.min, data->line);
      return (err_buf);
   }

   if (errno == EDOM || errno == ERANGE) {
      (void)snprintf (err_buf, max_data_ln_len, "Decoding minute line %s", start_point);
      return (err_buf);
   }

   start_point = end_point;
   errno = 0;
   data->time.sec = (int)(strtod(start_point, &end_point));
/*
   data->time.sec = strtod(start_point, &end_point, 10);
*/

   if (data->time.sec < 0 || data->time.min > 59) {
      (void)snprintf (err_buf, max_data_ln_len, "EXTRACT_time - second %ld out of range.\nline:%s",
         data->time.sec, data->line);
      return (err_buf);
   }

   if (errno == EDOM || errno == ERANGE) {
      (void)snprintf (err_buf, max_data_ln_len, "Decoding second line %s\n", start_point);
      return (err_buf);
   }

   data->start_of_data = data->line + (end_point - line);

   julday (&(data->time));

   return (NULL);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : CHECK_data
 | COMMENT      : See if there is an error with the data value.
 | PARAMETERS   :
 | RETURN VALUE : Error code
 | RESTRICTIONS : None
\*--------------------------------------------------------------------*/
int CHECK_data (int en, FILE_DATA *cur_fd) {
   if (en == EDOM || en == ERANGE) {
      (void)fprintf (stderr,"read_line");
      perror (" ");
      (void)fprintf (stderr, "Reading line %s\n", cur_fd->line);
      return (ENDOFDATA);
   }
   return (0);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : DATA_find_end
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void DATA_find_end (DATETIME *start_of_data, DATETIME *end_of_data) {
  FILE   *f_ptr;
  int      i, num_data_files;
  DATETIME check;
  static char *line = NULL;

   if (line == NULL) {
	   line = (char *) umalloc(max_data_ln_len * sizeof(char));
   }

  num_data_files = control_var_size ("data_file");

  /*
  **  Get start and end of first file.
  */
  f_ptr = fopen ((fd[0])->name, "r");
  
  fgets (line, max_data_ln_len, f_ptr);
  while (strncmp (line, "####", 4))
    fgets (line, max_data_ln_len, f_ptr);

  fgets (line, max_data_ln_len, f_ptr);
  INSERT_time (line, start_of_data);

  while (fgets (line, max_data_ln_len, f_ptr));

  INSERT_time (line, end_of_data);

  fclose (f_ptr);

  /*
  **  Loop through the other ones.
  */
  for (i = 1; i < num_data_files; i++) {
    f_ptr = fopen ((fd[i])->name, "r");

    fgets (line, max_data_ln_len, f_ptr);
    while (strncmp (line, "####", 4))
      fgets (line, max_data_ln_len, f_ptr);

    fgets (line, max_data_ln_len, f_ptr);
    INSERT_time (line, &check);

    if (check.jt < start_of_data->jt) {
      start_of_data->year = check.year;
      start_of_data->month = check.month;
      start_of_data->day = check.day;
      start_of_data->hour = check.hour;
      start_of_data->min = check.min;
      start_of_data->sec = check.sec;
      start_of_data->jd = check.jd;
      start_of_data->jt = check.jt;
    }

    while (fgets (line, max_data_ln_len, f_ptr));
    INSERT_time (line, &check);

    if (check.jt > end_of_data->jt) {
      end_of_data->year = check.year;
      end_of_data->month = check.month;
      end_of_data->day = check.day;
      end_of_data->hour = check.hour;
      end_of_data->min = check.min;
      end_of_data->sec = check.sec;
      end_of_data->jd = check.jd;
      end_of_data->jt = check.jt;
    }
    fclose (f_ptr);
  }
}

/**7****************** LOCAL FUNCTION DEFINITIONS *********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : INSERT_time
 | COMMENT      :
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
static void INSERT_time (char *line, DATETIME *ptr) {
   char   *start_point, *end_point;

   start_point = line;
   ptr->year = strtol (start_point, &end_point, 10);

   start_point = end_point;
   ptr->month = strtol (start_point, &end_point, 10);

   start_point = end_point;
   ptr->day = strtol (start_point, &end_point, 10);

   start_point = end_point;
   ptr->hour = strtol(start_point, &end_point, 10);

   start_point = end_point;
   ptr->min = strtol(start_point, &end_point, 10);

   start_point = end_point;
   ptr->sec = strtol(start_point, &end_point, 10);

   julday (ptr);
}
