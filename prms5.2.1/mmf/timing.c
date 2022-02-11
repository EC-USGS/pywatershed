/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : timing
 * COMMENT  : timing functions
 *            The routines with a _ suffix are called from Fortran
 *            The routines without the suffix are called from C
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define TIMING_C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**************************************************************************
 * dattim and dattim_ : get start, end or current data date and time
 *
 * args - dwhen : string, "start", "end", and "now"
 *        timearray: integer or long array which accepts the time
 *                   and date
 * 
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : dattim_
 | COMMENT		: called from Fortran, sorts out args and calls dattim()
 |                 get start, end or current data date and time
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void dattim_ (char *dwhen, ftnint *timearray, ftnlen dwhenlen) {

//  char *when;
  char when[80];
  long ta[6];

  /*
   * copy when and terminate
   */

//  when = (char *) umalloc(dwhenlen + 1);
//  strncpy(when, dwhen, dwhenlen);
//  when[dwhenlen] = '\0';

    strncpy (when, dwhen, dwhenlen);
    *(when + dwhenlen) = '\0';


  /*
   * call C version of dattim()
   */

  dattim(when, ta);
  timearray[0] = ta[0];
  timearray[1] = ta[1];
  timearray[2] = ta[2];
  timearray[3] = ta[3];
  timearray[4] = ta[4];
  timearray[5] = ta[5];

}

/**************************************************************************
 */


/*--------------------------------------------------------------------*\
 | FUNCTION     : dattim
 | COMMENT		: called from C
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void dattim (char *when, long *timearray) {

  DATETIME *time;

  /*
   * set time according to when argument
   */

  if (!strcmp(when, "start"))
    time = Mstrttime;
  else if(!strcmp(when, "end"))
    time = Mendtime;
  else if (!strcmp(when, "now"))
    time = Mnowtime;
  else {
    (void)fprintf(stderr,
	    "ERROR - dattim - illegal argument '%s'.\n", when);
    exit(1);
  }

  /*
   * load up time array
   */

  timearray[0] = time->year;
  timearray[1] = time->month;
  timearray[2] = time->day;
  timearray[3] = time->hour;
  timearray[4] = time->min;
  timearray[5] = time->sec;

}

/*--------------------------------------------------------------------*\
 | FUNCTION     : nowjt_
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double nowjt_ () {
   return Mnowtime->jt;
}

/**************************************************************************
 * julian_ and julian : returns the julian date of the data stream relative
 *                      to calendar, solar and water year start dates.
 *                         (1 JAN)   (22 DEC)  (1 OCT)
 *
 * args - when : string, "start", "end", "now"
 *        type : string, "calendar", "solar", "water", "absolute"
 * 
 * julian_() is called from Fortran, sorts out args and calls julian()
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : julian_
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long julian_ (char *jwhen, char *jtype, ftnlen jwhenlen, ftnlen jtypelen) {

//  char *when, *type;
  char when[80], type[80];
  long retval;

  /*
   * copy strings and terminate
   */

//  when = (char *) umalloc(jwhenlen + 1);
//  strncpy(when, jwhen, jwhenlen);
//  when[jwhenlen] = '\0';
    strncpy (when, jwhen, jwhenlen);
    *(when + jwhenlen) = '\0';

//  type = (char *) umalloc(jtypelen + 1);
//  strncpy(type, jtype, jtypelen);
//  type[jtypelen] = '\0';
    strncpy (type, jtype, jtypelen);
    *(type + jtypelen) = '\0';

  /*
   * call C version of julian()
   */

  retval = julian(when, type);

  return retval;

}

/**************************************************************************
 * julian() is called from C
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : julian
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long julian (char *when, char *type) {

  DATETIME *time, reftime;

  /*
   * set time according to when argument
   */

  if (!strcmp(when, "start"))
    time = Mstrttime;
  else if(!strcmp(when, "end"))
    time = Mendtime;
  else if (!strcmp(when, "now"))
    time = Mnowtime;
  else {
    (void)fprintf(stderr,
	    "ERROR - julian - illegal argument '%s'.\n", when);
    exit(1);
  }

  /*
   * set reftime depending on type arg
   */

  if (!strcmp(type, "calendar")) {
    reftime.year = time->year - 1;
    reftime.month = 12;
    reftime.day = 31;
  } else if(!strcmp(type, "solar")) {
    if ((time->month == 12) && (time->day > 21))
      reftime.year = time->year;
    else
      reftime.year = time->year - 1;
    reftime.month = 12;
    reftime.day = 21;
  } else if(!strcmp(type, "spring")) {
	  if ((time->month > 3) || (time->month == 3 && time->day > 20)) {
		reftime.year = time->year;
	  } else {
		reftime.year = time->year - 1;
	  }
	reftime.month = 3;
	reftime.day = 20;
  } else if (!strcmp(type, "water")) {
    if (time->month > 9)
      reftime.year = time->year;
    else
      reftime.year = time->year - 1;
    reftime.month = 9;
    reftime.day = 30;
  } else if(!strcmp(type, "absolute")) {
    julday(time);
    return (time->jd);
  } else {
    (void)fprintf(stderr,
	    "ERROR - julian - illegal argument '%s'.\n", type);
    exit(1);
  }

  reftime.hour = 0;
  reftime.min = 0;
  reftime.sec = 0;

  /*
   * compute the julian dates
   */

  julday(time);
  julday(&reftime);

  return (time->jd - reftime.jd);

}

/**************************************************************************
 * deltim_() is called from Fortran, deltim()
 */

double deltim_(void) {

/* printf ("from deltim:  %f\n", deltim()); */
  return deltim();

}

/**************************************************************************
 * deltim() is called from C
 */
/*--------------------------------------------------------------------*\
 | FUNCTION     : deltim
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double deltim (void) {
  return (double) Mdeltat * 24.0;
}

/**************************************************************************
 * getstep_() is called from Fortran, getstep()
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : getstep_
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long getstep_ (void) {
  return getstep();
}

/**************************************************************************
 * getstep() is called from C
 */
/*--------------------------------------------------------------------*\
 | FUNCTION     : getstep
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
long getstep (void) {
  return Mnsteps;
}

/**************************************************************************
 * djulian_ and djulian : returns the double julian date of the data stream
 *                      relative to calendar, solar and water year start dates.
 *                         (1 JAN)   (22 DEC)  (1 OCT)
 *
 * args - when : string, "start", "end", "now"
 *        type : string, "calendar", "solar", "water", "absolute"
 * 
 * julian_() is called from Fortran, sorts out args and calls julian()
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : djulian_
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double djulian_ (char *jwhen, char *jtype, ftnlen jwhenlen, ftnlen jtypelen) {

  char *when, *type;
  double retval;

  /*
   * copy strings and terminate
   */

  when = (char *) umalloc(jwhenlen + 1);
  strncpy(when, jwhen, jwhenlen);
  when[jwhenlen] = '\0';

  type = (char *) umalloc(jtypelen + 1);
  strncpy(type, jtype, jtypelen);
  type[jtypelen] = '\0';

  /*
   * call C version of djulian()
   */

  retval = djulian(when, type);

  /*
   * free up arrays
   */

//ufree(when);
//ufree(type);

  return retval;

}

/**************************************************************************
 * julian() is called from C
 */

/*--------------------------------------------------------------------*\
 | FUNCTION     : djulian
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double djulian (char *when, char *type) {

  DATETIME *time, reftime;

  /*
   * set time according to when argument
   */

  if (!strcmp(when, "start"))
    time = Mstrttime;
  else if(!strcmp(when, "end"))
    time = Mendtime;
  else if (!strcmp(when, "now"))
    time = Mnowtime;
  else {
    (void)fprintf(stderr,
	    "ERROR - julian - illegal argument '%s'.\n", when);
    exit(1);
  }

  /*
   * set reftime depending on type arg
   */

  if (!strcmp(type, "calendar")) {
    reftime.year = time->year - 1;
    reftime.month = 12;
    reftime.day = 31;
  } else if(!strcmp(type, "solar")) {
    if ((time->month == 12) && (time->day > 21))
      reftime.year = time->year;
    else
      reftime.year = time->year - 1;
    reftime.month = 12;
    reftime.day = 21;
  } else if (!strcmp(type, "water")) {
    if (time->month > 9)
      reftime.year = time->year;
    else
      reftime.year = time->year - 1;
    reftime.month = 9;
    reftime.day = 30;
  } else if(!strcmp(type, "absolute")) {
    julday(time);
    return (time->jt);
  } else {
    (void)fprintf(stderr,
	    "ERROR - julian - illegal argument '%s'.\n", type);
    exit(1);
  }

  reftime.hour = 0;
  reftime.min = 0;
  reftime.sec = 0;

  /*
   * compute the julian dates
   */

  julday(time);
  julday(&reftime);

  return (time->jt - reftime.jt);

}

/**************************************************************************
 * delnex_() is called from Fortran, delnex()
 */

double delnex_(void) {

/* printf ("from deltim:  %f\n", deltim()); */
  return delnex();

}

/**************************************************************************
 * delnex() is called from C
 */
/*--------------------------------------------------------------------*\
 | FUNCTION     : delnex
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : 
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
double delnex (void) {
  return (double) Mdeltanext * 24.0;
}

