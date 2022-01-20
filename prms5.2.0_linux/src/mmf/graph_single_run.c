/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : graph_single_run
 * COMMENT  : graph routines for mms run
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define GRAPH_SINGLE_RUN_C
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mms.h"

#define         MAXNUMBEROFGRAPHS               4

/**5*********************** LOCAL VARIABLES ***************************/
long NdispGraphs;
static double zero_time;
PUBVAR **disp_var;
int *disp_ele;
int numDispVars;

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : initializeRuntimeGraphs
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : int
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int initializeRuntimeGraphs (void) {
   CONTROL *control;
   int i;
   //long datetime[6];
   DATETIME starttime_copy;
   char *cptr, *cptr2;

   if (!runtime_graph_on) return (FALSE);

   //dattim("start", datetime);
   //zero_time = getjulday((int)datetime[1],(int)datetime[2],(int)datetime[0],
			//(int)datetime[3], (int)datetime[4],(double)datetime[5]);

      starttime_copy.year =Mstrttime->year;
   starttime_copy.month =Mstrttime->month;
   starttime_copy.day = Mstrttime->day;
   starttime_copy.hour =Mstrttime->hour;
   starttime_copy.min =Mstrttime->min;
   starttime_copy.sec =Mstrttime->sec;

   julday(&starttime_copy);

   //zero_time = zero_time - 1.0;

   zero_time = starttime_copy.jt;

/*
** Get the number of display vars
*/
   cptr = strdup ("dispVar_names");
   control = control_addr(cptr);
   if (control) {
      numDispVars = control->size;

      disp_var = (PUBVAR **)malloc (sizeof(PUBVAR *) * numDispVars);
      memset (disp_var, 0, sizeof (PUBVAR *) * numDispVars);

      disp_ele = (int *)malloc (sizeof (int) * numDispVars);
      memset (disp_ele, 0, sizeof (int) * numDispVars);

      for (i = 0; i < numDispVars; i++) {
/**
** Get address of each display variable for each graph
**/
         cptr = strdup ("dispVar_names");
		 cptr2 = (char *)control_sarray(cptr, i);

         disp_var[i] = var_addr (cptr2);

         cptr = strdup ("dispVar_element");
//         disp_ele[i] =  atoi (*control_sarray(cptr,i)) - 1;
         disp_ele[i] =  atoi (control_sarray(cptr,i)) - 1;
      }
   } else {
	   numDispVars = 0;
	   disp_var = NULL;
	   disp_ele = NULL;
	   runtime_graph_on = 0;
   }
   return (FALSE);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : plotRuntimeGraphValue
 | COMMENT	:
 | PARAMETERS   :
 | RETURN VALUE : int
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int plotRuntimeGraphValue (void) {
   double xval;
   float yval;
   int i;
   //long datetime[6];
   DATETIME nowtime_copy;

   if (!runtime_graph_on) return (FALSE);


   //dattim("now", datetime);
   nowtime_copy.year =Mnowtime->year;
   nowtime_copy.month =Mnowtime->month;
   nowtime_copy.day = Mnowtime->day;
   nowtime_copy.hour =Mnowtime->hour;
   nowtime_copy.min =Mnowtime->min;
   nowtime_copy.sec =Mnowtime->sec;

   julday(&nowtime_copy);
   //xval = getjulday(datetime[1],datetime[2],datetime[0],
	  //               datetime[3], datetime[4],(double)datetime[5]);

   xval = nowtime_copy.jt - zero_time;

   printf ("plotRuntimeGraphValue: xval = %f", xval);

   for (i = 0; i < numDispVars; i++) {
      yval = 0.0;
      switch ((disp_var[i])->type) {
         case M_LONG :
            yval = (float)(*(((long *)((disp_var[i])->value)) + disp_ele[i]));
            break;

         case M_DOUBLE :
            yval = (float)(*(((double *)((disp_var[i])->value)) + disp_ele[i]));
            break;

         case M_FLOAT :
            yval = *(((float *)((disp_var[i])->value)) + disp_ele[i]);
            break;
      }
      printf (" %f", yval);
   }
   printf ("\n");
   fflush (stdout);

   return (FALSE);
}

/*--------------------------------------------------------------------*\
 | FUNCTION     : closeRuntimeGraphs
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : int
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
int closeRuntimeGraphs (void) {
   if (!runtime_graph_on) return (FALSE);

   printf ("closeRuntimeGraph\n");
   return (FALSE);
}
