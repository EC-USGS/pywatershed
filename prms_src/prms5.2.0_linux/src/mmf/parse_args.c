/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : parse_args
 * COMMENT  : parses the command line arguments
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define PARSE_ARGS_C
#include <math.h> 
#include <string.h> 
#include <stdlib.h> 
#include "mms.h" 

/*--------------------------------------------------------------------*\
 | FUNCTION     : parse_args
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : void
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void parse_args (int argc, char **argv, int *set_count, char **set_name, char **set_value, int set_size) {

   int i;
   char *ptr;

   Mdebuglevel = 0;
   MAltContFile = strdup ("control");

/*
**  Get the model name.
*/
   ptr = strrchr (argv[0], '/');
   if (!ptr) ptr = strrchr (argv[0], '\\');
   if (ptr) ++ptr;
   else ptr = argv[0];

   model_name = strdup (ptr);

   executable_model = strdup (argv[0]);
   ptr = strstr (executable_model, ".exe");
   if (ptr) *ptr = '\0';

   if (argc >= 2) {
      for (i = 1; i < argc ; i++) {
		 if (!strcmp(argv[i], "-debug")) {
			 Mdebuglevel = atoi(argv[i+1]);
			 i++;

		 } else if (!strncmp(argv[i],"-C",2)) {
            MAltContFile = (char *)((argv[i]));
            MAltContFile+=2;

         } else if (!strncmp(argv[i],"-batch", 6)){
            batch_run_mode = TRUE;

         } else if (!strncmp(argv[i],"-print", 6)){
            print_mode = TRUE;

         } else if (!strncmp(argv[i],"-por", 4)){
            run_period_of_record = TRUE;

         } else if (!strncmp(argv[i],"-rtg", 4)){
            runtime_graph_on = TRUE;

		 } else if (!strncmp(argv[i],"-preprocess", 11)){
            preprocess_on = TRUE;

         } else if (!strncmp(argv[i],"-set",4)){

            if ((*set_count) >= set_size) {
               printf("parse_args: Overflow. Too many command line arguments set with -set flag.\n\n\n");
               exit(1);
            }

            i++;
            *(set_name + *set_count) = strdup ((char *)((argv[i])));
            i++;
            *(set_value + *set_count) = strdup ((char *)((argv[i])));
            (*set_count)++;

		} else if (!strncmp(argv[i],"-MAXDATALNLEN",13)){
            max_data_ln_len = atoi(argv[i+1]);
			i++;

		 } else { // Assume argument with no flag is control file name
			MAltContFile = (char *)((argv[i]));
		 }
      }
   }
}
