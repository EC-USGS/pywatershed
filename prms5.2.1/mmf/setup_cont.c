/*+
 * United States Geological Survey
 *
 * PROJECT  : Modular Modeling System (MMS)
 * FUNCTION : setup_cont
 * COMMENT  :
 *
 * $Id$
 *
-*/

/**1************************ INCLUDE FILES ****************************/
#define SETUP_CONT_C

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include "mms.h"

/**5*********************** LOCAL VARIABLES ***************************/
extern void decl_control_string (char *key, char *valstr);
extern void decl_control_int_array (char *key, long size, long *valstr);
extern void decl_control_float_array (char *key, long size, float *valstr);
extern void decl_control_string_array (char *key, long size, char *valstr);

/**6**************** EXPORTED FUNCTION DEFINITIONS ********************/
/*--------------------------------------------------------------------*\
 | FUNCTION     : setup_cont
 | COMMENT		:
 | PARAMETERS   :
 | RETURN VALUE : void
 | RESTRICTIONS :
\*--------------------------------------------------------------------*/
void setup_cont (void) {
        long *lval;
        float *fval;

        static long start_date[] = {2000,10,1,0,0,0};
        static long end_date[] = {2001,9,30,0,0,0};

/*
**	GSFLOW control variables
*/
        decl_control_string ("model_mode", "PRMS5");
        decl_control_string ("precip_module", "precip_1sta");
        decl_control_string ("temp_module", "temp_1sta");
        decl_control_string ("et_module", "potet_jh");
        decl_control_string ("srunoff_module", "srunoff_smidx");
        decl_control_string ("solrad_module", "ddsolrad");
		decl_control_string ("soilzone_module", "soilzone");
		decl_control_string ("capillary_module", "soilzone");
		decl_control_string ("strmflow_module", "strmflow");
        decl_control_string ("transp_module", "transp_tindex");
        decl_control_string ("gsflow_output_file", "gsflow.out");
        decl_control_string ("gsflow_csv_file", "gsflow.csv");
		//decl_control_string ("creator_email", "unknown");

/*
        cval = (char *)umalloc (sizeof (long));
        cval[0] = "recharge";
        decl_control_string_array ("mapOutVar_names", 20, cval);
*/

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 7;
        decl_control_int_array ("rpt_days", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
        decl_control_int_array ("gsf_rpt", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("print_debug", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
		decl_control_int_array ("cascade_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
		decl_control_int_array ("cascadegw_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
		decl_control_int_array ("subbasin_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("gwr_swale_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("frozen_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dprst_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
		decl_control_int_array ("parameter_check_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 1;
		decl_control_int_array ("cbh_check_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("cbh_binary_flag", 1, lval);		

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_imperv_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_intcp_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_covden_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_sro2dprst_perv_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_sro2dprst_imperv_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_covtype_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_transp_flag", 1, lval);

		lval = (long *)umalloc(sizeof (long));
		lval[0] = 0;
		decl_control_int_array("dyn_transp_on_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_fallfrost_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_springfrost_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_potet_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_soil_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_radtrncf_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array("dyn_snareathresh_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_sro_to_dprst_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_sro_to_imperv_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dyn_dprst_flag", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("dyn_ag_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("stream_temp_flag", 1, lval);

		lval = (long *)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("strmtemp_humidity_flag", 1, lval);
		
		lval = (long *)umalloc(sizeof (long));
		lval[0] = 0;
		decl_control_int_array("stream_temp_shade_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("segment_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("gwr_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("external_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("consumed_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("lake_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dprst_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dprst_transfer_water_use", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("dprst_add_water_use", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("soilzone_transferON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("canopy_transferON_OFF", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("outputSelectDatesON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("seg2hru_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("glacier_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("mbInit_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("musroute_flag", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("orad_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("snow_cbh_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("gwflow_cbh_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("humidity_cbh_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("windspeed_cbh_flag", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("albedo_cbh_flag", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("cloud_cover_cbh_flag", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("segmentOutON_OFF", 1, lval);

		lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("ignore_data_file_end", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("snarea_curve_flag", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("soilzone_aet_flag", 1, lval);

		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("snow_cloudcover_flag", 1, lval);

/*
**	file names
*/
        decl_control_string ("executable_desc", "MOWS executable");
        decl_control_string ("executable_model", "prmsIV");
        decl_control_string ("data_file", "prms.data");
        decl_control_string ("param_file", "prms.params");
        decl_control_string ("var_save_file", "prms_ic.out");
        decl_control_string ("var_init_file", "prms_ic.in");
        decl_control_string ("stat_var_file", "statvar.out");
        decl_control_string ("ani_output_file", "animation.out");
        decl_control_string ("model_output_file", "prms.out");
        decl_control_string ("stats_output_file", "stats.out");
		decl_control_string ("tmax_day", "tmax.day");
        decl_control_string ("tmin_day", "tmin.day");
        decl_control_string ("precip_day", "precip.day");
        decl_control_string ("swrad_day", "swrad.day");
        decl_control_string ("potet_day", "potet.day");
        decl_control_string ("transp_day", "transp.day");
        decl_control_string ("windspeed_day", "windspeed.day");
        decl_control_string ("humidity_day", "humidity.day");
		decl_control_string ("albedo_day", "albedo.day");
		decl_control_string ("cloud_cover_day", "cloudcover.day");
		decl_control_string ("pkwater_equiv_day", "pkwater_equiv.day");
        decl_control_string ("pk_depth_day", "pk_depth.day");
        decl_control_string ("snow_evap_day", "snow_evap.day");
        decl_control_string ("snowcov_area_day", "snowcov_area.day");
        decl_control_string ("snowmelt_day", "snowmelt.day");
        decl_control_string ("gwres_flow_day", "gwres_flow.day");
        decl_control_string ("dprst_area_dynamic", "dyndprst_area");
        decl_control_string ("dprst_depth_dynamic", "dyndprst_depth");
        decl_control_string ("dprst_frac_dynamic", "dyndprst_frac");
		decl_control_string ("snow_intcp_dynamic", "dynsnowintcp");
		decl_control_string ("srain_intcp_dynamic", "dynsrainintcp");
		decl_control_string ("wrain_intcp_dynamic", "dynwrainintcp");
		decl_control_string ("imperv_frac_dynamic", "dynimperv");
		decl_control_string ("imperv_stor_dynamic", "dynimperv");
		decl_control_string ("covtype_dynamic", "dyncovtype");
		decl_control_string ("covden_sum_dynamic", "dyncovden_sum");
		decl_control_string ("covden_win_dynamic", "dyncovden_win");
		decl_control_string ("jhcoef_dynamic", "dynjhcoef");
		decl_control_string ("potet_coef_dynamic", "dynpotetcoef");
		decl_control_string ("transpbeg_dynamic", "dyntranspbeg");
		decl_control_string ("transpend_dynamic", "dyntranspend");
		decl_control_string ("fallfrost_dynamic", "dynfallfrost");
		decl_control_string ("springfrost_dynamic", "dynspringfrost");
		decl_control_string ("soilrechr_dynamic", "dynsoilrechr");
		decl_control_string ("soilmoist_dynamic", "dynsoilmoist");
		decl_control_string ("radtrncf_dynamic", "dynradtrncf");
		decl_control_string ("sro2dprst_perv_dynamic", "dynsro2dprst_perv");
		decl_control_string ("sro2dprst_imperv_dynamic", "dynsro2dprst_imperv");
		decl_control_string ("transp_on_dynamic", "dyntranspon");
		decl_control_string ("csv_output_file", "prms_summary.csv");
        decl_control_string ("nhruOutBaseFileName", "nhruout_path");
		decl_control_string ("nsubOutBaseFileName", "nsubout_path");
		decl_control_string ("basinOutBaseFileName", "basinout_path");
		decl_control_string("nsegmentOutBaseFileName", "nsegmentout_path");
		decl_control_string("dynamic_param_log_file", "dynamic_parameter.out");
		decl_control_string("selectDatesFileName", "selectDates.in");
		decl_control_string("precip_map_file", "precip_map.dat");
		decl_control_string("temp_map_file", "temp_map.dat");
/*
**	run start and end times
*/
        decl_control_int_array("start_time", 6, start_date);
        decl_control_int_array("end_time", 6, end_date);

/*
**	flag for initializing vars from file
*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("init_vars_from_file", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("save_vars_to_file", 1, lval);

/*
**	initial delta-t - hours
*/
        fval = (float *)umalloc (sizeof (float));
		fval[0] = 24.0;
        decl_control_float_array ("initial_deltat", 1, fval);

/*
**	stats analysis
*/

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("statsON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("nstatVars", 1, lval);

/*
**	animation output
*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("aniOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("naniOutVars", 1, lval);

/*
**	map output
*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("mapOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("nmapOutVars", 1, lval);

/*
**	nhru_summary
*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nhruOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nhruOutVars", 1, lval);
        lval = (long *)umalloc(sizeof (long));
		lval[0] = 1;
		decl_control_int_array("nhruOut_freq", 1, lval);
        lval = (long*)umalloc(sizeof(long));
		lval[0] = 1;
		decl_control_int_array("nhruOut_format", 1, lval);
		lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("nhruOutNcol", 1, lval);

        lval = (long*)umalloc(sizeof(long));
		lval[0] = 0;
		decl_control_int_array("prms_warmup", 1, lval);


		/*
		**	nsub_summary
		*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nsubOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nsubOutVars", 1, lval);
        lval = (long *)umalloc(sizeof (long));
		lval[0] = 1;
		decl_control_int_array("nsubOut_freq", 1, lval);
        lval = (long*)umalloc(sizeof(long));
		lval[0] = 1;
		decl_control_int_array("nsubOut_format", 1, lval);

		/*
		**	basin_summary
		*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("basinOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("basinOutVars", 1, lval);
        lval = (long *)umalloc(sizeof (long));
		lval[0] = 1;
		decl_control_int_array("basinOut_freq", 1, lval);
        lval = (long*)umalloc(sizeof(long));
		lval[0] = 1;
		decl_control_int_array("basinOut_format", 1, lval);
		/*
		**	nsegment_summary
		*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nsegmentOutON_OFF", 1, lval);
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
		decl_control_int_array ("nsegmentOutVars", 1, lval);
        lval = (long *)umalloc(sizeof (long));
		lval[0] = 1;
		decl_control_int_array("nsegmentOut_freq", 1, lval);
        lval = (long*)umalloc(sizeof(long));
		lval[0] = 1;
		decl_control_int_array("nsegmentOut_format", 1, lval);

/*
**	graphics display
*/
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("ndispGraphs", 1, lval);

        lval = (long *)umalloc (sizeof (long));
		lval[0] = 50;
        decl_control_int_array ("dispGraphsBuffSize", 1, lval);

// CSV output
        lval = (long *)umalloc (sizeof (long));
		lval[0] = 0;
        decl_control_int_array ("csvON_OFF", 1, lval);

}
