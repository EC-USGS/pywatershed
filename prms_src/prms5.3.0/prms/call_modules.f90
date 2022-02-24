!***********************************************************************
! Defines the computational sequence, valid modules, and dimensions
!***********************************************************************
!***********************************************************************
      SUBROUTINE call_modules(Arg)
      USE PRMS_CONSTANTS, ONLY: ERROR_control
      use PRMS_CONTROL_FILE, only: read_control_file
      use PRMS_DATA_FILE, only: read_prms_data_file
      use PRMS_MMFAPI, only: Num_variables, Variable_data, MAXVARIABLES, declvar_int, declvar_real
      USE PRMS_MODULE
      use PRMS_READ_PARAM_FILE, only: declparam, check_parameters, getparam_int, getparam_real, &
                                      read_parameter_file_dimens, read_parameter_file_params, setup_params
      use PRMS_SET_TIME, only: prms_time
      use prms_utils, only: error_stop, module_error, numchars, print_module, PRMS_open_output_file, read_error
      IMPLICIT NONE
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Arg
! Functions
      INTRINSIC :: DATE_AND_TIME, INT
      INTEGER, EXTERNAL :: check_dims, basin, climateflow, setup
      INTEGER, EXTERNAL :: cascade, obs, soltab, transp_tindex
      INTEGER, EXTERNAL :: transp_frost, frost_date, routing
      INTEGER, EXTERNAL :: temp_1sta_laps, temp_dist2
      INTEGER, EXTERNAL :: precip_1sta_laps, climate_hru
      INTEGER, EXTERNAL :: precip_dist2, xyz_dist, ide_dist
      INTEGER, EXTERNAL :: ddsolrad, ccsolrad
      INTEGER, EXTERNAL :: potet_pan, potet_jh, potet_hamon, potet_hs, potet_pt, potet_pm
      INTEGER, EXTERNAL :: intcp, snowcomp, gwflow
      INTEGER, EXTERNAL :: srunoff, soilzone, soilzone_ag
      INTEGER, EXTERNAL :: strmflow, subbasin, basin_sum, map_results, write_climate_hru
      INTEGER, EXTERNAL :: strmflow_in_out, muskingum, muskingum_lake
      INTEGER, EXTERNAL :: water_use_read, dynamic_param_read, potet_pm_sta
      INTEGER, EXTERNAL :: stream_temp, glacr
      EXTERNAL :: precip_map, temp_map
      EXTERNAL :: call_modules_restart, water_balance, summary_output
      EXTERNAL :: prms_summary, module_doc, convert_params
! Local Variables
      INTEGER :: i, iret, nc, ierr
!***********************************************************************
      ierr = 0

      Process = Arg

      IF ( Process(:3)=='run' ) THEN
        Process_flag = RUN !(0=run, 1=declare, 2=init, 3=clean, 4=setdims)

      ELSEIF ( Process(:4)=='decl' ) THEN
        CALL DATE_AND_TIME(VALUES=Elapsed_time_start)
        Execution_time_start = Elapsed_time_start(5)*3600 + Elapsed_time_start(6)*60 + &
     &                         Elapsed_time_start(7) + Elapsed_time_start(8)*0.001

        Process_flag = DECL

        IF ( check_dims()/=0 ) ERROR STOP ERROR_dim

        IF ( Print_debug>DEBUG_minimum ) THEN
          PRINT 10, PRMS_VERSION
          WRITE ( PRMS_output_unit, 10 ) PRMS_VERSION
        ENDIF
  10  FORMAT (///, 25X, 'U.S. Geological Survey', /, 15X, &
     &        'Precipitation-Runoff Modeling System (PRMS)', /, 24X, A)
  15  FORMAT (/, 8X, 'Process',  12X, 'Available Modules', /, 68('-'), /, &
     &        '  Basin Definition: basin', /, &
     &        '    Cascading Flow: cascade', /, &
     &        '  Time Series Data: obs, water_use_read, dynamic_param_read', /, &
     &        '   Potet Solar Rad: soltab', /, &
     &        '  Temperature Dist: temp_1sta, temp_laps, temp_dist2, climate_hru,', /, &
     &        '                    temp_map', /, &
     &        '       Precip Dist: precip_1sta, precip_laps, precip_dist2,', /, &
     &        '                    climate_hru, precip_map', /, &
     &        'Temp & Precip Dist: xyz_dist, ide_dist', /, &
     &        '    Solar Rad Dist: ccsolrad, ddsolrad, climate_hru', /, &
     &        'Transpiration Dist: transp_tindex, climate_hru, transp_frost', /, &
     &        '      Potential ET: potet_hamon, potet_jh, potet_pan, climate_hru,', /, &
     &        '                    potet_hs, potet_pt, potet_pm, potet_pm_sta', /, &
     &        '      Interception: intcp', /, &
     &        'Snow & Glacr Dynam: snowcomp, glacr_melt', /, &
     &        '    Surface Runoff: srunoff_smidx, srunoff_carea', /, &
     &        '         Soil Zone: soilzone, soilzone_ag', /, &
     &        '       Groundwater: gwflow', /, &
     &        'Streamflow Routing: strmflow, strmflow_in_out, muskingum,', /, &
     &        '                    muskingum_lake, muskingum_mann', /, &
     &        'Stream Temperature: stream_temp', /, &
     &        '    Output Summary: basin_sum, subbasin, map_results, prms_summary,', /, &
     &        '                    nhru_summary, nsub_summary, water_balance', /, &
     &        '                    basin_summary, nsegment_summary, statvar_out', /, &
     &        '     Preprocessing: write_climate_hru, frost_date', /, 68('-'))
  16  FORMAT (//, 4X, 'Active modules listed in the order in which they are called', //, 8X, 'Process', 20X, &
     &        'Module', 9X, 'Version Date', /, A)

        IF ( Print_debug>DEBUG_minimum ) THEN
          PRINT 15
          PRINT 9002
          WRITE ( PRMS_output_unit, 15 )
          PRINT 16, EQULS(:62)
          WRITE ( PRMS_output_unit, 16 ) EQULS(:62)
        ENDIF
        CALL print_module(MODDESC, MODNAME, PRMS_versn)

        CALL print_module('Read Control File', 'read_control_file', Version_read_control_file)
        CALL print_module('Read Parameter File', 'read_parameter_file', Version_read_parameter_file)
        CALL print_module('Read Data File', 'read_data_file', Version_read_data_file)
        CALL read_prms_data_file()

        Num_variables = 0
        ALLOCATE ( Variable_data(MAXVARIABLES) ) ! don't know how many, need to read var_name file

        ALLOCATE ( Hru_type(Nhru) )
        Hru_type = 0
        IF ( declparam(MODNAME, 'hru_type', 'nhru', 'integer', &
     &       '1', '0', '4', &
     &       'HRU type', 'Type of each HRU (0=inactive; 1=land; 2=lake; 3=swale; 4=glacier)', &
     &       'none')/=0 ) CALL read_error(1, 'hru_type')

        Timestep = 0
        IF ( Init_vars_from_file>OFF ) CALL call_modules_restart(READ_INIT)

      ELSEIF ( Process(:4)=='init' ) THEN
        Process_flag = INIT

        IF ( getparam_int(MODNAME, 'hru_type', Nhru, Hru_type)/=0 ) CALL read_error(2, 'hru_type')

        nc = numchars(Model_control_file)
        IF ( Print_debug>DEBUG_less ) PRINT 9004, 'Using Control File: ', Model_control_file(:nc)
        IF ( Print_debug>DEBUG_minimum ) WRITE ( PRMS_output_unit, 9004 ) 'Using Control File: ', Model_control_file(:nc)

        nc = numchars(Param_file)
        IF ( Print_debug>DEBUG_less ) PRINT 9004, 'Using Parameter File: ', Param_file(:nc)
        IF ( Print_debug>DEBUG_minimum ) WRITE ( PRMS_output_unit, 9004 ) 'Using Parameter File: ', Param_file(:nc)

        IF ( Init_vars_from_file>0 ) THEN
          nc = numchars(Var_init_file)
          IF ( Print_debug>DEBUG_less ) PRINT 9004, 'Using var_init_file: ', Var_init_file(:nc)
        ENDIF
        IF ( Save_vars_to_file==ACTIVE ) THEN
          nc = numchars(Var_save_file)
          IF ( Print_debug>DEBUG_less ) PRINT 9004, 'Writing var_save_file: ', Var_save_file(:nc)
        ENDIF

        IF ( Print_debug>DEBUG_minimum ) THEN
          nc = numchars(Model_output_file)
          PRINT 9004, 'Writing PRMS Water Budget File: ', Model_output_file(:nc)
        ENDIF

      ELSEIF ( Process(:7)=='setdims' ) THEN
        Process_flag = SETDIMENS
        Kkiter = 1 ! set for PRMS-only mode
        Canopy_iter = 1
        Soilzone_add_water_use = OFF
        Dprst_add_water_use = OFF
        Dprst_transfer_water_use = OFF
        Gwr_add_water_use = OFF
        Gwr_transfer_water_use = OFF
        Lake_add_water_use = OFF
        Lake_transfer_water_use = OFF
        CALL setdims()
        IF ( ierr/=0 ) CALL module_error('setdims', Arg, ierr)
        CALL setup_params()
        CALL read_parameter_file_dimens()

      ELSE  !IF ( Process(:5)=='clean' ) THEN
        Process_flag = CLEAN

        IF ( Init_vars_from_file>OFF ) CLOSE ( Restart_inunit )
        IF ( Save_vars_to_file==ACTIVE ) THEN
          CALL PRMS_open_output_file(Restart_outunit, Var_save_file, 'var_save_file', 1, iret)
          IF ( iret/=0 ) ERROR STOP ERROR_open_out
          CALL call_modules_restart(SAVE_INIT)
        ENDIF
      ENDIF

      IF ( Model==DOCUMENTATION ) THEN
        IF ( Process_flag==SETDIMENS .OR. Process_flag==DECL ) THEN
          Init_vars_from_file = 0 ! make sure this is set so all variables and parameters are declared
          CALL module_doc()
          IF ( ierr/=0 ) CALL module_error('DOCUMENTATION', Arg, ierr)
          RETURN
        ELSE
          STOP
        ENDIF
      ENDIF

! All modules must be called for setdims, declare, initialize, and cleanup
      IF ( Process_flag/=RUN ) THEN
        ierr = basin()
        IF ( ierr/=0 ) CALL module_error('basin', Arg, ierr)

        IF ( Call_cascade==ACTIVE ) THEN
          ierr = cascade()
          IF ( ierr/=0 ) CALL module_error('cascade', Arg, ierr)
        ENDIF

        ierr = climateflow()
        IF ( ierr/=0 ) CALL module_error('climateflow', Arg, ierr)

        ierr = soltab()
        IF ( ierr/=0 ) CALL module_error('soltab', Arg, ierr)

        ierr = obs()
        IF ( ierr/=0 ) CALL module_error('obs', Arg, ierr)

        ierr = setup()
        IF ( ierr/=0 ) CALL module_error('setup', Arg, ierr)
      ENDIF

      ierr = prms_time()
      IF ( ierr/=0 ) CALL module_error('prms_time', Arg, ierr)

      IF ( Water_use_flag==ACTIVE ) THEN
        ierr = water_use_read()
        IF ( ierr/=0 ) CALL module_error('water_use_read', Arg, ierr)
      ENDIF

      IF ( Dynamic_flag==ACTIVE ) THEN
        ierr = dynamic_param_read()
        IF ( ierr/=0 ) CALL module_error('dynamic_param_read', Arg, ierr)
      ENDIF

      IF ( Climate_hru_flag==ACTIVE ) THEN
        ierr = climate_hru()
        IF ( ierr/=0 ) CALL module_error('climate_hru', Arg, ierr)
      ENDIF

      IF ( Climate_temp_flag==OFF ) THEN
        IF ( Temp_combined_flag==ACTIVE ) THEN
          ierr = temp_1sta_laps()
        ELSEIF ( Temp_flag==xyz_dist_module ) THEN
          ierr = xyz_dist()
        ELSEIF ( Temp_flag==temp_dist2_module ) THEN
          ierr = temp_dist2()
        ELSEIF ( Temp_flag==ide_dist_module ) THEN
          ierr = ide_dist()
        ELSE !IF ( Temp_flag==temp_map_module )
          CALL temp_map()
        ENDIF
        IF ( ierr/=0 ) CALL module_error(Temp_module, Arg, ierr)
      ENDIF

      IF ( Climate_precip_flag==OFF ) THEN
        IF ( Precip_combined_flag==ACTIVE ) THEN
          ierr = precip_1sta_laps()
        ELSEIF ( Precip_flag==precip_dist2_module ) THEN
          ierr = precip_dist2()
        ENDIF
        IF ( ierr/=0 ) CALL module_error(Precip_module, Arg, ierr)
      ENDIF

      IF ( Model==CLIMATE ) THEN
        IF ( Process_flag==RUN ) THEN
          CALL summary_output()
          RETURN
        ENDIF
      ENDIF

! frost_date is a pre-process module
      IF ( Model==FROST ) THEN
        ierr = frost_date()
        IF ( ierr/=0 ) CALL module_error('frost_date', Arg, ierr)
        IF ( Process_flag==RUN ) THEN
          CALL summary_output()
          RETURN
        ENDIF
        IF ( Process_flag==CLEAN ) STOP
      ENDIF

      IF ( Climate_swrad_flag==0 ) THEN
        IF ( Solrad_flag==ddsolrad_module ) THEN
          ierr = ddsolrad()
        ELSE !IF ( Solrad_flag==ccsolrad_module ) THEN
          ierr = ccsolrad()
        ENDIF
        IF ( ierr/=0 ) CALL module_error(Solrad_module, Arg, ierr)
      ENDIF

      IF ( Transp_flag==1 ) THEN
        ierr = transp_tindex()
      ELSEIF ( Transp_flag==2 ) THEN
        ierr = transp_frost()
      ENDIF
      IF ( ierr/=0 ) CALL module_error(Transp_module, Arg, ierr)

      IF ( Model==TRANSPIRE ) THEN
        IF ( Process_flag==RUN ) THEN
          CALL summary_output()
          RETURN
        ENDIF
      ENDIF

      IF ( Climate_potet_flag==OFF ) THEN
        IF ( Et_flag==potet_jh_module ) THEN
          ierr = potet_jh()
        ELSEIF ( Et_flag==potet_hamon_module ) THEN
          ierr = potet_hamon()
        ELSEIF ( Et_flag==potet_pan_module ) THEN
          ierr = potet_pan()
        ELSEIF ( Et_flag==potet_pt_module ) THEN
          ierr = potet_pt()
        ELSEIF ( Et_flag==potet_pm_sta_module ) THEN
          ierr = potet_pm_sta()
        ELSEIF ( Et_flag==potet_pm_module ) THEN
          ierr = potet_pm()
        ELSE !IF ( Et_flag==potet_hs_module ) THEN
          ierr = potet_hs()
        ENDIF
        IF ( ierr/=0 ) CALL module_error(Et_module, Arg, ierr)
      ENDIF

      IF ( Model==WRITE_CLIMATE ) THEN
        ierr = write_climate_hru()
        IF ( ierr/=0 ) CALL module_error('write_climate_hru', Arg, ierr)
        IF ( Process_flag==RUN ) RETURN
      ENDIF

      IF ( Model==POTET ) THEN
        IF ( Process_flag==RUN ) THEN
          CALL summary_output()
          RETURN
        ENDIF
      ENDIF

      ierr = intcp()
      IF ( ierr/=0 ) CALL module_error('intcp', Arg, ierr)

      IF ( no_snow_flag==OFF ) THEN
        ! rsr, need to do something if snow_cbh_flag=1
        ierr = snowcomp()
        IF ( ierr/=0 ) CALL module_error('snowcomp', Arg, ierr)

        IF ( Glacier_flag==ACTIVE ) THEN
          ierr = glacr()
          IF ( ierr/=0 ) CALL module_error('glacr', Arg, ierr)
        ENDIF
      ENDIF

      ierr = srunoff()
      IF ( ierr/=0 ) CALL module_error(Srunoff_module, Arg, ierr)

      IF ( AG_flag==ACTIVE ) THEN
        ierr = soilzone_ag()
      ELSE
        ierr = soilzone()
      ENDIF
      IF ( ierr/=0 ) CALL module_error(Soilzone_module, Arg, ierr)

      ! rsr, need to do something if gwflow_cbh_flag=1
      ierr = gwflow()
      IF ( ierr/=0 ) CALL module_error('gwflow', Arg, ierr)

      IF ( Stream_order_flag==ACTIVE ) THEN
        ierr = routing()
        IF ( ierr/=0 ) CALL module_error('routing', Arg, ierr)
      ENDIF

      IF ( Strmflow_flag==strmflow_noroute_module ) THEN
        ierr = strmflow()
      ELSEIF ( Muskingum_flag==ACTIVE ) THEN ! muskingum = 4; muskingum_mann = 7
        ierr = muskingum()
      ELSEIF ( Strmflow_flag==strmflow_in_out_module ) THEN
        ierr = strmflow_in_out()
      ELSEIF ( Strmflow_flag==strmflow_muskingum_lake_module ) THEN
        ierr = muskingum_lake()
      ENDIF
      IF ( ierr/=0 ) CALL module_error(Strmflow_module, Arg, ierr)

      IF ( Stream_temp_flag==ACTIVE ) ierr = stream_temp()

      IF ( Print_debug>DEBUG_minimum ) THEN
        ierr = basin_sum()
        IF ( ierr/=0 ) CALL module_error('basin_sum', Arg, ierr)
      ENDIF

      IF ( Print_debug==DEBUG_WB ) CALL water_balance()

      IF ( MapOutON_OFF>OFF ) THEN
        ierr = map_results()
        IF ( ierr/=0 ) CALL module_error('map_results', Arg, ierr)
      ENDIF

      IF ( Subbasin_flag==ACTIVE ) THEN
        ierr = subbasin()
        IF ( ierr/=0 ) CALL module_error('subbasin', Arg, ierr)
      ENDIF

      CALL summary_output()

      IF ( CsvON_OFF>OFF ) CALL prms_summary()

      IF ( ierr/=0 ) CALL module_error(MODNAME, Arg, ierr)
      IF ( Process_flag==RUN ) THEN
        RETURN
      ELSEIF ( Process_flag==CLEAN ) THEN
        CALL DATE_AND_TIME(VALUES=Elapsed_time_end)
        Execution_time_end = Elapsed_time_end(5)*3600 + Elapsed_time_end(6)*60 + &
     &                       Elapsed_time_end(7) + Elapsed_time_end(8)*0.001
        Elapsed_time = Execution_time_end - Execution_time_start
        Elapsed_time_minutes = INT(Elapsed_time/60.0)
        IF ( Print_debug>DEBUG_less ) THEN
          PRINT 9001
          PRINT 9003, 'start', (Elapsed_time_start(i),i=1,3), (Elapsed_time_start(i),i=5,7)
          PRINT 9003, 'end  ', (Elapsed_time_end(i),i=1,3), (Elapsed_time_end(i),i=5,7)
          PRINT '(A,I5,A,F6.2,A,/)', 'Execution elapsed time', Elapsed_time_minutes, ' minutes', &
     &                               Elapsed_time - Elapsed_time_minutes*60.0, ' seconds'
        ENDIF
        IF ( Print_debug>DEBUG_minimum ) &
     &       WRITE ( PRMS_output_unit,'(A,I5,A,F6.2,A,/)') 'Execution elapsed time', Elapsed_time_minutes, ' minutes', &
     &                                                     Elapsed_time - Elapsed_time_minutes*60.0, ' seconds'
        IF ( Print_debug>DEBUG_minimum ) CLOSE ( PRMS_output_unit )
        IF ( Save_vars_to_file==ACTIVE ) CLOSE ( Restart_outunit )
      ELSEIF ( Process_flag==DECL ) THEN
        CALL read_parameter_file_params()
        IF ( Print_debug>DEBUG_minimum ) THEN
          PRINT '(A)', EQULS(:62)
          WRITE ( PRMS_output_unit, '(A)' ) EQULS(:62)
        ENDIF
        IF ( Model==CONVERT ) CALL convert_params()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Inputerror_flag==1 ) THEN
          PRINT '(//,A,//,A,/,A,/,A)', '**Fix input errors in your Parameter File to continue**', &
     &          '  Set control parameter parameter_check_flag to 0 after', &
     &          '  all parameter values are valid.'
          PRINT '(/,A,/,A,/,A,/,A,/,A,/)', &
     &          'If input errors are related to paramters used for automated', &
     &          'calibration processes, with CAUTION, set control parameter', &
     &          'parameter_check_flag to 0. After calibration set the', &
     &          'parameter_check_flag to 1 to verify that those calibration', &
     &          'parameters have valid and compatible values.'
        ENDIF
        IF ( Parameter_check_flag==2 ) STOP
        IF ( Inputerror_flag==1 ) ERROR STOP ERROR_param
        IF ( Model==CONVERT ) THEN
          CALL convert_params()
          STOP
        ENDIF
        IF ( Print_debug>DEBUG_minimum ) &
     &       PRINT 4, 'Simulation time period:', Start_year, Start_month, Start_day, ' -', End_year, End_month, End_day, EQULS
      ENDIF

    4 FORMAT (/, 2(A, I5, 2('/',I2.2)), //, A, /)
 9001 FORMAT (/, 26X, 25('='), /, 26X, 'Normal completion of PRMS', /, 26X, 25('='), /)
 9002 FORMAT (//, 74('='), /, 'Please give careful consideration to fixing all ERROR and WARNING messages', /, 74('='))
 9003 FORMAT ('Execution ', A, ' date and time (yyyy/mm/dd hh:mm:ss)', I5, 2('/',I2.2), I3, 2(':',I2.2), /)
 9004 FORMAT (/, 2A)

      END SUBROUTINE call_modules

!***********************************************************************
!     declare the dimensions
!***********************************************************************
      SUBROUTINE setdims()
      USE PRMS_CONSTANTS, ONLY: ERROR_control
      use PRMS_CONTROL_FILE, only: get_control_arguments, read_control_file, control_integer, control_integer_array, control_string !, control_file_name
      USE PRMS_MODULE
      use PRMS_READ_PARAM_FILE, only: decldim, declfix, read_parameter_file_dimens, setup_dimens
      use prms_utils, only: compute_julday, error_stop, module_error, PRMS_open_input_file, PRMS_open_output_file, read_error
      IMPLICIT NONE
! Functions
      EXTERNAL :: check_module_names
! Local Variables
      ! Maximum values are no longer limits
! Local Variables
      INTEGER :: idim, iret, j
!***********************************************************************
      Inputerror_flag = 0

      CALL read_control_file()
      CALL get_control_arguments()

      ! debug print flag:
      ! -1=quiet - reduced screen output (DEBUG_less)
      ! 0=none; 1=water balances; 2=basin;
      ! 4=basin_sum; 5=soltab; 7=soil zone;
      ! 9=snowcomp; 13=cascade; 14=subbasin tree
      IF ( control_string(Data_file, 'data_file')/=0 ) CALL read_error(5, 'data_file')
      IF ( control_integer(Print_debug, 'print_debug')/=0 ) Print_debug = 0

      IF ( control_integer(Parameter_check_flag, 'parameter_check_flag')/=0 ) Parameter_check_flag = 1

      IF ( control_string(Model_mode, 'model_mode')/=0 ) CALL read_error(5, 'model_mode')
      IF ( Model_mode(:4)=='    ' ) Model_mode = 'PRMS5'
      PRMS4_flag = ACTIVE
      IF ( Model_mode(:5)=='PRMS5' .OR. Model_mode(:5)=='prms5' ) PRMS4_flag = OFF
      IF ( Model_mode(:4)=='PRMS' .OR. Model_mode(:4)=='prms' .OR. Model_mode(:5)=='DAILY' ) THEN
        Model = PRMS
      ELSEIF ( Model_mode(:5)=='FROST' ) THEN
        Model = FROST
      ELSEIF ( Model_mode(:13)=='WRITE_CLIMATE' ) THEN
        Model = WRITE_CLIMATE
      ELSEIF ( Model_mode(:7)=='CLIMATE' ) THEN
        Model = CLIMATE
      ELSEIF ( Model_mode(:5)=='POTET' ) THEN
        Model = POTET
      ELSEIF ( Model_mode(:9)=='TRANSPIRE' ) THEN
        Model = TRANSPIRE
      ELSEIF ( Model_mode(:7)=='CONVERT' ) THEN ! can be CONVERT4 or CONVERT5 or CONVERT (=CONVERT5)
        Model = CONVERT
      ELSEIF ( Model_mode(:13)=='DOCUMENTATION' ) THEN
        Model = DOCUMENTATION
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid model_mode value: ', Model_mode
        Inputerror_flag = 1
      ENDIF

      ! get simulation start_time and end_time
      Starttime = -1
      DO j = 1, 6
        IF ( control_integer_array(Starttime(j), j, 'start_time')/=0 ) THEN
          PRINT *, 'ERROR, start_time, index:', j, 'value: ', Starttime(j)
          Inputerror_flag = 1
        ENDIF
      ENDDO
      Start_year = Starttime(1)
      IF ( Start_year<0 ) THEN
        PRINT *, 'ERROR, control parameter start_time must be specified'
        Inputerror_flag = 1
      ENDIF
      Start_month = Starttime(2)
      Start_day = Starttime(3)
      Nowyear = Start_year
      Nowmonth = Start_month
      Nowday = Start_day
      Endtime = -1
      DO j = 1, 6
        IF ( control_integer_array(Endtime(j), j, 'end_time')/=0 ) THEN
          PRINT *, 'ERROR, end_time, index:', j, 'value: ', Endtime(j)
          Inputerror_flag = 1
        ENDIF
      ENDDO
      End_year = Endtime(1)
      IF ( End_year<0 ) THEN
        PRINT *, 'ERROR, control parameter start_time must be specified'
        Inputerror_flag = 1
      ENDIF
      End_month = Endtime(2)
      End_day = Endtime(3)

      startday = compute_julday(Start_year, Start_month, Start_day)
      endday = compute_julday(End_year, End_month, End_day)
      Number_timesteps = endday - startday + 1

      IF ( control_integer(Init_vars_from_file, 'init_vars_from_file')/=0 ) Init_vars_from_file = 0
      IF ( control_integer(Save_vars_to_file, 'save_vars_to_file')/=0 ) Save_vars_to_file = 0

      CALL setup_dimens()

      ! Open PRMS module output file
      IF ( control_string(Model_output_file, 'model_output_file')/=0 ) CALL read_error(5, 'model_output_file')
      IF ( Print_debug>DEBUG_minimum ) THEN
        CALL PRMS_open_output_file(PRMS_output_unit, Model_output_file, 'model_output_file', 0, iret)
        IF ( iret/=0 ) Inputerror_flag = 1
      ENDIF
!      IF ( control_file_name(Model_control_file)/=0 ) CALL read_error(5, 'control_file_name')
      IF ( control_string(Param_file, 'param_file')/=0 ) CALL read_error(5, 'param_file')

      ! Check for restart files
      IF ( Init_vars_from_file>0 ) THEN
        IF ( control_string(Var_init_file, 'var_init_file')/=0 ) CALL read_error(5, 'var_init_file')
        CALL PRMS_open_input_file(Restart_inunit, Var_init_file, 'var_init_file', 1, iret)
        IF ( iret/=0 ) Inputerror_flag = 1
      ENDIF
      IF ( Save_vars_to_file==ACTIVE ) THEN
        IF ( control_string(Var_save_file, 'var_save_file')/=0 ) CALL read_error(5, 'var_save_file')
      ENDIF

      Temp_module = 'temp_1sta'
      IF ( control_string(Temp_module, 'temp_module')/=0 ) CALL read_error(5, 'temp_module')
      Precip_module = 'precip_1sta'
      IF ( control_string(Precip_module, 'precip_module')/=0 ) CALL read_error(5, 'precip_module')
      Transp_module = 'transp_index'
      IF ( control_string(Transp_module, 'transp_module')/=0 ) CALL read_error(5, 'transp_module')
      Et_module = 'potet_jh'
      IF ( control_string(Et_module, 'et_module')/=0 ) CALL read_error(5, 'et_module')
      Srunoff_module = 'srunoff_smidx'
      IF ( control_string(Srunoff_module, 'srunoff_module')/=0 ) CALL read_error(5, 'srunoff_module')
      Solrad_module = 'ddsolrad'
      IF ( control_string(Solrad_module, 'solrad_module')/=0 ) CALL read_error(5, 'solrad_module')
      Soilzone_module = 'soilzone'
      IF ( control_string(Soilzone_module, 'soilzone_module')/=0 ) CALL read_error(5, 'soilzone_module')
      AG_flag = OFF
      IF ( Soilzone_module=='soilzone_ag' ) AG_flag = ACTIVE
      Strmflow_module = 'strmflow'
      IF ( control_string(Strmflow_module, 'strmflow_module')/=0 ) CALL read_error(5, 'strmflow_module')
      IF ( control_string(irrigated_area_module, 'irrigated_area_module')/=0 ) CALL read_error(5, 'irrigated_area_module')
      IF ( control_string(AET_module, 'AET_module')/=0 ) CALL read_error(5, 'AET_module')
      IF ( control_string(PET_ag_module, 'PET_ag_module')/=0 ) CALL read_error(5, 'PET_ag_module')
      IF ( irrigated_area_module(:11)=='climate_hru' ) irrigated_area_flag = ACTIVE
      IF ( AET_module(:11)=='climate_hru' ) AET_cbh_flag = ACTIVE
      IF ( PET_ag_module(:11)=='climate_hru' ) PET_cbh_flag = ACTIVE

      IF ( Parameter_check_flag>0 ) CALL check_module_names()

      Climate_precip_flag = OFF
      Climate_temp_flag = OFF
      Climate_transp_flag = OFF
      Climate_potet_flag = OFF
      Climate_swrad_flag = OFF

      IF ( Precip_module(:11)=='precip_1sta' .OR. Precip_module(:11)=='precip_prms') THEN
        Precip_flag = precip_1sta_module
      ELSEIF ( Precip_module(:11)=='precip_laps' ) THEN
        Precip_flag = precip_laps_module
      ELSEIF ( Precip_module(:12)=='precip_dist2' ) THEN
        Precip_flag = precip_dist2_module
      ELSEIF ( Precip_module(:8)=='ide_dist' ) THEN
        Precip_flag = ide_dist_module
      ELSEIF ( Precip_module(:11)=='climate_hru' ) THEN
        Precip_flag = climate_hru_module
        Climate_precip_flag = 1
      ELSEIF ( Precip_module(:8)=='xyz_dist' ) THEN
        Precip_flag = xyz_dist_module
      ELSEIF ( Precip_module(:15)=='precip_map' ) THEN
        Precip_flag = precip_map_module
      ELSE
        PRINT '(/,2A)', 'ERROR: invalid precip_module value: ', Precip_module
        Inputerror_flag = 1
      ENDIF
      Precip_combined_flag = OFF
      IF ( Precip_flag==precip_1sta_module .OR. Precip_flag==precip_laps_module ) Precip_combined_flag = 1

      IF ( Temp_module(:9)=='temp_1sta' ) THEN
        Temp_flag = temp_1sta_module
      ELSEIF ( Temp_module(:9)=='temp_laps' ) THEN
        Temp_flag = temp_laps_module
      ELSEIF ( Temp_module(:10)=='temp_dist2' ) THEN
        Temp_flag = temp_dist2_module
      ELSEIF ( Temp_module(:8)=='ide_dist' ) THEN
        Temp_flag = ide_dist_module
      ELSEIF ( Temp_module(:11)=='climate_hru' ) THEN
        Temp_flag = climate_hru_module
        Climate_temp_flag = ACTIVE
      ELSEIF ( Temp_module(:8)=='xyz_dist' ) THEN
        Temp_flag = xyz_dist_module
      ELSEIF ( Temp_module(:8)=='temp_sta' ) THEN
        Temp_flag = temp_sta_module
      ELSEIF ( Temp_module(:15)=='temp_map' ) THEN
        Temp_flag = temp_map_module
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid temp_module value: ', Temp_module
        Inputerror_flag = 1
      ENDIF
      Temp_combined_flag = OFF
      IF ( Temp_flag==temp_1sta_module .OR. Temp_flag==temp_laps_module .OR. Temp_flag==temp_sta_module ) Temp_combined_flag = 1

      IF ( Transp_module(:13)=='transp_tindex' ) THEN
        Transp_flag = 1
      ELSEIF ( Transp_module(:12)=='transp_frost' ) THEN
        Transp_flag = 2
      ELSEIF ( Transp_module(:11)=='climate_hru' ) THEN
        Transp_flag = 3
        Climate_transp_flag = 1
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid transp_module value: ', Transp_module
        Inputerror_flag = 1
      ENDIF

      IF ( Et_module(:8)=='potet_jh' ) THEN
        Et_flag = potet_jh_module
      ELSEIF ( Et_module(:11)=='potet_hamon' ) THEN
        Et_flag = potet_hamon_module
      ELSEIF ( Et_module(:11)=='climate_hru' ) THEN
        Et_flag = climate_hru_module
        Climate_potet_flag = ACTIVE
      ELSEIF ( Et_module(:8)=='potet_hs' ) THEN
        Et_flag = potet_hs_module
      ELSEIF ( Et_module(:12)=='potet_pm_sta' ) THEN
        Et_flag = potet_pm_sta_module
      ELSEIF ( Et_module(:8)=='potet_pm' ) THEN
        Et_flag = potet_pm_module
      ELSEIF ( Et_module(:8)=='potet_pt' ) THEN
        Et_flag = potet_pt_module
      ELSEIF ( Et_module(:9)=='potet_pan' ) THEN
        Et_flag = potet_pan_module
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid et_module value: ', Et_module
        Inputerror_flag = 1
      ENDIF

      ! stream_temp
      IF ( control_integer(Stream_temp_flag, 'stream_temp_flag')/=0 ) Stream_temp_flag = OFF
      ! 0 = CBH File; 1 = specified constant; 2 = Stations
      IF ( control_integer(Strmtemp_humidity_flag, 'strmtemp_humidity_flag')/=0 ) Strmtemp_humidity_flag = OFF

      IF ( control_integer(Snarea_curve_flag, 'snarea_curve_flag')/=0 ) Snarea_curve_flag = OFF
      IF ( control_integer(Soilzone_aet_flag, 'soilzone_aet_flag')/=0 ) Soilzone_aet_flag = OFF
      IF ( control_integer(Iter_aet_flag, 'iter_aet_flag')/=0 ) Iter_aet_flag = OFF
      IF ( control_integer(snow_cloudcover_flag, 'snow_cloudcover_flag')/=0 ) snow_cloudcover_flag = OFF

      IF ( control_integer(Humidity_cbh_flag, 'humidity_cbh_flag')/=0 ) Humidity_cbh_flag = OFF
      IF ( control_integer(Windspeed_cbh_flag, 'windspeed_cbh_flag')/=0 ) Windspeed_cbh_flag = OFF
      IF ( control_integer(Albedo_cbh_flag, 'albedo_cbh_flag')/=0 ) Albedo_cbh_flag = OFF
      IF ( control_integer(Cloud_cover_cbh_flag, 'cloud_cover_cbh_flag')/=0 ) Cloud_cover_cbh_flag = OFF
      IF ( Et_flag==potet_pm_module .OR. Et_flag==potet_pt_module .OR. &
     &     (Stream_temp_flag==ACTIVE .AND. Strmtemp_humidity_flag==OFF) ) Humidity_cbh_flag = ACTIVE
      IF ( Et_flag==potet_pm_module ) Windspeed_cbh_flag = ACTIVE

      IF ( Srunoff_module(:13)=='srunoff_smidx' ) THEN
        Sroff_flag = smidx_module
      ELSEIF ( Srunoff_module(:13)=='srunoff_carea' ) THEN
        Sroff_flag = carea_module
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid srunoff_module value: ', Srunoff_module
        Inputerror_flag = 1
      ENDIF

      IF ( control_integer(Orad_flag, 'orad_flag')/=0 ) Orad_flag = OFF
      IF ( Solrad_module(:8)=='ddsolrad' ) THEN
        Solrad_flag = ddsolrad_module
      ELSEIF ( Solrad_module(:11)=='climate_hru' ) THEN
        Solrad_flag = climate_hru_module
        Climate_swrad_flag = ACTIVE
      ELSEIF ( Solrad_module(:8)=='ccsolrad' ) THEN
        Solrad_flag = ccsolrad_module
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid solrad_module value: ', Solrad_module
        Inputerror_flag = 1
      ENDIF

      IF ( control_integer(Snow_cbh_flag, 'snow_cbh_flag')/=0 ) Snow_cbh_flag = OFF
      IF ( control_integer(Gwflow_cbh_flag, 'gwflow_cbh_flag')/=0 ) Gwflow_cbh_flag = OFF
      Snow_cbh_flag = OFF ! not implemented yet
      Gwflow_cbh_flag = OFF ! not implemented yet

      Climate_hru_flag = OFF
      IF ( Climate_temp_flag==ACTIVE .OR. Climate_precip_flag==ACTIVE .OR. Climate_potet_flag==ACTIVE .OR. &
     &     Climate_swrad_flag==ACTIVE .OR. Climate_transp_flag==ACTIVE .OR. &
     &     Humidity_cbh_flag==ACTIVE .OR. Windspeed_cbh_flag==ACTIVE .OR. &
     &     irrigated_area_flag==ACTIVE .OR. AET_cbh_flag==ACTIVE .OR. PET_cbh_flag==ACTIVE .OR. &
     &     Albedo_cbh_flag==ACTIVE .OR. Cloud_cover_cbh_flag==ACTIVE .OR. &
     &     Gwflow_cbh_flag==ACTIVE .OR. Snow_cbh_flag==ACTIVE ) Climate_hru_flag = ACTIVE

      Muskingum_flag = OFF
      IF ( Strmflow_module(:15)=='strmflow_in_out' ) THEN
        Strmflow_flag = strmflow_in_out_module
      ELSEIF ( Strmflow_module(:14)=='muskingum_lake' ) THEN
        Strmflow_flag = strmflow_muskingum_lake_module
      ELSEIF ( Strmflow_module(:13)=='strmflow_lake' ) THEN
        PRINT '(/,2A)', 'ERROR, invalid strmflow_module value: ', Strmflow_module
        Inputerror_flag = 1
      ELSEIF ( Strmflow_module(:8)=='strmflow' ) THEN
        Strmflow_flag = strmflow_noroute_module
      ELSEIF ( Strmflow_module(:14)=='muskingum_mann' ) THEN
        Strmflow_flag = strmflow_muskingum_mann_module
        Muskingum_flag = ACTIVE
      ELSEIF ( Strmflow_module(:9)=='muskingum' ) THEN
        Strmflow_flag = strmflow_muskingum_module
        Muskingum_flag = ACTIVE
      ELSE
        PRINT '(/,2A)', 'ERROR, invalid strmflow_module value: ', Strmflow_module
        Inputerror_flag = 1
      ENDIF

! cascade dimensions
      IF ( decldim('ncascade', 0, MAXDIM, &
     &     'Number of HRU links for cascading flow')/=0 ) CALL read_error(7, 'ncascade')
      IF ( decldim('ncascdgw', 0, MAXDIM, &
     &     'Number of GWR links for cascading flow')/=0 ) CALL read_error(7, 'ncascdgw')

! nsegment dimension
      IF ( decldim('nsegment', 0, MAXDIM, 'Number of stream-channel segments')/=0 ) CALL read_error(7, 'nsegment')

! subbasin dimensions
      IF ( control_integer(Subbasin_flag, 'subbasin_flag')/=0 ) Subbasin_flag = ACTIVE
      IF ( decldim('nsub', 0, MAXDIM, 'Number of internal subbasins')/=0 ) CALL read_error(7, 'nsub')

      IF ( control_integer(Dprst_flag, 'dprst_flag')/=0 ) Dprst_flag = OFF
      IF ( control_integer(Dprst_transfer_water_use, 'dprst_transfer_water_use')/=0 ) Dprst_transfer_water_use = OFF
      IF ( control_integer(Dprst_add_water_use, 'dprst_add_water_use')/=0 ) Dprst_add_water_use = OFF

      ! 0 = off, 1 = on, 2 = lauren version
      IF ( control_integer(CsvON_OFF, 'csvON_OFF')/=0 ) CsvON_OFF = OFF

! map results dimensions
      IF ( control_integer(MapOutON_OFF, 'mapOutON_OFF')/=0 ) MapOutON_OFF = OFF
      idim = 0
      IF ( MapOutON_OFF>OFF ) idim = 1
      IF ( decldim('nhrucell', idim, MAXDIM, &
     &     'Number of unique intersections between HRUs and spatial units of a target map for mapped results')/=0 ) &
     &     CALL read_error(7, 'nhrucell')
      IF ( decldim('ngwcell', 0, MAXDIM, &
     &     'Number of spatial units in the target map for mapped results')/=0 ) CALL read_error(7, 'ngwcell')

! declare precip_map and temp_map module specific dimensions
      IF ( decldim('nmap2hru', 0, MAXDIM, 'Number of intersections between HRUs and input climate map')/=0 ) &
     &     CALL read_error(7, 'nmap2hru')
      IF ( decldim('nmap', 0, MAXDIM, 'Number of mapped values')/=0 ) CALL read_error(7, 'nmap')
      IF ( control_integer(Glacier_flag, 'glacier_flag')/=0 ) Glacier_flag = OFF
      IF ( control_integer(no_snow_flag, 'no_snow_flag')/=0 ) no_snow_flag = OFF
      IF ( control_integer(Frozen_flag, 'frozen_flag')/=0 ) Frozen_flag = OFF
      IF ( control_integer(Dyn_imperv_flag, 'dyn_imperv_flag')/=0 ) Dyn_imperv_flag = OFF
      IF ( control_integer(Dyn_intcp_flag, 'dyn_intcp_flag')/=0 ) Dyn_intcp_flag = OFF
      IF ( control_integer(Dyn_covden_flag, 'dyn_covden_flag')/=0 ) Dyn_covden_flag = OFF
      IF ( control_integer(Dyn_dprst_flag, 'dyn_dprst_flag')/=0 ) Dyn_dprst_flag = OFF
      IF ( control_integer(Dyn_potet_flag, 'dyn_potet_flag')/=0 ) Dyn_potet_flag = OFF
      IF ( control_integer(Dyn_covtype_flag, 'dyn_covtype_flag')/=0 ) Dyn_covtype_flag = OFF
      IF ( control_integer(Dyn_transp_flag, 'dyn_transp_flag')/=0 ) Dyn_transp_flag = OFF
      IF ( control_integer(Dyn_soil_flag, 'dyn_soil_flag')/=0 ) Dyn_soil_flag = OFF
      IF ( control_integer(Dyn_radtrncf_flag, 'dyn_radtrncf_flag')/=0 ) Dyn_radtrncf_flag = OFF
      IF ( control_integer(Dyn_sro2dprst_perv_flag, 'dyn_sro2dprst_perv_flag')/=0 ) Dyn_sro2dprst_perv_flag = OFF
      IF ( control_integer(Dyn_sro2dprst_imperv_flag, 'dyn_sro2dprst_imperv_flag')/=0 ) Dyn_sro2dprst_imperv_flag = OFF
      IF ( control_integer(Dyn_fallfrost_flag, 'dyn_fallfrost_flag')/=0 ) Dyn_fallfrost_flag = OFF
      IF ( control_integer(Dyn_springfrost_flag, 'dyn_springfrost_flag')/=0 ) Dyn_springfrost_flag = OFF
      IF ( control_integer(Dyn_snareathresh_flag, 'dyn_snareathresh_flag')/=0 ) Dyn_snareathresh_flag = OFF
      IF ( control_integer(Dyn_transp_on_flag, 'dyn_transp_on_flag')/=0 ) Dyn_transp_on_flag = OFF
      IF ( control_integer(Dyn_ag_frac_flag, 'dyn_ag_frac_flag')/=0 ) Dyn_ag_frac_flag = OFF
      IF ( control_integer(Dyn_ag_soil_flag, 'dyn_ag_soil_flag')/=0 ) Dyn_ag_soil_flag = OFF
      Dynamic_flag = 0
      IF ( Dyn_imperv_flag/=OFF .OR. Dyn_intcp_flag/=0 .OR. Dyn_covden_flag/=0 .OR. Dyn_dprst_flag/=OFF .OR. &
     &     Dyn_potet_flag/=OFF .OR. Dyn_covtype_flag/=0 .OR. Dyn_transp_flag/=0 .OR. Dyn_soil_flag /=OFF .OR. &
     &     Dyn_radtrncf_flag/=OFF .OR. Dyn_sro2dprst_perv_flag/=0 .OR. Dyn_sro2dprst_imperv_flag/=OFF .OR. &
     &     Dyn_fallfrost_flag/=OFF .OR. Dyn_springfrost_flag/=0 .OR. Dyn_snareathresh_flag/=0 .OR. &
     &     Dyn_transp_on_flag/=OFF .OR. Dyn_ag_frac_flag==ACTIVE .OR. Dyn_ag_soil_flag==ACTIVE ) Dynamic_flag = ACTIVE
      IF ( control_integer(Gwr_transferON_OFF, 'gwr_transferON_OFF')/=0) Gwr_transferON_OFF = OFF
      IF ( control_integer(External_transferON_OFF, 'external_transferON_OFF')/=0 ) External_transferON_OFF = OFF
      IF ( control_integer(Dprst_transferON_OFF, 'dprst_transferON_OFF')/=0 ) Dprst_transferON_OFF = OFF
      IF ( control_integer(Segment_transferON_OFF, 'segment_transferON_OFF')/=0 ) Segment_transferON_OFF = OFF
      IF ( control_integer(Lake_transferON_OFF, 'lake_transferON_OFF')/=0 ) Lake_transferON_OFF = OFF
      IF ( control_integer(Gwr_swale_flag, 'gwr_swale_flag')/=0 ) Gwr_swale_flag = OFF

! nhru_summary
      IF ( control_integer(NhruOutON_OFF, 'nhruOutON_OFF')/=0 ) NhruOutON_OFF = 0

! nsub_summary
      IF ( control_integer(NsubOutON_OFF, 'nsubOutON_OFF')/=0 ) NsubOutON_OFF = 0

! basin_summary
      IF ( control_integer(BasinOutON_OFF, 'basinOutON_OFF')/=0 ) BasinOutON_OFF = 0

! nsegment_summary
      IF ( control_integer(NsegmentOutON_OFF, 'nsegmentOutON_OFF')/=0 ) NsegmentOutON_OFF = 0

! statvar_out
      IF ( control_integer(statsON_OFF, 'statsON_OFF')/=0 ) statsON_OFF = 0
      IF ( statsON_OFF==ACTIVE ) THEN
        IF ( control_string(stat_var_file, 'stat_var_file')/=0 ) CALL read_error(5, 'stat_var_file')
      ENDIF

      IF ( control_integer(Prms_warmup, 'prms_warmup')/=0 ) Prms_warmup = 0
      IF ( NhruOutON_OFF>0 .OR. NsubOutON_OFF>0 .OR. BasinOutON_OFF>0 .OR. NsegmentOutON_OFF>0 ) THEN
        IF ( Start_year+Prms_warmup>End_year ) THEN ! change to start full date ???
          PRINT *, 'ERROR, prms_warmup > than simulation time period:', Prms_warmup
          Inputerror_flag = 1
        ENDIF
      ENDIF

! cascade
      ! if cascade_flag = 2 (CASCADE_HRU_SEGMENT), use hru_segment parameter for cascades, ncascade=ncascdgw=nhru (typical polygon HRUs)
      IF ( control_integer(Cascade_flag, 'cascade_flag')/=0 ) Cascade_flag = CASCADE_NORMAL
      ! if cascadegw_flag = 2 (CASCADEGW_SAME), use same cascades as HRUs
      IF ( control_integer(Cascadegw_flag, 'cascadegw_flag')/=0 ) Cascadegw_flag = CASCADE_NORMAL

! spatial units
      IF ( decldim('ngw', 1, MAXDIM, 'Number of GWRs')/=0 ) CALL read_error(7, 'ngw')
      IF ( decldim('nhru', 1, MAXDIM, 'Number of HRUs')/=0 ) CALL read_error(7, 'nhru')
      IF ( decldim('nssr', 1, MAXDIM, 'Number of subsurface reservoirs')/=0 ) CALL read_error(7, 'nssr')
      IF ( decldim('nlake', 0, MAXDIM, 'Number of lakes')/=0 ) CALL read_error(7, 'nlake')
      IF ( decldim('nlake_hrus', 0, MAXDIM, 'Number of lake HRUs')/=0 ) CALL read_error(7, 'nlake_hrus')
      IF ( decldim('npoigages', 0, MAXDIM, 'Number of POI gages')/=0 ) CALL read_error(7, 'npoigages')

! Time-series data stations, need to know if in Data File
      IF ( decldim('nrain', 0, MAXDIM, 'Number of precipitation-measurement stations')/=0 ) CALL read_error(7, 'nrain')
      IF ( decldim('nsol', 0, MAXDIM, 'Number of solar-radiation measurement stations')/=0 ) CALL read_error(7, 'nsol')
      IF ( decldim('ntemp', 0, MAXDIM, 'Number of air-temperature-measurement stations')/=0 ) CALL read_error(7, 'ntemp')
      IF ( decldim('nobs', 0, MAXDIM, 'Number of streamflow-measurement stations')/=0 ) CALL read_error(7, 'nobs')
      IF ( decldim('nevap', 0, MAXDIM, 'Number of pan-evaporation data sets')/=0 ) CALL read_error(7, 'nevap')
      IF ( decldim('nratetbl', 0, MAXDIM, 'Number of rating-table data sets for lake elevations') &
     &     /=0 ) CALL read_error(7, 'nratetbl')
      IF ( decldim('nsnow', 0, MAXDIM, 'Number of snow-depth-measurement stations')/=0 ) CALL read_error(7, 'nsnow')

! depletion curves
      IF ( decldim('ndepl', 1, MAXDIM, 'Number of snow-depletion curves')/=0 ) CALL read_error(7, 'ndelp')
      IF ( decldim('ndeplval', 11, MAXDIM, 'Number of values in all snow-depletion curves (set to ndepl*11)')/=0 ) &
     &     CALL read_error(7, 'ndelplval')

! water-use
      IF ( decldim('nwateruse', 0, MAXDIM, 'Number of water-use data sets')/=0 ) CALL read_error(7, 'nwateruse')
      IF ( decldim('nexternal', 0, MAXDIM, &
     &      'Number of external water-use sources or destinations')/=0 ) CALL read_error(7, 'nexternal')
      IF ( decldim('nconsumed', 0, MAXDIM, 'Number of consumptive water-use destinations')/=0 ) CALL read_error(7, 'nconsumed')

! fixed dimensions
      IF ( declfix('ndays', MAX_DAYS_PER_YEAR, MAX_DAYS_PER_YEAR, 'Maximum number of days in a year ')/=0 ) &
     &     CALL read_error(7, 'ndays')
      IF ( declfix('nmonths', 12, 12, 'Number of months in a year')/=0 ) CALL read_error(7, 'nmonths')
      IF ( declfix('one', 1, 1, 'Number of values for scaler array')/=0 ) CALL read_error(7, 'one')

      IF ( Inputerror_flag==1 ) THEN
        PRINT '(//,A,/,A)', '**FIX input errors in your Control File to continue**', &
     &        'NOTE: some errors may be due to use of defalut values'
        STOP
      ENDIF

      END SUBROUTINE setdims

!***********************************************************************
!     Get and check consistency of dimensions with flags
!***********************************************************************
      INTEGER FUNCTION check_dims()
      use PRMS_READ_PARAM_FILE, only: getdim
      USE PRMS_MODULE
      use prms_utils, only: read_error
      IMPLICIT NONE
      EXTERNAL :: check_dimens
!***********************************************************************

      Nhru = getdim('nhru')
      IF ( Nhru==-1 ) CALL read_error(7, 'nhru')

      Nssr = getdim('nssr')
      IF ( Nssr==-1 ) CALL read_error(7, 'nssr')

      Ngw = getdim('ngw')
      IF ( Ngw==-1 ) CALL read_error(7, 'ngw')

      Ntemp = getdim('ntemp')
      IF ( Ntemp==-1 ) CALL read_error(6, 'ntemp')

      Nrain = getdim('nrain')
      IF ( Nrain==-1 ) CALL read_error(6, 'nrain')

      Nsol = getdim('nsol')
      IF ( Nsol==-1 ) CALL read_error(6, 'nsol')

      Nobs = getdim('nobs')
      IF ( Nobs==-1 ) CALL read_error(6, 'nobs')

      Nevap = getdim('nevap')
      IF ( Nevap==-1 ) CALL read_error(6, 'nevap')

      Ncascade = getdim('ncascade')
      IF ( Ncascade==-1 ) CALL read_error(7, 'ncascade')
      Ncascdgw = getdim('ncascdgw')
      IF ( Ncascdgw==-1 ) CALL read_error(7, 'ncascdgw')
      IF ( Cascade_flag==CASCADE_HRU_SEGMENT ) THEN
        Ncascade = Nhru
        Cascadegw_flag = CASCADEGW_SAME
      ENDIF
      IF ( Cascadegw_flag==CASCADEGW_SAME ) Ncascdgw = Ncascade
      IF ( Ncascade==0 ) Cascade_flag = CASCADE_OFF
      IF ( Ncascdgw==0 ) Cascadegw_flag = CASCADEGW_OFF
      IF ( (Cascade_flag>CASCADE_OFF .OR. Cascadegw_flag>CASCADEGW_OFF) .AND. Model/=CONVERT ) THEN ! don't call if model_mode = CONVERT
        Call_cascade = ACTIVE
      ELSE
        Call_cascade = OFF
      ENDIF

      Nwateruse = getdim('nwateruse')
      IF ( Nwateruse==-1 ) CALL read_error(7, 'nwateruse')

      Nexternal = getdim('nexternal')
      IF ( Nexternal==-1 ) CALL read_error(6, 'nexternal')

      Nconsumed = getdim('nconsumed')
      IF ( Nconsumed==-1 ) CALL read_error(6, 'nconsumed')

      Npoigages = getdim('npoigages')
      IF ( Npoigages==-1 ) CALL read_error(6, 'npoigages')

      Nlake = getdim('nlake')
      IF ( Nlake==-1 ) CALL read_error(7, 'nlake')

      Nlake_hrus = getdim('nlake_hrus')
      IF ( Nlake_hrus==-1 ) CALL read_error(7, 'nlake_hrus')
      IF ( Nlake>0 .AND. Nlake_hrus==0 ) Nlake_hrus = Nlake

      Ndepl = getdim('ndepl')
      IF ( Ndepl==-1 ) CALL read_error(7, 'ndepl')

      Ndeplval = getdim('ndeplval')
      IF ( Ndeplval==-1 ) CALL read_error(7, 'ndeplval')

      Nsub = getdim('nsub')
      IF ( Nsub==-1 ) CALL read_error(7, 'nsub')
      ! default = 1, turn off if no subbasins
      IF ( Subbasin_flag==1 .AND. Nsub==0 ) Subbasin_flag = 0

      Nsegment = getdim('nsegment')
      IF ( Nsegment==-1 ) CALL read_error(7, 'nsegment')

      Nhrucell = getdim('nhrucell')
      IF ( Nhrucell==-1 ) CALL read_error(6, 'nhrucell')

      Ngwcell = getdim('ngwcell')
      IF ( Ngwcell==-1 ) CALL read_error(6, 'ngwcell')

      Nratetbl = getdim('nratetbl')
      IF ( Nratetbl==-1 ) CALL read_error(6, 'nratetbl')

      Nmap2hru = getdim('nmap2hru')
      IF ( Nmap2hru==-1 ) CALL read_error(6, 'nmap2hru')
      Nmap = getdim('nmap')
      IF ( Nmap==-1 ) CALL read_error(6, 'nmap')

      Water_use_flag = OFF
      IF ( Nwateruse>0 ) THEN
        IF ( Segment_transferON_OFF==ACTIVE .OR. Gwr_transferON_OFF==ACTIVE .OR. External_transferON_OFF==ACTIVE .OR. &
     &       Dprst_transferON_OFF==ACTIVE .OR. Lake_transferON_OFF==ACTIVE ) Water_use_flag = ACTIVE
        IF ( Water_use_flag==OFF ) THEN
          PRINT *, 'WARNING, nwateruse specified > 0 without transfers active'
          Nwateruse = 0
        ENDIF
      ENDIF

      IF ( Segment_transferON_OFF==ACTIVE .OR. Gwr_transferON_OFF==ACTIVE .OR. External_transferON_OFF==ACTIVE .OR. &
     &     Dprst_transferON_OFF==ACTIVE .OR. Lake_transferON_OFF==ACTIVE .OR. Nconsumed>0 ) THEN
        IF ( Dprst_transferON_OFF==ACTIVE .AND. Dprst_flag==OFF ) THEN
          PRINT *, 'ERROR, specified water-use event based dprst input and have dprst inactive'
          Inputerror_flag = 1
        ENDIF
        IF ( Lake_transferON_OFF==ACTIVE .AND. Strmflow_flag==strmflow_muskingum_lake_module ) THEN
          PRINT *, 'ERROR, specified water-use event based lake input and have lake simulation inactive'
          Inputerror_flag = 1
        ENDIF
      ENDIF

      Stream_order_flag = 0
      IF ( Nsegment>0 .AND. Strmflow_flag>1 .AND. Model/=0 ) THEN
        Stream_order_flag = 1 ! strmflow_in_out, muskingum, muskingum_lake, muskingum_mann
      ENDIF

      IF ( Nsegment<1 .AND. Model/=DOCUMENTATION ) THEN
        IF ( Stream_order_flag==1 .OR. Call_cascade==1 ) THEN
          PRINT *, 'ERROR, streamflow and cascade routing require nsegment > 0, specified as:', Nsegment
          Inputerror_flag = 1
        ENDIF
      ENDIF

      Lake_route_flag = OFF
      IF ( Nlake>0 .AND. Strmflow_flag==3 ) Lake_route_flag = ACTIVE ! muskingum_lake

      IF ( Stream_temp_flag>0 .AND. Stream_order_flag==0 ) THEN
        PRINT *, 'ERROR, stream temperature computation requires streamflow routing, thus strmflow_module'
        PRINT *, '       must be set to strmflow_in_out, muskingum, muskingum_mann, or muskingum_lake'
        Inputerror_flag = 1
      ENDIF

      IF ( NsubOutON_OFF==1 .AND. Nsub==0 ) THEN
        NsubOutON_OFF = 0
        IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, nsubOutON_OFF = 1 and nsub = 0, thus nsub_summary not used'
      ENDIF

      IF ( NsegmentOutON_OFF==1 .AND. Nsegment==0 ) THEN
        NsegmentOutON_OFF = 0
        IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, nsegmentOutON_OFF = 1 and nsegment = 0, thus nsegment_summary not used'
      ENDIF

      IF ( Model==DOCUMENTATION .OR. Parameter_check_flag>0 ) CALL check_dimens()

      check_dims = Inputerror_flag
      END FUNCTION check_dims

!***********************************************************************
!     Check consistency of dimensions with flags
!***********************************************************************
      SUBROUTINE check_dimens()
      USE PRMS_MODULE
      IMPLICIT NONE
!***********************************************************************
      IF ( Nhru==0 .OR. Nssr==0 .OR. Ngw==0 ) THEN
        PRINT *, 'ERROR, nhru, nssr, and ngw must be > 0: nhru=', Nhru, ', nssr=', Nssr, ', ngw=', Ngw
        Inputerror_flag = 1
      ELSEIF ( Nssr/=Nhru .OR. Ngw/=Nhru ) THEN
        PRINT *, 'ERROR, nhru, nssr, and ngw must equal: nhru=', Nhru, ', nssr=', Nssr, ', ngw=', Ngw
        Inputerror_flag = 1
      ENDIF
      IF ( Ndepl==0 ) THEN
        PRINT *, 'ERROR, ndepl must be > 0: ndepl=', Ndepl
        Inputerror_flag = 1
      ENDIF
      IF ( Ndeplval/=Ndepl*11 ) THEN
        PRINT *, 'ERROR, ndeplval must be = ndepl*11: ndeplval:', Ndeplval, ', ndepl=', Ndepl
        Inputerror_flag = 1
      ENDIF

      IF ( Model==DOCUMENTATION ) THEN
        IF ( Ntemp==0 ) Ntemp = 1
        IF ( Nrain==0 ) Nrain = 1
        IF ( Nlake==0 ) Nlake = 1
        IF ( Nlake_hrus==0 ) Nlake_hrus = 1
        IF ( Nsol==0 ) Nsol = 1
        IF ( Nobs==0 ) Nobs = 1
        IF ( Ncascade==0 ) Ncascade = 1
        IF ( Ncascdgw==0 ) Ncascdgw = 1
        IF ( Nsub==0 ) Nsub = 1
        IF ( Nevap==0 ) Nevap = 1
        IF ( Nhrucell==0 ) Nhrucell = 1
        IF ( Ngwcell==0 ) Ngwcell = 1
        IF ( Nsegment==0 ) Nsegment = 1
        IF ( Nratetbl==0 ) Nratetbl = 4
        IF ( Nwateruse==0 ) Nwateruse = 1
        IF ( Nexternal==0 ) Nexternal = 1
        IF ( Nconsumed==0 ) Nconsumed = 1
        IF ( Npoigages==0 ) Npoigages = 1
        IF ( Nmap2hru==0 ) Nmap2hru = 1
        IF ( Nmap==0 ) Nmap = 1
        Subbasin_flag = ACTIVE
        Cascade_flag = CASCADE_NORMAL
        Cascadegw_flag = CASCADE_NORMAL
        Call_cascade = ACTIVE
        Stream_order_flag = ACTIVE
        Climate_hru_flag = ACTIVE
        Lake_route_flag = ACTIVE
        Water_use_flag = ACTIVE
        Segment_transferON_OFF = ACTIVE
        Gwr_transferON_OFF = ACTIVE
        External_transferON_OFF = ACTIVE
        Dprst_transferON_OFF = ACTIVE
        Lake_transferON_OFF = ACTIVE
      ENDIF

      END SUBROUTINE check_dimens

!***********************************************************************
!     Call output summary routines
!***********************************************************************
      SUBROUTINE summary_output()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF
      USE PRMS_MODULE, ONLY: NhruOutON_OFF, NsubOutON_OFF, BasinOutON_OFF, NsegmentOutON_OFF, statsON_OFF
      IMPLICIT NONE
      ! Functions
      EXTERNAL :: nhru_summary, nsub_summary, basin_summary, nsegment_summary, statvar_out
!***********************************************************************
      IF ( NhruOutON_OFF>OFF ) CALL nhru_summary()

      IF ( NsubOutON_OFF==ACTIVE ) CALL nsub_summary()

      IF ( BasinOutON_OFF==ACTIVE ) CALL basin_summary()

      IF ( NsegmentOutON_OFF>OFF ) CALL nsegment_summary()

      IF ( statsON_OFF==ACTIVE ) CALL statvar_out()

      END SUBROUTINE summary_output

!**********************************************************************
!     Module documentation
!**********************************************************************
      SUBROUTINE module_doc()
      USE PRMS_CONSTANTS, ONLY: DECL
      USE PRMS_MODULE, ONLY: Process_flag
      use PRMS_SET_TIME, only: prms_time
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: basin, climateflow
      INTEGER, EXTERNAL :: cascade, obs, soltab, transp_tindex
      INTEGER, EXTERNAL :: transp_frost, frost_date, routing
      INTEGER, EXTERNAL :: temp_1sta_laps, temp_dist2
      INTEGER, EXTERNAL :: precip_1sta_laps, climate_hru
      INTEGER, EXTERNAL :: precip_dist2, xyz_dist, ide_dist
      INTEGER, EXTERNAL :: ddsolrad, ccsolrad
      INTEGER, EXTERNAL :: potet_pan, potet_jh, potet_hamon, potet_hs, potet_pt, potet_pm
      INTEGER, EXTERNAL :: intcp, snowcomp, gwflow, srunoff, soilzone_ag
      INTEGER, EXTERNAL :: strmflow, subbasin, basin_sum, map_results, strmflow_in_out
      INTEGER, EXTERNAL :: write_climate_hru, muskingum, muskingum_lake
      INTEGER, EXTERNAL :: stream_temp
      EXTERNAL :: nhru_summary, prms_summary, water_balance, nsub_summary, basin_summary, nsegment_summary
      INTEGER, EXTERNAL :: dynamic_param_read, water_use_read, potet_pm_sta, glacr, setup
      EXTERNAL :: precip_map, temp_map
! Local variable
      INTEGER :: test
!**********************************************************************
      test = basin()
      test = cascade()
      test = climateflow()
      test = soltab()
      test = setup()
      test = prms_time()
      test = obs()
      test = water_use_read()
      test = dynamic_param_read()
      test = temp_1sta_laps()
      test = temp_dist2()
      test = xyz_dist()
      test = ide_dist()
      CALL temp_map()
      CALL precip_map()
      test = climate_hru()
      test = precip_1sta_laps()
      test = precip_dist2()
      test = ddsolrad()
      test = ccsolrad()
      test = transp_tindex()
      test = frost_date()
      test = transp_frost()
      test = potet_jh()
      test = potet_hamon()
      test = potet_pan()
      test = potet_hs()
      test = potet_pt()
      test = potet_pm()
      test = potet_pm_sta()
      test = write_climate_hru()
      test = intcp()
      test = snowcomp()
      test = srunoff()
      test = glacr()
      test = soilzone_ag()
      test = gwflow()
      test = routing()
      test = strmflow()
      test = strmflow_in_out()
      test = muskingum()
      test = muskingum_lake()
      test = stream_temp()
      test = basin_sum()
      test = map_results()
      CALL nhru_summary()
      CALL nsub_summary()
      CALL basin_summary()
      CALL nsegment_summary()
      CALL prms_summary()
      CALL water_balance()
      test = subbasin()

      IF ( Process_flag==DECL ) PRINT 9001

 9001 FORMAT (//, ' All available modules have been called.', /, &
     &        ' All parameters have been declared.', /, &
     &        ' Note, no simulation was computed.', /)

      END SUBROUTINE module_doc

!***********************************************************************
!     check module names
!***********************************************************************
      SUBROUTINE check_module_names()
      USE PRMS_CONSTANTS, ONLY: ERROR_control
      USE PRMS_MODULE, ONLY: Temp_module, Precip_module, Et_module, Solrad_module, &
     &    Transp_module, Srunoff_module, Strmflow_module
      IMPLICIT NONE
! Local Variables
      INTEGER :: ierr
!***********************************************************************
      ierr = 0
      IF ( Temp_module(:14)=='temp_1sta_prms' ) THEN
        PRINT *, 'WARNING, deprecated temp_module value, change temp_1sta_prms to temp_1sta'
        Temp_module = 'temp_1sta'
      ELSEIF ( Temp_module(:14)=='temp_laps_prms' ) THEN
        PRINT *, 'WARNING, deprecated temp_module value, change temp_laps_prms to temp_laps'
        Temp_module = 'temp_laps'
      ELSEIF ( Temp_module(:15)=='temp_dist2_prms' ) THEN
        PRINT *, 'WARNING, deprecated temp_module value, change temp_dist2_prms to temp_dist2'
        Temp_module = 'temp_dist2'
      ELSEIF ( Temp_module(:9)=='temp_2sta' ) THEN
        PRINT *, 'ERROR, module temp_2sta_prms not available, use a different temp_module'
        ierr = 1
      ENDIF

      IF ( Precip_module(:11)=='precip_prms' ) THEN
        PRINT *, 'WARNING, deprecated precip_module value, change precip_prms to precip_1sta'
        Precip_module = 'precip_1sta'
      ELSEIF ( Precip_module(:16)=='precip_laps_prms' ) THEN
        PRINT *, 'WARNING, deprecated precip_module value, change precip_laps_prms to precip_laps'
        Precip_module = 'precip_laps'
      ELSEIF ( Precip_module(:17)=='precip_dist2_prms' ) THEN
        PRINT *, 'WARNING, deprecated precip_module value, change precip_dist2_prms to precip_dist2'
        Precip_module = 'precip_dist2'
      ENDIF

      IF ( Temp_module(:8)=='ide_dist' .AND. Precip_module(:8)/='ide_dist') THEN
        PRINT '(/,A,/,2A)', 'ERROR, if ide_dist is specified for temp_module,', &
     &        'it also must be specified for precip_module: ', Precip_module
        ierr = 1
      ELSEIF ( Precip_module(:8)=='ide_dist' .AND. Temp_module(:8)/='ide_dist') THEN
        PRINT '(/,A,/,2A)', 'ERROR, if ide_dist is specified for precip_module,', &
     &        'it also must be specified for temp_module: ', Temp_module
        ierr = 1
      ELSEIF ( Temp_module(:8)=='xyz_dist' .AND. Precip_module(:8)/='xyz_dist') THEN
        PRINT '(/,A,/,2A)', 'ERROR, if xyz_dist is specified for temp_module,', &
     &        'it also must be specified for precip_module: ', Precip_module
        ierr = 1
      ELSEIF ( Precip_module(:8)=='xyz_dist' .AND. Temp_module(:8)/='xyz_dist') THEN
        PRINT '(/,A,/,2A)', 'ERROR, if xyz_dist is specified for precip_module,', &
     &        'it also must be specified for temp_module: ', Temp_module
        ierr = 1
      ENDIF

      IF ( Transp_module(:18)=='transp_tindex_prms' ) THEN
        PRINT *, 'WARNING, deprecated transp_module value, change transp_tindex_prms to transp_tindex'
        Transp_module = 'transp_tindex'
      ENDIF

      IF ( Et_module(:13)=='potet_jh_prms' ) THEN
        PRINT *, 'WARNING, deprecated et_module value, change potet_jh_prms to potet_jh'
        Et_module = 'potet_jh'
      ELSEIF ( Et_module(:14)=='potet_pan_prms' ) THEN
        PRINT *, 'WARNING, deprecated et_module value, change potet_pan_prms to potet_pan'
        Et_module = 'potet_pan'
      ELSEIF ( Et_module(:15)=='potet_epan_prms' ) THEN
        PRINT *, 'ERROR, deprecated et_module value, change potet_epan_prms to potet_pan'
        ierr = 1
      ELSEIF ( Et_module(:20)=='potet_hamon_hru_prms' ) THEN
        PRINT *, 'WARNING, deprecated et_module value, change potet_hamon_hru_prms to potet_hamon_hru'
        Et_module = 'potet_hamon'
      ELSEIF ( Et_module(:16)=='potet_hamon_prms' ) THEN
        PRINT *, 'WARNING, deprecated et_module value, change potet_hamon_prms to potet_hamon'
        Et_module = 'potet_hamon'
      ENDIF

      IF ( Solrad_module(:17)=='ddsolrad_hru_prms' ) THEN
        PRINT *, 'WARNING, deprecated solrad_module value, change ddsolrad_hru_prms to ddsolrad'
        Solrad_module = 'ddsolrad'
      ELSEIF ( Solrad_module(:17)=='ccsolrad_hru_prms' ) THEN
        PRINT *, 'WARNING, deprecated solrad_module value, change ccsolrad_hru_prms to ccsolrad'
        Solrad_module = 'ccsolrad'
      ELSEIF ( Solrad_module(:13)=='ddsolrad_prms' ) THEN
        PRINT *, 'WARNING, deprecated solrad_module value, change ddsolrad_prms to ddsolrad'
        Solrad_module = 'ddsolrad'
      ELSEIF ( Solrad_module(:13)=='ccsolrad_prms' ) THEN
        PRINT *, 'WARNING, deprecated solrad_module value, change ccsolrad_prms to ccsolrad'
        Solrad_module = 'ccsolrad'
      ENDIF

      IF ( Srunoff_module(:18)=='srunoff_carea_prms' ) THEN
        PRINT *, 'WARNING, deprecated srunoff_module value, change srunoff_carea_prms to srunoff_carea'
        Srunoff_module = 'srunoff_carea'
      ELSEIF ( Srunoff_module(:18)=='srunoff_smidx_prms' ) THEN
        PRINT *, 'WARNING, deprecated srunoff_module value, change srunoff_smidx_prms to srunoff_smidx'
        Srunoff_module = 'srunoff_smidx'
      ENDIF

      IF ( Strmflow_module(:13)=='strmflow_prms' ) THEN
        PRINT *, 'WARNING, deprecated strmflow_module value, change strmflow_prms to strmflow'
        Strmflow_module = 'strmflow'
      ELSEIF ( Strmflow_module(:13)=='strmflow_lake' ) THEN
        PRINT *, 'ERROR, module strmflow_lake not available, use a different strmflow_module, such as muskingum_lake'
        ierr = 1
      ENDIF
      IF ( ierr==1 ) ERROR STOP ERROR_control

      END SUBROUTINE check_module_names

!***********************************************************************
!     call_modules_restart - write or read restart file
!***********************************************************************
      SUBROUTINE call_modules_restart(In_out)
      USE PRMS_MODULE
      use prms_utils, only: check_restart, check_restart_dimen
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      INTRINSIC :: TRIM
      ! Local Variables
      INTEGER :: nhru_test, dprst_test, nsegment_test, temp_test, et_test, ierr, time_step
      INTEGER :: cascade_test, cascdgw_test, nhrucell_test, nlake_test, transp_test, start_time(6), end_time(6)
      CHARACTER(LEN=MAXCONTROL_LENGTH) :: model_test
      CHARACTER(LEN=12) :: module_name
!***********************************************************************
      IF ( In_out==0 ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Timestep, Nhru, Dprst_flag, Nsegment, Temp_flag, Et_flag, &
     &          Cascade_flag, Cascadegw_flag, Nhrucell, Nlake, Transp_flag, Model_mode
        WRITE ( Restart_outunit ) Starttime, Endtime
      ELSE
        ierr = 0
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) time_step, nhru_test, dprst_test, nsegment_test, temp_test, et_test, &
     &         cascade_test, cascdgw_test, nhrucell_test, nlake_test, transp_test, model_test
        READ ( Restart_inunit ) start_time, end_time
        IF ( Print_debug>DEBUG_minimum ) PRINT 4, EQULS, 'Simulation time period of Restart File:', &
     &       start_time(1), start_time(2), start_time(3), ' -', end_time(1), end_time(2), end_time(3), &
     &       'Last time step of simulation: ', time_step, EQULS
    4   FORMAT (/, A, /, 2(A, I5, 2('/',I2.2)), /, A, I0, /, A, /)
        IF ( TRIM(Model_mode)/=TRIM(model_test) ) THEN
          PRINT *, 'ERROR, Initial Conditions File saved for model_mode=', model_test
          PRINT *, '       Current model has model_mode=', Model_mode, ' they must be equal'
          ierr = 1
        ENDIF
        CALL check_restart_dimen('nhru', nhru_test, Nhru, ierr)
        CALL check_restart_dimen('nhrucell', nhrucell_test, Nhrucell, ierr)
        CALL check_restart_dimen('nlake', nlake_test, Nlake, ierr)
        IF ( Dprst_flag/=dprst_test ) THEN
          PRINT *, 'ERROR, Initial Conditions File saved for model with dprst_flag=', dprst_test
          PRINT *, '       Current model has dprst_flag=', Dprst_flag, ' they must be equal'
          ierr = 1
        ENDIF
        IF ( Cascade_flag/=cascade_test ) THEN
          PRINT *, 'ERROR, Initial Conditions File saved for model with cascade_flag=', cascade_test
          PRINT *, '       Current model has cascade_flag=', Cascade_flag, ' they must be equal'
          ierr = 1
        ENDIF
        IF ( Cascadegw_flag/=cascdgw_test ) THEN
          PRINT *, 'ERROR, Initial Conditions File saved for model with cascadegw_flag=', cascdgw_test
          PRINT *, '       Current model has cascadegw_flag=', Cascadegw_flag, ' they must be equal'
          ierr = 1
        ENDIF
        CALL check_restart_dimen('nsegment', nsegment_test, Nsegment, ierr)
        ! Temp_flag (1=temp_1sta; 2=temp_laps; 3=temp_dist2; 5=ide_dist; 6=xyz_dist; 7=climate_hru; 8=temp_sta
        IF ( Temp_flag/=temp_test ) THEN
          IF ( Temp_flag<4 .OR. temp_test<4 ) THEN
            PRINT *, 'ERROR, Initial Conditions File saved for model with different temperature module'
            PRINT *, '       than current model, cannot switch to/from temp_1sta, temp_laps, or temp_dist2'
            ierr = 1
          ENDIF
        ENDIF
        IF ( Et_flag/=et_test ) THEN
          IF ( Et_flag==4 .OR. et_test==4 ) THEN
            PRINT *, 'ERROR, Cannot switch to/from potet_pan module for restart simulations'
            ierr = 1
          ENDIF
        ENDIF
        IF ( Transp_flag/=transp_test ) THEN
          IF ( Transp_flag==1 .OR. transp_test==1 ) THEN
            PRINT *, 'ERROR, Cannot switch to/from transp_tindex module for restart simulations'
            ierr = 1
          ENDIF
        ENDIF
        IF ( ierr==1 ) ERROR STOP ERROR_restart
      ENDIF
      END SUBROUTINE call_modules_restart
