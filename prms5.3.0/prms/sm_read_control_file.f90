!***********************************************************************
! Read Control File
!***********************************************************************
submodule(PRMS_CONTROL_FILE) sm_read_control_file
  use PRMS_CONSTANTS, only: MAXCONTROL_LENGTH, MAXFILE_LENGTH, INT_TYPE, REAL_TYPE, CHAR_TYPE, ERROR_control
  use PRMS_MODULE
  use PRMS_MAP_RESULTS, only: NmapOutVars, MapOutVar_names
  use PRMS_STATVAR_OUT, only: statvarOut_format, nstatVars, statVar_element, statVar_names
  use PRMS_NHRU_SUMMARY, only: NhruOutVars, NhruOut_freq, NhruOutBaseFileName, NhruOutVar_names, NhruOut_format, NhruOutNcol
  use PRMS_NSUB_SUMMARY, only: NsubOutVars, NsubOut_freq, NsubOutBaseFileName, NsubOutVar_names, NsubOut_format
  use PRMS_BASIN_SUMMARY, only: BasinOutVars, BasinOut_freq, BasinOutBaseFileName, BasinOutVar_names
  use PRMS_NSEGMENT_SUMMARY, only: NsegmentOutVars, NsegmentOut_freq, NsegmentOutBaseFileName, &
                                   NsegmentOutVar_names, NsegmentOut_format
  use PRMS_WATER_USE, only: Segment_transfer_file, Gwr_transfer_file, Dprst_transfer_file, &
                            External_transfer_file, Lake_transfer_file
  use PRMS_GLACR, only: Mbinit_flag
  use PRMS_PRECIP_MAP, only: Precip_map_file
  use PRMS_TEMP_MAP, only: Tmax_map_file, Tmin_map_file

contains

  module subroutine read_control_file()
    use PRMS_CONSTANTS, only: MAXCONTROL_LENGTH, MAXFILE_LENGTH
    use PRMS_MODULE, only:  Model_control_file
    use prms_utils, only: numchars, PRMS_open_input_file, PRMS_open_output_file, read_error, write_outfile
    implicit none
    ! Functions
    intrinsic TRIM
    ! Local Variables
    character(LEN=MAXCONTROL_LENGTH) :: paramname
    character(LEN=4) :: string
    integer ios, numvalues, param_type, control_unit, j
    integer, allocatable :: int_parameter_values(:)
    character(LEN=MAXFILE_LENGTH), allocatable :: parameter_values(:)
    character(LEN=MAXCONTROL_LENGTH) :: paramstring
    real, allocatable :: real_parameter_values(:)
    !***********************************************************************

    ! control filename cannot include blanks
    call get_control_filename()
    call PRMS_open_input_file(control_unit, Model_control_file, 'model_control_file', 0, ios)
    if (ios /= 0) call read_error(10, trim(Model_control_file))
    ! read header
    read (control_unit, '(A)', IOSTAT=ios) Control_description
    if (ios /= 0) call read_error(12, trim(Model_control_file))

    call setup_cont() ! set default control parameter values

    ! Read all Control Parameters
    do
      read (control_unit, '(A)', IOSTAT=ios) string
      if (ios == -1) exit ! found end of Control File
      if (ios /= 0) call read_error(12, 'missing #### delimiter')
      if (string(:4) /= '####') cycle ! skip until delimiter found, such as blank of // comment lines
      read (control_unit, '(A)', IOSTAT=ios) paramname ! parameter name
      if (ios /= 0) call read_error(5, 'missing parameter name')
      read (control_unit, *, IOSTAT=ios) numvalues
      if (ios /= 0) call read_error(5, 'invalid number of values: '//trim(paramname))
      read (control_unit, *, IOSTAT=ios) param_type
      if (ios /= 0) call read_error(5, 'invalid parameter type: '//trim(paramstring))
      if (param_type<1 .or. param_type>4 .or. param_type==3) call read_error(5, 'invalid parameter type: '//trim(paramstring))
      allocate (int_parameter_values(numvalues), real_parameter_values(numvalues), parameter_values(numvalues))
      if (param_type == INT_TYPE) then
        read (Control_unit, *, IOSTAT=ios) (int_parameter_values(j), j=1, numvalues)
        if (ios /= 0) call read_error(5, 'invalid integer value: '//trim(paramname))
      elseif (param_type == CHAR_TYPE) then
        do j = 1, numvalues
          read (Control_unit, '(A)', IOSTAT=ios) parameter_values(j)
          if (ios /= 0) call read_error(5, 'invalid character value: '//trim(paramname))
        end do
      else
        read (Control_unit, *, IOSTAT=ios) (real_parameter_values(j), j=1, numvalues)
        if (ios /= 0) call read_error(5, 'invalid real value: '//trim(paramname))
      end if
      call set_control_parameter(paramname, numvalues, int_parameter_values, real_parameter_values, parameter_values)
      deallocate (int_parameter_values, real_parameter_values, parameter_values)
    end do
    ! reset control parameters based on command line
    close (control_unit)

  end subroutine read_control_file

  !***********************************************************************
  ! setup_cont - Set control parameter value defaults
  !***********************************************************************
  module subroutine setup_cont()
    use PRMS_CONSTANTS, only: DEBUG_normal, ACTIVE, OFF
    use PRMS_MODULE, only: Print_debug, Dprst_flag, Cascade_flag, Soilzone_aet_flag, &
                           Albedo_cbh_flag, Cloud_cover_cbh_flag, Csv_output_file, irrigated_area_module, AET_module, &
                           PET_ag_module, selectDatesFileName, outputSelectDatesON_OFF, snow_cloudcover_flag, &
                           Dyn_ag_frac_flag, Dyn_ag_soil_flag, AET_cbh_flag, PET_cbh_flag, Dprst_add_water_use, &
                           Dprst_transfer_water_use, Iter_aet_flag
    use PRMS_CLIMATE_HRU, only: Precip_day, Tmax_day, Tmin_day, Potet_day, Transp_day, Swrad_day, Albedo_day, Cloud_cover_day, &
                                Cbh_check_flag, Cbh_binary_flag, Windspeed_day, Humidity_day, &
                                AET_cbh_file, PET_cbh_file, irrigated_area_cbh_file
    use PRMS_DYNAMIC_PARAM_READ, only: imperv_frac_dynamic, imperv_stor_dynamic, dprst_depth_dynamic, dprst_frac_dynamic, &
                                       wrain_intcp_dynamic, srain_intcp_dynamic, snow_intcp_dynamic, covtype_dynamic, &
                                       potetcoef_dynamic, transpbeg_dynamic, transpend_dynamic, Dynamic_param_log_file, &
                                       soilmoist_dynamic, soilrechr_dynamic, radtrncf_dynamic, &
                                       fallfrost_dynamic, springfrost_dynamic, transp_on_dynamic, &
                                       covden_sum_dynamic, covden_win_dynamic, sro2dprst_perv_dyn, sro2dprst_imperv_dyn, &
                                       ag_soilmoist_dynamic, ag_soilrechr_dynamic, ag_frac_dynamic !, snareathresh_dynamic
    implicit none
    ! Local Variables
    integer i, numvalues
    !***********************************************************************
    Num_control_parameters = Max_num_control_parameters
    ! allocate and store parameter data
    allocate (Control_parameter_data(Num_control_parameters))
    do i = 1, Num_control_parameters
      Control_parameter_data(i) % read_flag = 2 ! 2 = set to default; 1 = read from Control File
      Control_parameter_data(i) % data_type = INT_TYPE ! 1 = integer, 2 = real, 4 = string
      Control_parameter_data(i) % allocate_flag = 0 ! set to 1 if allocatable
      Control_parameter_data(i) % numvals = 1
      Control_parameter_data(i) % name = ' '
      ! WARNING, parameter index is set based on order defaults defined
      Control_parameter_data(i) % index = i
      allocate (Control_parameter_data(i) % values_int(1))
      allocate (Control_parameter_data(i) % values_real(1))
      allocate (Control_parameter_data(i) % values_character(1))
      Control_parameter_data(i) % values_int(1) = 0
      Control_parameter_data(i) % values_real(1) = 0.0
      Control_parameter_data(i) % values_character(1) = ' '
    end do

    !      DO i = Num_control_parameters+1, Num_control_parameters+20
    !        Control_parameter_data(i)%read_flag = 0 ! 0 means not set
    !        Control_parameter_data(i)%data_type = 0 ! 1 = integer, 2 = real, 4 = string
    !        Control_parameter_data(i)%numvals = 0
    !        Control_parameter_data(i)%name = ' '
    ! WARNING, parameter index is set based on order defaults defined
    !        Control_parameter_data(i)%index = i
    !      ENDDO

    ! assign default value for integer flags
    ! note: default value for all parameters set to 0, only need to reset if other than 0
    numvalues = 1
    i = 1
    Control_parameter_data(i) % name = 'print_debug'
    Print_debug = DEBUG_normal
    i = i + 1
    Control_parameter_data(i) % name = 'parameter_check_flag'
    Parameter_check_flag = ACTIVE
    Control_parameter_data(i) % values_int(1) = Parameter_check_flag
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_flag'
    Dprst_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_add_water_use'
    Dprst_add_water_use = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_transfer_water_use'
    Dprst_transfer_water_use = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'cascade_flag'
    Cascade_flag = ACTIVE
    Control_parameter_data(i) % values_int(1) = Cascade_flag
    i = i + 1
    Control_parameter_data(i) % name = 'cascadegw_flag'
    Cascadegw_flag = ACTIVE
    Control_parameter_data(i) % values_int(1) = Cascadegw_flag
    i = i + 1
    Control_parameter_data(i) % name = 'save_vars_to_file'
    Save_vars_to_file = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'init_vars_from_file'
    Init_vars_from_file = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'frozen_flag'
    Frozen_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'glacier_flag'
    Glacier_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'no_snow_flag'
    no_snow_flag = ACTIVE
    i = i + 1
    Control_parameter_data(i) % name = 'stream_temp_flag'
    Stream_temp_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'strmtemp_humidity_flag'
    Strmtemp_humidity_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'stream_temp_shade_flag'
    Stream_temp_shade_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'snarea_curve_flag'
    Snarea_curve_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'snow_cloudcover_flag'
    snow_cloudcover_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'orad_flag'
    Orad_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'subbasin_flag'
    Subbasin_flag = ACTIVE
    Control_parameter_data(i) % values_int(1) = Subbasin_flag
    i = i + 1
    Control_parameter_data(i) % name = 'cbh_check_flag'
    Cbh_check_flag = ACTIVE
    Control_parameter_data(i) % values_int(1) = Cbh_check_flag
    i = i + 1
    Control_parameter_data(i) % name = 'cbh_binary_flag'
    Cbh_binary_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'gwr_swale_flag'
    Gwr_swale_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'snow_cbh_flag'
    Snow_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'gwflow_cbh_flag'
    Gwflow_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'humidity_cbh_flag'
    Humidity_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'windspeed_cbh_flag'
    Windspeed_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'mbinit_flag'
    Mbinit_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'outputSelectDatesON_OFF'
    outputSelectDatesON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOutON_OFF'
    NhruOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOut_freq'
    NhruOut_freq = 1
    Control_parameter_data(i) % values_int(1) = NhruOut_freq
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOutNcol'
    NhruOutNcol = 0
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOut_format'
    NhruOut_format = 1
    Control_parameter_data(i) % values_int(1) = NhruOut_format
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOutVars'
    NhruOutVars = 0
    i = i + 1
    Control_parameter_data(i) % name = 'statvarOut_format'
    statvarOut_format = 1
    Control_parameter_data(i) % values_int(1) = statvarOut_format
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOutON_OFF'
    NsubOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOutVars'
    NsubOutVars = 0
    i = i + 1
    Control_parameter_data(i) % name = 'basinOutON_OFF'
    BasinOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'basinOutVars'
    BasinOutVars = 0
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOutON_OFF'
    NsegmentOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOutVars'
    NsegmentOutVars = 0
    i = i + 1
    Control_parameter_data(i) % name = 'statsON_OFF'
    StatsON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'csvON_OFF'
    CsvON_OFF = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'aniOutON_OFF'
!    AniOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'mapOutON_OFF'
    MapOutON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'nstatVars'
    nstatVars = 0
    i = i + 1
    Control_parameter_data(i) % name = 'nmapOutVars'
    NmapOutVars = 0
!    i = i + 1
!    Control_parameter_data(i) % name = 'naniOutVars'
!    NaniOutVars = 0
!    i = i + 1
!    Control_parameter_data(i) % name = 'ndispGraphs'
!    NdispGraphs = 0
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOut_freq'
    NsubOut_freq = 1
    Control_parameter_data(i) % values_int(1) = NsubOut_freq
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOut_format'
    NsubOut_format = 1
    Control_parameter_data(i) % values_int(1) = NsubOut_format
    i = i + 1
    Control_parameter_data(i) % name = 'basinOut_freq'
    BasinOut_freq = 1
    Control_parameter_data(i) % values_int(1) = BasinOut_freq
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOut_freq'
    NsegmentOut_freq = 1
    Control_parameter_data(i) % values_int(1) = NsegmentOut_freq
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOut_format'
    NsegmentOut_format = 1
    Control_parameter_data(i) % values_int(1) = NsegmentOut_format
    i = i + 1
    Control_parameter_data(i) % name = 'prms_warmup'
    Prms_warmup = 0
    Control_parameter_data(i) % values_int(1) = Prms_warmup
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_imperv_flag'
    Dyn_imperv_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_intcp_flag'
    Dyn_intcp_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_covden_flag'
    Dyn_covden_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'soilzone_aet_flag'
    Soilzone_aet_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'iter_aet_flag'
    Iter_aet_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'albedo_cbh_flag'
    Albedo_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'cloud_cover_cbh_flag'
    Cloud_cover_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_ag_frac_flag'
    Dyn_ag_frac_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_ag_soil_flag'
    Dyn_ag_soil_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'AET_cbh_flag'
    AET_cbh_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'PET_cbh_flag'
    PET_cbh_flag = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'dispGraphsBuffSize'
!    DispGraphsBuffSize = 50
!    Control_parameter_data(i) % values_int(1) = DispGraphsBuffSize
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_transp_flag'
    Dyn_transp_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_transp_on_flag'
    Dyn_transp_on_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_sro2dprst_perv_flag'
    Dyn_sro2dprst_perv_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_sro2dprst_imperv_flag'
    Dyn_sro2dprst_imperv_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_covtype_flag'
    Dyn_covtype_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_fallfrost_flag'
    Dyn_fallfrost_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_springfrost_flag'
    Dyn_springfrost_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_potet_flag'
    Dyn_potet_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_soil_flag'
    Dyn_soil_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_radtrncf_flag'
    Dyn_radtrncf_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_snareathresh_flag'
    Dyn_snareathresh_flag = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'dyn_sro_to_dprst_flag'
!    Dyn_sro_to_dprst_flag = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'dyn_sro_to_imperv_flag'
!    Dyn_sro_to_imperv_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dyn_dprst_flag'
    Dyn_dprst_flag = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'segment_transferON_OFF'
    Segment_transferON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'gwr_transferON_OFF'
    Gwr_transferON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'external_transferON_OFF'
    External_transferON_OFF = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'consumed_transferON_OFF'
!    Consumed_transferON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'lake_transferON_OFF'
    Lake_transferON_OFF = OFF
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_transferON_OFF'
    Dprst_transferON_OFF = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'soilzone_transferON_OFF'
!    Soilzone_transferON_OFF = OFF
!    i = i + 1
!    Control_parameter_data(i) % name = 'canopy_transferON_OFF'
!    Canopy_transferON_OFF = OFF
    i = i + 1

    ! parameters that get allocated if in Control File
    Control_parameter_data(i) % name = 'mapOutVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'statVar_element'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'statVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
!    i = i + 1
!    Control_parameter_data(i) % name = 'aniOutVar_names'
!    Control_parameter_data(i) % data_type = CHAR_TYPE
!    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
!    i = i + 1
!    Control_parameter_data(i) % name = 'dispVar_names'
!    Control_parameter_data(i) % data_type = CHAR_TYPE
!    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
!    i = i + 1
!    Control_parameter_data(i) % name = 'dispVar_plot'
!    Control_parameter_data(i) % data_type = CHAR_TYPE
!    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
!    i = i + 1
!    Control_parameter_data(i) % name = 'dispVar_element'
!    Control_parameter_data(i) % data_type = CHAR_TYPE
!    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOutVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'basinOutVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOutVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOutVar_names'
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    i = i + 1

    ! floating point parameters
 !   Control_parameter_data(i) % name = 'initial_deltat'
 !   Initial_deltat = 24.0
 !   Control_parameter_data(i) % values_real(1) = Initial_deltat
 !   Control_parameter_data(i) % data_type = REAL_TYPE
 !   i = i + 1

    ! assign default value for character parameters
    Control_parameter_data(i) % name = 'selectDatesFileName'
    selectDatesFileName = 'dates'
    Control_parameter_data(i) % values_character(1) = selectDatesFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'model_mode'
    Model_mode = 'GSFLOW5'
    Control_parameter_data(i) % values_character(1) = Model_mode
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'executable_desc'
    Executable_desc = 'PRMS 5.3'
    Control_parameter_data(i) % values_character(1) = Executable_desc
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'executable_model'
    Executable_model = 'gsflow'
    Control_parameter_data(i) % values_character(1) = Executable_model
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'precip_module'
    Precip_module = 'precip_1sta'
    Control_parameter_data(i) % values_character(1) = Precip_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'temp_module'
    Temp_module = 'temp_1sta'
    Control_parameter_data(i) % values_character(1) = Temp_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'solrad_module'
    Solrad_module = 'ddsolrad'
    Control_parameter_data(i) % values_character(1) = Solrad_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'et_module'
    Et_module = 'potet_jh'
    Control_parameter_data(i) % values_character(1) = Et_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'srunoff_module'
    Srunoff_module = 'srunoff_smidx'
    Control_parameter_data(i) % values_character(1) = Srunoff_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'strmflow_module'
    Strmflow_module = 'strmflow'
    Control_parameter_data(i) % values_character(1) = Strmflow_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'transp_module'
    Transp_module = 'transp_tindex'
    Control_parameter_data(i) % values_character(1) = Transp_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'soilzone_module'
    Soilzone_module = 'soilzone'
    Control_parameter_data(i) % values_character(1) = Soilzone_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'data_file'
    Data_file = 'prms.data'
    Control_parameter_data(i) % values_character(1) = Data_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'param_file'
    Param_file = 'prms.params'
    Control_parameter_data(i) % values_character(1) = Param_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    Control_parameter_data(i) % allocate_flag = 1 ! need to allocate
    Param_file_control_parameter_id = i
    i = i + 1
    Control_parameter_data(i) % name = 'model_output_file'
    Model_output_file = 'prms.out'
    Control_parameter_data(i) % values_character(1) = Model_output_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'csv_output_file'
    Csv_output_file = 'prms_summary.csv'
    Control_parameter_data(i) % values_character(1) = Csv_output_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'var_save_file'
    Var_save_file = 'prms_ic.out'
    Control_parameter_data(i) % values_character(1) = Var_save_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'var_init_file'
    Var_init_file = 'prms_ic.in'
    Control_parameter_data(i) % values_character(1) = Var_init_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'stat_var_file'
    Stat_var_file = 'statvar.out'
    Control_parameter_data(i) % values_character(1) = Stat_var_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'segment_transfer_file'
    Segment_transfer_file = 'segment_transfer.in'
    Control_parameter_data(i) % values_character(1) = Segment_transfer_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'gwr_transfer_file'
    Gwr_transfer_file = 'gwr_transfer.in'
    Control_parameter_data(i) % values_character(1) = Gwr_transfer_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_transfer_file'
    Dprst_transfer_file = 'dprst_transfer.in'
    Control_parameter_data(i) % values_character(1) = Dprst_transfer_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'external_transfer_file'
    External_transfer_file = 'external_transfer.in'
    Control_parameter_data(i) % values_character(1) = External_transfer_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'lake_transfer_file'
    Lake_transfer_file = 'lake_transfer.in'
    Control_parameter_data(i) % values_character(1) = Lake_transfer_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'nhruOutBaseFileName'
    NhruOutBaseFileName = 'nhruout_path'
    Control_parameter_data(i) % values_character(1) = NhruOutBaseFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'nsubOutBaseFileName'
    NsubOutBaseFileName = 'nsubout_path'
    Control_parameter_data(i) % values_character(1) = NsubOutBaseFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'basinOutBaseFileName'
    BasinOutBaseFileName = 'basinout_path'
    Control_parameter_data(i) % values_character(1) = BasinOutBaseFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'nsegmentOutBaseFileName'
    NsubOutBaseFileName = 'nsegmentout_path'
    Control_parameter_data(i) % values_character(1) = NsegmentOutBaseFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'tmax_day'
    Tmax_day = 'tmax_day'
    Control_parameter_data(i) % values_character(1) = Tmax_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'tmin_day'
    Tmin_day = 'tmin_day'
    Control_parameter_data(i) % values_character(1) = Tmin_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'precip_day'
    Precip_day = 'precip_day'
    Control_parameter_data(i) % values_character(1) = Precip_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'swrad_day'
    Swrad_day = 'swrad_day'
    Control_parameter_data(i) % values_character(1) = Swrad_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'potet_day'
    Potet_day = 'potet_day'
    Control_parameter_data(i) % values_character(1) = Potet_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'transp_day'
    Transp_day = 'transp_day'
    Control_parameter_data(i) % values_character(1) = Transp_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'windspeed_day'
    Windspeed_day = 'windspeed_day'
    Control_parameter_data(i) % values_character(1) = Windspeed_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'humidity_day'
    Humidity_day = 'humidity_day'
    Control_parameter_data(i) % values_character(1) = Humidity_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'albedo_day'
    Albedo_day = 'albedo.day'
    Control_parameter_data(i) % values_character(1) = Albedo_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'cloud_cover_day'
    Cloud_cover_day = 'Cloud_cover_day'
    Control_parameter_data(i) % values_character(1) = Cloud_cover_day
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'irrigated_area_module'
    irrigated_area_module = 'irrigated_area_module'
    Control_parameter_data(i) % values_character(1) = irrigated_area_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'PET_ag_module'
    PET_ag_module = 'PET_ag_module'
    Control_parameter_data(i) % values_character(1) = PET_ag_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'AET_module'
    PET_ag_module = 'AET_module'
    Control_parameter_data(i) % values_character(1) = AET_module
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'AET_cbh_file'
    AET_cbh_file = 'AET_cbh_day'
    Control_parameter_data(i) % values_character(1) = AET_cbh_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'PET_cbh_file'
    PET_cbh_file = 'PET_cbh_day'
    Control_parameter_data(i) % values_character(1) = PET_cbh_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'irrigated_area_cbh_file'
    irrigated_area_cbh_file = 'irrigated_area_cbh_day'
    Control_parameter_data(i) % values_character(1) = irrigated_area_cbh_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'ag_soilmoist_dynamic'
    ag_soilmoist_dynamic = 'ag_soilmoist.dynamic'
    Control_parameter_data(i) % values_character(1) = ag_soilmoist_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'ag_soilrechr_dynamic'
    ag_soilrechr_dynamic = 'ag_soilrechr.dynamic'
    Control_parameter_data(i) % values_character(1) = ag_soilrechr_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'dynamic_param_log_file'
    Dynamic_param_log_file = 'dynamic_parameter.out'
    Control_parameter_data(i) % values_character(1) = Dynamic_param_log_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'ag_frac_dynamic'
    Dprst_depth_dynamic = 'dynag_frac'
    Control_parameter_data(i) % values_character(1) = ag_frac_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_depth_dynamic'
    Dprst_depth_dynamic = 'dyndprst_depth'
    Control_parameter_data(i) % values_character(1) = Dprst_depth_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'dprst_frac_dynamic'
    Dprst_frac_dynamic = 'dyndprst_frac'
    Control_parameter_data(i) % values_character(1) = Dprst_frac_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'snow_intcp_dynamic'
    Snow_intcp_dynamic = 'dynsnowintcp'
    Control_parameter_data(i) % values_character(1) = Snow_intcp_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'srain_intcp_dynamic'
    Srain_intcp_dynamic = 'dynsrainintcp'
    Control_parameter_data(i) % values_character(1) = Srain_intcp_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'wrain_intcp_dynamic'
    Wrain_intcp_dynamic = 'dynwrainintcp'
    Control_parameter_data(i) % values_character(1) = Wrain_intcp_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'imperv_frac_dynamic'
    Imperv_frac_dynamic = 'dynimperv_frac'
    Control_parameter_data(i) % values_character(1) = Imperv_frac_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'imperv_stor_dynamic'
    Imperv_stor_dynamic = 'dynimperv_stor'
    Control_parameter_data(i) % values_character(1) = Imperv_stor_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'covtype_dynamic'
    Covtype_dynamic = 'dyncovtype'
    Control_parameter_data(i) % values_character(1) = Covtype_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'covden_sum_dynamic'
    Covden_sum_dynamic = 'dyncovden_sum'
    Control_parameter_data(i) % values_character(1) = Covden_sum_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'covden_win_dynamic'
    Covden_win_dynamic = 'dyncovden_win'
    Control_parameter_data(i) % values_character(1) = Covden_win_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'potetcoef_dynamic'
    Potetcoef_dynamic = 'dynpotetcoef'
    Control_parameter_data(i) % values_character(1) = Potetcoef_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'transpbeg_dynamic'
    Transpbeg_dynamic = 'dyntranspbeg'
    Control_parameter_data(i) % values_character(1) = Transpbeg_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'transpend_dynamic'
    Transpend_dynamic = 'dyntranspend'
    Control_parameter_data(i) % values_character(1) = Transpend_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'fallfrost_dynamic'
    Fallfrost_dynamic = 'dynfallfrost'
    Control_parameter_data(i) % values_character(1) = Fallfrost_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'springfrost_dynamic'
    Springfrost_dynamic = 'dynspringfrost'
    Control_parameter_data(i) % values_character(1) = Springfrost_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'soilrechr_dynamic'
    Soilrechr_dynamic = 'dynsoilrechr'
    Control_parameter_data(i) % values_character(1) = Soilrechr_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'soilmoist_dynamic'
    Soilmoist_dynamic = 'dynsoilmoist'
    Control_parameter_data(i) % values_character(1) = Soilmoist_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'radtrncf_dynamic'
    Radtrncf_dynamic = 'dynradtrncf'
    Control_parameter_data(i) % values_character(1) = Radtrncf_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'sro2dprst_perv_dynamic'
    Sro2dprst_perv_dyn = 'dynsro2dprst_perv'
    Control_parameter_data(i) % values_character(1) = Sro2dprst_perv_dyn
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'sro2dprst_imperv_dynamic'
    Sro2dprst_imperv_dyn = 'dynsro2dprst_imperv'
    Control_parameter_data(i) % values_character(1) = Sro2dprst_imperv_dyn
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'transp_on_dynamic'
    Transp_on_dynamic = 'dyntranspon'
    Control_parameter_data(i) % values_character(1) = Transp_on_dynamic
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'dynamic_param_log_file'
    Dynamic_param_log_file = 'dynamic_param.log'
    Control_parameter_data(i) % values_character(1) = Dynamic_param_log_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'selectDatesFileName'
    SelectDatesFileName = 'print_dates.dat'
    Control_parameter_data(i) % values_character(1) = SelectDatesFileName
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'precip_map_file'
    Precip_map_file = 'precip.map'
    Control_parameter_data(i) % values_character(1) = Precip_map_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'tmax_map_file'
    Tmax_map_file = 'tmax.map'
    Control_parameter_data(i) % values_character(1) = Tmax_map_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    i = i + 1
    Control_parameter_data(i) % name = 'tmin_map_file'
    Tmin_map_file = 'tmin.map'
    Control_parameter_data(i) % values_character(1) = Tmin_map_file
    Control_parameter_data(i) % data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'pkwater_equiv_day'
    !      Pkwater_equiv_day = 'pkwater_equiv.day'
    !      Control_parameter_data(i)%values_character(1) = Pkwater_equiv_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'pk_depth_day'
    !      Pk_depth_day = 'pk_depth.day'
    !      Control_parameter_data(i)%values_character(1) = Pk_depth_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'snow_evap_day'
    !      Snow_evap_day = 'snow_evap.day'
    !      Control_parameter_data(i)%values_character(1) = Snow_evap_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'snowcov_area_day'
    !      Snowcov_area_day = 'snowcov_area.day'
    !      Control_parameter_data(i)%values_character(1) = Snowcov_area_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'snowmelt_day'
    !      Snowmelt_day = 'snowmelt.day'
    !      Control_parameter_data(i)%values_character(1) = Snowmelt_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE
    !      i = i + 1
    !      Control_parameter_data(i)%name = 'gwres_flow_day'
    !      Gwres_flow_day = 'gwres_flow.day'
    !      Control_parameter_data(i)%values_character(1) = Gwres_flow_day
    !      Control_parameter_data(i)%data_type = CHAR_TYPE

    ! time arrays
    i = i + 1
    Control_parameter_data(i) % name = 'start_time'
    deallocate (Control_parameter_data(i) % values_int)
    allocate (Control_parameter_data(i) % values_int(6))
    Starttime(1) = 2000
    Starttime(2) = 10
    Starttime(3) = 1
    Starttime(4) = 0
    Starttime(5) = 0
    Starttime(6) = 0
    Control_parameter_data(i) % values_int = Starttime
    Control_parameter_data(i) % numvals = 6
    i = i + 1
    Control_parameter_data(i) % name = 'end_time'
    deallocate (Control_parameter_data(i) % values_int)
    allocate (Control_parameter_data(i) % values_int(6))
    Endtime(1) = 2001
    Endtime(2) = 9
    Endtime(3) = 30
    Endtime(4) = 0
    Endtime(5) = 0
    Endtime(6) = 0
    Control_parameter_data(i) % values_int = Endtime
    Control_parameter_data(i) % numvals = 6

    Num_control_parameters = i

  end subroutine setup_cont

  !***********************************************************************
  ! Get Control File from command line or user interaction.
  !***********************************************************************
  module subroutine get_control_filename()
    use PRMS_CONSTANTS, only: MAXCMDLINE_LENGTH, ERROR_control
    use PRMS_MODULE, only: Print_debug, EQULS, Model_control_file
    use prms_utils, only: error_stop
    implicit none
    ! Functions
    intrinsic :: GET_COMMAND_ARGUMENT, COMMAND_ARGUMENT_COUNT, GET_COMMAND, TRIM
    ! Local Variables
    character(LEN=MAXCMDLINE_LENGTH) command_line_arg, command_line
    logical :: exists
    integer :: status, nchars, numargs
    !***********************************************************************
    ! Subroutine GET_COMMAND_ARGUMENT may not be available with all compilers-it is a Fortran 2003 routine
    ! This routine expects the Control File name to be the first argument, if present
    call GET_COMMAND(command_line)
    print *, 'Command line: ', trim(command_line)
    numargs = COMMAND_ARGUMENT_COUNT()
    if (Print_debug > -1) print '(/,A)', EQULS
    call GET_COMMAND_ARGUMENT(1, command_line_arg, nchars, status)
    if (status /= 0) then
      write (*, '(/,A)') 'Enter the name of the PRMS Control File or quit:'
      read (*, '(A)') Model_control_file
      if (Model_control_file(:4) == 'quit' .or. Model_control_file(:4) == 'QUIT') ERROR stop ERROR_control
    else
      if (trim(command_line_arg(:2)) == '-C') then
        Model_control_file = trim(command_line_arg(3:))
        !          CALL GET_COMMAND_ARGUMENT(2, Model_control_file, nchars, status)
        if (status /= 0) call error_stop('bad argment value after -C argument', ERROR_control)
      else
        Model_control_file = trim(command_line_arg)
      end if
    end if

    inquire (FILE=trim(Model_control_file), EXIST=exists)
    if (.not. exists) then
      write (*, '(/,A)') 'Control File does not exist, file name: '//trim(Model_control_file)
      print *, 'Note: Control File names cannot include spaces'
      ERROR stop - 2
    end if

  end subroutine get_control_filename

  !***********************************************************************
  ! Get Control File set arguments from command line.
  !***********************************************************************
  module subroutine get_control_arguments()
    use PRMS_CONSTANTS, only: DEBUG_less, MAXFILE_LENGTH, ERROR_control
    use PRMS_MODULE, only: Print_debug, EQULS
    use prms_utils, only: error_stop
    implicit none
    ! Functions
    intrinsic :: GET_COMMAND_ARGUMENT, COMMAND_ARGUMENT_COUNT, GET_COMMAND, trim
    ! Local Variables
    character(LEN=MAXFILE_LENGTH) command_line_arg, command_line
    integer :: status, i, j, ii, nchars, numargs, index, param_type, num_param_values
    !***********************************************************************
    ! Subroutine GET_COMMAND_ARGUMENT may not be available with all compilers-it is a Fortran 2003 routine
    ! This routine expects the Control File name to be the first argument, if present
    call GET_COMMAND(command_line)
    numargs = COMMAND_ARGUMENT_COUNT()
    PlotsON_OFF = 0
    i = 0
    do while (i < numargs)
      i = i + 1
      call GET_COMMAND_ARGUMENT(i, command_line_arg, nchars, status)
      if (status /= 0) call error_stop('setting control parameters from command line', ERROR_control)
      if (trim(command_line_arg(:2)) == '-C') then ! rsr, this doesn't work if no space
        i = i + 2
        cycle
      else
        if (Print_debug > -1) print *, 'PRMS command line argument,', i, ': ', trim(command_line_arg)
        if (i == 1) cycle
        if (trim(command_line_arg) == '-rtg') then
          PlotsON_OFF = 1
        elseif (trim(command_line_arg) == '-set') then
          ! find control file parameter and reset it, need type and number of values
          i = i + 1
          call GET_COMMAND_ARGUMENT(i, command_line_arg, nchars, status)
          if (status /= 0) call error_stop('bad argment value after -set argument', ERROR_control)
          if (Print_debug > -1) print *, 'PRMS command line argument,', i, ': ', trim(command_line_arg)
          index = 0
          do j = 1, Num_control_parameters
            if (trim(command_line_arg) == trim(Control_parameter_data(j) % name)) then
              param_type = Control_parameter_data(j) % data_type
              num_param_values = Control_parameter_data(j) % numvals
              index = j
              exit
            end if
          end do
          if (index == 0) call error_stop('control parameter argument not found', ERROR_control)
          do j = 1, num_param_values
            i = i + 1
            call GET_COMMAND_ARGUMENT(i, command_line_arg, nchars, status)
            if (status /= 0) call error_stop('bad value after -set argument', ERROR_control)
            if (Print_debug > -1) print *, 'PRMS command line argument,', i, ': ', trim(command_line_arg)
            if (param_type == 1) then
              if (Control_parameter_data(index) % name(:10) == 'start_time' .or. &
                  Control_parameter_data(index) % name(:8) == 'end_time') then
                read (command_line_arg, *, IOSTAT=status) &
                  (Control_parameter_data(index) % values_int(ii), ii=1, num_param_values)
                !                  print *, 'values ', Control_parameter_data(index)%values_int
                exit
              end if
              read (command_line_arg, *, IOSTAT=status) Control_parameter_data(index) % values_int(j)
              if (status /= 0) call error_stop('reading integer command line argument', ERROR_control)
            elseif (param_type == 4) then
              Control_parameter_data(index) % values_character(j) = command_line_arg
            elseif (param_type == 2) then
              read (command_line_arg, *) Control_parameter_data(index) % values_real(j)
            else
              call error_stop('control parameter type not implemented', ERROR_control)
            end if
          end do
        else
          call error_stop('command line argument invalid', ERROR_control)
        end if
      end if
    end do

    if (Print_debug > DEBUG_less) print '(A)', EQULS
  end subroutine get_control_arguments

  !***********************************************************************
  ! check control parameter if in Control File
  !***********************************************************************
  module subroutine set_control_parameter(Paramname, Numvalues, Paramval_int, Paramval_real, Paramval_char) ! allow arrays
    use PRMS_CONSTANTS, only: ERROR_control
    use prms_utils, only: error_stop
    implicit none
    ! Arguments
    character(LEN=MAXCONTROL_LENGTH), intent(IN) :: Paramname
    integer, intent(IN) :: Numvalues
    integer, intent(IN) :: Paramval_int(Numvalues)
    real, intent(IN) :: Paramval_real(Numvalues)
    character(LEN=MAXFILE_LENGTH), intent(IN) :: Paramval_char(Numvalues)
    ! Functions
    intrinsic :: trim
    ! Local Variables
    integer :: i, j, found, dtype
    !***********************************************************************
    found = 0
    do i = 1, Num_control_parameters
      if (trim(Paramname) == trim(Control_parameter_data(i) % name)) then
        found = i
        dtype = Control_parameter_data(i) % data_type
        Control_parameter_data(i) % numvals = Numvalues
        if (Control_parameter_data(i) % allocate_flag == 1) then ! one of variably sized parameters
          if (dtype == 1) then
            deallocate (Control_parameter_data(i) % values_int)
            allocate (Control_parameter_data(i) % values_int(Numvalues))
          elseif (dtype == 4) then
            deallocate (Control_parameter_data(i) % values_character)
            allocate (Control_parameter_data(i) % values_character(Numvalues))
          else
            call error_stop('allocatable control parameter that is real', ERROR_control)
          end if
        end if
        Control_parameter_data(i) % read_flag = 1
        if (dtype == 1) then
          do j = 1, Numvalues
            Control_parameter_data(i) % values_int(j) = Paramval_int(j)
          end do
        elseif (dtype == 4) then
          do j = 1, Numvalues
            Control_parameter_data(i) % values_character(j) = Paramval_char(j)
          end do
        else !IF ( dtype==2 ) THEN
          do j = 1, Numvalues
            Control_parameter_data(i) % values_real(j) = Paramval_real(j)
          end do
        end if
        exit
      end if
    end do

    if (found == 0) print *, 'WARNING, control parameter not used: ', trim(Paramname), ', ignored'
  end subroutine set_control_parameter

  !***********************************************************************
  ! control_integer
  ! control parameters are read, this sets integer values stored in the
  ! data base and checks to be sure a required parameter has a value (read or default)
  !***********************************************************************
  integer module function control_integer(Parmval, Paramname)
    use PRMS_CONSTANTS, only: ERROR_control
    use PRMS_CONTROL_FILE, only: Num_control_parameters, Control_parameter_data, Max_num_control_parameters
    use prms_utils, only: error_stop
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Paramname
    integer, intent(OUT) :: Parmval
    ! Functions
    intrinsic :: trim
    ! Local Variables
    integer :: i, found
    !***********************************************************************
    found = 0
    do i = 1, Num_control_parameters
      if (trim(Paramname) == trim(Control_parameter_data(i) % name)) then
        Parmval = Control_parameter_data(i) % values_int(1)
        found = i
        exit
      end if
    end do
    if (found == 0) then
      Num_control_parameters = Num_control_parameters + 1
      if (Num_control_parameters > Max_num_control_parameters) &
  &       call error_stop('exceeded maximum number of control parameters', ERROR_control)
      print *, 'WARNING, control parameter not in Control File: ', trim(Paramname), ', set to 0'
      Control_parameter_data(Num_control_parameters) % read_flag = 2 ! set to default
      Control_parameter_data(Num_control_parameters) % data_type = 1
      Control_parameter_data(Num_control_parameters) % numvals = 1
      Control_parameter_data(Num_control_parameters) % name = paramname
      Control_parameter_data(Num_control_parameters) % values_int(1) = 0 !???
    end if

    control_integer = 0
  end function control_integer

  !***********************************************************************
  ! control_integer_array
  ! control parameters are read are read and verified this
  ! function checks to be sure a required parameter has a value (read or default)
  !***********************************************************************
  integer module function control_integer_array(Parmval, Array_index, Paramname)
    use PRMS_CONSTANTS, only: ERROR_control
    use PRMS_CONTROL_FILE, only: Control_parameter_data, Num_control_parameters
    use prms_utils, only: error_stop
    implicit none
    ! Arguments
    ! Array_index and Parmval not used, only used with MMF
    integer, intent(IN) :: Array_index
    character(LEN=*), intent(IN) :: Paramname
    integer, intent(OUT) :: Parmval
    ! Functions
    intrinsic :: trim
    ! Local Variables
    integer :: found, i
    !***********************************************************************
    found = 0
    do i = 1, Num_control_parameters
      if (trim(Paramname) == trim(Control_parameter_data(i) % name)) then
        Parmval = Control_parameter_data(i) % values_int(Array_index)
        found = i
        exit
      end if
    end do
    if (found == 0) then
      print *, 'invalid array control parameter: ', trim(Paramname)
      call error_stop('execution terminated', ERROR_control)
    end if

    control_integer_array = 0
  end function control_integer_array

  !***********************************************************************
  ! control_string
  ! control parameters are read are read and verified this
  ! function checks to be sure a required parameter has a value (read or default)
  !***********************************************************************
  integer module function control_string(Parmval, Paramname)
    use PRMS_CONTROL_FILE, only: Num_control_parameters, Control_parameter_data, Max_num_control_parameters
    implicit none
    ! Functions
    intrinsic :: trim
    ! Arguments
    character(LEN=*), intent(IN) :: Paramname
    character(LEN=*), intent(OUT) :: Parmval
    ! Local Variables
    integer :: found, i
    !***********************************************************************
    found = 0
    do i = 1, Num_control_parameters
      if (trim(Paramname) == trim(Control_parameter_data(i) % name)) then
        Parmval = Control_parameter_data(i) % values_character(1)
        found = i
        exit
      end if
    end do
    if (found == 0) then
      Num_control_parameters = Num_control_parameters + 1
      print *, 'WARNING, control parameter not in Control File: ', trim(Paramname), ', set to blank'
      Control_parameter_data(Num_control_parameters) % read_flag = 2 ! set to default
      Control_parameter_data(Num_control_parameters) % data_type = 4
      Control_parameter_data(Num_control_parameters) % numvals = 1
      Control_parameter_data(Num_control_parameters) % name = paramname
      Control_parameter_data(Num_control_parameters) % values_character(1) = ' '
    end if

    control_string = 0
  end function control_string

  !***********************************************************************
  ! control_string_array
  ! control parameters are read are read and verified this
  ! function checks to be sure a required parameter has a value (read or default)
  !***********************************************************************
  integer module function control_string_array(Parmval, Paramname, Array_index)
    use PRMS_CONSTANTS, only: ERROR_control
    use PRMS_CONTROL_FILE, only: Control_parameter_data, Num_control_parameters
    use prms_utils, only: error_stop
    implicit none
    ! Arguments
    ! Array_index and Parmval not used, only used with MMF
    integer, intent(IN) :: Array_index
    character(LEN=*), intent(IN) :: Paramname
    character(LEN=*), intent(OUT) :: Parmval
    ! Functions
    intrinsic :: trim
    ! Local Variables
    integer :: found, i
    !***********************************************************************
    found = 0
    do i = 1, Num_control_parameters
      if (trim(Paramname) == trim(Control_parameter_data(i) % name)) then
        Parmval = Control_parameter_data(i) % values_character(Array_index)
        found = i
        exit
      end if
    end do
    if (found == 0) then
      print *, 'invalid array control parameter: ', trim(Paramname)
      call error_stop('execution terminated', ERROR_control)
    end if

    control_string_array = 0
  end function control_string_array
end submodule
