!***********************************************************************
! Declares and initializes climate and flow parameters and variables
!***********************************************************************
      MODULE PRMS_CLIMATEVARS
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ACTIVE, OFF, MONTHS_PER_YEAR, ERROR_dim, DEBUG_less, &
     &    potet_pt_module, potet_pm_module, potet_pm_sta_module, climate_hru_module, &
     &    precip_laps_module, xyz_dist_module, ide_dist_module, temp_1sta_module, &
     &    temp_laps_module, temp_sta_module, temp_dist2_module, potet_pan_module, &
     &    FEET, FEET2METERS, METERS2FEET, FAHRENHEIT, INACTIVE, LAKE, ERROR_PARAM, &
     &    ddsolrad_module, ccsolrad_module
      USE PRMS_MODULE, ONLY: Nhru, Nssr, Ngw, Nsegment, Nevap, Nlake, Ntemp, Nrain, Nsol, &
     &    Model, Print_debug, Init_vars_from_file, Temp_flag, Precip_flag, &
     &    Strmflow_module, Temp_module, Stream_order_flag, GSFLOW_flag, &
     &    Precip_module, Solrad_module, Transp_module, Et_module, PRMS4_flag, &
     &    Soilzone_module, Srunoff_module, Call_cascade, Et_flag, Dprst_flag, Solrad_flag, &
     &    Parameter_check_flag, Inputerror_flag, Humidity_cbh_flag, Glacier_flag
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Common States and Fluxes'
      character(len=11), parameter :: MODNAME = 'climateflow'
      character(len=*), parameter :: Version_climateflow = '2021-01-11'
      INTEGER, SAVE :: Use_pandata, Solsta_flag
      ! Tmax_hru and Tmin_hru are in temp_units
      REAL, SAVE, ALLOCATABLE :: Tmax_hru(:), Tmin_hru(:)
      REAL, SAVE, ALLOCATABLE :: Tsta_elev_feet(:), Tsta_elev_meters(:)
      REAL, SAVE, ALLOCATABLE :: Psta_elev_feet(:), Psta_elev_meters(:)
      REAL, SAVE, ALLOCATABLE :: Tmax_allsnow_f(:, :), Tmax_allsnow_c(:, :)
      REAL, SAVE, ALLOCATABLE :: Tmax_allrain_f(:, :)
!   Declared Variables - Precip
      INTEGER, SAVE, ALLOCATABLE :: Newsnow(:), Pptmix(:)
      DOUBLE PRECISION, SAVE :: Basin_ppt, Basin_rain, Basin_snow, Basin_obs_ppt
      REAL, SAVE, ALLOCATABLE :: Hru_ppt(:), Hru_rain(:), Hru_snow(:), Prmx(:)
!   Declared Variables - Temp
      DOUBLE PRECISION, SAVE :: Basin_temp, Basin_tmax, Basin_tmin
      REAL, SAVE :: Solrad_tmax, Solrad_tmin
      REAL, SAVE, ALLOCATABLE :: Tmaxf(:), Tminf(:), Tavgf(:)
      REAL, SAVE, ALLOCATABLE :: Tmaxc(:), Tminc(:), Tavgc(:)
!   Declared Variables - Transp
      INTEGER, SAVE :: Basin_transp_on
      INTEGER, SAVE, ALLOCATABLE :: Transp_on(:)
!   Declared Parameters and Variables - Potential ET
      DOUBLE PRECISION, SAVE :: Basin_potet, Basin_humidity
      REAL, SAVE, ALLOCATABLE :: Potet(:)
      ! for potet_pt, potet_pm and potet_pm_sta
      REAL, SAVE, ALLOCATABLE :: Tempc_dewpt(:), Vp_actual(:), Lwrad_net(:), Vp_slope(:)
      REAL, SAVE, ALLOCATABLE :: Vp_sat(:)
      REAL, SAVE, ALLOCATABLE :: Humidity_percent(:, :)
!   Declared Parameters and Variables - Solar Radiation
      INTEGER, SAVE :: Basin_solsta
      INTEGER, SAVE, ALLOCATABLE :: Hru_solsta(:), Hru_pansta(:)
      DOUBLE PRECISION, SAVE :: Basin_potsw, Basin_swrad, Basin_orad, Basin_horad
      REAL, SAVE :: Rad_conv, Orad
      REAL, SAVE, ALLOCATABLE :: Swrad(:), Orad_hru(:)
      REAL, SAVE, ALLOCATABLE :: Ppt_rad_adj(:, :), Radmax(:, :), Radj_sppt(:), Radj_wppt(:)
!   Declared Parameters - Temp
      INTEGER, SAVE :: Temp_units, Basin_tsta
      INTEGER, SAVE, ALLOCATABLE :: Hru_tsta(:)
      REAL, SAVE, ALLOCATABLE :: Tsta_elev(:), Tmax_aspect_adjust(:, :), Tmin_aspect_adjust(:, :)
!   Declared Parameters - Precip
      INTEGER, SAVE :: Precip_units
      REAL, SAVE :: Ppt_zero_thresh
      REAL, SAVE, ALLOCATABLE :: Tmax_allsnow(:, :), Adjmix_rain(:, :), Tmax_allrain(:, :), Tmax_allrain_offset(:, :)
      REAL, SAVE, ALLOCATABLE :: Psta_elev(:)
!   Declared Parameters - Intcp
      REAL, SAVE, ALLOCATABLE :: Epan_coef(:, :), Potet_sublim(:)
      END MODULE PRMS_CLIMATEVARS

!***********************************************************************
! Declares parameters and variables related to flows from soilzone,
! smbal, ssflow, srunoff_carea, srunoff_smidx
!***********************************************************************
      MODULE PRMS_FLOWVARS
      IMPLICIT NONE
!   Declared Variables
      ! snow
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Pkwater_equiv(:)
      ! soilzone
      DOUBLE PRECISION, SAVE :: Basin_ssflow, Basin_soil_to_gw
      DOUBLE PRECISION, SAVE :: Basin_actet, Basin_lakeevap
      DOUBLE PRECISION, SAVE :: Basin_swale_et, Basin_perv_et, Basin_sroff
      DOUBLE PRECISION, SAVE :: Basin_soil_moist, Basin_ssstor
      REAL, SAVE, ALLOCATABLE :: Hru_actet(:), Soil_moist(:)
      REAL, SAVE, ALLOCATABLE :: Soil_to_gw(:), Slow_flow(:)
      REAL, SAVE, ALLOCATABLE :: Soil_to_ssr(:), Ssres_in(:)
      REAL, SAVE, ALLOCATABLE :: Ssr_to_gw(:), Slow_stor(:)
      REAL, SAVE, ALLOCATABLE :: Ssres_stor(:), Ssres_flow(:), Soil_rechr(:)
      ! srunoff
      REAL, SAVE, ALLOCATABLE :: Sroff(:), Imperv_stor(:), Infil(:)
      ! Surface-Depression Storage
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Dprst_vol_open(:), Dprst_vol_clos(:)
      ! gwflow
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gwres_stor(:)
      ! lakes
      DOUBLE PRECISION, SAVE :: Basin_lake_stor
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_vol(:)
      ! streamflow
      DOUBLE PRECISION, SAVE :: Basin_cfs, Basin_cms, Basin_ssflow_cfs, Basin_sroff_cfs
      DOUBLE PRECISION, SAVE :: Basin_stflow_in, Basin_gwflow_cfs, Basin_stflow_out, Flow_out
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Seg_upstream_inflow(:), Seg_lateral_inflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Seg_outflow(:), Seg_inflow(:)
      ! glacr
      REAL, SAVE, ALLOCATABLE :: Glacier_frac(:), Alt_above_ela(:), Glrette_frac(:)
!   Declared Parameters
      REAL, SAVE, ALLOCATABLE :: Soil_moist_max(:), Soil_rechr_max(:), Sat_threshold(:)
      REAL, SAVE, ALLOCATABLE :: Snowinfil_max(:), Imperv_stor_max(:)
      END MODULE PRMS_FLOWVARS

!***********************************************************************
!     Main climateflow routine
!***********************************************************************
      INTEGER FUNCTION climateflow()
      USE PRMS_CONSTANTS, ONLY: DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: climateflow_decl, climateflow_init
      EXTERNAL :: climateflow_restart
!***********************************************************************
      climateflow = 0

      IF ( Process_flag==DECL ) THEN
        climateflow = climateflow_decl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL climateflow_restart(READ_INIT)
        climateflow = climateflow_init()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL climateflow_restart(SAVE_INIT)
      ENDIF

      END FUNCTION climateflow

!***********************************************************************
!     climateflow_decl - declare climate and flow variables and parameters
!***********************************************************************
      INTEGER FUNCTION climateflow_decl()
      USE PRMS_CLIMATEVARS
      USE PRMS_FLOWVARS
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declvar, declparam
      EXTERNAL :: read_error, print_module
!***********************************************************************
      climateflow_decl = 0

      CALL print_module(MODDESC, MODNAME, Version_climateflow)

      ALLOCATE ( Tmaxf(Nhru) )
      IF ( declvar(Temp_module, 'tmaxf', 'nhru', Nhru, 'real', &
     &     'Maximum air temperature distributed to each HRU', &
     &     'degrees Fahrenheit', Tmaxf)/=0 ) CALL read_error(3, 'tmaxf')

      ALLOCATE ( Tminf(Nhru) )
      IF ( declvar(Temp_module, 'tminf', 'nhru', Nhru, 'real', &
     &     'Minimum air temperature distributed to each HRU', &
     &     'degrees Fahrenheit', Tminf)/=0 ) CALL read_error(3, 'tminf')

      ALLOCATE ( Tavgf(Nhru) )
      IF ( declvar(Temp_module, 'tavgf', 'nhru', Nhru, 'real', &
     &     'Average air temperature distributed to each HRU', &
     &     'degrees Fahrenheit', Tavgf)/=0 ) CALL read_error(3, 'tavgf')

      ALLOCATE ( Tmaxc(Nhru) )
      IF ( declvar(Temp_module, 'tmaxc', 'nhru', Nhru, 'real', &
     &     'Maximum air temperature distributed to each HRU', &
     &     'degrees Celsius', Tmaxc)/=0 ) CALL read_error(3, 'tmaxc')

      ALLOCATE ( Tminc(Nhru) )
      IF ( declvar(Temp_module, 'tminc', 'nhru', Nhru, 'real', &
     &     'Minimum air temperature distributed to each HRU', &
     &     'degrees Celsius', Tminc)/=0 ) CALL read_error(3, 'tminc')

      ALLOCATE ( Tavgc(Nhru) )
      IF ( declvar(Temp_module, 'tavgc', 'nhru', Nhru, 'real', &
     &     'Average air temperature distributed to each HRU', &
     &     'degrees Celsius', Tavgc)/=0 ) CALL read_error(3, 'tavgc')

      IF ( declvar(Temp_module, 'basin_tmax', 'one', 1, 'double', &
     &     'Basin area-weighted average maximum air temperature', &
     &     'temp_units', Basin_tmax)/=0 ) CALL read_error(3, 'basin_tmax')

      IF ( declvar(Temp_module, 'basin_tmin', 'one', 1, 'double', &
     &     'Basin area-weighted average minimum air temperature', &
     &     'temp_units', Basin_tmin)/=0 ) CALL read_error(3, 'basin_tmin')

      IF ( declvar(Temp_module, 'basin_temp', 'one', 1, 'double', &
     &     'Basin area-weighted average air temperature', &
     &     'temp_units', Basin_temp)/=0 ) CALL read_error(3, 'basin_temp')

      IF ( declvar(Temp_module, 'solrad_tmax', 'one', 1, 'real', &
     &     'Basin daily maximum temperature for use with solar radiation calculations', &
     &     'temp_units', Solrad_tmax)/=0 ) CALL read_error(3, 'solrad_tmax')

      IF ( declvar(Temp_module, 'solrad_tmin', 'one', 1, 'real', &
     &     'Basin daily minimum temperature for use with solar radiation calculations', &
     &     'temp_units', Solrad_tmin)/=0 ) CALL read_error(3, 'solrad_tmin')

! PRECIPITATION VARIABLES AND PARAMETERS
      ALLOCATE ( Pptmix(Nhru) )
      IF ( declvar(Precip_module, 'pptmix', 'nhru', Nhru, 'integer', &
     &     'Flag to indicate if precipitation is a mixture of rain'// &
     &     ' and snow for each HRU (0=no; 1=yes)', &
     &     'none', Pptmix)/=0 ) CALL read_error(3, 'pptmix')

      ALLOCATE ( Newsnow(Nhru) )
      IF ( declvar(Precip_module, 'newsnow', 'nhru', Nhru, 'integer', &
     &    'Flag to indicate if new snow fell on each HRU (0=no; 1=yes)', &
     &    'none', Newsnow)/=0 ) CALL read_error(3, 'newsnow')

      ALLOCATE ( Prmx(Nhru) )
      IF ( declvar(Precip_module, 'prmx', 'nhru', Nhru, 'real', &
     &     'Fraction of rain in a mixed precipitation event for each HRU', &
     &     'decimal fraction', Prmx)/=0 ) CALL read_error(3, 'prmx')

      IF ( declvar(Precip_module, 'basin_rain', 'one', 1, 'double', &
     &     'Basin area-weighted average rainfall', &
     &     'inches', Basin_rain)/=0 ) CALL read_error(3, 'basin_rain')

      IF ( declvar(Precip_module, 'basin_snow', 'one', 1, 'double', &
     &     'Basin area-weighted average snowfall for basin', &
     &     'inches', Basin_snow)/=0 ) CALL read_error(3, 'basin_snow')

      IF ( declvar(Precip_module, 'basin_ppt', 'one', 1, 'double', &
     &     'Basin area-weighted average precipitation', &
     &     'inches', Basin_ppt)/=0 ) CALL read_error(3, 'basin_ppt')

! DANGER - Not sure what to do about this one.  For right now
!          I'm setting basin_ppt and basin_obs_ppt to the same
!          variable.  In the precip_1sta module, basin_obs_ppt
!          seems to be the area weighted precipitation average before
!          the correction factor is applied.  In other modules,
!          the correction "error" is applied to the station
!          precipitation rather than the hru precipitation.
      IF ( declvar(Precip_module, 'basin_obs_ppt', 'one', 1, 'double', &
     &     'Basin area-weighted average measured precipitation', &
     &     'inches', Basin_obs_ppt)/=0 ) CALL read_error(3, 'basin_obs_ppt')

      ALLOCATE ( Hru_ppt(Nhru) )
      IF ( declvar(Precip_module, 'hru_ppt', 'nhru', Nhru, 'real', &
     &     'Precipitation distributed to each HRU', &
     &     'inches', Hru_ppt)/=0 ) CALL read_error(3, 'hru_ppt')

      ALLOCATE ( Hru_rain(Nhru) )
      IF ( declvar(Precip_module, 'hru_rain', 'nhru', Nhru, 'real', &
     &     'Rain distributed to each HRU', &
     &     'inches', Hru_rain)/=0 ) CALL read_error(3, 'hru_rain')

      ALLOCATE ( Hru_snow(Nhru) )
      IF ( declvar(Precip_module, 'hru_snow', 'nhru', Nhru, 'real', &
     &     'Snow distributed to each HRU', &
     &     'inches', Hru_snow)/=0 ) CALL read_error(3, 'hru_snow')

! Solar Radiation variables
      ALLOCATE ( Swrad(Nhru) )
      IF ( declvar(Solrad_module, 'swrad', 'nhru', Nhru, 'real', &
     &     'Shortwave radiation distributed to each HRU', &
     &     'Langleys', Swrad)/=0 ) CALL read_error(3, 'swrad')

      IF ( declvar(Solrad_module, 'orad', 'one', 1, 'real', &
     &     'Measured or computed solar radiation on a horizontal surface', &
     &     'Langleys', Orad)/=0 ) CALL read_error(3, 'orad')

      IF ( declvar(Solrad_module, 'basin_horad', 'one', 1, 'double', &
     &     'Potential shortwave radiation for the basin centroid', &
     &     'Langleys', Basin_horad)/=0 ) CALL read_error(3, 'basin_horad')

      IF ( declvar(Solrad_module, 'basin_swrad', 'one', 1, 'double', &
     &     'Basin area-weighted average shortwave radiation', &
     &     'Langleys', Basin_swrad)/=0 ) CALL read_error(3, 'basin_swrad')

      IF ( declvar(Solrad_module, 'basin_potsw', 'one', 1, 'double', &
     &     'Basin area-weighted average shortwave radiation', &
     &     'Langleys', Basin_potsw)/=0 ) CALL read_error(3, 'basin_potsw')

      IF ( Solrad_flag==ddsolrad_module .OR. Solrad_flag==ccsolrad_module .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(Solrad_module, 'basin_orad', 'one', 1, 'double', &
     &       'Basin area-weighted average solar radiation on a horizontal surface', &
     &       'Langleys', Basin_orad)/=0 ) CALL read_error(3, 'basin_orad')

        ALLOCATE ( Orad_hru(Nhru) )
        IF ( declvar(Solrad_module, 'orad_hru', 'nhru', Nhru, 'real', &
     &       'Solar radiation on a horizontal surface for each HRU', &
     &       'Langleys', Orad_hru)/=0 ) CALL read_error(3, 'orad_hru')
      ENDIF

! Transpiration Variables
      ALLOCATE ( Transp_on(Nhru) )
      IF ( declvar(Transp_module, 'transp_on', 'nhru', Nhru, 'integer', &
     &     'Flag indicating whether transpiration is occurring (0=no; 1=yes)', &
     &     'none', Transp_on)/=0 ) CALL read_error(3, 'transp_on')

      IF ( declvar(Transp_module, 'basin_transp_on', 'one', 1,'integer', &
     &     'Flag indicating whether transpiration is occurring anywhere in the basin (0=no; 1=yes)', &
     &     'none', Basin_transp_on)/=0 ) CALL read_error(3, 'basin_transp_on')

! Potential ET Variables
      ALLOCATE ( Potet(Nhru) )
      IF ( declvar(Et_module, 'potet', 'nhru', Nhru, 'real', &
     &     'Potential ET for each HRU', &
     &     'inches', Potet)/=0 ) CALL read_error(3, 'potet')

      IF ( declvar(Et_module, 'basin_potet', 'one', 1, 'double', &
     &     'Basin area-weighted average potential ET', &
     &     'inches', Basin_potet)/=0 ) CALL read_error(3, 'basin_potet')

      IF ( Et_flag==potet_pt_module .OR. Et_flag==potet_pm_module .OR. &
     &     Et_flag==potet_pm_sta_module .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(Et_module, 'basin_humidity', 'one', 1, 'double', &
     &       'Basin area-weighted average humidity', &
     &       'percentage', Basin_humidity)/=0 ) CALL read_error(3, 'basin_humidity')
        ALLOCATE ( Tempc_dewpt(Nhru) )
        IF ( declvar(Et_module, 'tempc_dewpt', 'nhru', Nhru, 'real', &
     &       'Air temperature at dew point for each HRU', &
     &       'degrees Celsius', Tempc_dewpt)/=0 ) CALL read_error(3, 'tempc_dewpt')
        ALLOCATE ( Vp_actual(Nhru) )
        IF ( declvar(Et_module, 'vp_actual', 'nhru', Nhru, 'real', &
     &       'Actual vapor pressure for each HRU', &
     &       'kilopascals', Vp_actual)/=0 ) CALL read_error(3, 'vp_actual')
        ALLOCATE ( Lwrad_net(Nhru) )
        IF ( declvar(Et_module, 'lwrad_net', 'nhru', Nhru, 'real', &
     &       'Net long-wave radiation for each HRU', &
     &       'megajoules/m**2/day', Lwrad_net)/=0 ) CALL read_error(3, 'lwrad_net')
        ALLOCATE ( Vp_slope(Nhru) )
        IF ( declvar(Et_module, 'vp_slope', 'nhru', Nhru, 'real', &
     &       'Slope of saturation vapor pressure versus air temperature curve for each HRU', &
     &       'kilopascals/degrees Celsius', Vp_slope)/=0 ) CALL read_error(3, 'vp_slope')
        IF ( Et_flag==potet_pm_module .OR. Et_flag==potet_pm_sta_module .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Vp_sat(Nhru) )
          IF ( declvar(Et_module, 'vp_sat', 'nhru', Nhru, 'real', &
     &         'Saturation vapor pressure for each HRU', &
     &         'kilopascals', Vp_sat)/=0 ) CALL read_error(3, 'vp_sat')
        ENDIF
        ALLOCATE ( Humidity_percent(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(Et_module, 'humidity_percent', 'nhru,nmonths', 'real', &
     &       '0.0', '0.0', '100.0', &
     &       'Monthy humidity for each HRU', &
     &       'Monthy humidity for each HRU', &
     &       'percentage')/=0 ) CALL read_error(1, 'humidity_percent')
      ENDIF

! Soilzone variables
      ALLOCATE ( Soil_rechr(Nhru) )
      IF ( declvar(Soilzone_module, 'soil_rechr', 'nhru', Nhru, 'real', &
     &     'Storage for recharge zone (upper portion) of the'// &
     &     ' capillary reservoir that is available for both'// &
     &     ' evaporation and transpiration', &
     &     'inches', Soil_rechr)/=0 ) CALL read_error(3, 'soil_rechr')

      ALLOCATE ( Ssr_to_gw(Nssr) )
      IF ( declvar(Soilzone_module, 'ssr_to_gw', 'nssr', Nssr, 'real', &
     &     'Drainage from the gravity-reservoir to the associated GWR for each HRU', &
     &     'inches', Ssr_to_gw)/=0 ) CALL read_error(3, 'ssr_to_gw')

      ALLOCATE ( Ssres_stor(Nssr) )
      IF ( declvar(Soilzone_module, 'ssres_stor', 'nssr', Nssr, 'real', &
     &     'Storage in the gravity and preferential-flow reservoirs for each HRU', &
     &     'inches', Ssres_stor)/=0 ) CALL read_error(3, 'ssres_stor')

      ALLOCATE ( Slow_flow(Nhru) )
      IF ( declvar(Soilzone_module, 'slow_flow', 'nhru', Nhru, 'real', &
     &     'Interflow from gravity reservoir storage that flows to'// &
     &     ' the stream network for each HRU', &
     &     'inches', Slow_flow)/=0 ) CALL read_error(3, 'slow_flow')

      ALLOCATE ( Ssres_flow(Nssr) )
      IF ( declvar(Soilzone_module, 'ssres_flow', 'nssr', Nssr, 'real', &
     &     'Interflow from gravity and preferential-flow reservoirs'// &
     &     ' to the stream network for each HRU', &
     &     'inches', Ssres_flow)/=0 ) CALL read_error(3, 'ssres_flow')

      IF ( declvar(Soilzone_module, 'basin_ssflow', 'one', 1, 'double', &
     &     'Basin area-weighted average interflow from gravity and'// &
     &     ' preferential-flow reservoirs to the stream network', &
     &     'inches', Basin_ssflow)/=0 ) CALL read_error(3, 'basin_ssflow')

      IF ( declvar(Soilzone_module, 'basin_swale_et', 'one', 1, 'double', &
     &     'Basin area-weighted average ET from swale HRUs', &
     &     'inches', Basin_swale_et)/=0 ) CALL read_error(3, 'basin_swale_et')

      IF ( declvar(Soilzone_module, 'basin_soil_moist', 'one', 1, 'double', &
     &     'Basin area-weighted average capillary reservoir storage', &
     &     'inches', Basin_soil_moist)/=0 ) CALL read_error(3, 'basin_soil_moist')

      IF ( declvar(Soilzone_module, 'basin_ssstor', 'one', 1, 'double', &
     &     'Basin weighted average gravity and preferential-flow reservoir storage', &
     &     'inches', Basin_ssstor)/=0 ) CALL read_error(3, 'basin_ssstor')

      ALLOCATE ( Slow_stor(Nhru) )
      IF ( declvar(Soilzone_module, 'slow_stor', 'nhru', Nhru, 'real', &
     &     'Storage of gravity reservoir for each HRU', &
     &     'inches', Slow_stor)/=0 ) CALL read_error(3, 'slow_stor')

      ALLOCATE ( Soil_moist(Nhru) )
      IF ( declvar(Soilzone_module, 'soil_moist', 'nhru', Nhru, 'real', &
     &     'Storage of capillary reservoir for each HRU', &
     &     'inches', Soil_moist)/=0 ) CALL read_error(3, 'soil_moist')

      ALLOCATE ( Hru_actet(Nhru) )
      IF ( declvar(Soilzone_module, 'hru_actet', 'nhru', Nhru, 'real', &
     &     'Actual ET for each HRU', &
     &     'inches', Hru_actet)/=0 ) CALL read_error(3, 'hru_actet')

      IF ( declvar(Soilzone_module, 'basin_actet', 'one', 1, 'double', &
     &     'Basin area-weighted average actual ET', &
     &     'inches', Basin_actet)/=0 ) CALL read_error(3, 'basin_actet')

      IF ( declvar(Soilzone_module, 'basin_perv_et', 'one', 1, 'double', &
     &     'Basin area-weighted average ET from capillary reservoirs', &
     &     'inches', Basin_perv_et)/=0 ) CALL read_error(3, 'basin_perv_et')

      IF ( declvar(Soilzone_module, 'basin_lakeevap', 'one', 1, 'double', &
     &     'Basin area-weighted average lake evaporation', &
     &     'inches', Basin_lakeevap)/=0 ) CALL read_error(3, 'basin_lakeevap')

      ALLOCATE ( Ssres_in(Nssr) )
      IF ( declvar(Soilzone_module, 'ssres_in', 'nssr', Nssr, 'real', &
     &     'Inflow to the gravity and preferential-flow reservoirs for each HRU', &
     &     'inches', Ssres_in)/=0 ) CALL read_error(3, 'ssres_in')

      ALLOCATE ( Soil_to_gw(Nhru) )
      IF ( declvar(Soilzone_module, 'soil_to_gw', 'nhru', Nhru, 'real', &
     &     'Portion of excess flow to the capillary reservoir that'// &
     &     ' drains to the associated GWR for each HRU', &
     &     'inches', Soil_to_gw)/=0 ) CALL read_error(3, 'soil_to_gw')

      ALLOCATE ( Soil_to_ssr(Nhru) )
      IF ( declvar(Soilzone_module, 'soil_to_ssr', 'nhru', Nhru, 'real', &
     &     'Portion of excess flow to the capillary reservoir that'// &
     &     ' flows to the gravity reservoir for each HRU', &
     &     'inches', Soil_to_ssr)/=0 ) CALL read_error(3, 'soil_to_ssr')

      IF ( declvar(Soilzone_module, 'basin_soil_to_gw', 'one', 1, 'double', &
     &     'Basin average excess flow to capillary reservoirs that drains to GWRs', &
     &     'inches', Basin_soil_to_gw)/=0 ) CALL read_error(3, 'basin_soil_to_gw')

! gwflow
      IF ( GSFLOW_flag==OFF .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Gwres_stor(Nhru) )
        IF ( declvar('gwflow', 'gwres_stor', 'ngw', Nhru, 'double', &
     &       'Storage in each GWR', &
     &       'inches', Gwres_stor)/=0 ) CALL read_error(3, 'gwres_stor')
      ENDIF

! srunoff
      ALLOCATE ( Imperv_stor(Nhru) )
      IF ( declvar(Srunoff_module, 'imperv_stor', 'nhru', Nhru, 'real', &
     &     'Storage on impervious area for each HRU', &
     &     'inches', Imperv_stor)/=0 ) CALL read_error(3, 'imperv_stor')

      ALLOCATE ( Infil(Nhru) )
      IF ( declvar(Srunoff_module, 'infil', 'nhru', Nhru, 'real', &
     &     'Infiltration to the capillary and preferential-flow reservoirs for each HRU', &
     &     'inches', Infil)/=0 ) CALL read_error(3, 'infil')

      ALLOCATE ( Sroff(Nhru) )
      IF ( declvar(Srunoff_module, 'sroff', 'nhru', Nhru, 'real', &
     &     'Surface runoff to the stream network for each HRU', &
     &     'inches', Sroff)/=0 ) CALL read_error(3, 'sroff')

      IF ( declvar(Srunoff_module, 'basin_sroff', 'one', 1, 'double', &
     &     'Basin area-weighted average surface runoff to the stream network', &
     &     'inches', Basin_sroff)/=0 ) CALL read_error(3, 'basin_sroff')

! stream flow
      IF ( declvar(Strmflow_module, 'basin_cfs', 'one', 1, 'double', &
     &     'Streamflow leaving the basin through the stream network', &
     &     'cfs', Basin_cfs)/=0 ) CALL read_error(3, 'basin_cfs')

      IF ( declvar(Strmflow_module, 'basin_cms', 'one', 1, 'double', &
     &     'Streamflow leaving the basin through the stream network', &
     &     'cms', Basin_cms)/=0 ) CALL read_error(3, 'basin_cms')

      IF ( declvar(Strmflow_module, 'basin_stflow_in', 'one', 1, 'double', &
     &     'Basin area-weighted average lateral flow entering the stream network', &
     &     'inches', Basin_stflow_in)/=0 ) CALL read_error(3, 'basin_stflow_in')

      IF ( declvar(Strmflow_module, 'basin_stflow_out', 'one', 1, 'double', &
     &     'Basin area-weighted average streamflow leaving through the stream network', &
     &     'inches', Basin_stflow_out)/=0 ) CALL read_error(3, 'basin_stflow_out')

      IF ( declvar(Strmflow_module, 'basin_sroff_cfs', 'one', 1, 'double', &
     &     'Basin area-weighted average surface runoff to the stream network', &
     &     'cfs', Basin_sroff_cfs)/=0 ) CALL read_error(3, 'basin_sroff_cfs')

      IF ( declvar(Strmflow_module, 'basin_ssflow_cfs', 'one', 1, 'double', &
     &     'Interflow leaving the basin through the stream network', &
     &     'cfs', Basin_ssflow_cfs)/=0 ) CALL read_error(3, 'basin_ssflow')

      IF ( declvar(Strmflow_module, 'basin_gwflow_cfs', 'one', 1, 'double', &
     &     'Basin area-weighted average of groundwater flow to the stream network', &
     &     'cfs', Basin_gwflow_cfs)/=0 ) CALL read_error(3, 'basin_gwflow_cfs')

      IF ( Call_cascade==1 .OR. Stream_order_flag==ACTIVE ) THEN
        IF ( Nsegment==0 .AND. Model/=DOCUMENTATION ) THEN
          PRINT *, 'ERROR, nsegment=0, must be > 0 for selected module options'
          ERROR STOP ERROR_dim
        ENDIF
      ENDIF

      IF ( Stream_order_flag==ACTIVE ) THEN
        ALLOCATE ( Seg_outflow(Nsegment) )
        IF ( declvar(Strmflow_module, 'seg_outflow', 'nsegment', Nsegment, 'double', &
     &       'Streamflow leaving a segment', &
     &       'cfs', Seg_outflow)/=0 ) CALL read_error(3, 'seg_outflow')

        ALLOCATE ( Seg_inflow(Nsegment) )
        IF ( declvar(Strmflow_module, 'seg_inflow', 'nsegment', Nsegment, 'double', &
     &       'Total flow entering a segment', &
     &       'cfs', Seg_inflow)/=0 ) CALL read_error(3, 'seg_inflow')

        IF ( declvar(Strmflow_module, 'flow_out', 'one', 1, 'double', &
     &       'Total flow out of model domain', &
     &       'cfs', Flow_out)/=0 ) CALL read_error(3, 'flow_out')

        ALLOCATE ( Seg_lateral_inflow(Nsegment) )
        IF ( declvar(Strmflow_module, 'seg_lateral_inflow', 'nsegment', Nsegment, 'double', &
     &       'Lateral inflow entering a segment', &
     &       'cfs', Seg_lateral_inflow)/=0 ) CALL read_error(3, 'seg_lateral_inflow')

        ALLOCATE ( Seg_upstream_inflow(Nsegment) )
        IF ( declvar(Strmflow_module, 'seg_upstream_inflow', 'nsegment', Nsegment, 'double', &
     &       'Sum of inflow from upstream segments', &
     &       'cfs', Seg_upstream_inflow)/=0 ) CALL read_error(3, 'seg_upstream_inflow')
      ENDIF

      IF ( declvar(Strmflow_module, 'basin_lake_stor', 'one', 1, 'double', &
     &     'Basin volume-weighted average storage for all lakes using broad-crested weir or gate opening routing', &
     &     'inches', Basin_lake_stor)/=0 ) CALL read_error(3, 'basin_lake_stor')

      IF ( Nlake>0 ) THEN
        ALLOCATE ( Lake_vol(Nlake) )
        IF ( declvar(Strmflow_module, 'lake_vol', 'nlake', Nlake, 'double', &
     &       'Storage in each lake using broad-crested weir or gate opening routing', &
     &       'acre-feet', Lake_vol)/=0 ) CALL read_error(3, 'lake_vol')
      ENDIF

      IF ( Dprst_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Dprst_vol_open(Nhru) )
        IF ( declvar(Srunoff_module, 'dprst_vol_open', 'nhru', Nhru, 'double', &
     &       'Storage volume in open surface depressions for each HRU', &
     &       'acre-inches', Dprst_vol_open)/=0 ) CALL read_error(3, 'dprst_vol_open')
        ALLOCATE ( Dprst_vol_clos(Nhru) )
        IF ( declvar(Srunoff_module, 'dprst_vol_clos', 'nhru', Nhru, 'double', &
     &       'Storage volume in closed surface depressions for each HRU', &
     &       'acre-inches', Dprst_vol_clos)/=0 ) CALL read_error(3, 'dprst_vol_clos')
      ENDIF

      ALLOCATE ( Pkwater_equiv(Nhru) )
      IF ( declvar('snowcomp', 'pkwater_equiv', 'nhru', Nhru, 'double', &
     &     'Snowpack water equivalent on each HRU', &
     &     'inches', Pkwater_equiv)/=0 ) CALL read_error(3, 'pkwater_equiv')

! glacier variables
      IF ( Glacier_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Glacier_frac(Nhru) )
        IF ( declvar(MODNAME, 'glacier_frac', 'nhru', Nhru, 'real',       &
             'Fraction of glaciation (0=none; 1=100%)',                   &
             'decimal fraction', Glacier_frac)/=0 ) CALL read_error(3, 'glacier_frac')

        ALLOCATE ( Glrette_frac(Nhru) )
          IF ( declvar(MODNAME, 'glrette_frac', 'nhru', Nhru, 'real',     &
             'Fraction of snow field (too small for glacier dynamics)',   &
             'decimal fraction', Glrette_frac)/=0 ) CALL read_error(3, 'glrette_frac')

        ALLOCATE ( Alt_above_ela(Nhru) )
        IF ( declvar(MODNAME, 'alt_above_ela', 'nhru', Nhru, 'real',      &
             'Altitude above equilibrium line altitude (ELA)',            &
             'elev_units', Alt_above_ela)/=0 ) CALL read_error(3, 'alt_above_ela')
      ENDIF

      ! Allocate local variables
      IF ( Temp_flag<climate_hru_module .OR. Model==DOCUMENTATION ) ALLOCATE ( Tsta_elev_meters(Ntemp), Tsta_elev_feet(Ntemp) )
      IF ( Precip_flag==precip_laps_module .OR. Precip_flag==xyz_dist_module .OR. &
     &     Precip_flag==ide_dist_module .OR. Model==DOCUMENTATION ) &
     &     ALLOCATE ( Psta_elev_meters(Nrain), Psta_elev_feet(Nrain) )
      ALLOCATE ( Tmax_hru(Nhru), Tmin_hru(Nhru) )
      ALLOCATE ( Tmax_allsnow_f(Nhru,MONTHS_PER_YEAR), Tmax_allsnow_c(Nhru,MONTHS_PER_YEAR), Tmax_allrain_f(Nhru,MONTHS_PER_YEAR) )

! Declare Parameters
      IF ( Temp_flag<climate_hru_module .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Tsta_elev(Ntemp) )
        IF ( declparam(Temp_module, 'tsta_elev', 'ntemp', 'real', &
     &       '0.0', '-300.0', '30000.0', &
     &       'Temperature station elevation', &
     &       'Elevation of each air-temperature-measurement station', &
     &       'elev_units')/=0 ) CALL read_error(1, 'tsta_elev')
      ENDIF

      IF ( Temp_flag==temp_1sta_module .OR. Temp_flag==temp_laps_module .OR. Temp_flag==temp_sta_module &
     &     .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Hru_tsta(Nhru) )
        IF ( declparam(Temp_module, 'hru_tsta', 'nhru', 'integer', &
     &       '0', 'bounded', 'ntemp', &
     &       'Index of base temperature station for each HRU', &
     &       'Index of the base temperature station used for lapse rate calculations', &
     &       'none')/=0 ) CALL read_error(1, 'hru_tsta')
      ENDIF

      IF ( Temp_flag==temp_1sta_module .OR. Temp_flag==temp_laps_module .OR. Temp_flag==temp_dist2_module .OR. &
     &     Temp_flag==ide_dist_module .OR. Temp_flag==xyz_dist_module .OR. Temp_flag==temp_sta_module &
     &     .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Tmax_aspect_adjust(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(Temp_module, 'tmax_adj', 'nhru,nmonths', 'real', &
     &       '0.0', '-10.0', '10.0', &
     &       'HRU maximum temperature adjustment', &
     &       'Adjustment to maximum temperature for each HRU, estimated on the basis of slope and aspect', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmax_adj')

        ALLOCATE ( Tmin_aspect_adjust(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(Temp_module, 'tmin_adj', 'nhru,nmonths', 'real', &
     &       '0.0', '-10.0', '10.0', &
     &       'HRU minimum temperature adjustment', &
     &       'Adjustment to minimum temperature for each HRU, estimated on the basis of slope and aspect', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmin_adj')
      ENDIF

      ALLOCATE ( Potet_sublim(Nhru) )
      IF ( declparam('snowcomp', 'potet_sublim', 'nhru', 'real', &
     &     '0.5', '0.1', '0.75', &
     &     'Fraction of potential ET that is sublimated from snow for each HRU', &
     &     'Fraction of potential ET that is sublimated from snow in the canopy and snowpack for each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'potet_sublim')

      ALLOCATE ( Tmax_allrain_offset(Nhru,MONTHS_PER_YEAR), Tmax_allrain(Nhru,MONTHS_PER_YEAR) )
      IF ( PRMS4_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(Precip_module, 'tmax_allrain', 'nhru,nmonths', 'real', &
     &       '38.0', '-8.0', '75.0', &
     &       'Precipitation is rain if HRU max temperature >= this value', &
     &       'Monthly (January to December) maximum air temperature'// &
     &       ' when precipitation is assumed to be rain; if HRU air'// &
     &       ' temperature is greater than or equal to this value, precipitation is rain', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmax_allrain')
      ENDIF
      IF ( PRMS4_flag==OFF .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(Precip_module, 'tmax_allrain_offset', 'nhru,nmonths', 'real', &
     &       '1.0', '0.0', '50.0', &
     &       'Precipitation is rain if HRU max temperature >= tmax_allsnow + this value', &
     &       'Monthly (January to December) maximum air temperature'// &
     &       ' when precipitation is assumed to be rain; if HRU air'// &
     &       ' temperature is greater than or equal to tmax_allsnow plus this value, precipitation is rain', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmax_allrain_offset')
      ENDIF

      ALLOCATE ( Tmax_allsnow(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(Precip_module, 'tmax_allsnow', 'nhru,nmonths', 'real', &
     &     '32.0', '-10.0', '40.0', &
     &     'Maximum temperature when precipitation is all snow', &
     &     'Maximum air temperature when precipitation is assumed'// &
     &     ' to be snow; if HRU air temperature is less than or equal to this value, precipitation is snow', &
     &     'temp_units')/=0 ) CALL read_error(1, 'tmax_allsnow')

      ALLOCATE ( Adjmix_rain(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(Precip_module, 'adjmix_rain', 'nhru,nmonths', 'real', &
     &     '1.0', '0.0', '3.0', &
     &     'Adjustment factor for rain in a rain/snow mix', &
     &     'Monthly (January to December) factor to adjust rain proportion in a mixed rain/snow event', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'adjmix_rain')

      IF ( declparam(Precip_module, 'ppt_zero_thresh', 'one', 'real', &
     &     '0.0', '0.0', '0.1', &
     &     'Precipitation below this amount is set to 0.0', &
     &     'Precipitation below this amount is set to 0.0', &
     &     'precip_units')/=0 ) CALL read_error(1, 'ppt_zero_thresh')

      IF ( Precip_flag==precip_laps_module .OR. Precip_flag==xyz_dist_module .OR. Precip_flag==ide_dist_module &
     &     .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Psta_elev(Nrain) )
        IF ( declparam(Precip_module, 'psta_elev', 'nrain', 'real', &
     &       '0.0', '-300.0', '30000.0', &
     &       'Precipitation station elevation', &
     &       'Elevation of each precipitation measurement station', &
     &       'elev_units')/=0 ) CALL read_error(1, 'psta_elev')
      ENDIF

      IF ( declparam(Temp_module, 'temp_units', 'one', 'integer', &
     &     '0', '0', '1', &
     &     'Units flag for measured temperature', &
     &     'Flag to indicate the units of measured air-temperature values (0=Fahrenheit; 1=Celsius)', &
     &     'none')/=0 ) CALL read_error(1, 'temp_units')

      IF ( Temp_flag<ide_dist_module .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(Temp_module, 'basin_tsta', 'one', 'integer', &
     &       '0', 'bounded', 'ntemp', &
     &       'Index of main temperature station', &
     &       'Index of temperature station used to compute basin temperature values', &
     &       'none')/=0 ) CALL read_error(1, 'basin_tsta')
      ENDIF

      IF ( declparam(Precip_module, 'precip_units', 'one', 'integer', &
     &     '0', '0', '1', &
     &     'Units for measured precipitation', &
     &     'Units for measured precipitation (0=inches; 1=mm)', &
     &     'none')/=0 ) CALL read_error(1, 'precip_units')

      IF ( Solrad_flag==ddsolrad_module .OR. Solrad_flag==ccsolrad_module .OR. Model==DOCUMENTATION ) THEN
        IF ( Nsol>0 ) THEN
          IF ( declparam(Solrad_module, 'rad_conv', 'one', 'real', &
     &         '1.0', '0.1', '100.0', &
     &         'Conversion factor to Langleys for measured radiation', &
     &         'Conversion factor to Langleys for measured solar radiation', &
     &         'Langleys/radiation units')/=0 ) CALL read_error(1, 'rad_conv')

          IF ( declparam(Solrad_module, 'basin_solsta', 'one', 'integer', &
     &         '0', 'bounded', 'nsol', &
     &         'Index of main solar radiation station', &
     &         'Index of solar radiation station used to compute basin radiation values, used when dimension nsol>0', &
     &         'none')/=0 ) CALL read_error(1, 'basin_solsta')

          ALLOCATE ( Hru_solsta(Nhru) )
          IF ( declparam(Solrad_module, 'hru_solsta', 'nhru', 'integer', &
     &         '0', 'bounded', 'nsol', &
     &         'Index of solar radiation station associated with each HRU', &
     &         'Index of solar radiation station associated with each HRU', &
     &         'none')/=0 ) CALL read_error(1, 'hru_solsta')
        ENDIF

        ALLOCATE ( Ppt_rad_adj(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(Solrad_module, 'ppt_rad_adj', 'nhru,nmonths', 'real', &
     &       '0.02', '0.0', '0.5', &
     &       'Radiation reduced if HRU precipitation above this value', &
     &       'Monthly minimum precipitation, if HRU precipitation exceeds this value, radiation is'// &
     &       ' multiplied by radj_sppt or radj_wppt adjustment factor', &
     &       'inches')/=0 ) CALL read_error(1, 'ppt_rad_adj')
        ALLOCATE ( Radj_sppt(Nhru) )
        IF ( declparam(Solrad_module, 'radj_sppt', 'nhru', 'real', &
     &       '0.44', '0.0', '1.0', &
     &       'Adjustment to solar radiation on precipitation day - summer', &
     &       'Adjustment factor for computed solar radiation for summer day with greater than'// &
     &       ' ppt_rad_adj inches of precipitation for each HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'radj_sppt')
        ALLOCATE ( Radj_wppt(Nhru) )
        IF ( declparam(Solrad_module, 'radj_wppt', 'nhru', 'real', &
     &       '0.5', '0.0', '1.0', &
     &       'Adjustment to solar radiation on precipitation day - winter', &
     &       'Adjustment factor for computed solar radiation for winter day with greater than'// &
     &       ' ppt_rad_adj inches of precipitation for each HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'radj_wppt')
        ALLOCATE ( Radmax(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(Solrad_module, 'radmax', 'nhru,nmonths', 'real', &
     &       '0.8', '0.1', '1.0', &
     &       'Maximum fraction of potential solar radiation', &
     &       'Monthly (January to December) maximum fraction of the potential solar radiation'// &
     &       ' that may reach the ground due to haze, dust, smog, and so forth, for each HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'radmax')
      ENDIF

      ALLOCATE ( Epan_coef(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam('intcp', 'epan_coef', 'nhru,nmonths', 'real', &
     &     '1.0', '0.01', '3.0', &
     &     'Evaporation pan coefficient', &
     &     'Monthly (January to December) evaporation pan coefficient for each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'epan_coef')

      Use_pandata = OFF
      IF ( (Nevap>0 .AND. Et_flag==potet_pan_module) .OR. Model==DOCUMENTATION ) THEN
        Use_pandata = ACTIVE
        ALLOCATE ( Hru_pansta(Nhru) )
        IF ( declparam(Et_module, 'hru_pansta', 'nhru', 'integer', &
     &       '0', 'bounded', 'nevap', &
     &       'Index of pan evaporation station for each HRU', &
     &       'Index of pan evaporation station used to compute HRU potential ET', &
     &       'none')/=0 ) CALL read_error(1, 'hru_pansta')
      ENDIF

      ALLOCATE ( Sat_threshold(Nhru) )
      IF ( declparam(Soilzone_module, 'sat_threshold', 'nhru', 'real', &
     &     '999.0', '0.00001', '999.0', &
     &     'Soil saturation threshold, above field-capacity threshold', &
     &     'Water holding capacity of the gravity and preferential-'// &
     &     'flow reservoirs; difference between field capacity and'// &
     &     ' total soil saturation for each HRU', &
     &     'inches')/=0 ) CALL read_error(1, 'sat_threshold')

      ALLOCATE ( Soil_moist_max(Nhru) )
      IF ( declparam(Soilzone_module, 'soil_moist_max', 'nhru', 'real', &
     &     '2.0', '0.00001', '20.0', &
     &     'Maximum value of water for soil zone', &
     &     'Maximum available water holding capacity of capillary'// &
     &     ' reservoir from land surface to rooting depth of the'// &
     &     ' major vegetation type of each HRU', &
     &     'inches')/=0 ) CALL read_error(1, 'soil_moist_max')

      ALLOCATE ( Soil_rechr_max(Nhru) )
      IF ( PRMS4_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(Soilzone_module, 'soil_rechr_max', 'nhru', 'real', &
     &       '1.5', '0.00001', '20.0', &
     &       'Maximum storage for soil recharge zone', &
     &       'Maximum storage for soil recharge zone (upper portion of'// &
     &       ' capillary reservoir where losses occur as both'// &
     &       ' evaporation and transpiration); must be less than or equal to soil_moist_max', &
     &       'inches')/=0 ) CALL read_error(1, 'soil_rechr_max')
      ENDIF
      IF ( PRMS4_flag==OFF .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(Soilzone_module, 'soil_rechr_max_frac', 'nhru', 'real', &
     &       '1.0', '0.00001', '1.0', &
     &       'Fraction of capillary reservoir where losses occur as both evaporation and transpiration (soil recharge zone)', &
     &       'Fraction of the capillary reservoir water-holding capacity (soil_moist_max) where losses occur as both'// &
     &       ' evaporation and transpiration (upper zone of capillary reservoir) for each HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'soil_rechr_max_frac')
      ENDIF

      ALLOCATE ( Snowinfil_max(Nhru) )
      IF ( declparam(Srunoff_module, 'snowinfil_max', 'nhru', 'real', &
     &     '2.0', '0.0', '20.0', &
     &     'Maximum snow infiltration per day', &
     &     'Maximum snow infiltration per day for each HRU', &
     &     'inches/day')/=0 ) CALL read_error(1, 'snowinfil_max')

      ALLOCATE ( Imperv_stor_max(Nhru) )
      IF ( declparam(Srunoff_module, 'imperv_stor_max', 'nhru', 'real', &
     &     '0.05', '0.0', '0.5', &
     &     'HRU maximum impervious area retention storage', &
     &     'Maximum impervious area retention storage for each HRU', &
     &     'inches')/=0 ) CALL read_error(1, 'imperv_stor_max')

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==5 .OR. Model==DOCUMENTATION) THEN
        IF ( PRMS4_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'soil_rechr_init', 'nhru', 'real', &
     &         '1.0', '0.0', '20.0', &
     &         'Initial storage of water for soil recharge zone', &
     &         'Initial storage for soil recharge zone (upper part of'// &
     &         ' capillary reservoir where losses occur as both'// &
     &         ' evaporation and transpiration) for each HRU; must be'// &
     &         ' less than or equal to soil_moist_init', &
     &         'inches')/=0 ) CALL read_error(1, 'soil_rechr_init')
          IF ( declparam(MODNAME, 'soil_moist_init', 'nhru', 'real', &
     &         '3.0', '0.0', '20.0', &
     &         'Initial value of available water in capillary reservoir', &
     &         'Initial value of available water in capillary reservoir for each HRU', &
     &         'inches')/=0 ) CALL read_error(1, 'soil_moist_init')
          IF ( declparam(MODNAME, 'ssstor_init', 'nssr', 'real', &
     &         '0.0', '0.0', '10.0', &
     &         'Initial storage in each GVR and PFR', &
     &         'Initial storage of the gravity and preferential-flow reservoirs for each HRU', &
     &         'inches')/=0 ) CALL read_error(1, 'ssstor_init')
        ENDIF
        IF ( PRMS4_flag==OFF .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(Soilzone_module, 'soil_rechr_init_frac', 'nhru', 'real', &
     &         '0.0', '0.0', '1.0', &
     &         'Initial fraction of available water in the soil recharge zone within the capillary reservoir', &
     &         'Initial fraction of available water in the capillary reservoir where losses occur'// &
     &         ' as both evaporation and transpiration (upper zone of capillary reservoir) for each HRU', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'soil_rechr_init_frac')
          IF ( declparam(Soilzone_module, 'soil_moist_init_frac', 'nhru', 'real', &
     &         '0.0', '0.0', '1.0', &
     &         'Initial fraction available water in the capillary reservoir', &
     &         'Initial fraction of available water in the capillary reservoir (fraction of soil_moist_max for each HRU', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'soil_moist_init_frac')
          IF ( declparam(Soilzone_module, 'ssstor_init_frac', 'nssr', 'real', &
     &         '0.0', '0.0', '1.0', &
     &         'Initial fraction of available water in the gravity plus preferential-flow reservoirs', &
     &         'Initial fraction of available water in the gravity plus preferential-flow reservoirs'// &
     &         ' (fraction of sat_threshold) for each HRU', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'ssstor_init_frac')
        ENDIF
      ENDIF

      END FUNCTION climateflow_decl

!***********************************************************************
!     climateflow_init - Initialize module - get parameter values,
!                        set initial values and check parameter values
!***********************************************************************
      INTEGER FUNCTION climateflow_init()
      USE PRMS_CLIMATEVARS
      USE PRMS_FLOWVARS
      USE PRMS_BASIN, ONLY: Elev_units, Active_hrus, Hru_route_order, Hru_type
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: checkdim_param_limits, checkdim_bounded_limits
      REAL, EXTERNAL :: c_to_f, f_to_c
! Local variables
      INTEGER :: i, j, ierr
!***********************************************************************
      climateflow_init = 0

      IF ( Temp_flag<climate_hru_module ) THEN
        IF ( getparam(Temp_module, 'tsta_elev', Ntemp, 'real', Tsta_elev)/=0 ) CALL read_error(2, 'tsta_elev')
        DO i = 1, Ntemp
          IF ( Elev_units==FEET ) THEN
            Tsta_elev_feet(i) = Tsta_elev(i)
            Tsta_elev_meters(i) = Tsta_elev_feet(i)*FEET2METERS
          ELSE
            Tsta_elev_meters(i) = Tsta_elev(i)
            Tsta_elev_feet(i) = Tsta_elev_meters(i)*METERS2FEET
          ENDIF
        ENDDO
      ENDIF

      IF ( getparam('snowcomp', 'potet_sublim', Nhru, 'real', Potet_sublim)/=0 ) CALL read_error(2, 'potet_sublim')

      IF ( Temp_flag==temp_1sta_module .OR. Temp_flag==temp_laps_module .OR. Temp_flag==temp_dist2_module .OR. &
     &     Temp_flag==ide_dist_module .OR. Temp_flag==xyz_dist_module .OR. Temp_flag==temp_sta_module ) THEN
        IF ( getparam(Temp_module, 'tmax_adj', Nhru*MONTHS_PER_YEAR, 'real', Tmax_aspect_adjust)/=0 ) CALL read_error(2, 'tmax_adj')
        IF ( getparam(Temp_module, 'tmin_adj', Nhru*MONTHS_PER_YEAR, 'real', Tmin_aspect_adjust)/=0 ) CALL read_error(2, 'tmin_adj')
      ENDIF

      IF ( getparam(Temp_module, 'temp_units', 1, 'integer', Temp_units)/=0 ) CALL read_error(2, 'temp_units')

      IF ( Temp_flag<ide_dist_module ) THEN
        IF ( getparam(Temp_module, 'basin_tsta', 1, 'integer', Basin_tsta)/=0 ) CALL read_error(2, 'basin_tsta')
        CALL checkdim_param_limits(1, 'basin_tsta', 'ntemp', Basin_tsta, 1, Ntemp, Inputerror_flag)
      ELSE
        Basin_tsta = 0
      ENDIF

      IF ( Temp_flag==temp_1sta_module .OR. Temp_flag==temp_laps_module .OR. Temp_flag==temp_sta_module ) THEN
        IF ( getparam(Temp_module, 'hru_tsta', Nhru, 'integer', Hru_tsta)/=0 ) CALL read_error(2, 'hru_tsta')
        IF ( Parameter_check_flag>0 ) &
     &       CALL checkdim_bounded_limits('hru_tsta', 'ntemp', Hru_tsta, Nhru, 0, Ntemp, Inputerror_flag)
      ENDIF

      IF ( getparam(Precip_module, 'tmax_allsnow', Nhru*MONTHS_PER_YEAR, 'real', Tmax_allsnow)/=0 ) &
     &    CALL read_error(2, 'tmax_allsnow')

      IF ( PRMS4_flag==ACTIVE ) THEN
        IF ( getparam(Precip_module, 'tmax_allrain', Nhru*MONTHS_PER_YEAR, 'real', Tmax_allrain)/=0 ) &
     &       CALL read_error(2, 'tmax_allrain')
        DO j = 1, MONTHS_PER_YEAR
          DO i = 1, Nhru
            Tmax_allrain_offset(i, j) = Tmax_allrain(i, j) - Tmax_allsnow(i, j)
            IF ( Tmax_allrain_offset(i, j)<0.0 ) THEN
              IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, tmax_allsnow > tmax_allrain for HRU:', i, '; month:', j, &
     &                                               ' tmax_allrain set to tmax_allsnow'
              Tmax_allrain_offset(i, j) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ELSE
        IF ( getparam(Precip_module, 'tmax_allrain_offset', Nhru*MONTHS_PER_YEAR, 'real', Tmax_allrain_offset)/=0 ) &
     &                CALL read_error(2, 'tmax_allrain_offset')
      ENDIF

      ! Set tmax_allrain in units of the input values
      ! tmax_allsnow must be in the units of the input values
      IF ( Temp_units==FAHRENHEIT ) THEN
        Tmax_allsnow_f = Tmax_allsnow
        DO j = 1, MONTHS_PER_YEAR
          DO i = 1, Nhru
            Tmax_allrain_f(i, j) = Tmax_allsnow(i, j) + Tmax_allrain_offset(i, j)
            Tmax_allsnow_c(i, j) = f_to_c(Tmax_allsnow(i,j))
          ENDDO
        ENDDO
        Tmax_allrain = Tmax_allrain_f
      ELSE
        Tmax_allsnow_c = Tmax_allsnow
        DO i = 1, MONTHS_PER_YEAR
          DO j = 1, Nhru
            Tmax_allsnow_f(j, i) = c_to_f(Tmax_allsnow(j,i))
            Tmax_allrain(j, i) = Tmax_allsnow(j, i) + Tmax_allrain_offset(j, i)
            Tmax_allrain_f(j, i) = c_to_f(Tmax_allrain(j, i))
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE ( Tmax_allrain_offset )

      IF ( getparam(Precip_module, 'adjmix_rain', Nhru*MONTHS_PER_YEAR, 'real', Adjmix_rain)/=0 ) CALL read_error(2, 'adjmix_rain')

      IF ( getparam(Precip_module, 'precip_units', 1, 'integer', Precip_units)/=0 ) CALL read_error(2, 'precip_units')

      IF ( getparam(Precip_module, 'ppt_zero_thresh', 1, 'real', Ppt_zero_thresh)/=0 ) CALL read_error(2, 'ppt_zero_thresh')

      IF ( Precip_flag==precip_laps_module .OR. Precip_flag==xyz_dist_module .OR. Precip_flag==ide_dist_module ) THEN
        IF ( getparam(Precip_module, 'psta_elev', Nrain, 'real', Psta_elev)/=0 ) CALL read_error(2, 'psta_elev')
        DO i = 1, Nrain
          IF ( Elev_units==FEET ) THEN
            Psta_elev_feet(i) = Psta_elev(i)
            Psta_elev_meters(i) = Psta_elev_feet(i)*FEET2METERS
          ELSE
            Psta_elev_meters(i) = Psta_elev(i)
            Psta_elev_feet(i) = Psta_elev_meters(i)*METERS2FEET
          ENDIF
        ENDDO
      ENDIF

      IF ( Solrad_flag==ddsolrad_module .OR. Solrad_flag==ccsolrad_module ) THEN
        Solsta_flag = OFF
        IF ( Nsol>0 ) THEN
          IF ( getparam(Solrad_module, 'basin_solsta', 1, 'integer', Basin_solsta)/=0 ) CALL read_error(2, 'basin_solsta')

          IF ( getparam(Solrad_module, 'rad_conv', 1, 'real', Rad_conv)/=0 ) CALL read_error(2, 'rad_conv')

          IF ( getparam(Solrad_module, 'hru_solsta', Nhru, 'integer', Hru_solsta)/=0 ) CALL read_error(2, 'hru_solsta')

          IF ( Parameter_check_flag>0 ) THEN
            CALL checkdim_param_limits(1, 'basin_solsta', 'nsol', Basin_solsta, 1, Nsol, Inputerror_flag)
            CALL checkdim_bounded_limits('hru_solsta', 'nsol', Hru_solsta, Nhru, 0, Nsol, Inputerror_flag)
          ENDIF

          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            IF ( Hru_solsta(i)>0 ) THEN
              Solsta_flag = ACTIVE
              EXIT
            ENDIF
          ENDDO
        ENDIF
        IF ( getparam(Solrad_module, 'radj_sppt', Nhru, 'real', Radj_sppt)/=0 ) CALL read_error(2, 'radj_sppt')
        IF ( getparam(Solrad_module, 'radj_wppt', Nhru, 'real', Radj_wppt)/=0 ) CALL read_error(2, 'radj_wppt')
        IF ( getparam(Solrad_module, 'ppt_rad_adj', Nhru*MONTHS_PER_YEAR, 'real', Ppt_rad_adj)/=0 ) &
     &       CALL read_error(2, 'ppt_rad_adj')
        IF ( getparam(Solrad_module, 'radmax', Nhru*MONTHS_PER_YEAR, 'real', Radmax)/=0 ) CALL read_error(2, 'radmax')
      ELSE
        Basin_solsta = 0
        Rad_conv = 1.0
      ENDIF

      IF ( Use_pandata==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'hru_pansta', Nhru, 'integer', Hru_pansta)/=0 ) CALL read_error(2, 'hru_pansta')
        IF ( Parameter_check_flag>0 ) &
     &       CALL checkdim_bounded_limits('hru_pansta', 'nevap', Hru_pansta, Nhru, 0, Nevap, Inputerror_flag)
      ENDIF

      IF ( getparam('intcp', 'epan_coef', Nhru*MONTHS_PER_YEAR, 'real', Epan_coef)/=0 ) CALL read_error(2, 'epan_coef')

! FLOW VARIABLES AND PARAMETERS
      IF ( getparam(Soilzone_module, 'sat_threshold', Nhru, 'real', Sat_threshold)/=0 ) CALL read_error(2, 'sat_threshold')

      IF ( getparam(Soilzone_module, 'soil_moist_max', Nhru, 'real', Soil_moist_max)/=0 ) CALL read_error(2, 'soil_moist_max')

      IF ( PRMS4_flag==ACTIVE ) THEN
        IF ( getparam(Soilzone_module, 'soil_rechr_max', Nhru, 'real', Soil_rechr_max)/=0 ) CALL read_error(2, 'soil_rechr_max')
      ELSE
        IF ( getparam(Soilzone_module, 'soil_rechr_max_frac', Nhru, 'real', Soil_rechr_max)/=0 ) &
     &       CALL read_error(2, 'soil_rechr_max_frac')
        Soil_rechr_max = Soil_rechr_max*Soil_moist_max
      ENDIF

      ierr = 0
      IF ( Init_vars_from_file==OFF .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==5 ) THEN
        IF ( PRMS4_flag==ACTIVE ) THEN
          ! use PRMS4 parameters
          IF ( getparam(Soilzone_module, 'soil_moist_init', Nhru, 'real', Soil_moist)/=0 ) &
     &         CALL read_error(2, 'soil_moist_init')
          IF ( getparam(Soilzone_module, 'soil_rechr_init', Nhru, 'real', Soil_rechr)/=0 ) &
     &         CALL read_error(2, 'soil_rechr_init')
          IF ( getparam(Soilzone_module, 'ssstor_init', Nssr, 'real', Ssres_stor)/=0 ) &
     &         CALL read_error(2, 'ssstor_init')
! PRMS 5 parameters
        ELSE
          IF ( getparam(Soilzone_module, 'soil_moist_init_frac', Nhru, 'real', Soil_moist)/=0 ) &
     &         CALL read_error(2, 'soil_moist_init_frac')
          IF ( getparam(Soilzone_module, 'soil_rechr_init_frac', Nhru, 'real', Soil_rechr)/=0 ) &
     &         CALL read_error(2, 'soil_rechr_init_frac')
          IF ( getparam(Soilzone_module, 'ssstor_init_frac', Nssr, 'real', Ssres_stor)/=0 ) &
     &         CALL read_error(2, 'ssstor_init_frac')
          Soil_rechr = Soil_rechr*Soil_rechr_max
          Soil_moist = Soil_moist*Soil_moist_max
          Ssres_stor = Ssres_stor*Sat_threshold
        ENDIF
        Slow_stor = Ssres_stor
      ENDIF

      ! check parameters
      DO i = 1, Nhru
        IF ( Hru_type(i)==INACTIVE .OR. Hru_type(i)==LAKE ) CYCLE
        ! hru_type = land or swale or glacier
        IF ( Soil_moist_max(i)<0.00001 ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9006, i, Soil_moist_max(i)
            ierr = 1
          ELSE
            Soil_moist_max(i) = 0.00001
            IF ( Print_debug>DEBUG_less ) PRINT 9008, i
          ENDIF
        ENDIF
        IF ( Soil_rechr_max(i)<0.00001 ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9007, i, Soil_rechr_max(i)
            ierr = 1
          ELSE
            Soil_rechr_max(i) = 0.00001
            IF ( Print_debug>DEBUG_less ) PRINT 9009, i
          ENDIF
        ENDIF
        IF ( Soil_rechr_max(i)>Soil_moist_max(i) ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9002, i, Soil_rechr_max(i), Soil_moist_max(i)
            ierr = 1
          ELSE
            IF ( Print_debug>DEBUG_less ) PRINT 9012, i, Soil_rechr_max(i), Soil_moist_max(i)
            Soil_rechr_max(i) = Soil_moist_max(i)
          ENDIF
        ENDIF
        IF ( Soil_rechr(i)>Soil_rechr_max(i) ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9003, i, Soil_rechr(i), Soil_rechr_max(i)
            ierr = 1
          ELSE
            IF ( Print_debug>DEBUG_less ) PRINT 9013, i, Soil_rechr(i), Soil_rechr_max(i)
            Soil_rechr(i) = Soil_rechr_max(i)
          ENDIF
        ENDIF
        IF ( Soil_moist(i)>Soil_moist_max(i) ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9004, i, Soil_moist(i), Soil_moist_max(i)
            ierr = 1
          ELSE
            IF ( Print_debug>DEBUG_less ) PRINT 9014, i, Soil_moist(i), Soil_moist_max(i)
            Soil_moist(i) = Soil_moist_max(i)
          ENDIF
        ENDIF
        IF ( Soil_rechr(i)>Soil_moist(i) ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT 9005, i, Soil_rechr(i), Soil_moist(i)
            ierr = 1
          ELSE
            IF ( Print_debug>DEBUG_less ) PRINT 9015, i, Soil_rechr(i), Soil_moist(i)
            Soil_rechr(i) = Soil_moist(i)
          ENDIF
        ENDIF
        IF ( Ssres_stor(i)>Sat_threshold(i) ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT *, 'ERROR, HRU:', i, Ssres_stor(i), Sat_threshold(i), ' ssres_stor > sat_threshold'
            ierr = 1
          ELSE
            PRINT *, 'WARNING, HRU:', i, Ssres_stor(i), Sat_threshold(i), ' ssres_stor > sat_threshold, ssres_stor set to max'
            Ssres_stor(i) = Sat_threshold(i)
          ENDIF
        ENDIF
      ENDDO

      IF ( ierr>0 ) STOP ERROR_PARAM

      IF ( getparam(Srunoff_module, 'snowinfil_max', Nhru, 'real', Snowinfil_max)/=0 ) CALL read_error(2, 'snowinfil_max')

      IF ( getparam(Srunoff_module, 'imperv_stor_max', Nhru, 'real', Imperv_stor_max)/=0 ) CALL read_error(2, 'imperv_stor_max')

! initialize arrays (dimensioned Nhru)
      Tmaxf = 0.0
      Tminf = 0.0
      Tavgf = 0.0
      Tmaxc = 0.0
      Tminc = 0.0
      Tavgc = 0.0
      Tmax_hru = 0.0
      Tmin_hru = 0.0
      Pptmix = OFF
      Newsnow = OFF
      Prmx = 0.0
      Hru_ppt = 0.0
      Hru_rain = 0.0
      Hru_snow = 0.0
      Swrad = 0.0
      Potet = 0.0
      Slow_flow = 0.0
      Soil_to_gw = 0.0
      Soil_to_ssr = 0.0
      Hru_actet = 0.0
      Infil = 0.0
      Sroff = 0.0
! initialize arrays (dimensioned Nssr)
      Ssr_to_gw = 0.0
      Ssres_in = 0.0
      Ssres_flow = 0.0
      IF ( Solrad_flag==ddsolrad_module .OR. Solrad_flag==ccsolrad_module ) Orad_hru = 0.0
      IF ( Et_flag==potet_pt_module .OR. Et_flag==potet_pm_module .OR. Et_flag==potet_pm_sta_module ) THEN
        Tempc_dewpt = 0.0
        Vp_actual = 0.0
        Lwrad_net = 0.0
        Vp_slope = 0.0
        IF ( Et_flag==potet_pm_module .OR. Et_flag==potet_pm_sta_module ) Vp_sat = 0.0
        IF ( Humidity_cbh_flag==OFF ) THEN
          IF ( getparam(Et_module, 'humidity_percent', Nhru*12, 'real', Humidity_percent)/=0 ) &
     &         CALL read_error(2, 'humidity_percent')
        ELSE
          Humidity_percent = 1.0
        ENDIF
      ENDIF
! initialize arrays (dimensioned Nsegment)
      IF ( Stream_order_flag==ACTIVE ) THEN
        Seg_upstream_inflow = 0.0D0
        Seg_lateral_inflow = 0.0D0
      ENDIF

! initialize scalers
      Basin_temp = 0.0D0
      Basin_tmax = 0.0D0
      Basin_tmin = 0.0D0
      Basin_ppt = 0.0D0
      Basin_obs_ppt = 0.0D0
      Basin_rain = 0.0D0
      Basin_snow = 0.0D0
      Basin_horad = 0.0D0
      Basin_potsw = 0.0D0
      Basin_swrad = 0.0D0
      Basin_potet = 0.0D0
      Basin_humidity = 0.0D0
      Basin_orad = 0.0D0
      Basin_perv_et = 0.0D0
      Basin_actet = 0.0D0
      Basin_lakeevap = 0.0D0
      Basin_swale_et = 0.0D0
      Basin_soil_to_gw = 0.0D0
      Basin_ssflow = 0.0D0
      Basin_sroff = 0.0D0
      Solrad_tmax = 0.0
      Solrad_tmin = 0.0
      Basin_cfs = 0.0D0
      Basin_cms = 0.0D0
      Basin_stflow_in = 0.0D0
      Basin_stflow_out = 0.0D0
      Basin_ssflow_cfs = 0.0D0
      Basin_sroff_cfs = 0.0D0
      Basin_gwflow_cfs = 0.0D0
      Flow_out = 0.0D0
      Orad = 0.0

      IF ( Init_vars_from_file>0 .OR. ierr>0 ) RETURN

      Basin_soil_moist = 0.0D0
      Basin_ssstor = 0.0D0
      Basin_lake_stor = 0.0D0
      Basin_transp_on = OFF
! initialize arrays (dimensioned Nsegment)
      IF ( Stream_order_flag==ACTIVE ) THEN
        Seg_inflow = 0.0D0
        Seg_outflow = 0.0D0
      ENDIF
      Transp_on = OFF
! initialize storage variables
      Imperv_stor = 0.0
      Pkwater_equiv = 0.0D0
      Slow_stor = 0.0
      IF ( GSFLOW_flag==OFF ) Gwres_stor = 0.0D0 ! not needed for GSFLOW
      IF ( Dprst_flag==ACTIVE ) THEN
        Dprst_vol_open = 0.0D0
        Dprst_vol_clos = 0.0D0
      ENDIF
! initialize arrays (dimensioned nlake)
      IF ( Nlake>0 ) Lake_vol = 0.0D0

 9002 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_rechr_max > soil_moist_max', 2F10.5)
 9003 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_rechr_init > soil_rechr_max', 2F10.5)
 9004 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_moist_init > soil_moist_max', 2F10.5)
 9005 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_rechr > soil_moist based on init and max values', 2F10.5)
 9006 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_moist_max < 0.00001', F10.5)
 9007 FORMAT (/, 'ERROR, HRU: ', I0, ' soil_rechr_max < 0.00001', F10.5)
 9008 FORMAT (/, 'WARNING, HRU: ', I0, ' soil_moist_max < 0.00001, set to 0.00001')
 9009 FORMAT (/, 'WARNING, HRU: ', I0, ' soil_rechr_max < 0.00001, set to 0.00001')
 9012 FORMAT ('WARNING, HRU: ', I0, ' soil_rechr_max > soil_moist_max,', 2F10.5, /, 9X, &
     &        'soil_rechr_max set to soil_moist_max')
 9013 FORMAT ('WARNING, HRU: ', I0, ' soil_rechr_init > soil_rechr_max,', 2F10.5, /, 9X, &
     &        'soil_rechr set to soil_rechr_max')
 9014 FORMAT ('WARNING, HRU: ', I0, ' soil_moist_init > soil_moist_max,', 2F10.5, /, 9X, &
     &        'soil_moist set to soil_moist_max')
 9015 FORMAT ('WARNING, HRU: ', I0, ' soil_rechr_init > soil_moist_init,', 2F10.5, /, 9X, &
     &        'soil_rechr set to soil_moist based on init and max values')

      END FUNCTION climateflow_init

!***********************************************************************
!     Sets temperatures in both system of units for each HRU
!***********************************************************************
      SUBROUTINE temp_set(Ihru, Tmax, Tmin, Tmaxf, Tminf, Tavgf, Tmaxc, Tminc, Tavgc, Hru_area)
      USE PRMS_CLIMATEVARS, ONLY: Basin_temp, Basin_tmax, Basin_tmin, Temp_units, Tmax_hru, Tmin_hru
      USE PRMS_CONSTANTS, ONLY: MINTEMP, MAXTEMP, ERROR_temp
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Ihru
      REAL, INTENT(IN) :: Tmax, Tmin, Hru_area
      REAL, INTENT(OUT) :: Tmaxf, Tminf, Tavgf, Tmaxc, Tminc, Tavgc
! Functions
      INTRINSIC :: DBLE
      REAL, EXTERNAL :: c_to_f, f_to_c
      EXTERNAL :: print_date
!***********************************************************************
      IF ( Temp_units==0 ) THEN
!       degrees Fahrenheit
        Tmaxf = Tmax
        Tminf = Tmin
        Tavgf = (Tmax+Tmin)*0.5
        Tmaxc = f_to_c(Tmax)
        Tminc = f_to_c(Tmin)
        Tavgc = f_to_c(Tavgf)
        Basin_temp = Basin_temp + DBLE( Tavgf*Hru_area )
      ELSE
!       degrees Celsius
        Tmaxc = Tmax
        Tminc = Tmin
        Tavgc = (Tmax+Tmin)*0.5
        Tmaxf = c_to_f(Tmax)
        Tminf = c_to_f(Tmin)
        Tavgf = c_to_f(Tavgc)
        Basin_temp = Basin_temp + DBLE( Tavgc*Hru_area )
      ENDIF

      IF ( Tminf<MINTEMP .OR. Tmaxf>MAXTEMP ) THEN
        PRINT '(A,I0,1X,F0.4,1X,F0.4,/)', ' ERROR, invalid temperature value for HRU: ', Ihru, Tminf, Tmaxf
        CALL print_date(1)
        ERROR STOP ERROR_temp
      ENDIF
      Tmax_hru(Ihru) = Tmax ! in units temp_units
      Tmin_hru(Ihru) = Tmin ! in units temp_units

      Basin_tmax = Basin_tmax + DBLE( Tmax*Hru_area )
      Basin_tmin = Basin_tmin + DBLE( Tmin*Hru_area )

      END SUBROUTINE temp_set

!***********************************************************************
!     Computes precipitation form (rain, snow or mix) and depth for each HRU
!***********************************************************************
      SUBROUTINE precip_form(Precip, Hru_ppt, Hru_rain, Hru_snow, Tmaxf, &
     &           Tminf, Pptmix, Newsnow, Prmx, Tmax_allrain_f, Rain_adj, &
     &           Snow_adj, Adjmix_rain, Hru_area, Sum_obs, Tmax_allsnow_f)
      USE PRMS_CONSTANTS, ONLY: NEARZERO
      USE PRMS_CLIMATEVARS, ONLY: Basin_ppt, Basin_rain, Basin_snow
      IMPLICIT NONE
! Functions
      INTRINSIC :: ABS, DBLE
      EXTERNAL :: print_date
! Arguments
      REAL, INTENT(IN) :: Tmax_allrain_f, Tmax_allsnow_f, Rain_adj, Snow_adj
      REAL, INTENT(IN) :: Adjmix_rain, Tmaxf, Tminf, Hru_area
      DOUBLE PRECISION, INTENT(INOUT) :: Sum_obs
      INTEGER, INTENT(INOUT) :: Pptmix, Newsnow
      REAL, INTENT(INOUT) :: Precip, Hru_rain, Hru_snow, Prmx, Hru_ppt
! Local Variables
      REAL :: tdiff
!***********************************************************************
      ! basin precipitation before adjustments
      Sum_obs = Sum_obs + DBLE( Precip*Hru_area )

!******If maximum temperature is below or equal to the base temperature
!******for snow then precipitation is all snow
      IF ( Tmaxf<=Tmax_allsnow_f ) THEN
        Hru_ppt = Precip*Snow_adj
        Hru_snow = Hru_ppt
        Newsnow = 1

!******If minimum temperature is above base temperature for snow or
!******maximum temperature is above all_rain temperature then
!******precipitation is all rain
      ELSEIF ( Tminf>Tmax_allsnow_f .OR. Tmaxf>=Tmax_allrain_f ) THEN
        Hru_ppt = Precip*Rain_adj
        Hru_rain = Hru_ppt
        Prmx = 1.0

!******Otherwise precipitation is a mixture of rain and snow
      ELSE
        tdiff = Tmaxf - Tminf
        IF ( tdiff<0.0 ) THEN
          PRINT *, 'ERROR, tmax < tmin (degrees Fahrenheit), tmax:', Tmaxf, ' tmin:', TminF
          CALL print_date(1)
        ENDIF
        IF ( ABS(tdiff)<NEARZERO ) tdiff = 0.0001
        Prmx = ((Tmaxf-Tmax_allsnow_f)/tdiff)*Adjmix_rain
        IF ( Prmx<0.0 ) Prmx = 0.0

!******Unless mixture adjustment raises the proportion of rain to
!******greater than or equal to 1.0 in which case it all rain
!******If not, it is a rain/snow mixture
        IF ( Prmx<1.0 ) THEN
          Pptmix = 1
          Hru_ppt = Precip*Snow_adj
          Hru_rain = Prmx*Hru_ppt
          Hru_snow = Hru_ppt - Hru_rain
          Newsnow = 1
        ELSE
          Hru_ppt = Precip*Rain_adj
          Hru_rain = Hru_ppt
          Prmx = 1.0
        ENDIF
      ENDIF
      Basin_ppt = Basin_ppt + DBLE( Hru_ppt*Hru_area )
      Basin_rain = Basin_rain + DBLE( Hru_rain*Hru_area )
      Basin_snow = Basin_snow + DBLE( Hru_snow*Hru_area )

      END SUBROUTINE precip_form

!***********************************************************************
!     Write or read restart file
!***********************************************************************
      SUBROUTINE climateflow_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_CLIMATEVARS
      USE PRMS_FLOWVARS
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=11) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_transp_on, Basin_soil_moist, Basin_ssstor, Basin_lake_stor
        WRITE ( Restart_outunit ) Transp_on
        WRITE ( Restart_outunit ) Pkwater_equiv
        IF ( Glacier_flag==ACTIVE ) THEN
          WRITE ( Restart_outunit) Glacier_frac
          WRITE ( Restart_outunit) Glrette_frac
          WRITE ( Restart_outunit) Alt_above_ela
        ENDIF
        WRITE ( Restart_outunit ) Soil_moist
        WRITE ( Restart_outunit ) Slow_stor
        WRITE ( Restart_outunit ) Ssres_stor
        WRITE ( Restart_outunit ) Soil_rechr
        WRITE ( Restart_outunit ) Imperv_stor
        IF ( GSFLOW_flag==OFF ) WRITE ( Restart_outunit ) Gwres_stor
        IF ( Dprst_flag==ACTIVE ) THEN
          WRITE ( Restart_outunit ) Dprst_vol_open
          WRITE ( Restart_outunit ) Dprst_vol_clos
        ENDIF
        IF ( Stream_order_flag==ACTIVE ) THEN
          WRITE ( Restart_outunit ) Seg_inflow
          WRITE ( Restart_outunit ) Seg_outflow
        ENDIF
        IF ( Nlake>0 ) WRITE ( Restart_outunit ) Lake_vol
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_transp_on, Basin_soil_moist, Basin_ssstor, Basin_lake_stor
        READ ( Restart_inunit ) Transp_on
        READ ( Restart_inunit ) Pkwater_equiv
        IF ( Glacier_flag==ACTIVE ) THEN
          READ ( Restart_inunit) Glacier_frac
          READ ( Restart_inunit) Glrette_frac
          READ ( Restart_inunit) Alt_above_ela
        ENDIF
        READ ( Restart_inunit ) Soil_moist
        READ ( Restart_inunit ) Slow_stor
        READ ( Restart_inunit ) Ssres_stor
        READ ( Restart_inunit ) Soil_rechr
        READ ( Restart_inunit ) Imperv_stor
        IF ( GSFLOW_flag==OFF ) READ ( Restart_inunit ) Gwres_stor
        IF ( Dprst_flag==ACTIVE ) THEN
          READ ( Restart_inunit ) Dprst_vol_open
          READ ( Restart_inunit ) Dprst_vol_clos
        ENDIF
        IF ( Stream_order_flag==ACTIVE ) THEN
          READ ( Restart_inunit ) Seg_inflow
          READ ( Restart_inunit ) Seg_outflow
        ENDIF
        IF ( Nlake>0 ) READ ( Restart_inunit ) Lake_vol
      ENDIF
      END SUBROUTINE climateflow_restart
