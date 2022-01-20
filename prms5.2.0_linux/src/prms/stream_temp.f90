!***********************************************************************
! stream temperature module
!***********************************************************************
      MODULE PRMS_STRMTEMP
      USE PRMS_CONSTANTS, ONLY: MAX_DAYS_PER_YEAR, MONTHS_PER_YEAR, DOCUMENTATION, ACTIVE, OFF, &
     &    NEARZERO, ERROR_param, CFS2CMS_CONV, DAYS_YR, DAYS_PER_YEAR, DAYS_YR
      USE PRMS_MODULE, ONLY: Process_flag, Nsegment, Model, Init_vars_from_file, &
     &    Print_debug, Strmtemp_humidity_flag, Model, Inputerror_flag
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Stream Temperature'
      character(len=11), parameter :: MODNAME = 'stream_temp'
      character(len=*), parameter :: Version_stream_temp = '2020-12-02'
      INTEGER, SAVE, ALLOCATABLE :: Seg_hru_count(:), Seg_close(:)
      REAL, SAVE, ALLOCATABLE ::  seg_tave_ss(:), Seg_carea_inv(:), seg_tave_sroff(:), seg_tave_lat(:)
      REAL, SAVE, ALLOCATABLE :: seg_tave_gw(:), Flowsum(:)

      ! next variables only needed if strm_temp_shade_flag = 0
      REAL, SAVE, ALLOCATABLE :: Shade_jday(:, :), Svi_jday(:, :)
      REAL, SAVE, ALLOCATABLE :: Seg_lat(:), Seg_elev(:)
      REAL, SAVE, ALLOCATABLE :: Press(:)
      REAL, SAVE, ALLOCATABLE :: Cos_seg_lat(:), Sin_seg_lat(:), Horizontal_hour_angle(:, :), Total_shade(:, :)
      REAL, SAVE, ALLOCATABLE :: Sin_declination(:, :), Sin_lat_decl(:, :), Cos_lat_decl(:, :), Sin_alrs(:, :)
      REAL, SAVE, ALLOCATABLE :: Max_solar_altitude(:, :), Level_sunset_azimuth(:, :)
      REAL, SAVE, ALLOCATABLE :: Local_sunset_hour_angle(:, :), Local_sunrise_hour_angle(:, :)
      REAL, SAVE, ALLOCATABLE :: gw_sum(:), ss_sum(:)
      REAL, SAVE, ALLOCATABLE ::  gw_silo(:,:), ss_silo(:,:)
      REAL, SAVE, ALLOCATABLE :: hru_area_sum(:)
      INTEGER, SAVE, ALLOCATABLE :: upstream_count(:)
      INTEGER, SAVE, ALLOCATABLE :: upstream_idx(:,:)
      INTEGER, SAVE ::  gw_index, ss_index

!   Declared Variables
      REAL, SAVE, ALLOCATABLE :: Seg_tave_water(:), seg_tave_upstream(:), Seg_daylight(:)
      REAL, SAVE, ALLOCATABLE :: Seg_humid(:), Seg_width(:), Seg_ccov(:), seg_shade(:)
      REAL, SAVE, ALLOCATABLE :: Seg_tave_air(:), Seg_melt(:), Seg_rain(:)
      DOUBLE PRECISION, ALLOCATABLE :: Seg_potet(:)
!   Segment Parameters
      REAL, SAVE, ALLOCATABLE :: Seg_length(:) !, Mann_n(:)
      REAL, SAVE, ALLOCATABLE :: Seg_slope(:), Width_values(:, :)
      REAL, SAVE, ALLOCATABLE :: width_alpha(:), width_m(:), Stream_tave_init(:)
      INTEGER, SAVE:: Width_dim, Maxiter_sntemp
      REAL, SAVE, ALLOCATABLE :: Seg_humidity(:, :)
      REAL, SAVE, ALLOCATABLE :: lat_temp_adj(:, :)
      INTEGER, SAVE, ALLOCATABLE :: Seg_humidity_sta(:)
!   Shade Parameters needed if stream_temp_shade_flag = 0
      REAL, SAVE, ALLOCATABLE :: Azrh(:), Alte(:), Altw(:), Vce(:)
      REAL, SAVE, ALLOCATABLE :: Vdemx(:), Vhe(:), Voe(:), Vcw(:), Vdwmx(:), Vhw(:), Vow(:)
      REAL, SAVE, ALLOCATABLE :: Vdemn(:), Vdwmn(:)
      INTEGER, SAVE :: Spring_jday, Summer_jday, Autumn_jday, Winter_jday
!   Shade Parameters needed if stream_temp_shade_flag = 2
      REAL, SAVE, ALLOCATABLE :: Segshade_sum(:), Segshade_win(:)
      REAL, SAVE:: Albedo, Melt_temp
      ! INTEGER, SAVE :: Shadeflg, now using stream_temp_shade_flag
      INTEGER, SAVE, ALLOCATABLE :: Ss_tau(:), Gw_tau(:)
!   Control parameters
      INTEGER, SAVE :: Stream_temp_shade_flag
!   Conversions
      INTRINSIC :: ACOS
      REAL, PARAMETER :: HALF_PI = ACOS(0.0), ZERO_C = 273.16
      REAL, PARAMETER :: PI = ACOS(-1.0)
      REAL, PARAMETER :: DEG_TO_RAD = PI / 180.0, NOFLOW_TEMP = -98.9
      DOUBLE PRECISION :: MPS_CONVERT = 2.93981481D-07
      END MODULE PRMS_STRMTEMP

!***********************************************************************
!     Main stream temperature routine
!***********************************************************************
      INTEGER FUNCTION stream_temp()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: stream_temp_decl, stream_temp_init, stream_temp_run, stream_temp_setdims
      EXTERNAL :: stream_temp_restart
!***********************************************************************
      stream_temp = 0

      IF ( Process_flag==RUN ) THEN
         stream_temp  = stream_temp_run()
      ELSEIF ( Process_flag==DECL ) THEN
         stream_temp  = stream_temp_decl()
      ELSEIF ( Process_flag==INIT ) THEN
         IF ( Init_vars_from_file>0 ) CALL stream_temp_restart(1)
         stream_temp  = stream_temp_init()
      ELSEIF ( Process_flag==CLEAN ) THEN
         IF ( Save_vars_to_file==ACTIVE ) CALL stream_temp_restart(0)
      ENDIF

      END FUNCTION stream_temp

!***********************************************************************
!     stream_temp_decl - set up parameters and storage
!   Declared Parameters
!***********************************************************************
      INTEGER FUNCTION stream_temp_decl()
      USE PRMS_STRMTEMP
      IMPLICIT NONE
! Functions
      INTRINSIC :: INDEX
      INTEGER, EXTERNAL :: declparam, declvar, getdim, control_integer
      EXTERNAL :: read_error, print_module
!***********************************************************************
      stream_temp_decl = 0

      CALL print_module(MODDESC, MODNAME, Version_stream_temp)

      ! 0 = compute shade; 1 = specified constant
      IF ( control_integer(Stream_temp_shade_flag, 'stream_temp_shade_flag')/=0 ) Stream_temp_shade_flag = 0

! Declared Variables
      ALLOCATE ( Seg_width(Nsegment) )
      IF ( declvar( MODNAME, 'seg_width', 'nsegment', Nsegment, 'real', &
     &     'Width of each segment', &
     &     'meters', Seg_width)/=0 )  CALL read_error(3, 'seg_width')

      ALLOCATE (Seg_tave_water(Nsegment) ) ! previous ??
      IF ( declvar( MODNAME, 'seg_tave_water', 'nsegment', Nsegment, 'real', &
     &     'Computed daily mean stream temperature for each segment', &
     &     'degrees Celsius', Seg_tave_water)/=0 ) CALL read_error(3, 'seg_tave_water')

      ALLOCATE ( seg_tave_upstream(Nsegment) )
      IF ( declvar( MODNAME, 'seg_tave_upstream', 'nsegment', Nsegment, 'real', &
     &     'Temperature of streamflow entering each segment', &
     &     'degrees Celsius', seg_tave_upstream)/=0 )   CALL read_error(3,'seg_tave_upstream')

      ALLOCATE ( Seg_humid(Nsegment) )
      IF ( declvar( MODNAME, 'seg_humid', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average relative humidity for each segment from HRUs contributing flow to the segment', &
     &     'decimal fraction', Seg_humid)/=0 ) CALL read_error(3,'seg_humid')

      ALLOCATE ( Seg_melt(Nsegment) )
      IF ( declvar( MODNAME, 'seg_melt', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average snowmelt for each segment from HRUs contributing flow to the segment', &
     &     'inches', Seg_melt)/=0 ) CALL read_error(3, 'seg_melt')

      ALLOCATE ( Seg_rain(Nsegment) )
      IF ( declvar( MODNAME, 'seg_rain', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average rainfall for each segment from HRUs contributing flow to the segment', &
     &     'inches', Seg_rain)/=0 ) CALL read_error(3, 'seg_rain')

      ALLOCATE ( Seg_tave_air(Nsegment) )
      IF ( declvar( MODNAME, 'seg_tave_air', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average air temperature for each segment from HRUs contributing flow to the segment', &
     &     'degrees Celsius', Seg_tave_air)/=0 ) CALL read_error(3, 'seg_tave_air')

      ALLOCATE ( Seg_potet(Nsegment) )
      IF ( declvar( MODNAME, 'seg_potet', 'nsegment', Nsegment, 'double', &
     &     'HRU area-weighted average potential ET for each segment', &
     &     'inches', Seg_potet)/=0 ) CALL read_error(3, 'seg_potet')

      ALLOCATE ( Seg_ccov(Nsegment) )
      IF ( declvar( MODNAME, 'seg_ccov', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average cloud cover fraction for each segment from HRUs contributing flow to the segment', &
     &     'decimal fraction', Seg_ccov )/=0 ) CALL read_error(3, 'seg_ccov')

      ALLOCATE(Seg_shade(Nsegment))
      IF (declvar(MODNAME, 'seg_shade', 'nsegment', Nsegment, 'real', &
     &     'Area-weighted average shade fraction for each segment', &
     &     'decimal fraction', seg_shade)/=0 ) CALL read_error(3, 'seg_shade')

      ALLOCATE ( Seg_daylight(Nsegment) )
      IF ( declvar( MODNAME, 'seg_daylight', 'nsegment', Nsegment, 'real', &
     &     'Hours of daylight', &
     &     'hours', Seg_daylight)/=0 )   CALL read_error(3,'seg_daylight')

      ALLOCATE(seg_tave_gw(Nsegment))
      IF ( declvar( MODNAME, 'seg_tave_gw', 'nsegment', Nsegment, 'real', &
     &     'groundwater temperature', &
     &     'degrees Celsius', seg_tave_gw)/=0 )   CALL read_error(3,'seg_tave_gw')

      ALLOCATE(seg_tave_ss(Nsegment))
      IF ( declvar( MODNAME, 'seg_tave_ss', 'nsegment', Nsegment, 'real', &
     &     'subsurface temperature', &
     &     'degrees Celsius', seg_tave_ss)/=0 )   CALL read_error(3,'seg_tave_ss')

      ALLOCATE(seg_tave_sroff(Nsegment))
      IF ( declvar( MODNAME, 'seg_tave_sroff', 'nsegment', Nsegment, 'real', &
     &     'surface runoff temperature', &
     &     'degrees Celsius', seg_tave_sroff)/=0 )   CALL read_error(3,'seg_tave_sroff')

      ALLOCATE(seg_tave_lat(Nsegment))
      IF ( declvar( MODNAME, 'seg_tave_lat', 'nsegment', Nsegment, 'real', &
     &     'lateral flow temperature', &
     &     'degrees Celsius', seg_tave_lat)/=0 )   CALL read_error(3,'seg_tave_lat')

      ALLOCATE (Press(Nsegment) )
      ALLOCATE ( Seg_hru_count(Nsegment) )
      ALLOCATE (Seg_carea_inv(Nsegment) )
      ALLOCATE ( Seg_close(Nsegment) )
      ALLOCATE (gw_sum(Nsegment), ss_sum(Nsegment))
      ALLOCATE (gw_silo(nsegment,DAYS_PER_YEAR), ss_silo(nsegment,DAYS_PER_YEAR))
      ALLOCATE (hru_area_sum(nsegment))

      IF ( declparam( MODNAME, 'albedo', 'one', 'real', &
     &     '0.10', '0.0', '1.0', &
     &     'Short-wave solar radiation reflected by streams', &
     &     'Short-wave solar radiation reflected by streams', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'albedo')

      ALLOCATE(lat_temp_adj(Nsegment,MONTHS_PER_YEAR))
      IF ( declparam( MODNAME, 'lat_temp_adj', 'nsegment,nmonths', 'real', &
     &     '0.0', '-5.0', '5.0', &
     &     'Correction factor to adjust the bias of the temperature of the lateral inflow', &
     &     'Correction factor to adjust the bias of the temperature of the lateral inflow', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'lat_temp_adj')

      ALLOCATE ( Seg_length(Nsegment) )
      IF ( declparam( MODNAME, 'seg_length', 'nsegment', 'real', &
     &     '1000.0', '1.0', '100000.0', &
     &     'Length of each segment', &
     &     'Length of each segment', &
     &     'meters')/=0 ) CALL read_error(1, 'seg_length')

      ALLOCATE ( Seg_slope(Nsegment) )
      IF ( declparam( MODNAME, 'seg_slope', 'nsegment', 'real', &
     &     '0.015', '0.0001', '2.0', &
     &     'Bed slope of each segment', &
     &     'Bed slope of each segment', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'seg_slope')

      ALLOCATE ( width_alpha(Nsegment) )
      IF ( declparam( MODNAME, 'width_alpha', 'nsegment', 'real', &
     &     '0.015', '0.0001', '2.0', &
     &     'Alpha coefficient in power function for width calculation', &
     &     'Alpha coefficient in power function for width calculation', &
     &     'unknown')/=0 ) CALL read_error(1, 'width_alpha')

      ALLOCATE ( width_m(Nsegment) )
      IF ( declparam( MODNAME, 'width_m', 'nsegment', 'real', &
     &     '0.015', '0.0001', '2.0', &
     &     'M value in power function for width calculation', &
     &     'M value in power function for width calculation', &
     &     'unknown')/=0 ) CALL read_error(1, 'width_m')

      IF ( Stream_temp_shade_flag==OFF .OR. Model==DOCUMENTATION ) THEN
         ALLOCATE ( Azrh(Nsegment) )
         IF ( declparam( MODNAME, 'azrh', 'nsegment', 'real', &
     &       '0.0', '-1.5708', '1.5708', &
     &       'Azimuth angle of each segment', &
     &       'Azimuth angle of each segment', &
     &       'radians')/=0 ) CALL read_error(1, 'azrh')

         ALLOCATE ( Alte(Nsegment) )
         IF ( declparam( MODNAME, 'alte', 'nsegment', 'real', &
     &       '0.0', '0.0','1.57079633', &
     &       'East bank topographic altitude', &
     &       'East bank topographic altitude of each segment', &
     &       'radians')/=0 ) CALL read_error(1, 'alte')

         ALLOCATE ( Altw(Nsegment) )
         IF ( declparam( MODNAME, 'altw', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.57079633', &
     &       'West bank topographic altitude', &
     &       'West bank topographic altitude of each segment', &
     &       'radians')/=0 ) CALL read_error(1, 'altw')

         ALLOCATE ( Vce(Nsegment) )
         IF ( declparam( MODNAME, 'vce', 'nsegment', 'real', &
     &       '0.0', '0.0', '15.0', &
     &       'East bank average vegetation crown width', &
     &       'East bank average vegetation crown width for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'vce')

         ALLOCATE ( Vdemx(Nsegment) )
         IF ( declparam( MODNAME, 'vdemx', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Maximum east bank vegetation density', &
     &       'Maximum east bank vegetation density for each segment', &
     &       'decimal fraction')/=0 )  CALL read_error(1, 'vdemx')

         ALLOCATE ( Vdemn(Nsegment) )
         IF ( declparam( MODNAME, 'vdemn', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Minimum east bank vegetation density', &
     &       'Minimum east bank vegetation density for each segment', &
     &       'decimal fraction')/=0 )  CALL read_error(1, 'vdemn')

         ALLOCATE ( Vhe(Nsegment) )
         IF ( declparam( MODNAME, 'vhe', 'nsegment', 'real', &
     &       '0.0', '0.0', '30.0', &
     &       'East bank vegetation height', &
     &       'East bank average vegetation height for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'vhe')

         ALLOCATE ( Voe(Nsegment) )
         IF ( declparam( MODNAME, 'voe', 'nsegment', 'real', &
     &       '0.0', '0.0', '100.0',&
     &       'East bank vegetation offset', &
     &       'East bank vegetation offset for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'voe')

         ALLOCATE ( Vcw(Nsegment) )
         IF ( declparam( MODNAME, 'vcw', 'nsegment', 'real', &
     &       '0.0', '0.0', '15.0', &
     &       'West bank vegetation crown width', &
     &       'West bank average vegetation crown width for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'vcw')

         ALLOCATE ( Vdwmx(Nsegment) )
         IF ( declparam( MODNAME, 'vdwmx', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Maximum west bank vegetation density', &
     &       'Maximum west bank vegetation density for each segment', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'vdwmx')

         ALLOCATE ( Vdwmn(Nsegment) )
         IF ( declparam( MODNAME, 'vdwmn', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Minimum west bank vegetation density', &
     &       'Minimum west bank vegetation density for each segment', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'vdwmn')

         ALLOCATE ( Vhw(Nsegment) )
         IF ( declparam( MODNAME, 'vhw', 'nsegment', 'real', &
     &       '0.0', '0.0', '30.0', &
     &       'West bank vegetation height', &
     &       'West bank average vegetation height for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'vhw')

         ALLOCATE ( Vow(Nsegment) )
         IF ( declparam( MODNAME, 'vow', 'nsegment', 'real', &
     &       '0.0', '0.0', '100.0', &
     &       'West bank vegetation offset', &
     &       'West bank vegetation offset for each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'vow')
      ENDIF

      IF ( Stream_temp_shade_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
         ALLOCATE ( Segshade_sum(Nsegment) )
         IF ( declparam( MODNAME, 'segshade_sum', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0.', &
     &       'Total shade fraction for summer vegetation', &
     &       'Total shade fraction for summer vegetation; required when stream_temp_flag=1', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'segshade_sum')

         ALLOCATE ( Segshade_win(Nsegment) )
         IF ( declparam( MODNAME, 'segshade_win', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0.', &
     &       'Total shade fraction for winter vegetation', &
     &       'Total shade fraction for winter vegetation; required when stream_temp_flag=1', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'segshade_win')
      ENDIF

      ALLOCATE (ss_tau(Nsegment) )
      IF ( declparam( MODNAME, 'ss_tau', 'nsegment', 'integer', &
     &     '30', '1', '365', &
     &     'Average residence time of subsurface interflow', &
     &     'Average residence time of subsurface interflow', &
     &     'days')/=0 ) CALL read_error(1, 'ss_tau')

      ALLOCATE (gw_tau(Nsegment) )
      IF ( declparam( MODNAME, 'gw_tau', 'nsegment', 'integer', &
     &     '365', '1', '365', &
     &     'Average residence time in groundwater flow', &
     &     'Average residence time in groundwater flow', &
     &     'days')/=0 ) CALL read_error(1, 'gw_tau')

      IF ( declparam( MODNAME, 'melt_temp', 'one', 'real', &
     &     '1.5', '0.0', '10.0', &
     &     'Temperature at which snowmelt enters a stream', &
     &     'Temperature at which snowmelt enters a stream', &
     &     'degrees Celsius')/=0 ) CALL read_error(1, 'melt_temp')

      IF ( declparam( MODNAME, 'maxiter_sntemp', 'one', 'integer', &
     &     '1000', '10', '2000', &
     &     'Maximum number of Newton-Raphson iterations to compute stream temperature', &
     &     'Maximum number of Newton-Raphson iterations to compute stream temperature', &
     &     'none')/=0 ) CALL read_error(1, 'maxiter_sntemp')

      IF ( Strmtemp_humidity_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN  ! specified constant
         ALLOCATE ( Seg_humidity(Nsegment, MONTHS_PER_YEAR) )
         IF ( declparam( MODNAME, 'seg_humidity', 'nsegment,nmonths', 'real', &
     &       '0.7', '0.0', '1.0', &
     &       'Mean monthly humidity for each segment', &
     &       'Mean monthly humidity for each segment, used when values not input in CBH File', &
     &       'decimal fraction')/=0 )  CALL read_error(1, 'seg_humidity')
      ELSEIF ( Strmtemp_humidity_flag==2 .OR. Model==DOCUMENTATION ) THEN  ! use station data
         ALLOCATE ( Seg_humidity_sta(Nsegment) )
         IF ( declparam(MODNAME, 'seg_humidity_sta', 'nsegment', 'integer', &
     &       '0', 'bounded', 'nhumid', &
     &       'Index of humidity measurement station for each stream segment', &
     &       'Index of humidity measurement station for each stream segment', &
     &       'none')/=0 ) CALL read_error(1, 'seg_humidity_sta')
      ENDIF

      ALLOCATE (seg_lat(nsegment))
      IF ( declparam( MODNAME, 'seg_lat', 'nsegment', 'real', &
     &     '40.0', '-90.0', '90.0', &
     &     'Segment latitude', &
     &     'Latitude of each segment', &
     &     'degrees North')/=0 ) CALL read_error(1, 'seg_lat')

      ALLOCATE (seg_elev(nsegment))
      IF (declparam(MODNAME, 'seg_elev', 'nsegment', 'real', &
     &     '0.0', '-1000.0', '30000.0', &
     &     'Segment elevation at midpoint', 'Segment elevation at midpoint', &
     &     'meters')/=0 ) CALL read_error(1, 'seg_elev')

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==8 ) THEN
        ALLOCATE ( Stream_tave_init(Nsegment) )
        IF ( declparam(MODNAME, 'stream_tave_init', 'nsegment', 'real', &
     &       '0.0', '-10.0', '100.0', &
     &       'Initial average stream temperature in each segment', &
     &       'Initial average stream temperature in each segment at the beginning of a simulation', &
     &       'degrees Celsius')/=0 ) CALL read_error(1, 'stream_tave_init')
      ENDIF

      END FUNCTION stream_temp_decl

!***********************************************************************
!    stream_temp_init - Initialize module - get parameter values
!***********************************************************************
      INTEGER FUNCTION stream_temp_init()
      USE PRMS_STRMTEMP
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
      USE PRMS_OBS, ONLY: Nhumid
      USE PRMS_ROUTING, ONLY: Hru_segment, Tosegment, Segment_order, Segment_up
      IMPLICIT NONE
! Functions
      INTRINSIC :: COS, SIN, ABS, SIGN, ASIN, maxval
      INTEGER, EXTERNAL :: getparam
      REAL, EXTERNAL :: solalt
      EXTERNAL :: read_error, checkdim_param_limits, error_stop
! Local Variables
      INTEGER :: i, j, k, iseg, ierr, ii, this_seg
      REAL :: tan_d, tano, sinhro, temp, decl, cos_d, tanod, alrs
!***********************************************************************
      stream_temp_init = 0

      IF ( getparam( MODNAME, 'albedo', 1, 'real', Albedo)/=0 ) CALL read_error(2, 'albedo')
      IF ( getparam( MODNAME, 'lat_temp_adj', Nsegment*MONTHS_PER_YEAR, 'real', lat_temp_adj)/=0 ) &
     &     CALL read_error(2, 'lat_temp_adj')
      IF ( getparam( MODNAME, 'seg_length', Nsegment, 'real', Seg_length)/=0 ) CALL read_error(2, 'seg_length')

      IF (getparam(MODNAME, 'seg_lat', Nsegment, 'real', Seg_lat)/=0 ) CALL read_error(2, 'seg_lat')
!     Convert latitude from degrees to radians
      seg_lat = seg_lat * DEG_TO_RAD

      IF (getparam(MODNAME, 'seg_elev', Nsegment, 'real', Seg_elev)/=0 ) CALL read_error(2, 'seg_elev')

! convert stream length in meters to km
      Seg_length = Seg_length / 1000.0

      IF ( getparam( MODNAME, 'seg_slope', Nsegment, 'real', Seg_slope)/=0 ) CALL read_error(2, 'seg_slope')
      IF ( getparam( MODNAME, 'width_alpha', Nsegment, 'real', width_alpha)/=0 ) CALL read_error(2, 'width_alpha')
      IF ( getparam( MODNAME, 'width_m', Nsegment, 'real', width_m)/=0 ) CALL read_error(2, 'width_m')

      IF ( Stream_temp_shade_flag==OFF ) THEN
         IF ( getparam( MODNAME, 'azrh', Nsegment, 'real', Azrh)/=0 ) CALL read_error(2, 'azrh')
         IF ( getparam( MODNAME, 'alte', Nsegment, 'real', Alte)/=0 ) CALL read_error(2, 'alte')
         IF ( getparam( MODNAME, 'altw', Nsegment, 'real', Altw)/=0 ) CALL read_error(2, 'altw')
         IF ( getparam( MODNAME, 'vce', Nsegment, 'real', Vce)/=0 ) CALL read_error(2, 'vce')
         IF ( getparam( MODNAME, 'vdemx', Nsegment, 'real', Vdemx)/=0 ) CALL read_error(2, 'vdemx')
         IF ( getparam( MODNAME, 'vdemn', Nsegment, 'real', Vdemn)/=0 ) CALL read_error(2, 'vdemn')
         IF ( getparam( MODNAME, 'vhe', Nsegment, 'real', Vhe)/=0 ) CALL read_error(2, 'vhe')
         IF ( getparam( MODNAME, 'voe', Nsegment, 'real', Voe)/=0 ) CALL read_error(2, 'voe')
         IF ( getparam( MODNAME, 'vcw', Nsegment, 'real', Vcw)/=0 ) CALL read_error(2, 'vcw')
         IF ( getparam( MODNAME, 'vdwmx', Nsegment, 'real', Vdwmx)/=0 ) CALL read_error(2, 'vdwmx')
         IF ( getparam( MODNAME, 'vdwmn', Nsegment, 'real', Vdwmn)/=0 ) CALL read_error(2, 'vdwmn')
         IF ( getparam( MODNAME, 'vhw', Nsegment, 'real', Vhw)/=0 ) CALL read_error(2, 'vhw')
         IF ( getparam( MODNAME, 'vow', Nsegment, 'real', Vow)/=0 ) CALL read_error(2, 'vow')
      ELSE
         IF ( getparam( MODNAME, 'segshade_sum', Nsegment, 'real', Segshade_sum)/=0 ) CALL read_error(2, 'segshade_sum')
         IF ( getparam( MODNAME, 'segshade_win', Nsegment, 'real', Segshade_win)/=0 ) CALL read_error(2, 'segshade_win')
      ENDIF

      IF ( getparam( MODNAME, 'ss_tau', Nsegment, 'integer', Ss_tau)/=0 ) CALL read_error(2, 'ss_tau')
      IF ( getparam( MODNAME, 'gw_tau', Nsegment, 'integer', Gw_tau)/=0 ) CALL read_error(2, 'Gw_tau')
      IF ( getparam( MODNAME, 'melt_temp', 1, 'real', Melt_temp)/=0 ) CALL read_error(2, 'melt_temp')
      IF ( getparam( MODNAME, 'maxiter_sntemp', 1, 'real', Maxiter_sntemp)/=0 ) CALL read_error(2, 'maxiter_sntemp')

      ierr = 0
      IF ( Strmtemp_humidity_flag==1 ) THEN
         IF ( getparam( MODNAME, 'seg_humidity', Nsegment*MONTHS_PER_YEAR, 'real', Seg_humidity)/=0 ) &
     &      CALL read_error(2, 'seg_humidity')
      ELSEIF ( Strmtemp_humidity_flag==2 ) THEN ! use station data
         IF ( getparam(MODNAME, 'seg_humidity_sta', Nsegment, 'integer', Seg_humidity_sta)/=0 ) &
     &      CALL read_error(2, 'seg_humidity_sta')
         DO i = 1, Nsegment
            CALL checkdim_param_limits(i, 'seg_humidity_sta', 'nhumid', Seg_humidity_sta(i), 1, Nhumid, ierr)
         ENDDO
      ENDIF

! Initialize declared variables
      seg_tave_upstream = 0.0
      Seg_potet = 0.0D0
      Seg_humid = 0.0
      Seg_width = 0.0
      Seg_ccov = 0.0
      Seg_tave_air = 0.0
      seg_tave_gw = 0.0
      seg_tave_ss = 0.0
      seg_tave_sroff = 0.0

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==8 ) THEN
        IF ( getparam(MODNAME, 'stream_tave_init', Nsegment, 'real', Stream_tave_init)/=0 ) CALL read_error(2, 'stream_tave_init')
        Seg_tave_water = Stream_tave_init
      ENDIF
      IF ( Init_vars_from_file == 0 ) THEN
         gw_silo =  0.0
         ss_silo =  0.0
         gw_sum = 0.0
         ss_sum = 0.0
! these are set to zero because they will be incremented to 1 down in the run function
         gw_index = 0
         ss_index = 0
      ENDIF

      Seg_daylight = 12.0
      IF ( Stream_temp_shade_flag==OFF ) THEN
         ALLOCATE ( Cos_seg_lat(Nsegment), Sin_seg_lat(Nsegment), Horizontal_hour_angle(MAX_DAYS_PER_YEAR,Nsegment) )
         ALLOCATE ( Total_shade(MAX_DAYS_PER_YEAR,Nsegment), Sin_declination(MAX_DAYS_PER_YEAR,Nsegment) )
         ALLOCATE ( Sin_lat_decl(MAX_DAYS_PER_YEAR,Nsegment), Cos_lat_decl(MAX_DAYS_PER_YEAR,Nsegment) )
         ALLOCATE ( Max_solar_altitude(MAX_DAYS_PER_YEAR,Nsegment), Level_sunset_azimuth(MAX_DAYS_PER_YEAR,Nsegment) )
         ALLOCATE ( Local_sunset_hour_angle(MAX_DAYS_PER_YEAR,Nsegment), Local_sunrise_hour_angle(MAX_DAYS_PER_YEAR,Nsegment) )
         ALLOCATE ( Shade_jday(Nsegment, MAX_DAYS_PER_YEAR), Svi_jday(Nsegment, MAX_DAYS_PER_YEAR) )
         ALLOCATE ( Sin_alrs(MAX_DAYS_PER_YEAR,Nsegment) )
         Shade_jday = 0.0
         Svi_jday = 0.0
         Seg_lat = 0.0
      ENDIF

! Figure out how many HRUs are connected to each segment
      Seg_hru_count = 0
      DO k = 1, Active_hrus
         j = Hru_route_order(k)
         i = Hru_segment(j)
         IF ( i==0 ) CYCLE
         Seg_hru_count(i) = Seg_hru_count(i) + 1
      ENDDO

! find segments that are too short and print them out as they are found
      DO i = 1, Nsegment
         IF ( Seg_length(i)<NEARZERO ) THEN
            PRINT *, 'ERROR, seg_length too small for segment:', i, ', value:', Seg_length(i)
            ierr = 1
         ENDIF
      ENDDO

! exit if there are any segments that are too short
      IF ( ierr==1 ) THEN
         Inputerror_flag = ierr
         RETURN
      ENDIF

      Seg_close = Segment_up ! assign upstream values
      DO j = 1, Nsegment ! set values based on routing order for segments without associated HRUs
         i = Segment_order(j)

      ! If a segment does not have any HRUs, need to find the closest one for elevation and latitude info
      ! NOTE: seg_close variable can go upstream, downstream, or offstream looking for the "closest" segment with
      ! an HRU. This is not approprite to use in a situation where computed values are going to be taken from
      ! the closest HRU (i.e. flow).
      !
      ! This does work for NHM network (most comprehensive test).
      !
         IF ( Seg_hru_count(i)==0 ) THEN
            IF ( Segment_up(i)==0 ) THEN
               IF ( Tosegment(i)>0 ) THEN ! assign downstream values
                  Seg_close(i) = Tosegment(i) ! don't have a value yet, need to fix
               ELSE ! no upstream or downstream segment
                  IF ( j>1 ) THEN
                     Seg_close(i) = Segment_order(j-1) ! set to previous segment id
                  ELSE
                     Seg_close(i) = Segment_order(j+1) ! assume at least 2 segments
                  ENDIF
               ENDIF
            ENDIF
            IF ( Seg_elev(Seg_close(i))==30000.0 ) THEN ! need different segment
               iseg = -1
               DO k = j+1, Nsegment ! find first segment with valid values
                  ii = Segment_order(k)
                  IF ( Seg_hru_count(ii)>0 ) THEN
                     Seg_close(i) = ii
                     EXIT
                  ENDIF
               ENDDO
               IF ( iseg==-1 ) THEN
                  IF ( j>1 ) THEN
                     Seg_close(i) = Segment_order(j-1) ! set to previous segment id
                  ELSE ! this is a problem, shouldn't happen
                     CALL error_stop('segments do not have associated HRUs', ERROR_param)
                    ! Seg_close(i) = Segment_order(1) ! set to first segment id
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! Compute atmospheric pressure based on segment elevation.
         Press(i) = 1013.0 - (0.1055 * Seg_elev(i))

         IF ( Stream_temp_shade_flag==0 ) THEN
!  LATITUDE TRIGONOMETRIC PARAMETERS
            Cos_seg_lat(i) = COS(Seg_lat(i)) ! coso
            IF ( Cos_seg_lat(i) < NEARZERO ) Cos_Seg_lat(i) = NEARZERO
            Sin_seg_lat(i) = SIN(Seg_lat(i)) ! sino
            tano = Sin_seg_lat(i) / Cos_seg_lat(i)
            DO k = 1, MAX_DAYS_PER_YEAR
!  DECLINATION TRIGONOMETRIC PARAMETERS
               decl = 0.40928 * COS(((2.0 * PI) / DAYS_YR) * (172.0 - k))
               cos_d = COS(decl)
               Sin_declination(k, i) = SIN(decl) ! sin_d
               IF ( cos_d < NEARZERO ) cos_d = NEARZERO
               tan_d = Sin_declination(k, i) / cos_d
!
!  JOINT LATITUDE & DECLINATION TRIGONOMETRIC PARAMETERS
               Cos_lat_decl(k, i) = Cos_seg_lat(i) * cos_d ! cosod
               Sin_lat_decl(k, i) = Sin_seg_lat(i) * Sin_declination(k, i) ! sinod
               tanod = tano * tan_d
               IF ( ABS(tanod) > 1.0 ) tanod = SIGN(1.0,tanod)

!  LEVEL-PLAIN SUNRISE/SET HOUR ANGLE
               Horizontal_hour_angle(k, i) = ACOS(-tanod) ! hrso
               sinhro = SIN(Horizontal_hour_angle(k, i))
!
!  LEVEL-PLAIN SOLAR AZIMUTH
               temp = -Sin_declination(k, i)/Cos_Seg_lat(i)
               IF ( ABS(temp) > 1.0 ) temp = SIGN(1.0,temp)
               Level_sunset_azimuth(k, i) = ACOS(temp) ! azso
!
!  MAXIMUM POSSIBLE SOLAR ALTITUDE
               Max_solar_altitude(k, i) = ASIN( Sin_lat_decl(k,i) + Cos_lat_decl(k,i) ) ! alsmx
!
!  TOTAL POTENTIAL SHADE ON LEVEL-PLAIN ! totsh
               Total_shade(k, i) = 2.0 * ((Horizontal_hour_angle(k, i) * Sin_lat_decl(k, i)) + (sinhro * Cos_lat_decl(k, i)))
               IF ( Total_shade(k, i) < NEARZERO ) Total_shade(k, i) = NEARZERO
!
!  CHECK FOR REACH AZIMUTH LESS THAN SUNRISE
               IF ( Azrh(i) <= (-Level_sunset_azimuth(k, i)) ) THEN
                  alrs = 0.0
!
!  CHECK FOR REACH AZIMUTH GREATER THAN SUNSET
               ELSEIF ( Azrh(i) >= Level_sunset_azimuth(k, i) ) THEN
                  alrs = 0.0
!
!  REACH AZIMUTH IS BETWEEN SUNRISE & SUNSET
               ELSEIF ( Azrh(i) == 0.0 ) THEN
                  alrs = Max_Solar_altitude(k, i)
               ELSE
                  alrs = solalt(Cos_seg_lat(i), Sin_seg_lat(i), Sin_declination(k,i), Azrh(i), 0.0, Max_Solar_altitude(k,i))
                  Sin_alrs(k, i) = SIN(alrs)
!
!  END REACH & SOLAR AZIMUTH CHECK
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!     There may be headwater segments that do not have any HRUs and do not have any upstream segments to produce
!     streamflow. These segments will never have any streamflow, and consequently never be able to simulate
!     stream temperature. This block finds these and sets the stream temperature value to -99.9. Subsequent code
!     should be able to check if the temperature value is less than -99.0 and know that it doesn't need to do
!     any stream temperature calculation because there will never be any water in the segment.
!
!     This code is similar to the code above that computes latitude and elevation, but is different because it
!     must always look upstream because the downstream computations will not have been done when the current
!     segment is being calculated.
      do j = 1, nsegment
         this_seg = segment_order(j)

!         Check if this segment has any HRUs, keep moving up stream if not.
         do
            if (seg_hru_count(this_seg) .eq. 0) then
               ! Hit the headwater segment without finding any HRUs (i.e. sources of streamflow)
               ! Set the stream temp to -99.9 for this segment because there will never be any flow in this segment
               if (segment_up(this_seg) .eq. 0) then
                  Seg_tave_water(segment_order(j)) = -99.9
                  exit
               endif

               ! There is an upstream segment, check that segment for HRUs
               this_seg = segment_up(this_seg)
            else
               ! This segment has HRUs so there will be no streamflow
               exit
            endif
         enddo
      enddo

! For each segment, figure out how many upstream segments.
      ALLOCATE(upstream_count(Nsegment))
      upstream_count = 0
      do i = 1, nsegment
         do j = 1, nsegment
            if (tosegment(j) .eq. i) then
               upstream_count(i) = upstream_count(i) + 1
            endif
         end do
      end do

! For each segment, figure out the upstream segments. These will be looped over to determine inflows and temps to each segment
      ALLOCATE(upstream_idx(Nsegment, maxval(upstream_count)))
      upstream_idx = 0
      upstream_count = 0
      do i = 1, nsegment
         do j = 1, nsegment
            if (tosegment(j) .eq. i) then
               upstream_count(i) = upstream_count(i) + 1
               upstream_idx(i,upstream_count(i)) = j
            endif
         end do
      end do

!      do i = 1, nsegment
!         write(*, fmt="(1x,a,i0)", advance="no") "segment #", i
!         write(*, fmt="(1x,a,i0)", advance="no") " ", upstream_count(i)
!         do j = 1, upstream_count(i)
!            write(*, fmt="(1x,a,i0)", advance="no") " ", upstream_idx(i,j)
!         end do
!         write(*, fmt="(1x,a)",advance="yes") " done"
!      end do

      END FUNCTION stream_temp_init

!***********************************************************************
!     stream_temp_run - Computes stream temperatures
!***********************************************************************
      INTEGER FUNCTION stream_temp_run()
      USE PRMS_STRMTEMP
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area
      USE PRMS_SET_TIME, ONLY: Summer_flag, Nowmonth
      USE PRMS_CLIMATEVARS, ONLY: Tavgc, Potet, Hru_rain, Swrad
      USE PRMS_CLIMATE_HRU, ONLY: Humidity_hru
      USE PRMS_FLOWVARS, ONLY: Seg_outflow
      USE PRMS_SNOW, ONLY: Snowmelt
      USE PRMS_ROUTING, ONLY: Hru_segment, Segment_order, Seginc_swrad
      USE PRMS_OBS, ONLY: Humidity
      USE PRMS_SET_TIME, ONLY: Nowmonth, Jday !, Nowyear, Nowday
      USE PRMS_SOLTAB, ONLY: Soltab_potsw, Hru_cossl

      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
      REAL, EXTERNAL :: twavg, twmax, get_segwidth
      EXTERNAL :: equilb, lat_inflow, shday
! Local Variables
      REAL :: harea, svi, fs
      INTEGER :: i, j, k, kk, iseg
      REAL :: te, ak1, ak2, ccov
      DOUBLE PRECISION :: qlat
      REAL :: t_o, up_temp
!***********************************************************************
      stream_temp_run = 0
      Seg_tave_air = 0.0

! Humidity info come from parameter file when Strmtemp_humidity_flag==1
! Otherwise it comes as daily values per HRU from CBH. Code for this is
! down in the HRU loop.
      IF ( Strmtemp_humidity_flag==1 ) THEN
         DO i = 1, Nsegment
            Seg_humid(i) = Seg_humidity(i, Nowmonth)
         ENDDO
      ELSEIF ( Strmtemp_humidity_flag==2 ) THEN ! use station data
         DO i = 1, Nsegment
            Seg_humid(i) = Humidity(Seg_humidity_sta(i))
         ENDDO
      ELSE
         Seg_humid = 0.0
      ENDIF

      Seg_potet = 0.0D0
      Seg_ccov = 0.0
      Seg_melt = 0.0
      Seg_rain = 0.0
      hru_area_sum = 0.0

      ! Compute segment lateral inflow temperatures and segment meteorological values

      DO k = 1, Active_hrus
         j = Hru_route_order(k)

! DANGER HACK
! On restart, sometimes soltab_potsw comes in as zero. It should never be zero as 
! this results in divide by 0.0
         if (Soltab_potsw(jday, j) <= 10.0) then
            ccov = 1.0 - (Swrad(j) / sngl(10.0) * sngl(Hru_cossl(j)))
         else
            ccov = 1.0 - (Swrad(j) / sngl(Soltab_potsw(jday, j)) * sngl(Hru_cossl(j)))
         endif

         if (ccov .ne. ccov) then
             write(*,*) "ccov is nan", j, Swrad(j), Soltab_potsw(jday, j), Hru_cossl(j)
         endif

         IF ( ccov<NEARZERO ) THEN
            ccov = 0.0
         ELSEIF ( ccov>1.0 ) THEN
            ccov = 1.0
         ENDIF

         harea = Hru_area(j)
         i = Hru_segment(j)
         IF ( i==0 ) CYCLE

! Compute temperature of surface runoff here for HRU and stream segments
         Seg_tave_air(i) = Seg_tave_air(i) + Tavgc(j)*harea
         hru_area_sum(i) = hru_area_sum(i) + harea

! Compute segment humidity if info is specified in CBH as timeseries by HRU
         IF ( Strmtemp_humidity_flag==0 ) then
            Seg_humid(i) = Seg_humid(i) + Humidity_hru(j)*harea
         endif

! Figure out the contributions of the HRUs to each segment for these drivers.
         Seg_ccov(i) = Seg_ccov(i) + ccov*harea

         Seg_potet(i) = Seg_potet(i) + DBLE( Potet(j)*harea )
         Seg_melt(i) = Seg_melt(i) + Snowmelt(j)*harea
         Seg_rain(i) = Seg_rain(i) + Hru_rain(j)*harea
      ENDDO


      DO j = 1, Nsegment
         i = Segment_order(j)
         IF ( Seg_hru_count(i)>0 ) THEN
!            carea = Seg_carea_inv(i)
            Seg_ccov(i) = Seg_ccov(i) / hru_area_sum(i)
            Seg_potet(i) = Seg_potet(i) / dble(hru_area_sum(i))
            Seg_tave_air(i) = Seg_tave_air(i) / hru_area_sum(i)
            Seg_melt(i) = Seg_melt(i) / hru_area_sum(i)
            Seg_rain(i) = Seg_rain(i) / hru_area_sum(i)
            IF ( Strmtemp_humidity_flag==0 ) then
               Seg_humid(i) = Seg_humid(i) / hru_area_sum(i)

! DANGER potential hack here: Should CBH humidity data be converted to decimal fraction in
! the CBH file? Probably so. For now, convert it here.
! Humidity coming from CBH is in percent, not decimal fraction
               Seg_humid(i) = Seg_humid(i) * 0.01
            endif
         ELSE
! This block for segments that don't have contributing HRUs
            iseg = Seg_close(i) ! doesn't work if upstream segment
            Seg_tave_air(i) = Seg_tave_air(iseg)
            Seg_ccov(i) = Seg_ccov(iseg)
            Seg_potet(i) = Seg_potet(iseg)
            Seg_melt(i) = Seg_melt(iseg)
            Seg_rain(i) = Seg_rain(iseg)
            IF ( Strmtemp_humidity_flag==0 ) then
               Seg_humid(i) = Seg_humid(iseg)*Seg_carea_inv(iseg) ! ??
! DANGER Humidity coming from CBH is in percent, not decimal fraction
! Same as comment in above block
               Seg_humid(i) = Seg_humid(i) * 0.01
            endif
         ENDIF

!         if (Seg_tave_air(i) .ne. Seg_tave_air(i)) then
!            write(*,*) "Seg_tave_air is nan", i
!         endif
!
!         if (Seg_ccov(i) .ne. Seg_ccov(i)) then
!            write(*,*) "Seg_ccov is nan", i
!         endif
!
!         if (Seg_potet(i) .ne. Seg_potet(i)) then
!            write(*,*) "Seg_potet is nan", i
!         endif
!
!         if (Seg_melt(i) .ne. Seg_melt(i)) then
!            write(*,*) "Seg_melt is nan", i
!         endif
!             
!         if (Seg_rain(i) .ne. Seg_rain(i)) then
!            write(*,*) "Seg_rain is nan", i
!         endif
!             
!         if (Seg_humid(i) .ne. Seg_humid(i)) then
!            write(*,*) "Seg_humid is nan", i
!         endif
 
      ENDDO


! Compute the running averages for groundwater and subsurface temperatures.
      if (gw_index >= gw_tau(i)) then
         gw_index = 1
      else
         gw_index = gw_index + 1
      endif

      if (ss_index >= ss_tau(i)) then
         ss_index = 1
      else
         ss_index = ss_index + 1
      endif

      ! Mark all of the upstream segment temperatures as not having been computed yet.
      ! If the value is something other than -100.0, then I know that it has been computed.
      ! Trying to get at the differece between computed bad values and segments that have not been
      ! computed yet.
      seg_tave_upstream(i) = -100.0

! Big do loop
      DO j = 1, Nsegment
         i = Segment_order(j)

         ! !! LOOP BREAKS HERE !!
         !
         ! If the seg_tave_water value has been set to -99.9 (in init), then this is a segment that will
         ! never have streamflow because it does not have any HRUs connected to it and none of the
         ! upstream segments (if there are any) have HRUs connected. Because there can never be any
         ! flow, the temperature calculation will always fail, so don't bother with it.
         if (Seg_tave_water(i) < -99.0) then
            cycle
         endif

         ! !! LOOP BREAKS HERE !!
         !
         ! If the seginc_swrad value has been set to -99.9 (route_run), then this segment will
         ! never have solar radiation because it does not have any HRUs connected to it and none of the
         ! upstream or downstream segments have HRUs connected.
         if (seginc_swrad(i) < -99.0) then
            Seg_tave_water(i) = -99.9
            cycle
         endif

! GW moving average
         gw_sum(i) = gw_sum(i) - gw_silo(i, gw_index)
         gw_silo(i, gw_index) = Seg_tave_air(i)
         gw_sum(i) = gw_sum(i) + gw_silo(i, gw_index)
         seg_tave_gw(i) = gw_sum(i) / gw_tau(i)

! SS moving average
         ss_sum(i) = ss_sum(i) - ss_silo(i, ss_index)
         ss_silo(i, ss_index) = Seg_tave_air(i)
         ss_sum(i) = ss_sum(i) + ss_silo(i, ss_index)
         seg_tave_ss(i) = ss_sum(i) / ss_tau(i)

! Find upstream intitial inflow temperature for segment i
! i is the current segment
! kk is the upstream segment
         fs = 0.0
         up_temp = 0.0
         do k = 1, upstream_count(i)
            kk = upstream_idx(i,k)
            if (Seg_tave_water(kk) > -1.0) then
               up_temp = up_temp + (Seg_tave_water(kk) * SNGL(Seg_outflow(kk)))
               fs = fs + SNGL(Seg_outflow(kk))
            ENDIF
         ENDDO

         ! Finish computing seg_tave_upstream
         IF ( fs > NEARZERO) THEN
            seg_tave_upstream(i) = up_temp / fs
         ELSE
            ! -98.9 is the code for no flow on this timestep
            seg_tave_upstream(i) = NOFLOW_TEMP
         ENDIF

! debug
!         if (seg_tave_upstream(i) > 100.0) then
!            write(*,*) "upstream_temp: i = ", i, " seg_tave_upstream = ", seg_tave_upstream(i), " fs = ", &
!      &        fs, " seg_tave_water = ", Seg_tave_water(i), " troff = " , Seg_tave_air(i), " up_temp = ", up_temp
!         endif

         ! Compute flow-dependent water-in-segment width value
         if (seg_outflow(i) > NEARZERO) then
            Seg_width(i) = width_alpha(i) * SNGL(Seg_outflow(i)) ** width_m(i)
         else
            Seg_width(i) = 0.0
            if (Seg_tave_water(i) > -99.0) then
               ! This segment has upstream HRUs somewhere, but the current day's flow is zero
               Seg_tave_water(i) = NOFLOW_TEMP
            endif
         endif

         ! Compute the shade on the segment. Either set by value in the parameter file or computed
         IF ( Stream_temp_shade_flag==1 ) THEN
            IF ( Summer_flag==0 ) THEN
               seg_shade(i) = Segshade_win(i)
            ELSE
               seg_shade(i) = Segshade_sum(i)
            ENDIF

           ! Svi    = RIPARIAN VEGETATION SHADE
            svi = 0.0
         ELSE
            CALL shday(i, seg_shade(i), svi)
         ENDIF

         ! Start working towards the computation of the equilibrium temperature
         qlat = 0.0D0
         seg_tave_lat(i) = 0.0
         ak1 = 0.0
         ak2 = 0.0

         ! Inputs: seg_tave_gw, Seg_tave_air, seg_tave_ss, seg_tave_upstream, Seg_melt, Seg_rain
         ! Outputs: qlat (in CMS), seg_tave_lat
         CALL lat_inflow(qlat, seg_tave_lat(i), i, seg_tave_gw(i), Seg_tave_air(i), seg_tave_ss(i), &
     &                   Seg_melt(i), Seg_rain(i))


         ! This code does not handle thermodynamics of ice, so temperatures below 0 are not allowed.
         ! The question is when to set temperatures below 0 to 0. If, after computing the running averages
         ! and mixing the different sources of lateral flow, the temperature is less than 0, set the lateral
         ! flow temperature to 0 here.
         if (seg_tave_lat(i) .lt. NEARZERO) then
            seg_tave_lat(i) = 0.0
         endif

! Compute t_o
! t_o is the temperature of the water at the beginning of the time step (this is To in equation 32)
         if (Seg_tave_water(i) < -99.0) then
!            No flow in this segment and there never will be becuase there are no upstream HRUs.
            t_o = Seg_tave_water(i)

         elseif (Seg_tave_water(i) < -98.0) then
!            No flow in this segment on this time step, but could be on future time step
            t_o = Seg_tave_water(i)

         elseif ((fs .le. NEARZERO) .and. (qlat .le. NEARZERO)) then
             ! If there is no flow, set the temperature to -98.9
             ! -99.9 means that the segment never has any flow (determined up in init).
             ! -98.9 means that this a segment that could have flow, but doesn't
            Seg_tave_water(i) = -98.9
            t_o = Seg_tave_water(i)

         elseif (fs .le. NEARZERO) then
             ! if this is true, then there is no flow from upstream, but there is lateral inflow
            t_o = seg_tave_lat(i) + lat_temp_adj(i,Nowmonth)

         elseif (qlat .le. NEARZERO) then
             ! if this is true, then there is no lateral flow, but there is flow from upstream
            t_o = seg_tave_upstream(i)

         else
             ! if this is true, then there is both lateral flow and flow from upstream
             !  qlat is in CMS so fs needs to be converted
            t_o = sngl((seg_tave_upstream(i) * fs * CFS2CMS_CONV) + &
     &                   (sngl(qlat) * (seg_tave_lat(i) + lat_temp_adj(i,Nowmonth)))) / &
     &                   sngl((fs * CFS2CMS_CONV) + sngl(qlat))
         endif

! debug
!         if (t_o .ne. t_o) then
!             write(*,*) "t_o is Nan, seg_tave_upstream = ", seg_tave_upstream(i), " fs = ", fs, &
!     &                    " qlat = ", qlat, " seg_tave_lat = ", seg_tave_lat(i), " lat_temp_adj = ", lat_temp_adj(i,Nowmonth)
!             continue
!         endif

! debug
!         if (t_o .gt. 100.0) then
!             write(*,*) "this is the place: t_o = ", t_o, " ted = ", te, " seg_id = ", i
!             write(*,*) "   seg_tave_upstream = ", seg_tave_upstream(i), " fs = ", fs, &
!     &                    " qlat = ", qlat, " seg_tave_lat = ", seg_tave_lat(i), " lat_temp_adj = ", lat_temp_adj(i,Nowmonth)
!             write(*,*) "   width = ", Seg_width(i), Nowyear, Nowmonth, Nowday
!             continue
!             exit
!          endif

!         Need a good value of t_o
          if (t_o .gt. -98.0) then
!             This block computes the value for seg_tave_water

!             Compute the equilibrium temerature
              ! Out: te, ak1, ak2
              ! In: seg_shade, svi, i, t_o
              CALL equilb(te, ak1, ak2, seg_shade(i), svi, i, t_o)

!             Compute the daily mean water temperature
              ! In: t_o, qlat, seg_tave_lat(i), te, ak1, ak2, i, seg_width, seg_length
              Seg_tave_water(i) = twavg(fs, t_o, qlat, seg_tave_lat(i), te, ak1, ak2, seg_width(i), seg_length(i))

          else
              ! bad t_o value
              Seg_tave_water(i) = NOFLOW_TEMP
          endif

!          if (Seg_tave_water(i) .ne. Seg_tave_water(i)) then
!             write(*,*) "seg_tave_water is NaN", i, qlat, seg_tave_lat(i), te, ak1, ak2,seg_shade(i), svi, i, t_o
!          endif

      ENDDO
      END FUNCTION stream_temp_run
!
!*********************************************************************************
! Compute the flow-weighted average temperature and a total sum of lateral inflows
!*********************************************************************************
      SUBROUTINE lat_inflow(Qlat, Tl_avg, id, tave_gw, tave_air, tave_ss, melt, rain)
      USE PRMS_CONSTANTS, ONLY: CFS2CMS_CONV, NEARZERO
      USE PRMS_STRMTEMP, ONLY: Melt_temp
      USE PRMS_FLOWVARS, ONLY: Seg_lateral_inflow
      USE PRMS_ROUTING, ONLY: Seginc_sroff, Seginc_ssflow, Seginc_gwflow
      IMPLICIT NONE
! Functions
      INTRINSIC :: SNGL
! Arguments
      INTEGER, INTENT(IN) :: id
      REAL, INTENT(IN) :: tave_gw, tave_air, tave_ss, melt, rain
      REAL, INTENT(OUT) :: Tl_avg
      DOUBLE PRECISION, INTENT(OUT) :: Qlat
! Local Variables
      REAL :: weight_roff, weight_ss, weight_gw, melt_wt, rain_wt, troff, tss
!*****************************************************************************

      Qlat = Seg_lateral_inflow(id) * CFS2CMS_CONV
      Tl_avg = 0.0
      IF ( Qlat>0.0D0 ) THEN ! weights do not include water-use if active, not sure it works for cascades
         weight_roff = SNGL( (Seginc_sroff(id) / Qlat) * CFS2CMS_CONV )
         weight_ss = SNGL( (Seginc_ssflow(id) / Qlat) * CFS2CMS_CONV )
         weight_gw = SNGL( (Seginc_gwflow(id) / Qlat) * CFS2CMS_CONV )
      ELSE
         weight_roff = 0.0
         weight_ss = 0.0
         weight_gw = 0.0
      ENDIF

      IF (melt > 0.0) THEN
         melt_wt = melt/(melt + rain)
         IF (melt_wt < 0.0) melt_wt = 0.0
         IF (melt_wt > 1.0) melt_wt = 1.0
         rain_wt = 1.0 - melt_wt
         IF (rain == 0.0) THEN
            troff = Melt_temp
            tss = Melt_temp
         ELSE
            troff = Melt_temp * melt_wt + tave_air * rain_wt
            tss = Melt_temp * melt_wt + tave_ss * rain_wt
         ENDIF
      ELSE
         troff = tave_air
         tss = tave_ss
      ENDIF

      Tl_avg = weight_roff * troff + weight_ss * tss + weight_gw * tave_gw

      END SUBROUTINE lat_inflow

!***********************************************************************************************
      REAL FUNCTION twavg(qup, T0, Qlat, Tl_avg, Te, Ak1, Ak2, width, length)
!
!     PURPOSE:
!        1. TO PREDICT THE AVERAGE DAILY WATER TEMPERATURE USING A SECOND-ORDER
!           CLOSED-FORM SOLUTION TO THE STEADY-STATE HEAT TRANSPORT EQUATION.
      USE PRMS_CONSTANTS, ONLY: NEARZERO, CFS2CMS_CONV
      IMPLICIT NONE
! Functions
      INTRINSIC :: ABS, EXP, ALOG, SNGL, SIGN
! Arguments
      REAL, INTENT(IN) :: T0, Tl_avg, Te, Ak1, Ak2, width, length, qup
      DOUBLE PRECISION, INTENT(IN) :: Qlat
! Local Variables
      REAL :: tep, b, r, rexp, tw, delt, q_init, denom, Ql
!***************************************************************************************************
! DETERMINE EQUATION PARAMETERS
      q_init = SNGL(qup  * CFS2CMS_CONV)
      Ql = SNGL( Qlat )

! This is confused logic coment out here and compute the terms as needed below
!      b = (Ql / Seg_length) + ((Ak1 * Seg_width) / 4182.0E03)
!      IF ( b < NEARZERO ) b = NEARZERO ! rsr, don't know what value this should be to avoid divide by 0
!      r = 1.0 + (Ql / q_init)
!      IF ( r < NEARZERO ) r = NEARZERO

      IF (Ql <= NEARZERO ) THEN
!
! ZERO LATERAL FLOW
         tep = Te
         b = (Ak1 * width) / 4182.0E03
         rexp = -1.0*(b * length) / q_init
         r = EXP(rexp)

! LOSING STREAM
! No such thing as losing streams in PRMS
      ELSEIF ( Ql < 0.0 ) THEN
         write(*,*) "twavg: losing stream!!! Should be no such thing in PRMS!"
         tep  = Te
         b = (Ql / length) + ((Ak1 * width) / 4182.0E03)
         rexp = (Ql - (b * length)) / Ql
         r = 1.0 + (Ql / q_init)
         r = r**rexp
!
! This is a headwaters (i.e. no streamflow from above, but lateral flow from HRUs.
! Treat the lateral flow as upstream flow to avoid divide by zero
      ELSEIF ( Ql > NEARZERO .and. q_init <= NEARZERO ) THEN
         tep = Te
         b = (Ak1 * width) / 4182.0E03
!         rexp = -1.0*(b * length) / q_init
         rexp = -1.0*(b * length) / Ql
         r = EXP(rexp)
!
! GAINING STREAM (ie both ql and q_init have > zero values)
      ELSE
         b = (Ql / length) + ((Ak1 * width) / 4182.0E03)
         tep = (((Ql / length) * Tl_avg) + (((Ak1 * width) / (4182.0E03)) * Te)) / b

! shouldn't need to do this because Ql will always be greater than 0 if in here.
         IF ( Ql > 0.0 ) THEN
            rexp = -b / (Ql / length)
         ELSE
            rexp = 0.0
         ENDIF

! DANGER -- replaced this potential divide by zero with the logic below
!          r = 1.0 + (Ql / q_init)
         if (q_init < NEARZERO) then
            r = 2.0
         else
            r = 1.0 + (Ql / q_init)
         endif
         r = r**rexp

! END LATERAL FLOW TERM LOGIC
      ENDIF
!
! DETERMINE WATER TEMPERATURE
      delt  = tep - T0
      denom = (1.0 + (Ak2 / Ak1) * delt * (1.0 - r))
      IF ( ABS(denom) < NEARZERO ) denom = SIGN(NEARZERO, denom)
      tw    = tep - (delt * r / denom)
      IF ( tw < 0.0 ) tw = 0.0

      twavg = tw
      END FUNCTION twavg
!
!*******************************************************************************
!    "equilb"
!*******************************************************************************
      SUBROUTINE equilb (Ted, Ak1d, Ak2d, Sh, Svi, Seg_id, t_o)
!
!     PURPOSE:
!        1. DETERMINE THE AVERAGE DAILY EQUILIBRIUM WATER TEMPERATURE PARAMETERS
!        2. DETERMINE THE MAXIMUM DAILY EQUILIBRIUM WATER TEMPERATURE PARAMETERS

      USE PRMS_CONSTANTS, ONLY: NEARZERO, CFS2CMS_CONV
      USE PRMS_STRMTEMP, ONLY: ZERO_C, Seg_width, Seg_humid, Press, MPS_CONVERT, &
     &    Seg_ccov, Seg_potet, Albedo, seg_tave_gw, Seg_slope
      USE PRMS_FLOWVARS, ONLY: Seg_inflow
      USE PRMS_ROUTING, ONLY: Seginc_swrad
      IMPLICIT NONE
! Functions
      INTRINSIC :: EXP, SQRT, ABS, SNGL, DBLE
      EXTERNAL :: teak1
      REAL, EXTERNAL :: sat_vapor_press_poly
! Arguments:
      REAL, INTENT(OUT) :: Ted
      REAL, INTENT(OUT) :: Ak1d, Ak2d
      REAL, INTENT(IN) :: Sh, Svi
      INTEGER, INTENT(IN) :: Seg_id
      REAL, INTENT(IN) :: t_o
! Local Variables:  !RSR, maybe declare enegry balance fluxes
      DOUBLE PRECISION :: ha, hv, taabs
      REAL :: hf, hs, b, c, d, delt, del_ht, ltnt_ht, bow_coeff
      REAL :: hnet, vp_sat, sw_power, evap, q_init
      REAL, PARAMETER :: AKZ = 1.65, A = 5.40E-8, RAD_CONVERT = 41840.0/86400.0
      REAL :: foo
! *******************************************************************************

      taabs = DBLE( t_o + ZERO_C )
      vp_sat = 6.108 * EXP(17.26939 * t_o/(t_o + 237.3))

!
!  Convert units and set up parameters
      q_init = SNGL( Seg_inflow(Seg_id) * CFS2CMS_CONV )
      IF ( q_init < NEARZERO ) q_init = NEARZERO

      ! sw_power should be in watts / m2
      ! seginc_swrad is in Langley / day
      ! Used to use RAD_CONVERT, the conversion I'm using now is a slightly different number.
      sw_power = 11.63 / 24.0 * SNGL(seginc_swrad(seg_id))

      del_ht = 2.36E06   ! could multiple by 10E6 for this and other terms later to reduce round-off
      ltnt_ht = 2495.0E06

! If humidity is 1.0, there is a divide by zero below.
      if (Seg_humid(Seg_id) > 0.99) then
          foo = 0.99
      else
          foo = Seg_humid(Seg_id)
      endif

      bow_coeff = (0.00061 * Press(Seg_id))/(vp_sat * (1.0 - foo))
      evap = SNGL( Seg_potet(Seg_id) * MPS_CONVERT )
!
! HEAT FLUX COMPONENTS
      ! document - ha = (1-rl)(1-sh)(1+0.17Cl**2)(0.61+0.05*SQRT(vp_sat)*stefan(Ta+273.16)**4

      ha = ( (3.354939D-8 + 2.74995D-9 * DBLE(SQRT(Seg_humid(Seg_id) * vp_sat))) * DBLE((1.0 - Sh) &
     &     * (1.0 + (0.17*(Seg_ccov(Seg_id)**2)))) ) * (taabs**4)

! hf is heat from stream friction. See eqn. 14.  q_init is in CMS
      hf = 9805.0 * (q_init/Seg_width(Seg_id)) * Seg_slope(Seg_id)
      hs = (1.0 - sh) * sw_power * (1.0 - Albedo)
      hv = 5.24D-8 * DBLE(Svi) * (taabs**4)

! Stefan-Boltzmann constant = 5.670373D-08; emissivity of water = 0.9526, times each other: 5.4016D-08
! hw = water-emitted longwave radiation
! hw = 5.4016D-08 * (taabs**4)  hw is include in other computations
!
! DETERMINE EQUILIBIRIUM COEFFICIENTS
      b = bow_coeff * evap * (ltnt_ht + (del_ht * t_o)) + AKZ - (del_ht * evap)
      c = bow_coeff * del_ht * evap
      d = (SNGL(ha + hv) + hf + hs) + (ltnt_ht * evap * ((bow_coeff * t_o) - 1.0) + (seg_tave_gw(Seg_id) * AKZ))

!
! DETERMINE EQUILIBRIUM TEMPERATURE & 1ST ORDER THERMAL EXCHANGE COEF.
      Ted = t_o

      CALL teak1(A, b, c, d, Ted, Ak1d)

!      if (ted .ne. ted) then
!         write(*,*) "   ted is NaN in equilib ", seg_id, A, b, c, d, Ted, Ak1d
!         write(*,*) "   ", seg_id, ha, hv, hf, hs, ltnt_ht, evap, bow_coeff, seg_tave_gw(Seg_id), AKZ
!         write(*,*) "   ha components ",  Seg_humid(Seg_id), vp_sat, Sh, Seg_ccov(Seg_id), taabs
!      endif

!
! DETERMINE 2ND ORDER THERMAL EXCHANGE COEFFICIENT
      hnet = (A * ((t_o + ZERO_C)**4)) + (b * t_o) - (c * (t_o**2.0)) - d
      delt = t_o - Ted

      IF ( ABS(delt) < NEARZERO) THEN
        Ak2d = 0.0
      ELSE
        Ak2d = ((delt * Ak1d) - hnet) / (delt**2)
      ENDIF
!
! RETURN TO STREAMTEMP FUNCTION
      END SUBROUTINE equilb

!**********************************************************************************
!    "teak1"
!**********************************************************************************
      SUBROUTINE teak1(A, B, C, D, Teq, Ak1c)
!     PURPOSE:
!        1. TO DETERMINE THE EQUILIBRIUM WATER TEMPERATURE FROM THE ENERGY BALANCE
!           EQUATION BY ITERATING NEWTON'S METHOD
!        2. TO DETERMINE THE 1ST THERMAL EXCHANGE COEFFICIENT.
      USE PRMS_STRMTEMP, ONLY: ZERO_C, Maxiter_sntemp
!      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Nowday
      IMPLICIT NONE
      INTRINSIC :: ABS
! Arguments
      REAL, INTENT(IN) ::  A, B, C, D
      REAL, INTENT(INOUT) :: Teq
      REAL, INTENT(OUT) :: Ak1c
! Local variables
      REAL :: teabs, fte, fpte, delte
      INTEGER :: kount
! Parameters
      ! SOLUTION CONVERGENCE TOLERANCE
      REAL, PARAMETER :: TOLRN = 1.0E-4
!**********************************************************************************
      fte = 99999.0 ! rsr, fte was not set
      delte = 99999.0 ! rsr, delte was not set
      kount = 0

! BEGIN NEWTON ITERATION SOLUTION FOR TE
      DO kount = 1, Maxiter_sntemp
         IF ( ABS(fte) < TOLRN ) EXIT
         IF ( ABS(delte) < TOLRN ) EXIT
         teabs  = Teq + ZERO_C
         fte    = (A * (teabs**4.0)) + (B * Teq) - (C * (Teq**2.0)) - D
         fpte   = (4.0 * A * (teabs**3.0)) + B - (2.0 * C * Teq)
         delte  = fte / fpte
         Teq    = Teq - delte
      ENDDO

! DETERMINE 1ST THERMAL EXCHANGE COEFFICIENT
      Ak1c = (4.0 * A * ((Teq + ZERO_C)**3.0)) + B - (2.0 * C * Teq)
!
! RETURN TO 'EQUILB' SUBROUTINE
      END SUBROUTINE teak1

!     "shday" ***********************************************************
      SUBROUTINE shday(Seg_id, Shade, Svi)
!
!      THIS SUBPROGRAM IS TO CALCULATE THE TOTAL DAILY SHADE FOR A
!  GIVEN REACH. BOTH TOPOGRAPHIC AND RIPARIAN VEGETATION SHADE
!  IS INCLUDED.
!
!  VARIABLE NAME LIST
!
!      Als    = CURRENT SOLAR ALTITUDE
!      Alrs   = SOLAR ALTITUDE WHEN SOLAR & REACH AZIMUTHS ARE EQUAL
!      Alsmx  = MAXIMUM POSSIBLE SOLAR ALTITUDE
!      Alsr   = LOCAL SUNRISE SOLAR ALTITUDE ! rsr, not used
!      Alss   = LOCAL SUNSET SOLAR ALTITUDE ! rsr, not used
!      Alt    = CURRENT TOPOGRAPHIC ALTITUDE
!      Alte   = EAST SIDE MAXIMUM TOPOGRAPHIC ALTITUDE
!      Altmx  = CURRENT MAXIMUM TOPOGRAPHIC ALTITUDE LIMIT
!      Altop  = CURRENT TOPOGRAPHIC ALTITUDE
!      Altw   = WEST SIDE MAXIMUM TOPOGRAPHIC ALTITUDE
!      Azrh   = STREAM REACH AZIMUTH
!      Azs    = CURRENT SOLAR AZIMUTH
!      Azsr   = LOCAL SUNRISE SOLAR AZIMUTH ! rsr, not used
!      Azss   = LOCAL SUNSET SOLAR AZIMUTH ! rsr, not used
!      Azso   = LEVEL-PLAIN SUNSET AZIMUTH
!      Bavg   = AVERAGE STREAM WIDTH
!      Bs     = SHADED PART OF STREAM WIDTH
!      Cosas  = COS(AS)
!      Cosd   = COS(DECL)
!      Coshs  = COS(HS)
!      Coso   = COS(XLAT)
!      Cosod  = COS(XLAT)*COS(DECL)
!      Dayrad = CONVERSION RATIO FOR JULIAN DAYS TO RADIANS
!      Decl   = CURRENT SOLAR DECLINATION
!      Delhsr = SUNRISE SIDE HOUR ANGLE INCREMENT
!      Delhss = SUNSET SIDE HOUR ANGLE INCREMENT
!      Hrrs   = REACH HOUR ANGLE WHEN SOLAR & REACH AZIMUTHS ARE EQUAL
!      Hrs    = CURRENT SOLAR HOUR ANGLE
!      Hrsr   = LOCAL SUNRISE SOLAR HOUR ANGLE
!      Hrss   = LOCAL SUNSET SOLAR HOUR ANGLE
!      Hrso   = LEVEL-PLAIN SUNRISE/SET SOLAR HOUR ANGLE
!      Nbhs   = NUMBER OF SUNRISE/SET HOUR ANGLE INCREMENTS
!      Shday  = TOTAL DAILY SHADE
!      Sinal  = SIN(Al)
!      Sinar  = SIN(Ar)
!      Sin_d   = SIN(DECL)
!      Sinhsr = SIN(Hrsr)
!      Sinhss = SIN(Hrss)
!      Sinhro = SIN(Hrso)
!      Sino   = SIN(XLAT)
!      Sinod  = SIN(XLAT)*SIN(DECL)
!      Snflag = SOLAR NOON LIMIT FLAG
!      Sti    = TOPOGRAPHIC SHADE
!      Svi    = RIPARIAN VEGETATION SHADE
!      Svri   = SUNRISE VEGETATIVE SHADE
!      Svsi   = SUNSET VEGETATIVE SHADE
!      Tanasr = TAN(Alsr)
!      Tanass = TAN(Alss)
!      Tanalt = TAN(Alt)
!      Tano   = TAN(XLAT)
!      Tanod  = TAN(XLAT)*TAN(DECL)
!      Totsh  = LEVEL-PLAIN TOTAL SHADE POTENTIAL
!      Tolrn  = CONVERGENCE TOLERANCE CRITERIA
!      Flgrs  = SUNRISE FLAG; TRUE IF SUNRISE, FALSE IF SUNSET
!      Flgst  = SUNSET FLAG; TRUE IF SUNSET, FALSE IF SUNRISE
!      Vc = CROWN DIAMETER, CURRENT VEGETATION
!      Vce    = CROWN DIAMETER, EAST SIDE VEGETATION
!      Vco    = CURRENT VEGETATION OVERHANG
!      Vcw    = CROWN DIAMETER, WEST SIDE VEGETATION
!      Vd     = DENSITY, CURRENT VEGETATION
!      Vde    = DENSITY, EAST SIDE VEGETATION
!      Vdw    = DENSITY, WEST SIDE VEGETATION
!      Vh     = HEIGHT, CURRENT VEGETATION
!      Vhe    = HEIGHT, EAST SIDE VEGETATION
!      Vhw    = HEIGHT, WEST SIDE VEGETATION
!      Vo     = OFFSET, CURRENT VEGETATION
!      Voe    = OFFSET, EAST SIDE VEGETATION
!      Vow    = OFFSET, WEST SIDE VEGETATION
!
      USE PRMS_SET_TIME, ONLY: Jday
      USE PRMS_STRMTEMP, ONLY: Azrh, Alte, Altw, Seg_daylight, Seg_width, &
     &    PI, HALF_PI, Cos_seg_lat, Sin_seg_lat, Cos_lat_decl, Horizontal_hour_angle, &
     &    Level_sunset_azimuth, Max_solar_altitude, Sin_alrs, Sin_declination, Sin_lat_decl, Total_shade
      IMPLICIT NONE
! Functions
      INTRINSIC :: COS, SIN, TAN, ACOS, ASIN, ATAN, ABS, MAX, SNGL
      REAL, EXTERNAL:: solalt, rprnvg
      EXTERNAL :: snr_sst
! Arguments
      INTEGER, INTENT(IN) :: Seg_id
      REAL, INTENT(OUT):: Shade, Svi
! Local Variables
      REAL :: coso, cosod, sin_d, sino, sinod
      REAL :: altmx, alsmx, als, almn, almx
      REAL :: azso, azmn, azmx, azs, hrrs, hrsr, hrss, hrso, hrs, hrrh
      REAL :: temp, totsh, sti
      REAL :: altop(3), aztop(3)
! PARAMETER
      REAL, PARAMETER :: RADTOHOUR = 24.0/(2.0 * PI)
!*********************************************************************************

!  LATITUDE TRIGONOMETRIC PARAMETERS
      coso = Cos_seg_lat(Seg_id)
      sino = Sin_seg_lat(Seg_id)
      sin_d = Sin_declination(Jday, Seg_id)
      sinod = Sin_lat_decl(Jday, Seg_id)
      cosod = Cos_lat_decl(Jday, Seg_id)
!
!  INITIALIZE LOCAL SUNRISE/SET SOLAR PARAMETERS
      hrsr = 0.0
      hrss = 0.0
!
!  MAXIMUM POSSIBLE SOLAR ALTITUDE
      alsmx = Max_solar_altitude(Jday, Seg_id)
!
!  LEVEL-PLAIN SUNRISE/SET HOUR ANGLE
      hrso = Horizontal_hour_angle(Jday, Seg_id)
!
!  LEVEL-PLAIN SOLAR AZIMUTH
      azso = Level_sunset_azimuth(Jday, Seg_id)
!
!  TOTAL POTENTIAL SHADE ON LEVEL-PLAIN
      totsh = Total_shade(Jday, Seg_id)
!
!  CHECK FOR REACH AZIMUTH LESS THAN SUNRISE
      IF ( Azrh(Seg_id) <= (-azso) ) THEN
         hrrs = -hrso
!
!  CHECK FOR REACH AZIMUTH GREATER THAN SUNSET
      ELSEIF ( Azrh(Seg_id) >= azso ) THEN
         hrrs = hrso
!
!  REACH AZIMUTH IS BETWEEN SUNRISE & SUNSET
      ELSEIF ( Azrh(Seg_id) == 0.0 ) THEN
         hrrs = 0.0
      ELSE
         temp = (Sin_alrs(Jday, Seg_id) - sinod) / cosod
         IF ( ABS(temp) > 1.0 ) temp = SIGN(1.0,temp)
         hrrs = SIGN(ACOS(temp), Azrh(Seg_id))
!
!  END REACH & SOLAR AZIMUTH CHECK
      ENDIF
!
!  CHECK IF LEVEL-PLAIN
      IF ( (Alte(Seg_id) == 0.0 ) .AND. (Altw(Seg_id) == 0.0) ) THEN
!         azsr = -azso
         hrsr = -hrso
!         azss = azso
         hrss = hrso
         sti = 0.0
         Svi = (rprnvg(hrsr, hrrs, hrss, sino, coso, sin_d, cosod, sinod, Seg_id)) / (Seg_width(Seg_id) * totsh)

      ELSE
!  INITIALIZE SHADE VALUES
!
!  INSERT STARTING TOPOGRAPHIC AZIMUTH VALUES BETWEEN LEVEL PLAIN SUNRISE AND SUNSET
         aztop = 0.0
!
!  DETERMINE SUNRISE HOUR ANGLE.
         altop = 0.0
         IF ( -azso <= Azrh(Seg_id) ) THEN
            altop(1) = Alte(Seg_id)
            aztop(1) = azso*(Alte(Seg_id)/HALF_PI) - azso
         ELSE
            altop(1) = Altw(Seg_id)
            aztop(1) = azso*(Altw(Seg_id)/HALF_PI) - azso
         ENDIF
! LEVEL PLAIN
         IF (altop(1) == 0.0) THEN
            hrsr = -hrso
! NOT
         ELSE
! LOOK FOR SOLUTION BETWEEN LIMITS OF LEVEL PLAIN SUNRISE AND NOON
            azmn = -azso
            azmx = 0.0
            azs = aztop(1)
            altmx = altop(1)
            almn = 0.0
            almx = 1.5708
            als = solalt(coso, sino, sin_d, azs, almn, almx)
            CALL snr_sst(coso, sino, sin_d, altmx, almn, almx, azmn, azmx, azs, als, hrs, Seg_id)
!            azsr = azs
!            alsr = als
            hrsr = hrs
!            altr = altmx
         ENDIF
!
!  DETERMINE SUNSET HOUR ANGLE.
         IF ( azso <= Azrh(Seg_id) )THEN
            altop(2) = Alte(Seg_id)
            aztop(2) = azso - azso*(Alte(Seg_id)/HALF_PI)
         ELSE
            altop(2) = Altw(Seg_id)
            aztop(2) = azso - azso*(Altw(Seg_id)/HALF_PI)
         ENDIF
! LEVEL PLAIN
         IF (altop(2) == 0.0) THEN
            hrss = hrso
! NOT
         ELSE
! LOOK FOR SOLUTION BETWEEN LIMITS OF NOON AND LEVEL PLAIN SUNSET
            azmn = 0.0
            azmx = azso
            azs = aztop(2)
            altmx = altop(2)
            almn = 0.0
            almx = 1.5708
            als = solalt(coso, sino, sin_d, azs, almn, almx)
            CALL snr_sst(coso, sino, sin_d, altmx, almn, almx, azmn, azmx, azs, als, hrs, Seg_id)
!            azss = azs
!            alss = als
            hrss = hrs
!            alts = altmx
         ENDIF
!
!  SOLVE FOR SHADE INCREMENTS THIS SEGMENT
         IF ( hrrs < hrsr ) THEN
            hrrh = hrsr
         ELSEIF ( hrrs > hrss ) THEN
            hrrh = hrss
         ELSE
            hrrh = hrrs
         ENDIF

         Seg_daylight(Seg_id) = (hrss - hrsr) * RADTOHOUR
         sti = 1.0 - ((((hrss - hrsr) * sinod) + ((SIN(hrss) - SIN(hrsr)) * cosod)) / (totsh))
         Svi = ((rprnvg(hrsr, hrrh, hrss, sino, coso, sin_d, cosod, sinod, Seg_id)) / (Seg_width(Seg_id)*totsh))
!
!  END SUNRISE/SUNSET CALCULATION
      ENDIF
!
!  CHECK FOR ROUNDOFF ERRORS
      IF ( sti < 0.0 ) sti = 0.0
      IF ( sti > 1.0 ) sti = 1.0
      IF ( Svi < 0.0 ) Svi = 0.0
      IF ( Svi > 1.0 )  Svi = 1.0
!
!  RECORD TOTAL SHADE
      Shade = sti + Svi

      END SUBROUTINE shday
!
!**********************************************************************************************************
!    "snr_sst"
      SUBROUTINE snr_sst (Coso, Sino, Sin_d, Alt, Almn, Almx, Azmn, Azmx, Azs, Als, Hrs, Seg_id)
!
!     THIS SUBPROGRAM DETERMINES THE LOCAL SOLAR SUNRISE/SET
! AZIMUTH, ALTITUDE, AND HOUR ANGLE
!
      USE PRMS_CONSTANTS, ONLY: NEARZERO
      USE PRMS_STRMTEMP, ONLY: Azrh, PI, Maxiter_sntemp
      IMPLICIT NONE
! Functions
      INTRINSIC :: TAN, SIN, COS, ACOS, ASIN, ABS
!  Arguments
      INTEGER, INTENT(IN):: Seg_id
      REAL, INTENT(IN):: Coso, Sino, Sin_d, Alt, Almn, Almx, Azmn, Azmx
      REAL, INTENT(INOUT):: Azs, Als
      REAL, INTENT(OUT):: Hrs
!  Local Variables
      REAL :: cosazs, sinazs, sinazr, cosazr, cosals, f, g, fazs, fals, gazs, gals, xjacob
      REAL :: sinals, tanalt, tano, tanals, temp, delazs, delals
      INTEGER :: count
!***********************************************************************************************************
! TRIG FUNCTION FOR LOCAL ALTITUDE
      tanalt = TAN(Alt)
      tano   = Sino / Coso
      f = 999999.0 !rsr, these need values
      delazs = 9999999.0
      g = 99999999.0
      delals = 99999999.0

! BEGIN NEWTON-RAPHSON SOLUTION
      DO count = 1, Maxiter_sntemp
         IF ( ABS(delazs) < NEARZERO ) EXIT
         IF ( ABS(delals) < NEARZERO ) EXIT
         IF ( ABS(f) < NEARZERO ) EXIT
         IF ( ABS(g) < NEARZERO ) EXIT

         cosazs = COS(Azs)
         sinazs = SIN(Azs)

         sinazr = ABS(SIN(Azs - Azrh(Seg_id)))
         IF ( (((Azs-Azrh(Seg_id)) <= 0.0 ) .AND. ((Azs-Azrh(Seg_id)) <= (-PI))) .OR. &
     &       (((Azs-Azrh(Seg_id)) > 0.0 ) .AND. ((Azs-Azrh(Seg_id)) <= PI)) ) THEN
            cosazr = COS(Azs-Azrh(Seg_id))
         ELSE
            cosazr = -COS(Azs-Azrh(Seg_id))
         ENDIF

         cosals = COS(Als)
         IF ( cosals < NEARZERO ) cosals = NEARZERO
         sinals = SIN(Als)
         tanals = sinals / cosals
! FUNCTIONS OF AZS & ALS
         f = cosazs- (((Sino * sinals) - Sin_d) / (Coso * cosals))
         g = tanals - (tanalt * sinazr)
! FIRST PARTIALS DERIVATIVES OF F & G
         fazs = -sinazs
         fals = ((tanals * (Sin_d / Coso)) - (tano / cosals)) / cosals
         gazs = -tanalt * cosazr
         gals = 1.0 / (cosals * cosals)
! JACOBIAN
         xjacob = (fals * gazs) - (fazs * gals)
! DELTA CORRECTIONS
         delazs = ((f * gals) - (g * fals)) / xjacob
         delals = ((g * fazs) - (f * gazs)) / xjacob
! NEW VALUES OF AZS & ALS
         Azs = Azs + delazs
         Als = Als + delals
! CHECK FOR LIMITS
         IF ( Azs < (Azmn + NEARZERO) ) Azs = (Azmn + NEARZERO)
         IF ( Azs > (Azmx - NEARZERO) ) Azs = (Azmx - NEARZERO)
         IF ( Als < (Almn + NEARZERO) ) Als = (Almn + NEARZERO)
         IF ( Als > (Almx - NEARZERO) ) Als = (Almx - NEARZERO)
      ENDDO
!
! ENSURE AZIMUTH REMAINS BETWEEN -PI & PI
      IF ( Azs < (-PI) ) THEN
         Azs = Azs + PI
      ELSEIF ( Azs > PI) THEN
         Azs = Azs - PI
      ENDIF
!
! DETERMINE LOCAL SUNRISE/SET HOUR ANGLE
      sinals = SIN(Als)
      temp = (sinals - (Sino * Sin_d)) / (Coso * COS(ASIN(Sin_d)))
      IF ( ABS(temp) > 1.0 ) temp = SIGN(1.0,temp)
      Hrs = SIGN(ACOS(temp), Azs)

      END SUBROUTINE snr_sst

!*****************************************************************************
!     "solalt"
      REAL FUNCTION solalt (Coso, Sino, Sin_d, Az, Almn, Almx)
!
!      THIS SUBPROGRAM IS TO DETERMINE THE SOLAR ALTITUDE WHEN THE
!  TRIGONOMETRIC PARAMETERS FOR LATITUDE, DECLINATION, AND AZIMUTH
!  ARE GIVEN.
!
!  VARIABLE NAME LIST
!
!      Al     = TRIAL SOLAR ALTITUDE
!      AZ     = SOLAR AZIMUTH
!      COSAL  = COS(AL)
!      COSAZ  = COS(AZ)
!      Coso   = COS(XLAT)
!      DELAL  = INCREMENTAL CORRECTION TO AL
!      FAL    = FUNCTION OF AL
!      FPAL   = FIRST DERIVATIVE OF FAL
!      FPPAL  = SECOND DERIVATIVE OF FAL
!      Sin_d   = SIN(DECL)
!      Sino   = SIN(XLAT)
      USE PRMS_CONSTANTS, ONLY: NEARZERO
      USE PRMS_STRMTEMP, ONLY: HALF_PI, Maxiter_sntemp
      IMPLICIT NONE
! Functions
      INTRINSIC ASIN, ABS, COS, SIN
! Arguments
      REAL, INTENT(IN):: Coso, Sino, Sin_d, Az, Almn, Almx
! Local Variables
      REAL :: cosal, sinal, fal, fpal, fppal, al, alold, delal, a, b, cosaz, temp
      INTEGER :: kount
!*************************************************************************************
!
!  CHECK COS(AZ) EQUAL TO 0
      IF ( ABS(ABS(Az) - HALF_PI) < NEARZERO ) THEN
         temp = ABS(Sin_d / Sino)
         IF ( temp > 1.0 ) temp = 1.0
         Al = ASIN(temp)
      ELSE
!
!  DETERMINE SOLAR ALTITUDE FUNCTION COEFFICIENTS
         cosaz  = COS(Az)
         a      = Sino / (cosaz * Coso)
         b      = Sin_d / (cosaz * Coso)
!
!  INITIALIZE
         al     = (Almn + Almx) / 2.0
         kount = 0
         fal = COS(al) - (a * SIN(al)) + b
         delal = fal/(-SIN(al) - (a * COS(al)))
!
!  BEGIN NEWTON SECOND-ORDER SOLUTION
         DO kount = 1, Maxiter_sntemp
            IF ( ABS(fal) < NEARZERO ) EXIT
            IF ( ABS(delal) < NEARZERO ) EXIT
            alold  = al
            cosal  = COS(al)
            sinal  = SIN(al)
            fal    =  cosal - (a * sinal) + b
            fpal   = -sinal - (a * cosal)
            IF ( kount <= 3 ) THEN
               delal = fal / fpal
            ELSE
               fppal = b - fal
               delal = (2.0 * fal * fpal) / ((2.0 * fpal * fpal) - (fal * fppal))
            ENDIF
            al = al - delal
            IF (al < Almn) al = (alold + Almn) / 2.0
            IF (al > Almx) al = (alold + Almx) / 2.0
         ENDDO
      ENDIF
!
!  SOLUTION OBTAINED
      solalt = al

      END FUNCTION solalt

!***********************************************************************
      REAL FUNCTION rprnvg (Hrsr, Hrrs, Hrss, Sino, Coso, Sin_d, Cosod, Sinod, Seg_id)
!
!      THIS SUBPROGRAM IS TO COMPUTE THE RIPARIAN VEGETATION SHADE
!  SEGMENT BETWEEN THE TWO HOUR ANGLES HRSR & HRSS.
!
      USE PRMS_CONSTANTS, ONLY: NEARZERO
      USE PRMS_STRMTEMP, ONLY: Azrh, Vce, Vdemx, Vhe, Voe, Vcw, Vdwmx, Vhw, Vow, Seg_width, &
     &    Vdemn, Vdwmn, HALF_PI
      USE PRMS_SET_TIME, ONLY: Summer_flag
      IMPLICIT NONE
! Functions
      INTRINSIC :: COS, SIN, ASIN, ACOS, ABS
! Arguments
      REAL, INTENT(IN) :: Hrsr, Hrrs, Hrss, Sino, Coso, Sin_d, Cosod, Sinod
      INTEGER, INTENT(IN):: Seg_id
! Local Variables
      REAL :: svri, svsi, hrs, vco, delhsr, coshrs
      REAL :: sinhrs, temp, als, cosals, sinals, azs, bs, delhss
      INTEGER :: n
! Parameters
      INTEGER, PARAMETER :: NBHS = 15
      DOUBLE PRECISION, SAVE :: Epslon(15), Weight(15)
      DATA Epslon / .006003741, .031363304, .075896109, .137791135, .214513914, &
     &              .302924330, .399402954, .500000000, .600597047, .697075674, &
     &              .785486087, .862208866, .924103292, .968636696, .993996259 /
      DATA Weight / .015376621, .035183024, .053579610, .069785339, .083134603, &
     &              .093080500, .099215743, .101289120, .099215743, .093080500, &
     &              .083134603, .069785339, .053579610, .035183024, .015376621 /
!******************************************************************************
!  ****************** Determine seasonal shade
!
!  CKECK FOR NO SUNRISE
      IF ( Hrsr == Hrss ) THEN
         svri = 0.0
         svsi = 0.0

      ELSE
!
!  VEGETATIVE SHADE BETWEEN SUNRISE & REACH HOUR ANGLES
         svri = 0.0
         IF ( Hrsr < Hrrs ) THEN
            vco = ( Vce(Seg_id)/2.0 ) - Voe(Seg_id)
!
!  DETERMINE SUNRISE SIDE HOUR ANGLE INCREMENT PARAMETERS
            delhsr = Hrrs - Hrsr
!
!  PERFORM NUMERICAL INTEGRATION
            DO n = 1, NBHS
!  CURRENT SOLAR HOUR ANGLE
               hrs = SNGL(Hrsr + (Epslon(n) * delhsr))
               coshrs = COS(hrs)
               sinhrs = SIN(hrs)
!  CURRENT SOLAR ALTITUDE
               temp = Sinod + (Cosod * coshrs)
               IF ( temp > 1.0 ) temp = 1.0
               als = ASIN(temp)
               cosals = COS(als)
               sinals = SIN(als)
               IF ( sinals == 0.0 ) sinals = NEARZERO
!  CURRENT SOLAR AZIMUTH
               temp = ((Sino * sinals) - Sin_d) / (Coso * cosals)
               IF ( ABS(temp) > 1.0 ) temp = SIGN(1.0, temp)
               azs = ACOS(temp)
               IF ( azs < 0.0 ) azs = HALF_PI - azs
               IF ( hrs < 0.0 ) azs = -azs
!  DETERMINE AMOUNT OF STREAM WIDTH SHADED
               bs = ((Vhe(Seg_id) * (cosals/sinals)) * ABS(SIN(azs-Azrh(Seg_id)))) + vco
               IF ( bs < 0.0 ) bs = 0.0
               IF ( bs > Seg_width(Seg_id) ) bs = Seg_width(Seg_id)
!  INCREMENT SUNRISE SIDE VEGETATIVE SHADE
               IF ( Summer_flag == 1 ) THEN ! put back spring and autumn
                  svri = svri + SNGL(Vdemx(Seg_id) * bs * sinals * Weight(n))
               ELSE
                  svri = svri + SNGL(Vdemn(Seg_id) * bs * sinals * Weight(n))
               ENDIF
            ENDDO
!
            svri = svri * delhsr
         ENDIF
!
!  VEGETATIVE SHADE BETWEEN REACH & SUNSET HOUR ANGLES
         svsi = 0.0
         IF ( Hrss > Hrrs ) THEN
            vco = (Vcw(Seg_id)/2.0 ) - Vow(Seg_id)
!
!  DETERMINE SUNSET SIDE HOUR ANGLE INCREMENT PARAMETERS
            delhss = Hrss - Hrrs
!
!  PERFORM NUMERICAL INTEGRATION
            DO n = 1, Nbhs
!  CURRENT SOLAR HOUR ANGLE
               hrs = SNGL(Hrrs + (Epslon(n) * delhss))
               coshrs = COS(hrs)
               sinhrs = SIN(hrs)
!  CURRENT SOLAR ALTITUDE
               temp = Sinod + (Cosod * coshrs)
               IF ( temp > 1.0 ) temp = 1.0
               als = ASIN(temp)
               cosals = COS(als)
               sinals = SIN(als)
               IF ( sinals == 0.0 ) sinals = NEARZERO
!  CURRENT SOLAR AZIMUTH
               temp = ((Sino * sinals) - Sin_d) / (Coso * cosals)
               IF ( ABS(temp) > 1.0 ) temp = SIGN(1.0, temp)
               azs = ACOS(temp)
               IF ( azs < 0.0 ) azs = HALF_PI - azs
               IF ( hrs < 0.0 ) azs = -azs
!  DETERMINE AMOUNT OF STREAM WIDTH SHADED
               bs = ((Vhw(Seg_id) * (cosals/sinals)) * ABS(SIN(azs-Azrh(Seg_id)))) + vco
               IF ( bs < 0.0 ) bs = 0.0
               IF ( bs > Seg_width(Seg_id) ) bs = Seg_width(Seg_id)
!  INCREMENT SUNSET SIDE VEGETATIVE SHADE
               IF ( Summer_flag == 1 ) THEN ! fix for seasons
                  svsi = SNGL(svsi + (Vdwmx(Seg_id) * bs * sinals * Weight(n)))
               ELSE
                  svsi = SNGL(svsi + (Vdwmn(Seg_id) * bs * sinals * Weight(n)))
               ENDIF
            ENDDO
            svsi = svsi * delhss
         ENDIF
      ENDIF

!  COMBINE SUNRISE/SET VEGETATIVE SHADE VALUES
      rprnvg = svri + svsi

      END FUNCTION rprnvg

!***********************************************************************
!     stream_temp_restart - write or read stream_temp restart file
!***********************************************************************
      SUBROUTINE stream_temp_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_STRMTEMP
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=11) :: module_name
!***********************************************************************
      IF ( In_out==0 ) THEN
         WRITE ( Restart_outunit ) MODNAME
         WRITE ( Restart_outunit ) Seg_tave_water
         WRITE ( Restart_outunit ) gw_silo
         WRITE ( Restart_outunit ) ss_silo
         WRITE ( Restart_outunit ) gw_sum
         WRITE ( Restart_outunit ) ss_sum
         WRITE ( Restart_outunit ) gw_index
         WRITE ( Restart_outunit ) ss_index
      ELSE
         READ ( Restart_inunit ) module_name
         CALL check_restart(MODNAME, module_name)
         READ ( Restart_inunit ) Seg_tave_water
         READ ( Restart_inunit ) gw_silo
         READ ( Restart_inunit ) ss_silo
         READ ( Restart_inunit ) gw_sum
         READ ( Restart_inunit ) ss_sum
         READ ( Restart_inunit ) gw_index
         READ ( Restart_inunit ) ss_index
      ENDIF
      END SUBROUTINE stream_temp_restart
