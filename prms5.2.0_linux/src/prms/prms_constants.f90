MODULE PRMS_CONSTANTS
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    !! INT32, REAL32, REAL64
    integer, parameter :: sp = REAL32
    !! Define real precision and range
    integer, parameter :: dp = REAL64
    !! Define double precision and range
    real(REAL32), parameter :: CLOSEZERO = EPSILON(0.0)
    real(REAL32), parameter :: NEARZERO = 1.0E-6
    real(REAL64), parameter :: DNEARZERO = EPSILON(0.0D0)

    integer, parameter :: MAXFILE_LENGTH = 256
    integer, parameter :: MAXCONTROL_LENGTH = 32
    integer, parameter :: MAXDIM = 500
    integer, parameter :: MONTHS_PER_YEAR = 12
    integer, parameter :: MAX_DAYS_PER_YEAR = 366
    integer, parameter :: DAYS_PER_YEAR = 365
    real(REAL32), parameter :: DAYS_YR = 365.242
    real(REAL64), parameter :: DAYS_IN_YEAR = 365.242D0
    real(REAL64), parameter :: SECS_PER_DAY = 86400.0D0
    real(REAL64), parameter :: SECS_PER_HOUR = 3600.0D0
    real(REAL64), parameter :: FT2_PER_ACRE = 43560.0D0
    real(REAL64), parameter :: INCHES_PER_FOOT = 12.0D0
    real(REAL64), parameter :: CFS2CMS_CONV = 0.028316847D0

    real(REAL32), parameter :: INCH2CM = 2.54
    real(REAL32), parameter :: INCH2MM = 25.4
    real(REAL32), parameter :: INCH2M = 0.0254
    real(REAL32), parameter :: MM2INCH = 1.0 / INCH2MM
    real(REAL32), parameter :: FEET2METERS = 0.3048
    real(REAL32), parameter :: METERS2FEET = 1.0 / FEET2METERS

    real(REAL32), parameter :: MAXTEMP = 200.0
    real(REAL32), parameter :: MINTEMP = -150.0

    ! Frequency values, used for basinOut_freq, nhruOut_freq, and nsubOut_freq
    integer, parameter :: DAILY = 1
    integer, parameter :: MONTHLY = 2
    integer, parameter :: DAILY_MONTHLY = 3
    integer, parameter :: MEAN_MONTHLY = 4
    integer, parameter :: MEAN_YEARLY = 5
    integer, parameter :: YEARLY = 6

    ! nowtime index
    integer, parameter :: YEAR = 1
    integer, parameter :: MONTH = 2
    integer, parameter :: DAY = 3
    integer, parameter :: HOUR = 4
    integer, parameter :: MINUTE = 5
    integer, parameter :: SECOND = 6

    ! precip_units
    integer, parameter :: INCHES = 0
    integer, parameter :: MM = 1

    ! temp_units
    integer, parameter :: FAHRENHEIT = 0
    integer, parameter :: CELSIUS = 1

    ! elev_units
    integer, parameter :: FEET = 0
    integer, parameter :: METERS = 1

    ! runoff_units
    integer, parameter :: CFS = 0
    integer, parameter :: CMS = 1

    ! cov_type
    integer, parameter :: BARESOIL = 0
    integer, parameter :: GRASSES = 1
    integer, parameter :: SHRUBS = 2
    integer, parameter :: TREES = 3
    integer, parameter :: CONIFEROUS = 4

    ! soil_type
    integer, parameter :: SAND = 1
    integer, parameter :: LOAM = 2
    integer, parameter :: CLAY = 3

    ! Hemisphere
    integer, parameter :: NORTHERN = 0
    integer, parameter :: SOUTHERN = 1

    ! hru_type
    integer, parameter :: INACTIVE = 0
    integer, parameter :: LAND = 1
    integer, parameter :: LAKE = 2
    integer, parameter :: SWALE = 3
    integer, parameter :: GLACIER = 4

    ! seg_type
    integer, parameter :: OUTFLOW_SEGMENT = 0

    ! cascade flags
    integer, parameter :: CASCADE_OFF = 0
    integer, parameter :: CASCADE_NORMAL = 1
    integer, parameter :: CASCADE_HRU_SEGMENT = 2
    integer, parameter :: CASCADEGW_OFF = 0
    integer, parameter :: CASCADEGW_SAME = 2

    ! model_mode
    integer, parameter :: GSFLOW = 0
    integer, parameter :: PRMS = 1
    integer, parameter :: MODFLOW = 2
    integer, parameter :: DOCUMENTATION = 99
    integer, parameter :: RUN = 0
    integer, parameter :: DECL = 1
    integer, parameter :: INIT = 2
    integer, parameter :: CLEAN = 3
    integer, parameter :: SETDIMENS = 4
    integer, parameter :: WRITE_CLIMATE=24, CONVERT=25, CLIMATE=26, POTET=27, TRANSPIRE=28, FROST=29

    ! Error Codes
    integer, parameter :: ERROR_modflow = -5
    integer, parameter :: ERROR_read = -4
    integer, parameter :: ERROR_open_out = -3
    integer, parameter :: ERROR_open_in = -2
    integer, parameter :: ERROR_write = -1
    integer, parameter :: ERROR_control = 1
    integer, parameter :: ERROR_var = 2
    integer, parameter :: ERROR_dim = 3
    integer, parameter :: ERROR_param = 4
    integer, parameter :: ERROR_data = 5
    integer, parameter :: ERROR_time = 6
    integer, parameter :: ERROR_temp = 7
    integer, parameter :: ERROR_streamflow = 8
    integer, parameter :: ERROR_basin = 9
    integer, parameter :: ERROR_cbh = 10
    integer, parameter :: ERROR_cascades = 11
    integer, parameter :: ERROR_restart = 12
    integer, parameter :: ERROR_dynamic = 13
    integer, parameter :: ERROR_water_use = 14
    integer, parameter :: ERROR_decl_get = 15
    integer, parameter :: ERROR_module = 16
    integer, parameter :: ERROR_lake = 17
    integer, parameter :: ERROR_soilzone = 18

      ! debug print flag:
      ! -2=DEBUG_minimum
      ! -1=DEBUG_less
      ! 0=DEBUG_normal; 1=DEBUG_WB (water balances); 2=basin;
      ! 4=basin_sum; 5=DEBUG_SOLTAB; 7=soil zone;
      ! 9=snowcomp; 13=cascade; 14=subbasin tree
    integer, parameter :: DEBUG_minimum = -2
    integer, parameter :: DEBUG_less = -1
    integer, parameter :: DEBUG_normal = 0
    integer, parameter :: DEBUG_WB = 1
    integer, parameter :: DEBUG_SOLTAB = 5

! module flags
    integer, parameter :: precip_1sta_module = 1, precip_laps_module = 2, precip_dist2_module = 3
    integer, parameter :: ide_dist_module = 5, xyz_dist_module = 6, precip_map_module = 9, climate_hru_module = 7
    integer, parameter :: temp_1sta_module = 1, temp_laps_module = 2, temp_dist2_module = 3
    integer, parameter :: temp_map_module = 9, temp_sta_module = 8
    integer, parameter :: potet_jh_module = 1, potet_hamon_module = 2, potet_pan_module = 4
    integer, parameter :: potet_pt_module = 5, potet_pm_sta_module = 6, potet_hs_module = 10, potet_pm_module = 11
    integer, parameter :: strmflow_muskingum_module = 4, strmflow_muskingum_mann_module = 7
    integer, parameter :: strmflow_muskingum_lake_module = 3, strmflow_in_out_module = 5, strmflow_noroute_module = 1
    integer, parameter :: smidx_module = 1, carea_module = 2, ddsolrad_module = 1, ccsolrad_module = 2

    integer, parameter :: OFF = 0
    integer, parameter :: ACTIVE = 1
    integer, parameter :: INT_TYPE = 1
    integer, parameter :: REAL_TYPE = 2
    integer, parameter :: DBLE_TYPE = 3
    integer, parameter :: SAVE_INIT = 0
    integer, parameter :: READ_INIT = 1

END MODULE PRMS_CONSTANTS
