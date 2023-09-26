!***********************************************************************
! Initiates development of a snowpack and simulates snow accumulation
! and depletion processes using an energy-budget approach
!
! Modified glacier melt and glacier basal melt
! These modifications includes albedo info for saving between runs 2/00
!***********************************************************************

! PRMS_SNOW module for defining stateful variables

      MODULE PRMS_SNOW
      IMPLICIT NONE
      !****************************************************************
      !   Local Constants
      INTRINSIC :: ACOS
      INTEGER, PARAMETER :: MAXALB = 15
      REAL, PARAMETER :: PI = ACOS(-1.0)
      INTEGER, PARAMETER :: not_a_glacier_hru = -1

      !****************************************************************
      !   Local Variables
      character(len=*), parameter :: MODDESC = 'Snow Dynamics'
      character(len=8), parameter :: MODNAME = 'snowcomp'
      character(len=*), parameter :: Version_snowcomp = '2021-09-15'
      INTEGER, SAVE :: Active_glacier
      INTEGER, SAVE, ALLOCATABLE :: Int_alb(:)
      REAL, SAVE :: Acum(MAXALB), Amlt(MAXALB)
      REAL, SAVE, ALLOCATABLE :: Snowcov_areasv(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Scrv(:), Pss(:), Pksv(:), Pst(:)
      REAL, SAVE, ALLOCATABLE :: Salb(:), Slst(:)

      !****************************************************************
      !   Declared Variables

      INTEGER :: Yrdays5
      INTEGER, SAVE, ALLOCATABLE :: Pptmix_nopack(:), Lst(:)
      INTEGER, SAVE, ALLOCATABLE :: Iasw(:), Iso(:), Mso(:), Lso(:)
      DOUBLE PRECISION, SAVE :: Basin_snowmelt, Basin_pweqv, Basin_tcal
      DOUBLE PRECISION, SAVE :: Basin_snowcov, Basin_snowevap
      DOUBLE PRECISION, SAVE :: Basin_snowdepth, Basin_pk_precip
      REAL, SAVE, ALLOCATABLE :: Snowmelt(:), Snow_evap(:)
      REAL, SAVE, ALLOCATABLE :: Albedo(:), Pk_temp(:), Pk_den(:)
      REAL, SAVE, ALLOCATABLE :: Pk_def(:), Pk_ice(:), Freeh2o(:)
      REAL, SAVE, ALLOCATABLE :: Snowcov_area(:), Tcal(:)
      REAL, SAVE, ALLOCATABLE :: Snsv(:), Pk_precip(:)
      REAL, SAVE, ALLOCATABLE :: Frac_swe(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Pk_depth(:), Pkwater_ante(:), Ai(:)
      DOUBLE PRECISION, SAVE :: Basin_glacrevap, Basin_snowicecov, Basin_glacrb_melt
      REAL, SAVE, ALLOCATABLE :: Glacrmelt(:), Glacr_evap(:), Glacr_albedo(:), Glacr_pk_den(:)
      REAL, SAVE, ALLOCATABLE :: Glacr_pk_ice(:), Glacr_freeh2o(:), Glacrcov_area(:)
      REAL, SAVE, ALLOCATABLE :: Glacrb_melt(:), Glacr_pk_def(:), Glacr_pk_temp(:), Ann_tempc(:)
      REAL, SAVE, ALLOCATABLE :: Glacr_air_5avtemp1(:), Glacr_air_deltemp(:), Glacr_air_5avtemp(:)
      REAL, SAVE, ALLOCATABLE :: Glacr_5avsnow1(:), Glacr_5avsnow(:), Glacr_delsnow(:), Glacr_freeh2o_capm(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Glacr_pkwater_ante(:), Glacr_pkwater_equiv(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Glacr_pk_depth(:), Glacr_pss(:), Glacr_pst(:)
      !****************************************************************
      !   Declared Parameters

      INTEGER, SAVE, ALLOCATABLE :: Melt_look(:), Melt_force(:), Tstorm_mo(:, :)
      INTEGER, SAVE, ALLOCATABLE :: Hru_deplcrv(:)
      REAL, SAVE :: Albset_rnm, Albset_rna, Albset_snm, Albset_sna
      REAL, SAVE, ALLOCATABLE :: Emis_noppt(:), Freeh2o_cap(:), Cecn_coef(:, :)
      REAL, SAVE, ALLOCATABLE :: Den_init(:), Settle_const(:), Den_max(:)
      REAL, SAVE, ALLOCATABLE :: Rad_trncf(:), Snarea_thresh(:), Snowpack_init(:)
      REAL, SAVE, ALLOCATABLE :: Snarea_curve(:, :)
      REAL, SAVE, ALLOCATABLE :: Snarea_a(:), Snarea_b(:), Snarea_c(:), Snarea_d(:)
      REAL, SAVE, ALLOCATABLE :: Glacr_layer(:), Albedo_coef(:), Albedo_ice(:)
      REAL, SAVE, ALLOCATABLE :: Glacr_freeh2o_cap(:), Glacier_frac_init(:), Glrette_frac_init(:)

      END MODULE PRMS_SNOW

!***********************************************************************
!     Main snowcomp routine
!***********************************************************************
      INTEGER FUNCTION snowcomp()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: snodecl, snoinit, snorun
      EXTERNAL :: snowcomp_restart
!***********************************************************************
      snowcomp = 0

      IF ( Process_flag==RUN ) THEN
        snowcomp = snorun()
      ELSEIF ( Process_flag==DECL ) THEN
        snowcomp = snodecl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL snowcomp_restart(READ_INIT)
        snowcomp = snoinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL snowcomp_restart(SAVE_INIT)
      ENDIF

      END FUNCTION snowcomp

!***********************************************************************
!     snodecl - set up parameters for snowmelt computations
!   Declared Parameters
!     den_init, settle_const, den_max, melt_look
!     melt_force, rad_trncf, hru_deplcrv, snarea_curve, snarea_thresh
!     albset_rnm, albset_rna, albset_snm, albset_sna, potet_sublim
!     emis_noppt, cecn_coef, freeh2o_cap, tstorm_mo, tmax_allsnow
!     hru_area, cov_type, covden_win
!     glacr_freeh2o_cap, glacr_layer
!***********************************************************************
      INTEGER FUNCTION snodecl()
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ACTIVE, OFF, MONTHS_PER_YEAR
      USE PRMS_MODULE, ONLY: Model, Nhru, Ndepl, &
     &    Init_vars_from_file, Glacier_flag, Snarea_curve_flag
      USE PRMS_SNOW
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      snodecl = 0

      CALL print_module(MODDESC, MODNAME, Version_snowcomp)

! declare variables
      ALLOCATE ( Scrv(Nhru) )
      IF ( declvar(MODNAME, 'scrv', 'nhru', Nhru, 'double', &
     &     'Snowpack water equivalent plus a portion of new snow on each HRU', &
     &     'inches', Scrv)/=0 ) CALL read_error(3, 'scrv')

      ALLOCATE ( Pksv(Nhru) )
      IF ( declvar(MODNAME, 'pksv', 'nhru', Nhru, 'double', &
     &     'Snowpack water equivalent when there is new snow and in melt phase;'// &
     &     ' used to interpolate between depletion curve and 100 percent on each HRU', &
     &     'inches', Pksv)/=0 ) CALL read_error(3, 'pksv')

      ALLOCATE ( Snowcov_areasv(Nhru) )
      IF ( declvar(MODNAME, 'snowcov_areasv', 'nhru', Nhru, 'real', &
     &     'Snow cover fraction when there is new snow and in melt phase;'// &
     &     ' used to interpolate between depletion curve and 100 percent on each HRU', &
     &     'decimal fraction', Snowcov_areasv)/=0 ) CALL read_error(3, 'snowcov_areasv')

      ALLOCATE ( Salb(Nhru) )
      IF ( declvar(MODNAME, 'salb', 'nhru', Nhru, 'real', &
     &     'Days since last new snow to reset albedo for each HRU', &
     &     'days', Salb)/=0 ) CALL read_error(3, 'salb')

      ALLOCATE ( Slst(Nhru) )
      IF ( declvar(MODNAME, 'slst', 'nhru', Nhru, 'real', &
     &     'Days since last new snow for each HRU', &
     &     'days', Slst)/=0 ) CALL read_error(3, 'slst')

      ALLOCATE ( Int_alb(Nhru) )
      IF ( declvar(MODNAME, 'int_alb', 'nhru', Nhru, 'integer', &
     &     'Flag to indicate (1: accumulation season curve; 2: use of the melt season curve)', &
     &     'none', Int_alb)/=0 ) CALL read_error(3, 'int_alb')

! Glacier declares
      IF ( Glacier_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'yrdays5', 'one', 1, 'integer', &
     &     'Number of days since last 5-year mark', &
     &     'none', Yrdays5)/=0 ) CALL read_error(3, 'yrdays5')

        ALLOCATE ( Glacr_freeh2o_capm(Nhru) )
        IF ( declvar(MODNAME, 'glacr_freeh2o_capm', 'nhru', Nhru, 'real', &
     &       'Free-water holding capacity of glacier or glacierette ice, changes to 0 if active layer melts', &
     &       'decimal fraction', Glacr_freeh2o_capm)/=0 ) CALL read_error(3, 'glacr_freeh2o_capm')

        ALLOCATE ( Glacrb_melt(Nhru) )
        IF ( declvar(MODNAME, 'glacrb_melt', 'nhru', Nhru, 'real', &
             'Glacier or glacierette basal melt, goes to soil', &
             'inches/day', Glacrb_melt)/=0 ) CALL read_error(3, 'glacrb_melt')

        ALLOCATE ( Ann_tempc(Nhru) )
        IF ( declvar(MODNAME, 'ann_tempc', 'nhru', Nhru, 'real', &
     &       'Current average year air temperature over each HRU', &
     &       'degrees Celsius', Ann_tempc)/=0 ) CALL read_error(3, 'ann_tempc')

       ALLOCATE ( Glacr_air_5avtemp(Nhru) )
        IF ( declvar(MODNAME, 'glacr_air_5avtemp', 'nhru', Nhru, 'real', &
     &       'Current 5-yr average summer (June July Aug) air temperature over glacier or glacierette HRUs', &
     &       'degrees Celsius', Glacr_air_5avtemp)/=0 ) CALL read_error(3, 'glacr_air_5avtemp')

        ALLOCATE ( Glacr_air_5avtemp1(Nhru) )
        IF ( declvar(MODNAME, 'glacr_air_5avtemp1', 'nhru', Nhru, 'real', &
     &       'First 5-yr average summer temperature over glacier or glacierette HRUs', &
     &       'degrees Celsius', Glacr_air_5avtemp1)/=0 ) CALL read_error(3, 'glacr_air_5avtemp1')

        ALLOCATE ( Glacr_air_deltemp(Nhru) )
        IF ( declvar(MODNAME, 'glacr_air_deltemp', 'nhru', Nhru, 'real', &
     &       'Change in 5-yr average air temperature over glacier or glacierette HRUs from first', &
     &       'degrees Celsius', Glacr_air_deltemp)/=0 ) CALL read_error(3, 'glacr_air_deltemp')

       ALLOCATE ( Glacr_5avsnow(Nhru) )
        IF ( declvar(MODNAME, 'glacr_5avsnow', 'nhru', Nhru, 'real', &
     &       'Current 5-yr average snow over glacier or glacierette HRUs', &
     &       'inches/yr', Glacr_5avsnow)/=0 ) CALL read_error(3, 'glacr_5avsnow')

        ALLOCATE ( Glacr_5avsnow1(Nhru) )
        IF ( declvar(MODNAME, 'glacr_5avsnow1', 'nhru', Nhru, 'real', &
     &       'First 5-yr average snow over glacier or glacierette HRUs', &
     &       'inches/yr', Glacr_5avsnow1)/=0 ) CALL read_error(3, 'glacr_5avsnow1')

        ALLOCATE ( Glacr_delsnow(Nhru) )
        IF ( declvar(MODNAME, 'glacr_delsnow', 'nhru', Nhru, 'real', &
     &       'Change in 5-yr average snow over glacier or glacierette HRUs from first', &
     &       'inches/yr', Glacr_delsnow)/=0 ) CALL read_error(3, 'glacr_delsnow')

        ALLOCATE ( Glacr_pk_temp(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pk_temp', 'nhru', Nhru, 'real', &
     &       'Temperature of the glacier or glacierette on each HRU', &
     &       'degrees Celsius', Glacr_pk_temp)/=0 ) CALL read_error(3, 'glacr_pk_temp')

        ALLOCATE ( Glacr_pk_def(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pk_def', 'nhru', Nhru, 'real', &
     &       'Heat deficit, amount of heat necessary to make the glacier'// &
     &       ' or or glacierette snowpack isothermal at 0 degrees Celsius', &
     &       'Langleys', Glacr_pk_def)/=0 ) CALL read_error(3, 'glacr_pk_def')

        ALLOCATE ( Glacr_pk_den(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pk_den', 'nhru', Nhru, 'real', &
     &       'Density of the icepack on each glacier or glacierette HRU, hard coded to equal 0.917', &
     &       'gm/cm3', Glacr_pk_den)/=0 ) CALL read_error(3, 'glacr_pk_den')

        ALLOCATE ( Glacr_albedo(Nhru) )
        IF ( declvar(MODNAME, 'glacr_albedo', 'nhru', Nhru, 'real', &
     &       'Ice surface albedo or the fraction of radiation reflected from'// &
     &       ' the icepack surface for each glacier or glacierette HRU', &
     &       'decimal fraction', Glacr_albedo)/=0 ) CALL read_error(3, 'glacr_albedo')

        ALLOCATE ( Glacr_evap(Nhru) )
        IF ( declvar(MODNAME, 'glacr_evap', 'nhru', Nhru, 'real', &
     &       'Evaporation and sublimation from icepack on each glacier or glacierette HRU', &
     &       'inches', Glacr_evap)/=0 ) CALL read_error(3, 'glacr_evap')

        ALLOCATE ( Glacrmelt(Nhru) )
        IF ( declvar(MODNAME, 'glacrmelt', 'nhru', Nhru, 'real', &
     &       'Melt from icepack on each glacier or glacierette HRU, includes rain water that does not absorb', &
     &       'inches', Glacrmelt)/=0 ) CALL read_error(3, 'glacrmelt')

        ALLOCATE ( Glacr_pkwater_equiv(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pkwater_equiv', 'nhru', Nhru, 'double', &
     &       'Icepack water equivalent on each glacier or glacierette HRU', &
     &       'inches', Glacr_pkwater_equiv)/=0 ) CALL read_error(3, 'glacr_pkwater_equiv')

        ALLOCATE ( Glacr_pkwater_ante(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pkwater_ante', 'nhru', Nhru, 'double', &
     &       'Antecedent icepack water equivalent on each glacier or glacierette HRU', &
     &       'inches', Glacr_pkwater_ante)/=0 ) CALL read_error(3, 'glacr_pkwater_ante')

        ALLOCATE ( Glacrcov_area(Nhru) )
        IF ( declvar(MODNAME, 'glacrcov_area', 'nhru', Nhru, 'real', &
     &       'Ice-covered area (no snowpack) on each glacier HRU or HRU with glacierette at start of step', &
     &       'decimal fraction', Glacrcov_area)/=0 ) CALL read_error(3, 'glacrcov_area')

        ALLOCATE ( Glacr_pk_ice(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pk_ice', 'nhru', Nhru, 'real', &
     &       'Storage of frozen water in the icepack on each glacier or glacierette HRU', &
     &       'inches', Glacr_pk_ice)/=0 ) CALL read_error(3, 'glacr_pk_ice')

        ALLOCATE ( Glacr_freeh2o(Nhru) )
        IF ( declvar(MODNAME, 'glacr_freeh2o', 'nhru', Nhru, 'real', &
     &       'Storage of free liquid water in the icepack on each glacier or glacierette HRU', &
     &       'inches', Glacr_freeh2o)/=0 ) CALL read_error(3, 'glacr_freeh2o')

        ALLOCATE ( Glacr_pk_depth(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pk_depth', 'nhru', Nhru, 'double', &
     &       'Depth of icepack on each glacier or glacierette HRU, make essentially infinite', &
     &       'inches', Glacr_pk_depth)/=0 ) CALL read_error(3, 'glacr_pk_depth')

        ALLOCATE ( Glacr_pss(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pss', 'nhru', Nhru, 'double', &
     &       'Previous glacier or glacierette pack water equivalent plus new ice', &
     &       'inches', Glacr_pss)/=0 ) CALL read_error(3, 'glacr_pss')

        ALLOCATE ( Glacr_pst(Nhru) )
        IF ( declvar(MODNAME, 'glacr_pst', 'nhru', Nhru, 'double', &
     &       'While an icepack exists, glacr_pst tracks the maximum ice water equivalent of that icepack', &
     &       'inches', Glacr_pst)/=0 ) CALL read_error(3, 'glacr_pst')

        IF ( declvar(MODNAME, 'basin_snowicecov', 'one', 1, 'double', &
     &       'Basin area-weighted average snow and glacier and glacierette covered area', &
     &       'decimal fraction', Basin_snowicecov)/=0 ) CALL read_error(3, 'basin_snowicecov')

        IF ( declvar(MODNAME, 'basin_glacrb_melt', 'one', 1, 'double', &
     &       'Basin area-weighted average basal melt of glacier, goes to soil', &
     &       'inches', Basin_glacrb_melt)/=0 ) CALL read_error(3, 'basin_glacrb_melt')

        IF ( declvar(MODNAME, 'basin_glacrevap', 'one', 1, 'double', &
     &       'Basin area-weighted average glacier ice evaporation and sublimation', &
     &       'inches', Basin_glacrevap)/=0 ) CALL read_error(3, 'basin_glacrevap')
      ENDIF

      IF ( declvar(MODNAME, 'basin_snowdepth', 'one', 1, 'double', &
     &     'Basin area-weighted average snow depth', &
     &     'inches', Basin_snowdepth)/=0 ) CALL read_error(3, 'basin_snowdepth')

      ALLOCATE ( Pk_precip(Nhru) )
      IF ( declvar(MODNAME, 'pk_precip', 'nhru', Nhru, 'real', &
     &     'Precipitation added to snowpack for each HRU', &
     &     'inches', Pk_precip)/=0 ) CALL read_error(3, 'pk_precip')

      IF ( declvar(MODNAME, 'basin_pk_precip', 'one', 1, 'double', &
     &     'Basin area-weighted average precipitation added to snowpack', &
     &     'inches', Basin_pk_precip)/=0 ) CALL read_error(3, 'basin_pk_precip')

      ALLOCATE ( Albedo(Nhru) )
      IF ( declvar(MODNAME, 'albedo', 'nhru', Nhru, 'real', &
     &     'Snow surface albedo or the fraction of radiation reflected from the snowpack surface for each HRU', &
     &     'decimal fraction', Albedo)/=0 ) CALL read_error(3, 'albedo')

      ALLOCATE ( Pk_temp(Nhru) )
      IF ( declvar(MODNAME, 'pk_temp', 'nhru', Nhru, 'real', &
     &     'Temperature of the snowpack on each HRU', &
     &     'degrees Celsius', Pk_temp)/=0 ) CALL read_error(3, 'pk_temp')

      ALLOCATE ( Pk_den(Nhru) )
      IF ( declvar(MODNAME, 'pk_den', 'nhru', Nhru, 'real', &
     &     'Density of the snowpack on each HRU', &
     &     'gm/cm3', Pk_den)/=0 ) CALL read_error(3, 'pk_den')

      IF ( declvar(MODNAME, 'basin_tcal', 'one', 1, 'double', &
     &     'Basin area-weighted average net snowpack energy balance', &
     &     'Langleys', Basin_tcal)/=0 ) CALL read_error(3, 'basin_tcal')

      ALLOCATE ( Tcal(Nhru) )
      IF ( declvar(MODNAME, 'tcal', 'nhru', Nhru, 'real', &
     &     'Net snowpack energy balance on each HRU', &
     &     'Langleys', Tcal)/=0 ) CALL read_error(3, 'tcal')

      ALLOCATE ( Snow_evap(Nhru) )
      IF ( declvar(MODNAME, 'snow_evap', 'nhru', Nhru, 'real', &
     &     'Evaporation and sublimation from snowpack on each HRU', &
     &     'inches', Snow_evap)/=0 ) CALL read_error(3, 'snow_evap')

      ALLOCATE ( Snowmelt(Nhru) )
      IF ( declvar(MODNAME, 'snowmelt', 'nhru', Nhru, 'real', &
     &     'Snowmelt from snowpack on each HRU (not including snow on glacier)', &
     &     'inches', Snowmelt)/=0 ) CALL read_error(3, 'snowmelt')

      IF ( declvar(MODNAME, 'basin_snowmelt', 'one', 1, 'double', &
     &     'Basin area-weighted average snowmelt (not on including snow on glacier)', &
     &     'inches', Basin_snowmelt)/=0 ) CALL read_error(3, 'basin_snowmelt')

      IF ( declvar(MODNAME, 'basin_pweqv', 'one', 1, 'double', &
     &     'Basin area-weighted average snowpack water equivalent (not including glacier)', &
     &     'inches', Basin_pweqv)/=0 ) CALL read_error(3, 'basin_pweqv')

      ALLOCATE ( Pkwater_ante(Nhru) )
      IF ( declvar(MODNAME, 'pkwater_ante', 'nhru', Nhru, 'double', &
     &     'Antecedent snowpack water equivalent on each HRU', &
     &     'inches', Pkwater_ante)/=0 ) CALL read_error(3, 'pkwater_ante')

      ALLOCATE ( Snowcov_area(Nhru) )
      IF ( declvar(MODNAME, 'snowcov_area', 'nhru', Nhru, 'real', &
     &     'Snow-covered area on each HRU prior to melt and sublimation unless snowpack depleted', &
     &     'decimal fraction', Snowcov_area)/=0 ) CALL read_error(3, 'snowcov_area')

      IF ( declvar(MODNAME, 'basin_snowevap', 'one', 1, 'double', &
     &     'Basin area-weighted average evaporation and sublimation from snowpack (not including glacier)', &
     &     'inches', Basin_snowevap)/=0 ) CALL read_error(3, 'basin_snowevap')

      IF ( declvar(MODNAME, 'basin_snowcov', 'one', 1, 'double', &
     &     'Basin area-weighted average snow-covered area', &
     &     'decimal fraction', Basin_snowcov)/=0 ) CALL read_error(3, 'basin_snowcov')

      !rpayn commented
      ALLOCATE ( Pptmix_nopack(Nhru) )
      IF ( declvar(MODNAME, 'pptmix_nopack', 'nhru', Nhru, 'integer', &
     &     'Flag indicating that a mixed precipitation event has'// &
     &     ' occurred with no snowpack present on an HRU (1), otherwise (0)', &
     &     'none', Pptmix_nopack)/=0 ) CALL read_error(3, 'pptmix_nopack')

      !rpayn commented
      ALLOCATE ( Iasw(Nhru) )
      IF ( declvar(MODNAME, 'iasw', 'nhru', Nhru, 'integer', &
     &     'Flag indicating that snow covered area is'// &
     &     ' interpolated between previous location on curve and maximum (1), or is on the defined curve (0)', &
     &     'none', Iasw)/=0 ) CALL read_error(3, 'iasw')

      !rpayn commented
      ALLOCATE ( Iso(Nhru) )
      IF ( declvar(MODNAME, 'iso', 'nhru', Nhru, 'integer', &
     &     'Flag to indicate if time is before (1) or after (2) the day to force melt season (melt_force)', &
     &     'none', Iso)/=0 ) CALL read_error(3, 'iso')

      !rpayn commented
      ALLOCATE ( Mso(Nhru) )
      IF ( declvar(MODNAME, 'mso', 'nhru', Nhru, 'integer', &
     &     'Flag to indicate if time is before (1) or after (2) the first potential day for melt season (melt_look)', &
     &     'none', Mso)/=0 ) CALL read_error(3, 'mso')

      !rpayn commented
      ALLOCATE ( Lso(Nhru) )
      IF ( declvar(MODNAME, 'lso', 'nhru', Nhru, 'integer', &
     &     'Counter for tracking the number of days the snowpack is at or above 0 degrees Celsius', &
     &     'number of iterations', Lso)/=0 ) CALL read_error(3, 'lso')

      !rpayn commented
      ALLOCATE ( Lst(Nhru) )
      IF ( declvar(MODNAME, 'lst', 'nhru', Nhru, 'integer', &
     &     'Flag indicating whether there was new snow that'// &
     &     ' was insufficient to reset the albedo curve (1) (albset_snm or albset_sna), otherwise (0)', &
     &     'none', Lst)/=0 ) CALL read_error(3, 'lst')

      !rpayn commented
      ALLOCATE ( Pk_def(Nhru) )
      IF ( declvar(MODNAME, 'pk_def', 'nhru', Nhru, 'real', &
     &     'Heat deficit, amount of heat necessary to make the snowpack isothermal at 0 degrees Celsius', &
     &     'Langleys', Pk_def)/=0 ) CALL read_error(3, 'pk_def')

      !rpayn commented
      ALLOCATE ( Pk_ice(Nhru) )
      IF ( declvar(MODNAME, 'pk_ice', 'nhru', Nhru, 'real', &
     &     'Storage of frozen water in the snowpack on each HRU', &
     &     'inches', Pk_ice)/=0 ) CALL read_error(3, 'pk_ice')

      !rpayn commented
      ALLOCATE ( Freeh2o(Nhru) )
      IF ( declvar(MODNAME, 'freeh2o', 'nhru', Nhru, 'real', &
     &     'Storage of free liquid water in the snowpack on each HRU', &
     &     'inches', Freeh2o)/=0 ) CALL read_error(3, 'freeh2o')

      !rpayn commented
      ALLOCATE ( Pk_depth(Nhru) )
      IF ( declvar(MODNAME, 'pk_depth', 'nhru', Nhru, 'double', &
     &     'Depth of snowpack on each HRU', &
     &     'inches', Pk_depth)/=0 ) CALL read_error(3, 'pk_depth')

      !rpayn commented
      ALLOCATE ( Pss(Nhru) )
      IF ( declvar(MODNAME, 'pss', 'nhru', Nhru, 'double', &
     &     'Previous snowpack water equivalent plus new snow', &
     &     'inches', Pss)/=0 ) CALL read_error(3, 'pss')

      !rpayn commented
      ALLOCATE ( Pst(Nhru) )
      IF ( declvar(MODNAME, 'pst', 'nhru', Nhru, 'double', &
     &     'While a snowpack exists, pst tracks the maximum snow water equivalent of that snowpack', &
     &     'inches', Pst)/=0 ) CALL read_error(3, 'pst')

      !rpayn commented
      ALLOCATE ( Snsv(Nhru) )
      IF ( declvar(MODNAME, 'snsv', 'nhru', Nhru, 'real', &
     &     'Tracks the cumulative amount of new snow until'// &
     &     ' there is enough to reset the albedo curve (albset_snm or albset_sna)', &
     &     'inches', Snsv)/=0 ) CALL read_error(3, 'snsv')

      ALLOCATE ( Ai(Nhru) )
      IF ( declvar(MODNAME, 'ai', 'nhru', Nhru, 'double', &
     &     'Maximum snowpack for each HRU', &
     &     'inches', Ai)/=0 ) CALL read_error(3, 'ai')

      ALLOCATE ( Frac_swe(Nhru) )
      IF ( declvar(MODNAME, 'frac_swe', 'nhru', Nhru, 'real', &
     &     'Fraction of maximum snow-water equivalent (snarea_thresh) on each HRU', &
     &     'decimal fraction', Frac_swe)/=0 ) CALL read_error(3, 'frac_swe')

! declare parameters
      IF ( Glacier_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Albedo_coef(Nhru) )
        IF ( declparam(MODNAME, 'albedo_coef', 'nhru', 'real', &
     &       '0.137', '0.1', '0.3', &
     &       'Coefficient in calculation of ice albedo', &
     &       'Coefficient in calculation of ice albedo', &
     &       'none')/=0 ) CALL read_error(1, 'albedo_coef')

        ALLOCATE ( Albedo_ice(Nhru) )
        IF ( declparam(MODNAME, 'albedo_ice', 'nhru', 'real', &
     &       '0.344', '0.2', '0.6', &
     &       'Ice albedo 300 meters below ELA', &
     &       'Ice albedo 300 meters below ELA', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'albedo_ice')

        ALLOCATE ( Glacr_freeh2o_cap(Nhru) )
        IF ( declparam(MODNAME, 'glacr_freeh2o_cap', 'nhru', 'real', &
     &       '0.002', '0.0', '0.01', &
     &       'Free-water holding capacity of glacier ice', &
     &       'Free-water holding capacity of glacier ice expressed as a' // &
     &       ' decimal fraction of the frozen water content of the glacier ice (glacr_pk_ice)',  &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'glacr_freeh2o_cap')

        ALLOCATE ( Glacr_layer(Nhru) )
        IF ( declparam(MODNAME, 'glacr_layer', 'nhru', 'real', &
     &       '3.94', '0.0', '590.6', &
     &       'Active layer on glacier', &
     &       'Active layer is 0 to 15 m (590.6 inches) thick at start of year, when' // &
     &       ' melts will set daily glacr_pk_temp to 0',  &
     &       'inches')/=0 ) CALL read_error(1, 'glacr_layer')

        IF ( Init_vars_from_file==OFF ) THEN
          ALLOCATE ( Glacier_frac_init(Nhru) )
          IF ( declparam(MODNAME, 'glacier_frac_init', 'nhru', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Initial fraction of glaciation (0=none; 1=100%) in glacier-capable HRU', &
     &       'Initial fraction of glaciation (0=none; 1=100%) in glacier-capable HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'glacier_frac_init')

          ALLOCATE ( Glrette_frac_init(Nhru) )
          IF ( declparam(MODNAME, 'glrette_frac_init', 'nhru', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Initial fraction of glacierette (too small for glacier dynamics)', &
     &       'Initial fraction of glacierette  (too small for glacier dynamics)', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'glrette_frac_init')
        ENDIF
      ENDIF

      ALLOCATE ( Den_init(Nhru) )
      IF ( declparam(MODNAME, 'den_init', 'nhru', 'real', &
     &     '0.10', '0.01', '0.5', &
     &     'Initial density of new-fallen snow', &
     &     'Initial density of new-fallen snow', &
     &     'gm/cm3')/=0 ) CALL read_error(1, 'den_init')

      ALLOCATE ( Settle_const(Nhru) )
      IF ( declparam(MODNAME, 'settle_const', 'nhru', 'real', &
     &     '0.10', '0.01', '0.5', &
     &     'Snowpack settlement time constant', &
     &     'Snowpack settlement time constant', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'settle_const')

      ALLOCATE ( Den_max(Nhru) )
      IF ( declparam(MODNAME, 'den_max', 'nhru', 'real', &
     &     '0.6', '0.1', '0.8', &
     &     'Average maximum snowpack density', &
     &     'Average maximum snowpack density', &
     &     'gm/cm3')/=0 ) CALL read_error(1, 'den_max')

      ALLOCATE ( Melt_look(Nhru) )
      IF ( declparam(MODNAME, 'melt_look', 'nhru', 'integer', &
     &     '90', '1', '366', &
     &     'Julian date to start looking for spring snowmelt for each HRU', &
     &     'Julian date to start looking for spring snowmelt stage for each HRU;'// &
     &     ' varies with region depending on length of time that'// &
     &     ' permanent snowpack exists', &
     &     'Julian day')/=0 ) CALL read_error(1, 'melt_look')

      ALLOCATE ( Melt_force(Nhru) )
      IF ( declparam(MODNAME, 'melt_force', 'nhru', 'integer', &
     &     '140', '1', '366', &
     &     'Julian date to force snowpack to spring snowmelt stage for each HRU', &
     &     'Julian date to force snowpack to spring snowmelt stage for each HRU;'// &
     &     ' varies with region depending on length of time that'// &
     &     ' permanent snowpack exists', &
     &     'Julian day')/=0 ) CALL read_error(1, 'melt_force')

      ALLOCATE ( Rad_trncf(Nhru) )
      IF ( declparam(MODNAME, 'rad_trncf', 'nhru', 'real', &
     &     '0.5', '0.0', '1.0', &
     &     'Solar radiation transmission coefficient', &
     &     'Transmission coefficient for short-wave radiation through'// &
     &     ' the winter vegetation canopy', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'rad_trncf')

      ALLOCATE ( Hru_deplcrv(Nhru) )
      ALLOCATE ( Snarea_curve(11, Ndepl) )
      IF ( Snarea_curve_flag==0 .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(MODNAME, 'hru_deplcrv', 'nhru', 'integer', &
     &       '1', 'bounded', 'ndepl', &
     &       'Index number for snowpack areal depletion curve', &
     &       'Index number for the snowpack areal depletion curve associated with each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'hru_deplcrv')
        IF ( declparam(MODNAME, 'snarea_curve', 'ndeplval', 'real', &
     &       '1.0', '0.0', '1.0', &
     &       'Snow area depletion curve values', &
     &       'Snow area depletion curve values, 11 values for each'// &
     &       ' curve (0.0 to 1.0 in 0.1 increments)', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'snarea_curve')
      ENDIF
      IF ( Snarea_curve_flag==1 .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Snarea_a(Nhru) )
        IF ( declparam(MODNAME, 'snarea_a', 'nhru', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Snow area depletion curve minimum SWE value', &
     &       'Snow area depletion curve minimum SWE value for each HRU', &
     &       'inches')/=0 ) CALL read_error(1, 'snarea_a')
        ALLOCATE ( Snarea_b(Nhru) )
        IF ( declparam(MODNAME, 'snarea_b', 'nhru', 'real', &
     &       '2.0', '0.5', '20.0', &
     &       'Snow area depletion S-curve coefficient B', &
     &       'Snow area depletion S-curve coefficient B for each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'snarea_b')
        ALLOCATE ( Snarea_c(Nhru) )
        IF ( declparam(MODNAME, 'snarea_c', 'nhru', 'real', &
     &       '1.5', '0.001', '3.0', &
     &       'Snow area depletion S-curve coefficient C', &
     &       'Snow area depletion S-curve coefficient C for each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'snarea_c')
        ALLOCATE ( Snarea_d(Nhru) )
        IF ( declparam(MODNAME, 'snarea_d', 'nhru', 'real', &
     &       '0.975', '0.0', '3.0', &
     &       'Snow area depletion S-curve coefficient D', &
     &       'Snow area depletion S-curve coefficient D for each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'snarea_d')
      ENDIF

      ALLOCATE ( Snarea_thresh(Nhru) )
      IF ( declparam(MODNAME, 'snarea_thresh', 'nhru', 'real', &
     &     '50.0', '0.0', '200.0', &
     &     'Maximum threshold water equivalent for snow depletion', &
     &     'Maximum threshold snowpack water equivalent below'// &
     &     ' which the snow-covered-area curve is applied', &
     &     'inches')/=0 ) CALL read_error(1, 'snarea_thresh')

      IF ( declparam(MODNAME, 'albset_rnm', 'one', 'real', &
     &     '0.6', '0.4', '1.0', &
     &     'Albedo reset - rain,  melt stage', &
     &     'Fraction of rain in a mixed precipitation event'// &
     &     ' above which the snow albedo is not reset; applied during'// &
     &     ' the snowpack melt stage', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'albset_rnm')

      IF ( declparam(MODNAME, 'albset_rna', 'one', 'real', &
     &     '0.8', '0.5', '1.0', &
     &     'Albedo reset - rain, accumulation stage', &
     &     'Fraction of rain in a mixed precipitation event'// &
     &     ' above which the snow albedo is not reset; applied during'// &
     &     ' the snowpack accumulation stage', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'albset_rna')

      IF ( declparam(MODNAME, 'albset_snm', 'one', 'real', &
     &     '0.2', '0.1', '1.0', &
     &     'Albedo reset - snow, melt stage', &
     &     'Minimum snowfall, in water equivalent, needed to reset'// &
     &     ' snow albedo during the snowpack melt stage', &
     &     'inches')/=0 ) CALL read_error(1, 'albset_snm')

      IF ( declparam(MODNAME, 'albset_sna', 'one', 'real', &
     &     '0.05', '0.01', '1.0', &
     &     'Albedo reset - snow, accumulation stage', &
     &     'Minimum snowfall, in water equivalent, needed to reset'// &
     &     ' snow albedo during the snowpack accumulation stage', &
     &     'inches')/=0 ) CALL read_error(1, 'albset_sna')

      ALLOCATE ( Emis_noppt(Nhru) )
      IF ( declparam(MODNAME, 'emis_noppt', 'nhru', 'real', &
     &     '0.757', '0.757', '1.0', &
     &     'Emissivity of air on days without precipitation for each HRU', &
     &     'Average emissivity of air on days without precipitation for each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'emis_noppt')

      ALLOCATE ( Cecn_coef(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'cecn_coef', 'nhru,nmonths', 'real', &
     &     '5.0', '0.02', '20.0', &
     &     'Monthly convection condensation energy coefficient for each HRU', &
     &     'Monthly (January to December) convection condensation energy coefficient for each HRU', &
     &     'calories per degree Celsius above 0')/=0 ) CALL read_error(1, 'cecn_coef')

      ALLOCATE ( Freeh2o_cap(Nhru) )
      IF ( declparam(MODNAME, 'freeh2o_cap', 'nhru', 'real', &
     &     '0.05', '0.01', '0.2', &
     &     'Free-water holding capacity of snowpack for each HRU', &
     &     'Free-water holding capacity of snowpack for each HRU, expressed as a'// &
     &     ' decimal fraction of the frozen water content of the snowpack (pk_ice)', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'freeh2o_cap')

      ALLOCATE ( Tstorm_mo(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'tstorm_mo', 'nhru,nmonths', 'integer', &
     &     '0', '0', '1', &
     &     'Set to 1 if thunderstorms prevalent during month for each HRU', &
     &     'Monthly flag (January to December) for prevalent storm'// &
     &     ' type for each HRU (0=frontal storms; 1=convective storms)', &
     &     'none')/=0 ) CALL read_error(1, 'tstorm_mo')

      IF ( Init_vars_from_file==OFF .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==3 ) THEN
        ALLOCATE ( Snowpack_init(Nhru) )
        IF ( declparam(MODNAME, 'snowpack_init', 'nhru', 'real', &
     &       '0.0', '0.0', '5000.0', &
     &       'Initial snowpack water equivalent in each HRU', &
     &       'Storage of snowpack in each HRU at the beginning of a simulation', &
     &       'inches')/=0 ) CALL read_error(1, 'snowpack_init')
      ENDIF

      END FUNCTION snodecl

!***********************************************************************
!     snoinit - Initialize snowcomp module - get parameter values,
!               compute initial values
!***********************************************************************
      INTEGER FUNCTION snoinit()
      USE PRMS_CONSTANTS, ONLY: LAND, GLACIER, FEET, FEET2METERS, DNEARZERO, ACTIVE, OFF, MONTHS_PER_YEAR, DEBUG_less
      USE PRMS_MODULE, ONLY: Nhru, Ndepl, Print_debug, Init_vars_from_file, Glacier_flag, Snarea_curve_flag
      USE PRMS_SNOW
      USE PRMS_BASIN, ONLY: Basin_area_inv, Hru_route_order, Active_hrus, Hru_area_dble, Elev_units, Hru_type
      USE PRMS_FLOWVARS, ONLY: Pkwater_equiv, Glacier_frac, Glrette_frac, Alt_above_ela
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE, ATAN, SNGL, MIN
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: read_error, sca_deplcrv, glacr_states_to_zero
! Local Variables
      INTEGER :: i, j
      REAL :: x
! Save Variables
      REAL, SAVE :: acum_init(MAXALB), amlt_init(MAXALB)
      DATA acum_init/.80, .77, .75, .72, .70, .69, .68, .67, .66, .65, .64, .63, .62, .61, .60/
      DATA amlt_init/.72, .65, .60, .58, .56, .54, .52, .50, .48, .46, .44, .43, .42, .41, .40/
!***********************************************************************
      snoinit = 0

      IF ( Glacier_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'glacr_freeh2o_cap', Nhru, 'real', Glacr_freeh2o_cap)/=0 ) CALL read_error(2, 'glacr_freeh2o_cap')
        IF ( getparam(MODNAME, 'albedo_ice', Nhru, 'real', Albedo_ice)/=0 ) CALL read_error(2, 'albedo_ice')
        IF ( getparam(MODNAME, 'albedo_coef', Nhru, 'real', Albedo_coef)/=0 ) CALL read_error(2, 'albedo_coef')
        IF ( getparam(MODNAME, 'glacr_layer', Nhru, 'real', Glacr_layer)/=0 ) CALL read_error(2, 'glacr_layer')
      ENDIF

      IF ( getparam(MODNAME, 'den_init', Nhru, 'real', Den_init)/=0 ) CALL read_error(2, 'den_init')
      IF ( getparam(MODNAME, 'den_max', Nhru, 'real', Den_max)/=0 ) CALL read_error(2, 'den_max')
      IF ( getparam(MODNAME, 'settle_const', Nhru, 'real', Settle_const)/=0 ) CALL read_error(2, 'settle_const')

      IF ( getparam(MODNAME, 'melt_look', Nhru, 'integer', Melt_look)/=0 ) CALL read_error(2, 'melt_look')
      IF ( getparam(MODNAME, 'melt_force', Nhru, 'integer', Melt_force)/=0 ) CALL read_error(2, 'melt_force')
      IF ( getparam(MODNAME, 'rad_trncf', Nhru, 'real', Rad_trncf)/=0 ) CALL read_error(2, 'rad_trncf')
      IF ( Snarea_curve_flag==OFF ) THEN
        IF ( getparam(MODNAME, 'hru_deplcrv', Nhru, 'integer', Hru_deplcrv)/=0 ) CALL read_error(2, 'hru_deplcrv')
        IF ( getparam(MODNAME, 'snarea_curve', Ndepl*11, 'real', Snarea_curve)/=0 ) CALL read_error(2, 'snarea_curve')
      ELSE
        IF ( getparam(MODNAME, 'snarea_a', Nhru, 'real', Snarea_a)/=0 ) CALL read_error(2, 'snarea_a')
        IF ( getparam(MODNAME, 'snarea_b', Nhru, 'real', Snarea_b)/=0 ) CALL read_error(2, 'snarea_b')
        IF ( getparam(MODNAME, 'snarea_c', Nhru, 'real', Snarea_c)/=0 ) CALL read_error(2, 'snarea_c')
        IF ( getparam(MODNAME, 'snarea_d', Nhru, 'real', Snarea_d)/=0 ) CALL read_error(2, 'snarea_d')
        DO i = 1, Nhru
          Hru_deplcrv(i) = i
          x = 0.0
          DO j = 1, 11
            Snarea_curve(j,i) = (Snarea_a(i) - Snarea_d(i)) / (1 + (x**Snarea_b(i) / Snarea_c(i))) + Snarea_d(i)
            IF ( Snarea_curve(j,i)>1.0 ) THEN
              IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, snarea_curve computed > 1.0 for HRU:', i, &
     &                                               '; value:', Snarea_curve(j,i)
              Snarea_curve(j,i) = 1.0
            ENDIF
!            write (777,*) snarea_curve(j,i), x
            x = x + 0.1
          ENDDO
        ENDDO
      ENDIF
      IF ( getparam(MODNAME, 'snarea_thresh', Nhru, 'real', Snarea_thresh)/=0 ) CALL read_error(2, 'snarea_thresh')
      IF ( getparam(MODNAME, 'albset_rnm', 1, 'real', Albset_rnm)/=0 ) CALL read_error(2, 'albset_rnm')
      IF ( getparam(MODNAME, 'albset_rna', 1, 'real', Albset_rna)/=0 ) CALL read_error(2, 'albset_rna')
      IF ( getparam(MODNAME, 'albset_sna', 1, 'real', Albset_sna)/=0 ) CALL read_error(2, 'albset_sna')
      IF ( getparam(MODNAME, 'albset_snm', 1, 'real', Albset_snm)/=0 ) CALL read_error(2, 'albset_snm')
      IF ( getparam(MODNAME, 'emis_noppt', Nhru, 'real', Emis_noppt)/=0 ) CALL read_error(2, 'emis_noppt')
      IF ( getparam(MODNAME, 'cecn_coef', Nhru*MONTHS_PER_YEAR, 'real', Cecn_coef)/=0 ) CALL read_error(2, 'cecn_coef')
      IF ( getparam(MODNAME, 'freeh2o_cap', Nhru, 'real', Freeh2o_cap)/=0 ) CALL read_error(2, 'freeh2o_cap')
      IF ( getparam(MODNAME, 'tstorm_mo', Nhru*MONTHS_PER_YEAR, 'integer', Tstorm_mo)/=0 ) CALL read_error(2, 'tstorm_mo')

      Pk_precip = 0.0
      Snowmelt = 0.0
      Snow_evap = 0.0
      Pptmix_nopack = OFF
      Tcal = 0.0
      Frac_swe = 0.0
      Acum = acum_init
      Amlt = amlt_init
      Basin_tcal = 0.0D0
      Basin_snowmelt = 0.0D0
      Basin_snowevap = 0.0D0
      Basin_pk_precip = 0.0D0
      Basin_glacrb_melt = 0.0D0
      Basin_glacrevap = 0.0D0
      IF ( Glacier_flag==ACTIVE ) THEN
        Glacrb_melt = 0.0
        Glacrmelt = 0.0
        Glacr_evap = 0.0
      ENDIF

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==3 ) THEN
        IF ( getparam(MODNAME, 'snowpack_init', Nhru, 'real', Snowpack_init)/=0 ) CALL read_error(2, 'snowpack_init')
        Pkwater_equiv = 0.0D0
        Pk_depth = 0.0D0
        Pk_den = 0.0
        Pk_ice = 0.0
        Freeh2o = 0.0
        Ai = 0.0D0
        Snowcov_area = 0.0
        Basin_pweqv = 0.0D0
        Basin_snowdepth = 0.0D0
        Basin_snowcov = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          Pkwater_equiv(i) = DBLE( Snowpack_init(i) )
          IF ( Pkwater_equiv(i)>0.0D0 ) THEN
            Basin_pweqv = Basin_pweqv + Pkwater_equiv(i)*Hru_area_dble(i)
            Pk_depth(i) = Pkwater_equiv(i)/DBLE(Den_init(i))
            Pk_den(i) = SNGL( Pkwater_equiv(i)/Pk_depth(i) )
            Pk_ice(i) = SNGL( Pkwater_equiv(i) )
            Freeh2o(i) = Pk_ice(i)*Freeh2o_cap(i)
            Ai(i) = Pkwater_equiv(i) ! [inches]
            IF ( Ai(i)>Snarea_thresh(i) ) Ai(i) = DBLE( Snarea_thresh(i) ) ! [inches]
            IF ( Ai(i)>DNEARZERO ) THEN
              Frac_swe(i) = SNGL( Pkwater_equiv(i)/Ai(i) ) ! [fraction]
              Frac_swe(i) = MIN( 1.0, Frac_swe(i) )
            ENDIF
            CALL sca_deplcrv(Snowcov_area(i), Snarea_curve(1,Hru_deplcrv(i)), Frac_swe(i))
            Basin_snowcov = Basin_snowcov + DBLE(Snowcov_area(i))*Hru_area_dble(i)
            Basin_snowdepth = Basin_snowdepth + Pk_depth(i)*Hru_area_dble(i)
          ENDIF
        ENDDO
        Basin_pweqv = Basin_pweqv*Basin_area_inv
        Basin_snowcov = Basin_snowcov*Basin_area_inv
        Basin_snowdepth = Basin_snowdepth*Basin_area_inv
        DEALLOCATE ( Snowpack_init )
        Pss = Pkwater_equiv
        Pst = Pkwater_equiv
      ENDIF

      IF ( Init_vars_from_file>0 ) RETURN

      Yrdays5 = 0
      Basin_snowicecov = 0.0D0
      Iasw = 0
      Iso = 1
      Mso = 1
      Lso = 0
      Pk_def = 0.0
      Pk_temp = 0.0
      Albedo = 0.0
      Snsv = 0.0
      Lst = 0
      Int_alb = 1
      Salb = 0.0
      Slst = 0.0
      Snowcov_areasv = Snowcov_area
      Scrv = 0.0D0
      Pksv = 0.0D0

      IF ( Glacier_flag==ACTIVE ) THEN ! do here when not a restart simulation
        IF ( getparam(MODNAME, 'glacier_frac_init', Nhru, 'real', Glacier_frac_init)/=0 ) CALL read_error(2, 'glacier_frac_init')
        Glacr_albedo = 0.0
        Glacier_frac = Glacier_frac_init
        IF ( getparam(MODNAME, 'glrette_frac_init', Nhru, 'real', Glrette_frac_init)/=0 ) CALL read_error(2, 'glrette_frac_init')
        Glrette_frac = Glrette_frac_init
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Glacier_frac(i)>0.0 ) THEN
            IF ( Hru_type(i)==GLACIER ) THEN
              IF ( Elev_units==FEET ) THEN !from Oerlemans 1992
                Glacr_albedo(i) = Albedo_ice(i) +(Albedo_coef(i)/PI)*ATAN( (Alt_above_ela(i)*FEET2METERS+300.0)/200.0 )
              ELSE
                Glacr_albedo(i) = Albedo_ice(i) +(Albedo_coef(i)/PI)*ATAN( (Alt_above_ela(i)+300.0)/200.0 )
              ENDIF
            ELSE
              PRINT *, 'Warning, glacier_frac > 0, but hru_type not equal to 4, glacier_frac set to 0'
              PRINT *, 'in HRU ', i, 'glacier_frac_init = ', Glacier_frac_init(i)
              Glacier_frac(i) = 0.0
            ENDIF
          ENDIF
          IF ( Glrette_frac(i)>0.0 ) THEN
            IF ( Hru_type(i)==LAND ) THEN
              Glacr_albedo(i) = Albedo_ice(i)
            ELSE
              PRINT *, 'Warning, glrette_frac > 0, but hru_type not equal to 1, glrette_frac set to 0'
              PRINT *, 'in HRU ', i, 'glrette_frac_init = ', Glrette_frac_init(i)
              Glrette_frac(i) = 0.0
            ENDIF
          ENDIF
        ENDDO
        DEALLOCATE ( Glacier_frac_init )

        Alt_above_ela = 0.0
        Ann_tempc = 0.0
        Glacr_air_5avtemp = 0.0
        Glacr_air_5avtemp1 = 0.0
        Glacr_air_deltemp = 0.0
        Glacr_5avsnow = 0.0
        Glacr_5avsnow1 = 0.0
        Glacr_delsnow = 0.0
        Glacr_pk_den = 0.0
        Glacr_pk_temp = 0.0
        Glacr_pk_ice = 0.0
        Glacr_pk_def = 0.0
        Glacr_pkwater_equiv = 0.0D0
        Glacr_pkwater_ante = 0.0D0
        Glacr_freeh2o = 0.0
        Glacr_pk_depth = 0.0D0
        Glacr_pst = 0.0D0
        Glacr_pss = 0.0D0
        Glacrcov_area = 0.0
        Glacr_freeh2o_capm = Glacr_freeh2o_cap
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Glacier_frac(i)>0.0 .AND. Hru_type(i)==GLACIER ) CALL glacr_states_to_zero(i,1)
        ENDDO
      ENDIF

      END FUNCTION snoinit

!***********************************************************************
!     snorun - daily mode snow estimates
!***********************************************************************
      INTEGER FUNCTION snorun()
      USE PRMS_CONSTANTS, ONLY: LAKE, LAND, GLACIER, SHRUBS, FEET, &
     &    INCH2M, FEET2METERS, DNEARZERO, ACTIVE, OFF, DEBUG_less, DAYS_YR
      USE PRMS_MODULE, ONLY: Nhru, Print_debug, Glacier_flag, Start_year, &
     &    Nowyear, Nowmonth, Albedo_cbh_flag, snow_cloudcover_flag
      USE PRMS_SNOW
      USE PRMS_SOLTAB, ONLY: Soltab_horad_potsw, Soltab_potsw, Hru_cossl
      USE PRMS_CLIMATE_HRU, ONLY: Albedo_hru
      USE PRMS_BASIN, ONLY: Hru_area, Active_hrus, Hru_type, &
     &    Basin_area_inv, Hru_route_order, Cov_type, Elev_units
      USE PRMS_CLIMATEVARS, ONLY: Newsnow, Pptmix, Orad, Basin_horad, Potet_sublim, &
     &    Hru_ppt, Prmx, Tmaxc, Tminc, Tavgc, Swrad, Potet, Transp_on, Tmax_allsnow_c
      USE PRMS_FLOWVARS, ONLY: Pkwater_equiv, Glacier_frac, Glrette_frac, Alt_above_ela
      USE PRMS_SET_TIME, ONLY: Jday, Julwater
      USE PRMS_INTCP, ONLY: Net_rain, Net_snow, Net_ppt, Canopy_covden, Hru_intcpevap
      IMPLICIT NONE
! Functions
      EXTERNAL :: ppt_to_pack, snowcov, snalbedo, snowbal, snowevap, glacr_states_to_zero
      INTRINSIC :: ABS, SQRT, DBLE, SNGL, EXP, DABS, MOD, ATAN
! Local Variables
      INTEGER :: i, j, k, niteda, isglacier
      REAL :: trd, sw, effk, cst, temp, cals, emis, esv, swn, cec
      REAL :: ieffk, icst, icals, isw, iswn, frac, orad_local
      DOUBLE PRECISION :: dpt1, dpt_before_settle
!***********************************************************************
      snorun = 0

      ! Set the basin totals to 0
      ! (recalculated at the end of the time step)
      Basin_snowmelt = 0.0D0
      Basin_pweqv = 0.0D0
      Basin_snowevap = 0.0D0
      Basin_snowcov = 0.0D0
      Basin_pk_precip = 0.0D0
      Basin_snowdepth = 0.0D0
      Basin_tcal = 0.0D0
      IF ( Glacier_flag==ACTIVE ) THEN
        Basin_snowicecov = 0.0D0
        Basin_glacrb_melt = 0.0D0
        Basin_glacrevap = 0.0D0
        IF ( Julwater==1 .AND. MOD(Nowyear-Start_year,5)==0 ) Yrdays5 = 0
        Glacr_pkwater_ante = Glacr_pkwater_equiv
      ENDIF

      ! Calculate the ratio of measured radiation to potential radiation
      ! (used as a cumulative indicator of cloud cover)
      IF ( snow_cloudcover_flag==OFF ) trd = Orad/SNGL(Basin_horad) ! [dimensionless ratio] ! old basin-wide equation
      Frac_swe = 0.0
      ! By default, the precipitation added to snowpack, snowmelt,
      ! and snow evaporation are 0
      Pk_precip = 0.0 ! [inches]
      Snowmelt = 0.0 ! [inches]
      Snow_evap = 0.0 ! [inches]
      Tcal = 0.0
      Ai = 0.0D0
      ! Keep track of the pack water equivalent before it is changed
      ! by precipitation during this time step
      Pkwater_ante = Pkwater_equiv

      ! Loop through all the active HRUs, in routing order
      DO j = 1, Active_hrus
        i = Hru_route_order(j) ! [counter]

        ! Skip the HRU if it is a lake
        IF ( Hru_type(i)==LAKE ) CYCLE

        Active_glacier = OFF
        isglacier = OFF
        IF ( Glacier_flag==ACTIVE ) THEN
          IF ( Hru_type(i)==GLACIER .OR. Hru_type(i)==LAND ) THEN
            Glacrmelt(i) = 0.0 ! [inches]
            Glacrb_melt(i) = 0.0 ! [inches]
            Glacr_evap(i) = 0.0 ! [inches]
            IF ( Glacier_frac(i)>0.0 .OR. Glrette_frac(i)>0.0 ) THEN
              IF ( Glacier_frac(i)>0.0 ) Active_glacier = 1
              IF ( Glrette_frac(i)>0.0 ) Active_glacier = 2
              Glacr_pk_den(i) = 0.917
              ! if no active layer make 0 deg and no holding capacity at start of each day
              IF ( Glacr_layer(i)==0.0 .OR. Glacr_pk_depth(i)>1.0D3 ) THEN
                Glacr_pk_def(i) = 0.0
                Glacr_pk_temp(i) = 0.0
                Glacr_freeh2o_capm(i) = 0.0
              ENDIF
            ELSE !zero out states for glacier if gone (glacier state changes in glacier module, not here)
              Glacr_pkwater_equiv(i) = 0.D0
              Glacrcov_area(i) = 0.0
              Glacr_pk_def(i) = 0.0
              Glacr_pk_temp(i) = 0.0
              Glacr_pk_ice(i) = 0.0
              Glacr_freeh2o(i) = 0.0
              Glacr_pk_depth(i) = 0.D0
              Glacr_pss = 0.0D0
              Glacr_pst(i) = 0.0D0
              Glacr_pk_den(i) = 0.0
              Glacr_freeh2o_capm(i) = 0.0
              Glacr_albedo(i) = 0.0
            ENDIF
            isglacier = ACTIVE
          ENDIF
        ENDIF

        ! If it's the first julian day of the water year, several
        ! variables need to be reset
        ! - reset the previous snow water eqivalent plus new snow to 0
        ! - reset flags to indicate it is not melt season or potetential melt season
        ! - reset the counter for the number of days a snowpack is at 0 deg Celsius
        !rsr, do we want to reset all HRUs, what about Southern Hemisphere
        IF ( Julwater==1 ) THEN
          Pss(i) = 0.0D0 ! [inches]
          Iso(i) = 1 ! [flag]
          Mso(i) = 1 ! [flag]
          Lso(i) = 0 ! [counter]

          IF ( Active_glacier>OFF ) THEN
            CALL glacr_states_to_zero(i,1) !all snow on glacier becomes firn, reset active layer thickness
            IF ( Active_glacier==1 ) THEN
! If Active_glacier>OFF we are zeroing out snowpack if have glacierettes even though possibly a lot of HRU is not glacierized.
! If Active_glacier==1 do not zero out glacierettes, but then will maybe never melt ice on glacierettes. If the climate is
!	correct the snowpack will deplete quick because there is a lot of lower elevation than the glacierette included in the HRU.
! Choice does not effect runoff much, but will effect Basin_pweqv and things like that
              ! if terminus glacier, and has snow will disappear off glacier but that is likely anyhow
              Pkwater_equiv(i) = 0.0
              Pk_depth(i) = 0.0D0
              Pss(i) = 0.0D0
              Snsv(i) = 0.0
              Lst(i) = 0
              Pst(i) = 0.0D0
              Iasw(i) = 0
              Pk_den(i) = 0.0
              Snowcov_area(i) = 0.0
              Pk_def(i) = 0.0
              Pk_temp(i) = 0.0
              Pk_ice(i) = 0.0
              Freeh2o(i) = 0.0
              Snowcov_areasv(i) = 0.0 ! rsr, not in original code
              IF ( Elev_units==FEET ) THEN !from Oerlemans 1992
                Glacr_albedo(i) = Albedo_ice(i) +(Albedo_coef(i)/PI)*ATAN( (Alt_above_ela(i)*FEET2METERS+300.0)/200.0 )
              ELSE
                Glacr_albedo(i) = Albedo_ice(i) +(Albedo_coef(i)/PI)*ATAN( (Alt_above_ela(i)+300.0)/200.0 )
              ENDIF
            ELSE !IF ( Active_glacier==2 )
              Glacr_albedo(i) = Albedo_ice(i) !glacr_albedo doesn't change if glacierette but could get zeroed out
            ENDIF
          ENDIF
          IF ( isglacier==ACTIVE ) THEN
            IF (Nowyear >= Start_year+10 .AND. MOD(Nowyear-Start_year,5)==0 ) THEN
              Glacr_air_deltemp(i) = Glacr_air_5avtemp1(i) - Glacr_air_5avtemp(i) !need 5 years of data
              IF ( Glacr_5avsnow1(i)>0.0 ) THEN
                Glacr_delsnow(i) = 10.0*(Glacr_5avsnow1(i) - Glacr_5avsnow(i))/Glacr_5avsnow1(i) !number of 10 percent (*100.0/10.0) changes
              ELSE
                Glacr_delsnow(i) = 0.0
              ENDIF
            ENDIF
            !keep before restart
            IF ( MOD(Nowyear-Start_year,5)==0 ) THEN
              IF ( Nowyear-Start_year==5 ) THEN
                Glacr_air_5avtemp1(i) = Glacr_air_5avtemp(i)
                Glacr_5avsnow1(i) = Glacr_5avsnow(i)
              ENDIF
              Glacr_air_5avtemp(i) = 0.0 !zero out for new year restart
              Glacr_5avsnow(i) = 0.0 !zero out for new year restart
            ENDIF
            Ann_tempc(i) = 0.0 !zero out for new year restart
          ENDIF !end start of year calculations
        ENDIF

! Do for summer
        IF ( isglacier==ACTIVE ) THEN
          IF (Julwater>151 .AND. Julwater<244) THEN ! Now following McGrath et al 2017, temp June-August, 92 days
            Yrdays5 = Yrdays5 + 1
            Glacr_air_5avtemp(i) = ( Glacr_air_5avtemp(i)*(Yrdays5-1)+ Tavgc(i) )/Yrdays5
          ENDIF
! Do for every time step
          Ann_tempc(i) = ( Ann_tempc(i)*(Julwater-1)+ Tavgc(i) )/Julwater
          Glacr_5avsnow(i) = Glacr_5avsnow(i) + Net_snow(i)/5.0
        ENDIF

        ! HRU SET-UP - SET DEFAULT VALUES AND/OR BASE
        !              CONDITIONS FOR THIS TIME PERIOD
        !**************************************************************

        ! By default, there has not been a mixed event without a
        ! snowpack
        Pptmix_nopack(i) = OFF ! [flag]

        ! If the day of the water year is beyond the forced melt day
        ! indicated by the parameter, then set the flag indicating
        ! melt season
        !rsr, need to rethink this at some point
!rsr10  IF ( Iso(i)/=2 ) THEN
          IF ( Jday==Melt_force(i) ) Iso(i) = 2 ! [flag]
!rsr10  ENDIF

        ! If the day of the water year is beyond the first day to
        ! look for melt season indicated by the parameter,
        ! then set the flag indicating to watch for melt season
        !rsr, need to rethink this at some point
!rsr10  IF ( Mso(i)/=2 ) THEN
          IF ( Jday==Melt_look(i) ) Mso(i) = 2 ! [flag]
!rsr10  ENDIF

        ! Skip the HRU if there is no snowpack and no new snow and not a glacier
        IF ( Active_glacier==OFF ) THEN
          IF ( Pkwater_equiv(i)<DNEARZERO .AND. Newsnow(i)==0 ) THEN
            Snowcov_area(i) = 0.0 ! reset to be sure it is zero if snowpack melted on last timestep
            CYCLE
          ENDIF
        ENDIF

        ! If there is no existing snow pack and there is new snow, the
        ! initial snow covered area is complete (1)
        IF ( Newsnow(i)==ACTIVE .AND. Pkwater_equiv(i)<DNEARZERO ) Snowcov_area(i) = 1.0 ! [fraction of area]
        IF ( Active_glacier==1 ) Glacrcov_area(i) =(1.0-Snowcov_area(i))*Glacier_frac(i)
        IF ( Active_glacier==2 ) Glacrcov_area(i) =(1.0-Snowcov_area(i))*Glrette_frac(i)

        ! HRU STEP 1 - DEAL WITH PRECIPITATION AND ITS EFFECT ON THE WATER
        !              CONTENT AND HEAT CONTENT OF SNOW PACK
        !************************************************************

        ! If there is net precipitation on an existing snowpack, OR if
        ! there is any net snow, add the incoming water (or ice) and
        ! heat (or heat deficit) to the snowpack
        IF ( (Pkwater_equiv(i)>0.0D0.AND.Net_ppt(i)>0.0) .OR. Net_snow(i)>0.0 ) &
     &       CALL ppt_to_pack(Pptmix(i), Iasw(i), Tmaxc(i), Tminc(i), Tavgc(i), &
     &                        Pkwater_equiv(i), Net_rain(i), Pk_def(i), &
     &                        Pk_temp(i), Pk_ice(i), Freeh2o(i), Snowcov_area(i), &
     &                        Snowmelt(i), Pk_depth(i), Pss(i), Pst(i), Net_snow(i), &
     &                        Pk_den(i), Pptmix_nopack(i), Pk_precip(i), Tmax_allsnow_c(i,Nowmonth), &
     &                        Freeh2o_cap(i), Den_max(i), not_a_glacier_hru)
        IF ( Active_glacier>OFF ) THEN
          IF ( Glacrcov_area(i)>0.0.AND.Glacr_pkwater_ante(i)>0.0D0.AND.Net_ppt(i)>0.0 &
     &         .AND.Pptmix(i)==0.AND.Net_snow(i)==0.0 ) THEN
             CALL ppt_to_pack(0, Iasw(i), Tmaxc(i), Tminc(i), Tavgc(i), &
     &                        Glacr_Pkwater_equiv(i), Net_rain(i), Glacr_pk_def(i), &
     &                        Glacr_pk_temp(i), Glacr_pk_ice(i), Glacr_freeh2o(i), Glacrcov_area(i), &
     &                        Glacrmelt(i), Glacr_pk_depth(i), Glacr_pss(i), Glacr_pst(i), 0.0, &
     &                        Glacr_pk_den(i), Pptmix_nopack(i), Pk_precip(i), Tmax_allsnow_c(i,Nowmonth), &
     &                        Glacr_freeh2o_capm(i), Den_max(i), i)
          ENDIF
          ! FOLLOWING does basal melt on glacier
          ! Paterson 2010 says 12 mm/yr for friction and geothermal heating
          IF ( Active_glacier==1 ) Glacrb_melt(i) = 12.0*0.03937/DAYS_YR*Glacier_frac(i)
          IF ( Active_glacier==2 ) Glacrb_melt(i) = 12.0*0.03937/DAYS_YR*Glrette_frac(i) !since not moving much, maybe =0
        ENDIF

        ! If there is still a snowpack
        IF ( Pkwater_equiv(i)>0.0D0 ) THEN

          ! HRU STEP 2 - CALCULATE THE NEW SNOW COVERED AREA
          !**********************************************************
          ! Compute snow-covered area from depletion curve
          k = Hru_deplcrv(i)
          ! calculate the new snow covered area
          CALL snowcov(Iasw(i), Newsnow(i), Snowcov_area(i), &
     &                 Snarea_curve(1, k), Pkwater_equiv(i), Pst(i), &
     &                 Snarea_thresh(i), Net_snow(i), Scrv(i), &
     &                 Pksv(i), Snowcov_areasv(i), Ai(i), Frac_swe(i))

          ! HRU STEP 3 - COMPUTE THE NEW ALBEDO
          !**********************************************************

          ! Compute albedo if there is any snowpack
          CALL snalbedo(Newsnow(i), Iso(i), Lst(i), Snsv(i), &
     &                  Prmx(i), Pptmix(i), Albset_rnm, Net_snow(i), &
     &                  Albset_snm, Albset_rna, Albset_sna, Albedo(i), &
     &                  Int_alb(i), Salb(i), Slst(i))
        ENDIF

        IF ( Active_glacier>OFF ) THEN
          IF ( Active_glacier==1 ) Glacrcov_area(i) =(1.0-Snowcov_area(i))*Glacier_frac(i)
          IF ( Active_glacier==2 ) Glacrcov_area(i) =(1.0-Snowcov_area(i))*Glrette_frac(i)
! Albedo so transition snow to ice smooothly, see Oerlemans 1992, this is albedo if snowcovered ice too
! Albedo can be input in a CBH File when albedo_cbh_flag = ACTIVE
          IF ( Albedo_cbh_flag==OFF ) THEN
            Albedo(i) = Albedo(i) - (Albedo(i)-Glacr_albedo(i))*EXP(-5.0*SNGL(Pkwater_equiv(i))*INCH2M)
            IF ( Albedo(i)<0.08 ) Albedo(i)=0.08 !See Brock 2000
            IF ( Albedo(i)>0.92 ) Albedo(i)=0.92 !See Brock 2000
          ELSE
            Albedo(i) = Albedo_hru(i)
          ENDIF
        ENDIF

        ! If there is still a snowpack or glacier
        IF ( Pkwater_equiv(i)>0.0D0 .OR. Active_glacier>OFF ) THEN

          ! HRU STEP 4 - DETERMINE RADIATION FLUXES AND SNOWPACK
          !              STATES NECESSARY FOR ENERGY BALANCE
          !**********************************************************

          ! Set the emissivity of the air to the emissivity when there
          ! is no precipitation
          emis = Emis_noppt(i) ! [fraction of radiation]
          ! Could use equation from Swinbank 63 using Temp, a is -13.638, b is 6.148
          ! temparature is halfway between the minimum and average temperature for the day
          !temp = (Tminc(i)+Tavgc(i))*0.5
          !emis = ((temp+273.16)**(Emis_coefb-4.0))*(10.0**(Emis_coefa+1.0))/5.670373E−8 ! /by Stefan Boltzmann in SI units
          ! If there is any precipitation in the HRU, reset the
          ! emissivity to 1
          IF ( Hru_ppt(i)>0.0 ) emis = 1.0 ! [fraction of radiation]
          ! Save the current value of emissivity
          esv = emis ! [fraction of radiation]
          ! Set the convection-condensation for a half-day interval
          cec = Cecn_coef(i, Nowmonth)*0.5 ! [cal/(cm^2 degC)]
                                           ! or [Langleys / degC]
          ! If the land cover is trees, reduce the convection-
          ! condensation parameter by half
          IF ( Cov_type(i)>SHRUBS ) cec = cec*0.5 ! [cal/(cm^2 degC)] cov_type>2 is valid for trees and coniferous)
                                                  ! or [Langleys / degC]
          ! Check whether to force spring melt
          ! Spring melt is forced if time is before the melt-force
          ! day and after the melt-look day (parameters)
          ! If between these dates, the spring melt applies if the
          ! snowpack temperature is above or equal to 0
          ! for more than 4 cycles of the snorun function

          ! If before the first melt-force day
          IF ( Iso(i)==1 ) THEN
            ! If after the first melt-look day
            IF ( Mso(i)==2 ) THEN

              ! Melt season is determined by the number of days the
              ! snowpack is above 0 degrees C.  The first time that
              ! the snowpack is isothermal at 0 degrees C for more
              ! than 4 days is the beginning of snowmelt season.
              ! 2 options below (if-then, else)

              ! (1) The snowpack temperature is 0 degrees
              IF ( Pk_temp(i)>=0.0 ) THEN
                ! Increment the number of days that the snowpack
                ! has been isothermal at 0 degrees C
                Lso(i) = Lso(i) + 1 ! [days]
                ! If the snowpack temperature has been 0 or greater
                ! for more than 4 cycles
                IF ( Lso(i)>4 ) THEN
                  ! Set the melt-force flag and reset counter
                  Iso(i) = 2 ! [flag]
                  Lso(i) = 0 ! [days]
                ENDIF

              ! (2) The snowpack temperature is less than 0 degrees
              ELSE
                ! Reset the counter for days snowpack temperature is above 0
                Lso(i) = 0 ! [days]
              ENDIF
            ENDIF
          ENDIF

          ! Compute energy balance for night period
          ! niteda is a flag indicating nighttime (1) or daytime (2)
          ! set the flag indicating night time
          niteda = 1 ! [flag]
          ! temparature is halfway between the minimum and average temperature
          ! for the day
          temp = (Tminc(i)+Tavgc(i))*0.5

          IF ( Pkwater_equiv(i)>0.0D0 ) THEN
            ! The incoming shortwave radiation is the HRU radiation
            ! adjusted by the albedo (some is reflected back into the
            ! atmoshphere) and the transmission coefficient (some is
            ! intercepted by the winter vegetative canopy)
            swn = Swrad(i)*(1.0-Albedo(i))*Rad_trncf(i) ! [cal/cm^2]
                                                        ! or [Langleys]
            ! Calculate the new snow depth (Riley et al. 1973)
            ! RSR: the following 3 lines of code were developed by Rob Payn, 7/10/2013
            ! The snow depth depends on the previous snow pack water
            ! equivalent plus the new net snow
            Pss(i) = Pss(i) + DBLE( Net_snow(i) ) ! [inches]
            dpt_before_settle = Pk_depth(i) + DBLE(Net_snow(i))/DBLE(Den_init(i))
            dpt1 = dpt_before_settle + DBLE(Settle_const(i)) * ((Pss(i)/DBLE(Den_max(i))) - dpt_before_settle)
!            dpt1 = Pk_depth(i) + (Net_snow(i)/DBLE(Den_init(i))) + &
!                   DBLE(Settle_const(i)) * ((Pss(i)/DBLE(Den_max(i))) - Pk_depth(i))
!            dpt1 = ((Net_snow(i)/DBLE(Den_init(i)))+ (Settle_const(i)/Den_max(i)*Pss(i))+Pk_depth(i))*(1.0/(1.0+Settle_const(i))) ! [inches]
            ! RAPCOMMENT - CHANGED TO THE APPROPRIATE FINITE DIFFERENCE
            !             APPROXIMATION OF SNOW DEPTH
            Pk_depth(i) = dpt1 ! [inches]

            ! Calculate the snowpack density
            IF ( dpt1>0.0D0 ) THEN
              Pk_den(i) = SNGL( Pkwater_equiv(i)/dpt1 )
            ELSE
              Pk_den(i) = 0.0
            ENDIF
                                 ! [inch water equiv / inch depth]

            ! The effective thermal conductivity is approximated
            ! (empirically) as 0.0077 times (snowpack density)^2
            ! [cal / (sec g degC)] Therefore, the effective
            ! conductivity term (inside the square root) in the
            ! equation for conductive heat exchange can be
            ! calculated as follows (0.0077*pk_den^2)/(pk_den*0.5)
            ! where 0.5 is the specific heat of ice [cal / (g degC)]
            ! this simplifies to the following
            effk = 0.0154*Pk_den(i) ! [unitless]
            ! 13751 is the number of seconds in 12 hours over pi
            ! So for a half day, to calculate the conductive heat
            ! exchange per cm snow per cm^2 area per degree
            ! temperature difference is the following
            ! In effect, multiplying cst times the temperature
            ! gradient gives the heatexchange by heat conducted
            ! (calories) per square cm of snowpack
            cst = Pk_den(i)*(SQRT(effk*13751.0)) ! [cal/(cm^2 degC)]
                                                 ! or [Langleys / degC]

            ! no shortwave (solar) radiation at night
            sw = 0.0 ! [cal / cm^2] or [Langleys]
! new equation for trd
            ! Calculate the ratio of measured radiation to potential radiation
            ! (used as a cumulative indicator of cloud cover)
            IF ( snow_cloudcover_flag==ACTIVE ) THEN
              orad_local = SNGL( (DBLE(Swrad(i))*Hru_cossl(i)*Soltab_horad_potsw(Jday,i))/Soltab_potsw(Jday,i) )
              trd = orad_local/SNGL(Soltab_horad_potsw(Jday,i)) ! [dimensionless ratio]
            ENDIF
            ! calculate the night time energy balance
            CALL snowbal(niteda, Tstorm_mo(i,Nowmonth), Iasw(i), &
     &                   temp, esv, Hru_ppt(i), trd, Emis_noppt(i), &
     &                   Canopy_covden(i), cec, Pkwater_equiv(i), &
     &                   Pk_def(i), Pk_temp(i), Pk_ice(i), Freeh2o(i), &
     &                   Snowcov_area(i), Snowmelt(i), Pk_depth(i), &
     &                   Pss(i), Pst(i), Pk_den(i), cst, cals, sw, &
     &                   Freeh2o_cap(i), Den_max(i), not_a_glacier_hru)
            ! track total heat flux from both night and day periods
            Tcal(i) = cals ! [cal/cm^2] or [Langleys]
          ENDIF
          iswn  = 0.0
          IF ( Active_glacier>OFF ) THEN
            IF ( Glacrcov_area(i)>0.0 ) THEN
              iswn = Swrad(i)*(1.0-Glacr_albedo(i))*Rad_trncf(i) ! [cal/cm^2] !want bare ice albedo
                                                                 ! or [Langleys]
              ! Calculate the Glacier icepack density
              !
              ! The effective thermal conductivity is approximated
              ! (empirically) as 0.0077 times (snowpack density)^2 cal/(cm sec degC)
              ! from Oke 1987
              ! ice is 2.1 W/(m degC) = 0.021 W/(cm deg C) = 0.00502 cal/(cm sec degC)
              ! = 0.00597 times (0.917**2),
              ! firn (old snow density .5) is closer to 0.0042 W/(cm deg C) = 0.00401 times (0.5**2)
              !  Therefore, the effective
              ! conductivity term (inside the square root) in the
              ! equation for conductive heat exchange can be
              ! calculated as follows (0.0597*pk_den^2)/(pk_den*0.5)
              ! where 0.5 is the specific heat of ice [cal / (g degC)]
              ! this simplifies to the following
              ! might want to use 0.005*2 = 0.01 half way between if doing mix of firn and ice
              ieffk = 0.01194*Glacr_pk_den(i) ! [unitless]
              icst = Glacr_pk_den(i)*(SQRT(ieffk*13751.0)) ! [cal/(cm^2 degC)]
                                                           ! or [Langleys / degC]
              isw = 0.0 ! [cal / cm^2] or [Langleys]
              CALL snowbal(niteda, Tstorm_mo(i,Nowmonth), Iasw(i), &
     &                     temp, esv, Hru_ppt(i), trd, Emis_noppt(i), &
     &                     Canopy_covden(i), cec, Glacr_pkwater_equiv(i), &
     &                     Glacr_pk_def(i), Glacr_pk_temp(i), Glacr_pk_ice(i), Glacr_freeh2o(i), &
     &                     Glacrcov_area(i), Glacrmelt(i), Glacr_pk_depth(i), &
     &                     Glacr_pss(i), Glacr_pst(i), Glacr_pk_den(i), icst, icals, isw, &
     &                     Glacr_freeh2o_capm(i), Den_max(i), i)
            ENDIF
          ENDIF

          ! Compute energy balance for day period
          ! set the flag indicating daytime
          niteda = 2 ! [flag]
          ! temparature is halfway between the maximum and average
          ! temperature for the day
          temp = (Tmaxc(i)+Tavgc(i))*0.5 ! [degrees C]

          IF ( Pkwater_equiv(i)>0.0D0 ) THEN !(if the snowpack still exists)
            ! set shortwave radiation as calculated earlier
            sw = swn ! [cal/cm^2] or [Langleys]
            CALL snowbal(niteda, Tstorm_mo(i,Nowmonth), Iasw(i), &
     &                   temp, esv, Hru_ppt(i), trd, Emis_noppt(i), &
     &                   Canopy_covden(i), cec, Pkwater_equiv(i), &
     &                   Pk_def(i), Pk_temp(i), Pk_ice(i), Freeh2o(i), &
     &                   Snowcov_area(i), Snowmelt(i), Pk_depth(i), &
     &                   Pss(i), Pst(i), Pk_den(i), cst, cals, sw, &
     &                   Freeh2o_cap(i), Den_max(i), not_a_glacier_hru)
            ! track total heat flux from both night and day periods
            Tcal(i) = Tcal(i) + cals ! [cal/cm^2] or [Langleys]
          ENDIF
          ! Compute energy balance for day period (if glacier exists)
          IF ( Active_glacier>OFF ) THEN
            IF ( Glacrcov_area(i)>0.0 ) THEN
              ! set shortwave radiation as calculated earlier
              isw = iswn ! [cal/cm^2] or [Langleys]
              CALL snowbal(niteda, Tstorm_mo(i,Nowmonth), Iasw(i), &
     &                     temp, esv, Hru_ppt(i), trd, Emis_noppt(i), &
     &                     Canopy_covden(i), cec, Glacr_pkwater_equiv(i), &
     &                     Glacr_pk_def(i), Glacr_pk_temp(i), Glacr_pk_ice(i), Glacr_freeh2o(i), &
     &                     Glacrcov_area(i), Glacrmelt(i), Glacr_pk_depth(i), &
     &                     Glacr_pss(i), Glacr_pst(i), Glacr_pk_den(i), icst, icals, isw, &
     &                     Glacr_freeh2o_capm(i), Den_max(i), i)
            ENDIF
          ENDIF

          !  HRU STEP 5 - CALCULATE SNOWPACK LOSS TO EVAPORATION
          !********************************************************

          ! Compute snow evaporation (if there is still a snowpack)
          ! Some of the calculated evaporation can come from interception
          ! rather than the snowpack.  Therefore, the effects of
          ! interception must be evaluated.
          IF ( Pkwater_equiv(i)>0.0D0 ) THEN
            ! Snow can evaporate when transpiration is not occuring
            ! or when transpiration is occuring with cover types of
            ! bare soil or grass
            IF ( Transp_on(i)==0 .OR. (Transp_on(i)==1 .AND. Cov_type(i)<SHRUBS) ) & ! cov_type < 2
     &           CALL snowevap(Potet_sublim(i), Potet(i), Snowcov_area(i), &
     &                         Snow_evap(i), Pkwater_equiv(i), Pk_ice(i), &
     &                         Pk_def(i), Freeh2o(i), Pk_temp(i), Hru_intcpevap(i))
          ELSEIF ( Pkwater_equiv(i)<0.0D0 ) THEN
            IF ( Print_debug>DEBUG_less ) THEN
              IF ( Pkwater_equiv(i)<-DNEARZERO ) PRINT *, 'snowpack issue 3, negative pkwater_equiv, &
     &             HRU:', i, ' value:', Pkwater_equiv(i)
            ENDIF
            Pkwater_equiv(i) = 0.0D0 ! just to be sure negative values are ignored
          ENDIF
          IF ( Active_glacier>OFF ) THEN
            IF ( Glacrcov_area(i)>0.0 ) &
     &           CALL snowevap(Potet_sublim(i), Potet(i), Glacrcov_area(i), &
     &                         Glacr_evap(i), Glacr_pkwater_equiv(i), Glacr_pk_ice(i), &
     &                         Glacr_pk_def(i), Glacr_freeh2o(i), Glacr_pk_temp(i), Hru_intcpevap(i))
          ENDIF

          !  HRU CLEAN-UP - ADJUST FINAL HRU SNOWPACK STATES AND
          !                 INCREMENT THE BASIN TOTALS
          !*********************************************************

          ! Final state of the snowpack depends on whether it still
          ! exists after all the processing above
          ! 2 options below (if-then, else)

          ! (1) Snow pack still exists
          IF ( Pkwater_equiv(i)>0.0D0 ) THEN
            ! Snowpack still exists
            IF ( Pk_den(i)>0.0 ) THEN
              Pk_depth(i) = Pkwater_equiv(i)/DBLE(Pk_den(i))
            ELSE
              Pk_den(i) = Den_max(i)
              Pk_depth(i) = Pkwater_equiv(i)/DBLE(Den_max(i))
            ENDIF
            Pss(i) = Pkwater_equiv(i)
            ! If it is during the melt period and snowfall was
            ! insufficient to reset albedo, then reduce the cumulative
            ! new snow by the amount melted during the period
            ! (but don't let it be negative)
            IF ( Lst(i)>0 ) THEN
              Snsv(i) = Snsv(i) - Snowmelt(i)
              IF ( Snsv(i)<0.0 ) Snsv(i) = 0.0
            ENDIF
          ENDIF

        ENDIF  ! ends iff statement between steps 3 and 4 about 265 lines above


! LAST check to clear out all arrays if packwater is gone
        IF ( Pkwater_equiv(i)<=0.0D0 ) THEN
          IF ( Print_debug>DEBUG_less ) THEN
            IF ( Pkwater_equiv(i)<-DNEARZERO ) &
     &           PRINT *, 'Snowpack problem, pkwater_equiv negative, HRU:', i, ' value:', Pkwater_equiv(i)
          ENDIF
          Pkwater_equiv(i) = 0.0D0 ! just to be sure negative values are ignored
          ! Snowpack has been completely depleted, reset all states
          ! to no-snowpack values
          Pk_depth(i) = 0.0D0
          Pss(i) = 0.0D0
          Snsv(i) = 0.0
          Lst(i) = 0
          Pst(i) = 0.0D0
          Iasw(i) = 0
          Albedo(i) = 0.0
          Pk_den(i) = 0.0
          Snowcov_area(i) = 0.0
          Pk_def(i) = 0.0
          Pk_temp(i) = 0.0
          Pk_ice(i) = 0.0
          Freeh2o(i) = 0.0
          Snowcov_areasv(i) = 0.0 ! rsr, not in original code
          Ai(i) = 0.0D0
          Frac_swe(i) = 0.0
          Scrv(i) = 0.0D0
          Pksv(i) = 0.0D0
        ENDIF
        frac = 1.0
        IF ( Active_glacier>OFF ) THEN
          IF ( Glacr_pkwater_equiv(i)>0.0D0 ) THEN
            Glacr_pk_depth(i) = Glacr_pkwater_equiv(i)/DBLE(Glacr_pk_den(i))
          ELSE
            CALL glacr_states_to_zero(i,0)
          ENDIF
          Basin_glacrb_melt = Basin_glacrb_melt + DBLE( Glacrb_melt(i)*Hru_area(i) )
          Basin_glacrevap = Basin_glacrevap + DBLE( Glacr_evap(i)*Hru_area(i) )
          IF ( Active_glacier==1 ) frac = (1.0 - Glacier_frac(i))
          IF ( Active_glacier==2 ) frac = (1.0 - Glrette_frac(i))
        ENDIF

        ! Sum volumes for basin totals
        Basin_snowmelt = Basin_snowmelt + DBLE( Snowmelt(i)*Hru_area(i)*frac ) !don't include stuff melting into glacier
        Basin_pweqv = Basin_pweqv + Pkwater_equiv(i)*DBLE( Hru_area(i) )
        Basin_snowevap = Basin_snowevap + DBLE( Snow_evap(i)*Hru_area(i) )
        Basin_snowcov = Basin_snowcov + DBLE( Snowcov_area(i)*Hru_area(i) )
        Basin_pk_precip = Basin_pk_precip + DBLE( Pk_precip(i)*Hru_area(i) )
        Basin_snowdepth = Basin_snowdepth + Pk_depth(i)*DBLE( Hru_area(i) )
        Basin_tcal = Basin_tcal + DBLE( Tcal(i)*Hru_area(i) )

      ENDDO

      ! Area normalize basin totals
      Basin_snowmelt = Basin_snowmelt*Basin_area_inv
      Basin_pweqv = Basin_pweqv*Basin_area_inv
      Basin_snowevap = Basin_snowevap*Basin_area_inv
      Basin_snowcov = Basin_snowcov*Basin_area_inv
      Basin_pk_precip = Basin_pk_precip*Basin_area_inv
      Basin_snowdepth = Basin_snowdepth*Basin_area_inv
      Basin_tcal = Basin_tcal*Basin_area_inv
      IF ( Glacier_flag==ACTIVE ) THEN
        Basin_glacrb_melt = Basin_glacrb_melt*Basin_area_inv
        Basin_glacrevap = Basin_glacrevap*Basin_area_inv
        Basin_snowicecov = Basin_snowcov
      ENDIF

      IF ( Print_debug==9 ) THEN
        PRINT 9001, Jday, (Net_rain(i), i=1, Nhru)
        PRINT 9001, Jday, (Net_snow(i), i=1, Nhru)
        PRINT 9001, Jday, (Snowmelt(i), i=1, Nhru)
      ENDIF

 9001 FORMAT (I5, 177F6.3)

      END FUNCTION snorun

!***********************************************************************
!      Subroutine to add rain and/or snow to snowpack
!***********************************************************************
      SUBROUTINE ppt_to_pack(Pptmix, Iasw, Tmaxc, Tminc, Tavgc, &
     &           Pkwater_equiv, Net_rain, Pk_def, Pk_temp, Pk_ice, &
     &           Freeh2o, Snowcov_area, Snowmelt, Pk_depth, Pss, Pst, &
     &           Net_snow, Pk_den, Pptmix_nopack, Pk_precip, Tmax_allsnow_c, &
     &           Freeh2o_cap, Den_max, Ihru_gl)
      USE PRMS_CONSTANTS, ONLY: CLOSEZERO, INCH2CM, ACTIVE !, DNEARZERO
      IMPLICIT NONE
! Functions
      REAL, EXTERNAL :: f_to_c
      EXTERNAL :: calin
      INTRINSIC :: ABS, DBLE, SNGL
! Arguments
      INTEGER, INTENT(IN) :: Pptmix, Ihru_gl
      INTEGER, INTENT(INOUT) :: Iasw, Pptmix_nopack
      REAL, INTENT(IN) :: Tmaxc, Tminc, Tavgc, Net_rain, Net_snow
      REAL, INTENT(IN) :: Freeh2o_cap, Tmax_allsnow_c, Den_max
      REAL, INTENT(INOUT) :: Snowmelt, Freeh2o, Pk_precip
      REAL, INTENT(INOUT) :: Pk_def, Pk_ice, Pk_den, Snowcov_area, Pk_temp
      DOUBLE PRECISION, INTENT(INOUT) :: Pkwater_equiv, Pk_depth, Pst, Pss
! Local Variables
      REAL :: train, tsnow, caln, pndz, calpr, calps
!***********************************************************************

      ! The temperature of precipitation will be different if it is mixed or
      ! all rain or snow 2 options below (if-then, else)

      ! If there is any snow, the snow temperature is the average
      ! temperature
      tsnow = Tavgc ! [degrees C]
      ! (1) If precipitation is mixed...
      IF ( Pptmix==ACTIVE ) THEN
        ! If there is any rain, the rain temperature is halfway between the maximum
        ! temperature and the allsnow temperature
        train = (Tmaxc+Tmax_allsnow_c)*0.5 ! [degrees C]

        ! Temperatures will be different, depending on if there is an
        ! existing snowpack or not

        ! If there is a snowpack, snow temperature is halfway between
        ! the minimum daily temperature and maximum temperature for
        ! which all precipitation is snow
        IF ( Pkwater_equiv>0.0D0 ) THEN
          tsnow = (Tminc+Tmax_allsnow_c)*0.5 ! [degrees C]

        ! If there is no existing snowpack, snow temperature is the
        ! average temperature for the day
        ELSEIF ( Pkwater_equiv<0.0D0 ) THEN
!          IF ( Pkwater_equiv<-DNEARZERO ) &
!     &         PRINT *, 'snowpack issue in ppt_to_pack, negative pkwater_equiv', Pkwater_equiv
          Pkwater_equiv = 0.0D0 ! to be sure negative snowpack is ignored
        ENDIF

      ! (2) If precipitation is all snow or all rain...
      ELSE ! on glacier ice goes in here only
        ! If there is any rain, the rain temperature is the average
        ! temperature
        train = Tavgc ! [degrees C]
        ! If average temperature is close to freezing, the rain
        ! temperature is halfway between the maximum daily temperature
        ! and maximum temperature for which all precipitation is snow
        IF ( train<CLOSEZERO ) train = (Tmaxc+Tmax_allsnow_c)*0.5 ! [degrees C]
      ENDIF

      IF ( train<0.0 ) train = 0.0 ! [degrees C] ! train can't be < 0
      IF ( tsnow>0.0 ) tsnow = 0.0 ! [degrees C] ! tsnow can't be > 0

      ! Leavesley comments...
      ! If snowpack already exists, add rain first, then add
      ! snow.  If no antecedent snowpack, rain is already taken care
      ! of, so start snowpack with snow.  This SUBROUTINE assumes
      ! that in a mixed event, the rain will be first and turn to
      ! snow as the temperature drops.

      ! Rain can only add to the snowpack if a previous snowpack
      ! exists, so rain or a mixed event is processed differently
      ! when a snowpack exists
      ! 2 options below (if-then, elseif)

      ! (1) If there is net rain on an existing snowpack...
      IF ( Pkwater_equiv>0.0D0 ) THEN
        IF ( Net_rain>0.0 ) THEN ! on glacier ice goes in here only
          ! Add rain water to pack (rain on snow) and increment the
          ! precipitation on the snowpack by the rain water
          Pkwater_equiv = Pkwater_equiv + DBLE(Net_rain) ! [inches]
          Pk_precip = Pk_precip + Net_rain ! [inches]

          ! Incoming rain water carries heat that must be added to
          ! the snowpack.
          ! This heat could both warm the snowpack and melt snow.
          ! Handling of this heat depends on the current thermal
          ! condition of the snowpack.
          ! 2 options below (if-then, else)

          ! (1.1) If the snowpack is colder than freezing it has a
          ! heat deficit (requires heat to be brought to isothermal
          ! at 0 degC)...
          IF ( Pk_def>0.0 ) THEN
            ! Calculate the number of calories given up per inch of
            ! rain when cooling it from the current rain temperature
            ! to 0 deg C and then freezing it (liquid to solid state
            ! latent heat)
            ! This calculation assumes a volume of an inch of rain
            ! over a square cm of area
            ! 80 cal come from freezing 1 cm3 at 0 C
            ! (latent heat of fusion is 80 cal/cm^3),
            ! 1 cal from cooling 1cm3 for every degree C
            ! (specific heat of water is 1 cal/(cm^3 degC)),
            ! convert from 1 cm depth over 1 square cm to
            ! 1 inch depth over 1 square cm (INCH2CM = 2.54 cm/in)
            caln = (80.0+train)*INCH2CM ! [cal / (in cm^2)]
            ! calculate the amount of rain in inches
            ! (at the current rain temperature)
            ! needed to bring the snowpack to isothermal at 0
            pndz = Pk_def/caln ! [inches]

            ! The effect of rain on the snowpack depends on if there
            ! is not enough, enough, or more than enough heat in the
            ! rain to bring the snowpack to isothermal at 0 degC or not
            ! 3 options below (if-then, elseif, else)

            ! (1.1.1) Exactly enough rain to bring pack to isothermal...
            IF ( ABS(Net_rain-pndz)<CLOSEZERO ) THEN
              ! Heat deficit and temperature of the snowpack go to 0
              Pk_def = 0.0  ! [cal/cm^2]
              Pk_temp = 0.0 ! [degrees C]
              ! In the process of giving up its heat, all the net rain
              ! freezes and becomes pack ice
              Pk_ice = Pk_ice + Net_rain ! [inches]

            ! (1.1.2) Rain not sufficient to bring pack to isothermal...
            ELSEIF ( Net_rain<pndz ) THEN
              ! The snowpack heat deficit decreases by the heat provided
              ! by rain and a new snowpack temperature is calculated
              ! 1.27 is the specific heat of ice (0.5 cal/(cm^3 degC))
              ! times the conversion of cm to inches (2.54 cm/in)
              Pk_def = Pk_def - (caln*Net_rain) ! [cal/(in cm^3)]
              Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0)
              ! All the net rain freezes and becomes pack ice
              Pk_ice = Pk_ice + Net_rain

            ! (1.1.3) Rain in excess of amount required to bring pack
            !         to isothermal...
            ELSE
              ! Heat deficit and temperature of the snowpack go to 0
              Pk_def = 0.0
              Pk_temp = 0.0
              ! The portion of net rain that brings the snowpack to
              ! isothermal freezes
              Pk_ice = Pk_ice + pndz
              ! The rest of the net rain becomes free water in the
              ! snowpack
              ! Note that there cannot be previous Freeh2o because the
              ! snowpack had a heat deficit (all water was ice) before
              ! this condition was reached.
              Freeh2o = Net_rain - pndz
              ! Calculate the excess heat per area added by the portion
              ! of rain that does not bring the snowpack to isothermal
              ! (using specific heat of water)
              calpr = train*(Net_rain-pndz)*INCH2CM ! [cal/cm^2]
              ! Add the new heat to the snow pack
              ! (the heat in this excess rain will melt some of the
              ! pack ice when the water cools to 0 degC)
              CALL calin(calpr, Pkwater_equiv, Pk_def, Pk_temp, &
     &                   Pk_ice, Freeh2o, Snowcov_area, Snowmelt, &
     &                   Pk_depth, Pss, Pst, Iasw, Pk_den, Freeh2o_cap, Den_max, Ihru_gl)
            ENDIF

          ! (1.2) Rain on snowpack that is isothermal
          !       at 0 degC (no heat deficit)...
          ! on glacier ice not active_layer goes in here only, as Pk_def, pndz = 0,
          ELSE
            ! All net rain is added to free water in the snowpack
            Freeh2o = Freeh2o + Net_rain
            ! Calculate the heat per area added by the rain
            ! (using specific heat of water)
            calpr = train*Net_rain*INCH2CM ! [cal/cm^2]
            ! Add the new heat to the snow pack
            ! (the heat in rain will melt some of the pack ice when
            ! the water cools to 0 degC)
            CALL calin(calpr, Pkwater_equiv, Pk_def, Pk_temp, &
     &                 Pk_ice, Freeh2o, Snowcov_area, Snowmelt, &
     &                 Pk_depth, Pss, Pst, Iasw, Pk_den, Freeh2o_cap, Den_max, Ihru_gl)
          ENDIF
        ENDIF

      ! (2) If there is net rain but no snowpack, set flag for a mix
      !     on no snowpack.
      ELSEIF ( Net_rain>0.0 ) THEN
        ! Be careful with the code here.
        ! If this subroutine is called when there is an all-rain day
        ! on no existing snowpack (currently, it will not),
        ! then the flag here will be set inappropriately.
        Pptmix_nopack = ACTIVE ! [flag]
      ENDIF

      ! At this point, the subroutine has handled all conditions
      ! where there is net rain, so if there is net snow
      ! (doesn't matter if there is a pack or not)...
      IF ( Net_snow>0.0 ) THEN
        ! add the new snow to the pack water equivalent, precip, and ice
        Pkwater_equiv = Pkwater_equiv + DBLE(Net_snow)
        Pk_precip = Pk_precip + Net_snow
        Pk_ice = Pk_ice + Net_snow

        ! The temperature of the new snow will determine its effect on
        ! snowpack heat deficit
        ! 2 options below (if-then, else)

        ! (1) if the new snow is at 0 degC...
        IF ( tsnow>=0.0 ) THEN
          ! incoming snow does not change the overall heat content of
          ! the snowpack.
          ! However, the temperature will change, because the total heat
          ! content of the snowpack will be "spread out" among
          ! more snow.  Calculate the snow pack temperature from the
          ! heat deficit, specific heat of snow,
          ! and the new total snowpack water content
          Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]

        ! (2) if the new snow is colder than 0 degC...
        ELSE
          ! calculate the amount of heat the new snow will absorb if
          ! warming it to 0C (negative number).
          ! This is the negative of the heat deficit of the new snow.
          calps = tsnow*Net_snow*1.27 ! [cal/cm^2]

          ! The heat to warm the new snow can come from different
          ! sources depending on the state of the snowpack
          ! 2 options below (if-then, else)

          ! (2.1) if there is free water in the pack
          !       (at least some of it is going to freeze)...
          IF ( Freeh2o>0.0 ) THEN
            CALL caloss(calps, Pkwater_equiv, Pk_def, Pk_temp, Pk_ice, Freeh2o, Ihru_gl)

          ! (2.2) if there is no free water (snow pack has a
          !       heat deficit greater than or equal to 0)...
          ELSE
            ! heat deficit increases because snow is colder than
            ! pack (minus a negative number = plus)
            ! and calculate the new pack temperature
            Pk_def = Pk_def - calps ! [cal/cm^2]
            Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]
          ENDIF
        ENDIF
      ENDIF

      END SUBROUTINE ppt_to_pack

!***********************************************************************
!      Subroutine to compute change in snowpack when a net loss in
!        heat energy has occurred.
!***********************************************************************
      SUBROUTINE caloss(Cal, Pkwater_equiv, Pk_def, Pk_temp, Pk_ice, Freeh2o, Ihru_gl)
      USE PRMS_CONSTANTS, ONLY: CLOSEZERO !, DNEARZERO
      IMPLICIT NONE
! Functions
      INTRINSIC :: SNGL
! Arguments
      INTEGER, INTENT(IN) :: Ihru_gl
      REAL, INTENT(IN) :: Cal
      DOUBLE PRECISION, INTENT(INOUT) :: Pkwater_equiv
      REAL, INTENT(INOUT) :: Pk_def, Pk_ice, Freeh2o, Pk_temp
! Functions
      EXTERNAL :: glacr_states_to_zero
! Local Variables
      REAL :: calnd, dif
!***********************************************************************

      ! Loss of heat is handled differently if there is liquid water in
      ! the snowpack or not
      ! 2 options below (if-then, else)

      ! (1) No free water exists in pack
      IF ( Freeh2o<CLOSEZERO ) THEN
        ! heat deficit increases because snow is colder than pack
        ! (minus a negative number = plus)
        Pk_def = Pk_def - Cal ! [cal/cm^2]

      ! (2) Free water exists in pack
      ELSE
        ! calculate the total amount of heat per area that can be
        ! released by free water freezing
        calnd = Freeh2o*203.2 ! [cal/cm^2]
        ! determine the difference between heat in free water and the
        ! heat that can be absorbed by new snow (without melting)
        ! remember that cal is a negative number
        dif = Cal + calnd ! [cal/cm^2]

        ! The effect of freezing water depends on whether all or only
        ! part of the liquid water freezes
        ! 2 options below (if-then, else)

        ! (2) Only part of free water freezes
        IF ( dif>0.0 ) THEN
          ! the calories absorbed by the new snow freezes some
          ! of the free water
          ! (increase in ice, decrease in free water)
          Pk_ice = Pk_ice + (-Cal/203.2) ! [inches]
          Freeh2o = Freeh2o - (-Cal/203.2) ! [inches]
          RETURN
        ! (1) All free water freezes
        ELSE ! IF ( dif<=0.0 ) THEN
          ! if all the water freezes, then the remaining heat
          ! that can be absorbed by new snow (that which is not
          ! provided by freezing free water) becomes the new pack
          ! heat deficit
          IF ( dif<0.0 ) Pk_def = -dif ! [cal/cm^2]
          ! free pack water becomes ice
          Pk_ice = Pk_ice + Freeh2o ! [inches]
          Freeh2o = 0.0 ! [inches]

        ENDIF
      ENDIF

      ! if there is still a snowpack, calculate the new temperature
      IF ( Pkwater_equiv>0.0D0 ) THEN
        Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0)  ! [degrees C]
      ELSE
        IF ( Pkwater_equiv<0.0D0 ) THEN
!          IF ( Pkwater_equiv<-DNEARZERO ) &
!     &         PRINT *, 'snowpack issue 4, negative pkwater_equiv', Pkwater_equiv
          Pkwater_equiv = 0.0D0
        ENDIF
        ! If on melting glacier ice/firn, Ihru_gl >0, so melted active layer (won't melt infinite ice layer)
        If (Ihru_gl>0) CALL glacr_states_to_zero(Ihru_gl,0)
      ENDIF

      END SUBROUTINE caloss

!***********************************************************************
!      Subroutine to compute changes in snowpack when a net gain in
!        heat energy has occurred.
!***********************************************************************
      SUBROUTINE calin(Cal, Pkwater_equiv, Pk_def, Pk_temp, &
     &                 Pk_ice, Freeh2o, Snowcov_area, Snowmelt, &
     &                 Pk_depth, Pss, Pst, Iasw, Pk_den, Freeh2o_cap, Den_max, Ihru_gl)
      USE PRMS_CONSTANTS, ONLY: DEBUG_less, OFF
      USE PRMS_MODULE, ONLY: Print_debug
      USE PRMS_SNOW, ONLY: Active_glacier
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(INOUT) :: Iasw
      INTEGER, INTENT(IN) :: Ihru_gl
      REAL, INTENT(IN) :: Cal, Freeh2o_cap, Snowcov_area, Den_max
      REAL, INTENT(INOUT) :: Freeh2o
      DOUBLE PRECISION, INTENT(INOUT) :: Pkwater_equiv
      REAL, INTENT(INOUT) :: Pk_def, Pk_temp, Pk_ice, Pk_den, Snowmelt
      DOUBLE PRECISION, INTENT(INOUT) :: Pss, Pst, Pk_depth
! Functions
      INTRINSIC :: SNGL, DBLE
      EXTERNAL :: print_date, glacr_states_to_zero
! Local Variables
      REAL :: dif, pmlt, apmlt, apk_ice, pwcap
      DOUBLE PRECISION :: dif_dble
!***********************************************************************

      ! Calculate the difference between the incoming calories and the
      ! calories needed to bring the pack to isothermal
      ! at 0 (heat deficit)
      dif = Cal - Pk_def ! [cal/cm^2]

      ! The way incoming heat is handled depends on whether there is
      ! not enough, just enough, or more than enough heat to overcome
      ! the heat deficit of the snowpack.
      ! 3 choices below (if-then, elseif, else)

      ! (1) Not enough heat to overcome heat deficit...
      IF ( dif<0.0 ) THEN
        ! Reduce the heat deficit by the amount of incoming calories
        ! and adjust to the new temperature based on new heat deficit
        Pk_def = Pk_def - Cal ! [cal/cm^2]
        Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]

      ! (3) More than enough heat to overcome heat deficit
      !     (melt ice)...
      ELSEIF ( dif>0.0 ) THEN
        ! calculate the potential amount of snowmelt from excess
        ! heat in rain it takes 203.2 calories / (in cm^2) to melt snow
        ! (latent heat of fusion)
        ! convert from 1 cm depth over 1 square cm to
        ! 1 inch depth over 1 square cm 80.0*(INCH2CM = 2.54 cm/in) = 203.2
        pmlt = dif/203.2 ! [inches]
        ! Actual snowmelt can only come from snow covered area, so to
        ! calculate the actual potential snowmelt, the potential
        ! snowmelt from snowcovered area must be re-normalized to
        ! HRU area (rather than snowcover area)
        ! In effect, the potential snowmelt per area is reduced by the
        ! fraction of the watershed that is actually covered by snow
        apmlt = pmlt*Snowcov_area ! [inches]
        ! Set the heat deficit and temperature of the remaining
        ! snowpack to 0
        Pk_def = 0.0 ! [cal/cm^2]
        Pk_temp = 0.0 ! [degrees C]
        ! The only pack ice that is melted is in the snow covered area,
        ! so the pack ice needs to be re-normalized to the snowcovered
        ! area (rather than HRU area)
        ! In effect, the pack ice per area is increased by the fraction
        ! of the watershed that is actually covered by snow
        IF ( Snowcov_area>0.0 ) THEN
          apk_ice = Pk_ice/Snowcov_area ! [inches]
        ELSE
!          PRINT *, 'snowcov_area really small, melt all ice', snowcov_area, ' pmlt:', pmlt, ' dif:', dif, ' pk_ice:', Pk_ice
          apk_ice = 0.0
        ENDIF

        ! If snow is melting, the heat is handled based on whether all
        ! or only part of the pack ice melts
        ! 2 options below (if-then, else)

        ! (3.1) Heat applied to snow covered area is sufficient
        !       to melt all the ice in that snow pack...
        ! if on snow over glacier or active_layer and have excess energy from day over
        !        depth can melt from layer thickness, add depth to that layer
        IF ( Active_glacier>OFF ) THEN
          IF ( pmlt>apk_ice ) THEN
            !fractionate density with snow/active layer melting vs extra ice underneath melting
            Pk_den = Pk_den*SNGL(apk_ice/pmlt) + 0.917*((pmlt-apk_ice)/pmlt)
            apk_ice = pmlt
            Pk_ice =  apmlt
            Pkwater_equiv = apmlt
            Freeh2o = 0.0 ! [inches]
            Iasw = 0
            Pk_def = 0.0   ! [cal / cm^2]
            Pk_temp = 0.0  ! [degreees C]
            Pst = 0.0D0      ! [inches]
          ENDIF
        ENDIF

        IF ( pmlt>apk_ice ) THEN ! will not happen if Active_glacier>OFF because of above
          ! All pack water equivalent becomes meltwater
          Snowmelt = Snowmelt + SNGL( Pkwater_equiv ) ! [inches]
          Pkwater_equiv = 0.0D0 ! [inches]
          Iasw = 0 ! snow area does not change
          ! Set all snowpack states to 0
          ! Snowcov_area = 0.0 ! [fraction of area] ! shouldn't be changed with melt
          Pk_def = 0.0   ! [cal / cm^2]
          Pk_temp = 0.0  ! [degreees C]
          Pk_ice = 0.0   ! [inches]
          Freeh2o = 0.0  ! [inches]
          Pk_depth = 0.0D0 ! [inches]
          Pss = 0.0D0      ! [inches]
          Pst = 0.0D0      ! [inches]
          Pk_den = 0.0     ! [fraction of depth]

        ! (3.2) Heat only melts part of the ice in the snow pack...
        ELSE
          ! Remove actual melt from frozen water and add melt to
          ! free water
          Pk_ice = Pk_ice - apmlt ! [inches]
          Freeh2o = Freeh2o + apmlt ! [inches]
          ! Calculate the capacity of the snowpack to hold free water
          ! according to its current level of frozen water
          pwcap = Freeh2o_cap*Pk_ice ! [inches]
          ! Calculate the amount of free water in excess of the
          ! capacity to hold free water
          dif_dble = DBLE( Freeh2o - pwcap ) ! [inches]
          ! If there is more free water than the snowpack can hold,
          ! then there is going to be melt...
          IF ( dif_dble>0.0D0 ) THEN
            IF ( dif_dble>Pkwater_equiv ) dif_dble = Pkwater_equiv
            ! total packwater decreases by the excess and a new depth
            ! is calculated based on density
            Pkwater_equiv = Pkwater_equiv - dif_dble ! [inches]
            ! free water is at the current capacity
            Freeh2o = pwcap ! [inches]
            IF ( Pk_den>0.0 ) THEN
              Pk_depth = Pkwater_equiv/DBLE(Pk_den) ! [inches]
            ! RAPCOMMENT - added the conditional statement to make
            !   sure there is no division by zero (this can happen
            !   if there is a mixed event on no existing snowpack
            !   because a pack density has not been calculated, yet
            ELSE
            !rsr, this should not happen, remove later
              IF ( Print_debug>DEBUG_less ) THEN
                PRINT *, 'snow density problem', Pk_depth, Pk_den, Pss, Pkwater_equiv
                CALL print_date(1)
              ENDIF
              IF ( Active_glacier==OFF ) Pk_den = Den_max
              Pk_depth = Pkwater_equiv/DBLE(Den_max) ! [inches]
            ENDIF

            ! snowmelt increases by the excess free water
            Snowmelt = Snowmelt + SNGL( dif_dble ) ! [inches]
            ! reset the previous-snowpack-plus-new-snow to the
            ! current pack water equivalent
            Pss = Pkwater_equiv ! [inches]
          ENDIF
        ENDIF
      ! (2) Just enough heat to overcome heat deficit
      ELSE ! IF ( dif==0.0 ) THEN ! rsr 1/27/2016 why not set all snow states to 0 ???
        ! Set temperature and heat deficit to zero
        Pk_temp = 0.0 ! [degrees C]
        Pk_def = 0.0 ! [cal/cm^2]
      ENDIF
      IF ( .NOT.(Pkwater_equiv>0.0D0) ) Pk_den = 0.0
      ! If on melting glacier ice/firn, Ihru_gl >0, so melted active layer (won't melt infinite ice layer)
      IF ( Ihru_gl>0) THEN
        IF ( .NOT.(Pkwater_equiv>0.0D0) ) CALL glacr_states_to_zero(Ihru_gl,0)
      ENDIF

      END SUBROUTINE calin

!***********************************************************************
!      Subroutine to compute snowpack albedo
!***********************************************************************
      SUBROUTINE snalbedo(Newsnow, Iso, Lst, Snsv, Prmx, Pptmix, Albset_rnm, &
     &                    Net_snow, Albset_snm, Albset_rna, Albset_sna, Albedo, &
     &                    Int_alb, Salb, Slst)
      USE PRMS_CONSTANTS, ONLY: OFF
      USE PRMS_SNOW, ONLY: MAXALB, Acum, Amlt
      IMPLICIT NONE
! Functions
      INTRINSIC :: INT
! Arguments
      INTEGER, INTENT(IN) :: Newsnow, Iso, Pptmix
      INTEGER, INTENT(INOUT) :: Int_alb, Lst
      REAL, INTENT(IN) :: Albset_rnm, Albset_snm, Albset_rna, Albset_sna, Prmx, Net_snow
      REAL, INTENT(INOUT) :: Salb, Slst, Snsv
      REAL, INTENT(OUT) :: Albedo
! Local Variables
      INTEGER :: l
!***********************************************************************

      ! The albedo is always reset to a new initial (high) value when
      ! there is new snow above a threshold (parameter).  Albedo
      ! is then a function of the number of days since the last new snow
      ! Intermediate conditions apply when there is new snow
      ! below the threshold to reset the albedo to its highest value.
      ! The curve for albedo change (decreasing) is different for the
      ! snow accumulation season and the snow melt season.
      ! The albedo first depends on if there is no new snow during the
      ! current time step, if there is new snow during accumulation
      ! season, or if there is new snow during melt season.
      ! 3 options below (if-then, elseif, else)

      ! (1) There is no new snow
      IF ( Newsnow==OFF ) THEN
        ! If no new snow, check if there was previous new snow that
        ! was not sufficient to reset the albedo (Lst=1)
        ! Lst can only be greater than 0 during melt season (see below)
        IF ( Lst>0 ) THEN
          ! Slst is the number of days (float) since the last
          ! new snowfall
          ! Set the albedo curve back three days from the number
          ! of days since the previous snowfall
          ! (see Salb assignment below)
          ! (note that "shallow new snow" indicates new snow that
          ! is insufficient to completely reset the albedo curve)
          ! In effect, a shallow new snow sets the albedo curve back
          ! a few days, rather than resetting it entirely.
          Slst = Salb - 3.0 ! [days]
          ! Make sure the number of days since last new snow
          ! isn't less than 1
          IF ( Slst<1.0 ) Slst = 1.0 ! [days]
          ! If not in melt season
          IF ( Iso/=2 ) THEN
            ! Note that this code is unreachable in its current state.
            ! This code is only run during melt season due to the
            ! fact that Lst can only be set to 1 in the melt season.
            ! Therefore, Iso is always going to be equal to 2.
            ! Make sure the maximum point on the albedo curve is 5
            ! In effect, if there is any new snow, the albedo can
            ! only get so low in accumulation season, even if the
            ! new snow is insufficient to reset albedo entirely
            IF ( Slst>5.0 ) Slst = 5.0 ! [days]
          ENDIF
          ! Reset the shallow new snow flag and cumulative shallow
          ! snow variable (see below)
          Lst = 0 ! [flag]
          Snsv = 0.0 ! [inches]
        ENDIF

      ! (2) New snow during the melt season
      ELSEIF ( Iso==2 ) THEN
! RAPCOMMENT - CHANGED TO ISO FROM MSO

        ! If there is too much rain in a precipitation mix,
        ! albedo will not be reset
        ! New snow changes albedo only if the fraction rain
        ! is less than the threshold above which albedo is not reset
        IF ( Prmx<Albset_rnm ) THEN

          ! If the fraction rain doesn't prevent the albedo from
          ! being reset, then how the albedo changes depends on
          ! whether the snow amount is above or below the threshold
          ! for resetting albedo
          ! 2 options below (if-then, else)

          ! (2.1) If there is enough new snow to reset the albedo
          IF ( Net_snow>Albset_snm ) THEN
            ! Reset number of days since last new snow to 0
            Slst = 0.0 ! [days]
            Lst = 0 ! [flag]
            ! Reset the saved new snow to 0
            Snsv = 0.0 ! [inches]

          ! (2.2) If there is not enough new snow this time period
          ! to reset the albedo on its own
          ELSE
            ! Snsv tracks the amount of snow that has fallen as long
            ! as the total new snow is not
            ! enough to reset the albedo.
            Snsv = Snsv + Net_snow ! [inches]

            ! Even if the new snow during this time period is
            ! insufficient to reset the albedo, it may still reset the
            ! albedo if it adds enough to previous shallow snow
            ! accumulation.  The change in Albedo depends on if the
            ! total amount of accumulated shallow snow has become enough
            ! to reset the albedo or not.
            ! 2 options below (if-then, else)

            ! (2.2.1) If accumulated shallow snow is enough to reset
            !         the albedo
            IF ( Snsv>Albset_snm ) THEN
              ! Reset the albedo states.
              Slst = 0.0 ! [days]
              Lst = 0 ! [flag]
              Snsv = 0.0 ! [inches]

            ! (2.2.2) If the accumulated shallow snow is not enough to
            !         reset the albedo curve
            ELSE
              ! Salb records the number of days since the last new snow
              ! that reset albedo
              IF ( Lst==0 ) Salb = Slst ! [days]
              ! Reset the number of days since new snow
              Slst = 0.0 ! [days]
              ! set the flag indicating that there is shallow new snow
              ! (i.e. not enough new snow to reset albedo)
              Lst = 1 ! [flag]
            ENDIF
          ENDIF
        ENDIF

      ! (3) New snow during the accumulation season
      ELSE

        ! The change in albedo depends on if the precipitation is a mix,
        ! if the rain is above a threshold,  or if the snow is above
        ! a threshold.
        ! 4 options below (if-then, elseif, elseif, else)

        ! (3.1) If it is not a mixed event...
        IF ( Pptmix==OFF ) THEN
          ! During the accumulation season, the threshold for resetting
          ! the albedo does not apply if there is a snow-only event.
          ! Therefore, no matter how little snow there is, it will
          ! always reset the albedo curve the the maximum, if it
          ! occurs during the accumulation season.
          ! reset the time since last snow to 0
          Slst = 0.0 ! [days]
          ! there is no new shallow snow
          Lst = 0 ! [flag]

        ! (3.2) If it is a mixed event and the fraction rain is above
        !       the threshold above which albedo is not reset...
        ELSEIF ( Prmx>=Albset_rna ) THEN
          ! there is no new shallow snow
          Lst = 0 ! [flag]
          ! albedo continues to decrease on the curve

        ! (3.3) If it is a mixed event and there is enough new snow
        !       to reset albedo...
        ELSEIF ( Net_snow>=Albset_sna ) THEN
          ! reset the albedo
          Slst = 0.0 ! [days]
          ! there is no new shallow snow
          Lst = 0 ! [flag]

        ! (3.4) If it is a mixed event and the new snow was not
        !       enough to reset the albedo...
        ELSE
          ! set the albedo curve back 3 days (increasing the albedo)
          Slst = Slst - 3.0 ! [days]
          ! Make sure the number of days since last new snow is not
          ! less than 0
          IF ( Slst<0.0 ) Slst = 0.0 ! [days]
          ! Make sure the number of days since last new snow is not
          ! greater than 5
          ! In effect, if there is any new snow, the albedo can
          ! only get so low in accumulation season, even if the
          ! new snow is insufficient to reset albedo entirely
          IF ( Slst>5.0 ) Slst = 5.0 ! [days]
          Lst = 0 ! [flag]
        ENDIF
        Snsv = 0.0 ! [inches]
      ENDIF
      ! At this point, the subroutine knows where on the curve the
      ! albedo should be based on current conditions and the
      ! new snow (determined by value of Slst variable)

      ! Get the integer value for days (or effective days)
      ! since last snowfall
      l = INT(Slst+0.5) ! [days]

      ! Increment the state variable for days since the
      ! last snowfall
      Slst = Slst + 1.0 ! [days]

      !******Compute albedo
      ! Albedo will only be different from the max (default value)
      ! if it has been more than 0 days since the last new snow
      ! capable of resetting the albedo.  If albedo is at the
      ! maximum, the maximum is different for accumulation and
      ! melt season.
      ! 3 options below (if-then, elseif, else)

      ! (1) It has been more than 0 days since the last new snow
      IF ( l>0 ) THEN

        ! Albedo depends on whether it is currently on the
        ! accumulation season curve or on the melt season curve.
        ! 3 options below (if-then, elseif, else)

        ! (1.1) Currently using the melt season curve
        !       (Old snow - Spring melt period)...
        IF ( Int_alb==2 ) THEN
          ! Don't go past the last possible albedo value
          IF ( l>MAXALB ) l = MAXALB ! [days]
          ! Get the albedo number from the melt season curve
          Albedo = Amlt(l) ! [fraction of radiation]

        ! (1.2) Currently using the accumulation season curve
        !       (Old snow - Winter accumulation period)...
        ! and not past the maximum curve index
        ELSEIF ( l<=MAXALB ) THEN
          ! Get the albedo number from the accumulation season curve
          Albedo = Acum(l) ! [fraction of radiation]

        ! (1.3) Currently using the accumulation season curve and
        !       past the maximum curve index...
        ELSE
          ! start using the the MELT season curve at 12 days
          ! previous to the current number of days since the last
          ! new snow
          l = l - 12 ! [days]
          ! keep using the melt season curve until its minimum
          ! value (maximum index) is reached or until there is new snow
          IF ( l>MAXALB ) l = MAXALB ! [days]
          ! get the albedo value from the melt season curve
          Albedo = Amlt(l) ! [fraction of radiation]
        ENDIF

      ! (2) New snow has reset the albedo and it is melt season
      ELSEIF ( Iso==2 ) THEN
! RAPCOMMENT - CHANGED TO ISO FROM MSO
        ! Set albedo to initial value during melt season
        Albedo = 0.72 ! [fraction of radiation] value Rob suggested
!       Albedo = 0.81 ! [fraction of radiation] original value
        ! Int_alb is a flag to indicate use of the melt season curve (2)
        ! or accumulation season curve (1)
        ! Set flag to indicate melt season curve
        Int_alb = 2 ! [flag]

      ! (3) New snow has reset the albedo and it is accumulation season
      ELSE
        ! Set albedo to initial value during accumulation season
        Albedo = 0.91 ! [fraction of radiation]
        ! Set flag to indicate accumulation season curve
        Int_alb = 1 ! [flag]
      ENDIF

      END SUBROUTINE snalbedo

!***********************************************************************
!      Subroutine to compute energy balance of snowpack
!        1st call is for night period, 2nd call for day period
!***********************************************************************
      SUBROUTINE snowbal(Niteda, Tstorm_mo, Iasw, Temp, Esv, Hru_ppt, &
     &           Trd, Emis_noppt, Canopy_covden, Cec, Pkwater_equiv, &
     &           Pk_def, Pk_temp, Pk_ice, Freeh2o, Snowcov_area, &
     &           Snowmelt, Pk_depth, Pss, Pst, Pk_den, Cst, Cal, Sw, Freeh2o_cap, Den_max, Ihru_gl)
      USE PRMS_CONSTANTS, ONLY: CLOSEZERO
      IMPLICIT NONE
! Functions
      INTRINSIC :: SNGL
      EXTERNAL :: calin, caloss
! Arguments
      INTEGER, INTENT(IN) :: Niteda, Tstorm_mo, Ihru_gl
      INTEGER, INTENT(INOUT) :: Iasw
      REAL, INTENT(IN) :: Temp, Esv, Trd, Cec, Cst, Canopy_covden
      REAL, INTENT(IN) :: Emis_noppt, Sw, Freeh2o_cap
      REAL, INTENT(IN) :: Hru_ppt, Snowcov_area, Den_max
      DOUBLE PRECISION, INTENT(INOUT) :: Pst, Pss
      REAL, INTENT(OUT) :: Cal
      REAL, INTENT(INOUT) :: Pk_den, Pk_def, Pk_temp, Pk_ice
      REAL, INTENT(INOUT) :: Freeh2o, Snowmelt
      DOUBLE PRECISION, INTENT(INOUT) :: Pkwater_equiv, Pk_depth
! Local Variables
      REAL :: air, ts, emis, sno, sky, can, cecsub, qcond, pk_defsub, pkt, pks
      REAL, PARAMETER :: ONETHIRD = 1.0/3.0
!***********************************************************************
      ! Calculate the potential long wave energy from air based on
      ! temperature (assuming perfect black-body emission)
      ! Stefan Boltzmann/2 = (11.71E-8)/2 = 0.585E-7 because add for day and night
      air = 0.585E-7*((Temp+273.16)**4.0) ! [cal/cm^2] or [Langleys]
      ! set emissivity, which is the fraction of perfect black-body
      ! emission that is actually applied
      emis = Esv ! [fraction of radiation]

      ! The snowpack surface temperature and long-wave radiation
      ! FROM the snowpack depend on the air temperature (effectively,
      ! snowpack temperature cannot be larger than 0 degC)
      ! 2 options below (if-then, else)

      ! (1) If the temperature is below freezing, surface snow
      !     temperature and long wave energy are determined
      !     by temperature...
      IF ( Temp<0.0 ) THEN
        ts = Temp ! [degrees C]
        sno = air ! [cal/cm^2] or [Langleys]

      ! (2) If the temperature is at or above freezing, snow
      !     temperature and long wave energy are set to values
      !     corresponding to a temperature of 0 degC...
      ELSE
        ts = 0.0 ! [degrees C]
        sno = 325.7 ! [cal/cm^2] or [Langleys]
      ENDIF

      ! If precipitation over the time period was due to
      ! convective thunderstorms, then the emissivity should be reset
      IF ( Hru_ppt>0.0 ) THEN
        IF ( Tstorm_mo==1 ) THEN

          ! The emissivity of air depends on if it is day or night
          ! and the fraction of measured short wave radiation to
          ! potential short wave radiation is used as a surrogate
          ! to the duration of the convective storms
          ! 2 options below (if-then, else)

          ! (1) Night
          IF ( Niteda==1 ) THEN
            ! set the default emissivity
            emis = 0.85 ! [fraction of radiation]
            ! if measured radiation is greater than 1/3 potential
            ! radiation through the time period, then the emissivity
            ! is set to the "no precipitation" value
            IF ( Trd>ONETHIRD ) emis = Emis_noppt ![fraction of radiation]

          ! (2) Day
          ELSE
            ! if measured radiation is greater than 1/3 potential
            ! radiation but less than 1/2, then the emissivity is
            ! interpolated between 1.0 and 0.85
            ! if measured radiation is greater than 1/2 potential
            ! radiation, then the emissivity is interpolated between
            ! 0.85 and 0.75
            IF ( Trd>ONETHIRD ) emis = 1.29 - (0.882*Trd)
                                              ! [fraction of radiation]
            IF ( Trd>=0.5 ) emis = 0.95 - (0.2*Trd)
                                              ! [fraction of radiation]
          ENDIF
        ENDIF
      ENDIF

      ! Calculate the net incoming long wave radiation coming from the
      ! sky or canopy in the uncovered or covered portions of the
      ! snowpack, respectively.
      ! Note that the canopy is assumed to be a perfect blackbody
      ! (emissivity = 1) and the air has emissivity as determined
      ! from previous calculations
      sky = (1.0-Canopy_covden)*((emis*air)-sno) ! [cal/cm^2] or [Langleys]
      can = Canopy_covden*(air-sno) ! [cal/cm^2] or [Langleys]
!RAPCOMMENT  - CHECK THE INTERECEPT MODULE FOR CHANGE.  What if the land
! cover is grass? Is this automatically covered by canopy_covden being zero
! if the cover type is grass?

      ! If air temperature is above 0 degC then set the energy from
      ! condensation and convection, otherwise there is
      ! no energy from convection or condensation
      cecsub = 0.0 ! [cal/cm^2] or [Langleys]
      IF ( Temp>0.0 ) THEN
        IF ( Hru_ppt>0.0 ) cecsub = Cec*Temp ! [cal/cm^2]
                                                     ! or [Langleys]
      ENDIF

      ! Total energy potentially available from atmosphere: longwave,
      ! shortwave, and condensation/convection
      Cal = sky + can + cecsub + Sw ! [cal/cm^2] or [Langleys]

      ! If the surface temperature of the snow is 0 degC, and there
      ! is net incoming energy, then energy conduction has to be from
      ! the surface into the snowpack.
      ! Therefore, the energy from the atmosphere is applied to the
      ! snowpack and subroutine terminates
      IF ( ts>=0.0 ) THEN
        IF ( Cal>0.0 ) THEN
          CALL calin(Cal, Pkwater_equiv, Pk_def, Pk_temp, &
     &               Pk_ice, Freeh2o, Snowcov_area, Snowmelt, &
     &               Pk_depth, Pss, Pst, Iasw, Pk_den, Freeh2o_cap, Den_max, Ihru_gl)
          RETURN
        ENDIF
      ENDIF

      ! If the program gets to this point, then either the surface
      ! temperature is less than 0 degC, or the total energy from the
      ! atmosphere is not providing energy to the snowpack

      ! Because the temperature of the surface of the snowpack is
      ! assumed to be  controlled by air temperature, there is a
      ! potential heat flux due to conduction between the deeper
      ! snowpack and its surface.
      ! Calculate conductive heat flux as a function of the
      ! temperature gradient then set new snowpack conditions
      ! depending on the direction of heat flow
      qcond = Cst*(ts-Pk_temp) ! [cal/cm^2] or [Langleys]
!RAPCOMMENT - The original equation in the paper implies that the
! this equation should be relative to the temperature gradient
! in degF, not degC (Anderson 1968).  Which is correct?

      ! The energy flow depends on the direction of conduction and the
      ! temperature of the surface of the snowpack.  The total energy
      ! from the atmosphere can only penetrate into the snow pack if
      ! the temperature gradient allows conduction from the surface
      ! into the snowpack.
      ! 4 options below (if-then, elseif, elseif, else)

      ! (1) Heat is conducted from the snowpack to the surface
      !     (atmospheric energy is NOT applied to snowpack)...
      IF ( qcond<0.0 ) THEN
        ! If the temperature of the snowpack is below 0 degC,
        ! add to the heat deficit.  Otherwise, remove heat
        ! from the 0 degC isothermal snow pack.
        IF ( Pk_temp<0.0 ) THEN
          ! increase the heat deficit (minus a negative)
          ! and adjust temperature
          Pk_def = Pk_def - qcond ! [cal/cm^2] or [Langleys]
          Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]
        ELSE
          ! remove heat from the snowpack
          CALL caloss(qcond, Pkwater_equiv, Pk_def, Pk_temp, Pk_ice, Freeh2o, Ihru_gl)
        ENDIF
      ! Even though Cal is not applied to the snowpack under this
      ! condition, it maintains its value and the referencing code
      ! uses it to calculate the total energy balance of the snowpack.
      ! Right, now, Cal isn't used for anything outside this subroutine,
      ! but care should be taken if it is.

      ! (2)  There is no heat conduction, qcond = 0.0
      ELSEIF ( qcond<CLOSEZERO ) THEN

        ! if the pack temperature is isothermal at 0 degC, then apply
        ! any incoming radiation, condensation (latent heat),
        ! and convection heat to the snowpack
        IF ( Pk_temp>=0.0 ) THEN
          ! It does not appear that the interior of the following if
          ! statement is reachable in its current form, because if these
          ! conditions are true, then the code for surface temperature=0
          ! and cal=positive number would have run and the subroutine
          ! will have terminated
          IF ( Cal>0.0 ) CALL calin(Cal, Pkwater_equiv, Pk_def, Pk_temp, &
     &                              Pk_ice, Freeh2o, Snowcov_area, &
     &                              Snowmelt, Pk_depth, Pss, Pst, Iasw, Pk_den, &
     &                              Freeh2o_cap, Den_max, Ihru_gl)
        ENDIF

      ! (3) conduction is from the surface to the snowpack and the
      !     surface temperature is 0 degrees C...
      ELSEIF ( ts>=0.0 ) THEN
        ! note that Cal must be <= 0 for this condition to apply.
        ! Otherwise, the program wouldn't have gotten to this point.

        ! determine if the conductive heat is enough to overcome the
        ! current heat deficit
        pk_defsub = Pk_def - qcond
        IF ( pk_defsub<0.0 ) THEN
          ! deficit is overcome and snowpack becomes
          ! isothermal at 0 degC
          Pk_def = 0.0  ! [cal/cm^2] or [Langleys]
          Pk_temp = 0.0 ! [degrees C]
        ELSE
          ! deficit is decreased by conducted heat and temperature
          ! is recalculated
          Pk_def = pk_defsub ! [cal/cm^2] or [Langleys]
          Pk_temp = -pk_defsub/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]
        ENDIF

      ! (4) conduction is from the surface to the snowpack and the
      !     surface temperature is less than 0 degrees C...
      ELSE
        ! calculate the pack deficit if the snowpack was all at the
        ! surface temperature, then calculate how many calories to
        ! shift the pack to that deficit (pks will be a positive
        ! number because the conduction direction is from the surface
        ! into the snowpack)
        pkt = -ts*SNGL(Pkwater_equiv*1.27D0) ! [cal/cm^2] or [Langleys]
        pks = Pk_def - pkt ! [cal/cm^2] or [Langleys]
        ! determine if the conducted heat is enough to shift the
        ! pack to the deficit relative to the surface temperature
        pk_defsub = pks - qcond ! [cal/cm^2] or [Langleys]

        ! The effect of incoming conducted heat depends on whether
        ! it is enough to bring the snowpack to the same temperature
        ! as the surface or not
        ! 2 options below (if-then, else)

        ! (4.1) There is enough conducted heat to bring the deep
        !       snowpack to the surface temperature...
        IF ( pk_defsub<0.0 ) THEN
          ! there is enough conduction to change to the new pack deficit
          Pk_def = pkt ! [cal/cm^2] or [Langleys]
          Pk_temp = ts ! [degrees C]

        ! (4.2) There is not enough conducted heat to bring the deep
        !       snowpack to the surface temperature...
        ELSE
          ! the pack deficit doesn't make it all the way to the surface
          ! deficit, but is decreased relative to the conducted heat
          ! note that the next statement is equivalent to
          ! Pk_def = Pk_def - qcond
          Pk_def = pk_defsub + pkt ! [cal/cm^2] or [Langleys]
          Pk_temp = -Pk_def/SNGL(Pkwater_equiv*1.27D0) ! [degrees C]
        ENDIF
      ENDIF

      END SUBROUTINE snowbal

!***********************************************************************
!      Subroutine to compute evaporation from snowpack
!***********************************************************************
      SUBROUTINE snowevap(Potet_sublim, Potet, Snowcov_area, Snow_evap, &
     &                    Pkwater_equiv, Pk_ice, Pk_def, Freeh2o, Pk_temp, Hru_intcpevap)
      USE PRMS_CONSTANTS, ONLY: CLOSEZERO, DNEARZERO, DEBUG_less, OFF
      USE PRMS_MODULE, ONLY: Print_debug
      USE PRMS_SNOW, ONLY: Active_glacier
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE, SNGL
! Arguments
      REAL, INTENT(IN) :: Potet_sublim, Potet, Snowcov_area, Hru_intcpevap
      REAL, INTENT(INOUT) :: Pk_ice, Pk_def, Pk_temp
      DOUBLE PRECISION, INTENT(INOUT) :: Pkwater_equiv
      REAL, INTENT(OUT) :: Snow_evap, Freeh2o
! Local Variables
      REAL :: avail_et, cal, ez
!***********************************************************************
      ! the amount of evaporation affecting the snowpack is the
      ! total evaporation potential minus the evaporation from
      ! the interception storage
      ez = Potet_sublim*Potet*Snowcov_area - Hru_intcpevap ! [inches]

      ! The effects of evaporation depend on whether there is any
      ! potential for evaporation, and if the potential evapotation
      ! is enough to completely deplete the snow pack or not
      ! 3 options below (if-then, elseif, else)

      ! (1) There is no potential for evaporation...
      ! if on snow over glacier or active_layer and have excess energy from day over
      !        depth can evap from layer thickness, add depth to that layer
      IF ( Active_glacier>OFF ) THEN
        IF ( ez>Pkwater_equiv ) Pkwater_equiv = DBLE(ez)
      ENDIF
      IF ( ez<CLOSEZERO ) THEN
        Snow_evap = 0.0 ! [inches]

      ! (2) Enough potential evaporation to entirely deplete
      !     the snowpack...
      ELSEIF ( ez>=Pkwater_equiv ) THEN
        ! Set the evaporation to the pack water equivalent and set
        ! all snowpack variables to no-snowpack values
        Snow_evap = SNGL(Pkwater_equiv) ! [inches]
        Pkwater_equiv = 0.0D0 ! [inches]
        Pk_ice = 0.0 ! [inches]
        Pk_def = 0.0 ! [cal/cm^2]
        Freeh2o = 0.0 ! [inches]
        Pk_temp = 0.0 ! [degrees C]

      ! (3) Potential evaporation only partially depletes snowpack
      ELSE
        ! Evaporation depletes the amount of ice in the snowpack
        ! (sublimation)
        Pk_ice = Pk_ice - ez

        ! Change the pack conditions according to whether there is
        ! any ice left in the snowpack
        IF ( Pk_ice<0.0 ) THEN
!RAPCOMMENT - CHANGED TO CHECK FOR NEGATIVE PACK ICE
          ! If all pack ice is removed, then there cannot be a
          ! heat deficit
          ! JLM mass balance fix, if we take pkwater_equiv - ez below
          ! JLM we need to take ez from both pk_ice and freeh2o here
          freeh2o = freeh2o + pk_ice
          Pk_ice = 0.0
          Pk_def = 0.0
          Pk_temp = 0.0
        ELSE
          ! Calculate the amount of heat deficit that is removed
          ! by the sublimating ice
          ! Note that this only changes the heat deficit if the
          ! pack temperature is less than 0degC
          cal = Pk_temp*ez*1.27
          Pk_def = Pk_def + cal
        ENDIF
        ! Remove the evaporated water from the pack water equivalent
        Pkwater_equiv = Pkwater_equiv - ez
        Snow_evap = ez
      ENDIF

      IF ( Snow_evap<0.0 ) THEN
        Pkwater_equiv = Pkwater_equiv - DBLE(Snow_evap)

        IF ( Pkwater_equiv<0.0D0 ) THEN
          IF ( Print_debug>DEBUG_less ) THEN
            IF ( Pkwater_equiv<-DNEARZERO ) &
     &           PRINT *, 'snowpack issue, negative pkwater_equiv in snowevap', Pkwater_equiv
            Pkwater_equiv = 0.0D0  ! JLM: this is INSIDE a debug statement? will change the answers
            !  is this in the originial source
          ENDIF
        ENDIF

        Snow_evap = 0.0
      ENDIF

      avail_et = Potet - Hru_intcpevap - Snow_evap
      IF ( avail_et<0.0 ) THEN
!        PRINT *, 'snow evap', snow_evap, avail_et, pkwater_equiv
        Snow_evap = Snow_evap + avail_et
        Pkwater_equiv = Pkwater_equiv - DBLE(avail_et)
        IF ( Snow_evap<0.0 ) THEN
          Pkwater_equiv = Pkwater_equiv - Snow_evap
          IF ( Pkwater_equiv<0.0D0 ) THEN
            IF ( Print_debug>DEBUG_less ) THEN
              IF ( Pkwater_equiv<-DNEARZERO ) &
     &           PRINT *, 'snowpack issue 2, negative pkwater_equiv in snowevap', Pkwater_equiv
            ENDIF
            Pkwater_equiv = 0.0D0 ! to be sure negative snowpack is ignored
          ENDIF
          Snow_evap = 0.0
        ENDIF
      ENDIF

      END SUBROUTINE snowevap

!***********************************************************************
!      Subroutine to compute snow-covered area
!***********************************************************************
      SUBROUTINE snowcov(Iasw, Newsnow, Snowcov_area, Snarea_curve, &
     &                   Pkwater_equiv, Pst, Snarea_thresh, Net_snow, &
     &                   Scrv, Pksv, Snowcov_areasv, Ai, Frac_swe)
      USE PRMS_CONSTANTS, ONLY: DNEARZERO
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Newsnow
      INTEGER, INTENT(INOUT) :: Iasw
      REAL, INTENT(IN) :: Snarea_thresh, Net_snow, Snarea_curve(11)
      DOUBLE PRECISION, INTENT(IN) :: Pkwater_equiv
      REAL, INTENT(INOUT) :: Snowcov_area
      DOUBLE PRECISION, INTENT(OUT) :: Ai
      REAL, INTENT(INOUT) :: Snowcov_areasv
      DOUBLE PRECISION, INTENT(INOUT) :: Pst, Scrv, Pksv
      REAL, INTENT(OUT) :: Frac_swe
! Functions
      INTRINSIC :: DBLE, SNGL, MIN
      EXTERNAL :: sca_deplcrv
! Local Variables
      REAL :: snowcov_area_ante
      DOUBLE PRECISION :: fracy, difx, dify
!***********************************************************************
      snowcov_area_ante = Snowcov_area
      ! Reset snowcover area to the maximum
      Snowcov_area = Snarea_curve(11) ! [fraction of area]

      ! Track the maximum pack water equivalent for the current
      ! snow pack
      IF ( Pkwater_equiv>Pst ) Pst = Pkwater_equiv ! [inches]

      ! Set ai to the maximum packwater equivalent, but no higher than
      ! the threshold for complete snow cover
      Ai = Pst ! [inches]
      IF ( Ai>Snarea_thresh ) Ai = DBLE( Snarea_thresh ) ! [inches]

      ! calculate the ratio of the current packwater equivalent to
      ! the maximum packwater equivalent for the given snowpack
      IF ( Ai>DNEARZERO ) THEN
        Frac_swe = SNGL( Pkwater_equiv/Ai ) ! [fraction]
        Frac_swe = MIN( 1.0, Frac_swe )
      ELSE
!        print *, ai, snarea_thresh
        Frac_swe = 0.0
      ENDIF

      ! There are 3 potential conditions for the snow area curve:
      ! A. snow is accumulating and the pack is currently at its
      !    maximum level
      ! B. snow is depleting and the area is determined by the
      !    snow area curve
      ! C. new snow has occured on a depleting pack, temporarily
      !    resetting to 100% cover.
      ! For case (C), the snow covered area is linearly interpolated
      ! between 100% and the snow covered area before the new snow.
      ! In general, 1/4 of the new snow has to melt before the snow
      ! covered area goes below 100%, and then the remaining 3/4 has
      ! to melt to return to the previous snow covered area.

      ! First, the code decides whether snow is accumulating (A)
      ! or not (B/C).
      ! 2 options below (if-then, else)

      ! (1) The pack water equivalent is at the maximum
      IF ( Pkwater_equiv>=Ai ) THEN
        ! Stay on the snow area curve (it will be at the maximum
        ! because the pack water equivalent is equal to ai
        ! and it can't be higher)
        Iasw = 0

      ! (2) The pack water equivalent is less than the maximum
      ELSE

        ! If the snowpack isn't accumulating to a new maximum,
        ! it is either on the curve (condition B above) or being
        ! interpolated between the previous place on the curve and
        ! 100% (condition C above)
        ! 2 options below (if-then, elseif)

        ! (2.1) There was new snow...
        IF ( Newsnow/=0 ) THEN

          ! New snow will always reset the snow cover to 100%.
          ! However, different states changes depending  on whether
          ! the previous snow area condition was on the curve or
          ! being interpolated between the curve and 100%
          ! 2 options below (if-then, else)

          ! (2.1.1) The snow area is being interpolated between 100%
          !         and a previous location on the curve...
          IF ( Iasw>0 ) THEN
            ! The location on the interpolated line is based on how
            ! much of the new snow has melted.  Because the first 1/4
            ! of the new snow doesn't matter, it has to keep track of
            ! the current snow pack plus 3/4 of the new snow.
            Scrv = Scrv + (0.75D0*DBLE(Net_snow)) ! [inches]
            ! Scrv = Pkwater_equiv - (0.25D0*DBLE(Net_snow))) ! [inches]
!RAPCOMMENT - CHANGED TO INCREMENT THE SCRV VALUE IF ALREADY
!             INTERPOLATING BETWEEN CURVE AND 100%

          ! (2.1.2) The current snow area is on the curve...
          ELSE
            ! If switching from the snow area curve to interpolation
            ! between the curve and 100%, the current state of the snow
            ! pack has to be saved so that the interpolation can
            ! continue until back to the original conditions.
            ! First, set the flag to indicate interpolation between 100%
            ! and the previous area should be done
            Iasw = 1 ! [flag]
            ! Save the current snow covered area
            ! (before the new net snow)
            Snowcov_areasv = snowcov_area_ante ! [inches]
            ! Save the current pack water equivalent
            ! (before the new net snow)
            Pksv = Pkwater_equiv - DBLE( Net_snow ) ! [inches]
            ! The location on the interpolated line is based on how much
            ! of the new snow has melted.  Because the first 1/4
            ! of the new snow doesn't matter, it has to keep track of
            ! the current snow pack plus 3/4 of the new snow.
            Scrv = Pkwater_equiv - (0.25D0*DBLE(Net_snow)) ! [inches]
          ENDIF
          ! The subroutine terminates here because the snow covered area
          ! always starts at 100% if there is any new snow (no need to
          ! reset it from the maximum value set at the beginning of the
          ! subroutine).
          RETURN

        ! (2.2) There was no new snow, but the snow covered area is
        !       currently being interpolated between 100%
        !       from a previous new snow and the snow covered area
        !       before that previous new snow...
        ELSEIF ( Iasw/=0 ) THEN
          ! If the first 1/4 of the previous new snow has not melted,
          ! yet, then the snow covered area is still
          ! 100% and the subroutine can terminate.
          IF ( Pkwater_equiv>Scrv ) RETURN

          ! At this point, the program is almost sure it is
          ! interpolating between the previous snow covered area and
          ! 100%, but it is possible that enough snow has melted to
          ! return to the snow covered area curve instead.
          ! 2 options below (if-then, else)

          ! (2.2.1) The snow pack still has a larger water equivalent
          !         than before the previous new snow.  I.e., new snow
          !         has not melted back to original area...
          IF ( Pkwater_equiv>=Pksv ) THEN
            ! Do the interpolation between 100% and the snow covered
            ! area before the previous new snow.

            ! Calculate the difference between the maximum snow
            ! covered area (remember that Snowcov_area is always
            ! set to the maximum value at this point) and the snow
            ! covered area before the last new snow.
            difx = DBLE( Snowcov_area - Snowcov_areasv )
            ! Calculate the difference between the water equivalent
            ! before the last new snow and the previous water
            ! equivalent plus 3/4 of the last new snow.
            ! In effect, get the value of 3/4 of the previous
            ! new snow.
            dify = Scrv - Pksv ! [inches]                       !gl1098

            ! If 3/4 of the previous new snow is significantly
            ! different from zero, then calculate the ratio of the
            ! unmelted amount of previous new snow in the snow pack
            ! to the value of 3/4 of previous new snow.
            ! In effect, this is the fraction of the previous new snow
            ! that determines the current interpolation
            ! of snow covered area.
            fracy = 0.0D0 ! [fraction]                             !gl1098
            IF ( dify>0.0D0 ) fracy = (Pkwater_equiv-Pksv)/dify
                                                           ! [fraction]
            ! Linearly interpolate the new snow covered area.
            Snowcov_area = Snowcov_areasv + SNGL(fracy*difx)
                                                   ! [fraction of area]
            ! Terminate the subroutine
            RETURN

          ! (2.2.2) The snow pack has returned to the snow water
          ! equivalent before the previous new snow. I.e. back to
          ! original area before new snow.
          ELSE
            ! Reset the flag to use the snow area curve
            Iasw = 0 ! [flag]
          ENDIF

        ENDIF

        ! If this subroutine is still running at this point, then the
        ! program knows that the snow covered area needs to be
        ! adjusted according to the snow covered area curve.  So at
        ! this point it must interpolate between points on the snow
        ! covered area curve (not the same as interpolating between
        ! 100% and the previous spot on the snow area depletion curve).

        CALL sca_deplcrv(Snowcov_area, Snarea_curve, Frac_swe)

      ENDIF

      END SUBROUTINE snowcov

!***********************************************************************
!     Interpolate along snow covered area depletion curve
!***********************************************************************
      SUBROUTINE sca_deplcrv(Snowcov_area, Snarea_curve, Frac_swe)
      IMPLICIT NONE
! Functions
      INTRINSIC :: INT, FLOAT
! Arguments
      REAL, INTENT(OUT) :: Snowcov_area
      REAL, INTENT(IN) :: Snarea_curve(11), Frac_swe
! Local Variables
      INTEGER :: idx, jdx
      REAL :: af, dify, difx
!***********************************************************************
      IF ( Frac_swe>1.0 ) THEN
        Snowcov_area = Snarea_curve(11)
      ELSE

        ! get the indices (as integers) of the depletion curve that
        ! bracket the given Frac_swe (next highest and next lowest)
        idx = INT( 10.0*(Frac_swe+0.2) ) ! [index]
        jdx = idx - 1 ! [index]
        IF ( idx>11 ) idx = 11
        ! calculate the fraction of the distance (from the next lowest)
        ! the given Frac_swe is between the next highest and lowest
        ! curve values
        af = FLOAT( jdx-1 )
        dify = (Frac_swe*10.0) - af ! [fraction]
        ! calculate the difference in snow covered area represented
        ! by next highest and lowest curve values
        difx = Snarea_curve(idx) - Snarea_curve(jdx)
        ! linearly interpolate a snow covered area between those
        ! represented by the next highest and lowest curve values
        Snowcov_area = Snarea_curve(jdx) + dify*difx
      ENDIF
      END SUBROUTINE sca_deplcrv

!***********************************************************************
!     Set all glacier states to 0
!***********************************************************************
      SUBROUTINE glacr_states_to_zero(Ihru, active_layer_present)
      USE PRMS_SNOW, ONLY: Glacr_freeh2o_cap, Glacr_freeh2o_capm, Glacr_pk_def, Glacr_pk_depth, &
     &    Glacr_layer, Glacr_pk_temp, Ann_tempc, Glacr_pkwater_equiv, Glacr_pk_den, &
     &    Glacr_pk_ice, Glacr_pkwater_ante, Glacr_freeh2o, Glacr_pss, Glacr_pk_den
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Ihru, active_layer_present
! Functions
      INTRINSIC :: ATAN, SNGL
! Local Variables
      REAL :: reduce
!***********************************************************************
      IF ( Glacr_layer(Ihru)==0.0 .OR. active_layer_present==0) THEN
        Glacr_pk_depth(Ihru) = 1.0D5
        Glacr_pk_temp(Ihru) = 0.0
        Glacr_pk_def(Ihru) = 0.0
        Glacr_freeh2o_capm(Ihru) = 0.0
        reduce = 1.0
      ElSE
        Glacr_pk_depth(Ihru) = DBLE(Glacr_layer(Ihru))
        Glacr_pk_temp(Ihru) = Ann_tempc(Ihru) !start at average last year temp like Oerlemans 1992
        IF ( Glacr_pk_temp(Ihru) > 0.0) Glacr_pk_temp(Ihru) = 0.0
        Glacr_freeh2o_capm(Ihru) = Glacr_freeh2o_cap(Ihru)
        reduce = 0.8 !if start Glacr_pk_ice too close to Glacr_pk_depth can't grow with energy loss to free water gain
      ENDIF
      Glacr_pk_den(Ihru) = 0.917
      Glacr_pkwater_equiv(Ihru) = Glacr_pk_den(Ihru)*Glacr_pk_depth(Ihru)
      Glacr_pkwater_ante(Ihru) = Glacr_pkwater_equiv(Ihru)
      Glacr_pk_ice(Ihru) = reduce*SNGL(Glacr_pkwater_equiv(Ihru)-Glacr_freeh2o(Ihru))/0.9340 !density of pure ice
      Glacr_pss(Ihru) = Glacr_pkwater_equiv(Ihru)

      END SUBROUTINE glacr_states_to_zero

!***********************************************************************
!     snowcomp_restart - write or read snowcomp restart file
!***********************************************************************
      SUBROUTINE snowcomp_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT, ACTIVE
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit, Glacier_flag
      USE PRMS_SNOW
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=8) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_pweqv, Basin_snowcov, Basin_snowdepth, Basin_snowicecov
        WRITE ( Restart_outunit ) Int_alb
        WRITE ( Restart_outunit ) Scrv
        WRITE ( Restart_outunit ) Pksv
        WRITE ( Restart_outunit ) Snowcov_areasv
        WRITE ( Restart_outunit ) Salb
        WRITE ( Restart_outunit ) Slst
        WRITE ( Restart_outunit ) Lst
        WRITE ( Restart_outunit ) Iasw
        WRITE ( Restart_outunit ) Iso
        WRITE ( Restart_outunit ) Mso
        WRITE ( Restart_outunit ) Lso
        WRITE ( Restart_outunit ) Albedo
        WRITE ( Restart_outunit ) Pk_temp
        WRITE ( Restart_outunit ) Pk_den
        WRITE ( Restart_outunit ) Pk_def
        WRITE ( Restart_outunit ) Pk_ice
        WRITE ( Restart_outunit ) Freeh2o
        WRITE ( Restart_outunit ) Snowcov_area
        WRITE ( Restart_outunit ) Pss
        WRITE ( Restart_outunit ) Pst
        WRITE ( Restart_outunit ) Snsv
        WRITE ( Restart_outunit ) Pk_depth
        IF ( Glacier_flag==ACTIVE ) THEN
          WRITE ( Restart_outunit ) Glacr_albedo
          WRITE ( Restart_outunit ) Glacr_pk_den
          WRITE ( Restart_outunit ) Glacr_pk_ice
          WRITE ( Restart_outunit ) Glacr_freeh2o
          WRITE ( Restart_outunit ) Glacrcov_area
          WRITE ( Restart_outunit ) Glacr_pss
          WRITE ( Restart_outunit ) Glacr_pst
          WRITE ( Restart_outunit ) Glacr_pk_depth
          WRITE ( Restart_outunit ) Glacr_pkwater_equiv
          WRITE ( Restart_outunit ) Glacr_pk_temp
          WRITE ( Restart_outunit ) Ann_tempc, Yrdays5
          WRITE ( Restart_outunit ) Glacr_air_5avtemp, Glacr_air_5avtemp1, Glacr_air_deltemp
          WRITE ( Restart_outunit ) Glacr_5avsnow, Glacr_5avsnow1, Glacr_delsnow
          WRITE ( Restart_outunit ) Glacr_pk_def
          WRITE ( Restart_outunit ) Glacr_freeh2o_capm
        ENDIF
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_pweqv, Basin_snowcov, Basin_snowdepth, Basin_snowicecov
        READ ( Restart_inunit ) Int_alb
        READ ( Restart_inunit ) Scrv
        READ ( Restart_inunit ) Pksv
        READ ( Restart_inunit ) Snowcov_areasv
        READ ( Restart_inunit ) Salb
        READ ( Restart_inunit ) Slst
        READ ( Restart_inunit ) Lst
        READ ( Restart_inunit ) Iasw
        READ ( Restart_inunit ) Iso
        READ ( Restart_inunit ) Mso
        READ ( Restart_inunit ) Lso
        READ ( Restart_inunit ) Albedo
        READ ( Restart_inunit ) Pk_temp
        READ ( Restart_inunit ) Pk_den
        READ ( Restart_inunit ) Pk_def
        READ ( Restart_inunit ) Pk_ice
        READ ( Restart_inunit ) Freeh2o
        READ ( Restart_inunit ) Snowcov_area
        READ ( Restart_inunit ) Pss
        READ ( Restart_inunit ) Pst
        READ ( Restart_inunit ) Snsv
        READ ( Restart_inunit ) Pk_depth
        IF ( Glacier_flag==ACTIVE ) THEN
          READ ( Restart_inunit ) Glacr_albedo
          READ ( Restart_inunit ) Glacr_pk_den
          READ ( Restart_inunit ) Glacr_pk_ice
          READ ( Restart_inunit ) Glacr_freeh2o
          READ ( Restart_inunit ) Glacrcov_area
          READ ( Restart_inunit ) Glacr_pss
          READ ( Restart_inunit ) Glacr_pst
          READ ( Restart_inunit ) Glacr_pk_depth
          READ ( Restart_inunit ) Glacr_pkwater_equiv
          READ ( Restart_inunit ) Glacr_pk_temp
          READ ( Restart_inunit ) Ann_tempc, Yrdays5
          READ ( Restart_inunit ) Glacr_air_5avtemp, Glacr_air_5avtemp1, Glacr_air_deltemp
          READ ( Restart_inunit ) Glacr_5avsnow, Glacr_5avsnow1, Glacr_delsnow
          READ ( Restart_inunit ) Glacr_pk_def
          READ ( Restart_inunit ) Glacr_freeh2o_capm
        ENDIF
      ENDIF
      END SUBROUTINE snowcomp_restart
