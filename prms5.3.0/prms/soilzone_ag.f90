!***********************************************************************
! Computes inflows to and outflows from soil zone of each HRU and
! includes inflows from infiltration, groundwater, and upslope HRUs,
! and outflows to gravity drainage, interflow, and surface runoff to
! downslope HRUs; merge of smbal_prms and ssflow_prms with enhancements
!
! Daily accounting for soil zone;
!    adds infiltration
!    computes et
!    computes recharge of soil zone
!    computes interflow to stream or cascade
!    adjusts storage in soil zone
!    sends dunnian runoff to stream or cascade by adding to sroff
!    computes drainage to groundwater
!***********************************************************************
! FORTRAN module for soilzone_ag
!***********************************************************************
      MODULE PRMS_SOILZONE_AG

      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC_AG = 'Soilzone Computations'
      character(len=11), parameter :: MODNAME_AG = 'soilzone_ag'
      character(len=*), parameter :: Version_soilzone_ag = '2022-01-25'
      INTEGER, SAVE :: Soil_iter, Iter_aet
      DOUBLE PRECISION, SAVE :: Basin_ag_soil_to_gw, Basin_ag_up_max
      DOUBLE PRECISION, SAVE :: Basin_ag_actet, Last_ag_soil_moist, Basin_ag_soil_rechr, Last_ag_soil_rechr
      REAL, SAVE, ALLOCATABLE :: It0_ag_soil_rechr(:), It0_ag_soil_moist(:)
      REAL, SAVE, ALLOCATABLE :: It0_potet(:), It0_sroff(:), It0_strm_seg_in(:)
      !REAL, SAVE, ALLOCATABLE :: Ag_slow_flow(:), Ag_ssres_in(:), Ag_water_maxin(:)
!   Agriculture Declared Variables
      INTEGER, SAVE, ALLOCATABLE :: Ag_soil_saturated(:)
      DOUBLE PRECISION, SAVE :: Basin_ag_waterin
      REAL, SAVE, ALLOCATABLE :: Ag_hortonian(:), Unused_ag_et(:), Ag_soil_to_gvr(:), Ag_soilwater_deficit(:)
      REAL, SAVE, ALLOCATABLE :: Ag_actet(:), Ag_irrigation_add(:), Ag_irrigation_add_vol(:)
      REAL, SAVE, ALLOCATABLE :: Ag_soil_to_gw(:), hru_ag_actet(:)
      REAL, SAVE, ALLOCATABLE :: Ag_soil_lower(:), Ag_soil_lower_stor_max(:), Ag_potet_rechr(:), Ag_potet_lower(:)
!      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Ag_upslope_dunnian(:)
      INTEGER, SAVE :: total_iters, iter_nonconverge
      real, save :: unsatisfied_big
      ! parameters
! have covden a monthly, later
      INTEGER, SAVE, ALLOCATABLE :: Ag_soil_type(:) !, Ag_crop_type(:)
      REAL, SAVE, ALLOCATABLE :: Ag_soilwater_deficit_min(:), Ag_covden_sum(:,:), Ag_covden_win(:,:)
!      REAL, SAVE, ALLOCATABLE :: Ag_sat_threshold(:)
      REAL, SAVE, ALLOCATABLE :: Ag_soil_rechr_max_frac(:), Ag_soil2gw_max(:) ! Ag_crop_coef later, will specify PET
      INTEGER, SAVE :: max_soilzone_ag_iter
      REAL, SAVE :: soilzone_aet_converge

      END MODULE PRMS_SOILZONE_AG

!***********************************************************************
!     Main soilzone_ag routine
!***********************************************************************
      INTEGER FUNCTION soilzone_ag()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
! Functions
      INTEGER, EXTERNAL :: szdecl, szinit, szrun_ag, szdecl_ag, szinit_ag
      EXTERNAL :: soilzone_restart_ag
!***********************************************************************
      soilzone_ag = 0

      IF ( Process_flag==RUN ) THEN
        soilzone_ag = szrun_ag()
      ELSEIF ( Process_flag==DECL ) THEN
        soilzone_ag = szdecl()
        soilzone_ag = szdecl_ag()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL soilzone_restart_ag(READ_INIT)
        soilzone_ag = szinit()
        soilzone_ag = szinit_ag()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL soilzone_restart_ag(SAVE_INIT)
      ENDIF

      END FUNCTION soilzone_ag

!***********************************************************************
!     szdecl_ag - set up parameters for agriculture soil zone computations
!   Declared Parameters
!     sat_threshold, ssstor_init_frac fastcoef_lin, fastcoef_sq
!     ssr2gw_rate, ssr2gw_exp, soil2gw_max, soil_type, ag_soil_type, ag_soilwater_deficit_min
!     soil_rechr_max_frac, soil_rechr_init_frac, soil_moist_max, soil_moist_init_frac
!     pref_flow_den, slowcoef_lin, cov_type
!     hru_area, slowcoef_sq, gvr_hru_id
!***********************************************************************
      INTEGER FUNCTION szdecl_ag()
      USE PRMS_CONSTANTS, ONLY: OFF, ACTIVE, DOCUMENTATION, MONTHS_PER_YEAR
      use PRMS_MMFAPI, only: declvar_dble, declvar_int, declvar_real
      use PRMS_READ_PARAM_FILE, only: declparam, getdim
      USE PRMS_MODULE, ONLY: Nhru, Nlake, AG_flag, iter_aet_flag !, Cascade_flag
      USE PRMS_SOILZONE
      USE PRMS_SOILZONE_AG
      use prms_utils, only: error_stop, print_module, PRMS_open_module_file, read_error
      IMPLICIT NONE
!***********************************************************************
      szdecl_ag = 0

      total_iters = 0

      CALL print_module(MODDESC_AG, MODNAME_AG, Version_soilzone_ag)

      Iter_aet = OFF
      IF ( AG_flag==ACTIVE .OR. iter_aet_flag==ACTIVE ) Iter_aet = ACTIVE
      IF ( Iter_aet==ACTIVE ) THEN
        ALLOCATE ( It0_sroff(Nhru), It0_strm_seg_in(Nhru) )
        IF ( Nlake>0 ) ALLOCATE ( It0_potet(Nhru) )
      ENDIF

! Agriculture variables and parameters
      ALLOCATE ( Ag_soil_to_gw(Nhru) )
      CALL declvar_real(MODNAME, 'ag_soil_to_gw', 'nhru', Nhru, &
     &     'Direct recharge from agriculture capillary reservoir to groundwater reservior for each HRU', &
     &     'inches', Ag_soil_to_gw)

!      IF ( Cascade_flag>OFF ) THEN
!        ALLOCATE ( Ag_upslope_dunnian(Nhru) )
!        CALL declvar_dble(MODNAME, 'ag_upslope_dunnian', 'nhru', Nhru, &
!     &       'Cascading Dunnian surface runoff that'// &
!     &       ' flows to the agriculture capillary reservoir of each downslope HRU for each upslope HRU', &
!     &       'inches', Ag_upslope_dunnian)
!      ENDIF

      ALLOCATE ( Ag_actet(Nhru) )
      CALL declvar_real(MODNAME, 'ag_actet', 'nhru', Nhru, &
     &     'Actual ET for agriculture capillary reservoir for each HRU', &
     &     'inches', Ag_actet)
      ALLOCATE ( hru_ag_actet(Nhru) )
      CALL declvar_real(MODNAME, 'hru_ag_actet', 'nhru', Nhru, &
     &     'Actual ET for agriculture capillary reservoir averaged over each HRU', &
     &     'inches', hru_ag_actet)

      ALLOCATE ( Unused_ag_et(Nhru) )
      CALL declvar_real(MODNAME, 'unused_ag_et', 'nhru', Nhru, &
     &     'Actual ET for agriculture capillary reservoir for each HRU', &
     &     'inches', Unused_ag_et)

      ALLOCATE ( Ag_hortonian(Nhru) )
      CALL declvar_real(MODNAME, 'ag_hortonian', 'nhru', Nhru, &
     &     'Hortonian surface runoff that flows to the stream network from the agriculture fraction of each HRU', &
     &     'inches', Ag_hortonian)

      ALLOCATE ( Ag_soil_to_gvr(Nhru) )
      CALL declvar_real(MODNAME, 'ag_soil_to_gvr', 'nhru', Nhru, &
     &     'Excess capillary water that flows to the agriculture gravity reservoir from the agriculture fraction of each HRU', &
     &     'inches', Ag_soil_to_gvr)

      CALL declvar_dble(MODNAME, 'Basin_ag_waterin', 'one', 1, &
     &     'Basin area-weighted average infiltration,'// &
     &     ' cascading interflow and Dunnian flow added to agriculture reservoir storage', &
     &     'inches', Basin_ag_waterin)

      ALLOCATE ( Ag_irrigation_add(Nhru) )
      CALL declvar_real(MODNAME, 'ag_irrigation_add', 'nhru', Nhru, &
     &     'Irrigation water added to agriculture fraction when ag_actet < PET_external for each HRU', &
     &     'inches', Ag_irrigation_add)

      ALLOCATE ( Ag_soilwater_deficit(Nhru) )
      CALL declvar_real(MODNAME, 'ag_soilwater_deficit', 'nhru', Nhru, &
     &     'Soil-water deficit of agriculture fraction for each HRU', &
     &     'inches', Ag_soilwater_deficit)

      ALLOCATE ( Ag_irrigation_add_vol(Nhru) )
      CALL declvar_real(MODNAME, 'ag_irrigation_add_vol', 'nhru', Nhru, &
     &     'Irrigation water added to agriculture fraction when ag_actet < PET_external for each HRU', &
     &     'acre-inches', Ag_irrigation_add_vol)

      ALLOCATE ( Ag_soil_lower(Nhru), Ag_soil_lower_stor_max(Nhru) )
      CALL declvar_real(MODNAME, 'ag_soil_lower', 'nhru', Nhru, &
     &   'Storage in the lower zone of the agriculture'// &
     &   ' reservoir that is only available for transpiration for each HRU', &
     &   'inches', Ag_soil_lower)

      ALLOCATE ( Ag_potet_lower(Nhru) ) !, Ag_water_maxin(Nhru) )
      CALL declvar_real(MODNAME, 'ag_potet_lower', 'nhru', Nhru, &
     &     'Potential ET in the lower zone of the agriculture reservoir for each HRU', &
     &     'inches', Ag_potet_lower)

      ALLOCATE ( Ag_potet_rechr(Nhru) )
      CALL declvar_real(MODNAME, 'ag_potet_rechr', 'nhru', Nhru, &
     &     'Potential ET in the recharge zone of the agriculture reservoir for each HRU', &
     &     'inches', Ag_potet_rechr)

      ALLOCATE ( Ag_soil_saturated(Nhru) )
      CALL declvar_int(MODNAME, 'ag_soil_saturated', 'nhru', Nhru, &
     &     'Flag set if infiltration saturates capillary reservoir (0=no, 1=yes)', &
     &     'none', Ag_soil_saturated)

      IF ( declparam(MODNAME, 'max_soilzone_ag_iter', 'one', 'integer', &
     &     '10', '1', '9999', &
     &     'Maximum number of iterations to optimize computed AET and input AET', &
     &     'Maximum number of iterations to optimize computed AET and input AET', &
     &     'none')/=0 ) CALL read_error(1, 'max_soilzone_ag_iter')

      IF ( declparam(MODNAME, 'soilzone_aet_converge', 'one', 'real', &
     &     '0.0001', '0.0', '1.0', &
     &     'Convergence criteria to iterate computed AET compared to input AET', &
     &     'Convergence criteria to iterate computed AET compared to input AET', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'soilzone_aet_converge')

      ALLOCATE ( Ag_soil_type(Nhru) )
      IF ( declparam(MODNAME, 'ag_soil_type', 'nhru', 'integer', &
     &     '2', '1', '3', &
     &     'Agriculture soil type', 'Soil type of agriculture in each HRU (1=sand; 2=loam; 3=clay)', &
     &     'none')/=0 ) CALL read_error(1, 'ag_soil_type')

      ALLOCATE ( Ag_soilwater_deficit_min(Nhru) )
      IF ( declparam(MODNAME, 'ag_soilwater_deficit_min', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Minimum soil-water deficit to begin agriculture irrigaition', &
     &     'Minimum soil-water deficit fraction to begin agriculture irrigaition', &
     &     'fraction')/=0 ) CALL read_error(1, 'ag_soilwater_deficit_min')

     ! ALLOCATE ( Ag_sat_threshold(Nhru) )
     ! IF ( declparam(MODNAME, 'ag_sat_threshold', 'nhru', 'real', &
     !&     '999.0', '0.00001', '999.0', &
     !&     'Soil saturation threshold, above field-capacity threshold of agriculture reservoir', &
     !&     'Water holding capacity of the gravity and preferential-'// &
     !&     'flow reservoirs; difference between field capacity and'// &
     !&     ' total soil saturation for each HRU', &
     !&     'inches')/=0 ) CALL read_error(1, 'ag_sat_threshold')

!      ALLOCATE ( Ag_crop_type(Nhru) ) ! find Mastin's code on different crops
!      IF ( declparam(MODNAME, 'ag_crop_type', 'nhru', 'integer', &
!     &     '3', '0', '4', &
!     &     'Agriculture cover type designation for each HRU', &
!     &     'Vegetation cover type for agriculture in each HRU (0=none;'// &
!     &     ' 1=grasses; 2=grain; 3=trees; 4=vegetable)', &
!     &     'none')/=0 ) CALL read_error(1, 'ag_crop_type')

      ALLOCATE ( Ag_covden_sum(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'ag_covden_sum', 'nhru,nmonths', 'real', &
     &     '-1.0', '-1.0', '1.0', &
     &     'Summer vegetation cover density for agriculture crop type', &
     &     'Summer vegetation cover density for the agriculture crop type in each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'ag_covden_sum')

      ALLOCATE ( Ag_covden_win(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'ag_covden_win', 'nhru,nmonths', 'real', &
     &     '-1.0', '-1.0', '1.0', &
     &     'Winter vegetation cover density for crop type', &
     &     'Winter vegetation cover density for the crop type in each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'ag_covden_win')

      ALLOCATE ( Ag_soil2gw_max(Nhru) )
      IF ( declparam(MODNAME, 'ag_soil2gw_max', 'nhru', 'real', &
     &     '-1.0', '0.0', '5.0', &
     &     'Maximum value for agriculture capillary reservoir excess to groundwater storage', &
     &     'Maximum amount of the agriculture capillary reservoir excess that'// &
     &     ' is routed directly to the groundwater storage for each HRU', &
     &     'inches')/=0 ) CALL read_error(1, 'ag_soil2gw_max')

      END FUNCTION szdecl_ag

!***********************************************************************
!     szinit_ag - Initialize soilzone module - get parameter values,
!                 set initial values and check parameter values
!***********************************************************************
      INTEGER FUNCTION szinit_ag()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, LAKE, GLACIER, INACTIVE, OFF, MONTHS_PER_YEAR
      use PRMS_READ_PARAM_FILE, only: getparam_int, getparam_real
      USE PRMS_MODULE, ONLY: Init_vars_from_file, Nhru, Hru_type, Iter_aet_flag
      USE PRMS_SOILZONE, ONLY: MODNAME, Soil2gw_max
      USE PRMS_SOILZONE_AG
      USE PRMS_BASIN, ONLY: Basin_area_inv, Ag_area, Covden_win, Covden_sum
      USE PRMS_FLOWVARS, ONLY: Basin_ag_soil_moist, Ag_soil_moist, Ag_soil_rechr, Ag_soil_moist_max, Ag_soil_rechr_max
      use prms_utils, only: checkdim_bounded_limits, error_stop, read_error
      IMPLICIT NONE
! Functions
      EXTERNAL :: init_basin_vars
      INTRINSIC :: MIN, DBLE
! Local Variables
      INTEGER :: ihru
!***********************************************************************
      szinit_ag = 0

!??? figure out what to save in restart file ???
      IF ( getparam_int(MODNAME, 'max_soilzone_ag_iter', 1, max_soilzone_ag_iter)/=0 ) &
     &     CALL read_error(2, 'max_soilzone_ag_iter')
      IF ( getparam_real(MODNAME, 'soilzone_aet_converge', 1, soilzone_aet_converge)/=0 ) &
     &     CALL read_error(2, 'soilzone_aet_converge')
      IF ( getparam_int(MODNAME, 'ag_soil_type', Nhru, Ag_soil_type)/=0 ) CALL read_error(2, 'ag_soil_type')
      IF ( getparam_real(MODNAME, 'ag_soilwater_deficit_min', Nhru, Ag_soilwater_deficit_min)/=0 ) &
     &     CALL read_error(2, 'ag_soilwater_deficit_min')
!      IF ( getparam_int(MODNAME, 'ag_crop_type', Nhru, Ag_crop_type)/=0 ) CALL read_error(2, 'ag_crop_type')
      IF ( getparam_real(MODNAME, 'ag_covden_sum', Nhru*MONTHS_PER_YEAR, Ag_covden_sum)/=0 ) CALL read_error(2, 'ag_covden_sum')
      IF ( Ag_covden_sum(1,1)<0.0 ) Ag_covden_sum = Covden_sum
      IF ( getparam_real(MODNAME, 'ag_covden_win', Nhru*MONTHS_PER_YEAR, Ag_covden_win)/=0 ) CALL read_error(2, 'ag_covden_win')
      IF ( Ag_covden_win(1,1)<0.0 ) Ag_covden_win = Covden_win
      IF ( getparam_real(MODNAME, 'ag_soil2gw_max', Nhru, Ag_soil2gw_max)/=0 ) CALL read_error(2, 'ag_soil2gw_max')
      IF ( Ag_soil2gw_max(1)<0.0 ) Ag_soil2gw_max = Soil2gw_max
      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==5 ) Ag_soil_lower = 0.0
      ! dimensioned nhru
      Ag_soil_to_gw = 0.0
      Ag_hortonian = 0.0
      Ag_soil_to_gvr = 0.0
      Ag_soil_lower_stor_max = 0.0
      Ag_potet_lower = 0.0
      Ag_potet_rechr = 0.0
      Ag_actet = 0.0
      Unused_ag_et = 0.0
      Ag_soilwater_deficit = 0.0
      Basin_ag_soil_moist = 0.0D0
      Basin_ag_soil_rechr = 0.0D0
      Ag_soil_saturated = OFF
      DO ihru = 1, Nhru
        ! make sure LAKE, INACTIVE, GLACIER have agriculture values of 0
        IF ( Hru_type(ihru)==LAKE .OR. Hru_type(ihru)==INACTIVE .OR. Hru_type(ihru)==GLACIER ) Ag_area(ihru) = 0.0
        IF ( Ag_area(ihru)>0.0 ) THEN
          Basin_ag_soil_moist = Basin_ag_soil_moist + DBLE( Ag_soil_moist(ihru)*Ag_area(ihru) )
          Basin_ag_soil_rechr = Basin_ag_soil_rechr + DBLE( Ag_soil_rechr(ihru)*Ag_area(ihru) )
        ELSE
          Ag_soil_moist_max(ihru) = 0.0
          Ag_soil_rechr_max(ihru) = 0.0
          Ag_soil_moist(ihru) = 0.0
          Ag_soil_rechr(ihru) = 0.0
        ENDIF
      ENDDO
      Basin_ag_soil_moist = Basin_ag_soil_moist*Basin_area_inv
      Basin_ag_soil_rechr = Basin_ag_soil_rechr*Basin_area_inv

      IF ( Iter_aet_flag==ACTIVE ) ALLOCATE ( It0_ag_soil_rechr(Nhru), It0_ag_soil_moist(Nhru) )

      Soil_iter = 1
      iter_nonconverge = 0
      unsatisfied_big = 0.0

      END FUNCTION szinit_ag

!***********************************************************************
!     szrun_ag - Does soil water balance for each HRU, adds in infiltration
!                then computes actual et and apportions remainder between
!                recharge of soil moisture, soil storage available for
!                interflow, excess routed to stream,
!                and groundwater reservoirs
!***********************************************************************
      INTEGER FUNCTION szrun_ag()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, NEARZERO, LAND, LAKE, SWALE, GLACIER, &
     &    DEBUG_less, DEBUG_WB, ERROR_param, CASCADE_OFF
      USE PRMS_MODULE, ONLY: Nlake, Print_debug, Dprst_flag, Cascade_flag, &
     &    Frozen_flag, Soilzone_add_water_use, Call_cascade, &
     &    Nowmonth, Nowyear, Nowday, Iter_aet_flag, Hru_type
      USE PRMS_SOILZONE
      USE PRMS_SOILZONE_AG
      USE PRMS_BASIN, ONLY: Hru_perv, Hru_frac_perv, Hru_storage, &
     &    Hru_route_order, Active_hrus, Basin_area_inv, Hru_area, &
     &    Lake_hru_id, Cov_type, Numlake_hrus, Hru_area_dble, Ag_frac, Ag_area, Ag_cov_type
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt, Transp_on, Potet, Basin_potet, Basin_transp_on
! WARNING!!! Sroff, Basin_sroff, and Strm_seg_in can be updated
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_actet, Hru_actet, &
     &    Ssres_flow, Soil_to_gw, Basin_soil_to_gw, Ssr_to_gw, &
     &    Soil_to_ssr, Basin_lakeevap, Basin_perv_et, Basin_swale_et, &
     &    Sroff, Soil_moist_max, Infil, Soil_rechr_max, Ssres_in, Snowcov_area, Snow_evap, &
     &    Basin_soil_moist, Basin_ssstor, Slow_stor, Slow_flow, Pkwater_equiv, &
     &    Ssres_stor, Soil_moist, Sat_threshold, Soil_rechr, Basin_sroff, Basin_lake_stor, &
     &    Ag_soil_rechr, Ag_soil_moist, Ag_soil_rechr_max, Ag_soil_moist_max, Basin_ag_soil_moist, &
     &    Soil_moist_ante, Soil_rechr_ante, Ssres_stor_ante, Slow_stor_ante, Pref_flow_stor_ante
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SRUNOFF, ONLY: Hru_impervstor, Dprst_stor_hru
      USE PRMS_WATER_USE, ONLY: Soilzone_gain, Soilzone_gain_hru
      USE PRMS_CLIMATE_HRU, ONLY: AET_external, PET_external
      USE PRMS_CASCADE, ONLY: Ncascade_hru
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_INTCP, ONLY: Hru_intcpevap
      USE PRMS_SRUNOFF, ONLY: Hru_impervevap, Strm_seg_in, Dprst_evap_hru, Dprst_seep_hru, Frozen, Infil_ag
      use prms_utils, only: error_stop, print_date
      IMPLICIT NONE
! Functions
      INTRINSIC :: MIN, ABS, MAX, SNGL, DBLE
      EXTERNAL :: compute_soilmoist, compute_szactet, compute_cascades
      EXTERNAL :: compute_interflow, compute_gwflow, init_basin_vars
! Local Variables
      INTEGER :: i, k, update_potet, compute_lateral, perv_on_flag
      REAL :: dunnianflw, interflow, perv_area, harea
      REAL :: dnslowflow, dnpreflow, dndunn, availh2o, avail_potet, hruactet, ag_hruactet
      REAL :: topfr !, depth, tmp
      REAL :: dunnianflw_pfr, dunnianflw_gvr, pref_flow_maxin
      REAL :: perv_frac, capwater_maxin, ssresin
      REAL :: cap_upflow_max, unsatisfied_et, pervactet, prefflow, ag_water_maxin
      REAL :: ag_upflow_max, ag_capacity, excess, agfrac, ag_avail_potet, ag_potet
      REAL :: ag_AETtarget, ag_avail_targetAET, cap_ag_water_maxin, agactet
      INTEGER :: cfgi_frozen_hru
      INTEGER :: num_hrus_ag_iter, ag_on_flag, keep_iterating, add_estimated_irrigation
!***********************************************************************
      szrun_ag = 0

! It0 and _ante variables to save iteration states.
      Soil_moist_ante = Soil_moist
      Soil_rechr_ante = Soil_rechr
      Ssres_stor_ante = Ssres_stor
      Slow_stor_ante = Slow_stor
      IF ( Pref_flag==ACTIVE ) Pref_flow_stor_ante = Pref_flow_stor
      Last_soil_moist = Basin_soil_moist
      Last_ssstor = Basin_ssstor
      IF ( Iter_aet==ACTIVE ) THEN
        ! computed in srunoff
        It0_sroff = Sroff
        IF ( Call_cascade==ACTIVE ) It0_strm_seg_in = Strm_seg_in
        Ag_irrigation_add = 0.0
        It0_ag_soil_moist = Ag_soil_moist
        It0_ag_soil_rechr = Ag_soil_rechr
        IF ( Nlake>0 ) It0_potet = Potet
      ENDIF
      Last_ag_soil_moist = Basin_ag_soil_moist
      Last_ag_soil_rechr = Basin_ag_soil_rechr

      keep_iterating = ACTIVE
      Soil_iter = 1

! ***************************************
      DO WHILE ( keep_iterating==ACTIVE )
! ***************************************

      IF ( Soil_iter>1 ) THEN
        DO k = 1, Active_hrus
          i = Hru_route_order(k)
          Soil_moist(i) = Soil_moist_ante(i)
          Soil_rechr(i) = Soil_rechr_ante(i)
          Ssres_stor(i) = Ssres_stor_ante(i)
          Slow_stor(i) = Slow_stor_ante(i)
          IF ( Pref_flag==ACTIVE ) Pref_flow_stor(i) = Pref_flow_stor_ante(i)
          IF ( Nlake>0 ) Potet(i) = It0_potet(i)
        ENDDO
        IF ( Iter_aet_flag==ACTIVE ) THEN
          ! computed in srunoff
          Sroff = It0_sroff
          IF ( Call_cascade==ACTIVE ) Strm_seg_in = It0_strm_seg_in
        ENDIF
      ENDIF

      IF ( Cascade_flag>CASCADE_OFF ) THEN
        Upslope_interflow = 0.0D0
        Upslope_dunnianflow = 0.0D0
        IF ( Numlake_hrus>0 ) THEN
          Lakein_sz = 0.0D0
          Basin_lakeinsz = 0.0D0
        ENDIF
      ENDIF

      IF ( Soil_iter>1 ) THEN
        Ag_soil_moist = It0_ag_soil_moist
        Ag_soil_rechr = It0_ag_soil_rechr
      ENDIF
      Basin_ag_soil_moist = 0.0D0
      Basin_ag_soil_rechr = 0.0D0
      Basin_ag_soil_to_gw = 0.0D0
      Basin_ag_up_max = 0.0D0
      Basin_ag_actet = 0.0D0
      Basin_ag_waterin = 0.0D0
      CALL init_basin_vars()
      update_potet = OFF
      unsatisfied_big = 0.0
      IF ( Soilzone_add_water_use==ACTIVE ) Soilzone_gain_hru = 0.0
      add_estimated_irrigation = OFF
      num_hrus_ag_iter = 0

! ***************************************
      DO k = 1, Active_hrus
        i = Hru_route_order(k)
! ***************************************
        ! Soil_to_gw and Soil_to_ssr for whole HRU
        Soil_to_gw(i) = 0.0
        Soil_to_ssr(i) = 0.0
        ! gravity reservoir variables for whole HRU
        Ssr_to_gw(i) = 0.0
        Slow_flow(i) = 0.0
        Ssres_flow(i) = 0.0
        Cap_waterin(i) = 0.0
        Soil_saturated(i) = OFF
! initialize all HRU values in case dynamic ag frac
        Ag_soil_saturated(i) = OFF
        Ag_soil_to_gvr(i) = 0.0
        Ag_soil_to_gw(i) = 0.0
        Ag_hortonian(i) = 0.0
        Unused_ag_et(i) = 0.0
        Ag_soilwater_deficit(i) = 0.0
        hruactet = Hru_impervevap(i) + Hru_intcpevap(i) + Snow_evap(i)
        IF ( Dprst_flag==ACTIVE ) hruactet = hruactet + Dprst_evap_hru(i)
        harea = Hru_area(i)
        agfrac = Ag_frac(i)

        IF ( Hru_type(i)==LAKE ) THEN ! lake or reservoir
          !WARNING, RSR, if hru_actet>water in lake, then budget error
          hruactet = (Potet(i) - hruactet)*Lake_evap_adj(Nowmonth,Lake_hru_id(i))
          IF ( hruactet>Potet(i) ) THEN
            PRINT *, 'WARNING, lake evap > potet, for HRU:', i, ' potential ET increased to adjusted lake ET'
            PRINT *, hruactet, Potet(i), hruactet - Potet(i)
            Potet(i) = hruactet ! this could be a problem when it happens
            update_potet = ACTIVE
          ENDIF
          Unused_potet(i) = Potet(i) - hruactet
          Basin_actet = Basin_actet + DBLE( hruactet*harea )
          Basin_lakeevap = Basin_lakeevap + DBLE( hruactet*harea )
          Basin_lakeprecip = Basin_lakeprecip + DBLE( Hru_ppt(i)*harea )
          IF ( Cascade_flag>CASCADE_OFF ) THEN
            ! if lake HRU doesn't cascade, should we limit ET to
            !  water entering the HRU to this point (no gwflow yet)
            Lakein_sz(i) = Upslope_interflow(i) + Upslope_dunnianflow(i)
            Basin_lakeinsz = Basin_lakeinsz + Lakein_sz(i)*Hru_area_dble(i)
          ENDIF
          Hru_actet(i) = hruactet
          CYCLE
        ENDIF

        !Hru_type can be 1 (land) or 3 (swale) or 4 (glacier)
        compute_lateral = OFF ! swale
        IF ( Hru_type(i)==LAND .OR. Hru_type(i)==GLACIER ) compute_lateral = ACTIVE
        perv_area = Hru_perv(i)
        perv_frac = Hru_frac_perv(i)
        perv_on_flag = OFF
        IF ( perv_area>0.0 ) perv_on_flag = ACTIVE
        ag_water_maxin = 0.0
        ag_on_flag = OFF
        IF ( Ag_area(i)>0.0 ) THEN
          ag_on_flag = ACTIVE
          ag_water_maxin = Infil_ag(i)
        ENDIF

        avail_potet = Potet(i) - hruactet
        IF ( avail_potet<0.0 ) THEN
          print *, 'avail_potet<0', avail_potet, Potet(i), Hru_impervevap(i), Hru_intcpevap(i), Snow_evap(i), hruactet
          avail_potet = 0.0
          hruactet = Potet(i)
        ENDIF

        !Hru_type can be 1 (land) or 3 (swale) or 4 (glacier)

!******Add infiltration to soil and compute excess
        dunnianflw = 0.0
        dunnianflw_pfr = 0.0
        dunnianflw_gvr = 0.0
        interflow = 0.0

!******Add infiltration to soil and compute excess
        !infil_tot is the depth in whole HRU
        !capillary reservoir for pervious area
        !agriculture reservoir for irrigated area
        !preferential flow reservoir for whole HRU
        !gravity reservoir for whole HRU
        !upslope flow for whole HRU

!******if cascading flow available from upslope cascades
!****** add soil excess (Dunnian flow) to infiltration
        ! infil for pervious portion of HRU
        ! infil_ag for pervious and agriculture portion of HRU
        capwater_maxin = Infil(i)

        IF ( Iter_aet_flag==ACTIVE ) ag_water_maxin = ag_water_maxin + Ag_irrigation_add(i) ! units of inches over Ag_area
        IF ( Soilzone_add_water_use==ACTIVE ) THEN
          IF ( Soilzone_gain(i)>0.0 ) THEN
            IF ( perv_on_flag==OFF ) THEN
              PRINT *, 'perv_area=0.0 for HRU:', i
              CALL error_stop('soilzone gain specified and perv_area=0.0', ERROR_param)
            ENDIF
            Soilzone_gain_hru(i) = Soilzone_gain(i)/perv_area/SNGL(Cfs_conv) ! ??? is this harea
            cap_ag_water_maxin = Soilzone_gain_hru(i)
          ENDIF
          capwater_maxin = capwater_maxin + cap_ag_water_maxin
        ENDIF
        IF ( ag_on_flag==ACTIVE ) THEN
          Ag_soilwater_deficit(i) = Ag_soil_moist_max(i) - Ag_soil_moist(i)
          excess = MAX( 0.0, ag_water_maxin + Ag_soil_moist(i) - Ag_soil_moist_max(i) )
          IF ( excess>0.0 ) THEN
            Ag_hortonian(i) = excess
            Sroff(i) = Sroff(i) + excess
            ag_water_maxin = ag_water_maxin - excess
          ENDIF
        ENDIF

        cfgi_frozen_hru = OFF
        !Frozen is HRU variable that says if frozen gravity reservoir
        ! For CFGI all inflow is assumed to be Dunnian Flow when frozen
        IF ( Frozen_flag==ACTIVE ) THEN
          IF ( Frozen(i)==ACTIVE ) THEN
            IF ( Hru_type(i)==SWALE ) THEN
              PRINT *, 'ERROR, a swale HRU cannot be frozen for CFGI, HRU:', i
              ERROR STOP ERROR_param
            ENDIF
            cfgi_frozen_hru = ACTIVE
          ENDIF
        ENDIF

        ! compute preferential flow and storage, and any dunnian flow
        ! pref_flow for whole HRU
! ??? should cascading flow go to preferential flow fraction ???
        prefflow = 0.0
        dunnianflw_pfr = 0.0
        IF ( Pref_flow_infil_frac(i)>0.0 ) THEN
          pref_flow_maxin = 0.0
          Pref_flow_infil(i) = 0.0
          IF ( capwater_maxin>0.0 ) THEN
            ! pref_flow for whole HRU
            pref_flow_maxin = capwater_maxin*Pref_flow_infil_frac(i)
            capwater_maxin = capwater_maxin - pref_flow_maxin
            pref_flow_maxin = pref_flow_maxin*perv_frac
          ENDIF
          IF ( ag_water_maxin>0.0 .AND. ag_on_flag==ACTIVE ) THEN
            pref_flow_maxin = ag_water_maxin*Pref_flow_infil_frac(i)
            ag_water_maxin = ag_water_maxin - pref_flow_maxin
            pref_flow_maxin = pref_flow_maxin*agfrac
          ENDIF
          IF ( pref_flow_maxin>0.0 ) THEN
            IF ( cfgi_frozen_hru==ACTIVE ) THEN
              dunnianflw_pfr = pref_flow_maxin
            ELSE
              ! compute contribution to preferential-flow reservoir storage
              Pref_flow_stor(i) = Pref_flow_stor(i) + pref_flow_maxin
              dunnianflw_pfr = MAX( 0.0, Pref_flow_stor(i)-Pref_flow_max(i) )
            ENDIF
            IF ( dunnianflw_pfr>0.0 ) THEN
              Basin_dunnian_pfr = Basin_dunnian_pfr + DBLE( dunnianflw_pfr*harea )
              Pref_flow_stor(i) = Pref_flow_max(i)
            ENDIF
            Pref_flow_infil(i) = pref_flow_maxin - dunnianflw_pfr
            Basin_pref_flow_infil = Basin_pref_flow_infil + DBLE( Pref_flow_infil(i)*harea )
          ENDIF
          Pfr_dunnian_flow(i) = dunnianflw_pfr
        ENDIF

        IF ( Cascade_flag>CASCADE_OFF ) THEN
          IF ( perv_on_flag==ACTIVE ) THEN
            cap_upflow_max = SNGL(Upslope_dunnianflow(i)+Upslope_interflow(i))/perv_frac
            capwater_maxin = capwater_maxin + cap_upflow_max
            Basin_cap_up_max = Basin_cap_up_max + DBLE( cap_upflow_max*perv_area )
          ENDIF
          IF ( ag_on_flag==ACTIVE ) THEN
            ag_upflow_max = SNGL(Upslope_dunnianflow(i)+Upslope_interflow(i))/agfrac
            Basin_ag_up_max = Basin_ag_up_max + ag_upflow_max*Ag_area(i)
            ag_water_maxin = ag_water_maxin + ag_upflow_max
          ENDIF
        ENDIF
        IF ( perv_on_flag==ACTIVE ) THEN
          Cap_infil_tot(i) = capwater_maxin*perv_frac
          Basin_cap_infil_tot = Basin_cap_infil_tot + DBLE( Cap_infil_tot(i)*harea )
        ENDIF

!******Add infiltration to soil and compute excess
        Cap_waterin(i) = capwater_maxin

        IF ( cfgi_frozen_hru==OFF ) THEN
          IF ( perv_on_flag==ACTIVE ) THEN
            ! call even if capwater_maxin = 0, just in case soil_moist now > Soil_moist_max
            IF ( capwater_maxin+Soil_moist(i)>0.0 ) THEN
              CALL compute_soilmoist(Cap_waterin(i), Soil_moist_max(i), &
     &                               Soil_rechr_max(i), Soil2gw_max(i), Soil_to_ssr(i), &
     &                               Soil_moist(i), Soil_rechr(i), Soil_to_gw(i), perv_frac)
              Cap_waterin(i) = Cap_waterin(i)*perv_frac
              Basin_capwaterin = Basin_capwaterin + DBLE( Cap_waterin(i)*harea )
            ENDIF
          ENDIF
          IF ( ag_on_flag==ACTIVE ) THEN
            IF ( ag_water_maxin+Ag_soil_moist(i)>0.0 ) THEN
              if ( Ag_soil_moist(i)<Ag_soil_rechr(i) ) print *, 'AG1 soilrechr, before', i, &
     &             Ag_soil_moist(i)-Ag_soil_rechr(i), Ag_soil_moist(i), Ag_soil_rechr(i), Ag_soil_moist_max(i), Ag_soil_rechr_max(i)
              CALL compute_soilmoist(ag_water_maxin, Ag_soil_moist_max(i), &
     &                               Ag_soil_rechr_max(i), Soil2gw_max(i), Ag_soil_to_gvr(i), &
     &                               Ag_soil_moist(i), Ag_soil_rechr(i), Ag_soil_to_gw(i), agfrac)
              if ( Ag_soil_moist(i)< Ag_soil_rechr(i)) print *, 'AG1 soilrechr, after', i, &
                   Ag_soil_moist(i)-Ag_soil_rechr(i), Ag_soil_moist(i), Ag_soil_rechr(i), Ag_soil_moist_max(i), Ag_soil_rechr_max(i)
              ag_water_maxin = ag_water_maxin*agfrac
!              Ag_water_maxin(i) = ag_water_maxin
              Basin_ag_waterin = Basin_ag_waterin + DBLE( ag_water_maxin*harea )
              Soil_to_gw(i) = Soil_to_gw(i) + Ag_soil_to_gw(i)
            ENDIF
          ENDIF
          Basin_soil_to_gw = Basin_soil_to_gw + DBLE( Soil_to_gw(i)*harea )
          Soil_to_ssr(i) = Soil_to_ssr(i) + Ag_soil_to_gvr(i)
          Basin_sm2gvr_max = Basin_sm2gvr_max + DBLE( Soil_to_ssr(i)*harea )
        ENDIF

! compute slow interflow and ssr_to_gw
        topfr = 0.0
        ag_capacity = 0.0
        IF ( ag_on_flag==ACTIVE ) ag_capacity = (Ag_soil_moist_max(i) - Ag_soil_moist(i))*agfrac
        availh2o = Slow_stor(i) + Soil_to_ssr(i)
        IF ( compute_lateral==ACTIVE ) THEN
          topfr = MAX( 0.0, availh2o-Pref_flow_thrsh(i) )
          ssresin = Soil_to_ssr(i) - topfr
          Slow_stor(i) = availh2o - topfr
          ! compute slow contribution to interflow, if any
          IF ( Slow_stor(i)>0.0 ) &
     &         CALL compute_interflow(Slowcoef_lin(i), Slowcoef_sq(i), ssresin, Slow_stor(i), Slow_flow(i))
        ELSEIF ( Hru_type(i)==SWALE ) THEN
          Slow_stor(i) = availh2o
        ENDIF
        IF ( Slow_stor(i)>0.0 .AND. Ssr2gw_rate(i)>0.0 ) &
     &       CALL compute_gwflow(Ssr2gw_rate(i), Ssr2gw_exp(i), Ssr_to_gw(i), Slow_stor(i))

        ! compute contribution to Dunnian flow from PFR, if any
        IF ( Pref_flow_den(i)>0.0 ) THEN
          availh2o = Pref_flow_stor(i) + topfr
          dunnianflw_gvr = MAX( 0.0, availh2o-Pref_flow_max(i) )
          IF ( dunnianflw_gvr>0.0 ) THEN
            topfr = topfr - dunnianflw_gvr
            IF ( topfr<0.0 ) THEN
!              IF ( topfr<-NEARZERO .AND. Print_debug>DEBUG_less ) PRINT *, 'gvr2pfr<0', topfr, dunnianflw_gvr, &
!     &             Pref_flow_max(i), Pref_flow_stor(i), Soil_to_ssr(i)
              topfr = 0.0
            ENDIF
          ENDIF
          Pref_flow_in(i) = Pref_flow_infil(i) + topfr
          Pref_flow_stor(i) = Pref_flow_stor(i) + topfr
          IF ( Pref_flow_stor(i)>0.0 ) &
     &         CALL compute_interflow(Fastcoef_lin(i), Fastcoef_sq(i), &
     &                                Pref_flow_in(i), Pref_flow_stor(i), prefflow)
          Basin_pref_stor = Basin_pref_stor + DBLE( Pref_flow_stor(i)*harea )
          IF ( Pref_flow_max(i)>0.0 ) Basin_pfr_stor_frac = Basin_pfr_stor_frac + DBLE( Pref_flow_stor(i)/Pref_flow_max(i)*harea )
        ELSEIF ( compute_lateral==ACTIVE ) THEN
          dunnianflw_gvr = topfr  !?? is this right
        ENDIF
        Gvr2pfr(i) = topfr

        Basin_sm2gvr = Basin_sm2gvr + DBLE( Soil_to_ssr(i)*harea )
        Basin_dunnian_gvr = Basin_dunnian_gvr + DBLE( dunnianflw_gvr*harea )
        Basin_sz2gw = Basin_sz2gw + DBLE( Ssr_to_gw(i)*harea )

!******Compute actual evapotranspiration
        Snow_free(i) = 1.0 - Snowcov_area(i)
        Potet_rechr(i) = 0.0
        Potet_lower(i) = 0.0
        pervactet = 0.0
        agactet = 0.0
        ag_hruactet = 0.0
        unsatisfied_et = 0.0
        IF ( cfgi_frozen_hru==OFF ) THEN
          IF ( Soil_moist(i)>0.0 .AND. avail_potet>0.0 ) THEN
            CALL compute_szactet(Soil_moist_max(i), Soil_rechr_max(i), Transp_on(i), Cov_type(i), &
     &                           Soil_type(i), Soil_moist(i), Soil_rechr(i), pervactet, avail_potet, &
     &                           Snow_free(i), Potet_rechr(i), Potet_lower(i), &
     &                           Potet(i), perv_frac, Soil_saturated(i))
          ENDIF

          IF ( ag_on_flag==ACTIVE .AND. Iter_aet==ACTIVE .AND. Ag_soil_moist(i)>0.0 ) THEN
            IF ( Iter_aet_flag==ACTIVE ) THEN
              !soilwater_deficit = MAX( sz_deficit_param, ag_soil_moist(i)/ag_soil_moist_max(i) ) irrigation happens at a deficit threshold
              ag_AETtarget = AET_external(i) ! ??? rsr, should this be PET target
              IF ( AET_external(i)>PET_external(i) ) THEN
                PRINT *, 'external AET > PET ', AET_external(i), PET_external(i)
                PRINT *, PET_external(i), Potet(i), ag_avail_targetAET, AET_external(i), i, Ag_irrigation_add(i)
                AET_external(i) = PET_external(i)
              ENDIF
            ELSE
              ag_AETtarget = Potet(i)
            ENDIF

            ag_avail_targetAET = ag_AETtarget - hruactet

            IF ( ag_avail_targetAET<0.0 ) ag_avail_targetAET = 0.0
            IF ( ag_avail_targetAET>0.0 ) THEN
              if ( Ag_soil_moist(i)< Ag_soil_rechr(i)) print *, 'AG1 szactet, before', i, Ag_soil_moist(i)-Ag_soil_rechr(i), &
              Ag_soil_moist(i), Ag_soil_rechr(i), Ag_soil_moist_max(i), Ag_soil_rechr_max(i)
              CALL compute_szactet(Ag_soil_moist_max(i), Ag_soil_rechr_max(i), Transp_on(i), Ag_cov_type(i), &
     &                             Ag_soil_type(i), Ag_soil_moist(i), Ag_soil_rechr(i), agactet, ag_avail_targetAET, & !?? instead of ag_avail_potet use AET_external
     &                             Snow_free(i), Ag_potet_rechr(i), Ag_potet_lower(i), &
     &                             ag_AETtarget, agfrac, Ag_soil_saturated(i))
              if ( Ag_soil_moist(i)< Ag_soil_rechr(i)) print *, 'AG1 szactet, after', i, Ag_soil_moist(i)-Ag_soil_rechr(i), &
              Ag_soil_moist(i), Ag_soil_rechr(i), Ag_soil_moist_max(i), Ag_soil_rechr_max(i)
              ! sanity check
              IF ( Iter_aet_flag==ACTIVE ) THEN
                IF ( agactet-ag_AETtarget>NEARZERO ) THEN
                  PRINT *, 'ag_actet problem', agactet, ag_avail_potet, agfrac, AET_external(i), ag_potet, i, hruactet
                  PRINT *, ag_AETtarget-agactet, Ag_soil_moist(i), Ag_soil_rechr(i)
                ENDIF
              ENDIF
              ag_hruactet = agactet*agfrac

              avail_potet = ag_AETtarget - (hruactet + ag_hruactet)
              IF ( ag_AETtarget<agactet ) THEN
                PRINT *, 'WARNING, external agriculture available target AET from CBH File < computed AET', i, &
     &                   Nowyear, Nowmonth, Nowday, num_hrus_ag_iter
                PRINT '(4(A,F0.6))', '         AET_external: ', ag_AETtarget, '; ag_actet: ', &
     &                agactet, ' PET_external: ', PET_external(i), ' AET_external: ', AET_external(i)
                print *, Hru_impervevap(i), Hru_intcpevap(i), Snow_evap(i)
              ENDIF
              Unused_ag_et(i) = ag_AETtarget - agactet
            ENDIF
            unsatisfied_et = ag_AETtarget - agactet
            if ( unsatisfied_et<0.0 ) print *, unsatisfied_et, i, 'unsat', soilzone_aet_converge, ag_AETtarget, Ag_actet(i)
          ENDIF
          Unused_potet(i) = Unused_potet(i) - Unused_ag_et(i)
        ENDIF
        Ag_actet(i) = agactet
        hru_ag_actet(i) = ag_hruactet
        hruactet = hruactet + pervactet*perv_frac
        Hru_actet(i) = hruactet + ag_hruactet
        avail_potet = Potet(i) - hruactet
        ! sanity check
        IF ( avail_potet<0.0 ) THEN
          IF ( Print_debug>-1 ) THEN
            IF ( avail_potet<-NEARZERO ) THEN
              PRINT *, 'hru_actet>potet', i, hruactet, &
     &                 Nowmonth, Nowday, Hru_actet(i), Potet(i), avail_potet, ag_hruactet, &
     &                 pervactet*perv_frac, perv_frac, agfrac
              print *, 'potet', Potet(i), 'external', PET_external(i)
              PRINT *, 'hruactet', hruactet, Hru_impervevap(i) + Hru_intcpevap(i) + Snow_evap(i), &
     &                 Hru_impervevap(i), Hru_intcpevap(i), Snow_evap(i)
            ENDIF
            avail_potet = 0.0
          ENDIF
!          Hru_actet(i) = Potet(i)
!          IF ( perv_on_flag==ACTIVE ) THEN
!            tmp = avail_potet/perv_frac
!            pervactet = pervactet + tmp
!            Soil_moist(i) = Soil_moist(i) - tmp
!            Soil_rechr(i) = Soil_rechr(i) - tmp
!            IF ( Soil_rechr(i)<0.0 ) Soil_rechr(i) = 0.0
!            IF ( Soil_moist(i)<0.0 ) Soil_moist(i) = 0.0
!          ENDIF
        ENDIF
        Perv_actet(i) = pervactet
        hru_perv_actet(i) = pervactet * perv_frac

! soil_moist & soil_rechr multiplied by perv_area instead of harea
        Soil_lower(i) = Soil_moist(i) - Soil_rechr(i)
        Basin_soil_moist = Basin_soil_moist + DBLE( Soil_moist(i)*perv_area )
        Basin_soil_rechr = Basin_soil_rechr + DBLE( Soil_rechr(i)*perv_area )
        Basin_perv_et = Basin_perv_et + DBLE( Perv_actet(i)*perv_area )

! if HRU cascades,
! compute interflow and excess flow to each HRU or stream
        IF ( compute_lateral==ACTIVE ) THEN
          interflow = Slow_flow(i) + prefflow
          Basin_interflow_max = Basin_interflow_max + interflow*harea
          dunnianflw = dunnianflw_gvr + dunnianflw_pfr
          Dunnian_flow(i) = dunnianflw
          IF ( Cascade_flag>CASCADE_OFF ) THEN
            IF ( Ncascade_hru(i)>0 ) THEN
              dnslowflow = 0.0
              dnpreflow = 0.0
              dndunn = 0.0
              IF ( interflow+dunnianflw>0.0 ) THEN
                CALL compute_cascades(i, Ncascade_hru(i), Slow_flow(i), &
     &                                prefflow, Dunnian_flow(i), dnslowflow, &
     &                                dnpreflow, dndunn)
                Basin_dninterflow = Basin_dninterflow + DBLE( (dnslowflow+dnpreflow)*harea )
                Basin_dndunnianflow = Basin_dndunnianflow + DBLE( dndunn*harea )
              ENDIF
              Hru_sz_cascadeflow(i) = dnslowflow + dnpreflow + dndunn
              Basin_dncascadeflow = Basin_dncascadeflow + DBLE( Hru_sz_cascadeflow(i)*harea )
            ENDIF
          ENDIF

! treat pref_flow as interflow
          Ssres_flow(i) = Slow_flow(i)
          IF ( Pref_flow_den(i)>0.0 ) THEN
            Pref_flow(i) = prefflow
            Ssres_flow(i) = Ssres_flow(i) + prefflow
            Basin_prefflow = Basin_prefflow + DBLE( prefflow*harea )
            Basin_gvr2pfr = Basin_gvr2pfr + DBLE( Gvr2pfr(i)*harea )
          ENDIF
          Basin_ssflow = Basin_ssflow + DBLE( Ssres_flow(i)*harea )
          Basin_slowflow = Basin_slowflow + DBLE( Slow_flow(i)*harea )

! treat dunnianflw as surface runoff to streams
          Sroff(i) = Sroff(i) + Dunnian_flow(i)
          Basin_sroff = Basin_sroff + DBLE( Sroff(i)*harea )
          Basin_dunnian = Basin_dunnian + DBLE( Dunnian_flow(i)*harea )
          Ssres_stor(i) = Slow_stor(i) + Pref_flow_stor(i)

        ELSE ! for swales
          availh2o = Slow_stor(i) - Sat_threshold(i)
          Swale_actet(i) = 0.0
          IF ( availh2o>0.0 ) THEN ! if ponding, as storage > sat_threshold
            unsatisfied_et = Potet(i) - Hru_actet(i)
            IF ( unsatisfied_et>0.0 ) THEN
              availh2o = MIN ( availh2o, unsatisfied_et )
              Swale_actet(i) = availh2o
              Hru_actet(i) = Hru_actet(i) + Swale_actet(i)
              Slow_stor(i) = Slow_stor(i) - Swale_actet(i)
              Basin_swale_et = Basin_swale_et + DBLE( Swale_actet(i)*harea )
            ENDIF
            IF ( Print_debug==7 ) THEN
              IF ( Slow_stor(i)>Swale_limit(i) ) THEN
                WRITE ( DBGUNT, * ) 'Swale ponding, HRU:', i, &
     &                  ' gravity reservoir is 3*sat_threshold', Slow_stor(i), Sat_threshold(i)
                CALL print_date(DBGUNT)
              ENDIF
            ENDIF
          ENDIF
          Ssres_stor(i) = Slow_stor(i)
        ENDIF

        IF ( Soil_lower_stor_max(i)>0.0 ) Soil_lower_ratio(i) = Soil_lower(i)/Soil_lower_stor_max(i)
        Ssres_in(i) = Soil_to_ssr(i) + Pref_flow_infil(i)
        Basin_ssin = Basin_ssin + DBLE( Ssres_in(i)*harea )
        Basin_ssstor = Basin_ssstor + DBLE( Ssres_stor(i)*harea )
        Basin_slstor = Basin_slstor + DBLE( Slow_stor(i)*harea )
        Soil_moist_tot(i) = Ssres_stor(i) + Soil_moist(i)*perv_frac
        IF ( ag_on_flag==ACTIVE ) Soil_moist_tot(i) = Soil_moist_tot(i) + Ag_soil_moist(i)*agfrac
        Basin_soil_moist_tot = Basin_soil_moist_tot + DBLE( Soil_moist_tot(i)*harea )
        IF ( Soil_zone_max(i)>0.0 ) Basin_sz_stor_frac = Basin_sz_stor_frac + Soil_moist_tot(i)/Soil_zone_max(i)*harea
        IF ( perv_on_flag==ACTIVE ) THEN
          Basin_cpr_stor_frac = Basin_cpr_stor_frac + Soil_moist(i)/Soil_moist_max(i)*perv_area
          Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac + Soil_rechr(i)/Soil_rechr_max(i)*perv_area
          Basin_soil_lower_stor_frac = Basin_soil_lower_stor_frac + Soil_lower_ratio(i)*perv_area
        ENDIF
        IF ( Pref_flow_thrsh(i)>0.0 ) Basin_gvr_stor_frac = Basin_gvr_stor_frac + Slow_stor(i)/Pref_flow_thrsh(i)*harea
        Recharge(i) = Soil_to_gw(i) + Ssr_to_gw(i)
        IF ( Dprst_flag==1 ) Recharge(i) = Recharge(i) + SNGL( Dprst_seep_hru(i) )
        Basin_recharge = Basin_recharge + DBLE( Recharge(i)*harea )
        Grav_dunnian_flow(i) = dunnianflw_gvr
        Unused_potet(i) = Potet(i) - Hru_actet(i)
        IF ( ag_on_flag==ACTIVE ) THEN
          Ag_soil_lower(i) = Ag_soil_moist(i) - Ag_soil_rechr(i)
          Basin_ag_soil_moist = Basin_ag_soil_moist + DBLE( Ag_soil_moist(i)*Ag_area(i) )
          Basin_ag_soil_rechr = Basin_ag_soil_rechr + DBLE( Ag_soil_rechr(i)*Ag_area(i) )
          Basin_ag_actet = Basin_ag_actet + DBLE( Ag_actet(i)*Ag_area(i) )
          IF ( Iter_aet_flag==ACTIVE .AND.Transp_on(i)==ACTIVE ) THEN
          !IF ( Iter_aet_flag==ACTIVE ) THEN
            !agriculture_external(i)
            !IF ( Unused_potet(i)>0.0 ) THEN
            IF ( unsatisfied_et>soilzone_aet_converge ) THEN
              IF ( Ag_soilwater_deficit(i)>Ag_soilwater_deficit_min(i) ) &
                  Ag_irrigation_add(i) = Ag_irrigation_add(i) + unsatisfied_et
              keep_iterating = ACTIVE
              add_estimated_irrigation = ACTIVE
              num_hrus_ag_iter = num_hrus_ag_iter + 1
              IF ( unsatisfied_et>unsatisfied_big ) unsatisfied_big = unsatisfied_et
            ENDIF
          ENDIF
        ENDIF
        Basin_actet = Basin_actet + DBLE( Hru_actet(i)*harea )
        Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Pkwater_equiv(i)
        IF ( Dprst_flag==ACTIVE ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)

! ***************************************
      ENDDO ! end HRU loop
! ***************************************

      Soil_iter = Soil_iter + 1
      IF ( Iter_aet_flag==OFF .OR. Basin_transp_on==OFF ) keep_iterating = OFF
      IF ( Soil_iter>max_soilzone_ag_iter .OR. add_estimated_irrigation==OFF ) keep_iterating = OFF
!      if (ag_actet(i)>potet(i)) print *, 'ag_actet>potet', ag_actet(i), potet(i)
!      if (unused_potet(i)<0.00) print *, 'unused', unused_potet(i), potet(i)
!      if ( hru_actet(i)>potet(i)) print *, 'hru_actet', hru_actet(i), potet(i)

! ***************************************
      ENDDO ! end iteration while loop
! ***************************************

      Soil_iter = Soil_iter - 1
      IF ( Iter_aet_flag==ACTIVE ) THEN
        Ag_irrigation_add_vol = Ag_irrigation_add*Ag_area
        Basin_ag_soil_moist = Basin_ag_soil_moist*Basin_area_inv
        Basin_ag_soil_rechr = Basin_ag_soil_rechr*Basin_area_inv
        Basin_ag_waterin = Basin_ag_waterin*Basin_area_inv
      ENDIF
!      IF ( num_hrus_ag_iter>0 ) print '(2(A,I0))', 'number of hrus still iterating on AET: ', &
!     &     num_hrus_ag_iter
!      if ( soil_iter==max_soilzone_ag_iter ) iter_nonconverge = iter_nonconverge + 1
!      print *, 'iterations: ', Soil_iter, '; nonconverged', iter_nonconverge
      total_iters = total_iters + Soil_iter
!      print *, Nowtime, unsatisfied_big, unsatisfied_big/basin_potet, total_iters

      Basin_actet = Basin_actet*Basin_area_inv
      Basin_perv_et = Basin_perv_et*Basin_area_inv
      Basin_swale_et = Basin_swale_et*Basin_area_inv
      Basin_soil_rechr = Basin_soil_rechr*Basin_area_inv
      Basin_soil_to_gw = Basin_soil_to_gw*Basin_area_inv
      Basin_soil_moist = Basin_soil_moist*Basin_area_inv
      Basin_soil_moist_tot = Basin_soil_moist_tot*Basin_area_inv
      IF ( Nlake>0 ) THEN
        Basin_lakeevap = Basin_lakeevap*Basin_area_inv
        Basin_lakeprecip = Basin_lakeprecip*Basin_area_inv
        Basin_lakeinsz = Basin_lakeinsz*Basin_area_inv
        Basin_lake_stor = Basin_lake_stor + Basin_lakeprecip - Basin_lakeevap
      ENDIF
      IF ( Pref_flag==ACTIVE ) THEN
        Basin_pref_stor = Basin_pref_stor*Basin_area_inv
        Basin_pref_flow_infil = Basin_pref_flow_infil*Basin_area_inv
        Basin_prefflow = Basin_prefflow*Basin_area_inv
        Basin_dunnian_pfr = Basin_dunnian_pfr*Basin_area_inv
        Basin_pfr_stor_frac = Basin_pfr_stor_frac*Basin_area_inv
      ENDIF
      Basin_dunnian_gvr = Basin_dunnian_gvr*Basin_area_inv
      Basin_ssstor = Basin_ssstor*Basin_area_inv
      Basin_ssflow = Basin_ssflow*Basin_area_inv
      Basin_interflow_max = Basin_interflow_max*Basin_area_inv
      Basin_sz2gw = Basin_sz2gw*Basin_area_inv
      Basin_ssin = Basin_ssin*Basin_area_inv
      Basin_slstor = Basin_slstor*Basin_area_inv
      Basin_sroff = Basin_sroff*Basin_area_inv
      Basin_dunnian = Basin_dunnian*Basin_area_inv
      Basin_sm2gvr = Basin_sm2gvr*Basin_area_inv
      Basin_sm2gvr_max = Basin_sm2gvr_max*Basin_area_inv
      Basin_capwaterin = Basin_capwaterin*Basin_area_inv
      Basin_cap_infil_tot = Basin_cap_infil_tot*Basin_area_inv
      Basin_cap_up_max = Basin_cap_up_max*Basin_area_inv
      Basin_dninterflow = Basin_dninterflow*Basin_area_inv
      Basin_dndunnianflow = Basin_dndunnianflow*Basin_area_inv
      Basin_dncascadeflow = Basin_dncascadeflow*Basin_area_inv
      Basin_gvr2pfr = Basin_gvr2pfr*Basin_area_inv
      Basin_slowflow = Basin_slowflow*Basin_area_inv
      Basin_recharge = Basin_recharge*Basin_area_inv
      Basin_cpr_stor_frac = Basin_cpr_stor_frac*Basin_area_inv
      Basin_gvr_stor_frac = Basin_gvr_stor_frac*Basin_area_inv
      Basin_sz_stor_frac = Basin_sz_stor_frac*Basin_area_inv
      Basin_soil_lower_stor_frac = Basin_soil_lower_stor_frac*Basin_area_inv
      Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac*Basin_area_inv
      IF ( update_potet==ACTIVE ) THEN
        Basin_potet = 0.0D0
        DO k = 1, Active_hrus
          i = Hru_route_order(k)
          Basin_potet = Basin_potet + DBLE( Potet(i)*Hru_area(i) )
        ENDDO
        Basin_potet = Basin_potet*Basin_area_inv
      ENDIF

      END FUNCTION szrun_ag

!***********************************************************************
!     soilzone_restart_ag - write or read soilzone_ag restart file
!***********************************************************************
      SUBROUTINE soilzone_restart_ag(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT, ACTIVE
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_SOILZONE
      USE PRMS_SOILZONE_AG
      use prms_utils, only: check_restart
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Local Variable
      CHARACTER(LEN=8) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_soil_rechr, Basin_slstor, Basin_soil_moist_tot, Basin_pref_stor
        WRITE ( Restart_outunit ) Pref_flow_stor
        WRITE ( Restart_outunit ) Ag_soil_lower
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_soil_rechr, Basin_slstor, Basin_soil_moist_tot, Basin_pref_stor
        READ ( Restart_inunit ) Pref_flow_stor
        READ ( Restart_inunit ) Ag_soil_lower
      ENDIF
      END SUBROUTINE soilzone_restart_ag
