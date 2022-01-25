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
      MODULE PRMS_SOILZONE
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ERROR_dim, ACTIVE, OFF, DEBUG_WB, NEARZERO, &
     &    LAND, LAKE, SWALE, INACTIVE, GLACIER, BARESOIL, DEBUG_less, &
     &    ERROR_soilzone, SAND, CLAY, LOAM, CASCADE_OFF, ERROR_param, MONTHS_PER_YEAR
      USE PRMS_MODULE, ONLY: Model, Nhru, Nssr, Nsegment, Nlake, Nhrucell, Print_debug, Dprst_flag, &
     &    Init_vars_from_file, Cascade_flag, GSFLOW_flag, Parameter_check_flag, Inputerror_flag, &
     &    Kkiter, Soilzone_aet_flag, Frozen_flag, Soilzone_add_water_use, Call_cascade
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Soilzone Computations'
      character(len=8), parameter :: MODNAME = 'soilzone'
      character(len=*), parameter :: Version_soilzone = '2021-01-11'
      INTEGER, SAVE :: DBGUNT
      INTEGER, SAVE :: Max_gvrs, Et_type, Pref_flag
      REAL, SAVE, ALLOCATABLE :: Gvr2pfr(:), Swale_limit(:)
      REAL, SAVE, ALLOCATABLE :: Soil_lower_stor_max(:)
      REAL, SAVE, ALLOCATABLE :: Soil_moist_ante(:), Ssres_stor_ante(:)
      REAL, SAVE, ALLOCATABLE :: Grav_dunnian_flow(:), Pfr_dunnian_flow(:)
      DOUBLE PRECISION, SAVE :: Last_soil_moist, Last_ssstor
!   GSFLOW variables
      INTEGER, SAVE, ALLOCATABLE :: Hru_gvr_count(:), Hru_gvr_index(:, :), Hrucheck(:)
      REAL, SAVE, ALLOCATABLE :: Replenish_frac(:)
      REAL, SAVE, ALLOCATABLE :: It0_soil_rechr(:), It0_soil_moist(:)
      REAL, SAVE, ALLOCATABLE :: It0_pref_flow_stor(:), It0_ssres_stor(:)
      REAL, SAVE, ALLOCATABLE :: It0_gravity_stor_res(:), It0_sroff(:)
      REAL, SAVE, ALLOCATABLE :: It0_slow_stor(:), It0_potet(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: It0_strm_seg_in(:)
      DOUBLE PRECISION, SAVE :: It0_basin_soil_moist, It0_basin_ssstor, Basin_sz_gwin
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gvr_hru_pct_adjusted(:)
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_sz2gw, Basin_cap_infil_tot
      DOUBLE PRECISION, SAVE :: Basin_interflow_max, Basin_sm2gvr_max ! this is the same as basin_sm2gvr
      DOUBLE PRECISION, SAVE :: Basin_soil_rechr, Basin_dunnian_gvr
      DOUBLE PRECISION, SAVE :: Basin_recharge, Basin_pref_flow_infil
      DOUBLE PRECISION, SAVE :: Basin_ssin, Basin_dunnian_pfr
      DOUBLE PRECISION, SAVE :: Basin_sm2gvr, Basin_dninterflow
      DOUBLE PRECISION, SAVE :: Basin_dncascadeflow, Basin_dndunnianflow
      DOUBLE PRECISION, SAVE :: Basin_capwaterin, Basin_dunnian
      DOUBLE PRECISION, SAVE :: Basin_gvr2pfr, Basin_slowflow
      DOUBLE PRECISION, SAVE :: Basin_pref_stor, Basin_slstor, Basin_prefflow
      DOUBLE PRECISION, SAVE :: Basin_lakeinsz, Basin_lakeprecip
      DOUBLE PRECISION, SAVE :: Basin_cap_up_max
      DOUBLE PRECISION, SAVE :: Basin_soil_moist_tot
      DOUBLE PRECISION, SAVE :: Basin_soil_lower_stor_frac, Basin_soil_rechr_stor_frac, Basin_sz_stor_frac
      DOUBLE PRECISION, SAVE :: Basin_cpr_stor_frac, Basin_gvr_stor_frac, Basin_pfr_stor_frac
      REAL, SAVE, ALLOCATABLE :: Perv_actet(:), Pref_flow_thrsh(:)
      REAL, SAVE, ALLOCATABLE :: Soil_moist_tot(:), Recharge(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Upslope_interflow(:), Upslope_dunnianflow(:), Lakein_sz(:)
      REAL, SAVE, ALLOCATABLE :: Dunnian_flow(:), Cap_infil_tot(:)
      REAL, SAVE, ALLOCATABLE :: Pref_flow_stor(:), Pref_flow(:)
      REAL, SAVE, ALLOCATABLE :: Pref_flow_infil(:), Pref_flow_in(:)
      REAL, SAVE, ALLOCATABLE :: Hru_sz_cascadeflow(:), Swale_actet(:)
      REAL, SAVE, ALLOCATABLE :: Pref_flow_max(:), Snow_free(:)
      REAL, SAVE, ALLOCATABLE :: Cap_waterin(:), Soil_lower(:), Soil_zone_max(:)
      REAL, SAVE, ALLOCATABLE :: Potet_lower(:), Potet_rechr(:), Soil_lower_ratio(:)
      REAL, SAVE, ALLOCATABLE :: Unused_potet(:)
      INTEGER, SAVE, ALLOCATABLE :: Soil_saturated(:)
!      REAL, SAVE, ALLOCATABLE :: Cascade_interflow(:), Cascade_dunnianflow(:), Interflow_max(:)
!      REAL, SAVE, ALLOCATABLE :: Cpr_stor_frac(:), Pfr_stor_frac(:), Gvr_stor_frac(:), Soil_moist_frac(:)
!      REAL, SAVE, ALLOCATABLE :: Soil_rechr_ratio(:), Snowevap_aet_frac(:), Perv_avail_et(:), Cap_upflow_max(:)
!   GSFLOW Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_gvr2sm
      REAL, SAVE, ALLOCATABLE :: Sm2gw_grav(:), Gw2sm_grav(:)
      REAL, SAVE, ALLOCATABLE :: Gravity_stor_res(:), Gvr2sm(:), Grav_gwin(:)
!   Declared Parameters
      INTEGER, SAVE, ALLOCATABLE :: Soil_type(:), Gvr_hru_id(:)
      REAL, SAVE, ALLOCATABLE :: Pref_flow_den(:)
      REAL, SAVE, ALLOCATABLE :: Fastcoef_lin(:), Fastcoef_sq(:)
      REAL, SAVE, ALLOCATABLE :: Slowcoef_lin(:), Slowcoef_sq(:)
      REAL, SAVE, ALLOCATABLE :: Ssr2gw_rate(:), Ssr2gw_exp(:)
      REAL, SAVE, ALLOCATABLE :: Soil2gw_max(:)
      REAL, SAVE, ALLOCATABLE :: Lake_evap_adj(:, :)
      END MODULE PRMS_SOILZONE

!***********************************************************************
!     Main soilzone routine
!***********************************************************************
      INTEGER FUNCTION soilzone()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
! Functions
      INTEGER, EXTERNAL :: szdecl, szinit, szrun
      EXTERNAL :: soilzone_restart
!***********************************************************************
      soilzone = 0

      IF ( Process_flag==RUN ) THEN
        soilzone = szrun()
      ELSEIF ( Process_flag==DECL ) THEN
        soilzone = szdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL soilzone_restart(READ_INIT)
        soilzone = szinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL soilzone_restart(SAVE_INIT)
      ENDIF

      END FUNCTION soilzone

!***********************************************************************
!     szdecl - set up parameters for soil zone computations
!   Declared Parameters
!     sat_threshold, ssstor_init_frac fastcoef_lin, fastcoef_sq
!     ssr2gw_rate, ssr2gw_exp, soil2gw_max, soil_type
!     soil_rechr_max_frac, soil_rechr_init_frac, soil_moist_max, soil_moist_init_frac
!     pref_flow_den, slowcoef_lin, cov_type
!     hru_area, slowcoef_sq, gvr_hru_id
!***********************************************************************
      INTEGER FUNCTION szdecl()
      USE PRMS_SOILZONE
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar, getdim
      EXTERNAL :: read_error, print_module, PRMS_open_module_file, error_stop
!***********************************************************************
      szdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_soilzone)

! Declare Variables
      IF ( declvar(MODNAME, 'basin_capwaterin', 'one', 1, 'double', &
     &     'Basin area-weighted average infiltration,'// &
     &     ' cascading interflow and Dunnian flow added to capillary reservoir storage', &
     &     'inches', Basin_capwaterin)/=0 ) CALL read_error(3, 'basin_capwaterin')

      IF ( declvar(MODNAME, 'basin_cap_infil_tot', 'one', 1, 'double', &
     &     'Basin area-weighted average infiltration with cascading flow into capillary reservoirs', &
     &     'inches', Basin_cap_infil_tot)/=0 ) CALL read_error(3, 'basin_cap_infil_tot')

      IF ( declvar(MODNAME, 'basin_cap_up_max', 'one', 1, 'double', &
     &     'Basin area-weighted average maximum cascade flow that flows to capillary reservoirs', &
     &     'inches', Basin_cap_up_max)/=0 ) CALL read_error(3, 'basin_cap_up_max')

      IF ( declvar(MODNAME, 'basin_pref_flow_infil', 'one', 1, 'double', &
     &     'Basin area-weighted average infiltration to preferential-flow reservoir storage', &
     &     'inches', Basin_pref_flow_infil)/=0 ) CALL read_error(3, 'basin_pref_flow_infil')

      IF ( declvar(MODNAME, 'basin_dunnian_pfr', 'one', 1, 'double', &
     &     'Basin area-weighted average excess infiltration to'// &
     &     ' preferential-flow reservoirs from variable infil', &
     &     'inches', Basin_dunnian_pfr)/=0 ) CALL read_error(3, 'basin_dunnian_pfr')

      IF ( declvar(MODNAME, 'basin_dunnian_gvr', 'one', 1, 'double', &
     &     'Basin area-weighted average excess flow to preferential'// &
     &     '-flow reservoirs from gravity reservoirs', &
     &     'inches', Basin_dunnian_gvr)/=0 ) CALL read_error(3, 'basin_dunnian_gvr')

      ALLOCATE ( Cap_infil_tot(Nhru) )
      IF ( declvar(MODNAME, 'cap_infil_tot', 'nhru', Nhru, 'real', &
     &     'Infiltration and cascading interflow and Dunnian'// &
     &     ' flow added to capillary reservoir storage for each HRU', &
     &     'inches', Cap_infil_tot)/=0 ) CALL read_error(3, 'cap_infil_tot')

      IF ( declvar(MODNAME, 'basin_soil_moist_tot', 'one', 1, 'double', &
     &     'Basin area-weighted average total soil-zone water storage', &
     &     'inches', Basin_soil_moist_tot)/=0 ) CALL read_error(3, 'basin_soil_moist_tot')

      ALLOCATE ( Soil_moist_tot(Nhru) )
      IF ( declvar(MODNAME, 'soil_moist_tot', 'nhru', Nhru, 'real', &
     &     'Total soil-zone water storage (soil_moist + ssres_stor)', &
     &     'inches', Soil_moist_tot)/=0 ) CALL read_error(3, 'soil_moist_tot')

      IF ( declvar(MODNAME, 'basin_cpr_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of capillary reservoir storage of the maximum storage', &
     &     'decimal fraction', Basin_cpr_stor_frac)/=0 ) CALL read_error(3, 'basin_cpr_stor_frac')

      IF ( declvar(MODNAME, 'basin_gvr_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of gravity reservoir storage of the maximum storage', &
     &     'decimal fraction', Basin_gvr_stor_frac)/=0 ) CALL read_error(3, 'basin_gvr_stor_frac')

      IF ( declvar(MODNAME, 'basin_pfr_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of preferential-flow reservoir storage of the maximum storage', &
     &     'decimal fraction', Basin_pfr_stor_frac)/=0 ) CALL read_error(3, 'basin_pfr_stor_frac')

      IF ( declvar(MODNAME, 'basin_soil_lower_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of soil lower zone storage of the maximum storage', &
     &     'decimal fraction', Basin_soil_lower_stor_frac)/=0 ) CALL read_error(3, 'basin_soil_lower_stor_frac')

      IF ( declvar(MODNAME, 'basin_soil_rechr_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of soil recharge zone storage of the maximum storage', &
     &     'decimal fraction', Basin_soil_rechr_stor_frac)/=0 ) CALL read_error(3, 'basin_soil_rechr_stor_frac')

      IF ( declvar(MODNAME, 'basin_sz_stor_frac', 'one', 1, 'double', &
     &     'Basin area-weighted average fraction of soil zone storage of the maximum storage', &
     &     'decimal fraction', Basin_sz_stor_frac)/=0 ) CALL read_error(3, 'basin_sz_stor_frac')

!      ALLOCATE ( Cpr_stor_frac(Nhru) )
!      IF ( declvar(MODNAME, 'cpr_stor_frac', 'nhru', Nhru, 'real', &
!     &     'Fraction of capillary reservoir storage of the maximum storage for each HRU', &
!     &     'decimal fraction', Cpr_stor_frac)/=0 ) CALL read_error(3, 'cpr_stor_frac')

!      ALLOCATE ( Pfr_stor_frac(Nhru) )
!      IF ( declvar(MODNAME, 'pfr_stor_frac', 'nhru', Nhru, 'real', &
!     &     'Fraction of preferential flow reservoir storage of the maximum storage for each HRU', &
!     &     'decimal fraction', Pfr_stor_frac)/=0 ) CALL read_error(3, 'pfr_stor_frac')

!      ALLOCATE ( Gvr_stor_frac(Nhru) )
!      IF ( declvar(MODNAME, 'gvr_stor_frac', 'nhru', Nhru, 'real', &
!     &     'Fraction of gravity reservoir storage of the maximum storage for each HRU', &
!     &     'decimal fraction', Gvr_stor_frac)/=0 ) CALL read_error(3, 'gvr_stor_frac')

!      ALLOCATE ( Soil_moist_frac(Nhru) )
!      IF ( declvar(MODNAME, 'soil_moist_frac', 'nhru', Nhru, 'real', &
!     &     'Fraction of soil zone storage of the maximum storage for each HRU', &
!     &     'decimal fraction', Soil_moist_frac)/=0 ) CALL read_error(3, 'soil_moist_frac')

      IF ( declvar(MODNAME, 'basin_sm2gvr', 'one', 1, 'double', &
     &     'Basin area-weighted average excess flow from'// &
     &     ' capillary reservoirs to gravity reservoir storage', &
     &     'inches', Basin_sm2gvr)/=0 ) CALL read_error(3, 'basin_sm2gvr')

      IF ( declvar(MODNAME, 'basin_gvr2pfr', 'one', 1, 'double', &
     &     'Basin area-weighted average excess flow to'// &
     &     ' preferential-flow reservoir storage from gravity reservoirs', &
     &     'inches', Basin_gvr2pfr)/=0 ) CALL read_error(3, 'basin_gvr2pfr')

      IF ( declvar(MODNAME, 'basin_slowflow', 'one', 1, 'double', &
     &     'Basin area-weighted average interflow from gravity reservoirs to the stream network', &
     &     'inches', Basin_slowflow)/=0 ) CALL read_error(3, 'basin_slowflow')

      IF ( declvar(MODNAME, 'basin_prefflow', 'one', 1, 'double', &
     &     'Basin area-weighted average interflow from'// &
     &     ' preferential-flow reservoirs to the stream network', &
     &     'inches', Basin_prefflow)/=0 ) CALL read_error(3, 'basin_prefflow')

      IF ( declvar(MODNAME, 'basin_slstor', 'one', 1, 'double', &
     &     'Basin area-weighted average storage of gravity reservoirs', &
     &     'inches', Basin_slstor)/=0 ) CALL read_error(3, 'basin_slstor')

      ALLOCATE ( Dunnian_flow(Nhru) )
      IF ( declvar(MODNAME, 'dunnian_flow', 'nhru', Nhru, 'real', &
     &     'Dunnian surface runoff that flows to the stream network for each HRU', &
     &     'inches', Dunnian_flow)/=0 ) CALL read_error(3, 'dunnian_flow')

      IF ( declvar(MODNAME, 'basin_dunnian', 'one', 1, 'double', &
     &     'Basin area-weighted average Dunnian surface runoff that flows to the stream network', &
     &     'inches', Basin_dunnian)/=0 ) CALL read_error(3, 'basin_dunnian')

      IF ( declvar(MODNAME, 'basin_soil_rechr', 'one', 1, 'double', &
     &     'Basin area-weighted average storage for recharge zone;'// &
     &     ' upper portion of capillary reservoir where both'// &
     &     ' evaporation and transpiration occurs', &
     &     'inches', Basin_soil_rechr)/=0 ) CALL read_error(3, 'basin_soil_rechr')

      IF ( declvar(MODNAME, 'basin_sz2gw', 'one', 1, 'double', &
     &     'Basin area-weighted average drainage from gravity reservoirs to GWRs', &
     &     'inches', Basin_sz2gw)/=0 ) CALL read_error(3, 'basin_sz2gw')

      ALLOCATE ( Pref_flow_in(Nhru) )
      IF ( declvar('soilzone', 'pref_flow_in', 'nhru', Nhru, 'real', &
     &     'Infiltration and flow from gravity reservoir to the preferential-flow reservoir', &
     &     'inches', Pref_flow_in)/=0 ) CALL read_error(3, 'pref_flow_in')

      IF ( declvar(MODNAME, 'basin_sm2gvr_maxin', 'one', 1, 'double', &
     &     'Basin area-weighted average maximum excess flow from'// &
     &     ' capillary reservoirs that flows to gravity reservoirs', &
     &     'inches', Basin_sm2gvr_max)/=0 ) CALL read_error(3, 'basin_sm2gvr_max')

      IF ( declvar(MODNAME, 'basin_interflow_max', 'one', 1, 'double', &
     &     'Basin area-weighted average maximum interflow that flows from gravity reservoirs', &
     &     'inches', Basin_interflow_max)/=0 ) CALL read_error(3, 'basin_interflow_max')

      ALLOCATE ( Perv_actet(Nhru) )
      IF ( declvar(MODNAME, 'perv_actet', 'nhru', Nhru, 'real', &
     &     'Actual ET from the capillary reservoir of each HRU', &
     &     'inches', Perv_actet)/=0 ) CALL read_error(3, 'perv_actet')

!      ALLOCATE ( Perv_avail_et(Nhru) )
!      IF ( declvar(MODNAME, 'perv_avail_et', 'nhru', Nhru, 'real', &
!     &     'Unsatisfied ET available to the capillary reservoir of each HRU', &
!     &     'inches', Perv_avail_et)/=0 ) CALL read_error(3, 'perv_avail_et')

      ! added to be compatible with ssflow_prms
      IF ( declvar(MODNAME, 'basin_ssin', 'one', 1, 'double', &
     &     'Basin area-weighted average inflow to gravity and preferential-flow reservoir storage', &
     &     'inches', Basin_ssin)/=0 ) CALL read_error(3, 'basin_ssin')

!      ALLOCATE ( Interflow_max(Nhru) )
!      IF ( declvar(MODNAME, 'interflow_max', 'nhru', Nhru, 'real', &
!     &     'Maximum interflow for each HRU', &
!     &     'inches', Interflow_max)/=0 ) CALL read_error(3, 'interflow_max')

      IF ( Cascade_flag>CASCADE_OFF .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'basin_dndunnianflow', 'one', 1, 'double', &
     &       'Basin area-weighted average cascading Dunnian flow', &
     &       'inches', Basin_dndunnianflow)/=0 ) CALL read_error(3, 'basin_dndunnianflow')

        IF ( declvar(MODNAME, 'basin_dninterflow', 'one', 1, 'double', &
     &       'Basin area-weighted average cascading interflow', &
     &       'inches', Basin_dninterflow)/=0 ) CALL read_error(3, 'basin_dninterflow')

        IF ( declvar(MODNAME, 'basin_dncascadeflow', 'one', 1, 'double', &
     &       'Basin area-weighted average cascading interflow and Dunnian surface runoff', &
     &       'inches', Basin_dncascadeflow)/=0 ) CALL read_error(3, 'basin_dncascadeflow')

        ALLOCATE ( Upslope_interflow(Nhru) )
        IF ( declvar(MODNAME, 'upslope_interflow', 'nhru', Nhru, 'double', &
     &       'Cascading interflow runoff that flows to'// &
     &       ' the capillary reservoir of each downslope HRU for each upslope HRU', &
     &       'inches', Upslope_interflow)/=0 ) CALL read_error(3, 'upslope_interflow')

        ALLOCATE ( Upslope_dunnianflow(Nhru) )
        IF ( declvar(MODNAME, 'upslope_dunnianflow', 'nhru', Nhru, 'double', &
     &       'Cascading Dunnian surface runoff that'// &
     &       ' flows to the capillary reservoir of each downslope HRU for each upslope HRU', &
     &       'inches', Upslope_dunnianflow)/=0 ) CALL read_error(3, 'upslope_dunnianflow')

        ALLOCATE ( Hru_sz_cascadeflow(Nhru) )
        IF ( declvar(MODNAME, 'hru_sz_cascadeflow', 'nhru', Nhru, 'real', &
     &       'Cascading interflow and Dunnian surface runoff from each HRU', &
     &       'inches', Hru_sz_cascadeflow)/=0 ) CALL read_error(3, 'hru_sz_cascadeflow')

!        ALLOCATE ( Cap_upflow_max(Nhru) )
!        IF ( declvar(MODNAME, 'cap_upflow_max', 'nhru', Nhru, 'real', &
!     &       'Maximum infiltration and any cascading interflow and'// &
!     &       ' Dunnian surface runoff that can be added to capillary reservoir storage for each HRU', &
!     &       'inches', Cap_upflow_max)/=0 ) CALL read_error(3, 'cap_upflow_max')

!        ALLOCATE ( Cascade_interflow(Nhru) )
!        IF ( declvar(MODNAME, 'cascade_interflow', 'nhru', Nhru, 'real', &
!     &       'Cascading interflow for each HRU', &
!     &       'inches', Cascade_interflow)/=0 ) CALL read_error(3, 'cascade_interflow')

!        ALLOCATE ( Cascade_dunnianflow(Nhru) )
!        IF ( declvar(MODNAME, 'cascade_dunnianflow', 'nhru', Nhru, 'real', &
!     &       'Cascading Dunnian flow for each HRU', &
!     &       'inches', Cascade_dunnianflow)/=0 ) CALL read_error(3, 'cascade_dunnianflow')

        IF ( Nlake>0 ) THEN
          ALLOCATE ( Lakein_sz(Nhru) )
          IF ( declvar(MODNAME, 'lakein_sz', 'nhru', Nhru, 'double', &
     &         'Cascading interflow and Dunnian surface runoff to lake HRUs from each upslope HRU', &
     &         'inches', Lakein_sz)/=0 ) CALL read_error(3, 'lakein_sz')

          IF ( declvar(MODNAME, 'basin_lakeinsz', 'one', 1, 'double', &
     &         'Basin area-weighted average lake inflow from land HRUs', &
     &         'inches', Basin_lakeinsz)/=0 ) CALL read_error(3, 'basin_lakeinsz')
        ENDIF
      ENDIF

      IF ( declvar(MODNAME, 'basin_pref_stor', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in preferential-flow reservoirs', &
     &     'inches', Basin_pref_stor)/=0 ) CALL read_error(3, 'basin_pref_stor')

      ALLOCATE ( Pref_flow_infil(Nhru) )
      IF ( declvar(MODNAME, 'pref_flow_infil', 'nhru', Nhru, 'real', &
     &     'Infiltration to the preferential-flow reservoir storage for each HRU', &
     &     'inches', Pref_flow_infil)/=0 ) CALL read_error(3, 'pref_flow_infil')

      ALLOCATE ( Pref_flow_stor(Nhru) )
      IF ( declvar(MODNAME, 'pref_flow_stor', 'nhru', Nhru, 'real', &
     &     'Storage in preferential-flow reservoir for each HRU', &
     &     'inches', Pref_flow_stor)/=0 ) CALL read_error(3, 'pref_flow_stor')

      ALLOCATE ( Pref_flow(Nhru) )
      IF ( declvar(MODNAME, 'pref_flow', 'nhru', Nhru, 'real', &
     &     'Interflow from the preferential-flow reservoir that'// &
     &     ' flows to the stream network for each HRU', &
     &     'inches', Pref_flow)/=0 ) CALL read_error(3, 'pref_flow')

      ALLOCATE ( Pref_flow_thrsh(Nhru) )
      IF ( declvar(MODNAME, 'pref_flow_thrsh', 'nhru', Nhru, 'real', &
     &     'Soil storage threshold defining storage between field'// &
     &     ' capacity and maximum soil saturation minus preferential-flow storage', &
     &     'inches', Pref_flow_thrsh)/=0 ) CALL read_error(3, 'pref_flow_thrsh')

      ALLOCATE ( Pref_flow_max(Nhru) )
      IF ( declvar(MODNAME, 'pref_flow_max', 'nhru', Nhru, 'real', &
     &     'Maximum storage of the preferential-flow reservoir for each HRU', &
     &     'inches', Pref_flow_max)/=0 ) CALL read_error(3, 'pref_flow_max')

      ALLOCATE ( Soil_zone_max(Nhru) )
!      IF ( declvar(MODNAME, 'soil_zone_max', 'nhru', Nhru, 'real', &
!     &     'Maximum storage of all soil zone reservoirs', &
!     &     'inches', Soil_zone_max)/=0 ) CALL read_error(3, 'soil_zone_max')

      IF ( declvar(MODNAME, 'basin_lakeprecip', 'one', 1, 'double', &
     &     'Basin area-weighted average precipitation on lake HRUs', &
     &     'inches', Basin_lakeprecip)/=0 ) CALL read_error(3, 'basin_lakeprecip')

      ALLOCATE ( Swale_actet(Nhru) )
      IF ( declvar(MODNAME, 'swale_actet', 'nhru', Nhru, 'real', &
     &     'Evaporation from the gravity and preferential-flow reservoirs that exceeds sat_threshold', &
     &     'inches', Swale_actet)/=0 ) CALL read_error(3, 'swale_actet')

      IF ( declvar(MODNAME, 'basin_recharge', 'one', 1, 'double', &
     &     'Basin area-weighted average recharge to GWRs', &
     &     'inches', Basin_recharge)/=0 ) CALL read_error(3, 'basin_recharge')

      ALLOCATE ( Recharge(Nhru) )
      IF ( declvar(MODNAME, 'recharge', 'nhru', Nhru, 'real', &
     &     'Recharge to the associated GWR as sum of soil_to_gw and ssr_to_gw for each HRU', &
     &     'inches', Recharge)/=0 ) CALL read_error(3, 'recharge')

      ALLOCATE ( Cap_waterin(Nhru) )
      IF ( declvar(MODNAME, 'cap_waterin', 'nhru', Nhru, 'real', &
     &     'Infiltration and any cascading interflow and'// &
     &     ' Dunnian surface runoff added to capillary reservoir storage for each HRU', &
     &     'inches', Cap_waterin)/=0 ) CALL read_error(3, 'cap_waterin')

      ALLOCATE ( Soil_lower(Nhru) )
      IF ( declvar(MODNAME, 'soil_lower', 'nhru', Nhru, 'real', &
     &     'Storage in the lower zone of the capillary'// &
     &     ' reservoir that is only available for transpiration for each HRU', &
     &     'inches', Soil_lower)/=0 ) CALL read_error(3, 'soil_lower')

      ALLOCATE ( Potet_lower(Nhru) )
      IF ( declvar(MODNAME, 'potet_lower', 'nhru', Nhru, 'real', &
     &     'Potential ET in the lower zone of the capillary reservoir for each HRU', &
     &     'inches', Potet_lower)/=0 ) CALL read_error(3, 'potet_lower')

      ALLOCATE ( Potet_rechr(Nhru) )
      IF ( declvar(MODNAME, 'potet_rechr', 'nhru', Nhru, 'real', &
     &     'Potential ET in the recharge zone of the capillary reservoir for each HRU', &
     &     'inches', Potet_rechr)/=0 ) CALL read_error(3, 'potet_rechr')

      ALLOCATE ( Soil_lower_ratio(Nhru), Soil_lower_stor_max(Nhru) )
      IF ( declvar(MODNAME, 'soil_lower_ratio', 'nhru', Nhru, 'real', &
     &     'Water content ratio in the lower zone of the capillary reservoir for each HRU', &
     &     'decimal fraction', Soil_lower_ratio)/=0 ) CALL read_error(3, 'soil_lower_ratio')

!      ALLOCATE ( Soil_rechr_ratio(Nhru) )
!      IF ( declvar(MODNAME, 'soil_rechr_ratio', 'nhru', Nhru, 'real', &
!     &     'Water content ratio in the recharge zone of the capillary reservoir for each HRU', &
!     &     'decimal fraction', Soil_rechr_ratio)/=0 ) CALL read_error(3, 'soil_rechr_ratio')

      ALLOCATE ( Snow_free(Nhru) )
      IF ( declvar(MODNAME, 'snow_free', 'nhru', Nhru, 'real', &
     &     'Fraction of snow-free surface for each HRU', &
     &     'decimal fraction', Snow_free)/=0 ) CALL read_error(3, 'snow_free')

      ALLOCATE ( Unused_potet(Nhru) )
      IF ( declvar(MODNAME, 'unused_potet', 'nhru', Nhru, 'real', &
     &     'Unsatisfied potential evapotranspiration', &
     &     'inches', Unused_potet)/=0 ) CALL read_error(3, 'unused_potet')

      ALLOCATE ( Soil_saturated(Nhru) )
      IF ( declvar(MODNAME, 'soil_saturated', 'nhru', Nhru, 'integer', &
     &     'Flag set if infiltration saturates capillary reservoir (0=no, 1=yes)', &
     &     'none', Soil_saturated)/=0 ) CALL read_error(3, 'soil_saturated')

!      ALLOCATE ( Snowevap_aet_frac(Nhru) )
!      IF ( declvar(MODNAME, 'snowevap_aet_frac', 'nhru', Nhru, 'double', &
!     &     'Fraction of sublimation of AET for each HRU', &
!     &     'decimal fraction', Snowevap_aet_frac)/=0 ) CALL read_error(3, 'snowevap_aet_frac')

      IF ( GSFLOW_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( Nhrucell<-1 ) CALL error_stop('dimension nhrucell not specified > 0', ERROR_dim)
        ALLOCATE ( Gravity_stor_res(Nhrucell) )
        IF ( declvar(MODNAME, 'gravity_stor_res', 'nhrucell', Nhrucell, 'real', &
     &       'Storage in each gravity-flow reservoir', &
     &       'inches', Gravity_stor_res)/=0 ) CALL read_error(3, 'gravity_stor_res')

        ALLOCATE ( Sm2gw_grav(Nhrucell) )
        IF ( declvar(MODNAME, 'sm2gw_grav', 'nhrucell', Nhrucell, 'real', &
     &       'Drainage from each gravity reservoir to each MODFLOW cell', &
     &       'inches', Sm2gw_grav)/=0 ) CALL read_error(3, 'sm2gw_grav')

        IF ( declvar(MODNAME, 'basin_gvr2sm', 'one', 1, 'double', &
     &       'Basin area-weighted average gravity flow to capillary reservoirs', &
     &       'inches', Basin_gvr2sm)/=0 ) CALL read_error(3, 'basin_gvr2sm')

        ALLOCATE ( Gvr2sm(Nhru) )
        IF ( declvar(MODNAME, 'gvr2sm', 'nhru', Nhru, 'real', &
     &       'Gravity flow to soil moist replenishment for each HRU', &
     &       'inches', Gvr2sm)/=0 ) CALL read_error(3, 'gvr2sm')

        ALLOCATE ( Gw2sm_grav(Nhrucell) )
        IF ( declvar(MODNAME, 'gw2sm_grav', 'nhrucell', Nhrucell, 'real', &
     &       'Groundwater discharge to gravity-flow reservoirs', &
     &       'inches', Gw2sm_grav)/=0 ) CALL read_error(3, 'gw2sm_grav')

        ALLOCATE ( Grav_gwin(Nhru) ) ! ???
        IF ( declvar(MODNAME, 'grav_gwin', 'nhru', Nhru, 'real', &
     &       'Groundwater discharge to gravity-flow reservoirs for each HRU', &
     &       'inches', Grav_gwin)/=0 ) CALL read_error(3, 'grav_gwin')

        ALLOCATE ( Gvr_hru_pct_adjusted(Nhrucell), Replenish_frac(Nhru) )
        ALLOCATE ( Hru_gvr_count(Nhru), Hrucheck(Nhru) )
        ALLOCATE ( It0_ssres_stor(Nhru), It0_soil_rechr(Nhru), It0_soil_moist(Nhru) )
        ALLOCATE ( It0_gravity_stor_res(Nhrucell), It0_sroff(Nhru), It0_slow_stor(Nhru) )
        IF ( Nlake>0 ) ALLOCATE ( It0_potet(Nhru) )
        IF ( Call_cascade==ACTIVE ) ALLOCATE ( It0_strm_seg_in(Nsegment) )
! Declare Parameters
        ALLOCATE ( Gvr_hru_id(Nhrucell) )
        IF ( Nhru/=Nhrucell .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'gvr_hru_id', 'nhrucell', 'integer', &
     &         '0', 'bounded', 'nhru', &
     &         'Corresponding HRU id of each GVR', &
     &         'Index of the HRU associated with each gravity reservoir', &
     &         'none')/=0 ) CALL read_error(1, 'gvr_hru_id')
        ENDIF
      ENDIF

! Allocate arrays for local and variables from other modules
      ALLOCATE ( Gvr2pfr(Nhru), Swale_limit(Nhru) )
      ALLOCATE ( Pfr_dunnian_flow(Nhru), Grav_dunnian_flow(Nhru) )
      IF ( Print_debug==DEBUG_WB ) ALLOCATE( Soil_moist_ante(Nhru), Ssres_stor_ante(Nhru) )

      IF ( Print_debug==7 ) CALL PRMS_open_module_file(DBGUNT, 'soilzone.dbg')

! Declare Parameters
      IF ( Nlake>0 ) THEN
        ALLOCATE ( Lake_evap_adj(MONTHS_PER_YEAR,Nlake) )
        IF ( declparam(MODNAME, 'lake_evap_adj', 'nmonths,nlake', 'real', &
     &       '1.0', '0.5', '1.5', &
     &       'Monthly potet factor to adjust potet on lakes', &
     &       'Monthly (January to December) adjustment factor for potential ET for each lake', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'lake_evap_adj')
      ENDIF

      ALLOCATE ( Slowcoef_lin(Nhru) )
      IF ( declparam(MODNAME, 'slowcoef_lin', 'nhru', 'real', &
     &     '0.015', '0.0', '1.0', &
     &     'Linear gravity-flow reservoir routing coefficient', &
     &     'Linear coefficient in equation to route gravity-reservoir storage downslope for each HRU', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'slowcoef_lin')

      ALLOCATE ( Slowcoef_sq(Nhru) )
      IF ( declparam(MODNAME, 'slowcoef_sq', 'nhru', 'real', &
     &     '0.1', '0.0', '1.0', &
     &     'Non-linear gravity-flow reservoir routing coefficient', &
     &     'Non-linear coefficient in equation to route'// &
     &     ' gravity-reservoir storage downslope for each HRU', &
     &     'none')/=0 ) CALL read_error(1, 'slowcoef_sq')

      ALLOCATE ( Pref_flow_den(Nhru) )
      IF ( declparam(MODNAME, 'pref_flow_den', 'nhru', 'real', &
     &     '0.0', '0.0', '0.5', &
     &     'Fraction of the gravity reservoir in which preferential flow occurs for each HRU', &
     &     'Fraction of the gravity reservoir in which preferential flow occurs for each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1,'pref_flow_den')

      ALLOCATE ( Soil2gw_max(Nhru) )
      IF ( declparam(MODNAME, 'soil2gw_max', 'nhru', 'real', &
     &     '0.0', '0.0', '5.0', &
     &     'Maximum value for capillary reservoir excess to GWR', &
     &     'Maximum amount of the capillary reservoir excess that'// &
     &     ' is routed directly to the GWR for each HRU', &
     &     'inches')/=0 ) CALL read_error(1, 'soil2gw_max')

      ALLOCATE ( Soil_type(Nhru) )
      IF ( declparam(MODNAME, 'soil_type', 'nhru', 'integer', &
     &     '2', '1', '3', &
     &     'HRU soil type', 'Soil type of each HRU (1=sand; 2=loam; 3=clay)', &
     &     'none')/=0 ) CALL read_error(1, 'soil_type')

      ALLOCATE ( Fastcoef_lin(Nhru) )
      IF ( declparam(MODNAME, 'fastcoef_lin', 'nhru', 'real', &
     &     '0.1', '0.0', '1.0', &
     &     'Linear preferential-flow routing coefficient', &
     &     'Linear coefficient in equation to route preferential-flow storage downslope for each HRU', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'fastcoef_lin')

      ALLOCATE ( Fastcoef_sq(Nhru) )
      IF ( declparam(MODNAME, 'fastcoef_sq', 'nhru', 'real', &
     &     '0.8', '0.0', '1.0', &
     &     'Non-linear preferential-flow routing coefficient', &
     &     'Non-linear coefficient in equation used to route'// &
     &     ' preferential-flow storage downslope for each HRU', &
     &     'none')/=0 ) CALL read_error(1, 'fastcoef_sq')

      ALLOCATE ( Ssr2gw_rate(Nhru) )
      IF ( declparam(MODNAME, 'ssr2gw_rate', 'nssr', 'real', &
     &     '0.1', '0.0001', '999.0', &
     &     'Coefficient to route water from gravity reservoir to GWR', &
     &     'Linear coefficient in equation used to route water from'// &
     &     ' the gravity reservoir to the GWR for each HRU', &
     &     'inches/day')/=0 ) CALL read_error(1, 'ssr2gw_rate')

      ALLOCATE ( Ssr2gw_exp(Nhru) )
      IF ( declparam(MODNAME, 'ssr2gw_exp', 'nssr', 'real', &
     &     '1.0', '0.0', '3.0', &
     &     'Coefficient to route water from subsurface to groundwater', &
     &     'Non-linear coefficient in equation used to route water'// &
     &     ' from the gravity reservoir to the GWR for each HRU', &
     &     'none')/=0 ) CALL read_error(1, 'ssr2gw_exp')

      END FUNCTION szdecl

!***********************************************************************
!     szinit - Initialize soilzone module - get parameter values,
!              set initial values and check parameter values
!***********************************************************************
      INTEGER FUNCTION szinit()
      USE PRMS_SOILZONE
      USE PRMS_BASIN, ONLY: Hru_type, Hru_perv, Active_hrus, Hru_route_order, &
     &    Basin_area_inv, Hru_area, Hru_frac_perv, Numlake_hrus
      USE PRMS_FLOWVARS, ONLY: Soil_moist_max, Soil_rechr_max, &
     &    Ssres_stor, Basin_ssstor, Basin_soil_moist, Slow_stor, &
     &    Soil_moist, Sat_threshold, Soil_rechr
      USE PRMS_SNOW, ONLY: Snowcov_area
      IMPLICIT NONE
! Functions
      EXTERNAL :: init_basin_vars, checkdim_bounded_limits, error_stop
      INTEGER, EXTERNAL :: getparam
      INTRINSIC :: MIN, DBLE
! Local Variables
      INTEGER :: i, ii, ihru, icnt, j
      REAL :: hruarea, perv_area
!***********************************************************************
      szinit = 0

      IF ( getparam(MODNAME, 'slowcoef_lin', Nhru, 'real', Slowcoef_lin)/=0 ) CALL read_error(2, 'slowcoef_lin')
      IF ( getparam(MODNAME, 'slowcoef_sq', Nhru, 'real', Slowcoef_sq)/=0 ) CALL read_error(2, 'slowcoef_sq')
      IF ( getparam(MODNAME, 'pref_flow_den', Nhru, 'real', Pref_flow_den)/=0 ) CALL read_error(2, 'pref_flow_den')
      IF ( getparam(MODNAME, 'fastcoef_lin', Nhru, 'real', Fastcoef_lin)/=0 ) CALL read_error(2, 'fastcoef_lin')
      IF ( getparam(MODNAME, 'fastcoef_sq', Nhru, 'real', Fastcoef_sq)/=0 ) CALL read_error(2, 'fastcoef_sq')
      IF ( getparam(MODNAME, 'ssr2gw_rate', Nssr, 'real', Ssr2gw_rate)/=0 ) CALL read_error(2, 'ssr2gw_rate')
      IF ( getparam(MODNAME, 'ssr2gw_exp', Nssr, 'real', Ssr2gw_exp)/=0 ) CALL read_error(2, 'ssr2gw_exp')
      IF ( getparam(MODNAME, 'soil_type', Nhru, 'integer', Soil_type)/=0 ) CALL read_error(2, 'soil_type')
      IF ( getparam(MODNAME, 'soil2gw_max', Nhru, 'real', Soil2gw_max)/=0 ) CALL read_error(2, 'soil2gw_max')
      IF ( Nlake>0 ) THEN
        IF ( getparam(MODNAME, 'lake_evap_adj', MONTHS_PER_YEAR*Nlake, 'real', Lake_evap_adj)/=0 ) &
     &       CALL read_error(2, 'lake_evap_adj')
      ENDIF

      IF ( GSFLOW_flag==ACTIVE ) THEN
        IF ( Nhru/=Nhrucell ) THEN
          IF ( getparam(MODNAME, 'gvr_hru_id', Nhrucell, 'integer', Gvr_hru_id)/=0 ) CALL read_error(2, 'gvr_hru_id')
          IF ( Parameter_check_flag==1 ) &
     &         CALL checkdim_bounded_limits('gvr_hru_id', 'nhru', Gvr_hru_id, Nhrucell, 1, Nhru, Inputerror_flag)
        ELSE
          DO i = 1, Nhru
            Gvr_hru_id(i) = i
          ENDDO
        ENDIF
        Grav_gwin = 0.0 ! dimension nhru
        Gw2sm_grav = 0.0
      ENDIF

      Swale_limit = 0.0
      Pref_flag = OFF
      Pfr_dunnian_flow = 0.0
      Grav_dunnian_flow = 0.0
      Soil_lower_ratio = 0.0
      Pref_flow_thrsh = 0.0
      Pref_flow_max = 0.0
      Soil_moist_tot = 0.0
      Soil_lower = 0.0
      Soil_zone_max = 0.0
      Soil_lower_stor_max = 0.0

      ! initialize scalers
      CALL init_basin_vars()

!      Pfr_stor_frac = 0.0
!      Gvr_stor_frac = 0.0
!      Cpr_stor_frac = 0.0
!      Soil_moist_frac = 0.0
      Snow_free = 1.0 - Snowcov_area

      DO i = 1, Nhru

! WARNING changing parameter value as invalid for these HRU types, need to be 0 for computations below
        IF ( Hru_type(i)==INACTIVE .OR. Hru_type(i)==LAKE .OR. Hru_type(i)==SWALE) Pref_flow_den(i) = 0.0

        IF ( Hru_type(i)==INACTIVE .OR. Hru_type(i)==LAKE ) THEN
          Sat_threshold(i) = 0.0 ! allow modification of parameter value for lake and inactive HRUs
          Soil_rechr(i) = 0.0
          Soil_moist(i) = 0.0
          Ssres_stor(i) = 0.0
          Slow_stor(i) = 0.0
          Pref_flow_stor(i) = 0.0
!          Soil_rechr_ratio(i) = 0.0
          CYCLE
        ENDIF

        IF ( Hru_type(i)==SWALE ) THEN ! swale
          Swale_limit(i) = 3.0*Sat_threshold(i)
          Pref_flow_thrsh(i) = Sat_threshold(i)
        ELSE ! land or glacier
          Pref_flow_thrsh(i) = Sat_threshold(i)*(1.0-Pref_flow_den(i))
          Pref_flow_max(i) = Sat_threshold(i) - Pref_flow_thrsh(i)
          IF ( Pref_flow_den(i)>0.0 ) Pref_flag = ACTIVE
        ENDIF

        ! hru_type = land, swale, or glacier
        IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==5 ) THEN
          Slow_stor(i) = MIN( Ssres_stor(i), Pref_flow_thrsh(i) )
          IF ( Pref_flow_den(i)>0.0 ) THEN
            Pref_flow_stor(i) = Ssres_stor(i) - Slow_stor(i)
          ELSE
            Pref_flow_stor(i) = 0.0
          ENDIF
        ENDIF

        hruarea = Hru_area(i)
        perv_area = Hru_perv(i)
        Soil_zone_max(i) = Sat_threshold(i) + Soil_moist_max(i)*Hru_frac_perv(i)
        Soil_moist_tot(i) = Ssres_stor(i) + Soil_moist(i)*Hru_frac_perv(i)
!        Soil_moist_frac(i) = Soil_moist_tot(i)/Soil_zone_max(i)
!        Cpr_stor_frac(i) = Soil_moist(i)/Soil_moist_max(i)
!        Basin_cpr_stor_frac = Basin_cpr_stor_frac + DBLE( Cpr_stor_frac(i)*perv_area )
        Basin_cpr_stor_frac = Basin_cpr_stor_frac + DBLE( Soil_moist(i)/Soil_moist_max(i)*perv_area )
        Soil_lower(i) = Soil_moist(i) - Soil_rechr(i)
        Soil_lower_stor_max(i) = Soil_moist_max(i) - Soil_rechr_max(i)
        IF ( Soil_lower_stor_max(i)>0.0 ) Soil_lower_ratio(i) = Soil_lower(i)/Soil_lower_stor_max(i)
!        Soil_rechr_ratio(i) = Soil_rechr(i)/Soil_rechr_max(i)
!        Basin_sz_stor_frac = Basin_sz_stor_frac + DBLE( Soil_moist_frac(i)*hruarea )
        Basin_sz_stor_frac = Basin_sz_stor_frac + DBLE( Soil_moist_tot(i)/Soil_zone_max(i)*hruarea )
        Basin_soil_lower_stor_frac = Basin_soil_lower_stor_frac + DBLE( Soil_lower_ratio(i)*perv_area )
!        Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac + DBLE( Soil_rechr_ratio(i)*perv_area )
        Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac + DBLE( Soil_rechr(i)/Soil_rechr_max(i)*perv_area )
        Basin_soil_moist = Basin_soil_moist + DBLE( Soil_moist(i)*Hru_perv(i) )
        Basin_soil_moist_tot = Basin_soil_moist_tot + DBLE( Soil_moist_tot(i)*hruarea )
        ! rsr, 6/12/2014 potential problem for GSFLOW if sum of slow_stor /= gravity_stor_res
        Basin_slstor = Basin_slstor + DBLE( Slow_stor(i)*hruarea )
        Basin_ssstor = Basin_ssstor + DBLE( Ssres_stor(i)*hruarea )
        Basin_soil_rechr = Basin_soil_rechr + DBLE( Soil_rechr(i)*perv_area )
        IF ( Pref_flow_den(i)>0.0 ) THEN
          Basin_pref_stor = Basin_pref_stor + DBLE( Pref_flow_stor(i)*hruarea )
!          Pfr_stor_frac(i) = Pref_flow_stor(i)/Pref_flow_max(i)
!          Basin_pfr_stor_frac = Basin_pfr_stor_frac + DBLE( Pfr_stor_frac(i)*hruarea )
          Basin_pfr_stor_frac = Basin_pfr_stor_frac + DBLE( Pref_flow_stor(i)/Pref_flow_max(i)*hruarea )
          IF ( Pref_flow_thrsh(i)>0.0 ) THEN
!            Gvr_stor_frac(i) = Slow_stor(i)/Pref_flow_thrsh(i)
            Basin_gvr_stor_frac = Basin_gvr_stor_frac + DBLE( Slow_stor(i)/Pref_flow_thrsh(i)*hruarea )
          ENDIF
        ENDIF
      ENDDO
      Basin_soil_rechr = Basin_soil_rechr*Basin_area_inv
      Basin_ssstor = Basin_ssstor*Basin_area_inv
      Basin_slstor = Basin_slstor*Basin_area_inv
      Basin_soil_moist = Basin_soil_moist*Basin_area_inv
      Basin_soil_moist_tot = Basin_soil_moist_tot*Basin_area_inv
      Basin_pref_stor = Basin_pref_stor*Basin_area_inv
      Last_soil_moist = Basin_soil_moist
      Last_ssstor = Basin_ssstor
      Basin_cpr_stor_frac = Basin_cpr_stor_frac*Basin_area_inv
      Basin_gvr_stor_frac = Basin_gvr_stor_frac*Basin_area_inv
      Basin_pfr_stor_frac = Basin_pfr_stor_frac*Basin_area_inv
      Basin_sz_stor_frac = Basin_sz_stor_frac*Basin_area_inv
      Basin_soil_lower_stor_frac = Basin_soil_lower_stor_frac*Basin_area_inv
      Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac*Basin_area_inv

! initialize arrays (dimensioned Nhru)
      Dunnian_flow = 0.0
      IF ( Cascade_flag>CASCADE_OFF ) THEN
        Upslope_interflow = 0.0D0
        Upslope_dunnianflow = 0.0D0
        Hru_sz_cascadeflow = 0.0
!        Cap_upflow_max = 0.0
!        Cascade_interflow = 0.0
!        Cascade_dunnianflow = 0.0
        IF ( Numlake_hrus>0 ) Lakein_sz = 0.0D0
      ENDIF
      Cap_infil_tot = 0.0
      Pref_flow_infil = 0.0
      Pref_flow_in = 0.0
      Pref_flow = 0.0
      Gvr2pfr = 0.0
      Swale_actet = 0.0
      Perv_actet = 0.0
!      Perv_avail_et = 0.0
      Recharge = 0.0
      Cap_waterin = 0.0
      Potet_lower = 0.0
      Potet_rechr = 0.0
      Unused_potet = 0.0 ! dimension nhru
      Soil_saturated = OFF
      IF ( Pref_flag==ACTIVE .AND. GSFLOW_flag==ACTIVE ) ALLOCATE ( It0_pref_flow_stor(Nhru) )
!      Interflow_max = 0.0
!      Snowevap_aet_frac = 0.0

! initialize GSFLOW arrays
      IF ( GSFLOW_flag==ACTIVE ) THEN
        Gvr2sm = 0.0 ! dimension nhru
        Sm2gw_grav = 0.0 ! dimension nhrucell

        Max_gvrs = 1
        Hrucheck = 1
        Hru_gvr_count = 0
        Replenish_frac = 0.0
        DO i = 1, Nhrucell
          ihru = Gvr_hru_id(i)
          IF ( Hru_type(ihru)==INACTIVE .OR. Hru_type(ihru)==LAKE ) THEN
            Gravity_stor_res(i) = 0.0
            Hrucheck(ihru) = 0
          ELSE
            ! set only for cold start simulations
            IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==5 ) &
      &          Gravity_stor_res(i) = Ssres_stor(ihru)
            Hru_gvr_count(ihru) = Hru_gvr_count(ihru) + 1
            IF ( Hru_gvr_count(ihru)>Max_gvrs ) Max_gvrs = Hru_gvr_count(ihru)
            Replenish_frac(ihru) = Soil_rechr_max(ihru)/Soil_moist_max(ihru)
          ENDIF
        ENDDO
        ALLOCATE ( Hru_gvr_index(Max_gvrs, Nhru) )
        Hru_gvr_index = 0
        IF ( Nhru==Nhrucell ) THEN
          IF ( Max_gvrs/=1 ) &
     &         CALL error_stop('nhru=nhrucell, but, gvr_hru_id array specifies more than one GVR for an HRU', ERROR_dim)
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            IF ( Hru_type(i)==LAKE ) CYCLE
            Hru_gvr_index(1, i) = i
          ENDDO
        ELSE
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            IF ( Hru_type(i)==LAKE ) CYCLE
            icnt = 0
            DO ii = 1, Nhrucell
              IF ( Gvr_hru_id(ii)==i ) THEN
                icnt = icnt + 1
                Hru_gvr_index(icnt, i) = ii
                IF ( icnt==Hru_gvr_count(i) ) EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      END FUNCTION szinit

!***********************************************************************
!     szrun - Does soil water balance for each HRU, adds in infiltration
!             then computes actual et and apportions remainder between
!             recharge of soil moisture, soil storage available for
!             interflow, excess routed to stream,
!             and groundwater reservoirs
!***********************************************************************
      INTEGER FUNCTION szrun()
      USE PRMS_SOILZONE
      USE PRMS_BASIN, ONLY: Hru_type, Hru_perv, Hru_frac_perv, &
     &    Hru_route_order, Active_hrus, Basin_area_inv, Hru_area, &
     &    Lake_hru_id, Cov_type, Numlake_hrus, Hru_area_dble
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt, Transp_on, Potet, Basin_potet
! WARNING!!! Sroff, Basin_sroff, and Strm_seg_in can be updated
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_actet, Hru_actet, &
     &    Ssres_flow, Soil_to_gw, Basin_soil_to_gw, Ssr_to_gw, &
     &    Soil_to_ssr, Basin_lakeevap, Basin_perv_et, Basin_swale_et, &
     &    Sroff, Soil_moist_max, Infil, Soil_rechr_max, Ssres_in, &
     &    Basin_soil_moist, Basin_ssstor, Slow_stor, Slow_flow, &
     &    Ssres_stor, Soil_moist, Sat_threshold, Soil_rechr, Basin_sroff, Basin_lake_stor
      USE PRMS_WATER_USE, ONLY: Soilzone_gain
      USE PRMS_CASCADE, ONLY: Ncascade_hru
      USE PRMS_SET_TIME, ONLY: Nowmonth, Cfs_conv !, Nowday
      USE PRMS_INTCP, ONLY: Hru_intcpevap
      USE PRMS_SNOW, ONLY: Snowcov_area, Snow_evap
      USE PRMS_SRUNOFF, ONLY: Hru_impervevap, Strm_seg_in, Dprst_evap_hru, Dprst_seep_hru, Frozen
      IMPLICIT NONE
! Functions
      INTRINSIC :: MIN, ABS, MAX, SNGL, DBLE
      EXTERNAL :: compute_soilmoist, compute_szactet, compute_cascades, compute_gravflow
      EXTERNAL :: compute_interflow, compute_gwflow, init_basin_vars, print_date
! Local Variables
      INTEGER :: i, k, update_potet, compute_lateral
      REAL :: dunnianflw, interflow, perv_area, harea
      REAL :: dnslowflow, dnpreflow, dndunn, availh2o, avail_potet
      REAL :: gvr_maxin, topfr !, tmp
      REAL :: dunnianflw_pfr, dunnianflw_gvr, pref_flow_maxin
      REAL :: perv_frac, capacity, capwater_maxin, ssresin
      REAL :: cap_upflow_max, unsatisfied_et, pervactet, prefflow
      DOUBLE PRECISION :: gwin
      INTEGER :: cfgi_frozen_hru
!***********************************************************************
      szrun = 0

! It0 variables used with MODFLOW integration to save iteration states.
      IF ( GSFLOW_flag==ACTIVE ) THEN
        IF ( Kkiter>1 ) THEN
          Soil_rechr = It0_soil_rechr
          Soil_moist = It0_soil_moist
          Ssres_stor = It0_ssres_stor
          Slow_stor = It0_slow_stor
          Sroff = It0_sroff
          Gravity_stor_res = It0_gravity_stor_res
          IF ( Pref_flag==ACTIVE ) Pref_flow_stor = It0_pref_flow_stor
          IF ( Nlake>0 ) Potet = It0_potet ! possible that potet is increased if lake_evap_adj > 1.0
          Basin_soil_moist = It0_basin_soil_moist
          Basin_ssstor = It0_basin_ssstor
          IF ( Call_cascade==ACTIVE ) Strm_seg_in = It0_strm_seg_in
        ELSE
          It0_soil_rechr = Soil_rechr
          It0_soil_moist = Soil_moist
          It0_ssres_stor = Ssres_stor
          It0_slow_stor = Slow_stor
          It0_sroff = Sroff
          It0_gravity_stor_res = Gravity_stor_res
          IF ( Pref_flag==ACTIVE ) It0_pref_flow_stor = Pref_flow_stor
          IF ( Nlake>0 ) It0_potet = Potet
          It0_basin_soil_moist = Basin_soil_moist
          It0_basin_ssstor = Basin_ssstor
          IF ( Call_cascade==ACTIVE ) It0_strm_seg_in = Strm_seg_in
          Gw2sm_grav = 0.0
        ENDIF
        Sm2gw_grav = 0.0
      ENDIF

      IF ( Cascade_flag>CASCADE_OFF ) THEN
        Upslope_interflow = 0.0D0
        Upslope_dunnianflow = 0.0D0
        IF ( Numlake_hrus>0 ) THEN
          Lakein_sz = 0.0D0
          Basin_lakeinsz = 0.0D0
        ENDIF
      ENDIF

      IF ( Print_debug==DEBUG_WB ) THEN
        Soil_moist_ante = Soil_moist
        Ssres_stor_ante = Ssres_stor
        Last_soil_moist = Basin_soil_moist
        Last_ssstor = Basin_ssstor
      ENDIF
      CALL init_basin_vars()
      gwin = 0.0D0
      ! Soil_to_gw and Soil_to_ssr for whole HRU
      Soil_to_gw = 0.0
      Soil_to_ssr = 0.0
!      Snowevap_aet_frac = 0.0
      ! gravity reservoir variables for whole HRU
      Ssr_to_gw = 0.0
      Slow_flow = 0.0
      Ssres_flow = 0.0
      Cap_waterin = 0.0
      Soil_saturated = OFF
      update_potet = OFF
      DO k = 1, Active_hrus
        i = Hru_route_order(k)

        Hru_actet(i) = Hru_impervevap(i) + Hru_intcpevap(i) + Snow_evap(i)
        IF ( Dprst_flag==ACTIVE ) Hru_actet(i) = Hru_actet(i) + Dprst_evap_hru(i)
        harea = Hru_area(i)

        IF ( Hru_type(i)==LAKE ) THEN ! lake or reservoir
          !WARNING, RSR, if hru_actet>water in lake, then budget error
          Hru_actet(i) = (Potet(i) - Hru_actet(i))*Lake_evap_adj(Nowmonth,Lake_hru_id(i))
          IF ( Hru_actet(i)>Potet(i) ) THEN
            PRINT *, 'WARNING, lake evap > potet, for HRU:', i, ' potential ET increased to adjusted lake ET'
            PRINT *, Hru_actet(i), Potet(i), Hru_actet(i) - Potet(i)
            Potet(i) = Hru_actet(i) ! this could be a problem when it happens
            update_potet = ACTIVE
          ENDIF
          Unused_potet(i) = Potet(i) - Hru_actet(i)
          Basin_actet = Basin_actet + DBLE( Hru_actet(i)*harea )
          Basin_lakeevap = Basin_lakeevap + DBLE( Hru_actet(i)*harea )
          Basin_lakeprecip = Basin_lakeprecip + DBLE( Hru_ppt(i)*harea )
          IF ( Cascade_flag>CASCADE_OFF ) THEN
            ! if lake HRU doesn't cascade, should we limit ET to
            !  water entering the HRU to this point (no gwflow yet)
            Lakein_sz(i) = Upslope_interflow(i) + Upslope_dunnianflow(i)
            Basin_lakeinsz = Basin_lakeinsz + Lakein_sz(i)*Hru_area_dble(i)
          ENDIF
          CYCLE
        ENDIF

        !Hru_type can be 1 (land) or 3 (swale) or 4 (glacier)
        compute_lateral = OFF ! swale
        IF ( Hru_type(i)==LAND .OR. Hru_type(i)==GLACIER ) compute_lateral = ACTIVE
        ! perv_frac has to be > 0.001
        perv_area = Hru_perv(i)
        perv_frac = Hru_frac_perv(i)

        avail_potet = Potet(i) - Hru_actet(i)
        IF ( avail_potet<0.0 ) avail_potet = 0.0
!        Snowevap_aet_frac(i) = 0.0

        !Hru_type can be 1 (land) or 3 (swale) or 4 (glacier)

!******Add infiltration to soil and compute excess
        dunnianflw = 0.0
        dunnianflw_pfr = 0.0
        dunnianflw_gvr = 0.0
        interflow = 0.0

!******Add infiltration to soil and compute excess
        !infil_tot is the depth in whole HRU
        !capillary reservoir for pervious area
        !preferential flow reservoir for whole HRU
        !gravity reservoir for whole HRU
        !upslope flow for whole HRU

!******if cascading flow available from upslope cascades
!****** add soil excess (Dunnian flow) to infiltration
        ! perv_frac has to be > 0.001
        ! infil for pervious portion of HRU
        capwater_maxin = Infil(i)

        IF ( Soilzone_add_water_use==ACTIVE ) THEN
          IF ( Soilzone_gain(i)>0.0 ) capwater_maxin = capwater_maxin + Soilzone_gain(i)/perv_area/SNGL(Cfs_conv) ! ??? is this harea
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
        prefflow = 0.0
        IF ( Pref_flow_den(i)>0.0 ) THEN
          pref_flow_maxin = 0.0
          Pref_flow_infil(i) = 0.0
          IF ( capwater_maxin>0.0 ) THEN
            ! pref_flow for whole HRU
            pref_flow_maxin = capwater_maxin*Pref_flow_den(i)
            capwater_maxin = capwater_maxin - pref_flow_maxin
            pref_flow_maxin = pref_flow_maxin*perv_frac
            IF ( cfgi_frozen_hru==ACTIVE ) THEN
              dunnianflw_pfr = pref_flow_maxin
            ELSE
              ! compute contribution to preferential-flow reservoir storage
              Pref_flow_stor(i) = Pref_flow_stor(i) + pref_flow_maxin
              dunnianflw_pfr = MAX( 0.0, Pref_flow_stor(i)-Pref_flow_max(i) )
            ENDIF
            IF ( dunnianflw_pfr>0.0 ) THEN
              Basin_dunnian_pfr = Basin_dunnian_pfr + dunnianflw_pfr*harea
              Pref_flow_stor(i) = Pref_flow_max(i)
            ENDIF
            Pref_flow_infil(i) = pref_flow_maxin - dunnianflw_pfr
            Basin_pref_flow_infil = Basin_pref_flow_infil + Pref_flow_infil(i)*harea
          ENDIF
          Pfr_dunnian_flow(i) = dunnianflw_pfr
        ENDIF

        IF ( Cascade_flag>CASCADE_OFF ) THEN
!          Cap_upflow_max(i) = SNGL(Upslope_dunnianflow(i)+Upslope_interflow(i))/perv_frac
!          capwater_maxin = capwater_maxin + Cap_upflow_max(i)
!          Basin_cap_up_max = Basin_cap_up_max + Cap_upflow_max(i)*perv_area
          cap_upflow_max = SNGL(Upslope_dunnianflow(i)+Upslope_interflow(i))/perv_frac
          capwater_maxin = capwater_maxin + cap_upflow_max
          Basin_cap_up_max = Basin_cap_up_max + cap_upflow_max*perv_area
        ENDIF
        Cap_infil_tot(i) = capwater_maxin*perv_frac
        Basin_cap_infil_tot = Basin_cap_infil_tot + DBLE( Cap_infil_tot(i)*harea )

!******Add infiltration to soil and compute excess
        gvr_maxin = 0.0
        Cap_waterin(i) = capwater_maxin

        IF ( cfgi_frozen_hru==OFF ) THEN
          ! call even if capwater_maxin = 0, just in case soil_moist now > Soil_moist_max
          IF ( capwater_maxin+Soil_moist(i)>0.0 ) THEN
            CALL compute_soilmoist(Cap_waterin(i), Soil_moist_max(i), &
     &           Soil_rechr_max(i), Soil2gw_max(i), gvr_maxin, &
     &           Soil_moist(i), Soil_rechr(i), Soil_to_gw(i), perv_frac)
            Cap_waterin(i) = Cap_waterin(i)*perv_frac
            Basin_capwaterin = Basin_capwaterin + DBLE( Cap_waterin(i)*harea )
            Basin_soil_to_gw = Basin_soil_to_gw + DBLE( Soil_to_gw(i)*harea )
            Basin_sm2gvr_max = Basin_sm2gvr_max + DBLE( gvr_maxin*harea )
            Soil_to_ssr(i) = gvr_maxin
          ENDIF
        ENDIF

! compute slow interflow and ssr_to_gw
        topfr = 0.0
        IF ( GSFLOW_flag==ACTIVE ) THEN
          ! capacity for whole HRU
          capacity = (Soil_moist_max(i) - Soil_moist(i))*perv_frac
          CALL compute_gravflow(i, capacity, Slowcoef_lin(i), &
     &                          Slowcoef_sq(i), Ssr2gw_rate(i), Ssr2gw_exp(i), &
     &                          gvr_maxin, Pref_flow_thrsh(i), topfr, &
     &                          Ssr_to_gw(i), Slow_flow(i), Slow_stor(i), &
     &                          Gvr2sm(i), Soil_to_gw(i), gwin, compute_lateral)
          ! adjust soil moisture with replenish amount
          IF ( Gvr2sm(i)>0.0 ) THEN
            Soil_moist(i) = Soil_moist(i) + Gvr2sm(i)/perv_frac ! ??? could this be bigger than soil_moist_max ??? (add to Dunnian)
!            IF ( Soil_moist(i)>Soil_moist_max(i) ) PRINT *, 'CAP sm>max', Soil_moist(i), Soil_moist_max(i), i
            IF ( Soilzone_aet_flag==ACTIVE ) THEN
              Soil_lower(i) = MIN( Soil_lower_stor_max(i), Soil_moist(i) - Soil_rechr(i) + Gvr2sm(i)/perv_frac )
              Soil_rechr(i) = Soil_moist(i) - Soil_lower(i)
!              excess = MAX( 0.0, Soil_lower(i) - Soil_lower_stor_max(i) )
!              if ( abs(soil_lower(i) + soil_rechr(i) - soil_moist(i))>NEARZERO ) THEN
!                print *, 'excess', excess, Soil_lower_stor_max(i), soil_lower(i) + soil_rechr(i)-soil_moist(i)
!                print *, soil_lower(i), soil_rechr(i), soil_moist(i)
!              endif
            ELSE
              Soil_rechr(i) = Soil_rechr(i) + Gvr2sm(i)/perv_frac*Replenish_frac(i)
              Soil_rechr(i) = MIN( Soil_rechr_max(i), Soil_rechr(i) )
            ENDIF
            Basin_gvr2sm = Basin_gvr2sm + DBLE( Gvr2sm(i)*harea )
!          ELSEIF ( Gvr2sm(i)<-NEARZERO ) THEN
!            PRINT *, 'negative gvr2sm, HRU:', i, Gvr2sm(i)
!            Gvr2sm(i) = 0.0
          ENDIF
          Grav_gwin(i) = SNGL( gwin )
          Basin_sz_gwin = Basin_sz_gwin + gwin*DBLE( harea )
        ELSE
          availh2o = Slow_stor(i) + gvr_maxin
          IF ( compute_lateral==ACTIVE ) THEN
            topfr = MAX( 0.0, availh2o-Pref_flow_thrsh(i) )
            ssresin = gvr_maxin - topfr
            Slow_stor(i) = availh2o - topfr
            ! compute slow contribution to interflow, if any
            IF ( Slow_stor(i)>0.0 ) &
     &           CALL compute_interflow(Slowcoef_lin(i), Slowcoef_sq(i), &
     &                                  ssresin, Slow_stor(i), Slow_flow(i))
          ELSEIF ( Hru_type(i)==SWALE ) THEN
            Slow_stor(i) = availh2o
          ENDIF
          IF ( Slow_stor(i)>0.0 .AND. Ssr2gw_rate(i)>0.0 ) &
       &       CALL compute_gwflow(Ssr2gw_rate(i), Ssr2gw_exp(i), Ssr_to_gw(i), Slow_stor(i))
        ENDIF

        ! compute contribution to Dunnian flow from PFR, if any
        IF ( Pref_flow_den(i)>0.0 ) THEN
          availh2o = Pref_flow_stor(i) + topfr
          dunnianflw_gvr = MAX( 0.0, availh2o-Pref_flow_max(i) )
          IF ( dunnianflw_gvr>0.0 ) THEN
            topfr = topfr - dunnianflw_gvr
            IF ( topfr<0.0 ) THEN
!              IF ( topfr<-NEARZERO .AND. Print_debug>DEBUG_less ) PRINT *, 'gvr2pfr<0', topfr, dunnianflw_gvr, &
!     &             Pref_flow_max(i), Pref_flow_stor(i), gvr_maxin
              topfr = 0.0
            ENDIF
          ENDIF
          Pref_flow_in(i) = Pref_flow_infil(i) + topfr
          Pref_flow_stor(i) = Pref_flow_stor(i) + topfr
          IF ( Pref_flow_stor(i)>0.0 ) &
     &         CALL compute_interflow(Fastcoef_lin(i), Fastcoef_sq(i), &
     &                                Pref_flow_in(i), Pref_flow_stor(i), prefflow)
          Basin_pref_stor = Basin_pref_stor + DBLE( Pref_flow_stor(i)*harea )
!          Pfr_stor_frac(i) = Pref_flow_stor(i)/Pref_flow_max(i)
!          Basin_pfr_stor_frac = Basin_pfr_stor_frac + Pfr_stor_frac(i)*harea
          Basin_pfr_stor_frac = Basin_pfr_stor_frac + Pref_flow_stor(i)/Pref_flow_max(i)*harea
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
        IF ( Soil_moist(i)>0.0 .AND. cfgi_frozen_hru==OFF ) THEN
          CALL compute_szactet(Soil_moist_max(i), Soil_rechr_max(i), Transp_on(i), Cov_type(i), &
     &                         Soil_type(i), Soil_moist(i), Soil_rechr(i), pervactet, avail_potet, &
     &                         Snow_free(i), Potet_rechr(i), Potet_lower(i), &
     &                         Potet(i), perv_frac, Soil_saturated(i))
          ! sanity check
!          IF ( pervactet>avail_potet ) THEN
!            Soil_moist(i) = Soil_moist(i) + pervactet - avail_potet
!            pervactet = avail_potet
!            PRINT *, 'perv_et problem', pervactet, Avail_potet
!          ENDIF
        ENDIF
!        Perv_avail_et(i) = avail_potet

        ! sanity check
!        IF ( Soil_moist(i)<0.0 ) THEN
!          IF ( Print_debug>-1 ) PRINT *, i, Soil_moist(i), ' negative'
!          IF ( pervactet>=ABS(Soil_moist(i)) ) THEN
!            pervactet = pervactet + Soil_moist(i)
!            Soil_moist(i) = 0.0
!          ENDIF
!          IF ( Soil_moist(i)<-NEARZERO ) THEN
!            IF ( Print_debug>-1 ) PRINT *, 'HRU:', i, ' soil_moist<0.0', Soil_moist(i)
!          ENDIF
!          Soil_moist(i) = 0.0
!        ENDIF

        Hru_actet(i) = Hru_actet(i) + pervactet*perv_frac
        avail_potet = Potet(i) - Hru_actet(i)
        ! sanity check
!        IF ( avail_potet<0.0 ) THEN
!          IF ( Print_debug>-1 ) THEN
!            IF ( avail_potet<-NEARZERO ) PRINT *, 'hru_actet>potet', i, &
!     &           Nowmonth, Nowday, Hru_actet(i), Potet(i), avail_potet
!          ENDIF
!          Hru_actet(i) = Potet(i)
!          tmp = avail_potet/perv_frac
!          pervactet = pervactet + tmp
!          Soil_moist(i) = Soil_moist(i) - tmp
!          Soil_rechr(i) = Soil_rechr(i) - tmp
!          IF ( Soil_rechr(i)<0.0 ) Soil_rechr(i) = 0.0
!          IF ( Soil_moist(i)<0.0 ) Soil_moist(i) = 0.0
!        ENDIF
        Perv_actet(i) = pervactet

! soil_moist & soil_rechr multiplied by perv_area instead of harea
        Soil_lower(i) = Soil_moist(i) - Soil_rechr(i)
        Basin_soil_moist = Basin_soil_moist + DBLE( Soil_moist(i)*perv_area )
        Basin_soil_rechr = Basin_soil_rechr + DBLE( Soil_rechr(i)*perv_area )
        Basin_perv_et = Basin_perv_et + DBLE( Perv_actet(i)*perv_area )

! if HRU cascades,
! compute interflow and excess flow to each HRU or stream
        IF ( compute_lateral==ACTIVE ) THEN
          interflow = Slow_flow(i) + prefflow
!          Interflow_max(i) = interflow
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
!              Cascade_interflow(i) = dnslowflow + dnpreflow
!              Cascade_dunnianflow(i) = dndunn
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
!        Soil_rechr_ratio(i) = Soil_rechr(i)/Soil_rechr_max(i)
        Ssres_in(i) = Soil_to_ssr(i) + Pref_flow_infil(i) + SNGL( gwin )
        Basin_ssin = Basin_ssin + DBLE( Ssres_in(i)*harea )
        Basin_ssstor = Basin_ssstor + DBLE( Ssres_stor(i)*harea )
        Basin_slstor = Basin_slstor + DBLE( Slow_stor(i)*harea )
        Soil_moist_tot(i) = Ssres_stor(i) + Soil_moist(i)*perv_frac
        Basin_soil_moist_tot = Basin_soil_moist_tot + DBLE( Soil_moist_tot(i)*harea )
!        Soil_moist_frac(i) = Soil_moist_tot(i)/Soil_zone_max(i)
!        Cpr_stor_frac(i) = Soil_moist(i)/Soil_moist_max(i)
!        IF ( Pref_flow_thrsh(i)>0.0 ) Gvr_stor_frac(i) = Slow_stor(i)/Pref_flow_thrsh(i)
!        Basin_cpr_stor_frac = Basin_cpr_stor_frac + Cpr_stor_frac(i)*perv_area
!        Basin_gvr_stor_frac = Basin_gvr_stor_frac + Gvr_stor_frac(i)*harea
!        Basin_sz_stor_frac = Basin_sz_stor_frac + Soil_moist_frac(i)*harea
        Basin_cpr_stor_frac = Basin_cpr_stor_frac + Soil_moist(i)/Soil_moist_max(i)*perv_area
        IF ( Pref_flow_thrsh(i)>0.0 ) Basin_gvr_stor_frac = Basin_gvr_stor_frac + Slow_stor(i)/Pref_flow_thrsh(i)*harea
        Basin_sz_stor_frac = Basin_sz_stor_frac + Soil_moist_tot(i)/Soil_zone_max(i)*harea
        Basin_soil_lower_stor_frac = Basin_soil_lower_stor_frac + Soil_lower_ratio(i)*perv_area
!        Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac + Soil_rechr_ratio(i)*perv_area
        Basin_soil_rechr_stor_frac = Basin_soil_rechr_stor_frac + Soil_rechr(i)/Soil_rechr_max(i)*perv_area
        Recharge(i) = Soil_to_gw(i) + Ssr_to_gw(i)
        IF ( Dprst_flag==1 ) Recharge(i) = Recharge(i) + SNGL( Dprst_seep_hru(i) )
        Basin_recharge = Basin_recharge + DBLE( Recharge(i)*harea )
        Grav_dunnian_flow(i) = dunnianflw_gvr
        Unused_potet(i) = Potet(i) - Hru_actet(i)
        Basin_actet = Basin_actet + DBLE( Hru_actet(i)*harea )
!        IF ( Hru_actet(i)>0.0 ) Snowevap_aet_frac(i) = Snow_evap(i)/Hru_actet(i)

      ENDDO
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
      Basin_gvr2sm = Basin_gvr2sm*Basin_area_inv
      Basin_sz_gwin = Basin_sz_gwin*Basin_area_inv
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

      END FUNCTION szrun

!***********************************************************************
!     Add infiltration to soil and compute excess
!     Soil_to_gw and Soil_to_ssr for whole HRU
!***********************************************************************
      SUBROUTINE compute_soilmoist(Infil, Soil_moist_max, &
     &           Soil_rechr_max, Soil2gw_max, Soil_to_ssr, Soil_moist, &
     &           Soil_rechr, Soil_to_gw, Perv_frac)
      IMPLICIT NONE
! Function
      INTRINSIC :: MIN
! Arguments
      REAL, INTENT(IN) :: Perv_frac, Soil_moist_max, Soil_rechr_max, Soil2gw_max
      REAL, INTENT(INOUT) :: Infil, Soil_moist, Soil_rechr, Soil_to_gw, Soil_to_ssr
! Local Variables
      REAL :: excs
!***********************************************************************
      Soil_rechr = MIN( (Soil_rechr+Infil), Soil_rechr_max )
      ! soil_moist_max from previous time step or soil_moist_max has
      ! changed for a restart simulation
      excs = Soil_moist + Infil
      Soil_moist = MIN( excs, Soil_moist_max )
      excs = (excs - Soil_moist_max)*Perv_frac
      IF ( excs>0.0 ) THEN
        IF ( Soil2gw_max>0.0 ) THEN
          Soil_to_gw = MIN( Soil2gw_max, excs )
          excs = excs - Soil_to_gw
        ENDIF
        IF ( excs>Infil*Perv_frac ) THEN !probably dynamic
          Infil = 0.0
        ELSE
          Infil = Infil - excs/Perv_frac         !???? what if Infil<0 ??? might happen with dynamic and small values, maybe ABS < NEARZERO = 0.0
!          IF ( Infil<0.0 ) THEN
!            IF ( Infil<-0.0001 ) THEN
!              PRINT *, 'negative infil', infil, soil_moist, excs
!              Soil_moist = Soil_moist + Infil
!            ENDIF
!            Infil = 0.0
!          ENDIF
        ENDIF

        Soil_to_ssr = excs
        IF ( Soil_to_ssr<0.0 ) Soil_to_ssr = 0.0
      ENDIF

      END SUBROUTINE compute_soilmoist

!***********************************************************************
!     Compute actual evapotranspiration
!***********************************************************************
      SUBROUTINE compute_szactet(Soil_moist_max, Soil_rechr_max, &
     &           Transp_on, Cov_type, Soil_type, &
     &           Soil_moist, Soil_rechr, Perv_actet, Avail_potet, &
     &           Snow_free, Potet_rechr, Potet_lower, Potet, Perv_frac, Soil_saturated)
      USE PRMS_CONSTANTS, ONLY: NEARZERO, BARESOIL, SAND, LOAM, CLAY, ACTIVE, OFF
      USE PRMS_MODULE, ONLY: Soilzone_aet_flag
      USE PRMS_SOILZONE, ONLY: Et_type
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Transp_on, Cov_type, Soil_type
      INTEGER, INTENT(INOUT) :: Soil_saturated
      REAL, INTENT(IN) :: Soil_moist_max, Soil_rechr_max, Snow_free, Potet, Perv_frac
      REAL, INTENT(INOUT) :: Soil_moist, Soil_rechr, Avail_potet, Potet_rechr, Potet_lower
      REAL, INTENT(OUT) :: Perv_actet
! Local Variables
      REAL, PARAMETER :: ONETHIRD = 1.0/3.0, TWOTHIRDS = 2.0/3.0
      REAL :: et, pcts, pctr, pet
!***********************************************************************
!******Determine if evaporation(Et_type = 2) or transpiration plus
!******evaporation(Et_type = 3) are active.  if not, Et_type = 1

      IF ( Soilzone_aet_flag==ACTIVE ) THEN
        pet = Potet
      ELSE
        pet = Avail_potet
      ENDIF
      IF ( Avail_potet<NEARZERO ) THEN
        Et_type = 1
        pet = 0.0
      ELSEIF ( Transp_on==OFF ) THEN
        IF ( Snow_free<0.01 ) THEN
          Et_type = 1
        ELSE
          Et_type = 2
        ENDIF
      ELSEIF ( Cov_type>BARESOIL ) THEN
        Et_type = 3
      ELSEIF ( Snow_free<0.01 ) THEN
        Et_type = 1
      ELSE
        Et_type = 2
      ENDIF

      IF ( Et_type>1 ) THEN
        pcts = Soil_moist/Soil_moist_max
        IF ( pcts>0.9999 ) Soil_saturated = 1
        pctr = Soil_rechr/Soil_rechr_max
        Potet_lower = pet
        Potet_rechr = pet

!******sandy soil
        IF ( Soil_type==SAND ) THEN
          IF ( pcts<0.25 ) Potet_lower = 0.5*pcts*pet
          IF ( pctr<0.25 ) Potet_rechr = 0.5*pctr*pet
!******loam soil
        ELSEIF ( Soil_type==LOAM ) THEN
          IF ( pcts<0.5 ) Potet_lower = pcts*pet
          IF ( pctr<0.5 ) Potet_rechr = pctr*pet
!******clay soil
        ELSEIF ( Soil_type==CLAY ) THEN
          IF ( pcts<TWOTHIRDS .AND. pcts>ONETHIRD ) THEN
            Potet_lower = pcts*pet
          ELSEIF ( pcts<=ONETHIRD ) THEN
            Potet_lower = 0.5*pcts*pet
          ENDIF
          IF ( pctr<TWOTHIRDS .AND. pctr>ONETHIRD ) THEN
            Potet_rechr = pctr*pet
          ELSEIF ( pctr<=ONETHIRD ) THEN
            Potet_rechr = 0.5*pctr*pet
          ENDIF
        ENDIF

!******Soil moisture accounting
        IF ( Et_type==2 ) Potet_rechr = Potet_rechr*Snow_free

!!!! need to limit Potet_rechr and Potet_lower (pervious area) by Avail_potet (whole HRU)
        IF ( Potet_rechr*Perv_frac > Avail_potet ) THEN
!          print *, 'reducing potet_rechr in cap reservoir', perv_frac, pet
!          print *, potet_rechr, potet_lower, avail_potet, (Potet_rechr+Potet_lower)*Perv_frac
          Potet_rechr = Avail_potet / Perv_frac
        ENDIF
        IF ( Potet_lower*Perv_frac > Avail_potet ) THEN
!          print *, 'reducing potet_lower in cap reservoir', perv_frac, pet
!          print *, potet_rechr, potet_lower, avail_potet, (Potet_rechr+Potet_lower)*Perv_frac
          Potet_lower = Avail_potet / Perv_frac
        ENDIF

        IF ( Potet_rechr>Soil_rechr ) THEN
          Potet_rechr = Soil_rechr
          Soil_rechr = 0.0
        ELSE
          Soil_rechr = Soil_rechr - Potet_rechr
        ENDIF
        IF ( Et_type==2 .OR. Potet_rechr>=Potet_lower ) THEN
          IF ( Potet_rechr>Soil_moist ) THEN
            Potet_rechr = Soil_moist
            Soil_moist = 0.0
          ELSE
            Soil_moist = Soil_moist - Potet_rechr
          ENDIF
          et = Potet_rechr
        ELSEIF ( Potet_lower>Soil_moist ) THEN
          et = Soil_moist
          Soil_moist = 0.0
        ELSE
          Soil_moist = Soil_moist - Potet_lower
          et = Potet_lower
        ENDIF
        IF ( Soil_rechr>Soil_moist ) Soil_rechr = Soil_moist
      ELSE
        et = 0.0
      ENDIF
      Perv_actet = et
      ! sanity check
      IF ( Perv_actet*Perv_frac-Avail_potet > NEARZERO ) THEN
        PRINT *, 'perv_et problem', Perv_actet*Perv_frac, Avail_potet, Perv_frac, Perv_actet
!        Soil_moist = Soil_moist + Perv_actet - Avail_potet
!        Perv_actet = Avail_potet
      ENDIF
      IF ( Perv_actet*Perv_frac>Potet ) THEN
        PRINT *, 'perv_et PET problem', Perv_actet*Perv_frac, Avail_potet, Perv_frac, Potet
      ENDIF

      END SUBROUTINE compute_szactet

!***********************************************************************
!     compute interflow and flow to groundwater reservoir
!***********************************************************************
      SUBROUTINE compute_gwflow(Ssr2gw_rate, Ssr2gw_exp, Ssr_to_gw, Slow_stor)
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Ssr2gw_rate, Ssr2gw_exp
      REAL, INTENT(INOUT) :: Slow_stor, Ssr_to_gw
!***********************************************************************
!******compute flow to groundwater
      Ssr_to_gw = Ssr2gw_rate*(Slow_stor**Ssr2gw_exp)
      IF ( Ssr_to_gw<0.0 ) THEN
        Ssr_to_gw = 0.0
      ELSEIF ( Ssr_to_gw>Slow_stor ) THEN
        Ssr_to_gw = Slow_stor
      ENDIF
      Slow_stor = Slow_stor - Ssr_to_gw

      END SUBROUTINE compute_gwflow

!***********************************************************************
!     Compute subsurface lateral flow
!***********************************************************************
      SUBROUTINE compute_interflow(Coef_lin, Coef_sq, Ssres_in, Storage, Inter_flow)
      USE PRMS_CONSTANTS, ONLY: ERROR_soilzone !, NEARZERO, CLOSEZERO
      IMPLICIT NONE
! Functions
      INTRINSIC :: EXP, SQRT
      EXTERNAL :: error_stop
! Arguments
      REAL, INTENT(IN) :: Coef_lin, Coef_sq, Ssres_in
      REAL, INTENT(INOUT) :: Storage, Inter_flow
! Local Variables
      REAL :: c1, c2, c3, sos
!***********************************************************************
! Inter_flow is in inches for the timestep
!******compute interflow
      IF ( Coef_lin<=0.0 .AND. Ssres_in<=0.0 ) THEN
        c1 = Coef_sq*Storage
        Inter_flow = Storage*(c1/(1.0+c1))
      ELSEIF ( Coef_lin>0.0 .AND. Coef_sq<=0.0 ) THEN
        c2 = 1.0 - EXP(-Coef_lin)
        Inter_flow = Ssres_in*(1.0-c2/Coef_lin) + Storage*c2
      ELSEIF ( Coef_sq>0.0 ) THEN
        c3 = SQRT(Coef_lin**2.0+4.0*Coef_sq*Ssres_in)
        sos = Storage - ((c3-Coef_lin)/(2.0*Coef_sq))
        IF ( c3==0.0 ) CALL error_stop('in compute_interflow sos=0, please contact code developers', ERROR_soilzone)
        c1 = Coef_sq*sos/c3
        c2 = 1.0 - EXP(-c3)
        IF ( 1.0+c1*c2>0.0 ) THEN
          Inter_flow = Ssres_in + (sos*(1.0+c1)*c2)/(1.0+c1*c2)
        ELSE
          Inter_flow = Ssres_in
        ENDIF
      ELSE
        Inter_flow = 0.0
      ENDIF

! sanity check
      IF ( Inter_flow<0.0 ) THEN
!        IF ( Inter_flow<-NEARZERO ) PRINT *, 'interflow<0', Inter_flow, Ssres_in, Storage
        Inter_flow = 0.0
      ELSEIF ( Inter_flow>Storage ) THEN
        Inter_flow = Storage
      ENDIF
      Storage = Storage - Inter_flow
!      IF ( Storage<0.0 ) THEN
!        IF ( Storage<-CLOSEZERO ) PRINT *, 'Sanity check, ssres_stor<0.0', Storage
!        Storage = 0.0
! rsr, if very small storage, add it to interflow
!      ELSEIF ( Storage>0.0 .AND. Storage<NEARZERO ) THEN
!        print *, 'small storage', storage, inter_flow
!        Inter_flow = Inter_flow + Storage
!        Storage = 0.0
!      ENDIF

      END SUBROUTINE compute_interflow

!***********************************************************************
!     Compute cascading interflow and excess flow
!***********************************************************************
      SUBROUTINE compute_cascades(Ihru, Ncascade_hru, Slowflow, Preflow, &
     &           Dunnian, Dnslowflow, Dnpreflow, Dndunnflow)
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_SOILZONE, ONLY: Upslope_dunnianflow, Upslope_interflow
      USE PRMS_CASCADE, ONLY: Hru_down, Hru_down_frac, Hru_down_fracwt, Cascade_area
      USE PRMS_SRUNOFF, ONLY: Strm_seg_in
      IMPLICIT NONE
! Functions
      INTRINSIC :: IABS, DBLE
! Arguments
      INTEGER, INTENT(IN) :: Ihru, Ncascade_hru
      REAL, INTENT(INOUT) :: Dunnian, Slowflow, Preflow
      REAL, INTENT(INOUT) :: Dnslowflow, Dnpreflow, Dndunnflow
! Local Variables
      INTEGER :: j, k
      REAL :: frac, fracwt
!***********************************************************************
      DO k = 1, Ncascade_hru
        j = Hru_down(k, Ihru)
        frac = Hru_down_frac(k, Ihru)
! if hru_down(k, Ihru) > 0, cascade contributes to a downslope HRU
        IF ( j>0 ) THEN
          fracwt = Hru_down_fracwt(k, Ihru)
          Upslope_interflow(j) = Upslope_interflow(j) + DBLE( (Slowflow+Preflow)*fracwt )
          Upslope_dunnianflow(j) = Upslope_dunnianflow(j) + DBLE( Dunnian*fracwt )
          Dnslowflow = Dnslowflow + Slowflow*frac
          Dnpreflow = Dnpreflow + Preflow*frac
          Dndunnflow = Dndunnflow + Dunnian*frac
! if hru_down(k, Ihru) < 0, cascade contributes to a stream
        ELSEIF ( j<0 ) THEN
          j = IABS(j)
          Strm_seg_in(j) = Strm_seg_in(j) + DBLE( (Slowflow+Preflow+Dunnian)*Cascade_area(k, Ihru) )*Cfs_conv
        ENDIF
      ENDDO

! reset Slowflow, Preflow, and Dunnian_flow as they accumulate flow to streams
      Slowflow = Slowflow - Dnslowflow
      Preflow = Preflow - Dnpreflow
      Dunnian = Dunnian - Dndunnflow

      END SUBROUTINE compute_cascades

!***********************************************************************
!     compute interflow and flow to groundwater reservoir
!***********************************************************************
      SUBROUTINE compute_gravflow(Ihru, Capacity, Slowcoef_lin, &
     &           Slowcoef_sq, Ssr2gw_rate, Ssr2gw_exp, Gvr_maxin, &
     &           Pref_flow_thrsh, Gvr2pfr, Ssr_to_gw, &
     &           Slow_flow, Slow_stor, Gvr2sm, Soil_to_gw, Gwin, Compute_lateral)
      USE PRMS_CONSTANTS, ONLY: DEBUG_less, ACTIVE
      USE PRMS_MODULE, ONLY: Dprst_flag, Print_debug
      USE PRMS_SOILZONE, ONLY: Gravity_stor_res, Sm2gw_grav, Hru_gvr_count, Hru_gvr_index, &
     &    Gw2sm_grav, Gvr_hru_pct_adjusted
      USE PRMS_SRUNOFF, ONLY: Dprst_seep_hru
      IMPLICIT NONE
! Functions
      INTRINSIC :: MAX, DBLE, SNGL
      EXTERNAL :: check_gvr_sm, compute_interflow
! Arguments
      INTEGER, INTENT(IN) :: Ihru, Compute_lateral
      REAL, INTENT(IN) :: Slowcoef_lin, Slowcoef_sq, Ssr2gw_rate, Ssr2gw_exp
      REAL, INTENT(IN) :: Pref_flow_thrsh, Soil_to_gw, Gvr_maxin
      REAL, INTENT(INOUT) :: Capacity
      REAL, INTENT(OUT) :: Ssr_to_gw, Slow_stor, Slow_flow, Gvr2pfr, Gvr2sm
      DOUBLE PRECISION, INTENT(OUT) :: Gwin
! Local Variables
      INTEGER :: j, igvr
      REAL :: perc, slowflow, extra_water, gvrin_actual, depth, input
      DOUBLE PRECISION :: topfr, slflow, togw, slowstor, frac
!***********************************************************************
      !Capacity is for whole HRU
      !Soil_to_gw is for whole HRU
      !TO DO
! use VKS as a function of slope (vector analysis) instead of coef_lin
! coef_lin for pref_flow needs to be VKS lateral times a factor
! change slow to interflow
! in init, set an array dimensioned by nhrucell to vks*mfl_to_inch

      Gwin = 0.0D0
      Gvr2sm = 0.0
      topfr = 0.0D0
      slflow = 0.0D0
      togw = 0.0D0
      slowstor = 0.0D0
      DO j = 1, Hru_gvr_count(Ihru)
        igvr = Hru_gvr_index(j, Ihru)
        frac = Gvr_hru_pct_adjusted(igvr)
        Gwin = Gwin + DBLE( Gw2sm_grav(igvr) )*frac
        input = Gvr_maxin + Gw2sm_grav(igvr)
        depth = Gravity_stor_res(igvr) + input
        IF ( depth>0.0 .AND. Capacity>0.0 ) CALL check_gvr_sm(Capacity, depth, frac, Gvr2sm, input)

        IF ( Compute_lateral==ACTIVE ) THEN
          extra_water = MAX( 0.0, depth-Pref_flow_thrsh )
          IF ( extra_water>0.0 ) THEN
            !compute contribution to preferential-flow reservoir storage
            topfr = topfr + DBLE( extra_water )*frac
            depth = Pref_flow_thrsh
          ENDIF
          gvrin_actual = MAX(0.0, input-extra_water)

! compute contribution to slow interflow, if any
          IF ( depth>0.0 ) THEN
            CALL compute_interflow(Slowcoef_lin, Slowcoef_sq, gvrin_actual, depth, slowflow)
            slflow = slflow + DBLE( slowflow )*frac
          ENDIF
        ENDIF

! compute flow to groundwater, if any
        IF ( depth>0.0 ) THEN
          IF ( Ssr2gw_rate>0.0 ) THEN
! use VKS instead of rate  ???????????????
            perc = Ssr2gw_rate*(depth**Ssr2gw_exp)
            IF ( perc<0.0 ) THEN
              perc = 0.0
            ELSEIF ( perc>depth ) THEN
              perc = depth
            ENDIF
            depth = depth - perc
            Sm2gw_grav(igvr) = perc
            togw = togw + DBLE( perc )*frac
          ENDIF
!        ELSE ! GVRs can go negative if flux change in MODFLOW final iteration decreases, so don't set to 0
!          if(depth<0.0) print *, 'depth<0', depth, ihru
!          depth = 0.0
        ENDIF

        Gravity_stor_res(igvr) = depth
        slowstor = slowstor + DBLE(depth)*frac

! add any direct recharge from soil infiltration
        Sm2gw_grav(igvr) = Sm2gw_grav(igvr) + Soil_to_gw
        IF ( Dprst_flag==1 ) Sm2gw_grav(igvr) = Sm2gw_grav(igvr) + SNGL( Dprst_seep_hru(Ihru) )

      ENDDO ! end loop of GVRs in the HRU

      Gvr2pfr = SNGL( topfr )
      Slow_flow = SNGL( slflow )
      Ssr_to_gw = SNGL( togw )
      Slow_stor = SNGL( slowstor )
      IF ( Slow_stor>Pref_flow_thrsh ) THEN
        IF ( Print_debug>DEBUG_less .AND. Compute_lateral==ACTIVE ) &
     &       PRINT *, 'slow_stor > thrsh', Slow_stor, Pref_flow_thrsh, ' HRU:', Ihru
      ENDIF

      END SUBROUTINE compute_gravflow

!***********************************************************************
!     adjust soil moist based on being below field capacity (capacity)
!     and preferential-flow threshold (Pref_flow_thrsh)
!***********************************************************************
      SUBROUTINE check_gvr_sm(Capacity, Depth, Frac, Gvr2sm, Input)
!      USE PRMS_CONSTANTS, ONLY: CLOSEZERO
      IMPLICIT NONE
! Functions
      INTRINSIC :: MAX, ABS, SNGL
! Arguments
      DOUBLE PRECISION, INTENT(IN) :: Frac
      REAL, INTENT(INOUT) :: Capacity, Gvr2sm, Depth, Input
! Local Variables
      REAL :: to_sm, frac_sngl
!***********************************************************************
! check to see if soil is below capacity, if so add up to field capacity
! Capacity is for whole HRU
! to_sm and Gvr2sm are for whole HRU

      frac_sngl = SNGL( Frac )
      ! fill up capillary with part of gravity water
      to_sm = Capacity
      ! take all gravity water and put in capillary
      IF ( to_sm>Depth ) to_sm = Depth

! compute adjusmtent to soil moist to get to field capacity
      Capacity = Capacity - to_sm*frac_sngl
      IF ( Capacity<0.0 ) THEN
        to_sm = to_sm - Capacity*frac_sngl
        Capacity = 0.0
      ENDIF
      Gvr2sm = Gvr2sm + to_sm*frac_sngl
      Depth = Depth - to_sm
      !IF ( Depth<0.0 ) PRINT *, 'depth<0', depth
!      IF ( Depth<CLOSEZERO ) Depth = 0.0
      Input = Input - to_sm*frac_sngl

      END SUBROUTINE check_gvr_sm

!***********************************************************************
!     Initialize basin variables
!***********************************************************************
      SUBROUTINE init_basin_vars()
      USE PRMS_SOILZONE
      USE PRMS_FLOWVARS, ONLY: Basin_actet, Basin_perv_et, &
     &    Basin_swale_et, Basin_lakeevap, Basin_soil_to_gw, Basin_ssflow, Basin_sroff, &
     &    Basin_soil_moist, Basin_ssstor
      IMPLICIT NONE
!***********************************************************************
      Basin_lakeinsz = 0.0D0
      Basin_recharge = 0.0D0
      Basin_gvr2sm = 0.0D0
      Basin_sz_gwin = 0.0D0
      Basin_ssin = 0.0D0
      Basin_sm2gvr = 0.0D0
      Basin_dninterflow = 0.0D0
      Basin_dndunnianflow = 0.0D0
      Basin_dncascadeflow = 0.0D0
      Basin_sz2gw = 0.0D0
      Basin_sm2gvr_max = 0.0D0
      Basin_interflow_max = 0.0D0
      Basin_dunnian = 0.0D0
      Basin_capwaterin = 0.0D0
      Basin_cap_infil_tot = 0.0D0
      Basin_cap_up_max = 0.0D0
      Basin_pref_flow_infil = 0.0D0
      Basin_dunnian_pfr = 0.0D0
      Basin_dunnian_gvr = 0.0D0
      Basin_gvr2pfr = 0.0D0
      Basin_slowflow = 0.0D0
      Basin_prefflow = 0.0D0
      Basin_lakeprecip = 0.0D0
      Basin_actet = 0.0D0
      Basin_perv_et = 0.0D0
      Basin_swale_et = 0.0D0
      Basin_lakeevap = 0.0D0
      Basin_soil_to_gw = 0.0D0
      Basin_ssflow = 0.0D0
      Basin_sroff = 0.0D0
      Basin_soil_moist = 0.0D0
      Basin_slstor = 0.0D0
      Basin_ssstor = 0.0D0
      Basin_pref_stor = 0.0D0
      Basin_soil_rechr = 0.0D0
      Basin_soil_moist_tot = 0.0D0
      Basin_soil_lower_stor_frac = 0.0D0
      Basin_soil_rechr_stor_frac = 0.0D0
      Basin_sz_stor_frac = 0.0D0
      Basin_cpr_stor_frac = 0.0D0
      Basin_gvr_stor_frac = 0.0D0
      Basin_pfr_stor_frac = 0.0D0

      END SUBROUTINE init_basin_vars

!***********************************************************************
!     soilzone_restart - write or read soilzone restart file
!***********************************************************************
      SUBROUTINE soilzone_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit, GSFLOW_flag
      USE PRMS_SOILZONE
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Function
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=8) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_soil_rechr, Basin_slstor, Basin_soil_moist_tot, Basin_pref_stor
        WRITE ( Restart_outunit ) Pref_flow_stor
        IF ( GSFLOW_flag==ACTIVE ) WRITE ( Restart_outunit ) Gravity_stor_res
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_soil_rechr, Basin_slstor, Basin_soil_moist_tot, Basin_pref_stor
        READ ( Restart_inunit ) Pref_flow_stor
        IF ( GSFLOW_flag==ACTIVE ) READ ( Restart_inunit ) Gravity_stor_res
      ENDIF
      END SUBROUTINE soilzone_restart
