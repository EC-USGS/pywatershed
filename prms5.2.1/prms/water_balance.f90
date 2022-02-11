!***********************************************************************
! Computes water balance for components of a PRMS model
!***********************************************************************
      MODULE PRMS_WATER_BALANCE
        IMPLICIT NONE
!   Local Variables
        character(len=*), parameter :: MODDESC = 'Water Balance Computations'
        character(len=*), parameter :: MODNAME_WB = 'water_balance'
        character(len=*), parameter :: Version_water_balance = '2021-08-13'
        INTEGER, SAVE :: BALUNT, SZUNIT, GWUNIT, INTCPUNT, SROUNIT, SNOWUNIT
        REAL, PARAMETER :: TOOSMALL = 3.1E-05, SMALL = 1.0E-04, BAD = 1.0E-03
        DOUBLE PRECISION, PARAMETER :: DSMALL = 1.0D-04, DTOOSMALL = 1.0D-05
        DOUBLE PRECISION, SAVE :: Last_basin_gwstor, Basin_dprst_wb
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Hru_storage_ante(:), Gwstor_ante(:)
!   Declared Variables
        DOUBLE PRECISION, SAVE :: Basin_capillary_wb, Basin_gravity_wb, Basin_soilzone_wb
!        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Hru_runoff(:)
      END MODULE PRMS_WATER_BALANCE

!***********************************************************************
!***********************************************************************
      SUBROUTINE water_balance()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN
      USE PRMS_MODULE, ONLY: Process_flag
      USE PRMS_WATER_BALANCE
      IMPLICIT NONE
! Functions
      EXTERNAL :: water_balance_decl, water_balance_init, water_balance_run
!***********************************************************************
      IF ( Process_flag==RUN ) THEN
        CALL water_balance_run()
      ELSEIF ( Process_flag==DECL ) THEN
        CALL water_balance_decl()
      ELSEIF ( Process_flag==INIT ) THEN
        CALL water_balance_init()
      ELSEIF ( Process_flag==CLEAN ) THEN
        CLOSE ( SZUNIT )
        CLOSE ( BALUNT )
        CLOSE ( INTCPUNT )
        CLOSE ( GWUNIT )
        CLOSE ( SROUNIT )
        CLOSE ( SNOWUNIT )
      ENDIF

      END SUBROUTINE water_balance

!***********************************************************************
!***********************************************************************
      SUBROUTINE water_balance_decl()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, DOCUMENTATION, CASCADE_OFF
      USE PRMS_MODULE, ONLY: Model, Nhru, Cascade_flag, Dprst_flag
      USE PRMS_WATER_BALANCE
      USE PRMS_SRUNOFF, ONLY: MODNAME
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declvar
      EXTERNAL :: read_error, print_module, PRMS_open_module_file
!***********************************************************************
      CALL print_module(MODDESC, MODNAME_WB, Version_water_balance)

! Declare Variables
      IF ( declvar(MODNAME_WB, 'basin_capillary_wb', 'one', 1, 'double', &
     &     'Basin area-weighted average capillary reservoir storage', &
     &     'inches', Basin_capillary_wb)/=0 ) CALL read_error(3, 'basin_capillary_wb')

      IF ( declvar(MODNAME_WB, 'basin_gravity_wb', 'one', 1, 'double', &
     &     'Basin area-weighted average gravity reservoir storage', &
     &     'inches', Basin_gravity_wb)/=0 ) CALL read_error(3, 'basin_gravity_wb')

      IF ( declvar(MODNAME_WB, 'basin_soilzone_wb', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in soilzone reservoirs', &
     &     'inches', Basin_soilzone_wb)/=0 ) CALL read_error(3, 'basin_soilzone_wb')

!      ALLOCATE ( Hru_runoff(Nhru) )
!      IF ( declvar(MODNAME, 'hru_runoff', 'nhru', Nhru, 'double', &
!     &     'Total lateral flow leaving each HRU (includes cascading flow)', &
!     &     'inches', Hru_runoff)/=0 ) CALL read_error(3, 'hru_runoff')

      ALLOCATE ( Hru_storage_ante(Nhru) )

      IF ( Dprst_flag==ACTIVE ) THEN
        IF ( declvar(MODNAME_WB, 'basin_dprst_wb', 'one', 1, 'double', &
     &       'Basin area-weighted average surface-depression storage water balance', &
     &       'inches', Basin_dprst_wb)/=0 ) CALL read_error(3, 'basin_dprst_wb')
      ENDIF

      ALLOCATE ( Gwstor_ante(Nhru) )

      IF ( Model/=DOCUMENTATION ) THEN
        CALL PRMS_open_module_file(INTCPUNT, 'intcp.wbal')
        WRITE ( INTCPUNT, 9003 )

        CALL PRMS_open_module_file(SNOWUNIT, 'snowcomp.wbal')
        WRITE ( SNOWUNIT, 9007 )

        CALL PRMS_open_module_file(SROUNIT, MODNAME//'.wbal')
        IF ( Cascade_flag>CASCADE_OFF ) THEN
          WRITE ( SROUNIT, 9006 )
        ELSE
          WRITE ( SROUNIT, 9005 )
        ENDIF

        CALL PRMS_open_module_file(SZUNIT, 'soilzone.wbal')
        WRITE ( SZUNIT, 9001 )

        CALL PRMS_open_module_file(GWUNIT, 'gwflow.wbal')
        WRITE ( GWUNIT, 9004 )

        CALL PRMS_open_module_file(BALUNT, 'wbal.msgs')
      ENDIF

 9001 FORMAT ('    Date     Water Bal     bsmbal    last SM  soilmoist  last stor    SS stor    perv ET      sz2gw  interflow', &
     &        '    soil2gw    dunnian    soil in   lakeinsz   downflow   swale ET  pref flow   pfr dunn   pfr stor', &
     &        '  slow stor dunnian gvr lake evap')
 9003 FORMAT ('    Date     Water Bal     Precip     Netppt  Intcpevap  Intcpstor  last_stor changeover  net apply     apply')
 9004 FORMAT ('    Date     Water Bal last store  GWR store', &
              '   GW input    GW flow    GW sink GW upslope minarea_in   downflow')
 9005 FORMAT ('    Date     Water Bal     Robal      Sroff      Infil  Impervevap Impervstor Dprst_evap Dprst_seep', &
     &        '   Perv Sro Imperv Sro  Dprst Sro  CFGI Sro  Glacrmelt')
 9006 FORMAT ('    Date     Water Bal     Robal      Sroff      Infil  Impervevap Impervstor Dprst_evap Dprst_seep', &
     &        '   Perv Sro Imperv Sro  Dprst Sro  Sroffdown  Srofflake   CFGI Sro  Glacrmelt')
 9007 FORMAT ('    Date     Water Bal  Snowpack    Snowmelt   Snowevap  Snowcover' )

      END SUBROUTINE water_balance_decl

!***********************************************************************
!     szinit - Initialize soilzone module - get parameter values,
!              set initial values and check parameter values
!***********************************************************************
      SUBROUTINE water_balance_init()
      USE PRMS_WATER_BALANCE
      USE PRMS_FLOWVARS, ONLY: Gwres_stor
      USE PRMS_GWFLOW, ONLY: Basin_gwstor
      USE PRMS_BASIN, ONLY: Hru_storage
!***********************************************************************
      Basin_capillary_wb = 0.0D0
      Basin_gravity_wb = 0.0D0
      Basin_soilzone_wb = 0.0D0
      Basin_dprst_wb = 0.0D0
!      Hru_runoff = 0.0D0
      Last_basin_gwstor = Basin_gwstor
      Gwstor_ante = Gwres_stor
      Hru_storage_ante = Hru_storage

      END SUBROUTINE water_balance_init

!***********************************************************************
!     water_balance_run - Computes balance for each HRU and model domain
!***********************************************************************
      SUBROUTINE water_balance_run()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, NEARZERO, DNEARZERO, LAKE, CASCADE_OFF, CASCADEGW_OFF
      USE PRMS_MODULE, ONLY: Cascade_flag, Cascadegw_flag, Dprst_flag, Glacier_flag, Nowyear, Nowmonth, Nowday
      USE PRMS_WATER_BALANCE
      USE PRMS_BASIN, ONLY: Hru_route_order, Active_hrus, Hru_frac_perv, Hru_area_dble, Hru_perv, &
     &    Hru_type, Basin_area_inv, Dprst_area_max, Hru_percent_imperv, Dprst_frac, Cov_type, Hru_storage
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt, Basin_ppt, Hru_rain, Hru_snow, Newsnow, Pptmix
      USE PRMS_FLOWVARS, ONLY: Basin_soil_moist, Basin_ssstor, Soil_to_gw, Soil_to_ssr, &
     &    Infil, Soil_moist_max, Ssr_to_gw, Ssres_flow, Basin_soil_to_gw, Soil_moist, Ssres_stor, &
     &    Slow_flow, Basin_perv_et, Basin_ssflow, Basin_swale_et, Slow_stor, Ssres_in, Soil_rechr, &
     &    Basin_lakeevap, Sroff, Hru_actet, Pkwater_equiv, Gwres_stor, Dprst_vol_open, Dprst_vol_clos, Basin_sroff
      USE PRMS_SET_TIME, ONLY: Nowtime
      USE PRMS_CASCADE, ONLY: Ncascade_hru
      USE PRMS_INTCP, ONLY: Hru_intcpstor, Basin_net_ppt, Basin_intcp_evap, Basin_changeover, &
     &    Basin_intcp_stor, Net_rain, Net_snow, Hru_intcpevap, &
     &    Srain_intcp, Wrain_intcp, Snow_intcp, Intcp_stor, Intcp_evap, &
     &    Canopy_covden, Intcp_changeover, Net_ppt, Intcp_stor_ante, Last_intcp_stor, &
     &    Net_apply, Gain_inches, Use_transfer_intcp, Basin_hru_apply, Basin_net_apply
      USE PRMS_SNOW, ONLY: Snowmelt, Pptmix_nopack, Snow_evap, Snowcov_area, Basin_pweqv, &
     &    Basin_snowmelt, Basin_snowevap, Basin_snowcov, Pkwater_ante, Glacrb_melt
      USE PRMS_GLACR, ONLY: Glacr_flow
      USE PRMS_SRUNOFF, ONLY: Basin_infil, Hru_hortn_cascflow, Upslope_hortonian, Hru_impervstor, &
     &    Dprst_stor_hru, Basin_sroffp, Basin_sroffi, Basin_dprst_sroff, Basin_sroff_down, &
     &    Basin_hortonian_lakes, Basin_imperv_evap, Basin_imperv_stor, &
     &    Basin_dprst_evap, Basin_dprst_seep, Hru_impervevap, Dprst_seep_hru, &
     &    Dprst_evap_hru, Dprst_sroff_hru, Dprst_insroff_hru, Sro_to_dprst_perv, &
     &    Dprst_area_clos, Hortonian_flow, Dprst_in, Hru_sroffp, Hru_sroffi, Imperv_stor_ante, &
     &    Dprst_stor_ante, Basin_cfgi_sroff
      USE PRMS_SOILZONE, ONLY: Swale_actet, Dunnian_flow, Basin_sz2gw, &
     &    Perv_actet, Cap_infil_tot, Pref_flow_infil, Cap_waterin, Upslope_interflow, &
     &    Upslope_dunnianflow, Pref_flow, Pref_flow_stor, Soil_lower, Gvr2pfr, Basin_ssin, &
     &    Basin_lakeinsz, Basin_dunnian, Pref_flow_max, Pref_flow_den, Pref_flow_thrsh, &
     &    Basin_sm2gvr_max, Basin_cap_infil_tot, Basin_slowflow, &
     &    Basin_dunnian_gvr, Basin_pref_flow_infil, Basin_dninterflow, Basin_pref_stor, Basin_dunnian_pfr, &
     &    Basin_dncascadeflow, Basin_capwaterin, Basin_sm2gvr, Basin_prefflow, Basin_slstor, Basin_gvr2pfr, &
     &    Hru_sz_cascadeflow, Pfr_dunnian_flow, Grav_dunnian_flow, Basin_dndunnianflow, &
     &    Soil_moist_tot, Soil_moist_ante, Ssres_stor_ante, Last_soil_moist, Last_ssstor
      USE PRMS_GWFLOW, ONLY: Basin_dnflow, Basin_gwsink, Basin_gwstor_minarea_wb, Gwres_flow, &
     &    Basin_gwstor, Basin_gwflow, Basin_gw_upslope, Basin_gwin, Gwres_sink, Hru_gw_cascadeflow, Gw_upslope, &
     &    Gwminarea_flag, Gwstor_minarea_wb, Gwin_dprst, Gwres_in
      IMPLICIT NONE
! Functions
      INTRINSIC ABS, DBLE, SNGL, DABS
      EXTERNAL print_date
! Local Variables
      INTEGER :: i, k
      REAL :: last_sm, last_ss, soilbal, perv_frac, gvrbal, test, waterin, waterout, hrubal
      REAL :: delstor, robal, gmelt
      DOUBLE PRECISION :: basin_bal, bsmbal, soil_in, gwbal, gwup, basin_robal, bsnobal
      DOUBLE PRECISION :: hru_out, hru_in, wbal, delta_stor, pptbal, brobal, dprst_hru_wb, harea
      CHARACTER(LEN=20), PARAMETER :: fmt1 = '(A, I5, 2("/",I2.2))'
!***********************************************************************
      Basin_capillary_wb = 0.0D0
      Basin_gravity_wb = 0.0D0
      basin_soilzone_wb = 0.0D0
      basin_bal = 0.0D0
      soil_in = 0.0D0
      basin_robal = 0.0D0
      bsnobal = 0.0D0
!      Hru_runoff = 0.0D0
      DO k = 1, Active_hrus
        i = Hru_route_order(k)

        IF ( Hru_type(i)==LAKE ) CYCLE ! no water balance for lakes

        harea = Hru_area_dble(i)
        perv_frac = Hru_frac_perv(i)

        ! intcp
        delstor = Hru_intcpstor(i) - Intcp_stor_ante(i)
        hrubal = Hru_rain(i) + Hru_snow(i) - Net_rain(i) - Net_snow(i) &
     &           - delstor - Hru_intcpevap(i) - Intcp_changeover(i)
        IF ( Use_transfer_intcp==ACTIVE ) hrubal = hrubal + Gain_inches(i) - Net_apply(i)
        IF ( ABS(hrubal)>TOOSMALL ) THEN
          IF ( ABS(hrubal)>SMALL ) THEN
            WRITE ( BALUNT, * ) 'Possible HRU interception water balance error'
          ELSE
            WRITE ( BALUNT, * ) 'Interception HRU rounding issue'
          ENDIF
          WRITE ( BALUNT,'(I7,6I5,15F10.5,I5)' ) i, Nowtime, hrubal, &
     &            Net_rain(i), Net_snow(i), Hru_rain(i), Hru_snow(i), &
     &            Intcp_stor(i), Intcp_stor_ante(i), Intcp_evap(i), Srain_intcp(i), &
     &            Wrain_intcp(i), Snow_intcp(i), Canopy_covden(i), delstor, &
     &            Hru_intcpstor(i), Intcp_changeover(i), Cov_type(i)
          IF ( Use_transfer_intcp==1 ) WRITE ( BALUNT, * ) Gain_inches(i), Net_apply(i)
        ENDIF

        ! Skip the HRU if there is no snowpack and no new snow
        IF ( Pkwater_ante(i)>DNEARZERO .OR. Newsnow(i)==ACTIVE ) THEN
          hrubal = SNGL( Pkwater_ante(i) - Pkwater_equiv(i) ) - Snow_evap(i) - Snowmelt(i)
          IF ( Pptmix_nopack(i)==ACTIVE ) THEN
            hrubal = hrubal + Net_snow(i)
          ELSE
            hrubal = hrubal + Net_ppt(i)
          ENDIF
          IF ( ABS(hrubal)>TOOSMALL ) THEN
            IF ( ABS(hrubal)>SMALL ) THEN
              WRITE ( BALUNT, * ) 'Possible snow water balance error'
            ELSE
              WRITE ( BALUNT, * ) 'Possible HRU snow rounding issue'
            ENDIF
            WRITE ( BALUNT,* ) i, hrubal, Nowyear, Nowmonth, Nowday, &
     &              Pkwater_ante(i), Pkwater_equiv(i), Snow_evap(i), &
     &              Snowmelt(i), Net_ppt(i), Net_snow(i), Net_rain(i), &
     &              Newsnow(i), Pptmix(i), Pptmix_nopack(i)
          ENDIF
          bsnobal = bsnobal + DBLE(hrubal)*harea
        ENDIF

        robal = Snowmelt(i) - Hortonian_flow(i) & !includes dprst runoff, if any
     &          - Infil(i)*perv_frac - Hru_impervevap(i) + Imperv_stor_ante(i) - Hru_impervstor(i) + Intcp_changeover(i)
        ! need to account for AG in water balance
        IF ( Use_transfer_intcp==ACTIVE ) robal = robal + Net_apply(i)*perv_frac !??? is net_apply for whole HRU (also for ag, impervious, dprst)
        IF ( Net_ppt(i)>0.0 ) THEN
          IF ( Pptmix_nopack(i)==ACTIVE ) THEN
            robal = robal + Net_rain(i)
          ELSEIF ( Snowmelt(i)<NEARZERO .AND. Pkwater_equiv(i)<DNEARZERO) THEN
            IF ( Snow_evap(i)<NEARZERO ) THEN
              robal = robal + Net_ppt(i)
            ELSEIF ( Net_snow(i)<NEARZERO ) THEN
              robal = robal + Net_rain(i)
            ENDIF
            !IF ( Net_snow(i)<NEARZERO ) robal = robal + Net_rain(i)
            !??  IF ( frzen==1 ) robal = robal + Net_rain(i)
          ENDIF
        ENDIF
        IF ( Cascade_flag>CASCADE_OFF ) robal = robal + SNGL( Upslope_hortonian(i) - Hru_hortn_cascflow(i) )
        IF ( Dprst_flag==ACTIVE ) robal = robal - Dprst_evap_hru(i) + &
     &                               SNGL( Dprst_stor_ante(i) - Dprst_stor_hru(i) - Dprst_seep_hru(i) ) !- Dprst_in(i) - Dprst_insroff_hru(i)
        gmelt = 0.0
        basin_robal = basin_robal + DBLE( robal )
        IF ( ABS(robal)>TOOSMALL ) THEN
          IF ( Dprst_flag==ACTIVE ) THEN
            dprst_hru_wb = Dprst_stor_ante(i) - Dprst_stor_hru(i) - Dprst_seep_hru(i) - Dprst_sroff_hru(i) + Dprst_in(i) &
     &                     - DBLE( Dprst_evap_hru(i) ) + DBLE( Dprst_insroff_hru(i) )
            Basin_dprst_wb = Basin_dprst_wb + dprst_hru_wb*harea
            WRITE ( BALUNT, * ) 'dprst', dprst_hru_wb, Dprst_stor_hru(i), Dprst_stor_hru(i), &
     &              Dprst_seep_hru(i), Dprst_evap_hru(i), Dprst_sroff_hru(i), Snowmelt(i), Net_rain(i), &
     &              Dprst_insroff_hru(i)
            WRITE ( BALUNT, * ) Dprst_vol_open(i), Dprst_vol_clos(i), &
     &              (Dprst_vol_open(i)+Dprst_vol_clos(i))/harea, Dprst_area_max(i), Pkwater_equiv(i), &
     &              Dprst_area_clos(i), Snowcov_area(i), Dprst_in(i), Sro_to_dprst_perv(i), Hru_sroffp(i), Hru_sroffi(i)
            WRITE ( BALUNT, * ) robal, Net_rain(i), Net_ppt(i), Net_rain(i)*Dprst_frac(i), &
     &              Dprst_frac(i), Pptmix_nopack(i)
            WRITE ( BALUNT, * ) Infil(i), perv_frac, Hru_impervevap(i), Imperv_stor_ante(i), Hru_impervstor(i), &
     &                          Hru_percent_imperv(i), Dprst_sroff_hru(i)
          ENDIF
          IF ( ABS(robal)>SMALL ) THEN
            WRITE ( BALUNT, '(A,I0,A,1X,I0)') 'Possible HRU surface runoff water balance ERROR, HRU: ', i, &
     &                                         ' hru_type:', Hru_type(i)
          ELSE
            WRITE ( BALUNT, '(A,I0,A,1X,I0)' ) 'HRU surface runoff rounding issue, HRU:', i, ' hru_type:', Hru_type(i)
          ENDIF
          IF ( Glacier_flag==ACTIVE ) gmelt = Glacrb_melt(i) + Glacr_flow(i)
          IF ( Cascade_flag>CASCADE_OFF ) THEN
            WRITE ( BALUNT, '(4I4,1X,F0.6,18(1X,F0.6))' ) Nowyear, Nowmonth, Nowday, Pptmix_nopack(i), robal, Snowmelt(i), &
     &              Upslope_hortonian(i), Imperv_stor_ante(i), Hru_hortn_cascflow(i), Infil(i), Hortonian_flow(i), &
     &              Hru_impervstor(i), Hru_impervevap(i), Net_ppt(i), &
     &              Pkwater_equiv(i), Snow_evap(i), Net_snow(i), Net_rain(i), Hru_sroffp(i), Hru_sroffi(i), gmelt, harea
          ELSE
            WRITE ( BALUNT,'(4I4,1X,F0.6,18(1X,F0.6))' ) Nowyear, Nowmonth, Nowday, Pptmix_nopack(i), &
     &              robal, Snowmelt(i), Imperv_stor_ante(i), Infil(i), Hortonian_flow(i), Hru_impervstor(i), &
     &              Hru_impervevap(i), Hru_percent_imperv(i), Net_ppt(i), Pkwater_equiv(i), Snow_evap(i), &
     &              Net_snow(i), Net_rain(i), Hru_sroffp(i), Hru_sroffi(i), Intcp_changeover(i), gmelt, harea
          ENDIF
        ENDIF

        last_sm = Soil_moist_ante(i)
        last_ss = Ssres_stor_ante(i)

        soilbal = (last_sm - Soil_moist(i) - Perv_actet(i))*perv_frac &
     &            - Soil_to_ssr(i) - Soil_to_gw(i) + Cap_infil_tot(i)

        IF ( ABS(soilbal)>TOOSMALL ) THEN
          WRITE ( BALUNT, * ) 'HRU capillary problem'
          WRITE ( BALUNT, * ) soilbal, Cap_infil_tot(i), last_sm, Soil_moist(i), Perv_actet(i), Soil_to_ssr(i), &
     &                        Soil_to_gw(i), i, Infil(i), Pref_flow_infil(i), perv_frac, &
     &                        Soil_moist_max(i), Cap_waterin(i)
          IF ( Cascade_flag>CASCADE_OFF ) WRITE ( BALUNT, * ) 'UP cascade', Upslope_interflow(i), Upslope_dunnianflow(i)
        ENDIF
        gvrbal = last_ss - Ssres_stor(i) + Soil_to_ssr(i) - Ssr_to_gw(i) - Swale_actet(i) - Dunnian_flow(i) &
     &           - Ssres_flow(i) + Pfr_dunnian_flow(i) + Pref_flow_infil(i)
        IF ( Cascade_flag>CASCADE_OFF ) gvrbal = gvrbal - Hru_sz_cascadeflow(i)
        test = ABS( gvrbal )
        IF ( test>TOOSMALL ) THEN
          WRITE ( BALUNT, * ) 'Bad GVR balance, HRU:', i, ' hru_type:', Hru_type(i)
          WRITE ( BALUNT, * ) gvrbal, last_ss, Ssres_stor(i), Ssr_to_gw(i), Swale_actet(i), &
     &            Dunnian_flow(i), Ssres_flow(i), Pfr_dunnian_flow(i), Pref_flow_thrsh(i), Ssres_in(i), &
     &            Pref_flow_infil(i), Grav_dunnian_flow(i), Slow_flow(i), Pref_flow(i), Soil_to_ssr(i), Gvr2pfr(i), &
     &            perv_frac, Slow_stor(i), Pref_flow_stor(i), Infil(i), Pref_flow_max(i), Pref_flow_den(i)
          IF ( Cascade_flag>CASCADE_OFF ) WRITE ( BALUNT, * ) 'sz cascade', Hru_sz_cascadeflow(i)
        ENDIF

        waterin = Cap_infil_tot(i) + Pref_flow_infil(i) + Pfr_dunnian_flow(i)
        waterout = Ssr_to_gw(i) + Ssres_flow(i) + Soil_to_gw(i) + Swale_actet(i) + Perv_actet(i)*perv_frac &
     &             + Dunnian_flow(i)
        IF ( Cascade_flag>CASCADE_OFF ) waterout = waterout + Hru_sz_cascadeflow(i)
        soil_in = soil_in + DBLE(Infil(i)*perv_frac)*harea
        soilbal = waterin - waterout + last_ss - Ssres_stor(i) + (last_sm-Soil_moist(i))*perv_frac
        basin_bal = basin_bal + DBLE(soilbal)*harea
        test = ABS( soilbal )
        IF ( test>TOOSMALL ) THEN
!          IF ( test>Ssres_stor(i)*TOOSMALL ) THEN
          WRITE ( BALUNT, * ) 'HRU:', i, ' Hru_type:', Hru_type(i)
          IF ( test>BAD ) THEN
            WRITE ( BALUNT, * ) 'HRU soilzone water balance ***ERROR***'
          ELSEIF ( test>SMALL ) THEN
            WRITE ( BALUNT, * ) 'Possible soilzone HRU water balance ERROR'
          ELSE
            WRITE ( BALUNT, * ) 'Possible soilzone HRU water balance rounding issue'
          ENDIF
          WRITE ( BALUNT, 9001 ) Nowyear, Nowmonth, Nowday, i, soilbal, Infil(i), last_sm, last_ss, &
     &            Soil_moist(i), Ssres_stor(i), Perv_actet(i), Ssr_to_gw(i), Slow_flow(i), Pref_flow(i), Ssres_flow(i), &
     &            Soil_to_gw(i), Pref_flow_infil(i), Pref_flow_stor(i), Slow_stor(i), Soil_rechr(i), &
     &            Soil_lower(i), Soil_to_ssr(i), Ssres_flow(i), waterin, Swale_actet(i)
          IF ( Cascade_flag>CASCADE_OFF ) WRITE ( BALUNT, * ) 'cascade', Upslope_dunnianflow(i), Upslope_interflow(i), &
     &                                                        Hru_sz_cascadeflow(i), Ncascade_hru(i)
          WRITE ( BALUNT, * ) Hru_perv(i), perv_frac, Pref_flow_den(i), (Infil(i)*perv_frac), Cap_infil_tot(i)
          WRITE ( BALUNT, * ) Dunnian_flow(i), Pfr_dunnian_flow(i)
!          ENDIF
        ENDIF

        hru_out = DBLE( Sroff(i) + Gwres_flow(i) + Ssres_flow(i) + Hru_actet(i) + Gwres_sink(i) )
        hru_in = DBLE( Hru_ppt(i) ) !+ Intcp_changeover(i)
        IF ( Cascade_flag>CASCADE_OFF ) THEN
          hru_out = hru_out + DBLE( Hru_sz_cascadeflow(i) ) + Hru_hortn_cascflow(i)
          hru_in = hru_in + Upslope_dunnianflow(i) + Upslope_interflow(i) + Upslope_hortonian(i)
        ENDIF
        IF ( Cascadegw_flag>CASCADEGW_OFF ) THEN
          hru_out = hru_out + Hru_gw_cascadeflow(i)
          hru_in = hru_in + Gw_upslope(i)/DBLE(harea)
        ENDIF
!        Hru_runoff(i) = hru_out - DBLE( Hru_actet(i) )
        wbal = hru_in - hru_out + Hru_storage_ante(i) - Hru_storage(i) !- Intcp_changeover(i)
        IF ( Gwminarea_flag==ACTIVE ) wbal = wbal + Gwstor_minarea_wb(i)
        IF ( DABS(wbal)>DTOOSMALL ) THEN
          WRITE ( BALUNT, * ) 'Possible HRU water balance issue:', wbal, '; HRU:', i, ' hru_type:', Hru_type(i), '; area:', harea
          WRITE ( BALUNT, * ) Sroff(i), Gwres_flow(i), Ssres_flow(i), Hru_actet(i), Gwres_sink(i), &
     &            Hru_storage_ante(i), Hru_storage(i), Hru_type(i),  Pfr_dunnian_flow(i), Intcp_changeover(i), &
     &            Snowmelt(i), Hru_ppt(i)
          !WRITE ( BALUNT, * ) Gwstor_minarea_wb(i)
          WRITE ( BALUNT, * ) Soil_moist_tot(i), Hru_intcpstor(i), Gwres_stor(i), Pkwater_equiv(i), Hru_impervstor(i)
          IF ( Cascade_flag>CASCADE_OFF ) THEN
            WRITE ( BALUNT, * ) Hru_sz_cascadeflow(i), Upslope_dunnianflow(i), Upslope_interflow(i)
            WRITE ( BALUNT, * ) Upslope_hortonian(i), Hru_hortn_cascflow(i)
          ENDIF
          IF ( Cascadegw_flag>CASCADEGW_OFF ) WRITE ( BALUNT, * ) Gw_upslope(i)/DBLE(harea), Hru_gw_cascadeflow(i)
          !CALL print_date(0)
          WRITE ( BALUNT, FMT1 ) '   Date:', Nowyear, Nowmonth, Nowday
        ENDIF

        wbal = Gwstor_ante(i) + Gwres_in(i)/harea - Gwres_stor(i) - DBLE( Gwres_sink(i) + Gwres_flow(i) )
        IF ( Cascadegw_flag>CASCADEGW_OFF ) wbal = wbal - Hru_gw_cascadeflow(i)
        IF ( Gwminarea_flag==ACTIVE ) wbal = wbal + Gwstor_minarea_wb(i)
        gwup = 0.0D0
        IF ( Cascadegw_flag>CASCADEGW_OFF ) gwup = Gw_upslope(i)
        IF ( DABS(wbal)>DTOOSMALL ) THEN
          WRITE ( BALUNT, * ) 'Possible GWR water balance issue', &
     &                        i, wbal, Gwstor_ante(i), Gwres_in(i)/harea, Gwres_stor(i), Gwres_flow(i), &
     &                        Gwres_sink(i), Soil_to_gw(i), Ssr_to_gw(i), gwup, harea
          IF ( Cascadegw_flag>CASCADEGW_OFF ) WRITE ( BALUNT, * ) 'gw cascade', Hru_gw_cascadeflow(i)
          IF ( Gwminarea_flag==ACTIVE ) WRITE ( BALUNT, * ) 'gwstor_minarea_wb', Gwstor_minarea_wb(i)
          IF ( Dprst_flag==ACTIVE ) WRITE ( BALUNT, * ) 'gwin_dprst', Gwin_dprst(i)
        ENDIF
        Hru_storage_ante(i) = Hru_storage(i)
        Gwstor_ante(i) = Gwres_stor(i)
      ENDDO
      Basin_dprst_wb = Basin_dprst_wb*Basin_area_inv

! intcp
      delta_stor = Basin_intcp_stor - Last_intcp_stor
      pptbal = Basin_ppt - Basin_net_ppt - delta_stor - Basin_intcp_evap - Basin_changeover
      IF ( Use_transfer_intcp==ACTIVE ) pptbal = pptbal + Basin_net_apply
      IF ( DABS(pptbal)>DSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'Possible basin interception water balance error', &
     &                         Nowyear, Nowmonth, Nowday, pptbal, Basin_changeover
      ELSEIF ( DABS(pptbal)>DTOOSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'Interception basin rounding issue', &
     &                         Nowyear, Nowmonth, Nowday, pptbal
      ENDIF
      WRITE ( INTCPUNT, 9002 ) Nowyear, Nowmonth, Nowday, pptbal, &
     &        Basin_ppt, Basin_net_ppt, Basin_intcp_evap, &
     &        Basin_intcp_stor, Last_intcp_stor, Basin_changeover, Basin_net_apply, Basin_hru_apply

! snowcomp
      bsnobal = bsnobal*Basin_area_inv
      IF ( DABS(bsnobal)>DSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'Possible basin snow water balance error', &
     &                         Nowyear, Nowmonth, Nowday, bsnobal
      ELSEIF ( DABS(bsnobal)>DTOOSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'Possible basin snow rounding issue', &
     &                         Nowyear, Nowmonth, Nowday, bsnobal
      ENDIF
      WRITE ( SNOWUNIT, 9002 ) Nowyear, Nowmonth, Nowday, bsnobal, Basin_pweqv, &
     &                         Basin_snowmelt, Basin_snowevap, Basin_snowcov

! srunoff
      brobal = Basin_sroff - Basin_sroffp - Basin_sroffi - Basin_dprst_sroff - Basin_cfgi_sroff
      IF ( Cascade_flag>CASCADE_OFF ) THEN
        brobal = brobal + Basin_sroff_down
        WRITE ( SROUNIT, 9002 ) Nowyear, Nowmonth, Nowday, basin_robal, &
     &          brobal, Basin_sroff, Basin_infil, Basin_imperv_evap, &
     &          Basin_imperv_stor, Basin_dprst_evap, Basin_dprst_seep, &
     &          Basin_sroffp, Basin_sroffi, Basin_dprst_sroff, &
     &          Basin_sroff_down, Basin_hortonian_lakes, Basin_cfgi_sroff
      ELSE
        WRITE ( SROUNIT, 9002 ) Nowyear, Nowmonth, Nowday, basin_robal, &
     &          brobal, Basin_sroff, Basin_infil, Basin_imperv_evap, &
     &          Basin_imperv_stor, Basin_dprst_evap, Basin_dprst_seep, &
     &          Basin_sroffp, Basin_sroffi, Basin_dprst_sroff, Basin_cfgi_sroff
      ENDIF
      IF ( DABS(basin_robal)>DSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'possible srunoff basin water balance ERROR', &
     &                         Nowyear, Nowmonth, Nowday, basin_robal
      ELSEIF ( DABS(basin_robal)>DTOOSMALL ) THEN
        WRITE ( BALUNT, 9003 ) 'possible srunoff basin water balance rounding issue', &
     &                         Nowyear, Nowmonth, Nowday, basin_robal
      ENDIF

! soilzone
      Basin_capillary_wb = Last_soil_moist - Basin_soil_moist - &
     &                     Basin_perv_et - Basin_sm2gvr_max + Basin_cap_infil_tot - Basin_soil_to_gw
      Basin_gravity_wb = Last_ssstor - Basin_ssstor + Basin_sm2gvr - Basin_dncascadeflow - &
     &                   Basin_ssflow - Basin_sz2gw - Basin_dunnian + Basin_dunnian_pfr - &
     &                   Basin_swale_et + Basin_pref_flow_infil
      Basin_soilzone_wb = Basin_infil + Last_ssstor - &
     &                    Basin_ssstor + Last_soil_moist - Basin_soil_moist - &
     &                    Basin_perv_et - Basin_swale_et - Basin_sz2gw - &
     &                    Basin_soil_to_gw - Basin_ssflow - Basin_dunnian - &
     &                    Basin_lakeinsz

      IF ( DABS(Basin_gravity_wb)>DTOOSMALL ) WRITE ( BALUNT, * ) 'basin gvrbal issue', Basin_gravity_wb, &
     &     Last_ssstor, Basin_ssstor, Basin_sm2gvr, Basin_ssflow, Basin_sz2gw, Basin_dunnian, &
     &     Basin_swale_et, Basin_pref_flow_infil, Basin_dninterflow, Basin_pref_stor, &
     &     Basin_dunnian_pfr, Basin_lakeinsz, Basin_dncascadeflow, &
     &     Basin_dunnian_gvr, Basin_slowflow, Basin_prefflow, Basin_gvr2pfr, Nowtime
      IF ( DABS(Basin_capillary_wb)>DTOOSMALL ) WRITE( BALUNT, * ) 'possible basin capillary balance issue', &
     &     Basin_capillary_wb, Last_soil_moist, Basin_soil_moist, Basin_perv_et, &
     &     Basin_sm2gvr, Basin_cap_infil_tot, Basin_soil_to_gw, Basin_sm2gvr_max, Basin_capwaterin, Nowtime
      IF ( DABS(Basin_soilzone_wb)>DTOOSMALL ) WRITE ( BALUNT, * ) 'possible basin soil zone rounding issue', &
     &     Basin_soilzone_wb, Basin_capwaterin, Basin_pref_flow_infil, Basin_infil, &
     &     Last_ssstor, Basin_ssstor, Last_soil_moist, Basin_soil_moist, Basin_perv_et, Basin_swale_et, &
     &     Basin_sz2gw, Basin_soil_to_gw, Basin_ssflow, Basin_dunnian, Basin_dncascadeflow, &
     &     Basin_sm2gvr, Basin_lakeinsz, Basin_dunnian_pfr, Nowtime

      soil_in = soil_in*Basin_area_inv
      basin_bal = basin_bal*Basin_area_inv
      bsmbal = Last_soil_moist - Basin_soil_moist + Last_ssstor - Basin_ssstor - Basin_perv_et - Basin_sz2gw + soil_in - &
     &         Basin_ssflow - Basin_soil_to_gw - Basin_dunnian - Basin_swale_et - Basin_lakeinsz

      WRITE ( SZUNIT, 9002 ) Nowyear, Nowmonth, Nowday, basin_bal, &
     &        bsmbal, Last_soil_moist, Basin_soil_moist, Last_ssstor, &
     &        Basin_ssstor, Basin_perv_et, Basin_sz2gw, Basin_ssflow, &
     &        Basin_soil_to_gw, Basin_dunnian, soil_in, Basin_lakeinsz, &
     &        Basin_dncascadeflow, Basin_swale_et, Basin_prefflow, Basin_dunnian_pfr, &
     &        Basin_pref_stor, Basin_slstor, Basin_dunnian_gvr, Basin_lakeevap

      IF ( DABS(bsmbal)>0.05D0 .OR. DABS(basin_bal)>0.001D0 ) THEN
        WRITE ( BALUNT, * ) '*ERROR, soilzone basin water balance', bsmbal, &
     &          basin_bal, Last_soil_moist, Basin_soil_moist, &
     &          Last_ssstor, Basin_ssstor, Basin_perv_et, Basin_sz2gw, &
     &          soil_in, Basin_ssflow, Basin_soil_to_gw, Basin_dunnian, &
     &          Basin_swale_et, Basin_lakeinsz
        WRITE ( BALUNT, * ) Basin_pref_stor, Basin_slstor
      ELSEIF ( DABS(bsmbal)>0.005D0 .OR. DABS(basin_bal)>DTOOSMALL ) THEN
        WRITE ( BALUNT, * ) 'Possible soilzone basin water balance ERROR', &
     &          bsmbal, basin_bal, Last_soil_moist, Basin_soil_moist, &
     &          Last_ssstor, Basin_ssstor, Basin_perv_et, Basin_sz2gw, &
     &          soil_in, Basin_ssflow, Basin_soil_to_gw, Basin_dunnian, &
     &          Basin_swale_et, Basin_lakeinsz
        WRITE ( BALUNT, * ) Basin_pref_stor, Basin_slstor
      ELSEIF ( DABS(bsmbal)>0.0005D0 .OR. DABS(basin_bal)>DTOOSMALL ) THEN
        WRITE ( BALUNT, '(A,2F12.7)' ) 'Basin soilzone rounding issue', bsmbal, basin_bal
        WRITE ( BALUNT, * ) Basin_soilzone_wb, Basin_ssin, &
     &          Basin_dninterflow, Basin_sm2gvr, Basin_capwaterin, &
     &          soil_in, Basin_gvr2pfr, Basin_dndunnianflow, (soil_in - Basin_infil)
      ENDIF

! gwflow

      ! not going to balance because gwstor under lakes is computed each time step fix for lakes
      ! basin_gwin includes upslope flow, gwin_dprst, soil_to_gw, ssr_to_gw
      gwbal = Basin_gwin + Last_basin_gwstor - Basin_gwstor - Basin_gwsink &
     &        - Basin_gwflow - Basin_dnflow + Basin_gwstor_minarea_wb
      IF ( DABS(gwbal)>DSMALL ) WRITE ( BALUNT, 9003 ) 'Possible GWR basin water balance issue', &
     &                                                 Nowyear, Nowmonth, Nowday, gwbal
      WRITE ( GWUNIT, 9002 ) Nowyear, Nowmonth, Nowday, &
     &        gwbal, Last_basin_gwstor, Basin_gwstor, Basin_gwin, &
     &        Basin_gwflow, Basin_gwsink, &
     &        Basin_gw_upslope, Basin_gwstor_minarea_wb, Basin_dnflow
      Last_basin_gwstor = Basin_gwstor

 9001 FORMAT (I5, 2('/', I2.2), I7, 26F11.5)
 9002 FORMAT (I5, 2('/', I2.2), 23F11.5)
 9003 FORMAT (A, I5, 2('/', I2.2), 2F12.5)

      END SUBROUTINE water_balance_run
