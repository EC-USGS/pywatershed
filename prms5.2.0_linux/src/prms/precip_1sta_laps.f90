!***********************************************************************
! Distributes precipitation and determines the form (rain, snow, mixed)
! to each HRU from one or more stations HRU using monthly correction
! factors to account for differences in altitude, spatial variation,
! topography, and measurement gage efficiency (precip_1sta)
! or by computing a daily lapse rate using elevation data and using
! specified adjustment factors for precipitation from two measurement
! stations(precip_laps)
!   Declared Parameters
!     tmax_allrain, tmax_allsnow, hru_psta, adjmix_rain, hru_area
!   Declared Parameters for precip_1sta
!     rain_adj = Rain_adj_lapse; snow_adj = Snow_adj_lapse
!   Declared Parameters for precip_laps
!     padj_rn, padj_sn
!     hru_plaps, psta_elev, pmn_mo, hru_elev
! Needs variable "precip" in the DATA FILE
! Needs computed variables tmaxf and tminf set in the temperature module
!***********************************************************************
      MODULE PRMS_PRECIP_1STA_LAPS
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, ACTIVE, OFF, GLACIER, &
     &      DEBUG_less, MM, MM2INCH, MONTHS_PER_YEAR, DOCUMENTATION, precip_1sta_module, precip_laps_module
        USE PRMS_MODULE, ONLY: Nhru, Nrain, Model, Process_flag, Inputerror_flag, Precip_flag, &
     &      Print_debug, Glacier_flag
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Precipitation Distribution'
        character(len=11) :: MODNAME
        character(len=*), parameter :: Version_precip = '2020-12-02'
        INTEGER, SAVE, ALLOCATABLE :: Psta_nuse(:)
        REAL, SAVE, ALLOCATABLE :: Rain_adj_lapse(:, :), Snow_adj_lapse(:, :), Precip_local(:)
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Hru_psta(:)
        ! Declared Parameters for precip_laps
        INTEGER, SAVE, ALLOCATABLE :: Hru_plaps(:)
        REAL, SAVE, ALLOCATABLE :: Padj_rn(:, :), Padj_sn(:, :), Pmn_mo(:, :)
      END MODULE PRMS_PRECIP_1STA_LAPS

      INTEGER FUNCTION precip_1sta_laps()
      USE PRMS_PRECIP_1STA_LAPS
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_area, Hru_route_order, Basin_area_inv, Hru_elev_ts, Hru_type
      USE PRMS_CLIMATEVARS, ONLY: Newsnow, Pptmix, Prmx, Basin_ppt, &
     &    Basin_rain, Basin_snow, Hru_ppt, Hru_rain, Hru_snow, &
     &    Basin_obs_ppt, Tmaxf, Tminf, Tmax_allrain_f, Tmax_allsnow_f, &
     &    Adjmix_rain, Precip_units
      USE PRMS_SET_TIME, ONLY: Nowmonth
      USE PRMS_OBS, ONLY: Precip
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, precip_form, print_module, compute_precip_laps
      EXTERNAL :: print_date, checkdim_param_limits
! Local Variables
      INTEGER :: i, ii, ierr
      REAL :: ppt
      DOUBLE PRECISION :: sum_obs
!***********************************************************************
      precip_1sta_laps = 0

      IF ( Process_flag==RUN ) THEN
        Precip_local = Precip
        DO i = 1, Nrain
          IF ( Psta_nuse(i)==ACTIVE ) THEN
            IF ( Precip_local(i)<0.0 ) THEN
              IF ( Print_debug>DEBUG_less ) THEN
                PRINT 9002, Precip_local(i), MODNAME, i
                CALL print_date(1)
              ENDIF
              Precip_local(i) = 0.0
            ELSEIF ( Precip_units==MM ) THEN
              Precip_local(i) = Precip_local(i)*MM2INCH
            ENDIF
          ENDIF
        ENDDO

        Basin_ppt = 0.0D0
        Basin_rain = 0.0D0
        Basin_snow = 0.0D0
        sum_obs = 0.0D0
        DO ii = 1, Active_hrus
          i = Hru_route_order(ii)
!******Zero precipitation on HRU
          Hru_ppt(i) = 0.0
          Hru_rain(i) = 0.0
          Hru_snow(i) = 0.0
          Prmx(i) = 0.0
          Newsnow(i) = OFF
          Pptmix(i) = OFF
          ppt = Precip_local(Hru_psta(i))
          IF ( Glacier_flag==ACTIVE ) THEN
            IF ( Hru_type(i)==GLACIER ) THEN
              ! Hru_elev_ts is the antecedent glacier elevation
              IF ( Precip_flag==precip_laps_module ) CALL compute_precip_laps(i, Hru_plaps(i), Hru_psta(i), Hru_elev_ts(i))
            ENDIF
          ENDIF
          IF ( ppt>0.0 ) &
     &         CALL precip_form(ppt, Hru_ppt(i), Hru_rain(i), Hru_snow(i), Tmaxf(i), &
     &                          Tminf(i), Pptmix(i), Newsnow(i), Prmx(i), &
     &                          Tmax_allrain_f(i,Nowmonth), Rain_adj_lapse(i,Nowmonth), Snow_adj_lapse(i,Nowmonth), &
     &                          Adjmix_rain(i,Nowmonth), Hru_area(i), sum_obs, Tmax_allsnow_f(i,Nowmonth))
        ENDDO
        Basin_ppt = Basin_ppt*Basin_area_inv
        Basin_obs_ppt = sum_obs*Basin_area_inv
        Basin_rain = Basin_rain*Basin_area_inv
        Basin_snow = Basin_snow*Basin_area_inv

      ELSEIF ( Process_flag==DECL ) THEN
        IF ( Precip_flag==precip_1sta_module ) THEN
          MODNAME = 'precip_1sta'
        ELSE
          MODNAME = 'precip_laps'
        ENDIF
        CALL print_module(MODDESC, MODNAME, Version_precip)

        ALLOCATE ( Psta_nuse(Nrain), Precip_local(Nrain) )

! Declare parameters
        ALLOCATE ( Hru_psta(Nhru) )
        IF ( declparam(MODNAME, 'hru_psta', 'nhru', 'integer', &
     &       '0', 'bounded', 'nrain', &
     &       'Index of base precipitation station for each HRU', &
     &       'Index of the base precipitation station used for lapse'// &
     &       ' rate calculations for each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'hru_psta')

        ALLOCATE ( Rain_adj_lapse(Nhru, MONTHS_PER_YEAR), Snow_adj_lapse(Nhru, MONTHS_PER_YEAR) )
        IF ( Precip_flag==precip_1sta_module .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'rain_adj', 'nhru,nmonths', 'real', &
     &         '1.0', '0.5', '10.0', &
     &         'Monthly rain adjustment factor for each HRU', &
     &         'Monthly (January to December) factor to adjust measured'// &
     &         ' precipitation on each HRU to account for'// &
     &         ' differences in elevation, and so forth', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'rain_adj')

          IF ( declparam(MODNAME, 'snow_adj', 'nhru,nmonths', 'real', &
     &         '1.0', '0.5', '2.5', &
     &         'Monthly snow adjustment factor for each HRU', &
     &         'Monthly (January to December) factor to adjust measured'// &
     &         ' precipitation on each HRU to account for'// &
     &         ' differences in elevation, and so forth', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'snow_adj')
        ENDIF

        IF ( Precip_flag==precip_laps_module .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Padj_rn(Nrain, MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'padj_rn', 'nrain,nmonths', 'real', &
     &         '1.0', '-2.0', '10.0', &
     &         'Rain adjustment factor, by month for each precipitation station', &
     &         'Monthly (January to December) factor to adjust'// &
     &         ' precipitation lapse rate computed between station hru_psta'// &
     &         ' and station hru_plaps; positive factors are mutiplied'// &
     &         ' times the lapse rate and negative factors are made'// &
     &         ' positive and substituted for the computed lapse rate', &
     &         'precip_units')/=0 ) CALL read_error(1, 'padj_rn')

          ALLOCATE ( Padj_sn(Nrain, MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'padj_sn', 'nrain,nmonths', 'real', &
     &         '1.0', '-2.0', '10.0', &
     &         'Snow adjustment factor, by month for each precipitation station', &
     &         'Monthly (January to December) factor to adjust'// &
     &         ' precipitation lapse rate computed between station hru_psta'// &
     &         ' and station hru_plaps; positive factors are mutiplied'// &
     &         ' times the lapse rate and negative factors are made'// &
     &         ' positive and substituted for the computed lapse rate', &
     &         'precip_units')/=0 ) CALL read_error(1, 'padj_sn')

          ALLOCATE ( Pmn_mo(Nrain, MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'pmn_mo', 'nrain,nmonths', 'real', &
     &         '1.0', '0.00001', '100.0', &
     &         'Mean monthly precipitation for each lapse precipitation station', &
     &         'Mean monthly (January to December) precipitation for'// &
     &         ' each lapse precipitation measurement station', &
     &         'precip_units')/=0 ) CALL read_error(1, 'pmn_mo')

          ALLOCATE ( Hru_plaps(Nhru) )
          IF ( declparam(MODNAME, 'hru_plaps', 'nhru', 'integer', &
     &         '0', 'bounded', 'nrain', &
     &         'Index of precipitation station to lapse against hru_psta', &
     &         'Index of the lapse precipitation measurement station used for lapse'// &
     &         ' rate calculations for each HRU', &
     &         'none')/=0 ) CALL read_error(1, 'hru_plaps')
        ENDIF

! Get parameters
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'hru_psta', Nhru, 'integer', Hru_psta)/=0 ) CALL read_error(2, 'hru_psta')

        IF ( Precip_flag==precip_1sta_module ) THEN
          IF ( getparam(MODNAME, 'rain_adj', Nhru*MONTHS_PER_YEAR, 'real', Rain_adj_lapse)/=0 ) CALL read_error(2, 'rain_adj')
          IF ( getparam(MODNAME, 'snow_adj', Nhru*MONTHS_PER_YEAR, 'real', Snow_adj_lapse)/=0 ) CALL read_error(2, 'snow_adj')
        ELSE
          IF ( getparam(MODNAME, 'padj_rn', Nrain*MONTHS_PER_YEAR, 'real', Padj_rn)/=0 ) CALL read_error(2, 'padj_rn')
          IF ( getparam(MODNAME, 'padj_sn', Nrain*MONTHS_PER_YEAR, 'real', Padj_sn)/=0 ) CALL read_error(2, 'padj_sn')
          IF ( getparam(MODNAME, 'hru_plaps', Nhru, 'integer', Hru_plaps)/=0 ) CALL read_error(2, 'hru_plaps')
          IF ( getparam(MODNAME, 'pmn_mo', Nrain*MONTHS_PER_YEAR, 'real', Pmn_mo)/=0 ) CALL read_error(2, 'pmn_mo')
        ENDIF

        Psta_nuse = 0
        DO ii = 1, Active_hrus
          i = Hru_route_order(ii)
          ierr = 0
          CALL checkdim_param_limits(i, 'hru_psta', 'nrain', Hru_psta(i), 1, Nrain, ierr)
          IF ( ierr==0 ) THEN
            Psta_nuse(Hru_psta(i)) = 1
          ELSE
            Inputerror_flag = 1
          ENDIF
          IF ( Precip_flag==precip_laps_module ) THEN
            ierr = 0
            CALL checkdim_param_limits(i, 'hru_plaps', 'nrain', Hru_plaps(i), 1, Nrain, ierr)
            IF ( ierr==0 ) THEN
              ! Hru_elev_ts is the current elevation, either hru_elev or for restart Hru_elev_ts
              CALL compute_precip_laps(i, Hru_plaps(i), Hru_psta(i), Hru_elev_ts(i))
            ELSE
              Inputerror_flag = 1
            ENDIF
          ENDIF
        ENDDO

      ENDIF

 9002 FORMAT (/, 'WARNING: negative precipitation value:', F10.3, /, 'specified for module ', A, /, &
     &        'precipitation station: ', I0, '; value set to 0.0')

      END FUNCTION precip_1sta_laps

!***********************************************************************
!     Compute lapse rate for an HRU
!***********************************************************************
      SUBROUTINE compute_precip_laps(Ihru, Hru_plaps, Hru_psta, Hru_elev)
      USE PRMS_CONSTANTS, ONLY: NEARZERO, MONTHS_PER_YEAR
      USE PRMS_PRECIP_1STA_LAPS, ONLY: Pmn_mo, Padj_sn, Padj_rn, Snow_adj_lapse, Rain_adj_lapse
      USE PRMS_CLIMATEVARS, ONLY: Psta_elev
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Ihru, Hru_psta, Hru_plaps
      REAL, INTENT(IN) :: Hru_elev
! Functions
      INTRINSIC ABS
! Local Variables
      INTEGER :: j
      REAL :: elp_diff, elh_diff, pmo_diff, pmo_rate, adj_p
!***********************************************************************
      elp_diff = Psta_elev(Hru_plaps) - Psta_elev(Hru_psta)
      IF ( ABS(elp_diff)<NEARZERO ) elp_diff = 1.0
      elh_diff = Hru_elev - Psta_elev(Hru_psta)
      DO j = 1, MONTHS_PER_YEAR
        pmo_diff = Pmn_mo(Hru_plaps, j) - Pmn_mo(Hru_psta, j)
        pmo_rate = pmo_diff / elp_diff
        adj_p = (pmo_rate*elh_diff)/Pmn_mo(Hru_psta, j)
        IF ( Padj_sn(Hru_psta, j)>=0.0 ) THEN
          Snow_adj_lapse(Ihru, j) = 1.0 + Padj_sn(Hru_psta, j)*adj_p
        ELSE
          Snow_adj_lapse(Ihru, j) = -Padj_sn(Hru_psta, j)
        ENDIF
        IF ( Padj_rn(Hru_psta,j)<0.0 ) THEN
          Rain_adj_lapse(Ihru, j) = -Padj_rn(Hru_psta, j)
        ELSE
          Rain_adj_lapse(Ihru, j) = 1.0 + Padj_rn(Hru_psta, j)*adj_p
        ENDIF
      ENDDO
      END SUBROUTINE compute_precip_laps
