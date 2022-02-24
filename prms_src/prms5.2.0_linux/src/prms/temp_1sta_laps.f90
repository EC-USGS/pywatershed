!***********************************************************************
! Distributes maximum, minimum, and average temperatures to each HRU
! using temperature data measured at one station and an estimated monthly
! lapse rate (temp_1sta) or by computing a daily lapse rate based on
! elevations with temperature data measured at two stations (temp_laps)
!
! Variables needed from DATA FILE: tmax, tmin
! Declared Parameters
!     tsta_elev, tmax_adj, tmin_adj, hru_type
!     tmax_adj = tmax_aspect_adjust; tmin_adj = tmin_aspect_adjust
!     hru_tsta, hru_elev, hru_area, temp_units, basin_tsta
! Declared Parameters for temp_1sta
!     tmax_lapse, tmin_lapse
! Declared Parameters for temp_laps
!     hru_tlaps
!***********************************************************************
      MODULE PRMS_TEMP_1STA_LAPS
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, &
     &      GLACIER, DEBUG_less, MONTHS_PER_YEAR, ERROR_temp, DOCUMENTATION, &
     &      MINTEMP, MAXTEMP, NEARZERO, temp_1sta_module, temp_laps_module, READ_INIT, SAVE_INIT
        USE PRMS_MODULE, ONLY: Process_flag, Nhru, Ntemp, Model, &
     &      Print_debug, Init_vars_from_file, Save_vars_to_file, &
     &      Temp_flag, Inputerror_flag, Start_month, Glacier_flag
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Temperature Distribution'
        character(len=9), SAVE :: MODNAME
        character(len=*), parameter :: Version_temp = '2020-12-02'
        INTEGER, SAVE, ALLOCATABLE :: Tmax_cnt(:), Tmin_cnt(:), Nuse_tsta(:)
        REAL, SAVE, ALLOCATABLE :: Elfac(:), Tmax_prev(:), Tmin_prev(:)
        REAL, SAVE, ALLOCATABLE :: Tcrn(:), Tcrx(:) ! temp_1sta
        REAL, SAVE :: Solrad_tmax_good, Solrad_tmin_good
        ! Declared Parameters
        INTEGER, SAVE :: Max_missing
        REAL, SAVE, ALLOCATABLE :: Tmax_lapse(:, :), Tmin_lapse(:, :)
        INTEGER, SAVE, ALLOCATABLE :: Hru_tlaps(:)
      END MODULE PRMS_TEMP_1STA_LAPS

      INTEGER FUNCTION temp_1sta_laps()
      USE PRMS_TEMP_1STA_LAPS
      USE PRMS_BASIN, ONLY: Hru_elev_ts, Hru_area, Active_hrus, Hru_route_order, Basin_area_inv, Hru_type
      USE PRMS_CLIMATEVARS, ONLY: Tmax_aspect_adjust, Tmin_aspect_adjust, Tsta_elev, &
     &    Hru_tsta, Solrad_tmax, Solrad_tmin, Basin_temp, Basin_tmax, &
     &    Basin_tmin, Tmaxf, Tminf, Tminc, Tmaxc, Tavgf, Tavgc, Basin_tsta, Tmax_allrain
      USE PRMS_SET_TIME, ONLY: Nowmonth, Nowday
      USE PRMS_OBS, ONLY: Tmax, Tmin
      IMPLICIT NONE
! Functions
      INTRINSIC :: INDEX, ABS
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, temp_set, print_module, temp_1sta_laps_restart, print_date, checkdim_param_limits
      EXTERNAL :: compute_temp_laps
! Local Variables
      INTEGER :: j, k, jj, i, kk, kkk, l, ierr
      REAL :: tmx, tmn
!***********************************************************************
      temp_1sta_laps = 0

      IF ( Process_flag==RUN ) THEN
        kk = 0
        kkk = 0
        DO i = 1, Ntemp
          IF ( Nuse_tsta(i)>0 ) THEN
            IF ( Tmax(i)<MINTEMP .OR. Tmax(i)>MAXTEMP ) THEN
              Tmax_cnt(i) = Tmax_cnt(i) + 1
              IF ( Tmax_cnt(i)<Max_missing ) THEN
                IF ( Print_debug>DEBUG_less ) THEN
                  PRINT 9001, 'tmax', Tmax(i), i, Tmax_prev(i)
                  CALL print_date(0)
                ENDIF
                Tmax(i) = Tmax_prev(i)
                kk = 1
              ELSE
                PRINT 9002, 'tmax', Tmax(i), i
                CALL print_date(0)
                ERROR STOP ERROR_temp
              ENDIF
            ELSE
              Tmax_prev(i) = Tmax(i)
              Tmax_cnt(i) = 0
            ENDIF
            IF ( Tmin(i)<MINTEMP .OR. Tmin(i)>MAXTEMP ) THEN
              Tmin_cnt(i) = Tmin_cnt(i) + 1
              IF ( Tmin_cnt(i)<Max_missing ) THEN
                IF ( Print_debug>DEBUG_less ) THEN
                  PRINT 9001, 'tmin', Tmin(i), i, Tmin_prev(i)
                  CALL print_date(0)
                ENDIF
                Tmin(i) = Tmin_prev(i)
                kkk = 1
              ELSE
                PRINT 9002, 'tmin', Tmin(i), i
                CALL print_date(0)
                ERROR STOP ERROR_temp
              ENDIF
            ELSE
              Tmin_prev(i) = Tmin(i)
              Tmin_cnt(i) = 0
            ENDIF
          ENDIF
        ENDDO
        ! if all values good, reset _cnt variable
        IF ( kk==0 ) Tmax_cnt = 0
        IF ( kkk==0 ) Tmin_cnt = 0

        Basin_tmax = 0.0D0
        Basin_tmin = 0.0D0
        Basin_temp = 0.0D0
        IF ( Temp_flag==temp_1sta_module ) THEN
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            k = Hru_tsta(j)
            IF ( Nowday==1 ) THEN
              IF ( Glacier_flag==ACTIVE ) THEN
                ! Hru_elev_ts is the antecedent glacier elevation
                IF ( Hru_type(j)==GLACIER ) Elfac(j) = (Hru_elev_ts(j) - Tsta_elev(k))/1000.0
              ENDIF
              Tcrx(j) = Tmax_lapse(j, Nowmonth)*Elfac(j) - Tmax_aspect_adjust(j, Nowmonth)
              Tcrn(j) = Tmin_lapse(j, Nowmonth)*Elfac(j) - Tmin_aspect_adjust(j, Nowmonth)
            ENDIF
            tmx = Tmax(k) - Tcrx(j)
            tmn = Tmin(k) - Tcrn(j)
            CALL temp_set(j, tmx, tmn, Tmaxf(j), Tminf(j), Tavgf(j), &
     &                    Tmaxc(j), Tminc(j), Tavgc(j), Hru_area(j))
          ENDDO
        ELSEIF ( Temp_flag==temp_laps_module ) THEN
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            k = Hru_tsta(j)
            l = Hru_tlaps(j)
            IF ( Glacier_flag==ACTIVE ) THEN
              ! Hru_elev_ts is the antecedent glacier elevation
              IF ( Hru_type(j)==GLACIER ) CALL compute_temp_laps(Elfac(j), Hru_elev_ts(j), Tsta_elev(l), Tsta_elev(k))
            ENDIF
            tmx = Tmax(k) + (Tmax(l) - Tmax(k))*Elfac(j) + Tmax_aspect_adjust(j, Nowmonth)
            tmn = Tmin(k) + (Tmin(l) - Tmin(k))*Elfac(j) + Tmin_aspect_adjust(j, Nowmonth)
            CALL temp_set(j, tmx, tmn, Tmaxf(j), Tminf(j), Tavgf(j), &
     &                    Tmaxc(j), Tminc(j), Tavgc(j), Hru_area(j))
          ENDDO
        ELSE ! Temp_flag = temp_sta_module
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            k = Hru_tsta(j)
            tmx = Tmax(k) + Tmax_aspect_adjust(j, Nowmonth)
            tmn = Tmin(k) + Tmin_aspect_adjust(j, Nowmonth)
            CALL temp_set(j, tmx, tmn, Tmaxf(j), Tminf(j), Tavgf(j), &
     &                    Tmaxc(j), Tminc(j), Tavgc(j), Hru_area(j))
          ENDDO
        ENDIF
        Basin_tmax = Basin_tmax*Basin_area_inv
        Basin_tmin = Basin_tmin*Basin_area_inv
        Basin_temp = Basin_temp*Basin_area_inv
        Solrad_tmax = Tmax(Basin_tsta)
        Solrad_tmin = Tmin(Basin_tsta)
        IF ( Solrad_tmax<MINTEMP .OR. Solrad_tmax>MAXTEMP ) THEN
          IF ( Print_debug>DEBUG_less ) THEN
            PRINT *, 'Bad temperature data to set solrad_tmax:', Solrad_tmax, ' using last valid value:', Solrad_tmax_good
            CALL print_date(0)
          ENDIF
          Solrad_tmax = Solrad_tmax_good
        ELSE
          Solrad_tmax_good = Solrad_tmax
        ENDIF
        IF ( Solrad_tmin<MINTEMP .OR. Solrad_tmin>MAXTEMP ) THEN
          IF ( Print_debug>DEBUG_less ) THEN
            PRINT *, 'Bad temperature data to set solrad_tmin:', Solrad_tmin, ' using last valid value:', Solrad_tmin_good
            CALL print_date(0)
          ENDIF
          Solrad_tmin = Solrad_tmin_good
        ELSE
          Solrad_tmin_good = Solrad_tmin
        ENDIF

      ELSEIF ( Process_flag==DECL ) THEN
        IF ( Temp_flag==temp_1sta_module ) THEN
          MODNAME = 'temp_1sta'
        ELSEIF ( Temp_flag==temp_laps_module ) THEN
          MODNAME = 'temp_laps'
        ELSE ! Temp_flag = temp_sta_module
          MODNAME = 'temp_sta '
        ENDIF
        CALL print_module(MODDESC, MODNAME, Version_temp)

        ALLOCATE ( Elfac(Nhru), Nuse_tsta(Ntemp) )
        ALLOCATE ( Tmin_cnt(Ntemp), Tmax_cnt(Ntemp), Tmax_prev(Ntemp), Tmin_prev(Ntemp) )

        IF ( Temp_flag==temp_1sta_module .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Tcrn(Nhru), Tcrx(Nhru) )
          ALLOCATE ( Tmax_lapse(Nhru, MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'tmax_lapse', 'nhru,nmonths', 'real', &
     &         '3.0', '-20.0', '20.0', &
     &         'Monthly maximum temperature lapse rate for each HRU', &
     &         'Monthly (January to December) values representing the change in maximum air temperature per 1000 elev_units of'// &
     &         ' elevation change for each HRU', &
     &         'temp_units/elev_units')/=0 ) CALL read_error(1, 'tmax_lapse')
! 3 degC/ 1000 ft is adiabatic, or 9.8 degC/ 1000 m, or 5.4 degF/ 1000 ft, or 17.64 degF /1000 m

          ALLOCATE ( Tmin_lapse(Nhru, MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'tmin_lapse', 'nhru,nmonths', 'real', &
     &         '3.0', '-20.0', '20.0', &
     &         'Monthly minimum temperature lapse rate for each HRU', &
     &         'Monthly (January to December) values representing the change in minimum air temperture per 1000 elev_units of'// &
     &         ' elevation change for each HRU', &
     &         'temp_units/elev_units')/=0 ) CALL read_error(1, 'tmin_lapse')
        ENDIF

        IF ( Temp_flag==temp_laps_module .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Hru_tlaps(Nhru) )
          IF ( declparam(MODNAME, 'hru_tlaps', 'nhru', 'integer', &
     &         '0', 'bounded', 'ntemp', &
     &         'Index of lapse temperature station for each HRU', &
     &         'Index of the lapse temperature station used for lapse rate calculations', &
     &         'none')/=0 ) CALL read_error(1, 'hru_tlaps')
        ENDIF

        IF ( declparam(MODNAME, 'max_missing', 'one', 'integer', &
     &       '3', '0', '10', &
     &       'Maximum number of consecutive missing values allowed for'// &
     &       ' any air-temperature-measurement station; 0=unlimited', &
     &       'Maximum number of consecutive missing values allowed for'// &
     &       ' any air-temperature-measurement station; missing value set'// &
     &       ' to last valid value; 0=unlimited', &
     &       'none')/=0 ) CALL read_error(1, 'max_missing')

      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL temp_1sta_laps_restart(READ_INIT)

        ! Initialize variables, get parameter values, compute Elfac
        IF ( Temp_flag==temp_1sta_module ) THEN
          IF ( getparam(MODNAME, 'tmin_lapse', Nhru*MONTHS_PER_YEAR, 'real', Tmin_lapse)/=0 ) CALL read_error(2, 'tmin_lapse')
          IF ( getparam(MODNAME, 'tmax_lapse', Nhru*MONTHS_PER_YEAR, 'real', Tmax_lapse)/=0 ) CALL read_error(2, 'tmax_lapse')
        ELSEIF ( Temp_flag==temp_laps_module ) THEN
          IF ( getparam(MODNAME, 'hru_tlaps', Nhru, 'integer', Hru_tlaps)/=0 ) CALL read_error(2, 'hru_tlaps')
        ENDIF
        IF ( getparam(MODNAME, 'max_missing', 1, 'integer', Max_missing)/=0 ) CALL read_error(2, 'max_missing')
        Max_missing = Max_missing + 1

        Nuse_tsta = 0
        Elfac = 0.0
        IF ( Temp_flag==temp_1sta_module ) THEN
          Tcrx = 0.0
          Tcrn = 0.0
          DO i = 1, Active_hrus
            j = Hru_route_order(i)
            k = Hru_tsta(j)
            Nuse_tsta(k) = 1
            ! Hru_elev_ts is the current elevation, either hru_elev or for restart Hru_elev_ts
            Elfac(j) = (Hru_elev_ts(j)-Tsta_elev(k))/1000.0
            Tcrx(j) = Tmax_lapse(j, Start_month)*Elfac(j) - Tmax_aspect_adjust(j, Start_month)
            Tcrn(j) = Tmin_lapse(j, Start_month)*Elfac(j) - Tmin_aspect_adjust(j, Start_month)
          ENDDO
        ELSEIF ( Temp_flag==temp_laps_module ) THEN
          ierr = 0
          DO i = 1, Active_hrus
            j = Hru_route_order(i)
            CALL checkdim_param_limits(j, 'hru_tlaps', 'ntemp', Hru_tlaps(j), 1, Ntemp, ierr)
            IF ( ierr==1 ) CYCLE ! if one error found no need to compute values
            k = Hru_tsta(j)
            Nuse_tsta(k) = ACTIVE
            l = Hru_tlaps(j)
            ! Hru_elev_ts is the current glacier elevation, either hru_elev or for restart Hru_elev_ts
            CALL compute_temp_laps(Elfac(j), Hru_elev_ts(j), Tsta_elev(l), Tsta_elev(k))
          ENDDO
          IF ( ierr==1 ) THEN
            Inputerror_flag = 1
            RETURN
          ENDIF
        ENDIF

        IF ( Init_vars_from_file==0 ) THEN
          Solrad_tmax_good = Solrad_tmax
          Solrad_tmin_good = Solrad_tmin
          Tmax_cnt = 0
          Tmin_cnt = 0

          DO i = 1, Ntemp
            Tmax_prev(i) = Tmax_allrain(1, Start_month)
          ENDDO
          Tmin_prev = Tmax_prev
        ENDIF

      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL temp_1sta_laps_restart(SAVE_INIT)

      ENDIF

 9001 FORMAT ('WARNING, bad temperature, ', A, ':', F10.3, &
     &        '; temperature station: ', I0, /, 'Value set to last valid value:', F10.3)
 9002 FORMAT (/, 'ERROR, too many consecutive bad temperatures, ', A, ':', F10.3, /, &
     &        'temperature station: ', I0, /, &
     &        'Fix Data File or increase parameter max_missing' )

      END FUNCTION temp_1sta_laps

!***********************************************************************
!     Compute lapse rate for an HRU
!***********************************************************************
      SUBROUTINE compute_temp_laps(Elfac, Hru_elev, Tsta_elev_laps, Tsta_elev_base)
      USE PRMS_CONSTANTS, ONLY: NEARZERO
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Hru_elev, Tsta_elev_laps, Tsta_elev_base
      REAL, INTENT(OUT) :: Elfac
! Functions
      INTRINSIC ABS
! Local Variables
      REAL :: tdiff
!***********************************************************************
      tdiff = Tsta_elev_laps - Tsta_elev_base
      IF ( ABS(tdiff)<NEARZERO ) tdiff = 1.0
      Elfac = (Hru_elev-Tsta_elev_base)/tdiff
      END SUBROUTINE compute_temp_laps

!***********************************************************************
!     Write to or read from restart file
!***********************************************************************
      SUBROUTINE temp_1sta_laps_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_TEMP_1STA_LAPS
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=9) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Solrad_tmax_good, Solrad_tmin_good
        WRITE ( Restart_outunit ) Tmax_cnt
        WRITE ( Restart_outunit ) Tmin_cnt
        WRITE ( Restart_outunit ) Tmax_prev
        WRITE ( Restart_outunit ) Tmin_prev
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Solrad_tmax_good, Solrad_tmin_good
        READ ( Restart_inunit ) Tmax_cnt
        READ ( Restart_inunit ) Tmin_cnt
        READ ( Restart_inunit ) Tmax_prev
        READ ( Restart_inunit ) Tmin_prev
      ENDIF
      END SUBROUTINE temp_1sta_laps_restart
