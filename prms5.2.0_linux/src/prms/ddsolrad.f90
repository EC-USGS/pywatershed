!***********************************************************************
! Distributes solar radiation to each HRU and estimates missing solar
! radiation data using a maximum temperature per degree-day relation;
! Declared Parameters
!     dday_slope, dday_intcp, radj_sppt, radj_wppt, basin_solsta
!     radadj_slope, radadj_intcp, radmax, ppt_rad_adj, rad_conv
!     tmax_index, tmax_allrain, hru_solsta
!RSR: 03/31/2008
!RSR: Warning, summer is based on equinox of Julian days 79 to 265 in
!RSR:          Northern hemisphere and Julian day 265 to 79 in Southern
!***********************************************************************
      MODULE PRMS_DDSOLRAD
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, DEBUG_less, MONTHS_PER_YEAR, ACTIVE, OFF
        USE PRMS_MODULE, ONLY: Process_flag, Print_debug, Nhru, Nsol
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Solar Radiation Distribution'
        character(len=*), parameter :: MODNAME = 'ddsolrad'
        character(len=*), parameter :: Version_ddsolrad = '2020-12-02'
        INTEGER, SAVE :: Observed_flag
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Radadj_slope(:, :), Radadj_intcp(:, :)
        REAL, SAVE, ALLOCATABLE :: Dday_slope(:, :), Dday_intcp(:, :), Tmax_index(:, :)
      END MODULE PRMS_DDSOLRAD

      INTEGER FUNCTION ddsolrad()
      USE PRMS_DDSOLRAD
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area, Basin_area_inv
      USE PRMS_CLIMATEVARS, ONLY: Swrad, Tmax_hru, Basin_orad, Orad_hru, &
     &    Rad_conv, Hru_solsta, Basin_horad, &
     &    Basin_potsw, Basin_swrad, Basin_solsta, Orad, Hru_ppt, Tmax_allrain, &
     &    Solsta_flag, Radj_sppt, Radj_wppt, Ppt_rad_adj, Radmax
      USE PRMS_SOLTAB, ONLY: Soltab_potsw, Soltab_basinpotsw, Hru_cossl, Soltab_horad_potsw
      USE PRMS_SET_TIME, ONLY: Jday, Nowmonth, Summer_flag
      USE PRMS_OBS, ONLY: Solrad
      IMPLICIT NONE
! Functions
      INTRINSIC :: INT, FLOAT, DBLE, SNGL
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module, print_date
! Local Variables
      INTEGER :: j, jj, k, kp, kp1
      REAL :: pptadj, radadj, dday, ddayi
      REAL, SAVE, DIMENSION(26) :: solf
      DATA solf/.20, .35, .45, .51, .56, .59, .62, .64, .655, .67, .682, &
     &          .69, .70, .71, .715, .72, .722, .724, .726, .728, .73, &
     &          .734, .738, .742, .746, .75/
!***********************************************************************
      ddsolrad = 0

      IF ( Process_flag==RUN ) THEN
!rsr using julian day as the soltab arrays are filled by julian day
        Basin_horad = Soltab_basinpotsw(Jday)
        Basin_swrad = 0.0D0
        Basin_orad = 0.0D0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)

          ! set degree day and radiation adjustment limited by radmax
          dday = Dday_slope(j, Nowmonth)*Tmax_hru(j) + Dday_intcp(j, Nowmonth) + 1.0
          IF ( dday<1.0 ) dday = 1.0
          IF ( dday<26.0 ) THEN
            kp = INT(dday)
            ddayi = FLOAT(kp)
            kp1 = kp + 1
            radadj = solf(kp) + ((solf(kp1)-solf(kp))*(dday-ddayi))
            IF ( radadj>Radmax(j,Nowmonth) ) radadj = Radmax(j, Nowmonth)
          ELSE
            radadj = Radmax(j, Nowmonth)
          ENDIF

          ! Set precipitation adjument factor based on temperature
          ! and amount of precipitation
          pptadj = 1.0
          IF ( Hru_ppt(j)>Ppt_rad_adj(j,Nowmonth) ) THEN
            IF ( Tmax_hru(j)<Tmax_index(j,Nowmonth) ) THEN
              pptadj = Radj_sppt(j)
              IF ( Tmax_hru(j)>=Tmax_allrain(j,Nowmonth) ) THEN
                IF ( Summer_flag==OFF ) pptadj = Radj_wppt(j) ! Winter
              ELSE
                pptadj = Radj_wppt(j)
              ENDIF
            ELSE
              pptadj = Radadj_intcp(j, Nowmonth) + &
     &                 Radadj_slope(j, Nowmonth)*(Tmax_hru(j)-Tmax_index(j,Nowmonth))
              IF ( pptadj>1.0 ) pptadj = 1.0
            ENDIF
          ENDIF

          radadj = radadj*pptadj
          IF ( radadj<0.2 ) radadj = 0.2
          Orad_hru(j) = radadj*SNGL( Soltab_horad_potsw(Jday,j) )
          Basin_orad = Basin_orad + DBLE( Orad_hru(j)*Hru_area(j) )

          IF ( Solsta_flag==1 ) THEN
            k = Hru_solsta(j)
            IF ( k>0 ) THEN
              IF ( Solrad(k)<0.0 .OR. Solrad(k)>10000.0 ) THEN
                IF ( Print_debug>DEBUG_less ) THEN
                  PRINT *, 'WARNING, measured solar radiation missing, HRU:', j, '; station:', k, '; value computed'
                  CALL print_date(1)
                ENDIF
              ELSE
                Swrad(j) = Solrad(k)*Rad_conv
                Basin_swrad = Basin_swrad + DBLE( Swrad(j)*Hru_area(j) )
                CYCLE
              ENDIF
            ENDIF
          ENDIF
          Swrad(j) = SNGL( Soltab_potsw(Jday, j)/Soltab_horad_potsw(Jday, j)*DBLE(Orad_hru(j))/Hru_cossl(j) )
          Basin_swrad = Basin_swrad + DBLE( Swrad(j)*Hru_area(j) )
        ENDDO
        Basin_orad = Basin_orad*Basin_area_inv
        IF ( Observed_flag==1 ) THEN
          Orad = Solrad(Basin_solsta)*Rad_conv
        ELSE
          Orad = SNGL( Basin_orad )
        ENDIF
        Basin_swrad = Basin_swrad*Basin_area_inv
        Basin_potsw = Basin_swrad

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_ddsolrad)

        ! Declare Parameters
        ALLOCATE ( Dday_slope(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'dday_slope', 'nhru,nmonths', 'real', &
     &       '0.4', '0.1', '1.4', &
     &       'Slope in temperature degree-day relationship', &
     &       'Monthly (January to December) slope in degree-day equation for each HRU', &
     &       'dday/temp_units')/=0 ) CALL read_error(1, 'dday_slope')
        ALLOCATE ( Dday_intcp(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'dday_intcp', 'nhru,nmonths', 'real', &
     &       '-40.0', '-60.0', '10.0', &
     &       'Intercept in temperature degree-day relationship', &
     &       'Monthly (January to December) intercept in degree-day equation for each HRU', &
     &       'dday')/=0 ) CALL read_error(1, 'dday_intcp')
        ALLOCATE ( Radadj_slope(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'radadj_slope', 'nhru,nmonths', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Slope in air temperature range adjustment to degree-day equation', &
     &       'Monthly (January to December) slope in air temperature range adjustment to degree-day equation for each HRU', &
     &       'dday/temp_units')/=0 ) CALL read_error(1, 'radadj_slope')
        ALLOCATE ( Radadj_intcp(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'radadj_intcp', 'nhru,nmonths', 'real', &
     &       '1.0', '0.0', '1.0', &
     &       'Intercept in air temperature range adjustment to degree-day equation', &
     &       'Monthly (January to December) intercept in air temperature range adjustment to degree-day equation for each HRU', &
     &       'dday')/=0 ) CALL read_error(1, 'radadj_intcp')
        ALLOCATE ( Tmax_index(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'tmax_index', 'nhru,nmonths', 'real', &
     &       '50.0', '-10.0', '110.0', &
     &       'Monthly index temperature', &
     &       'Monthly (January to December) index temperature used'// &
     &       ' to determine precipitation adjustments to solar radiation for each HRU', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmax_index')

      ELSEIF ( Process_flag==INIT ) THEN
! Get parameters
        IF ( getparam(MODNAME, 'dday_slope', Nhru*MONTHS_PER_YEAR, 'real', Dday_slope)/=0 ) CALL read_error(2, 'dday_slope')
        IF ( getparam(MODNAME, 'dday_intcp', Nhru*MONTHS_PER_YEAR, 'real', Dday_intcp)/=0 ) CALL read_error(2, 'dday_intcp')
        IF ( getparam(MODNAME, 'radadj_slope', Nhru*MONTHS_PER_YEAR, 'real', Radadj_slope)/=0 ) CALL read_error(2, 'radadj_slope')
        IF ( getparam(MODNAME, 'radadj_intcp', Nhru*MONTHS_PER_YEAR, 'real', Radadj_intcp)/=0 ) CALL read_error(2, 'radadj_intcp')
        IF ( getparam(MODNAME, 'tmax_index', Nhru*MONTHS_PER_YEAR, 'real', Tmax_index)/=0 ) CALL read_error(2, 'tmax_index')
        Observed_flag = OFF
        IF ( Nsol>0 .AND. Basin_solsta>0 ) Observed_flag = ACTIVE

      ENDIF

      END FUNCTION ddsolrad
