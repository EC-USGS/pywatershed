!***********************************************************************
! Distributes solar radiation to each HRU and estimates missing solar
! radiation data using a relation between solar radiation and cloud cover.
! Declared Parameters
!     ccov_slope, ccov_intcp, radj_sppt, radj_wppt, basin_solsta
!     crad_coef, crad_exp, radmax, ppt_rad_adj, rad_conv, hru_solsta
!RSR: 03/31/2008
!RSR: Warning, summer is based on equinox of Julian days 79 to 265 in
!RSR:          Northern hemisphere and Julian day 265 to 79 in Southern
!***********************************************************************
      MODULE PRMS_CCSOLRAD
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, DEBUG_less, MONTHS_PER_YEAR, ACTIVE, OFF
        USE PRMS_MODULE, ONLY: Process_flag, Print_debug, Nhru, Nsol
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Solar Radiation Distribution'
        character(len=*), parameter :: MODNAME = 'ccsolrad'
        character(len=*), parameter :: Version_ccsolrad = '2020-12-02'
        INTEGER, SAVE :: Observed_flag
        ! Declared Variables
        DOUBLE PRECISION, SAVE :: Basin_radadj, Basin_cloud_cover
        REAL, SAVE, ALLOCATABLE :: Cloud_radadj(:), Cloud_cover_hru(:)
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Crad_coef(:, :), Crad_exp(:, :)
        REAL, SAVE, ALLOCATABLE :: Ccov_slope(:, :), Ccov_intcp(:, :)
      END MODULE PRMS_CCSOLRAD
!***********************************************************************
      INTEGER FUNCTION ccsolrad()
      USE PRMS_CCSOLRAD
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area, Basin_area_inv
      USE PRMS_CLIMATEVARS, ONLY: Swrad, Basin_orad, Orad_hru, &
     &    Rad_conv, Hru_solsta, Basin_horad, Basin_potsw, Basin_swrad, Basin_solsta, Orad, Hru_ppt, &
     &    Tmax_hru, Tmin_hru, Solsta_flag, Radj_sppt, Radj_wppt, Ppt_rad_adj, Radmax
      USE PRMS_SOLTAB, ONLY: Soltab_potsw, Soltab_basinpotsw, Hru_cossl, Soltab_horad_potsw
      USE PRMS_SET_TIME, ONLY: Jday, Nowmonth, Summer_flag
      USE PRMS_OBS, ONLY: Solrad
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE, SNGL
      INTEGER, EXTERNAL :: declparam, getparam, declvar
      EXTERNAL :: read_error, print_module, print_date
! Local Variables
      INTEGER :: j, jj, k
      REAL :: pptadj, radadj, ccov
!***********************************************************************
      ccsolrad = 0

      IF ( Process_flag==RUN ) THEN
!rsr using julian day as the soltab arrays are filled by julian day
        Basin_horad = Soltab_basinpotsw(Jday)
        Basin_swrad = 0.0D0
        Basin_orad = 0.0D0
        Basin_radadj = 0.0D0
        Basin_cloud_cover = 0.0D0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)

          ! determine radiation adjustment due to precipitation
          IF ( Hru_ppt(j)>Ppt_rad_adj(j,Nowmonth) ) THEN
            IF ( Summer_flag==ACTIVE ) THEN
              pptadj = Radj_sppt(j)
            ELSE
              pptadj = Radj_wppt(j) ! Winter
            ENDIF
          ELSE
            pptadj = 1.0
          ENDIF

          ccov = Ccov_slope(j, Nowmonth)*(Tmax_hru(j)-Tmin_hru(j)) + Ccov_intcp(j, Nowmonth)
          IF ( ccov<0.0 ) THEN
            ccov = 0.0
          ELSEIF ( ccov>1.0 ) THEN
            ccov = 1.0
          ENDIF
          Cloud_cover_hru(j) = ccov
          Basin_cloud_cover = Basin_cloud_cover + DBLE( ccov*Hru_area(j) )

          radadj = Crad_coef(j, Nowmonth) + &
     &             (1.0-Crad_coef(j,Nowmonth))*((1.0-Cloud_cover_hru(j))**Crad_exp(j,Nowmonth))
          IF ( radadj>Radmax(j,Nowmonth) ) radadj = Radmax(j, Nowmonth)
          Cloud_radadj(j) = radadj*pptadj
          Basin_radadj = Basin_radadj + DBLE( Cloud_radadj(j)*Hru_area(j) )

          Orad_hru(j) = Cloud_radadj(j)*SNGL( Soltab_horad_potsw(Jday,j) )
          Basin_orad = Basin_orad + DBLE( Orad_hru(j)*Hru_area(j) )
          IF ( Solsta_flag==ACTIVE ) THEN
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
          Swrad(j) = SNGL( Soltab_potsw(Jday, j)*DBLE( Cloud_radadj(j))/Hru_cossl(j) )
          Basin_swrad = Basin_swrad + DBLE( Swrad(j)*Hru_area(j) )
        ENDDO
        Basin_orad = Basin_orad*Basin_area_inv
        Basin_radadj = Basin_radadj*Basin_area_inv
        IF ( Observed_flag==ACTIVE ) THEN
          Orad = Solrad(Basin_solsta)*Rad_conv
        ELSE
          Orad = SNGL( Basin_orad )
        ENDIF
        Basin_swrad = Basin_swrad*Basin_area_inv
        Basin_potsw = Basin_swrad
        Basin_cloud_cover = Basin_cloud_cover*Basin_area_inv

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_ccsolrad)

        ALLOCATE ( Cloud_radadj(Nhru) )
        IF ( declvar(MODNAME, 'cloud_radadj', 'nhru', Nhru, 'real', &
     &       'Radiation adjustment for cloud cover of each HRU', &
     &       'decimal fraction', Cloud_radadj)/=0 ) CALL read_error(3, 'cloud_radadj')

        IF ( declvar(MODNAME, 'basin_radadj', 'one', 1, 'double', &
     &       'Basin area-weighted average radiation adjustment for cloud cover', &
     &       'decimal fraction', Basin_radadj)/=0 ) CALL read_error(3, 'basin_radadj')

        ALLOCATE ( Cloud_cover_hru(Nhru) )
        IF ( declvar(MODNAME, 'cloud_cover_hru', 'nhru', Nhru, 'real', &
     &       'Cloud cover proportion of each HRU', &
     &       'decimal fraction', Cloud_cover_hru)/=0 ) CALL read_error(3, 'cloud_cover_hru')

        IF ( declvar(MODNAME, 'basin_cloud_cover', 'one', 1, 'double', &
     &       'Basin area-weighted average cloud cover proportion', &
     &       'decimal fraction', Basin_cloud_cover)/=0 ) CALL read_error(3, 'basin_cloud_cover')

        ! Declare Parameters
        ALLOCATE ( Crad_coef(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'crad_coef', 'nhru,nmonths', 'real', &
     &       '0.4', '0.1', '0.7', &
     &       'Coefficient in cloud cover-solar radiation relationship', &
     &       'Coefficient(B) in Thompson(1976) equation;' // &
     &       ' varies by region, contour map of values in reference', &
     &       'none')/=0 ) CALL read_error(1, 'crad_coef')
        ALLOCATE ( Crad_exp(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'crad_exp', 'nhru,nmonths', 'real', &
     &       '0.61', '0.2', '0.8', &
     &       'Exponent in cloud cover-solar radiation relationship', &
     &       'Exponent(P) in Thompson(1976) equation', &
     &       'none')/=0 ) CALL read_error(1, 'crad_exp')

        ALLOCATE ( Ccov_slope(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'ccov_slope', 'nhru,nmonths', 'real', &
     &       '-0.13', '-0.5', '-0.01', &
     &       'Slope in temperature cloud cover relationship', &
     &       'Monthly (January to December) coefficient in cloud-cover relationship', &
     &       'none')/=0 ) CALL read_error(1, 'ccov_slope')

        ALLOCATE ( Ccov_intcp(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'ccov_intcp', 'nhru,nmonths', 'real', &
     &       '1.83', '0.0', '5.0', &
     &       'Intercept in temperature cloud cover relationship', &
     &       'Monthly (January to December) intercept in cloud-cover relationship', &
     &       'none')/=0 ) CALL read_error(1, 'ccov_intcp')

      ELSEIF ( Process_flag==INIT ) THEN
! Get parameters
        IF ( getparam(MODNAME, 'crad_coef', Nhru*MONTHS_PER_YEAR, 'real', Crad_coef)/=0 ) CALL read_error(2, 'crad_coef')
        IF ( getparam(MODNAME, 'crad_exp', Nhru*MONTHS_PER_YEAR, 'real', Crad_exp)/=0 ) CALL read_error(2, 'crad_exp')
        IF ( getparam(MODNAME, 'ccov_slope', Nhru*MONTHS_PER_YEAR, 'real', Ccov_slope)/=0 ) CALL read_error(2, 'ccov_slope')
        IF ( getparam(MODNAME, 'ccov_intcp', Nhru*MONTHS_PER_YEAR, 'real', Ccov_intcp)/=0 ) CALL read_error(2, 'ccov_intcp')

        Cloud_radadj = 0.0
        Basin_radadj = 0.0D0
        Basin_cloud_cover = 0.0D0
        Cloud_cover_hru = 0.0

        Observed_flag = OFF
        IF ( Nsol>0 .AND. Basin_solsta>0 ) Observed_flag = ACTIVE

      ENDIF

      END FUNCTION ccsolrad
