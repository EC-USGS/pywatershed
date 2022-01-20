!***********************************************************************
! Distributes maximum and minimum temperatures to each HRU using a
! basin wide lapse rate applied to the temperature data, adjusted for
! distance, measured at each station
!
!     Revised 5/8/98 by Mark Mastin, J Vaccaro
!         --Declared variables basin_lapse_max and basin_lapse_min
!           They are computed in function t2dist2run
!       calculations now use all the stations and distance weight for
!       interpolating values, evens out distribution of temperature and can
!       smooth out if bad data, and still accounts for local effects
!
! Variables needed from DATA FILE: tmax, tmin
!***********************************************************************
      MODULE PRMS_TEMP_DIST2
      USE PRMS_CONSTANTS, ONLY: MONTHS_PER_YEAR, ACTIVE, DOCUMENTATION, READ_INIT, SAVE_INIT, &
     &    DNEARZERO, NEARZERO, MAXTEMP, MINTEMP, ERROR_data, GLACIER, ERROR_dim
      USE PRMS_MODULE, ONLY: Model, Nhru, Ntemp, Init_vars_from_file, Glacier_flag
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Temperature Distribution'
      character(len=10), parameter :: MODNAME = 'temp_dist2'
      character(len=*), parameter :: Version_temp = '2020-12-02'
      INTEGER, SAVE, ALLOCATABLE :: N_tsta(:), Nuse_tsta(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Dist(:, :)
      REAL, SAVE, ALLOCATABLE :: Delv(:, :), Elfac(:, :)
      REAL, SAVE :: Solrad_tmax_good, Solrad_tmin_good
!   Declared Variables
      REAL, SAVE :: Basin_lapse_max, Basin_lapse_min
!   Declared Parameters
      INTEGER, SAVE :: Max_tsta
      REAL, SAVE :: Dist_max
      REAL, SAVE :: Monmin(MONTHS_PER_YEAR), Monmax(MONTHS_PER_YEAR)
      REAL, SAVE :: Lapsemin_min(MONTHS_PER_YEAR), Lapsemin_max(MONTHS_PER_YEAR)
      REAL, SAVE :: Lapsemax_min(MONTHS_PER_YEAR), Lapsemax_max(MONTHS_PER_YEAR)
      REAL, SAVE, ALLOCATABLE :: Tsta_xlong(:), Tsta_ylat(:)
      REAL, SAVE, ALLOCATABLE :: Hru_xlong(:), Hru_ylat(:)
      END MODULE PRMS_TEMP_DIST2

!***********************************************************************
!     Main temp_dist2 routine
!***********************************************************************
      INTEGER FUNCTION temp_dist2()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, RUN, DECL, INIT, CLEAN, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Init_vars_from_file, Save_vars_to_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: t2dist2decl, t2dist2init, t2dist2run
      EXTERNAL :: temp_dist2_restart
!***********************************************************************
      temp_dist2 = 0

      IF ( Process_flag==RUN ) THEN
        temp_dist2 = t2dist2run()
      ELSEIF ( Process_flag==DECL ) THEN
        temp_dist2 = t2dist2decl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL temp_dist2_restart(READ_INIT)
        temp_dist2 = t2dist2init()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL temp_dist2_restart(SAVE_INIT)
      ENDIF

      END FUNCTION temp_dist2

!***********************************************************************
!     t2dist2decl - set up parameters for temperature computations
!   Declared Parameters
!     tsta_elev, tmax_adj, tmin_adj
!     hru_elev, hru_area, temp_units, basin_tsta, max_tsta
!     monmin, monmax, lapsemin_min, lapsemin_max, lapsemax_min
!     lapsemax_max, tsta_xlong, tsta_ylat, hru_ylat, hru_xlong, dist_max
!***********************************************************************
      INTEGER FUNCTION t2dist2decl()
      USE PRMS_TEMP_DIST2
      IMPLICIT NONE
! Functions
      INTRINSIC :: INDEX
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module, error_stop
!***********************************************************************
      t2dist2decl = 0

      CALL print_module(MODDESC, MODNAME, Version_temp)

      IF ( Ntemp<2 .AND. Model/=DOCUMENTATION ) &
     &     CALL error_stop('temp_dist2 requires at least 2 air-temperature-measurement stations', ERROR_dim)

! added by Mastin 5/8/98
      ALLOCATE ( Elfac(Nhru,Ntemp), Delv(Ntemp,Ntemp), Dist(Nhru,Ntemp), N_tsta(Nhru) )

      IF ( declvar(MODNAME, 'basin_lapse_max', 'one', 1, 'real', &
     &     'Basin area-weighted average maximum air temperature lapse rate per 1000 feet', &
     &     'degrees', Basin_lapse_max)/=0 ) CALL read_error(3, 'basin_lapse_max')

      IF ( declvar(MODNAME, 'basin_lapse_min', 'one', 1, 'real', &
     &     'Basin area-weighted average minimum air temperature lapse rate per 1000 feet', &
     &     'degrees', Basin_lapse_min)/=0 ) CALL read_error(3, 'basin_lapse_min')

      IF ( declparam(MODNAME, 'dist_max', 'one', 'real', &
     &     '1.0E9', '0.0', '1.0E9', &
     &     'Maximum distance from HRU to include a climate station', &
     &     'Maximum distance from an HRU to a measurement station for use in calcuations', &
     &     'feet')/=0 ) CALL read_error(1, 'dist_max')

      IF ( declparam(MODNAME, 'max_tsta', 'one', 'integer', &
     &     '0', 'bounded', 'ntemp', &
     &     'Maximum number of temperature stations to use for'// &
     &     ' distributing temperature to any HRU', &
     &     'Maximum number of air-temperature measurement stations to use for'// &
     &     ' distributing temperature to any HRU', &
     &     'none')/=0 ) CALL read_error(1, 'max_tsta')

! added THE FOLLOWING NEW PARAMETERS by J Vaccaro 7.98,
!       various parameters to interpolate and constrain lapse rates for temperature

      IF ( declparam(MODNAME, 'monmin', 'nmonths', 'real', &
     &     '-60.0', '-60.0', '65.0', &
     &     'Daily minimum temperature', &
     &     'Monthly minimum air temperature to constrain lowest'// &
     &     ' minimum measured air temperatures for bad values based'// &
     &     ' on historical temperature for all measurement stations', &
     &     'temp_units')/=0 ) CALL read_error(1, 'monmin')

      IF ( declparam(MODNAME, 'monmax', 'nmonths', 'real', &
     &     '100.0', '0.0', '115.0', &
     &     'Daily maximum temperature', &
     &     'Monthly maximum air temperature to constrain lowest'// &
     &     ' minimum measured air temperatures for bad values based'// &
     &     ' on historical temperature for all measurement stations', &
     &     'temp_units')/=0 ) CALL read_error(1, 'monmax')

      IF ( declparam(MODNAME, 'lapsemin_min', 'nmonths', 'real', &
     &     '-4.0', '-7.0', '-3.0', &
     &     'Monthly minimum lapse rate for minimum temperature', &
     &     'Monthly (January to December) minimum lapse rate to'// &
     &     ' constrain lowest minimum lapse rate on the basis of historical'// &
     &     ' daily air temperatures for all air-temperature measurement stations', &
     &     'temp_units/feet')/=0 ) CALL read_error(1, 'lapsemin_min')

      IF ( declparam(MODNAME, 'lapsemin_max', 'nmonths', 'real', &
     &     '3.0', '-2.0', '4.0', &
     &     'Monthly maximum lapse rate for minimum temperature', &
     &     'Monthly (January to December) minimum lapse rate to'// &
     &     ' constrain lowest maximum lapse rate on the basis of historical'// &
     &     ' daily air temperatures for all air-temperature measurement stations', &
     &     'temp_units/feet')/=0 ) CALL read_error(1, 'lapsemin_max')

      IF ( declparam(MODNAME, 'lapsemax_min', 'nmonths', 'real', &
     &     '-6.5', '-7.0', '-3.0', &
     &     'Monthly minimum lapse rate for maximum temperature', &
     &     'Monthly (January to December) maximum lapse rate to'// &
     &     ' constrain lowest minimum lapse rate on the basis of historical'// &
     &     ' daily air temperatures for all air-temperature measurement stations', &
     &     'temp_units/feet')/=0 ) CALL read_error(1, 'lapsemax_min')

      IF ( declparam(MODNAME, 'lapsemax_max', 'nmonths', 'real', &
     &     '2.0', '-3.0', '3.0', &
     &     'Monthly maximum lapse rate for maximum temperature', &
     &     'Monthly (January to December) maximum lapse rate to'// &
     &     ' constrain lowest maximum lapse rate on the basis of historical'// &
     &     ' daily air temperatures for all air-temperature measurement stations', &
     &     'temp_units/feet')/=0 ) CALL read_error(1, 'lapsemax_max')

      ALLOCATE ( Tsta_xlong(Ntemp) )
      IF ( declparam(MODNAME, 'tsta_xlong', 'ntemp', 'real', &
     &     '0.0', '-1.0E9', '1.0E9', &
     &     'Temperature station longitude, State Plane', &
     &     'Longitude of each air-temperature-measurement station,'// &
     &     ' State Plane Coordinate System', &
     &     'feet')/=0 ) CALL read_error(1, 'tsta_xlong')

      ALLOCATE ( Tsta_ylat(Ntemp) )
      IF ( declparam(MODNAME, 'tsta_ylat', 'ntemp', 'real', &
     &     '0.0', '-1.0E9', '1.0E9', &
     &     'Temperature station latitude, State Plane', &
     &     'Latitude of each air-temperature-measurement station,'// &
     &     ' State Plane Coordinate System', &
     &     'feet')/=0 ) CALL read_error(1, 'tsta_ylat')

      ALLOCATE ( Hru_ylat(Nhru) )
      IF ( declparam(MODNAME, 'hru_ylat', 'nhru', 'real', &
     &     '0.0', '-1.0E9', '1.0E9', &
     &     'HRU latitude of centroid, State Plane', &
     &     'Latitude of each HRU for the centroid, State Plane Coordinate System', &
     &     'feet')/=0 ) CALL read_error(1, 'hru_ylat')

      ALLOCATE ( Hru_xlong(Nhru) )
      IF ( declparam(MODNAME, 'hru_xlong', 'nhru', 'real', &
     &     '0.0', '-1.0E9', '1.0E9', &
     &     'HRU longitude of centroid, State Plane', &
     &     'Longitude of each HRU for the centroid, State Plane Coordinate System', &
     &     'feet')/=0 ) CALL read_error(1, 'hru_xlong')

! END NEW PARAMETERS

      END FUNCTION t2dist2decl

!***********************************************************************
!     t2dist2init - Initialize temp_dist2 module
!                 - get parameter values, compute elfac, dist
!***********************************************************************
      INTEGER FUNCTION t2dist2init()
      USE PRMS_TEMP_DIST2
      USE PRMS_BASIN, ONLY: Hru_elev
      USE PRMS_CLIMATEVARS, ONLY: Tsta_elev
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: read_error
      INTRINSIC :: DSQRT, ABS, DABS, DBLE
! Local Variables
      INTEGER :: i, j, k, n, kk, kkbig
      DOUBLE PRECISION :: distx, disty, distance, big_dist, dist2
      DOUBLE PRECISION, ALLOCATABLE :: nuse_tsta_dist(:, :)
!***********************************************************************
      t2dist2init = 0

      IF ( getparam(MODNAME, 'dist_max', 1, 'real', Dist_max)/=0 ) CALL read_error(2, 'dist_max')

      IF ( getparam(MODNAME, 'max_tsta', 1, 'real', Max_tsta)/=0 ) CALL read_error(2, 'max_tsta')
      IF ( Max_tsta==0 ) Max_tsta = Ntemp

      IF ( getparam(MODNAME, 'monmin', MONTHS_PER_YEAR, 'real', Monmin)/=0 ) CALL read_error(2, 'monmin')

      IF ( getparam(MODNAME, 'monmax', MONTHS_PER_YEAR, 'real', Monmax)/=0 ) CALL read_error(2, 'monmax')

      IF ( getparam(MODNAME, 'lapsemin_min', MONTHS_PER_YEAR, 'real', Lapsemin_min) &
     &     /=0 ) CALL read_error(2, 'lapsemin_min')

      IF ( getparam(MODNAME, 'lapsemin_max', MONTHS_PER_YEAR, 'real', Lapsemin_max) &
     &     /=0 ) CALL read_error(2, 'lapsemin_max')

      IF ( getparam(MODNAME, 'lapsemax_min', MONTHS_PER_YEAR, 'real', Lapsemax_min) &
     &     /=0 ) CALL read_error(2, 'lapsemax_min')

      IF ( getparam(MODNAME, 'lapsemax_max', MONTHS_PER_YEAR, 'real', Lapsemax_max) &
     &     /=0 ) CALL read_error(2, 'lapsemax_max')

      IF ( getparam(MODNAME, 'tsta_xlong', Ntemp, 'real', Tsta_xlong) &
     &     /=0 ) CALL read_error(2, 'tsta_xlong')

      IF ( getparam(MODNAME, 'tsta_ylat', Ntemp, 'real', Tsta_ylat) &
     &     /=0 ) CALL read_error(2, 'tsta_ylat')

      IF ( getparam(MODNAME, 'hru_xlong', Nhru, 'real', Hru_xlong) &
     &     /=0 ) CALL read_error(2, 'hru_xlong')

      IF ( getparam(MODNAME, 'hru_ylat', Nhru, 'real', Hru_ylat) &
     &     /=0 ) CALL read_error(2, 'hru_ylat')

      Basin_lapse_max = 0.0
      Basin_lapse_min = 0.0
      IF ( Init_vars_from_file==0 ) THEN
        Solrad_tmax_good = 0.0
        Solrad_tmin_good = 0.0
      ENDIF

! CALCULATE:  DISTANCE FROM EACH MRU TO EACH TEMPERATURE GAGE
!          :  ELEVATION FACTOR FOR EACH MRU TO EACH TEMPERATURE GAGE
      ALLOCATE ( Nuse_tsta(Max_tsta,Nhru), nuse_tsta_dist(Max_tsta,Nhru) )
      N_tsta = 0
      Nuse_tsta = 0
      nuse_tsta_dist = 0.0D0
      DO i = 1, Nhru
        DO k = 1, Ntemp
          Elfac(i, k) = (Hru_elev(i)-Tsta_elev(k))/1000.0
          distx = DBLE( (Hru_xlong(i)-Tsta_xlong(k))**2 )
          disty = DBLE( (Hru_ylat(i)-Tsta_ylat(k))**2 )
          distance = DSQRT(distx+disty)
          IF ( DABS(distance)<DNEARZERO ) distance = 1.0D0
          dist2 = 1.0D0/(distance/5280.0D0)
          Dist(i, k) = dist2*dist2
          IF ( distance<DBLE(Dist_max) ) THEN
            n = N_tsta(i)
            IF ( n<Max_tsta ) THEN
              n = n + 1
              Nuse_tsta(n, i) = k
              nuse_tsta_dist(n, i) = distance
              N_tsta(i) = n
            ELSE ! have max_tsta, don't use biggest distance
              big_dist = 0.0D0
              kkbig = 1
              DO kk = 1, Max_tsta
                IF ( big_dist<nuse_tsta_dist(kk,i) ) THEN
                  big_dist = nuse_tsta_dist(kk, i)
                  kkbig = kk
                ENDIF
              ENDDO
              IF ( distance<big_dist ) THEN ! if equal use first one
                Nuse_tsta(kkbig, i) = k
                nuse_tsta_dist(kkbig, i) = distance
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE ( nuse_tsta_dist )

      DO j = 1, Ntemp - 1
        DO k = j + 1, Ntemp
          Delv(j, k) = (Tsta_elev(j)-Tsta_elev(k))/1000.0
          IF ( ABS(Delv(j,k))<NEARZERO ) Delv(j, k) = 1.0
        ENDDO
      ENDDO
      ! DEALLOCATE ( Tsta_xlong, Tsta_ylat, Hru_xlong, Hru_ylat )

      END FUNCTION t2dist2init

!***********************************************************************
!     t2dist2run - Computes maximum, minumum and average temperature
!                  for each HRU based on average lapse rate for all
!                  stations. Average is constrained by maximum and
!                  minimum lapse rates for minimum and maximum
!                  temperatures (each has an upper and lower limit).
!                  Limits can be calculated using all data for the
!                  available period of record
!***********************************************************************
      INTEGER FUNCTION t2dist2run()
      USE PRMS_TEMP_DIST2
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area, Basin_area_inv, &
     &    Hru_elev_ts, Hru_type
      USE PRMS_CLIMATEVARS, ONLY: Solrad_tmax, Solrad_tmin, Basin_temp, Tmax_aspect_adjust, Tmin_aspect_adjust, &
     &    Basin_tmax, Basin_tmin, Tmaxf, Tminf, Tminc, Tmaxc, Tavgf, Tavgc, Basin_tsta, Tsta_elev
      USE PRMS_SET_TIME, ONLY: Nowmonth
      USE PRMS_OBS, ONLY: Tmax, Tmin
      IMPLICIT NONE
! Functions
      EXTERNAL :: temp_set, print_date, error_stop
      INTRINSIC :: FLOAT, DBLE, SNGL
! Local Variables
      INTEGER :: j, k, ntotx, ntotn, jj, kk, allmissing
      REAL :: tcrx, tcrn, diffn, diffx, mx, mn, tmx_sngl, tmn_sngl
      DOUBLE PRECISION :: sumtx, sumtn, tmx, tmn
      REAL :: lapsemaxmax, lapsemaxmin, lapseminmax, lapseminmin
      DOUBLE PRECISION :: sumdist
!***********************************************************************
      t2dist2run = 0

      mn = Monmin(Nowmonth)
      mx = Monmax(Nowmonth)

      Basin_tmax = 0.0D0
      Basin_tmin = 0.0D0
      Basin_temp = 0.0D0

! Calculate basin-average lapse rate using all temperature stations

      lapsemaxmax = Lapsemax_max(Nowmonth)
      lapsemaxmin = Lapsemax_min(Nowmonth)
      lapseminmax = Lapsemin_max(Nowmonth)
      lapseminmin = Lapsemin_min(Nowmonth)

      sumtx = 0.0D0
      sumtn = 0.0D0
      ntotx = 0
      ntotn = 0
      allmissing = 0
      DO j = 1, Ntemp - 1

! check for missing or bad temps based on min and max daily values
! observed for each month.

! the value of  -9999 = missing in HDB, and rdb

        IF ( Tmax(j)<mn ) CYCLE
        IF ( Tmin(j)<mn ) CYCLE
        IF ( Tmax(j)>mx ) CYCLE
        IF ( Tmin(j)>mx ) CYCLE

        DO k = j + 1, Ntemp

          IF ( Tmax(k)<mn ) CYCLE
          IF ( Tmin(k)<mn ) CYCLE
          IF ( Tmax(k)>mx ) CYCLE
          IF ( Tmin(k)>mx ) CYCLE
          allmissing = 1

          diffx = (Tmax(j)-Tmax(k))/Delv(j, k)
          diffn = (Tmin(j)-Tmin(k))/Delv(j, k)
          IF ( diffx>lapsemaxmax ) diffx = lapsemaxmax
          IF ( diffx<lapsemaxmin ) diffx = lapsemaxmin
          IF ( diffn>lapseminmax ) diffn = lapseminmax
          IF ( diffn<lapseminmin ) diffn = lapseminmin
          sumtx = sumtx + DBLE( diffx )
          ntotx = ntotx + 1
          sumtn = sumtn + DBLE( diffn )
          ntotn = ntotn + 1
        ENDDO
      ENDDO
      IF ( allmissing==0 ) THEN
        CALL print_date(1)
        CALL error_stop('all temperature stations have missing data', ERROR_data)
      ENDIF

      IF ( ntotx>0 ) THEN
        Basin_lapse_max = SNGL(sumtx)/FLOAT(ntotx)
      ELSE
        Basin_lapse_max = (lapsemaxmax+lapsemaxmin)*0.5
      ENDIF
      IF ( ntotn>0 ) THEN
        Basin_lapse_min = SNGL(sumtn)/FLOAT(ntotn)
      ELSE
        Basin_lapse_min = (lapseminmax+lapseminmin)*0.5
      ENDIF

! HRU loop for this day or timestep

      DO jj = 1, Active_hrus
        j = Hru_route_order(jj)

        tmx = 0.0D0
        tmn = 0.0D0
        sumdist = 0.0D0

        DO kk = 1, N_tsta(j)
          k = Nuse_tsta(kk, j)
          IF ( Hru_type(j)==GLACIER .AND. Glacier_flag==ACTIVE ) Elfac(j, k) = (Hru_elev_ts(j)-Tsta_elev(k))/1000.0

! check for missing or bad temps
          IF ( Tmax(k)<mn ) CYCLE
          IF ( Tmin(k)<mn ) CYCLE
          IF ( Tmax(k)>mx ) CYCLE
          IF ( Tmin(k)>mx ) CYCLE

          sumdist = sumdist + Dist(j, k)
          tcrx = Basin_lapse_max*Elfac(j, k)
          tcrn = Basin_lapse_min*Elfac(j, k)
          tmx = tmx + DBLE( (Tmax(k)+tcrx) )*Dist(j, k)
          tmn = tmn + DBLE( (Tmin(k)+tcrn) )*Dist(j, k)
        ENDDO

        IF ( sumdist>DNEARZERO ) THEN
          tmn = tmn/sumdist - DBLE( Tmin_aspect_adjust(j, Nowmonth) )
          tmx = tmx/sumdist - DBLE( Tmax_aspect_adjust(j, Nowmonth) )
        ELSE
          tmn = DBLE( (mn+mx)*0.5 )
          tmx = tmn
          PRINT *, 'WARNING, HRU:', j, '; no valid data available to set temperatures,'
          PRINT *, ' set values using monmax and monmin', tmx, tmn
          CALL print_date(1)
        ENDIF
        IF ( tmx<=tmn ) tmx = tmn + 0.01D0

        tmx_sngl = SNGL( tmx )
        tmn_sngl = SNGL( tmn )
        CALL temp_set(j, tmx_sngl, tmn_sngl, Tmaxf(j), Tminf(j), Tavgf(j), &
     &                Tmaxc(j), Tminc(j), Tavgc(j), Hru_area(j))
      ENDDO

      Basin_tmax = Basin_tmax*Basin_area_inv
      Basin_tmin = Basin_tmin*Basin_area_inv
      Basin_temp = Basin_temp*Basin_area_inv

      Solrad_tmax = Tmax(Basin_tsta)
      Solrad_tmin = Tmin(Basin_tsta)
      IF ( Solrad_tmax<MINTEMP .OR. Solrad_tmax>MAXTEMP ) THEN
        PRINT *, 'Bad temperature data to set solrad_tmax:', &
     &           Solrad_tmax, ' using last valid value:', solrad_tmax_good
        CALL print_date(0)
        Solrad_tmax = Solrad_tmax_good
      ELSE
        Solrad_tmax_good = Solrad_tmax
      ENDIF
      IF ( Solrad_tmin<MINTEMP .OR. Solrad_tmin>MAXTEMP ) THEN
        PRINT *, 'Bad temperature data to set solrad_tmin:', &
     &           Solrad_tmin, ' using last valid value:', solrad_tmin_good
        CALL print_date(0)
        Solrad_tmin = Solrad_tmin_good
      ELSE
        Solrad_tmin_good = Solrad_tmin
      ENDIF

      END FUNCTION t2dist2run

!***********************************************************************
!     temp_dist2_restart - write or read temp_dist2 restart file
!***********************************************************************
      SUBROUTINE temp_dist2_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_TEMP_DIST2
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Function
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=10) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Solrad_tmax_good, Solrad_tmin_good
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Solrad_tmax_good, Solrad_tmin_good
      ENDIF
      END SUBROUTINE temp_dist2_restart
