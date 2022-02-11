!***********************************************************************
! Compute potential solar radiation and sunlight hours for each HRU for
! each day of year; modification of soltab_prms
!
! References -- you *will* need these to figure out what is going on:
!   Swift, L.W., Jr., 1976, Algorithm for solar radiation on mountain
!   slopes: Water Resources Research, v. 12, no. 1, p. 108-112.
!
!   Lee, R., 1963, Evaluation of solar beam irradiation as a climatic parameter
!   of mountain watersheds, Colorado State University Hydrology Papers, 2,
!   50 pp.
!***********************************************************************
      MODULE PRMS_SOLTAB
      USE PRMS_CONSTANTS, ONLY: DAYS_IN_YEAR, MAX_DAYS_PER_YEAR
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Potential Solar Radiation'
      character(len=*), parameter :: MODNAME = 'soltab'
      character(len=*), parameter :: Version_soltab = '2021-08-13'
      DOUBLE PRECISION, PARAMETER :: PI=3.1415926535898D0
      DOUBLE PRECISION, PARAMETER :: RADIANS=PI/180.0D0, TWOPI=2.0D0*PI
      DOUBLE PRECISION, PARAMETER :: PI_12=12.0D0/PI
! TWOPI = 6.2831853071786
! RADIANS = 0.017453292519943
! PI_12 = 3.8197186342055
      DOUBLE PRECISION, PARAMETER :: ECCENTRICY = 0.01671D0
      ! 0.016723401  daily change -1.115E-09, eccen = 0.016723401 + (julhour-julhour(1966,1,0,18))+dmin/60)/24*-1.115E-09
      ! julday(1966,1,0.75 UT) = 2439126.25
      ! eccen = 0.01675104-0.00004180*T-0.000000126*T**2  T is julian centuries (days time from epoch, is GMT from Jan 0.0
      DOUBLE PRECISION, PARAMETER :: DEGDAY = 360.0D0/DAYS_IN_YEAR
      DOUBLE PRECISION, PARAMETER :: DEGDAYRAD = DEGDAY*RADIANS ! about 0.00143356672
! DEGDAY = 360 degrees/days in year
      DOUBLE PRECISION, SAVE :: Solar_declination(MAX_DAYS_PER_YEAR), Soltab_basinpotsw(MAX_DAYS_PER_YEAR)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_cossl(:), Soltab_sunhrs(:, :)
!   Declared Variables
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Soltab_potsw(:, :), Soltab_horad_potsw(:, :)
!   Declared Parameters
      REAL, SAVE, ALLOCATABLE :: Hru_aspect(:), Hru_slope(:)
      END MODULE PRMS_SOLTAB

!***********************************************************************
!     Main soltab routine
!***********************************************************************
      INTEGER FUNCTION soltab()
      USE PRMS_CONSTANTS, ONLY: DECL, INIT
      USE PRMS_MODULE, ONLY: Process_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: sthdecl, sthinit
!***********************************************************************
      soltab = 0

      IF ( Process_flag==DECL ) THEN
        soltab = sthdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        soltab = sthinit()
      ENDIF

      END FUNCTION soltab

!***********************************************************************
!     sthdecl - set up parameters for solar radiation computations
!   Declared Parameters
!     hru_aspect, hru_lat, hru_slope
!***********************************************************************
      INTEGER FUNCTION sthdecl()
      USE PRMS_CONSTANTS, ONLY: MAX_DAYS_PER_YEAR
      USE PRMS_MODULE, ONLY: Nhru
      USE PRMS_SOLTAB
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam !, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      sthdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_soltab)

      ALLOCATE ( Soltab_potsw(MAX_DAYS_PER_YEAR, Nhru) )
!      IF ( declvar(MODNAME, 'soltab_potsw', 'ndays,nhru', MAX_DAYS_PER_YEAR*Nhru, 'double', &
!     &     'Potential solar radiation for each Julian Day, for each HRU', &
!     &     'Langleys', Soltab_potsw)/=0 ) CALL read_error(3, 'soltab_potsw')

      ALLOCATE ( Soltab_horad_potsw(MAX_DAYS_PER_YEAR, Nhru) )
!      IF ( declvar(MODNAME, 'soltab_horad_potsw', 'ndays,nhru', MAX_DAYS_PER_YEAR*Nhru, 'double', &
!     &     'Potential solar radiation on a horizontal plane for each Julian Day, for each HRU', &
!     &     'Langleys', Soltab_horad_potsw)/=0 ) CALL read_error(3, 'soltab_horad_potsw')

      ALLOCATE ( Hru_cossl(Nhru), Soltab_sunhrs(MAX_DAYS_PER_YEAR, Nhru) )

!   Declared Parameters
      ALLOCATE ( Hru_slope(Nhru) )
      IF ( declparam(MODNAME, 'hru_slope', 'nhru', 'real', &
     &     '0.0', '0.0', '10.0', &
     &     'HRU slope', &
     &     'Slope of each HRU, specified as change in vertical length divided by change in horizontal length', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'hru_slope')

      ALLOCATE ( Hru_aspect(Nhru) )
      IF ( declparam(MODNAME, 'hru_aspect', 'nhru', 'real', &
     &     '0.0', '0.0', '360.0', &
     &     'HRU aspect', 'Aspect of each HRU', &
     &     'angular degrees')/=0 ) CALL read_error(1, 'hru_aspect')

      END FUNCTION sthdecl

!***********************************************************************
!     sthinit - Initialize soltab module - get parameter values,
!               compute soltab_potsw (potential shortwave radiation)
!               and soltab_sunhrs (hours between sunrise and sunset)
!               for each HRU for each day of the year.
!***********************************************************************
      INTEGER FUNCTION sthinit()
      USE PRMS_CONSTANTS, ONLY: MAX_DAYS_PER_YEAR, DEBUG_SOLTAB, OFF
      USE PRMS_MODULE, ONLY: Nhru, Print_debug, Glacier_flag
      USE PRMS_SOLTAB
      USE PRMS_BASIN, ONLY: Hru_type, Active_hrus, Hru_route_order, Basin_lat, Hru_lat
      IMPLICIT NONE
! Functions
      INTRINSIC :: SIN, COS, DBLE
!     INTRINSIC :: ASIN
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: compute_soltab, read_error, PRMS_open_module_file
! Local Variables
      CHARACTER(LEN=12) :: output_path
      INTEGER :: jd, j, n, file_unit, nn
      REAL :: lat
      DOUBLE PRECISION :: basin_cossl
      DOUBLE PRECISION :: basin_sunhrs(MAX_DAYS_PER_YEAR), obliquity(MAX_DAYS_PER_YEAR)
      DOUBLE PRECISION :: y, y2, y3, jddbl
!***********************************************************************
      sthinit = 0

      IF ( getparam(MODNAME, 'hru_slope', Nhru, 'real', Hru_slope)/=0 ) CALL read_error(2, 'hru_slope')
      IF ( getparam(MODNAME, 'hru_aspect', Nhru, 'real', Hru_aspect)/=0 ) CALL read_error(2, 'hru_aspect')

      DO jd = 1, MAX_DAYS_PER_YEAR
        jddbl = DBLE(jd)

        ! cosine of the solar zenith angle: http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesmBbrowser/html_code/csm_share/shr_orb_mod.F90.html
        !    jday   ! Julian cal day (1.xx to 365.xx)
        !    lat    ! Centered latitude (radians)
        !    lon    ! Centered longitude (radians)
        !    declin ! Solar declination (radians)
        ! shr_orb_cosz = sin(lat)*sin(declin) - &
   !&              cos(lat)*cos(declin)*cos(jday*2.0_SHR_KIND_R8*pi + lon)

        ! Eccentricity from equation E-2 (Dingman, S. L., 1994, Physical Hydrology. Englewood Cliffs, NJ: Prentice Hall, 575 p.)
        ! dayangle = (2*PI*(Jday-1))/365 = DEGDAYRAD*(jddbl-1.0D0) = day angle in radians
        ! eccentricity = 1.00011D0 + 0.034221D0*COS(dayangle) + 0.00128D0*SIN(dayangle) + 0.000719D0*COS(2.0D0*dayangle) + 0.000077D0*SIN(2.0D0*dayangle)
!rsr .0172 = 2PI/365 = RADIAN_YEAR = DEGDAYRAD
!rsr01/2006 commented out equations from Llowd W. Swift paper 2/1976
!       obliquity(jd) = 1.0D0 - (0.0167D0*COS((jd-3)*0.0172D0))
        obliquity(jd) = 1.0D0 - (ECCENTRICY*COS((jddbl-3.0D0)*DEGDAYRAD))
!       Solar_declination(jd) = 0.007D0 - (0.4067D0*COS((jd+10)*0.0172D0))
!       Solar_declination(jd) = ASIN(0.39785D0 * SIN( (278.9709D0+DEGDAY*jd)*RADIANS + 1.9163D0*RADIANS * SIN((356.6153D0+DEGDAY*jd)*RADIANS )) )
        ! hour = 12.0D0
!       y = DEGDAYRAD*(jddbl-1.0D0 +(hour-12.0D0)/24.0D0)
        y = DEGDAYRAD*(jddbl-1.0D0) ! assume noon
        y2 = 2.0D0*y
        y3 = 3.0D0*y
        Solar_declination(jd) = 0.006918D0 - 0.399912D0*COS(y) + 0.070257D0*SIN(y) &
     &                          - 0.006758D0*COS(y2) + 0.000907D0*SIN(y2) &
     &                          - 0.002697D0*COS(y3) + 0.00148D0*SIN(y3)
      ENDDO

!   Module Variables
      Soltab_sunhrs = 0.0D0
      Soltab_potsw = 0.0D0
      Soltab_horad_potsw = 0.0D0
      Hru_cossl = 0.0D0
      DO nn = 1, Active_hrus
        n = Hru_route_order(nn)
        CALL compute_soltab(obliquity, Solar_declination, 0.0, 0.0, Hru_lat(n), &
     &                      Hru_cossl(n), Soltab_horad_potsw(1, n), &
     &                      Soltab_sunhrs(1, n), Hru_type(n), n)
        CALL compute_soltab(obliquity, Solar_declination, Hru_slope(n), Hru_aspect(n), &
     &                      Hru_lat(n), Hru_cossl(n), Soltab_potsw(1, n), &
     &                      Soltab_sunhrs(1, n), Hru_type(n), n)
      ENDDO

      lat = SNGL( Basin_lat )
      CALL compute_soltab(obliquity, Solar_declination, 0.0, 0.0, lat, basin_cossl, &
     &                    Soltab_basinpotsw, basin_sunhrs, 0, 0)

      IF ( Print_debug==DEBUG_SOLTAB ) THEN
        output_path = 'soltab_debug'
        PRINT *, ''
        PRINT *, 'soltab debug data written to: ', output_path
        CALL PRMS_open_module_file(file_unit, output_path)
        DO n = 1, Nhru
          WRITE ( file_unit, * ) 'HRU:', n
          WRITE ( file_unit, * ) '***Soltab_sunhrs***'
          WRITE ( file_unit, '(13F8.3)' ) (Soltab_sunhrs(j,n), j=1,MAX_DAYS_PER_YEAR)
          WRITE ( file_unit, * ) '***Soltab_potsw***'
          WRITE ( file_unit, '(13F8.3)' ) (Soltab_potsw(j,n), j=1,MAX_DAYS_PER_YEAR)
        ENDDO
!       WRITE ( file_unit, * ) obliquity, Solar_declination
        WRITE ( file_unit, * ) 2.0D0/(obliquity(356)*obliquity(356)), 2.0D0/(obliquity(10)*obliquity(10)), &
     &                         2.0D0/(obliquity(23)*obliquity(23)), 2.0D0/(obliquity(38)*obliquity(38)), &
     &                         2.0D0/(obliquity(51)*obliquity(51)), 2.0D0/(obliquity(66)*obliquity(66)), &
     &                         2.0D0/(obliquity(80)*obliquity(80)), 2.0D0/(obliquity(94)*obliquity(94)), &
     &                         2.0D0/(obliquity(109)*obliquity(109)), 2.0D0/(obliquity(123)*obliquity(123)), &
     &                         2.0D0/(obliquity(138)*obliquity(138)), 2.0D0/(obliquity(152)*obliquity(152)), &
     &                         2.0D0/(obliquity(173)*obliquity(173))
        WRITE ( file_unit, * ) Solar_declination(356), Solar_declination(10), Solar_declination(23), &
     &                         Solar_declination(38), Solar_declination(51), Solar_declination(66), &
     &                         Solar_declination(80), Solar_declination(94), Solar_declination(109), &
     &                         Solar_declination(123), Solar_declination(138), Solar_declination(152), &
     &                         Solar_declination(173)
        CLOSE ( file_unit)
! from original soltab
!     data obliquity/2.06699,2.06317,2.05582,2.04520,2.03243,2.01706,2.00080,
!    +1.98553,1.96990,1.95714,1.94689,1.94005,1.93616/

!     data Solar_declination/-.410152,-.383391,-.337430,-.27198,-.190532,-.09832,0.,
!    +.09832,.190532,.27198,.33743,.383391,.410152/

!     data jday/356,10,23,38,51,66,80,94,109,123,138,152,173/
      ENDIF

      DEALLOCATE ( Hru_slope )
      IF ( Glacier_flag==OFF ) DEALLOCATE ( Hru_aspect )

      END FUNCTION sthinit

!***********************************************************************
!  compute soltab_potsw (potential shortwave radiation)
!  and soltab_sunhrs (hours between sunrise and sunset)
!  for each HRU for each day of the year.
!***********************************************************************
      SUBROUTINE compute_soltab(Obliquity, Solar_declination, Slope, Aspect, &
     &                          Latitude, Cossl, Soltab, Sunhrs, Hru_type, Id)
      USE PRMS_CONSTANTS, ONLY: MAX_DAYS_PER_YEAR, DNEARZERO
      USE PRMS_SOLTAB, ONLY: PI, TWOPI, RADIANS, PI_12
      IMPLICIT NONE
      EXTERNAL compute_t
!     Functions
      DOUBLE PRECISION, EXTERNAL :: func3
      INTRINSIC ASIN, SIN, COS, ATAN, ABS
!     Arguments
      INTEGER, INTENT(IN) :: Hru_type, Id
      DOUBLE PRECISION, INTENT(IN), DIMENSION(MAX_DAYS_PER_YEAR) :: Obliquity, Solar_declination
      REAL, INTENT(IN) :: Slope, Aspect, Latitude
      DOUBLE PRECISION, INTENT(OUT) :: Cossl
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(MAX_DAYS_PER_YEAR) :: Soltab, Sunhrs
!     Local Variables
      INTEGER :: jd
      DOUBLE PRECISION :: a, x0, x1, x2, r0, r1, d1, t, sunh, solt
      DOUBLE PRECISION :: t0, t1, t2, t3, t6, t7, d, sl
!***********************************************************************
! from SWIFT (1976)
! x0, x1, x2 = l0, l1, l2
! sl = i

      sl = ATAN(Slope)
      Cossl = COS(sl)
      a = Aspect*RADIANS

! x0 latitude of HRU
      x0 = Latitude*RADIANS

! x1 latitude of equivalent slope
! This is equation 13 from Lee, 1963
      x1 = ASIN(Cossl*SIN(x0)+SIN(sl)*COS(x0)*COS(a))

! d1 is the denominator of equation 12, Lee, 1963
      d1 = Cossl*COS(x0) - SIN(sl)*SIN(x0)*COS(a)
      IF ( ABS(d1)<DNEARZERO ) d1 = DNEARZERO

! x2 is the difference in longitude between the location of
! the HRU and the equivalent horizontal surface expressed in angle hour
! This is equation 12 from Lee, 1963
      x2 = ATAN(SIN(sl)*SIN(a)/d1)
      IF ( d1<0.0D0 ) x2 = x2 + PI

! r0 is the minute solar constant cal/cm2/min
      r0 = 2.0D0
! r0 could be 1.95 (Drummond, et al 1968)
      DO jd = 1, MAX_DAYS_PER_YEAR
        d = Solar_declination(jd)

! This is adjusted to express the variability of available insolation as
! a function of the earth-sun distance.  Lee, 1963, p 16.
! r1 is the hour solar constant cal/cm2/hour
! r0 is the minute solar constant cal/cm2/min
! 60.0D0 is minutes in an hour
! Obliquity is the obliquity of the ellipse of the earth's orbit around the sun. E
! is also called the radius vector of the sun (or earth) and is the ratio of
! the earth-sun distance on a day to the mean earth-sun distance.
! obliquity = ~23.439 (obliquity of sun)
        r1 = 60.0D0*r0/(Obliquity(jd)*Obliquity(jd))

!  compute_t is the sunrise equation.
!  t7 is the hour angle of sunset on the equivalent slope
!  t6 is the hour angle of sunrise on the equivalent slope
        CALL compute_t(x1, d, t)
        t7 = t - x2
        t6 = -t - x2

!  compute_t is the sunrise equation.
!  t1 is the hour angle of sunset on a hroizontal surface at the HRU
!  t0 is the hour angle of sunrise on a hroizontal surface at the HRU
        CALL compute_t(x0, d, t)
        t1 = t
        t0 = -t

! For HRUs that have an east or west direction component to their aspect, the
! longitude adjustment (moving the effective slope east or west) will cause either:
! (1) sunrise to be earlier than at the horizontal plane at the HRU
! (2) sunset to be later than at the horizontal plane at the HRU
! This is not possible. The if statements below check for this and adjust the
! sunrise/sunset angle hours on the equivalent slopes as necessary.
!
! t3 is the hour angle of sunrise on the slope at the HRU
! t2 is the hour angle of sunset on the slope at the HRU
        IF ( t7>t1 ) THEN
          t3 = t1
        ELSE
          t3 = t7
        ENDIF
        IF ( t6<t0 ) THEN
          t2 = t0
        ELSE
          t2 = t6
        ENDIF

        IF ( ABS(sl)<DNEARZERO ) THEN
!  solt is Swift's R4 (potential solar radiation on a sloping surface cal/cm2/day)
!  Swift, 1976, equation 6
          solt = func3(0.0D0, x0, t1, t0, r1, d)
!  sunh is the number of hours of direct sunlight (sunset minus sunrise) converted
!  from angle hours in radians to hours (24 hours in a day divided by 2 pi radians
!  in a day).
          sunh = (t1-t0)*PI_12
        ELSE
          IF ( t3<t2 ) THEN
            t2 = 0.0D0
            t3 = 0.0D0
          ENDIF
          t6 = t6 + TWOPI
          IF ( t6<t1 ) THEN
            solt = func3(x2, x1, t3, t2, r1, d) + func3(x2, x1, t1, t6, r1, d)
            sunh = (t3-t2+t1-t6)*PI_12
          ELSE
            t7 = t7 - TWOPI
            IF ( t7>t0 ) THEN
              solt = func3(x2, x1, t3, t2, r1, d) + func3(x2, x1, t7, t0, r1, d)
              sunh = (t3-t2+t7-t0)*PI_12
            ELSE
              solt = func3(x2, x1, t3, t2, r1, d)
              sunh = (t3-t2)*PI_12
            ENDIF
          ENDIF
        ENDIF
        IF ( solt<0.0D0 ) THEN
          PRINT *, 'WARNING: solar table value for day:', jd, &
     &             ' computed as:', solt, ' set to', 0.0, &
     &             ' for HRU:', Id, ' hru_type:', Hru_type
          PRINT *, 'slope, aspect, latitude, cossl', Slope, Aspect, Latitude, Cossl
          solt = 0.0D0
          PRINT *, Slope, Aspect, Latitude, Cossl, sunh
          PRINT *, t0, t1, t2, t3, t6, t7, d
        ENDIF
        IF ( sunh<DNEARZERO ) sunh = 0.0D0
        Sunhrs(jd) = sunh
        Soltab(jd) = solt

      ENDDO

      END SUBROUTINE compute_soltab

!***********************************************************************
!***********************************************************************
      SUBROUTINE compute_t(Lat, Solar_declination, T)
      USE PRMS_SOLTAB, ONLY: PI
      IMPLICIT NONE
      INTRINSIC TAN, ACOS
! Arguments
      DOUBLE PRECISION, INTENT(IN) :: Lat, Solar_declination
      DOUBLE PRECISION, INTENT(OUT) :: T
! Local Variables
      DOUBLE PRECISION :: tx
!***********************************************************************

!  This is the sunrise equation
!  Lat is the latitude
!  Solar_declination is the declination of the sun on a day
!  T is the angle hour from the local meridian (local solar noon) to the
!  sunrise (negative) or sunset (positive).  The Earth rotates at the angular
!  speed of 15 degrees/hour (2 pi / 24 hour in radians) and, therefore, T/15 degress (T*24/pi
!  in radians) gives the time of sunrise as the number of hours before the local
!  noon, or the time of sunset as the number of hours after the local noon.
!  Here the term local noon indicates the local time when the sun is exactly to
!  the south or north or exactly overhead.
      tx = -TAN(Lat)*TAN(Solar_declination)
      IF ( tx<-1.0D0 ) THEN
        T = PI
!rsr bug fix, old code would set t=acos(0.0) for tx>1 12/05
      ELSEIF ( tx>1.0D0 ) THEN
        T = 0.0D0
      ELSE
        T = ACOS(tx)
      ENDIF

      END SUBROUTINE compute_t

!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION func3(V, W, X, Y, R1, Solar_declination)
      USE PRMS_SOLTAB, ONLY: PI_12
      IMPLICIT NONE
      INTRINSIC SIN, COS
! Arguments
      DOUBLE PRECISION, INTENT(IN) :: V, W, X, Y, R1, Solar_declination
!***********************************************************************
!  This is the radian angle version of FUNC3 (eqn 6) from Swift, 1976
!  or Lee, 1963 equation 5.
!  func3 (R4) is potential solar radiation on the surface cal/cm2/day
!  V (L2) latitude angle hour offset between actual and equivalent slope
!  W (L1) latitude of the equivalent slope
!  X (T3) hour angle of sunset on equivalent slope
!  Y (T2) hour angle of sunrise on equivalent slope
!  R1 solar constant for 60 minutes
!  Solar_declination declination of sun
      func3 = R1*PI_12*(SIN(Solar_declination)*SIN(W)*(X-Y) + COS(Solar_declination)*COS(W)*(SIN(X+V)-SIN(Y+V)))

      END FUNCTION func3
