!***********************************************************************
! Computes the potential evapotranspiration using the Hamon
! formulation (Hamon, 1961)
!***********************************************************************
      MODULE PRMS_POTET_HAMON
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, MONTHS_PER_YEAR
        USE PRMS_MODULE, ONLY: Process_flag, Nhru
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Potential Evapotranspiration'
        character(len=*), parameter :: MODNAME = 'potet_hamon'
        character(len=*), parameter :: Version_potet = '2020-08-03'
        DOUBLE PRECISION, PARAMETER :: ONE_12TH = 1.0D0/12.0D0
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Hamon_coef(:, :)
      END MODULE PRMS_POTET_HAMON

      INTEGER FUNCTION potet_hamon()
      USE PRMS_POTET_HAMON
      USE PRMS_BASIN, ONLY: Basin_area_inv, Active_hrus, Hru_area, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Basin_potet, Potet, Tavgc
      USE PRMS_SOLTAB, ONLY: Soltab_sunhrs
      USE PRMS_SET_TIME, ONLY: Nowmonth, Jday
      IMPLICIT NONE
! Functions
      INTRINSIC :: EXP, DBLE, SNGL
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module
! Local Variables
      INTEGER :: i, j
      REAL :: dyl, vpsat, vdsat
!***********************************************************************
      potet_hamon = 0

      IF ( Process_flag==RUN ) THEN
!******Compute potential et for each HRU using Hamon formulation
        Basin_potet = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
! Convert daylength from hours to 12 hour multiple (equal day and night period)
          dyl = SNGL( Soltab_sunhrs(Jday, i)*ONE_12TH )
          vpsat = 6.108*EXP(17.26939*Tavgc(i)/(Tavgc(i)+237.3)) ! in Hamon 1963, eqn. used 273.3??
          vdsat = 216.7*vpsat/(Tavgc(i)+273.3)
          Potet(i) = Hamon_coef(i, Nowmonth)*dyl*dyl*vdsat !??? why day length squared??? Hamon 1963 did not square dyl, it was squared in 1961 original
          ! pet = 0.1651*dyl*vdsat*HC  (default HC would be 1.2 for units mm/day
          ! hamon_coef includes conversion to inches and 0.1651. hamon_coef = x*0.1651/25.4 = x*0.0065 so potet = inches/day
          ! Potet(i) = hamoncoef*0.1651/25.4*dyl*vdsat  1963 version
          IF ( Potet(i)<0.0 ) Potet(i) = 0.0
          Basin_potet = Basin_potet + DBLE( Potet(i)*Hru_area(i) )
        ENDDO
        Basin_potet = Basin_potet*Basin_area_inv

!******Declare parameters
      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_potet)

        ALLOCATE ( Hamon_coef(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'hamon_coef', 'nhru,nmonths', 'real', &
     &       '0.0055', '0.004', '0.008', &
     &       'Monthly air temperature coefficient - Hamon', &
     &       'Monthly (January to December) air temperature coefficient used in Hamon potential ET computations for each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'hamon_coef')

!******Get parameters
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'hamon_coef', Nhru*MONTHS_PER_YEAR, 'real', Hamon_coef)/=0 ) CALL read_error(2, 'hamon_coef')
      ENDIF

      END FUNCTION potet_hamon
