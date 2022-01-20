!***********************************************************************
! Computes the potential evapotranspiration using the Jensen-Haise
! formulation (Jensen and others, 1970)
!     Potet = Coef_t_mean*(Tavgf-Temp_x_mean)*Swrad/elh
!***********************************************************************
      MODULE PRMS_POTET_JH
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, MONTHS_PER_YEAR, INCH2CM
        USE PRMS_MODULE, ONLY: Process_flag, Nhru
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Potential Evapotranspiration'
        character(len=*), parameter :: MODNAME = 'potet_jh'
        character(len=*), parameter :: Version_potet = '2020-08-03'
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Jh_coef(:, :), Jh_coef_hru(:)
      END MODULE PRMS_POTET_JH

      INTEGER FUNCTION potet_jh()
      USE PRMS_POTET_JH
      USE PRMS_BASIN, ONLY: Basin_area_inv, Active_hrus, Hru_area, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Basin_potet, Potet, Tavgc, Tavgf, Swrad
      USE PRMS_SET_TIME, ONLY: Nowmonth
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module
! Local Variables
      INTEGER :: i, j
      REAL :: elh
!***********************************************************************
      potet_jh = 0

      IF ( Process_flag==RUN ) THEN
!***********************************************************************
! 597.3 cal/gm at 0 C is the energy required to change the state of 
! water to vapor
! elh is the latent heat of vaporization (not including the *2.54)
        Basin_potet = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          elh = (597.3-(0.5653*Tavgc(i)))*INCH2CM
          Potet(i) = Jh_coef(i, Nowmonth)*(Tavgf(i)-Jh_coef_hru(i))*Swrad(i)/elh
          IF ( Potet(i)<0.0 ) Potet(i) = 0.0
          Basin_potet = Basin_potet + DBLE( Potet(i)*Hru_area(i) )
        ENDDO
        Basin_potet = Basin_potet*Basin_area_inv

!******Declare parameters
      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_potet)

        ALLOCATE ( Jh_coef(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'jh_coef', 'nhru,nmonths', 'real', &
     &       '0.014', '-0.5', '1.5', &
     &       'Monthly air temperature coefficient for each HRU - Jensen-Haise', &
     &       'Monthly (January to December) air temperature coefficient used in Jensen-Haise potential ET computations'// &
     &       ' for each HRU', &
     &       'per degrees Fahrenheit')/=0 ) CALL read_error(1, 'jh_coef')

        ALLOCATE ( Jh_coef_hru(Nhru) )
        IF ( declparam(MODNAME, 'jh_coef_hru', 'nhru', 'real', &
     &       '13.0', '-99.0', '150.0', &
     &       'HRU air temperature coefficient - Jensen-Haise', &
     &       'Air temperature coefficient used in Jensen-Haise potential ET computations for each HRU', &
     &       'degrees Fahrenheit')/=0 ) CALL read_error(1, 'jh_coef_hru')

!******Get parameters
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'jh_coef', Nhru*MONTHS_PER_YEAR, 'real', Jh_coef)/=0 ) CALL read_error(2, 'jh_coef')
        IF ( getparam(MODNAME, 'jh_coef_hru', Nhru, 'real', Jh_coef_hru)/=0 ) CALL read_error(2, 'jh_coef_hru')

      ENDIF

      END FUNCTION potet_jh
