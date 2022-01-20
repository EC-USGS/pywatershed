!***********************************************************************
! Computes the potential evapotranspiration using the Hargreaves and
! Samani formulation
! Hargreaves, G.H. and Z.A. Samani, 1985. Reference crop
! evapotranspiration from temperature. Transaction of ASAE 1(2):96-99.
!***********************************************************************
      MODULE PRMS_POTET_HS
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, MONTHS_PER_YEAR
        USE PRMS_MODULE, ONLY: Process_flag, Nhru
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Potential Evapotranspiration'
        character(len=*), parameter :: MODNAME = 'potet_hs'
        character(len=*), parameter :: Version_potet = '2020-08-03'
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Hs_krs(:, :)
      END MODULE PRMS_POTET_HS

      INTEGER FUNCTION potet_hs()
      USE PRMS_POTET_HS
      USE PRMS_BASIN, ONLY: Basin_area_inv, Active_hrus, Hru_area, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Basin_potet, Potet, Tavgc, Tminc, Tmaxc, Swrad
      USE PRMS_SET_TIME, ONLY: Nowmonth
      IMPLICIT NONE
! Functions
      INTRINSIC :: SQRT, DBLE, ABS
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module
! Local Variables
      INTEGER :: i, j
      REAL :: temp_diff, swrad_inch_day !, coef_kt
!***********************************************************************
      potet_hs = 0

      IF ( Process_flag==RUN ) THEN
        Basin_potet = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          temp_diff = Tmaxc(i) - Tminc(i) ! should be mean monthlys???
!          swrad_mm_day = Swrad(i)/23.89/2.45
!          swrad_mm_day = Swrad(i)*0.04184/2.45
          swrad_inch_day = Swrad(i)*0.000673 ! Langleys->in/day
! http://www.zohrabsamani.com/research_material/files/Hargreaves-samani.pdf
!          coef_kt = 0.00185*(temp_diff**2) - 0.0433*temp_diff + 0.4023
!          Potet(i) = Hs_krs(i, Nowmonth)*coef_kt*swrad_inch_day*SQRT(temp_diff)*(Tavgc(i)+17.8)

! Danger markstro
!          Potet(i) = Hs_krs(i, Nowmonth)*swrad_inch_day*SQRT(temp_diff)*(Tavgc(i)+17.8)
          Potet(i) = Hs_krs(i, Nowmonth)*swrad_inch_day*SQRT(ABS(temp_diff))*(Tavgc(i)+17.8)
! End Danger
          IF ( Potet(i)<0.0 ) Potet(i) = 0.0
          Basin_potet = Basin_potet + DBLE( Potet(i)*Hru_area(i) )
        ENDDO
        Basin_potet = Basin_potet*Basin_area_inv

!******Declare parameters
      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_potet)

        ALLOCATE ( Hs_krs(Nhru,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'hs_krs', 'nhru,nmonths', 'real', &
     &       '0.0135', '0.01', '0.24', &
     &       'Potential ET adjustment factor - Hargreaves-Samani', &
     &       'Monthly (January to December) adjustment factor used in Hargreaves-Samani potential ET computations for each HRU', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'hs_krs')

!******Get parameters
      ELSEIF ( Process_FLAG==INIT ) THEN
        IF ( getparam(MODNAME, 'hs_krs', Nhru*MONTHS_PER_YEAR, 'real', Hs_krs)/=0 ) CALL read_error(2, 'hs_krs')
      ENDIF

      END FUNCTION potet_hs
