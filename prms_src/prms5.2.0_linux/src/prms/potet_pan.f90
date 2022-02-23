!***********************************************************************
! Computes the potential evapotranspiration for each HRU using
! pan-evaporation data
!   Declared Parameters: hru_pansta, epan_coef
!***********************************************************************
      MODULE PRMS_POTET_PAN
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, MONTHS_PER_YEAR, DEBUG_less, ERROR_dim, &
       &    READ_INIT, SAVE_INIT     
        USE PRMS_MODULE, ONLY: Process_flag, Nevap, Print_debug, Save_vars_to_file, Init_vars_from_file
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Potential Evapotranspiration'
        character(len=9), parameter :: MODNAME = 'potet_pan'
        character(len=*), parameter :: Version_potet = '2020-12-02'
        REAL, SAVE, ALLOCATABLE :: Last_pan_evap(:)
      END MODULE PRMS_POTET_PAN

      INTEGER FUNCTION potet_pan()
      USE PRMS_POTET_PAN
      USE PRMS_BASIN, ONLY: Basin_area_inv, Active_hrus, Hru_area, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Basin_potet, Potet, Hru_pansta, Epan_coef
      USE PRMS_SET_TIME, ONLY: Nowmonth
      USE PRMS_OBS, ONLY: Pan_evap
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module, potet_pan_restart, print_date, error_stop
! Local Variables
      INTEGER :: i, j, k
!***********************************************************************
      potet_pan = 0

      IF ( Process_flag==RUN ) THEN
        DO i = 1, Nevap
          IF ( Pan_evap(i)<0.0 ) THEN
            IF ( Print_debug>DEBUG_less ) THEN
              PRINT *, 'Pan_evap<0, set to last value, station:', i, '; value:', Pan_evap(i)
              CALL print_date(1)
            ENDIF
            Pan_evap(i) = Last_pan_evap(i)
          ENDIF
        ENDDO

        Basin_potet = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          k = Hru_pansta(i)
          Potet(i) = Pan_evap(k)*Epan_coef(i, Nowmonth)
          IF ( Potet(i)<0.0 ) Potet(i) = 0.0
          Basin_potet = Basin_potet + DBLE( Potet(i)*Hru_area(i) )
        ENDDO
        Basin_potet = Basin_potet*Basin_area_inv
        Last_pan_evap = Pan_evap

!******Declare parameters
      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_potet)

        IF ( Nevap==0 ) CALL error_stop('potet_pan module selected, but nevap=0', ERROR_dim)
        ALLOCATE ( Last_pan_evap(Nevap) )

      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) THEN
          CALL potet_pan_restart(READ_INIT)
        ELSE
          Last_pan_evap = 0.0
        ENDIF

      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL potet_pan_restart(SAVE_INIT)

      ENDIF

      END FUNCTION potet_pan

!***********************************************************************
!     Write to or read from restart file
!***********************************************************************
      SUBROUTINE potet_pan_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_POTET_PAN, ONLY: MODNAME, Last_pan_evap
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL check_restart
      ! Local Variable
      CHARACTER(LEN=9) :: module_name
!***********************************************************************
      IF ( In_out==0 ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Last_pan_evap
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Last_pan_evap
      ENDIF
      END SUBROUTINE potet_pan_restart
