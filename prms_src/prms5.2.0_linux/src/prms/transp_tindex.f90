!***********************************************************************
! Determines whether current time period is one of active transpiration
! based on a temperature index method.
!***********************************************************************
      MODULE PRMS_TRANSP_TINDEX
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, FAHRENHEIT, MONTHS_PER_YEAR, &
     &      READ_INIT, SAVE_INIT       
        USE PRMS_MODULE, ONLY: Process_flag, Nhru, Save_vars_to_file, Init_vars_from_file, &
     &      Start_month, Start_day
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Transpiration Distribution'
        character(len=13), parameter :: MODNAME = 'transp_tindex'
        character(len=*), parameter :: Version_transp = '2020-12-02'
        INTEGER, SAVE, ALLOCATABLE :: Transp_check(:)
        REAL, SAVE, ALLOCATABLE :: Tmax_sum(:), Transp_tmax_f(:)
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Transp_beg(:), Transp_end(:)
        REAL, SAVE, ALLOCATABLE :: Transp_tmax(:)
      END MODULE PRMS_TRANSP_TINDEX

      INTEGER FUNCTION transp_tindex()
      USE PRMS_TRANSP_TINDEX
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Tmaxf, Temp_units, Transp_on, Basin_transp_on 
      USE PRMS_SET_TIME, ONLY: Nowmonth, Nowday
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, getparam
      REAL, EXTERNAL :: c_to_f
      EXTERNAL :: read_error, print_module, transp_tindex_restart
! Local Variables
      INTEGER :: i, j, motmp
!***********************************************************************
      transp_tindex = 0

      IF ( Process_flag==RUN ) THEN
!******Set switch for active transpiration period
        Basin_transp_on = OFF
        DO j = 1, Active_hrus
          i = Hru_route_order(j)

!******check for month to turn check switch on or
!******transpiration switch off
          IF ( Nowday==1 ) THEN
            !******check for end of period
            IF ( Nowmonth==Transp_end(i) ) THEN
              Transp_on(i) = OFF
              Transp_check(i) = OFF
              Tmax_sum(i) = 0.0
            ENDIF
!******check for month to turn transpiration switch on or off
            IF ( Nowmonth==Transp_beg(i) ) THEN
              Transp_check(i) = ACTIVE
              Tmax_sum(i) = 0.0
            ENDIF
          ENDIF

!******If in checking period, then for each day
!******sum maximum temperature until greater than temperature index parameter,
!******at which time, turn transpiration switch on, check switch off
          ! freezing temperature assumed to be 32 degrees Fahrenheit
          IF ( Transp_check(i)==ACTIVE ) THEN
            IF ( Tmaxf(i)>32.0 ) Tmax_sum(i) = Tmax_sum(i) + Tmaxf(i)
            IF ( Tmax_sum(i)>Transp_tmax_f(i) ) THEN
              Transp_on(i) = ACTIVE
              Transp_check(i) = OFF
              Tmax_sum(i) = 0.0
            ENDIF
          ENDIF

          IF ( Basin_transp_on==OFF ) THEN
            IF ( Transp_on(i)==ACTIVE ) Basin_transp_on = ACTIVE
          ENDIF
        ENDDO

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_transp)

        ALLOCATE ( Tmax_sum(Nhru), Transp_check(Nhru), Transp_tmax_f(Nhru) )

        ALLOCATE ( Transp_beg(Nhru) )
        IF ( declparam(MODNAME, 'transp_beg', 'nhru', 'integer', &
     &       '1', '1', '12', &
     &       'Month to begin testing for transpiration', &
     &       'Month to begin summing the maximum air temperature for each HRU; when sum is greater than or'// &
     &       ' equal to transp_tmax, transpiration begins', &
     &       'month')/=0 ) CALL read_error(1, 'transp_beg')

        ALLOCATE ( Transp_end(Nhru) )
        IF ( declparam(MODNAME, 'transp_end', 'nhru', 'integer', &
     &       '13', '1', '13', &
     &       'Month to stop transpiration period', &
     &       'Month to stop transpiration computations; transpiration is computed through the end of previous month', &
     &       'month')/=0 ) CALL read_error(1, 'transp_end')

        ALLOCATE ( Transp_tmax(Nhru) )
        IF ( declparam(MODNAME, 'transp_tmax', 'nhru', 'real', &
     &       '1.0', '0.0', '1000.0', &
     &       'Tmax index to determine start of transpiration', &
     &       'Temperature index to determine the specific date of the start of the transpiration period;'// &
     &       ' the maximum air temperature for each HRU is summed starting with the first day of month transp_beg;'// &
     &       ' when the sum exceeds this index, transpiration begins', &
     &       'temp_units')/=0 ) CALL read_error(1, 'transp_tmax')

      ELSEIF ( Process_flag==INIT ) THEN

        IF ( getparam(MODNAME, 'transp_beg', Nhru, 'integer', Transp_beg)/=0 ) CALL read_error(2, 'transp_beg')
        IF ( getparam(MODNAME, 'transp_end', Nhru, 'integer', Transp_end)/=0 ) CALL read_error(2, 'transp_end')
        IF ( getparam(MODNAME, 'transp_tmax', Nhru, 'real', Transp_tmax)/=0 ) CALL read_error(2, 'transp_tmax')

        IF ( Init_vars_from_file>OFF ) CALL transp_tindex_restart(READ_INIT)
        IF ( Temp_units==FAHRENHEIT ) THEN
          Transp_tmax_f = Transp_tmax
        ELSE
          DO i = 1, Nhru
            Transp_tmax_f(i) = c_to_f(Transp_tmax(i))
          ENDDO
        ENDIF
        !DEALLOCATE ( Transp_tmax )

        IF ( Init_vars_from_file==0 ) Tmax_sum = 0.0

        motmp = Start_month + MONTHS_PER_YEAR
        Transp_check = OFF
        Basin_transp_on = OFF
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Start_month==Transp_beg(i) ) THEN
            IF ( Start_day>10 ) THEN ! rsr, why 10? if transp_tmax < 300, should be < 10
              Transp_on(i) = ACTIVE
            ELSE
              Transp_check(i) = ACTIVE
            ENDIF
          ELSEIF ( Transp_end(i)>Transp_beg(i) ) THEN
            IF ( Start_month>Transp_beg(i) .AND. Start_month<Transp_end(i) ) Transp_on(i) = ACTIVE
          ELSE
            IF ( Start_month>Transp_beg(i) .OR. motmp<Transp_end(i)+MONTHS_PER_YEAR ) Transp_on(i) = ACTIVE
          ENDIF
          IF ( Basin_transp_on==OFF ) THEN
            IF ( Transp_on(i)==ACTIVE ) Basin_transp_on = ACTIVE
          ENDIF
        ENDDO

      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL transp_tindex_restart(SAVE_INIT)

      ENDIF

      END FUNCTION transp_tindex

!***********************************************************************
!     Write to or read from restart file
!***********************************************************************
      SUBROUTINE transp_tindex_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_TRANSP_TINDEX
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL check_restart
      ! Local Variable
      CHARACTER(LEN=13) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Tmax_sum
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Tmax_sum
      ENDIF
      END SUBROUTINE transp_tindex_restart
