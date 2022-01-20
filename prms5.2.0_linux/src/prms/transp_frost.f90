!***********************************************************************
! Determine whether transpiration is occurring. Transpiration is based
! on time between the last spring and the first fall killing frost.
!***********************************************************************
      MODULE PRMS_TRANSP_FROST
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, ACTIVE, OFF
        USE PRMS_MODULE, ONLY: Process_flag, Nhru
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Transpiration Distribution'
        character(len=*), parameter :: MODNAME = 'transp_frost'
        character(len=*), parameter :: Version_transp = '2020-12-02'
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Fall_frost(:), Spring_frost(:)
      END MODULE PRMS_TRANSP_FROST

      INTEGER FUNCTION transp_frost()
      USE PRMS_TRANSP_FROST
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Transp_on, Basin_transp_on
      USE PRMS_SET_TIME, ONLY: Jsol
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, getparam
      EXTERNAL :: read_error, print_module
! Local Variables
      INTEGER :: i, j
!***********************************************************************
      transp_frost = 0

      IF ( Process_flag==RUN ) THEN
!******Set switch for active transpiration period
! If the current solar day is between the last frost of the
! spring and the first frost of the fall, then transpiration
! is on for the HRU. If any HRU is transpiring, then
! Basin_transp_on is set to 1 (ACTIVE).
        Basin_transp_on = OFF
        Transp_on = OFF
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Jsol>=Spring_frost(i) .AND. Jsol<=Fall_frost(i) ) THEN
            Transp_on(i) = ACTIVE
            Basin_transp_on = ACTIVE
          ENDIF
        ENDDO

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_transp)

        ALLOCATE ( Spring_frost(Nhru) )
        IF ( declparam(MODNAME, 'spring_frost', 'nhru', 'integer', &
     &       '111', '1', '366', &
     &       'The solar date (number of days after winter solstice) of the last killing frost of the spring', &
     &       'The solar date (number of days after winter solstice) of the last killing frost of the spring', &
     &       'Solar date')/=0 ) CALL read_error(1, 'spring_frost')

        ALLOCATE ( Fall_frost(Nhru) )
        IF ( declparam(MODNAME, 'fall_frost', 'nhru', 'integer', &
     &       '264', '1', '366', &
     &       'The solar date (number of days after winter solstice) of the first killing frost of the fall', &
     &       'The solar date (number of days after winter solstice) of the first killing frost of the fall', &
     &       'Solar date')/=0 ) CALL read_error(1, 'fall_frost')

      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'spring_frost', Nhru, 'integer', Spring_frost)/=0 ) CALL read_error(2, 'spring_frost')
        IF ( getparam(MODNAME, 'fall_frost', Nhru, 'integer', Fall_frost)/=0 ) CALL read_error(2, 'fall_frost')

        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Jsol>=Spring_frost(i) .AND. Jsol<=Fall_frost(i) ) THEN
            Transp_on(i) = ACTIVE
            Basin_transp_on = ACTIVE
          ENDIF
        ENDDO
      ENDIF

      END FUNCTION transp_frost
 