!***********************************************************************
! This is a preprocess module.
! Determine the latest "killing" frost date in the spring and the
! earliest date in the fall.
! Declared Parameters: frost_temp
!***********************************************************************
      INTEGER FUNCTION frost_date()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, DAYS_PER_YEAR, NORTHERN
      USE PRMS_MODULE, ONLY: Process_flag, Nhru
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area, Basin_area_inv, Hemisphere
      USE PRMS_CLIMATEVARS, ONLY: Tmin_hru
      USE PRMS_SET_TIME, ONLY: Jsol
      IMPLICIT NONE
      character(len=*), parameter :: MODDESC = 'Preprocessing'
      character(len=*), parameter :: MODNAME = 'frost_date'
      character(len=*), parameter :: Version_frost_date = '2020-12-02'
! Functions
      INTRINSIC :: NINT, DBLE
      INTEGER, EXTERNAL :: declparam, getparam, get_season
      EXTERNAL :: read_error, write_integer_param, PRMS_open_module_file, print_module
! Declared Parameters
      REAL, SAVE, ALLOCATABLE :: Frost_temp(:)
! Local Variables
      ! fall_frost: The number of solar days after winter solstice of
      !             the first killing frost of the fall
      ! spring_frost: The number of solar days after winter solstice of
      !               the last killing frost of the spring
      INTEGER, SAVE, ALLOCATABLE :: fall_frost(:), spring_frost(:)
      ! basin_fall_frost: The basin average solar date of the first
      !                   killing frost of the fall
      ! basin_spring_frost: The basin average solar date of the last
      !                     killing frost of the fall
      DOUBLE PRECISION, SAVE :: basin_fall_frost, basin_spring_frost
      INTEGER, SAVE :: oldSeason, fallFrostCount, springFrostCount, fall1, spring1
      INTEGER, SAVE :: switchToSpringToday, switchToFallToday, Iunit
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: fallFrostSum(:), springFrostSum(:)
      INTEGER, SAVE, ALLOCATABLE :: currentFallFrost(:), currentSpringFrost(:)
      INTEGER :: season, j, jj, basin_fall(1), basin_spring(1)
!***********************************************************************
      frost_date = 0

      IF ( Process_flag==RUN ) THEN
        season = get_season()

! Figure out if the season changes on this timestep. Putting this
! check here makes the blocks below easier to understand.
        IF ( oldSeason/=season ) THEN
          IF ( season==1 ) THEN
            switchToSpringToday = ACTIVE
          ELSE
            switchToFallToday = ACTIVE
          ENDIF
        ELSE
          switchToSpringToday = OFF
          switchToFallToday = OFF
        ENDIF
        oldSeason = season

! If this is the first timestep of fall, unset the CurrentFallFrost
! variable. Also since we are finished looking for spring frosts,
! add the CurrentSpringFrost dates to the spring_frost variable
! (average date of the spring frost for each HRU).
        IF ( switchToFallToday==ACTIVE ) THEN
          fallFrostCount = fallFrostCount + 1
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            currentFallFrost(j) = 0
            IF ( currentSpringFrost(j)==0 ) currentSpringFrost(j) = spring1
            springFrostSum(j) = springFrostSum(j) + DBLE( currentSpringFrost(j) )
          ENDDO

! If this is the first timestep of spring, unset the
! CurrentSpringFrost variable. Also since we are finished looking
! for fall frosts, add the CurrentFallFrost dates to the fall_frost
! variable (average date of the fall frost for each HRU).
        ELSEIF ( switchToSpringToday==ACTIVE ) THEN
          springFrostCount = springFrostCount + 1
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            currentSpringFrost(j) = 0
            IF ( currentFallFrost(j)==0 ) currentFallFrost(j) = fall1
            fallFrostSum(j) = fallFrostSum(j) + DBLE( currentFallFrost(j) )
          ENDDO
        ENDIF

! If the season is fall, look for the earliest fall frost.
! should probably look for x number of consecutive days or frosts before setting killing frost
        IF ( season==2 ) THEN
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            IF ( Tmin_hru(j)<=Frost_temp(j) .AND. currentFallFrost(j)==0 ) currentFallFrost(j) = Jsol
          ENDDO
! This is the spring phase, look for the latest spring frost.
        ELSE
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            IF ( Tmin_hru(j)<=Frost_temp(j) ) currentSpringFrost(j) = Jsol
          ENDDO
        ENDIF

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_frost_date)

        ALLOCATE ( Frost_temp(Nhru) )
        IF ( declparam(MODNAME, 'frost_temp', 'nhru', 'real', &
     &       '28.0', '-10.0', '32.0', &
     &       'Temperature of killing frost', 'Temperature of killing frost', &
     &       'temp_units')/=0 ) CALL read_error(1, 'frost_temp')

! Allocate arrays for local variables
        ALLOCATE ( fall_frost(Nhru), spring_frost(Nhru) )
        ALLOCATE ( fallFrostSum(Nhru), springFrostSum(Nhru) )
        ALLOCATE ( currentFallFrost(Nhru), currentSpringFrost(Nhru) )

      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'frost_temp', Nhru, 'real', Frost_temp)/=0 ) CALL read_error(2, 'frost_temp')
        fall_frost = 0
        spring_frost = 0
        currentFallFrost = 0
        currentSpringFrost = 0
        fallFrostSum = 0.0D0
        springFrostSum = 0.0D0
        fallFrostCount = 0
        springFrostCount = 0
        CALL PRMS_open_module_file(Iunit, 'frost_date.param')
        oldSeason = get_season()
        IF ( Hemisphere==NORTHERN ) THEN
          spring1 = 1
          fall1 = DAYS_PER_YEAR
        ELSE
          spring1 = DAYS_PER_YEAR
          fall1 = 1
        ENDIF

      ELSEIF ( Process_flag==CLEAN ) THEN
        basin_fall_frost = 0.0D0
        basin_spring_frost = 0.0D0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)
          IF ( fallFrostCount==0 ) fallFrostCount = 1
          fall_frost(j) = NINT( fallFrostSum(j)/DBLE( fallFrostCount ) )
          IF ( fallFrostCount==0 ) fallFrostCount = 1
          spring_frost(j) = NINT( springFrostSum(j)/DBLE( springFrostCount ) )
          fall_frost(j) = fall_frost(j) + 10
          IF ( fall_frost(j)>DAYS_PER_YEAR ) fall_frost(j) = DAYS_PER_YEAR
          spring_frost(j) = spring_frost(j) + 10
          IF ( spring_frost(j)>DAYS_PER_YEAR ) spring_frost(j) = spring_frost(j) - DAYS_PER_YEAR
          basin_fall_frost = basin_fall_frost + fall_frost(j)*Hru_area(j)
          basin_spring_frost = basin_spring_frost + spring_frost(j)*Hru_area(j)
        ENDDO
        basin_fall_frost = basin_fall_frost*Basin_area_inv
        basin_spring_frost = basin_spring_frost*Basin_area_inv

        CALL write_integer_param(Iunit, 'fall_frost', 'nhru', Nhru, fall_frost)
        CALL write_integer_param(Iunit, 'spring_frost', 'nhru', Nhru, spring_frost)
        basin_fall(1) = NINT( basin_fall_frost )
        basin_spring(1) = NINT( basin_spring_frost )
        CALL write_integer_param(Iunit, 'basin_fall_frost', 'one', 1, basin_fall)
        CALL write_integer_param(Iunit, 'basin_spring_frost', 'one', 1, basin_spring)
      ENDIF

      END FUNCTION frost_date

!*************************************************************
! Figure out if the current solar day is in "spring" or "fall"
!*************************************************************
      INTEGER FUNCTION get_season()
      USE PRMS_CONSTANTS, ONLY: NORTHERN
      USE PRMS_BASIN, ONLY: Hemisphere
      USE PRMS_SET_TIME, ONLY: Jsol
!*************************************************************
      get_season = 2 ! default is fall frost
      IF ( Hemisphere==NORTHERN ) THEN
        IF ( Jsol>0 .AND. Jsol<183 ) get_season = 1 ! This is the spring phase
      ELSE ! Southern Hemisphere
        IF ( Jsol>182 .AND. Jsol<367 ) get_season = 1 ! This is the spring phase
      ENDIF
      END FUNCTION get_season
