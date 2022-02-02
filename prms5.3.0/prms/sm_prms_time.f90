submodule(PRMS_SET_TIME) sm_prms_time

contains
!***********************************************************************
!***********************************************************************
  integer module function prms_time()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, YEAR, MONTH, DAY, HOUR, MINUTE, MAX_DAYS_PER_YEAR, DAYS_PER_YEAR, &
     &    ACTIVE, OFF, NORTHERN, FT2_PER_ACRE, SECS_PER_HOUR, INCHES_PER_FOOT, SECS_PER_DAY, ERROR_time
      use PRMS_MMFAPI, only: dattim, deltim
      USE PRMS_MODULE, ONLY: Process_flag, Timestep, Starttime, Nowyear, Nowmonth, Nowday
      USE PRMS_BASIN, ONLY: Hemisphere, Basin_area_inv
      use PRMS_DATA_FILE, only: read_data_line
      use prms_utils, only: leap_day, julian_day, compute_julday, print_module
      IMPLICIT NONE
      ! Functions
      INTRINSIC :: SNGL
      ! Local Variables
      INTEGER :: startday
      DOUBLE PRECISION :: dt
!***********************************************************************
      prms_time = 0

      IF ( Process_flag==RUN .OR. Process_flag==INIT ) THEN

        IF ( Process_flag==RUN ) THEN
          Timestep = Timestep + 1

          CALL dattim('now', Nowtime)
          Jday = julian_day('now', 'calendar')
          Jsol = julian_day('now', 'solar')
          Julwater = julian_day('now', 'water')
          Julian_day_absolute = Julian_day_absolute + 1
          CALL read_data_line()

        ELSE ! initialize
          Modays(1) = 31
          Modays(3) = 31
          Modays(4) = 30
          Modays(5) = 31
          Modays(6) = 30
          Modays(7) = 31
          Modays(8) = 31
          Modays(9) = 30
          Modays(10) = 31
          Modays(11) = 30
          Modays(12) = 31

          Nowtime = Starttime
          Jday = julian_day('start', 'calendar')
          Jsol = julian_day('start', 'solar')
          Julwater = julian_day('start', 'water')
          startday = compute_julday(Starttime(1), Starttime(2), Starttime(3))
          Julian_day_absolute = startday
        ENDIF

        Nowyear = Nowtime(YEAR)
        Nowmonth = Nowtime(MONTH)
        Nowday = Nowtime(DAY)
        Nowhour = Nowtime(HOUR)
        Nowminute = Nowtime(MINUTE)

        IF ( leap_day(Nowyear)==1 ) THEN
          Yrdays = MAX_DAYS_PER_YEAR
          Modays(2) = 29
        ELSE
          Yrdays = DAYS_PER_YEAR
          Modays(2) = 28
        ENDIF

        ! Summer is based on equinox:
        !   Julian days 79 to 265 for Northern hemisphere
        !   Julian day 265 to 79 in Southern hemisphere
        Summer_flag = ACTIVE ! 1 = summer, 0 = winter
        IF ( Hemisphere==NORTHERN ) THEN
          IF ( Jday<79 .OR. Jday>265 ) Summer_flag = OFF ! Equinox
        ELSE ! Southern Hemisphere
          IF ( Jday>79 .AND. Jday<265 ) Summer_flag = OFF ! Equinox
        ENDIF

        dt = deltim()
        Timestep_hours = SNGL( dt )
        Timestep_days = Timestep_hours/24.0
        Timestep_minutes = Timestep_hours*60.0
        Timestep_seconds = dt*SECS_PER_HOUR
        Cfs_conv = FT2_PER_ACRE/INCHES_PER_FOOT/Timestep_seconds
        Cfs2inches = Basin_area_inv*INCHES_PER_FOOT*Timestep_seconds/FT2_PER_ACRE

        ! Check to see if in a daily or subdaily time step
        IF ( Timestep_hours>24.0 ) THEN
          PRINT *, 'ERROR, timestep > daily, fix Data File, timestep:', Timestep_hours
          ERROR STOP ERROR_time
        ELSEIF ( Timestep_hours<24.0 ) THEN
          PRINT *, 'ERROR, timestep < daily for daily model, fix Data File', Timestep_hours
          ERROR STOP ERROR_time
        ENDIF

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_prms_time)
        Timestep_seconds = SECS_PER_DAY
        Cfs_conv = FT2_PER_ACRE/INCHES_PER_FOOT/Timestep_seconds
        Cfs2inches = Basin_area_inv*INCHES_PER_FOOT*Timestep_seconds/FT2_PER_ACRE
      ENDIF

      END FUNCTION prms_time
end submodule
