!***********************************************************************
! Read and makes available climate data (tmin, tmax, precip, potential
! solar radiation, potential evapotranspieration) and/or transpiration
! on, by HRU from files pre-processed Data Files available for other
! PRMS modules
!***********************************************************************
      MODULE PRMS_CLIMATE_HRU
        USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Climate Input'
        character(len=*), parameter :: MODNAME = 'climate_hru'
        character(len=*), parameter :: Version_climate_hru = '2021-11-19'
        INTEGER, SAVE :: Precip_unit, Tmax_unit, Tmin_unit, Et_unit, Swrad_unit, Transp_unit
        INTEGER, SAVE :: Humidity_unit, Windspeed_unit, AET_unit, PET_unit, Irrigated_area_unit
        INTEGER, SAVE :: Albedo_unit, Cloud_cover_unit
        ! Control Parameters
        CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: Tmin_day, Tmax_day, Precip_day, Potet_day, Swrad_day, Transp_day
        CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: Humidity_day, Windspeed_day, AET_cbh_file, PET_cbh_file, irrigated_area_cbh_file
        CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: Albedo_day, Cloud_cover_day
        INTEGER, SAVE :: Cbh_check_flag, Cbh_binary_flag
        ! Declared Variables
        DOUBLE PRECISION, SAVE :: Basin_windspeed, Basin_aet_external, Basin_pet_external, Basin_irrigated_area
        REAL, ALLOCATABLE :: Humidity_hru(:), Windspeed_hru(:), AET_external(:), PET_external(:), Irrigated_area(:)
        REAL, ALLOCATABLE :: Albedo_hru(:), Cloud_cover_cbh(:)
        ! Declared Parameters
        REAL, SAVE, ALLOCATABLE :: Rain_cbh_adj(:, :), Snow_cbh_adj(:, :), Potet_cbh_adj(:, :)
        REAL, SAVE, ALLOCATABLE :: Tmax_cbh_adj(:, :), Tmin_cbh_adj(:, :)
      END MODULE PRMS_CLIMATE_HRU

      INTEGER FUNCTION climate_hru()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, RUN, DECL, INIT, DOCUMENTATION, &
     &    MM2INCH, MINTEMP, MAXTEMP, ERROR_cbh, CELSIUS, MONTHS_PER_YEAR
      use PRMS_MMFAPI, only: declvar_dble, declvar_real
      use PRMS_READ_PARAM_FILE, only: declparam, getparam_real
      use PRMS_CONTROL_FILE, only: control_integer, control_string
      USE PRMS_MODULE, ONLY: Process_flag, Model, Nhru, Climate_transp_flag, Orad_flag, &
     &    Climate_precip_flag, Climate_temp_flag, Climate_potet_flag, Climate_swrad_flag, &
     &    Start_year, Start_month, Start_day, Humidity_cbh_flag, Windspeed_cbh_flag, &
     &    irrigated_area_flag, AET_cbh_flag, PET_cbh_flag, Albedo_cbh_flag, Cloud_cover_cbh_flag, Nowmonth
      USE PRMS_CLIMATE_HRU
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area, Basin_area_inv
      USE PRMS_CLIMATEVARS, ONLY: Solrad_tmax, Solrad_tmin, Basin_temp, &
     &    Basin_tmax, Basin_tmin, Tmaxf, Tminf, Tminc, Tmaxc, Tavgf, &
     &    Tavgc, Hru_ppt, Hru_rain, Hru_snow, Prmx, Pptmix, Newsnow, &
     &    Precip_units, Tmax_allrain_f, Adjmix_rain, &
     &    Basin_ppt, Basin_potet, Potet, Basin_snow, Basin_rain, &
     &    Basin_horad, Orad, Swrad, Basin_potsw, Basin_swrad, Basin_obs_ppt, &
     &    Transp_on, Basin_transp_on, Tmax_allsnow_f, Basin_humidity, Ppt_zero_thresh
      USE PRMS_SET_TIME, ONLY: Jday
      USE PRMS_SOLTAB, ONLY: Soltab_basinpotsw, Hru_cossl, Soltab_potsw
      use prms_utils, only: find_current_time, find_header_end, print_date, print_module, read_error
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE, SNGL
      EXTERNAL :: precip_form, read_cbh_date, check_cbh_value, check_cbh_intvalue
! Local Variables
      INTEGER :: yr, mo, dy, i, hr, mn, sec, jj, ierr, istop, missing, ios
      DOUBLE PRECISION :: sum_obs
      REAL :: tmax_hru, tmin_hru, ppt, harea
!***********************************************************************
      climate_hru = 0
      ierr = 0
      IF ( Process_flag==RUN ) THEN
        IF ( Climate_temp_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Tmax_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tmaxf(i), i=1,Nhru)
          ELSE
            READ ( Tmax_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tmaxf(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'tmaxf', ios, ierr)
          ENDIF
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Tmin_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tminf(i), i=1,Nhru)
          ELSE
            READ ( Tmin_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tminf(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'tminf', ios, ierr)
          ENDIF
          Basin_tmax = 0.0D0
          Basin_tmin = 0.0D0
          Basin_temp = 0.0D0
        ENDIF

        IF ( Climate_precip_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Precip_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Hru_ppt(i), i=1,Nhru)
          ELSE
            READ ( Precip_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Hru_ppt(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSE
            IF ( Ppt_zero_thresh>0.0 ) THEN
              DO i = 1, Nhru
                IF ( Hru_ppt(i)<Ppt_zero_thresh ) Hru_ppt(i) = 0.0
              ENDDO
            ENDIF
            IF ( Cbh_check_flag==ACTIVE ) CALL read_cbh_date(yr, mo, dy, 'Hru_ppt', ios, ierr)
          ENDIF
          Basin_ppt = 0.0D0
          Basin_rain = 0.0D0
          Basin_snow = 0.0D0
          sum_obs = 0.0D0
        ENDIF

        IF ( Climate_potet_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Et_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Potet(i), i=1,Nhru)
          ELSE
            READ ( Et_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Potet(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'potet', ios, ierr)
          ENDIF
          Basin_potet = 0.0D0
        ENDIF

        IF ( AET_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( AET_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (AET_external(i), i=1,Nhru)
          ELSE
            READ ( AET_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (AET_external(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'AET_external', ios, ierr)
          ENDIF
          Basin_aet_external = 0.0D0
        ENDIF

        IF ( PET_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( PET_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (PET_external(i), i=1,Nhru)
          ELSE
            READ ( PET_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (PET_external(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'PET_external', ios, ierr)
          ENDIF
          Basin_pet_external = 0.0D0
        ENDIF

        IF ( irrigated_area_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Irrigated_area_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Irrigated_area(i), i=1,Nhru)
          ELSE
            READ ( Irrigated_area_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Irrigated_area(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'irrigated_area', ios, ierr)
          ENDIF
          Basin_irrigated_area = 0.0D0
        ENDIF

        IF ( Climate_swrad_flag==ACTIVE ) THEN
          IF ( Orad_flag==OFF ) THEN
            IF ( Cbh_binary_flag==OFF ) THEN
              READ ( Swrad_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Swrad(i), i=1,Nhru)
            ELSE
              READ ( Swrad_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Swrad(i), i=1,Nhru)
            ENDIF
          ELSE
            IF ( Cbh_binary_flag==OFF ) THEN
              READ ( Swrad_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Swrad(i), i=1,Nhru), Orad
            ELSE
              READ ( Swrad_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Swrad(i), i=1,Nhru), Orad
            ENDIF
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'swrad', ios, ierr)
          ENDIF
          Basin_swrad = 0.0D0
        ENDIF

        IF ( Climate_transp_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Transp_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Transp_on(i), i=1,Nhru)
          ELSE
            READ ( Transp_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Transp_on(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'transp_on', ios, ierr)
          ENDIF
          Basin_transp_on = OFF
        ENDIF

        IF ( Humidity_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Humidity_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Humidity_hru(i), i=1,Nhru)
          ELSE
            READ ( Humidity_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Humidity_hru(i), i=1,Nhru)
          ENDIF
          IF ( Cbh_check_flag==ACTIVE ) CALL read_cbh_date(yr, mo, dy, 'humidity_hru', ios, ierr)
          Basin_humidity = 0.0D0
        ENDIF

        IF ( Albedo_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Albedo_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Albedo_hru(i), i=1,Nhru)
          ELSE
            READ ( Albedo_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Albedo_hru(i), i=1,Nhru)
          ENDIF
          IF ( Cbh_check_flag==ACTIVE ) CALL read_cbh_date(yr, mo, dy, 'albedo_hru', ios, ierr)
        ENDIF

        IF ( Cloud_cover_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Cloud_cover_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Cloud_cover_cbh(i), i=1,Nhru)
          ELSE
            READ ( Cloud_cover_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Cloud_cover_cbh(i), i=1,Nhru)
          ENDIF
          IF ( Cbh_check_flag==ACTIVE ) CALL read_cbh_date(yr, mo, dy, 'cloud_cover_cbh', ios, ierr)
        ENDIF

        IF ( Windspeed_cbh_flag==ACTIVE ) THEN
          IF ( Cbh_binary_flag==OFF ) THEN
            READ ( Windspeed_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Windspeed_hru(i), i=1,Nhru)
          ELSE
            READ ( Windspeed_unit, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Windspeed_hru(i), i=1,Nhru)
          ENDIF
          IF ( ios/=0 ) THEN
            ierr = 1
          ELSEIF ( Cbh_check_flag==ACTIVE ) THEN
            CALL read_cbh_date(yr, mo, dy, 'windspeed_hru', ios, ierr)
          ENDIF
          Basin_windspeed = 0.0D0
        ENDIF

        IF ( ierr/=0 ) ERROR STOP ERROR_cbh

        IF ( Cbh_check_flag==ACTIVE ) THEN
          missing = 0
          IF ( Climate_temp_flag==ACTIVE ) THEN
            CALL check_cbh_value('tmaxf', Tmaxf, MINTEMP, MAXTEMP, missing)
            CALL check_cbh_value('tminf', Tminf, MINTEMP, MAXTEMP, missing)
          ENDIF
          IF ( Climate_potet_flag==ACTIVE ) CALL check_cbh_value('potet', Potet, 0.0, 50.0, missing)
          IF ( Climate_swrad_flag==ACTIVE ) CALL check_cbh_value('swrad', Swrad, 0.0, 1000.0, missing)
          IF ( Climate_transp_flag==ACTIVE ) CALL check_cbh_intvalue('transp_on', Transp_on, 0, 1, missing)
          IF ( Climate_precip_flag==ACTIVE ) CALL check_cbh_value('hru_ppt', Hru_ppt, 0.0, 300.0, missing)
          IF ( Humidity_cbh_flag==ACTIVE ) CALL check_cbh_value('humidity_hru', Humidity_hru, 0.0, 100.0, missing)
          IF ( Albedo_cbh_flag==ACTIVE ) CALL check_cbh_value('albedo_hru', Albedo_hru, 0.0, 100.0, missing)
          IF ( Cloud_cover_cbh_flag==ACTIVE ) CALL check_cbh_value('cloud_cover_cbh', Cloud_cover_cbh, 0.0, 100.0, missing)
          IF ( Windspeed_cbh_flag==ACTIVE ) CALL check_cbh_value('windspeed_hru', Windspeed_hru, 0.0, 400.0, missing)
          IF ( AET_cbh_flag==ACTIVE ) CALL check_cbh_value('AET_external', AET_external, 0.0, 50.0, missing)
          IF ( PET_cbh_flag==ACTIVE ) CALL check_cbh_value('PET_external', PET_external, 0.0, 50.0, missing)
          IF ( irrigated_area_flag==ACTIVE ) CALL check_cbh_value('irrigated_area', Irrigated_area, 0.0, 9999.0, missing)

          IF ( missing==1 ) THEN
            CALL print_date(0)
            ERROR STOP ERROR_cbh
          ENDIF
        ENDIF

        IF ( Climate_precip_flag==ACTIVE ) THEN
!******Initialize HRU variables
          Pptmix = OFF
          Newsnow = OFF
          Prmx = 0.0
          Hru_rain = 0.0
          Hru_snow = 0.0
        ENDIF

        DO jj = 1, Active_hrus
          i = Hru_route_order(jj)
          harea = Hru_area(i)

          IF ( Climate_temp_flag==ACTIVE ) THEN
            tmax_hru = Tmaxf(i) + Tmax_cbh_adj(i, Nowmonth)
            tmin_hru = Tminf(i) + Tmin_cbh_adj(i, Nowmonth)
            CALL temp_set(i, tmax_hru, tmin_hru, Tmaxf(i), Tminf(i), &
     &                    Tavgf(i), Tmaxc(i), Tminc(i), Tavgc(i), harea)
          ENDIF

          IF ( Climate_potet_flag==ACTIVE ) THEN
            Potet(i) = Potet(i)*Potet_cbh_adj(i, Nowmonth)
            Basin_potet = Basin_potet + DBLE( Potet(i)*harea )
          ENDIF

          IF ( Climate_precip_flag==ACTIVE ) THEN
            IF ( Hru_ppt(i)>0.0 ) THEN
              IF ( Precip_units==CELSIUS ) Hru_ppt(i) = Hru_ppt(i)*MM2INCH
              ppt = Hru_ppt(i)
              CALL precip_form(ppt, Hru_ppt(i), Hru_rain(i), Hru_snow(i), &
     &                         Tmaxf(i), Tminf(i), Pptmix(i), Newsnow(i), &
     &                         Prmx(i), Tmax_allrain_f(i,Nowmonth), &
     &                         Rain_cbh_adj(i,Nowmonth), Snow_cbh_adj(i,Nowmonth), &
     &                         Adjmix_rain(i,Nowmonth), harea, sum_obs, Tmax_allsnow_f(i,Nowmonth))
            ELSEIF ( Hru_ppt(i)<0.0 ) THEN
              PRINT *, 'ERROR, negative precipitation value entered in CBH File, HRU:', i
              ierr = 1
!              Hru_ppt(i) = 0.0
            ENDIF
          ENDIF
          IF ( Climate_transp_flag==ACTIVE ) THEN
            IF ( Transp_on(i)==ACTIVE ) Basin_transp_on = ACTIVE
          ENDIF
          IF ( Climate_swrad_flag==ACTIVE ) Basin_swrad = Basin_swrad + DBLE( Swrad(i)*harea )
          IF ( Humidity_cbh_flag==ACTIVE ) Basin_humidity = Basin_humidity + DBLE( Humidity_hru(i)*harea )
          IF ( Windspeed_cbh_flag==ACTIVE ) Basin_windspeed = Basin_windspeed + DBLE( Windspeed_hru(i)*harea )
          IF ( AET_cbh_flag==ACTIVE ) Basin_aet_external = Basin_aet_external + DBLE( AET_external(i)*harea )
          IF ( PET_cbh_flag==ACTIVE ) Basin_pet_external = Basin_pet_external + DBLE( PET_external(i)*harea )
          IF ( irrigated_area_flag==ACTIVE ) Basin_irrigated_area = Basin_irrigated_area + DBLE( Irrigated_area(i)*harea )
        ENDDO

        IF ( ierr==1 ) THEN
          CALL print_date(0)
          ERROR STOP ERROR_cbh
        ENDIF

        IF ( Climate_temp_flag==ACTIVE ) THEN
          Basin_tmax = Basin_tmax*Basin_area_inv
          Basin_tmin = Basin_tmin*Basin_area_inv
          Basin_temp = Basin_temp*Basin_area_inv
          Solrad_tmax = SNGL( Basin_tmax )
          Solrad_tmin = SNGL( Basin_tmin )
        ENDIF

        IF ( Climate_precip_flag==ACTIVE ) THEN
          Basin_ppt = Basin_ppt*Basin_area_inv
          Basin_obs_ppt = sum_obs*Basin_area_inv
          Basin_rain = Basin_rain*Basin_area_inv
          Basin_snow = Basin_snow*Basin_area_inv
        ENDIF
        IF ( Climate_potet_flag==ACTIVE ) Basin_potet = Basin_potet*Basin_area_inv
        IF ( Climate_swrad_flag==ACTIVE ) THEN
          Basin_horad = Soltab_basinpotsw(Jday)
          IF ( Orad_flag==OFF ) Orad = SNGL( (DBLE(Swrad(1))*Hru_cossl(1)*Basin_horad)/Soltab_potsw(Jday,1) ) ! ??bad assumption using HRU 1
          Basin_swrad = Basin_swrad*Basin_area_inv
          Basin_potsw = Basin_swrad
        ENDIF
        IF ( Humidity_cbh_flag==ACTIVE ) Basin_humidity = Basin_humidity*Basin_area_inv
        IF ( Windspeed_cbh_flag==ACTIVE ) Basin_windspeed = Basin_windspeed*Basin_area_inv
        IF ( AET_cbh_flag==ACTIVE ) Basin_aet_external = Basin_aet_external*Basin_area_inv
        IF ( PET_cbh_flag==ACTIVE ) Basin_pet_external = Basin_pet_external*Basin_area_inv
        IF ( irrigated_area_flag==ACTIVE ) Basin_irrigated_area = Basin_irrigated_area*Basin_area_inv

      ELSEIF ( Process_flag==DECL ) THEN

        IF ( control_integer(Cbh_check_flag, 'cbh_check_flag')/=0 ) Cbh_check_flag = ACTIVE
        IF ( control_integer(Cbh_binary_flag, 'cbh_binary_flag')/=0 ) Cbh_binary_flag = OFF

        IF ( Climate_temp_flag==ACTIVE .OR. Model==DOCUMENTATION ) &
     &       CALL print_module('Temperature Distribution', MODNAME, Version_climate_hru)
        IF ( Climate_precip_flag==ACTIVE .OR. Model==DOCUMENTATION ) &
     &       CALL print_module('Precipitation Distribution', MODNAME, Version_climate_hru)
        IF ( Climate_swrad_flag==ACTIVE .OR. Model==DOCUMENTATION ) &
     &       CALL print_module('Solar Radiation Distribution', MODNAME, Version_climate_hru)
        IF ( Climate_potet_flag==ACTIVE .OR. Model==DOCUMENTATION ) &
     &       CALL print_module('Potential Evapotranspiration', MODNAME, Version_climate_hru)
        IF ( Climate_transp_flag==ACTIVE .OR. Model==DOCUMENTATION ) &
     &       CALL print_module('Transpiration Distribution', MODNAME, Version_climate_hru)

        IF ( Humidity_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Humidity_hru(Nhru) )
          CALL declvar_real(MODNAME, 'humidity_hru', 'nhru', Nhru, &
     &         'Relative humidity of each HRU read from CBH File', &
     &         'percentage', Humidity_hru)
        ENDIF
        IF ( Albedo_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Albedo_hru(Nhru) )
          CALL declvar_real(MODNAME, 'albedo_hru', 'nhru', Nhru, &
     &         'Snowpack albedo of each HRU read from CBH File', &
     &         'decimal fraction', Albedo_hru)
        ENDIF
        IF ( Cloud_cover_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Cloud_cover_cbh(Nhru) )
          CALL declvar_real(MODNAME, 'cloud_cover_cbh', 'nhru', Nhru, &
     &         'Cloud_cover of each HRU read from CBH File', &
     &         'decimal fraction', Cloud_cover_cbh)
        ENDIF
        IF ( Windspeed_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          CALL declvar_dble(MODNAME, 'basin_windspeed', 'one', 1, &
     &         'Basin area-weighted average wind speed', &
     &         'meters/second', Basin_windspeed)
          ALLOCATE ( Windspeed_hru(Nhru) )
          CALL declvar_real(MODNAME, 'windspeed_hru', 'nhru', Nhru, &
     &         'Wind speed for each HRU read from CBH File', &
     &         'meters/second', Windspeed_hru)
        ENDIF

        IF ( AET_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          CALL declvar_dble(MODNAME, 'basin_aet_external', 'one', 1, &
     &         'Basin area-weighted average PET read from CBH File', &
     &         'inches', Basin_aet_external)
          ALLOCATE ( AET_external(Nhru) )
          CALL declvar_real(MODNAME, 'AET_external', 'nhru', Nhru, &
     &         'Actual evapotranspiration for each HRU that is read from a CBH File', &
     &         'inches', AET_external)
        ENDIF

        IF ( PET_cbh_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          CALL declvar_dble(MODNAME, 'basin_pet_external', 'one', 1, &
     &         'Basin area-weighted average PET read from CBH File', &
     &         'inches', Basin_pet_external)
          ALLOCATE ( PET_external(Nhru) )
          CALL declvar_real(MODNAME, 'PET_external', 'nhru', Nhru, &
     &         'Potential evapotranspiration for each HRU that is read from a CBH File', &
     &         'inches', PET_external)
        ENDIF

        IF ( irrigated_area_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          CALL declvar_dble(MODNAME, 'basin_irrigated_area', 'one', 1, &
     &         'Basin area-weighted average irrigation area read from CBH File', &
     &         'acres', Basin_irrigated_area)
          ALLOCATE ( Irrigated_area(Nhru) )
          CALL declvar_real(MODNAME, 'irrigated_area', 'nhru', Nhru, &
     &         'Irrigated area for each HRU that is read from a CBH File', &
     &         'acres', Irrigated_area)
        ENDIF

!   Declared Parameters
        IF ( Climate_temp_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Tmax_cbh_adj(Nhru,MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'tmax_cbh_adj', 'nhru,nmonths', 'real', &
     &         '0.0', '-10.0', '10.0', &
     &         'Monthly maximum temperature adjustment factor for each HRU', &
     &         'Monthly (January to December) adjustment factor to maximum air temperature for each HRU,'// &
     &         ' estimated on the basis of slope and aspect', &
     &         'temp_units')/=0 ) CALL read_error(1, 'tmax_cbh_adj')

          ALLOCATE ( Tmin_cbh_adj(Nhru,MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'tmin_cbh_adj', 'nhru,nmonths', 'real', &
     &         '0.0', '-10.0', '10.0', &
     &         'Monthly minimum temperature adjustment factor for each HRU', &
     &         'Monthly (January to December) adjustment factor to minimum air temperature for each HRU,'// &
     &         ' estimated on the basis of slope and aspect', &
     &         'temp_units')/=0 ) CALL read_error(1, 'tmin_cbh_adj')
        ENDIF

        IF ( Climate_precip_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Rain_cbh_adj(Nhru,MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'rain_cbh_adj', 'nhru,nmonths', 'real', &
     &         '1.0', '0.5', '2.0', &
     &         'Rain adjustment factor, by month for each HRU', &
     &         'Monthly (January to December) adjustment factor to'// &
     &         ' measured precipitation determined to be rain on'// &
     &         ' each HRU to account for differences in elevation, and so forth', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'rain_cbh_adj')

          ALLOCATE ( Snow_cbh_adj(Nhru,MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'snow_cbh_adj', 'nhru,nmonths', 'real', &
     &         '1.0', '0.5', '2.0', &
     &         'Snow adjustment factor, by month for each HRU', &
     &         'Monthly (January to December) adjustment factor to'// &
     &         ' measured precipitation determined to be snow on'// &
     &         ' each HRU to account for differences in elevation, and so forth', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'snow_cbh_adj')
        ENDIF

        IF ( Climate_potet_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Potet_cbh_adj(Nhru,MONTHS_PER_YEAR) )
          IF ( declparam(MODNAME, 'potet_cbh_adj', 'nhru,nmonths', 'real', &
     &         '1.0', '0.5', '1.5', &
     &         'Potential ET adjustment factor, by month for each HRU', &
     &         'Monthly (January to December) adjustment factor to'// &
     &         ' potential evapotranspiration specified in CBH Files for each HRU', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'potet_cbh_adj')
        ENDIF

      ELSEIF ( Process_flag==INIT ) THEN
        Basin_humidity = 0.0D0
        Basin_windspeed = 0.0D0
        IF ( Humidity_cbh_flag==ACTIVE ) Humidity_hru = 0.0
        IF ( Windspeed_cbh_flag==ACTIVE ) Windspeed_hru = 0.0

        istop = 0
        ierr = 0

        IF ( Climate_precip_flag==ACTIVE ) THEN
          IF ( getparam_real(MODNAME, 'rain_cbh_adj', Nhru*MONTHS_PER_YEAR, Rain_cbh_adj)/=0 ) CALL read_error(2, 'rain_cbh_adj')
          IF ( getparam_real(MODNAME, 'snow_cbh_adj', Nhru*MONTHS_PER_YEAR, Snow_cbh_adj)/=0 ) CALL read_error(2, 'snow_cbh_adj')

          IF ( control_string(Precip_day, 'precip_day')/=0 ) CALL read_error(5, 'precip_day')
          CALL find_header_end(Precip_unit, Precip_day, 'precip_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Precip_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Precip_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Climate_temp_flag==ACTIVE ) THEN
          IF ( getparam_real(MODNAME, 'tmax_cbh_adj', Nhru*MONTHS_PER_YEAR, Tmax_cbh_adj)/=0 ) CALL read_error(2, 'tmax_cbh_adj')
          IF ( getparam_real(MODNAME, 'tmin_cbh_adj', Nhru*MONTHS_PER_YEAR, Tmin_cbh_adj)/=0 ) CALL read_error(2, 'tmin_cbh_adj')

          IF ( control_string(Tmax_day, 'tmax_day')/=0 ) CALL read_error(5, 'tmax_day')
          IF ( control_string(Tmin_day, 'tmin_day')/=0 ) CALL read_error(5, 'tmin_day')
          CALL find_header_end(Tmax_unit, Tmax_day, 'tmax_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Tmax_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Tmax_day
              istop = 1
            ENDIF
          ENDIF
          CALL find_header_end(Tmin_unit, Tmin_day, 'tmin_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Tmin_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Tmin_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Climate_potet_flag==ACTIVE ) THEN
          IF ( getparam_real(MODNAME, 'potet_cbh_adj', Nhru*MONTHS_PER_YEAR, Potet_cbh_adj)/=0 ) &
     &         CALL read_error(2, 'potet_cbh_adj')
          IF ( control_string(Potet_day, 'potet_day')/=0 ) CALL read_error(5, 'potet_day')
          CALL find_header_end(Et_unit, Potet_day, 'potet_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Et_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Potet_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( AET_cbh_flag==ACTIVE ) THEN
          IF ( control_string(AET_cbh_file, 'AET_cbh_file')/=0 ) CALL read_error(5, 'AET_cbh_file')
          CALL find_header_end(AET_unit, AET_cbh_file, 'AET_cbh_file', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(AET_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', AET_cbh_file
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( PET_cbh_flag==ACTIVE ) THEN
          IF ( control_string(PET_cbh_file, 'PET_cbh_file')/=0 ) CALL read_error(5, 'PET_cbh_file')
          CALL find_header_end(PET_unit, PET_cbh_file, 'PET_cbh_file', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(PET_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', PET_cbh_file
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( irrigated_area_flag==ACTIVE ) THEN
          IF ( control_string(irrigated_area_cbh_file, 'irrigated_cbh_file')/=0 ) CALL read_error(5, 'irrigated_area_cbh_file')
          CALL find_header_end(irrigated_area_unit, irrigated_area_cbh_file, 'irrigated_area_cbh_file', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(irrigated_area_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', irrigated_area_cbh_file
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Climate_transp_flag==ACTIVE ) THEN
          IF ( control_string(Transp_day, 'transp_day')/=0 ) CALL read_error(5, 'transp_day')
          CALL find_header_end(Transp_unit, Transp_day, 'transp_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Transp_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Transp_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Climate_swrad_flag==ACTIVE ) THEN
          IF ( control_string(Swrad_day, 'swrad_day')/=0 ) CALL read_error(5, 'swrad_day')
          CALL find_header_end(Swrad_unit, Swrad_day, 'swrad_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Swrad_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Swrad_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Humidity_cbh_flag==ACTIVE ) THEN
          IF ( control_string(Humidity_day, 'humidity_day')/=0 ) CALL read_error(5, 'humidity_day')
          ierr = 2 ! signals routine to ignore CBH file requirement and use a parameter
          CALL find_header_end(Humidity_unit, Humidity_day, 'humidity_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSEIF ( ierr==2 ) THEN
            Humidity_cbh_flag = 0
          ELSE
            CALL find_current_time(Humidity_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Humidity_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Windspeed_cbh_flag==ACTIVE ) THEN
          IF ( control_string(Windspeed_day, 'windspeed_day')/=0 ) CALL read_error(5, 'windspeed_day')
          CALL find_header_end(Windspeed_unit, Windspeed_day, 'windspeed_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Windspeed_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Windspeed_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Albedo_cbh_flag==ACTIVE ) THEN
          IF ( control_string(Albedo_day, 'albedo_day')/=0 ) CALL read_error(5, 'albedo_day')
          CALL find_header_end(Albedo_unit, Albedo_day, 'albedo_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Albedo_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Albedo_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( Cloud_cover_cbh_flag==ACTIVE ) THEN
          IF ( control_string(Cloud_cover_day, 'cloud_cover_day')/=0 ) CALL read_error(5, 'cloud_cover_day')
          CALL find_header_end(Cloud_cover_unit, Cloud_cover_day, 'cloud_cover_day', ierr, 1, Cbh_binary_flag)
          IF ( ierr==1 ) THEN
            istop = 1
          ELSE
            CALL find_current_time(Cloud_cover_unit, Start_year, Start_month, Start_day, ierr, Cbh_binary_flag)
            IF ( ierr==-1 ) THEN
              PRINT *, 'for first time step, CBH File: ', Cloud_cover_day
              istop = 1
            ENDIF
          ENDIF
        ENDIF

        IF ( istop==ACTIVE ) THEN
          PRINT *, 'ERROR in climate_hru'
          ERROR STOP ERROR_cbh
        ENDIF

      ENDIF

      END FUNCTION climate_hru

!***********************************************************************
!     Read a day in the CBH File
!***********************************************************************
      SUBROUTINE read_cbh_date(Year, Month, Day, Var, Ios, Iret)
      USE PRMS_MODULE, ONLY: Nowyear, Nowmonth, Nowday
      use prms_utils, only: print_date
! Argument
      INTEGER, INTENT(IN) :: Year, Month, Day, Ios
      CHARACTER(LEN=*), INTENT(IN) :: Var
      INTEGER, INTENT(INOUT) :: Iret
! Local Variables
      INTEGER :: right_day
!***********************************************************************
      right_day = 1
      IF ( Year/=Nowyear .OR. Month/=Nowmonth .OR. Day/=Nowday ) right_day = 0
      IF ( Ios/=0 .OR. right_day==0 ) THEN
        PRINT *, 'ERROR, reading CBH File, variable: ', Var, ' IOSTAT=', Ios
        IF ( Ios==-1 ) THEN
          PRINT *, '       End-of-File found'
        ELSEIF ( right_day==0 ) THEN
          PRINT *, '       Wrong day found'
        ELSE
          PRINT *, '       Invalid data value found'
        ENDIF
        CALL print_date(0)
        Iret = 1
      ENDIF
      END SUBROUTINE read_cbh_date

!***********************************************************************
!     Check CBH value limits
!***********************************************************************
      SUBROUTINE check_cbh_value(Var, Var_value, Lower_val, Upper_val, Missing)
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
! Argument
      REAL, INTENT(IN) :: Var_value(*), Lower_val, Upper_val
      CHARACTER(LEN=*), INTENT(IN) :: Var
      INTEGER, INTENT(INOUT) :: Missing
! Functions
      !INTRINSIC :: ISNAN
      !EXTERNAL :: print_date
! Local Variables
      INTEGER :: i, j
!***********************************************************************
      DO j = 1, Active_hrus
        i = Hru_route_order(j)
        !IF ( ISNAN(Var_value(i)) ) THEN
        !  PRINT *, 'ERROR, NaN value found for variable: ', Var
        !  Var_value(i) = 0.0
        !  Missing = 1
        !  CALL print_date(0)
        !ELSEIF ( Var_value(i)<Lower_val .OR. Var_value(i)>Upper_val ) THEN
        IF ( Var_value(i)<Lower_val .OR. Var_value(i)>Upper_val ) THEN
          PRINT '(3A,F0.8, A, I0)', 'ERROR, bad value, variable: ', Var, ' Value:', Var_value(i), ' HRU: ', i
          PRINT *,                  '       lower bound:', Lower_val, ' upper bound:', Upper_val
          Missing = 1
        ENDIF
      ENDDO
      END SUBROUTINE check_cbh_value

!***********************************************************************
!     Check CBH integer value limits
!***********************************************************************
      SUBROUTINE check_cbh_intvalue(Var, Var_value, Lower_val, Upper_val, Missing)
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
! Argument
      INTEGER, INTENT(IN) :: Var_value(*), Lower_val, Upper_val
      CHARACTER(LEN=*), INTENT(IN) :: Var
      INTEGER, INTENT(INOUT) :: Missing
      ! Local Variables
      INTEGER :: i, j
!***********************************************************************
      DO j = 1, Active_hrus
        i = Hru_route_order(j)
        IF ( Var_value(i)<Lower_val .OR. Var_value(i)>Upper_val ) THEN
          PRINT '(3A,I0, A, I0)', 'ERROR, bad value, variable: ', Var, ' Value:', Var_value(i), ', HRU: ', i
          PRINT *,                '       lower bound:', Lower_val, ' upper bound:', Upper_val
          Missing = 1
        ENDIF
      ENDDO
      END SUBROUTINE check_cbh_intvalue
