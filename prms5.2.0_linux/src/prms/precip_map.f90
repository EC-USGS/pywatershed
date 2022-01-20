!***********************************************************************
! Maximum, minimum, average temperatures, and precipitation are distributed to each HRU
! using temperature and/or precipitation data input as a time series of gridded
! or other spatial units using an area-weighted method and correction factors to
! account for differences in altitude, spatial variation, topography, and 
! measurement gage efficiency
!***********************************************************************
      MODULE PRMS_PRECIP_MAP
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, OFF, MAXFILE_LENGTH, &
     &      MM, MM2INCH, MONTHS_PER_YEAR, precip_map_module
        USE PRMS_MODULE, ONLY: Model, Process_flag, Start_year, Start_month, Start_day, Nmap2hru, Nmap
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Precipitation Distribution'
        character(len=*), parameter :: MODNAME = 'precip_map'
        character(len=*), parameter :: Version_precip_map = '2020-11-20'
        INTEGER, SAVE :: Precip_unit
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Hru2map_id(:), Map2hru_id(:)
        REAL, SAVE, ALLOCATABLE :: Hru2map_pct(:), Precip_map_values(:)
        REAL, SAVE, ALLOCATABLE :: Precip_map_adj(:, :)
        ! parameters in basin:
        !    hru_area
        ! Control Parameters
        CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: Precip_map_file
      END MODULE PRMS_PRECIP_MAP

      SUBROUTINE precip_map()
      USE PRMS_PRECIP_MAP
      USE PRMS_BASIN, ONLY: Hru_area, Basin_area_inv, Active_hrus, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt, Hru_rain, Hru_snow, Prmx, Pptmix, Newsnow, &
     &    Precip_units, Tmax_allrain_f, Adjmix_rain, Tmaxf, Tminf, &
     &    Basin_ppt, Basin_snow, Basin_rain, Basin_obs_ppt, Tmax_allsnow_f
      USE PRMS_SET_TIME, ONLY: Nowmonth
! Functions
      INTRINSIC :: SNGL
      INTEGER, EXTERNAL :: declparam, getparam, control_string
      EXTERNAL :: read_error, precip_form, find_header_end, find_current_time
      EXTERNAL :: read_cbh_date, print_module, print_date
! Local Variables
      INTEGER :: yr, mo, dy, i, hr, mn, sec, ierr, ios, j, kg, kh, istop
      REAL :: ppt, harea
!***********************************************************************
       IF ( Process_flag==RUN ) THEN
        READ ( Precip_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Precip_map_values(i), i=1,Nmap)
        Basin_ppt = 0.0D0
        Basin_rain = 0.0D0
        Basin_snow = 0.0D0
        Basin_obs_ppt = 0.0D0
        Hru_ppt = 0.0

        DO j = 1, Nmap2hru
          kg = Map2hru_id(j)
          kh = Hru2map_id(j)
          Hru_ppt(kh) = Hru_ppt(kh) + Precip_map_values(kg)*Hru2map_pct(j)*Precip_map_adj(kg, Nowmonth)
        ENDDO

        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          harea = Hru_area(i)
!******Initialize HRU variables
          Pptmix(i) = OFF
          Newsnow(i) = OFF
          Prmx(i) = 0.0
          Hru_rain(i) = 0.0
          Hru_snow(i) = 0.0
          IF ( Hru_ppt(i)>0.0 ) THEN
            IF ( Precip_units==MM ) Hru_ppt(i) = Hru_ppt(i)*MM2INCH
            ppt = Hru_ppt(i)
            CALL precip_form(ppt, Hru_ppt(i), Hru_rain(i), Hru_snow(i), &
     &                       Tmaxf(i), Tminf(i), Pptmix(i), Newsnow(i), &
     &                       Prmx(i), Tmax_allrain_f(i,Nowmonth), 1.0, 1.0, &
     &                       Adjmix_rain(i,Nowmonth), harea, Basin_obs_ppt, Tmax_allsnow_f(i,Nowmonth))
          ELSEIF ( Hru_ppt(i)<0.0 ) THEN
            PRINT *, 'WARNING, negative precipitation value entered in precipitation map file and set to 0.0, HRU:', i
            CALL print_date(0)
            Hru_ppt(i) = 0.0
          ENDIF
        ENDDO
        Basin_ppt = Basin_ppt*Basin_area_inv
        Basin_obs_ppt = Basin_obs_ppt*Basin_area_inv
        Basin_rain = Basin_rain*Basin_area_inv
        Basin_snow = Basin_snow*Basin_area_inv

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_precip_map)

        ALLOCATE ( Precip_map_values(Nmap) )

! Declare parameters
        ALLOCATE ( Precip_map_adj(Nmap,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'precip_map_adj', 'nmap,nmonths', 'real', &
     &     '1.0', '0.5', '2.0', &
     &     'Monthly rain adjustment factor for each mapped spatial unit', &
     &     'Monthly (January to December) multiplicative adjustment factor to mapped precipitation to account for'// &
     &     ' differences in elevation, and so forth', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'precip_map_adj')

        ALLOCATE ( Hru2map_id(Nmap2hru) )
        IF ( declparam(MODNAME, 'hru2map_id', 'nmap2hru', 'integer', &
     &       '1', 'bounded', 'nhru', &
     &       'HRU identification number for each HRU to mapped spatial units intersection', &
     &       'HRU identification number for each HRU to mapped spatial units intersection', &
     &       'none')/=0 ) CALL read_error(1, 'hru2map_id')

        !rsr, bounded value could be a problem if number of mapped spatial units > nhru
        ALLOCATE ( Map2hru_id(Nmap2hru) )
        IF ( declparam(MODNAME, 'map2hru_id', 'nmap2hru', 'integer', &
     &       '0', 'bounded', 'nmap', &
     &       'Mapped spatial unit identification number for each HRU to map intersection', &
     &       'Mapped spatial unit identification number for each HRU to map intersection', &
     &       'none')/=0 ) CALL read_error(1, 'map2hru_id')

        ALLOCATE ( Hru2map_pct(Nmap2hru) )
        IF ( declparam(MODNAME, 'hru2map_pct', 'nmap2hru', 'real', &
     &       '0.0', '0.0', '1.0', &
     &       'Portion of HRU associated with each HRU to map intersection', &
     &       'Portion of HRU associated with each HRU to map intersection', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'hru2map_pct')

! Get parameters
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( getparam(MODNAME, 'map2hru_id', Nmap2hru, 'integer', Map2hru_id)/=0 ) CALL read_error(2, 'map2hru_id')
        IF ( getparam(MODNAME, 'hru2map_id', Nmap2hru, 'integer', Hru2map_id)/=0 ) CALL read_error(2, 'hru2map_id')
        IF ( getparam(MODNAME, 'hru2map_pct', Nmap2hru, 'real', Hru2map_pct)/=0 ) CALL read_error(2, 'hru2map_pct')

        istop = 0
        ierr = 0
        IF ( getparam(MODNAME, 'precip_map_adj', Nmap*MONTHS_PER_YEAR, 'real', Precip_map_adj)/=0 ) &
     &       CALL read_error(2, 'precip_map_adj')
        IF ( control_string(Precip_map_file, 'precip_map_file')/=0 ) CALL read_error(5, 'precip_map_file')
        CALL find_header_end(Precip_unit, Precip_map_file, 'precip_map_file', ierr, 1, 0)
        IF ( ierr==1 ) THEN
          istop = 1
        ELSE
          CALL find_current_time(Precip_unit, Start_year, Start_month, Start_day, ierr, 0)
          IF ( ierr==-1 ) THEN
            PRINT *, 'for first time step, Precip Map File: ', Precip_map_file
            istop = 1
          ENDIF
        ENDIF
        IF ( istop==1 ) STOP 'ERROR in precip_map module'
      ENDIF

      END SUBROUTINE precip_map
