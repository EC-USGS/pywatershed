!***********************************************************************
! Maximum, minimum, average temperatures, and precipitation are distributed to each HRU
! using temperature and/or precipitation data input as a time series of gridded
! or other spatial units using an area-weighted method and correction factors to
! account for differences in altitude, spatial variation, topography, and 
! measurement gage efficiency
!***********************************************************************
      MODULE PRMS_TEMP_MAP
        USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, MAXFILE_LENGTH, MONTHS_PER_YEAR, temp_map_module
        USE PRMS_MODULE, ONLY: Model, Process_flag, Start_year, Start_month, Start_day, Nmap2hru, Nmap
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Temperature Distribution'
        character(len=*), parameter :: MODNAME = 'temp_map'
        character(len=*), parameter :: Version_temp_map = '2020-11-20'
        INTEGER, SAVE :: Tmax_unit, Tmin_unit
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Hru2map_id(:), Map2hru_id(:)
        REAL, SAVE, ALLOCATABLE :: Hru2map_pct(:), Tmax_map_values(:), Tmin_map_values(:), Precip_map_values(:)
        REAL, SAVE, ALLOCATABLE :: Tmax_map_adj(:, :), Tmin_map_adj(:, :), Precip_map_adj(:, :)
        ! parameters in basin:
        !    hru_area
        ! Control Parameters
        CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: Tmin_map_file, Tmax_map_file
      END MODULE PRMS_TEMP_MAP

      SUBROUTINE temp_map()
      USE PRMS_TEMP_MAP
      USE PRMS_BASIN, ONLY: Hru_area, Basin_area_inv, Active_hrus, Hru_route_order
      USE PRMS_CLIMATEVARS, ONLY: Solrad_tmax, Solrad_tmin, Basin_temp, &
     &    Basin_tmax, Basin_tmin, Tmaxf, Tminf, Tminc, Tmaxc, Tavgf, Tavgc
      USE PRMS_SET_TIME, ONLY: Nowmonth
! Functions
      INTRINSIC :: SNGL
      INTEGER, EXTERNAL :: declparam, getparam, getdim, decldim, control_string
      EXTERNAL :: read_error, precip_form, temp_set, find_header_end, find_current_time
      EXTERNAL :: read_cbh_date, print_module, print_date
! Local Variables
      INTEGER :: yr, mo, dy, i, hr, mn, sec, ierr, ios, j, kg, kh, istop
      REAL :: tmax_hru, tmin_hru, harea
!***********************************************************************
       IF ( Process_flag==RUN ) THEN
        READ ( Tmax_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tmax_map_values(i), i=1,Nmap)
        READ ( Tmin_unit, *, IOSTAT=ios ) yr, mo, dy, hr, mn, sec, (Tmin_map_values(i), i=1,Nmap)
        Basin_tmax = 0.0D0
        Basin_tmin = 0.0D0
        Basin_temp = 0.0D0
        Tmaxf = 0.0
        Tminf = 0.0

        DO j = 1, Nmap2hru
          kg = Map2hru_id(j)
          kh = Hru2map_id(j)
          Tmaxf(kh) = Tmaxf(kh) + Tmax_map_values(kg)*Hru2map_pct(j) + Tmax_map_adj(kg, Nowmonth)
          Tminf(kh) = Tminf(kh) + Tmin_map_values(kg)*Hru2map_pct(j) + Tmin_map_adj(kg, Nowmonth)
        ENDDO

        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          harea = Hru_area(i)
          tmax_hru = Tmaxf(i)
          tmin_hru = Tminf(i)
          CALL temp_set(i, tmax_hru, tmin_hru, Tmaxf(i), Tminf(i), &
     &                    Tavgf(i), Tmaxc(i), Tminc(i), Tavgc(i), harea)
        ENDDO
        Basin_tmax = Basin_tmax*Basin_area_inv
        Basin_tmin = Basin_tmin*Basin_area_inv
        Basin_temp = Basin_temp*Basin_area_inv
        Solrad_tmax = SNGL( Basin_tmax )
        Solrad_tmin = SNGL( Basin_tmin )

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_temp_map)

        ALLOCATE ( Tmax_map_values(Nmap), Tmin_map_values(Nmap) )

! Declare parameters
        ALLOCATE ( Tmax_map_adj(Nmap,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'tmax_map_adj', 'nmap,nmonths', 'real', &
     &       '0.0', '-10.0', '10.0', &
     &       'Monthly maximum temperature adjustment factor for each mapped spatial unit', &
     &       'Monthly (January to December) additive adjustment factor to maximum air temperature for each mapped,'// &
     &       ' spatial unit estimated on the basis of slope and aspect', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmax_map_adj')
        ALLOCATE ( Tmin_map_adj(Nmap,MONTHS_PER_YEAR) )
        IF ( declparam(MODNAME, 'tmin_map_adj', 'nmap,nmonths', 'real', &
     &       '0.0', '-10.0', '10.0', &
     &       'Monthly minimum temperature adjustment factor for each mapped spatial unit', &
     &       'Monthly (January to December) additive adjustment factor to minimum air temperature for each'// &
     &       ' mapped spatial unit, estimated on the basis of slope and aspect', &
     &       'temp_units')/=0 ) CALL read_error(1, 'tmin_map_adj')

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
        IF ( getparam(MODNAME, 'tmax_map_adj', Nmap*MONTHS_PER_YEAR, 'real', Tmax_map_adj)/=0 ) &
     &       CALL read_error(2, 'tmax_map_adj')
        IF ( getparam(MODNAME, 'tmin_map_adj', Nmap*MONTHS_PER_YEAR, 'real', Tmin_map_adj)/=0 ) &
     &       CALL read_error(2, 'tmin_map_adj')
        IF ( control_string(Tmax_map_file, 'tmax_map_file')/=0 ) CALL read_error(5, 'tmax_map_file')
        IF ( control_string(Tmin_map_file, 'tmin_map_file')/=0 ) CALL read_error(5, 'tmin_map_file')
        CALL find_header_end(Tmax_unit, Tmax_map_file, 'tmax_map_file', ierr, 1, 0)
        IF ( ierr==1 ) THEN
          istop = 1
        ELSE
          CALL find_current_time(Tmax_unit, Start_year, Start_month, Start_day, ierr, 0)
          IF ( ierr==-1 ) THEN
            PRINT *, 'for first time step, Tmax Map File: ', Tmax_map_file
            istop = 1
          ENDIF
        ENDIF
        CALL find_header_end(Tmin_unit, Tmin_map_file, 'tmin_map_file', ierr, 1, 0)
        IF ( ierr==1 ) THEN
          istop = 1
        ELSE
          CALL find_current_time(Tmin_unit, Start_year, Start_month, Start_day, ierr, 0)
          IF ( ierr==-1 ) THEN
            PRINT *, 'for first time step, Tmin Map File: ', Tmin_map_file
            istop = 1
          ENDIF
        ENDIF
        IF ( istop==1 ) STOP 'ERROR in temp_map module'
      ENDIF

      END SUBROUTINE temp_map
