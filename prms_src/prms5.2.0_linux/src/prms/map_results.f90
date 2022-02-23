!***********************************************************************
!     Output declared variables by HRU as a grid mapped to a specified
!     spatial resolution
!***********************************************************************
      MODULE PRMS_MAP_RESULTS
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, NEARZERO, REAL_TYPE, DBLE_TYPE, &
     &    DOCUMENTATION, ERROR_dim, DEBUG_less
      USE PRMS_MODULE, ONLY: Model, Nhru, Nhrucell, Ngwcell, MapOutON_OFF, &
     &    Print_debug, Inputerror_flag, Start_year, Start_month, Start_day, Parameter_check_flag, &
     &    Prms_warmup, End_year, End_month, End_day
      IMPLICIT NONE
! Module Variables
      character(len=*), parameter :: MODDESC = 'Output Summary'
      character(len=*), parameter :: MODNAME = 'map_results'
      character(len=*), parameter :: Version_map_results = '2020-12-02'
      INTEGER, SAVE :: Mapflg, Numvalues, Lastyear, Totdays
      INTEGER, SAVE :: Yrdays, Yrresults, Totresults, Monresults, Mondays
      INTEGER, SAVE :: Begin_results, Begyr, Dailyresults
      INTEGER, SAVE :: Prevyr, Prevmo, Prevday, Weekresults, Weekdays
      INTEGER, SAVE, ALLOCATABLE :: Totunit(:), Yrunit(:), Monunit(:), Weekunit(:), Dailyunit(:)
      INTEGER, SAVE, ALLOCATABLE :: Nc_vars(:), Map_var_type(:)
      DOUBLE PRECISION, SAVE :: Conv_fac
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Basin_var_tot(:), Basin_var_yr(:), Basin_var_mon(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gvr_map_frac_dble(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Basin_var_week(:), Map_var_id(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Map_var_mon(:, :), Map_var_week(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Map_var_yr(:, :), Map_var_tot(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Basin_var_daily(:), Map_var_daily(:, :)
      REAL, SAVE, ALLOCATABLE :: Map_var(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Map_var_dble(:, :)
      CHARACTER(LEN=15), SAVE :: Mapfmt
! Declared Parameters
      INTEGER, SAVE :: Ncol, Mapvars_freq, Mapvars_units
      INTEGER, SAVE, ALLOCATABLE :: Gvr_map_id(:), Gvr_hru_id(:)
      REAL, SAVE, ALLOCATABLE :: Gvr_map_frac(:)
! Control Parameters
      INTEGER, SAVE :: NmapOutVars
      CHARACTER(LEN=36), SAVE, ALLOCATABLE :: MapOutVar_names(:)
      END MODULE PRMS_MAP_RESULTS

!     ******************************************************************
!     Mapped results module
!     ******************************************************************
      INTEGER FUNCTION map_results()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN
      USE PRMS_MODULE, ONLY: Process_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: map_resultsdecl, map_resultsinit
      INTEGER, EXTERNAL :: map_resultsclean, map_resultsrun
!***********************************************************************
      map_results = 0

      IF ( Process_flag==RUN ) THEN
        map_results = map_resultsrun()
      ELSEIF ( Process_flag==DECL ) THEN
        map_results = map_resultsdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        map_results = map_resultsinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        map_results = map_resultsclean()
      ENDIF

      END FUNCTION map_results

!***********************************************************************
!     map_resultsdecl - declare parameters and variables
!***********************************************************************
      INTEGER FUNCTION map_resultsdecl()
      USE PRMS_MAP_RESULTS
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, control_string_array, control_integer
      EXTERNAL :: read_error, print_module, error_stop
! Local Variables
      INTEGER :: i
!***********************************************************************
      map_resultsdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_map_results)

      IF ( control_integer(NmapOutVars, 'nmapOutVars')/=0 ) NmapOutVars = 0
      IF ( NmapOutVars==0 ) THEN
        IF ( Model/=DOCUMENTATION ) THEN
          PRINT *, 'WARNING, map_results requested with nmapOutVars equal 0'
          PRINT *, 'no map_results output is produced'
          MapOutON_OFF = OFF
          RETURN
        ENDIF
      ELSE
        ALLOCATE ( MapOutVar_names(NmapOutVars), Map_var_type(NmapOutVars), Nc_vars(NmapOutVars) )
        ALLOCATE ( Map_var(Nhru, NmapOutVars), Map_var_dble(Nhru, NmapOutVars) )
        MapOutVar_names = ' '
        DO i = 1, NmapOutVars
          IF ( control_string_array(MapOutVar_names(i), 'mapOutVar_names', i)/=0 ) CALL read_error(5, 'mapOutVar_names')
        ENDDO
      ENDIF

      Mapflg = ACTIVE
      Numvalues = Nhru
      IF ( (Nhru/=Ngwcell .AND. Ngwcell/=0) .OR. mapOutON_OFF==2 ) THEN
        IF ( Ngwcell==0 ) CALL error_stop('in map_results, ngwcell must be specified > 0', ERROR_dim)
        Mapflg = OFF
        Numvalues = Ngwcell
      ENDIF

      IF ( Mapflg==OFF .OR. Model==DOCUMENTATION ) THEN
        IF ( Nhrucell<1 ) THEN
          IF ( Model==DOCUMENTATION ) THEN
            Nhrucell = 1
          ELSE
            CALL error_stop('in map_results, nhrucell = 0 and must be > 0', ERROR_dim)
          ENDIF
        ENDIF
        ALLOCATE ( Gvr_map_id(Nhrucell), Gvr_map_frac(Nhrucell), Gvr_hru_id(Nhrucell), Map_var_id(Ngwcell) )
        ALLOCATE ( Gvr_map_frac_dble(Nhrucell) )
      ENDIF

! Declared Parameters
      IF ( declparam(MODNAME, 'mapvars_freq', 'one', 'integer', &
     &     '0', '0', '7', &
     &     'Mapped results frequency', &
     &     'Flag to specify the output frequency (0=none; 1=monthly; 2=yearly; 3=total; 4=monthly and yearly;'// &
     &     ' 5=monthly, yearly, and total; 6=weekly; 7=daily)', &
     &     'none')/=0 ) CALL read_error(1, 'mapvars_freq')

      IF ( declparam(MODNAME, 'mapvars_units', 'one', 'integer', &
     &     '0', '0', '3', &
     &     'Flag to specify the output units of mapped results', &
     &     'Flag to specify the output units of mapped results (0=units of the variable;'// &
     &     ' 1=inches to feet; 2=inches to centimeters; 3=inches to meters; as states or fluxes)', &
     &     'none')/=0 ) CALL read_error(1, 'mapvars_units')

      IF ( declparam(MODNAME, 'ncol', 'one', 'integer', &
     &     '1', '1', '50000', &
     &     'Number of columns for each row of the mapped results', &
     &     'Number of columns for each row of the mapped results', &
     &     'none')/=0 ) CALL read_error(1, 'ncol')

      IF ( Mapflg==OFF .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(MODNAME, 'gvr_cell_id', 'nhrucell', 'integer', &
     &       '0', 'bounded', 'ngwcell', &
     &       'Corresponding grid cell id associated with each GVR', &
     &       'Index of the grid cell associated with each gravity reservoir', &
     &       'none')/=0 ) CALL read_error(1, 'gvr_cell_id')
        IF ( Nhrucell/=Ngwcell .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'gvr_cell_pct', 'nhrucell', 'real', &
     &         '0.0', '0.0', '1.0', &
     &         'Proportion of the grid cell associated with each GVR', &
     &         'Proportion of the grid cell area associated with each gravity reservoir', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'gvr_cell_pct')
        ENDIF
        IF ( Nhru/=Nhrucell .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'gvr_hru_id', 'nhrucell', 'integer', &
     &         '0', 'bounded', 'nhru', &
     &         'Corresponding HRU id of each GVR', &
     &         'Index of the HRU associated with each gravity reservoir', &
     &         'none')/=0 ) CALL read_error(1, 'gvr_hru_id')
        ENDIF
      ENDIF

      END FUNCTION map_resultsdecl

!***********************************************************************
!     map_resultsinit - Initialize map_results module
!***********************************************************************
      INTEGER FUNCTION map_resultsinit()
      USE PRMS_MAP_RESULTS
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
      INTEGER, EXTERNAL :: getparam, getvartype, numchars, getvarsize
      EXTERNAL :: read_error, PRMS_open_output_file, checkdim_bounded_limits
! Local Variables
      INTEGER :: i, jj, is, ios, ierr, size, dim
      REAL, ALLOCATABLE, DIMENSION(:) :: map_frac
!***********************************************************************
      map_resultsinit = 0

      Begin_results = ACTIVE
      IF ( Prms_warmup>0 ) Begin_results = OFF
      Begyr = Start_year + Prms_warmup
      Lastyear = Begyr

      IF ( getparam(MODNAME, 'mapvars_freq', 1, 'integer', Mapvars_freq)/=0 ) CALL read_error(1, 'mapvars_freq')
      IF ( Mapvars_freq==0 ) THEN
        PRINT *, 'WARNING, map_results requested with mapvars_freq equal 0'
        PRINT *, 'no map_resultsults output is produced'
        MapOutON_OFF = OFF
        RETURN
      ENDIF
      IF ( getparam(MODNAME, 'ncol', 1, 'integer', Ncol)/=0 ) CALL read_error(2, 'ncol')
      IF ( getparam(MODNAME, 'mapvars_units', 1, 'integer', Mapvars_units)/=0 ) CALL read_error(2, 'Mapvars_units')

      WRITE ( Mapfmt, 9001 ) Ncol

      IF ( Mapvars_units==0 ) THEN
        Conv_fac = 1.0D0
      ELSEIF ( Mapvars_units==1 ) THEN
        Conv_fac = 1.0D0/12.0D0
      ELSEIF ( Mapvars_units==2 ) THEN
        Conv_fac = 2.54D0
      ELSE
        Conv_fac = 0.0254D0
      ENDIF
      Lastyear = Begyr
      Totdays = 0
      Mondays = 0
      Yrdays = 0
      Weekdays = 0
      Prevyr = Begyr
      Prevmo = Start_month
      Prevday = Start_day
      Map_var_type = REAL_TYPE
      ierr = 0
      DO jj = 1, NmapOutVars
        Nc_vars(jj) = numchars(MapOutVar_names(jj))
        Map_var_type(jj) = getvartype(MapOutVar_names(jj)(:Nc_vars(jj)), Map_var_type(jj) )
        IF ( Map_var_type(jj)/=REAL_TYPE .AND. Map_var_type(jj)/=DBLE_TYPE ) THEN
          PRINT *, 'ERROR, invalid map_results variable:', MapOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only real or double variables allowed'
          ierr = 1
        ENDIF
        size = getvarsize(MapOutVar_names(jj)(:Nc_vars(jj)), dim )
        IF ( size/=Nhru ) THEN
          PRINT *, 'ERROR, invalid map_results variable:', MapOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only variables with the number of values equal to nhru allowed'
          PRINT *, '       other possible dimensions are nssr, ngw, nhrucell, and ngwcell'
          ierr = 1
        ENDIF
      ENDDO

      Dailyresults = OFF
      Weekresults = OFF
      Monresults = OFF
      Yrresults = OFF
      Totresults = OFF

      IF ( Mapvars_freq==7 ) THEN
        Dailyresults = ACTIVE
        ALLOCATE ( Map_var_daily(Nhru, NmapOutVars), Basin_var_daily(NmapOutVars) )
        Map_var_daily = 0.0D0
        Basin_var_daily = 0.0D0
        ALLOCATE ( Dailyunit(NmapOutVars) )
        DO jj = 1, NmapOutVars
          CALL PRMS_open_output_file(Dailyunit(jj), MapOutVar_names(jj)(:Nc_vars(jj))//'.daily', 'xxx', 0, ios)
          IF ( ios/=0 ) ierr = 1
        ENDDO
      ENDIF

      IF ( Mapvars_freq==6 ) THEN
        Weekresults = ACTIVE
        ALLOCATE ( Map_var_week(Nhru, NmapOutVars), Basin_var_week(NmapOutVars) )
        Map_var_week = 0.0D0
        Basin_var_week = 0.0D0
        ALLOCATE ( Weekunit(NmapOutVars) )
        DO jj = 1, NmapOutVars
          CALL PRMS_open_output_file(Weekunit(jj), MapOutVar_names(jj)(:Nc_vars(jj))//'.weekly', 'xxx', 0, ios)
          IF ( ios/=0 ) ierr = 1
        ENDDO
      ENDIF

      IF ( Mapvars_freq==1 .OR. Mapvars_freq==4 .OR. Mapvars_freq==5 ) THEN
        Monresults = ACTIVE
        ALLOCATE ( Map_var_mon(Nhru, NmapOutVars), Basin_var_mon(NmapOutVars) )
        Map_var_mon = 0.0D0
        Basin_var_mon = 0.0D0
        ALLOCATE ( Monunit(NmapOutVars) )
        DO jj = 1, NmapOutVars
          CALL PRMS_open_output_file(Monunit(jj), MapOutVar_names(jj)(:Nc_vars(jj))//'.monthly', 'xxx', 0, ios)
          IF ( ios/=0 ) ierr = 1
        ENDDO
      ENDIF

      IF ( Mapvars_freq==2 .OR. Mapvars_freq==4 .OR. Mapvars_freq==5 ) THEN
        Yrresults = ACTIVE
        ALLOCATE ( Map_var_yr(Nhru, NmapOutVars) )
        ALLOCATE ( Basin_var_yr(NmapOutVars) )
        Map_var_yr = 0.0D0
        Basin_var_yr = 0.0D0
        ALLOCATE ( Yrunit(NmapOutVars) )
        DO jj = 1, NmapOutVars
          CALL PRMS_open_output_file(Yrunit(jj), MapOutVar_names(jj)(:Nc_vars(jj))//'.yearly', 'xxx', 0, ios)
          IF ( ios/=0 ) ierr = 1
        ENDDO
      ENDIF

      IF ( Mapvars_freq==3 .OR. Mapvars_freq==5 ) THEN
        Totresults = ACTIVE
        ALLOCATE ( Map_var_tot(Nhru, NmapOutVars) )
        ALLOCATE ( Basin_var_tot(NmapOutVars) )
        Map_var_tot = 0.0D0
        Basin_var_tot = 0.0D0
        ALLOCATE ( Totunit(NmapOutVars) )
        DO jj = 1, NmapOutVars
          CALL PRMS_open_output_file(Totunit(jj), MapOutVar_names(jj)(:Nc_vars(jj))//'.total', 'xxx', 0, ios)
          IF ( ios/=0 ) ierr = 1
        ENDDO
      ENDIF
      IF ( ierr==1 ) Inputerror_flag = 1

      IF ( Mapflg==OFF ) THEN
        IF ( getparam(MODNAME, 'gvr_cell_id', Nhrucell, 'integer', Gvr_map_id)/=0 ) CALL read_error(2, 'gvr_cell_id')
        IF ( Nhru/=Nhrucell ) THEN
          IF ( getparam(MODNAME, 'gvr_hru_id', Nhrucell, 'integer', Gvr_hru_id)/=0 ) CALL read_error(2, 'gvr_hru_id')
          IF ( Parameter_check_flag>0 ) &
     &         CALL checkdim_bounded_limits('gvr_hru_id', 'nhru', Gvr_hru_id, Nhrucell, 1, Nhru, Inputerror_flag)
        ELSE
          DO i = 1, Nhrucell
            Gvr_hru_id(i) = i
          ENDDO
        ENDIF
        IF ( Nhrucell/=Ngwcell ) THEN
          IF ( getparam(MODNAME, 'gvr_cell_pct', Nhrucell, 'real', Gvr_map_frac)/=0 ) CALL read_error(2, 'gvr_cell_pct')
        ELSE
          Gvr_map_frac = 1.0
        ENDIF

        IF ( Parameter_check_flag>0 ) THEN
          ALLOCATE ( map_frac(Ngwcell) )
          map_frac = 0.0
          DO i = 1, Nhrucell
            IF ( Gvr_map_id(i)>0 .AND. Gvr_hru_id(i)>0 .AND. Gvr_map_frac(i)>0.0 ) THEN
              is = Gvr_map_id(i)
              map_frac(is) = map_frac(is) + Gvr_map_frac(i)
            ELSEIF ( Print_debug>DEBUG_less ) THEN
              PRINT *, 'WARNING, map intersection:', i, ' ignored as specification is incomplete'
              PRINT *, '         HRU:', Gvr_hru_id(i), ' Map unit:', Gvr_map_id(i), 'Fraction:', Gvr_map_frac(i)
            ENDIF
          ENDDO

          DO i = 1, Ngwcell
            IF ( map_frac(i)<0.0 ) THEN
              PRINT *, 'ERROR, accounting for area of mapped spatial unit is negative:'
              PRINT *, '       Map id:', i, ' Fraction:', map_frac(i)
              Inputerror_flag = 1
            ELSEIF ( map_frac(i)<NEARZERO ) THEN
              CYCLE
            ELSEIF ( map_frac(i)>1.0001 ) THEN
              IF ( Print_debug>DEBUG_less ) THEN
                PRINT *, 'WARNING, excess accounting for area of mapped spatial unit:'
                PRINT *, '         Map id:', i, ' Fraction:', map_frac(i)
              ENDIF
            ELSEIF ( map_frac(i)<0.9999 ) THEN
              IF ( Print_debug>DEBUG_less ) THEN
                PRINT *, 'WARNING, incomplete accounting for area of mapped spatial unit'
                PRINT *, '         Map unit:', i, 'Fraction:', map_frac(i)
              ENDIF
            ENDIF
          ENDDO
          DEALLOCATE ( map_frac )
        ENDIF

        DO i = 1, Nhrucell
          Gvr_map_frac_dble(i) = DBLE( Gvr_map_frac(i) )
        ENDDO
      ENDIF

 9001 FORMAT ( '(',I8,'E10.3)' )

      END FUNCTION map_resultsinit

!***********************************************************************
!     map_resultsrun - Outputs various declared variables in grid format
!                      mapped to a specified spatial resolution
!***********************************************************************
      INTEGER FUNCTION map_resultsrun()
      USE PRMS_MAP_RESULTS
      USE PRMS_BASIN, ONLY: Hru_area_dble, Active_hrus, Hru_route_order, Basin_area_inv
      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Nowday, Modays
      IMPLICIT NONE
! FUNCTIONS AND SUBROUTINES
      INTRINSIC :: DBLE
      INTEGER, EXTERNAL :: getvar
      EXTERNAL :: read_error, write_results
! Local Variables
      INTEGER :: j, i, k, jj, last_day
      DOUBLE PRECISION :: factor, map_var_double
!***********************************************************************
      map_resultsrun = 0

      IF ( Begin_results==OFF ) THEN
        IF ( Nowyear==Begyr .AND. Nowmonth==Start_month .AND. Nowday==Start_day ) THEN
          Begin_results = ACTIVE
        ELSE
          RETURN
        ENDIF
      ENDIF

! check for last day of simulation
      last_day = 0
      IF ( Nowyear==End_year ) THEN
        IF ( Nowmonth==End_month ) THEN
          IF ( Nowday==End_day ) THEN
            last_day = 1
            Prevyr = Nowyear
            Prevmo = Nowmonth
            Prevday = Nowday
          ENDIF
        ENDIF
      ENDIF

      IF ( Yrresults==ACTIVE ) THEN
! check for first time step of the next year
        IF ( Lastyear/=Nowyear .OR. last_day==1 ) THEN
          IF ( (Nowmonth==Start_month .AND. Nowday==Start_day) .OR. last_day==1 ) THEN
            Lastyear = Nowyear
            factor = Conv_fac/DBLE(Yrdays)
            Basin_var_yr = 0.0D0
            DO jj = 1, NmapOutVars
              DO j = 1, Active_hrus
                i = Hru_route_order(j)
                Map_var_yr(i, jj) = Map_var_yr(i, jj)*factor
                Basin_var_yr(jj) = Basin_var_yr(jj) + Map_var_yr(i, jj)*Hru_area_dble(i)
              ENDDO
              Basin_var_yr(jj) = Basin_var_yr(jj)*Basin_area_inv

              WRITE ( Yrunit(jj), 9002 ) Prevyr, Prevmo, Prevday, ' Basin yearly mean:', Basin_var_yr(jj)

              IF ( Mapflg==ACTIVE ) THEN
                CALL write_results(Yrunit(jj), Map_var_yr(1,jj))
              ELSE
                Map_var_id = 0.0D0
                DO k = 1, Nhrucell
                  IF ( Gvr_map_id(k)==0 ) CYCLE
                  IF ( Gvr_map_frac(k)<NEARZERO ) CYCLE
                  Map_var_id(Gvr_map_id(k)) = Map_var_id(Gvr_map_id(k)) + &
     &                                        Map_var_yr(Gvr_hru_id(k), jj)*Gvr_map_frac_dble(k)
                ENDDO
                CALL write_results(Yrunit(jj), Map_var_id)
              ENDIF
              WRITE ( Yrunit(jj), 9003 )
            ENDDO
            Map_var_yr = 0.0D0
            Yrdays = 0
          ENDIF
        ENDIF
        Yrdays = Yrdays + 1
        Prevyr = Nowyear
        Prevmo = Nowmonth
        Prevday = Nowday
      ENDIF

!-----------------------------------------------------------------------
! need getvars for each variable (only can have short string)
      DO jj = 1, NmapOutVars
        IF ( Map_var_type(jj)==REAL_TYPE ) THEN
          IF ( getvar(MODNAME, MapOutVar_names(jj)(:Nc_vars(jj)), &
     &         Nhru, 'real', Map_var(1, jj))/=0 ) &
     &         CALL read_error(4, MapOutVar_names(jj)(:Nc_vars(jj)))
        ELSEIF ( Map_var_type(jj)==DBLE_TYPE ) THEN
          IF ( getvar(MODNAME, MapOutVar_names(jj)(:Nc_vars(jj)), &
     &         Nhru, 'double', Map_var_dble(1, jj))/=0 ) &
     &         CALL read_error(4, MapOutVar_names(jj)(:Nc_vars(jj)))
        ENDIF
      ENDDO

      DO jj = 1, NmapOutVars
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          IF ( Map_var_type(jj)==REAL_TYPE ) THEN
            map_var_double = DBLE( Map_var(i, jj) )
            IF ( Dailyresults==ACTIVE ) THEN
              Map_var_daily(i, jj) = map_var_double
            ELSEIF ( Weekresults==ACTIVE ) THEN
              Map_var_week(i, jj) = Map_var_week(i, jj) + map_var_double
            ELSE
              IF ( Totresults==ACTIVE ) Map_var_tot(i, jj) = Map_var_tot(i, jj) + map_var_double
              IF ( Yrresults==ACTIVE ) Map_var_yr(i, jj) = Map_var_yr(i, jj) + map_var_double
              IF ( Monresults==ACTIVE ) Map_var_mon(i, jj) = Map_var_mon(i, jj) + map_var_double
            ENDIF
          ELSEIF ( Map_var_type(jj)==DBLE_TYPE ) THEN
            IF ( Dailyresults==ACTIVE ) THEN
              Map_var_daily(i, jj) = Map_var_dble(i, jj)
            ELSEIF ( Weekresults==ACTIVE ) THEN
              Map_var_week(i, jj) = Map_var_week(i, jj) + Map_var_dble(i, jj)
            ELSE
              IF ( Totresults==ACTIVE ) Map_var_tot(i, jj) = Map_var_tot(i, jj) + Map_var_dble(i, jj)
              IF ( Yrresults==ACTIVE ) Map_var_yr(i, jj) = Map_var_yr(i, jj) + Map_var_dble(i, jj)
              IF ( Monresults==ACTIVE ) Map_var_mon(i, jj) = Map_var_mon(i, jj) + Map_var_dble(i, jj)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF ( Dailyresults==ACTIVE ) THEN
        Basin_var_daily = 0.0D0
        DO jj = 1, NmapOutVars
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            Basin_var_daily(jj) = Basin_var_daily(jj) + Map_var_daily(i, jj)*Hru_area_dble(i)
          ENDDO
          Basin_var_daily(jj) = Basin_var_daily(jj)*Basin_area_inv

          WRITE ( Dailyunit(jj), 9002 ) Nowyear, Nowmonth, Nowday, ' Basin daily mean:', Basin_var_daily(jj)
          IF ( Mapflg==1 ) THEN
            CALL write_results(Dailyunit(jj), Map_var_daily(1,jj))
          ELSE
            Map_var_id = 0.0D0
            DO k = 1, Nhrucell
              IF ( Gvr_map_id(k)==0 ) CYCLE
              IF ( Gvr_map_frac(k)<NEARZERO ) CYCLE
              Map_var_id(Gvr_map_id(k)) = Map_var_id(Gvr_map_id(k)) + &
     &                                    Map_var_daily(Gvr_hru_id(k), jj)*Gvr_map_frac_dble(k)
            ENDDO
            CALL write_results(Dailyunit(jj), Map_var_id)
          ENDIF
          WRITE ( Dailyunit(jj), 9003 )
        ENDDO
      ENDIF

      IF ( Weekresults==ACTIVE ) THEN
        Weekdays = Weekdays + 1
! check for seventh day
        IF ( Weekdays==7 ) THEN
          factor = Conv_fac/7.0D0
          Basin_var_week = 0.0D0
          DO jj = 1, NmapOutVars
            DO j = 1, Active_hrus
              i = Hru_route_order(j)
              Map_var_week(i, jj) = Map_var_week(i, jj)*factor
              Basin_var_week(jj) = Basin_var_week(jj) + Map_var_week(i, jj)*Hru_area_dble(i)
            ENDDO
            Basin_var_week(jj) = Basin_var_week(jj)*Basin_area_inv

            WRITE ( Weekunit(jj), 9002 ) Nowyear, Nowmonth, Nowday, ' Basin weekly mean:', Basin_var_week(jj)
            IF ( Mapflg==ACTIVE ) THEN
              CALL write_results(Weekunit(jj), Map_var_week(1,jj))
            ELSE
              Map_var_id = 0.0D0
              DO k = 1, Nhrucell
                IF ( Gvr_map_id(k)==0 ) CYCLE
                IF ( Gvr_map_frac(k)<NEARZERO ) CYCLE
                Map_var_id(Gvr_map_id(k)) = Map_var_id(Gvr_map_id(k)) + &
     &                                      Map_var_week(Gvr_hru_id(k),jj)*Gvr_map_frac_dble(k)
              ENDDO
              CALL write_results(Weekunit(jj), Map_var_id)
            ENDIF
            WRITE ( Weekunit(jj), 9003 )
          ENDDO
          Map_var_week = 0.0D0
          Weekdays = 0
        ENDIF
      ENDIF

      IF ( Monresults==ACTIVE ) THEN
        Mondays = Mondays + 1
! check for last day of current month
        IF ( Nowday==Modays(Nowmonth) .OR. last_day==1 ) THEN
          factor = Conv_fac/DBLE(Mondays)
          Basin_var_mon = 0.0D0
          DO jj = 1, NmapOutVars
            DO j = 1, Active_hrus
              i = Hru_route_order(j)
              Map_var_mon(i, jj) = Map_var_mon(i, jj)*factor
              Basin_var_mon(jj) = Basin_var_mon(jj) + Map_var_mon(i, jj)*Hru_area_dble(i)
            ENDDO
            Basin_var_mon(jj) = Basin_var_mon(jj)*Basin_area_inv

            WRITE ( Monunit(jj), 9002 ) Nowyear, Nowmonth, Nowday, &
     &                                  ' Basin monthly mean:', Basin_var_mon(jj)
            IF ( Mapflg==ACTIVE ) THEN
              CALL write_results(Monunit(jj), Map_var_mon(1,jj))
            ELSE
              Map_var_id = 0.0D0
              DO k = 1, Nhrucell
                IF ( Gvr_map_id(k)==0 ) CYCLE
                IF ( Gvr_map_frac(k)<NEARZERO ) CYCLE
                Map_var_id(Gvr_map_id(k)) = Map_var_id(Gvr_map_id(k)) + &
     &                                      Map_var_mon(Gvr_hru_id(k), jj)*Gvr_map_frac_dble(k)
              ENDDO
              CALL write_results(Monunit(jj), Map_var_id)
            ENDIF
            WRITE ( Monunit(jj), 9003 )
          ENDDO
          Map_var_mon = 0.0D0
          Mondays = 0
        ENDIF
      ENDIF

      IF ( Totresults==ACTIVE ) THEN
        Totdays = Totdays + 1
! check for last day of simulation
        IF ( last_day==1 ) THEN
          factor = Conv_fac/DBLE(Totdays)
          Basin_var_tot = 0.0D0
          DO jj = 1, NmapOutVars
            DO j = 1, Active_hrus
              i = Hru_route_order(j)
              Map_var_tot(i, jj) = Map_var_tot(i, jj)*factor
              Basin_var_tot(jj) = Basin_var_tot(jj) + Map_var_tot(i, jj)*Hru_area_dble(i)
            ENDDO
            Basin_var_tot(jj) = Basin_var_tot(jj)*Basin_area_inv

            WRITE ( Totunit(jj), 9004 ) 'Time period: ', Begyr, Start_month, Start_day, &
     &              Nowyear, Nowmonth, Nowday, ' Basin simulation mean:', Basin_var_tot(jj)

            IF ( Mapflg==ACTIVE ) THEN
              CALL write_results(Totunit(jj), Map_var_tot(1,jj))
            ELSE
              Map_var_id = 0.0D0
              DO k = 1, Nhrucell
                IF ( Gvr_map_id(k)==0 ) CYCLE
                IF ( Gvr_map_frac(k)<NEARZERO ) CYCLE
                Map_var_id(Gvr_map_id(k)) = Map_var_id(Gvr_map_id(k)) + &
     &                                      Map_var_tot(Gvr_hru_id(k), jj)*Gvr_map_frac_dble(k)
              ENDDO
              CALL write_results(Totunit(jj), Map_var_id)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

 9002 FORMAT ( I5, 2('/',I2.2), A, E12.3 )
 9003 FORMAT ( '######################', / )
 9004 FORMAT ( A, 2(I5, 2('/',I2.2)), A, E12.4 )

      END FUNCTION map_resultsrun

!***********************************************************************
!     map_resultsclean - close files
!***********************************************************************
      INTEGER FUNCTION map_resultsclean()
      USE PRMS_MAP_RESULTS
      IMPLICIT NONE
      INTEGER :: jj
!***********************************************************************
      map_resultsclean = 0

      IF ( Totresults==ACTIVE ) THEN
        DO jj = 1, NmapOutVars
          CLOSE ( Totunit(jj) )
        ENDDO
      ENDIF
      IF ( Monresults==ACTIVE ) THEN
        DO jj = 1, NmapOutVars
          CLOSE ( Monunit(jj) )
        ENDDO
      ENDIF
      IF ( Yrresults==ACTIVE ) THEN
        DO jj = 1, NmapOutVars
          CLOSE ( Yrunit(jj) )
        ENDDO
      ENDIF
      IF ( Weekresults==ACTIVE ) THEN
        DO jj = 1, NmapOutVars
          CLOSE ( Weekunit(jj) )
        ENDDO
      ENDIF

      END FUNCTION map_resultsclean

!***********************************************************************
!     write_results - write results to map
!***********************************************************************
      SUBROUTINE write_results(Iunit, Values)
      USE PRMS_MAP_RESULTS, ONLY: Numvalues, Mapfmt
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit
      DOUBLE PRECISION, INTENT(IN) :: Values(Numvalues)
! Local Variables
      INTEGER :: j
!***********************************************************************
      WRITE ( Iunit, Mapfmt ) (Values(j), j=1,Numvalues)
      END SUBROUTINE write_results
