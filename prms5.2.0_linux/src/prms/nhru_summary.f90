!***********************************************************************
!     Output a set of declared variables by HRU in CSV format
!***********************************************************************
      MODULE PRMS_NHRU_SUMMARY
      USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH, ERROR_control, ERROR_open_out, &
     &    DAILY, MONTHLY, DAILY_MONTHLY, MEAN_MONTHLY, MEAN_YEARLY, YEARLY, ACTIVE, OFF, &
     &    REAL_TYPE, DBLE_TYPE, INT_TYPE, DNEARZERO, RUN, DECL, INIT, CLEAN, DOCUMENTATION
      USE PRMS_MODULE, ONLY: Process_flag, Model, Nhru, Nsub, NhruOutON_OFF, Prms_warmup, &
     &    Start_year, Start_month, Start_day, End_year, End_month, End_day
      IMPLICIT NONE
! Module Variables
      character(len=*), parameter :: MODDESC = 'Output Summary'
      character(len=*), parameter :: MODNAME = 'nhru_summary'
      character(len=*), parameter :: Version_nhru_summary = '2020-12-02'
      INTEGER, SAVE :: Begin_results, Begyr, Lastyear
      INTEGER, SAVE, ALLOCATABLE :: Dailyunit(:), Nc_vars(:), Nhru_var_type(:), Nhru_var_int(:, :)
      REAL, SAVE, ALLOCATABLE :: Nhru_var_daily(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Nhru_var_dble(:, :)
      CHARACTER(LEN=48), SAVE :: Output_fmt, Output_fmt2, Output_fmt3, Output_fmtint
      CHARACTER(LEN=48), SAVE :: Output_grid_fmt, Output_grid_fmtint, Output_date_fmt, Output_date_fmt3, Output_fmt3int
      INTEGER, SAVE :: Daily_flag, Double_vars, Yeardays, Monthly_flag, Integer_vars
      INTEGER, SAVE :: dates_next_year, dates_next_month, dates_next_day, selectDates_unit
      DOUBLE PRECISION, SAVE :: Monthdays
      INTEGER, SAVE, ALLOCATABLE :: Monthlyunit(:), Yearlyunit(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Nhru_var_monthly(:, :), Nhru_var_yearly(:, :)
! Paramters
      INTEGER, SAVE, ALLOCATABLE :: Nhm_id(:)
! Control Parameters
      INTEGER, SAVE :: NhruOutVars, NhruOut_freq, NhruOut_format, NhruOutNcol, outputSelectDatesON_OFF 
      CHARACTER(LEN=36), SAVE, ALLOCATABLE :: NhruOutVar_names(:)
      CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: NhruOutBaseFileName, selectDatesFileName
      END MODULE PRMS_NHRU_SUMMARY

!     ******************************************************************
!     nhru results module
!     ******************************************************************
      SUBROUTINE nhru_summary()
      USE PRMS_NHRU_SUMMARY
      IMPLICIT NONE
! Functions
      EXTERNAL :: nhru_summarydecl, nhru_summaryinit, nhru_summaryrun
! Local Variables
      INTEGER :: i
!***********************************************************************
      IF ( Process_flag==RUN ) THEN
        CALL nhru_summaryrun()
      ELSEIF ( Process_flag==DECL ) THEN
        CALL nhru_summarydecl()
      ELSEIF ( Process_flag==INIT ) THEN
        CALL nhru_summaryinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        DO i = 1, NhruOutVars
          IF ( Daily_flag==ACTIVE ) THEN
            IF ( Dailyunit(i)>0 ) CLOSE ( Dailyunit(i) )
          ENDIF
          IF ( NhruOut_freq>MEAN_MONTHLY ) THEN
            IF ( Yearlyunit(i)>0 ) CLOSE ( Yearlyunit(i) )
          ENDIF
          IF ( Monthly_flag==ACTIVE ) THEN
            IF ( Monthlyunit(i)>0 ) CLOSE ( Monthlyunit(i) )
          ENDIF
        ENDDO
      ENDIF

      END SUBROUTINE nhru_summary

!***********************************************************************
!     declare parameters and variables
!***********************************************************************
      SUBROUTINE nhru_summarydecl()
      USE PRMS_NHRU_SUMMARY
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: control_string_array, control_integer, control_string, declparam
      EXTERNAL :: read_error, print_module, error_stop
! Local Variables
      INTEGER :: i
!***********************************************************************
      CALL print_module(MODDESC, MODNAME, Version_nhru_summary)

      IF ( control_integer(NhruOutVars, 'nhruOutVars')/=0 ) NhruOutVars = 0
      ! 1 = daily, 2 = monthly, 3 = both, 4 = mean monthly, 5 = mean yearly, 6 = yearly total
      IF ( control_integer(NhruOut_freq, 'nhruOut_freq')/=0 ) NhruOut_freq = 0
      IF ( NhruOut_freq<DAILY .OR. NhruOut_freq>YEARLY ) CALL error_stop('invalid nhruOut_freq value', ERROR_control)
      ! 1 = ES10.3; 2 = F0.2; 3 = F0.3; 4 = F0.4; 5 = F0.5
      IF ( control_integer(NhruOut_format, 'nhruOut_format')/=0 ) NhruOut_format = 1
      IF ( NhruOut_format<1 .OR. NhruOut_format>5 ) CALL error_stop('invalid nhruOut_format value', ERROR_control)
      IF ( control_integer(NhruOutNcol, 'nhruOutNcol')/=0 ) NhruOutNcol = 0

      IF ( control_integer(outputSelectDatesON_OFF, 'outputSelectDatesON_OFF')/=0 ) outputSelectDatesON_OFF = OFF
      IF ( outputSelectDatesON_OFF==ACTIVE ) THEN
        IF ( control_string(selectDatesFileName, 'selectDatesFileName')/=0 ) CALL read_error(5, 'selectDatesFileName')
      ENDIF

      IF ( NhruOutVars==0 ) THEN
        IF ( Model/=DOCUMENTATION ) CALL error_stop('nhru_summary requested with nhruOutVars equal 0', ERROR_control)
      ELSE
        ALLOCATE ( NhruOutVar_names(NhruOutVars), Nhru_var_type(NhruOutVars), Nc_vars(NhruOutVars) )
        NhruOutVar_names = ' '
        DO i = 1, NhruOutVars
          IF ( control_string_array(NhruOutVar_names(i), 'nhruOutVar_names', i)/=0 ) CALL read_error(5, 'nhruOutVar_names')
        ENDDO
        IF ( control_string(NhruOutBaseFileName, 'nhruOutBaseFileName')/=0 ) CALL read_error(5, 'nhruOutBaseFileName')
      ENDIF

! Declared Parameters
      IF ( NhruOutON_OFF==2 .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Nhm_id(Nhru) )
        IF ( declparam(MODNAME, 'nhm_id', 'nhru', 'integer', &
     &       '1', '1', '9999999', &
     &       'National Hydrologic Model HRU ID', 'National Hydrologic Model HRU ID', &
     &       'none') /= 0 ) CALL read_error(1, 'nhm_id')
      ENDIF

      END SUBROUTINE nhru_summarydecl

!***********************************************************************
!     Initialize module values
!***********************************************************************
      SUBROUTINE nhru_summaryinit()
      USE PRMS_NHRU_SUMMARY
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getvartype, numchars, getvarsize, getparam
      EXTERNAL read_error, PRMS_open_output_file, find_header_end, find_current_file_time
! Local Variables
      INTEGER :: ios, ierr, size, dim, jj, j
      CHARACTER(LEN=MAXFILE_LENGTH) :: fileName
!***********************************************************************
      Begin_results = ACTIVE
      IF ( Prms_warmup>0 ) Begin_results = OFF
      Begyr = Start_year + Prms_warmup
      Lastyear = Begyr

      IF ( outputSelectDatesON_OFF==ACTIVE ) THEN
        CALL find_header_end(selectDates_unit, selectDatesFileName, 'selectDatesFileName', ierr, 0, 0)
        IF ( ierr==0 ) THEN
          CALL find_current_file_time(selectDates_unit, Start_year, Start_month, Start_day, &
     &                                dates_next_year, dates_next_month, dates_next_day)
        ELSE
          ierr = 1
        ENDIF
      ENDIF

      IF ( NhruOutNcol<1 ) NhruOutNcol = Nhru

      IF ( NhruOut_format==1 ) THEN
        WRITE ( Output_fmt, 9001 ) Nhru
        WRITE ( Output_grid_fmt, 9014 ) NhruOutNcol - 1
      ELSEIF ( NhruOut_format==2 ) THEN
        WRITE ( Output_fmt, 9007 ) Nhru
        WRITE ( Output_grid_fmt, 9017 ) NhruOutNcol - 1
      ELSEIF ( NhruOut_format==3 ) THEN
        WRITE ( Output_fmt, 9006 ) Nhru
        WRITE ( Output_grid_fmt, 9016 ) NhruOutNcol - 1
      ELSEIF ( NhruOut_format==4 ) THEN
        WRITE ( Output_fmt, 9005 ) Nhru
        WRITE ( Output_grid_fmt, 9015 ) NhruOutNcol - 1
      ELSEIF ( NhruOut_format==5 ) THEN
        WRITE ( Output_fmt, 9012 ) Nhru
        WRITE ( Output_grid_fmt, 9014 ) NhruOutNcol - 1
      ENDIF
      WRITE ( Output_fmtint, 9004 ) Nhru
      WRITE ( Output_grid_fmtint, 9018 ) NhruOutNcol - 1
      WRITE ( Output_date_fmt, 9013 )

      Double_vars = OFF
      Integer_vars = OFF
      ierr = 0
      DO jj = 1, NhruOutVars
        Nc_vars(jj) = numchars(NhruOutVar_names(jj))
        Nhru_var_type(jj) = getvartype(NhruOutVar_names(jj)(:Nc_vars(jj)), Nhru_var_type(jj) )
        IF ( Nhru_var_type(jj)==DBLE_TYPE ) Double_vars = ACTIVE
        IF ( Nhru_var_type(jj)==INT_TYPE ) Integer_vars = ACTIVE
        IF ( Nhru_var_type(jj)/=REAL_TYPE .AND. Nhru_var_type(jj)/=DBLE_TYPE .AND. Nhru_var_type(jj)/=INT_TYPE ) THEN
          PRINT *, 'ERROR, invalid nhru_summary variable:', NhruOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only integer, real or double variables allowed'
          ierr = 1
        ENDIF
        size = getvarsize(NhruOutVar_names(jj)(:Nc_vars(jj)), dim )
        IF ( size/=Nhru ) THEN
          PRINT *, 'ERROR, invalid nhru_summary variable:', NhruOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only variables dimensioned by nhru, nssr, or ngw allowed'
          ierr = 1
        ENDIF
      ENDDO
      IF ( ierr==1 ) ERROR STOP ERROR_control
      IF ( Double_vars==ACTIVE ) THEN
        ALLOCATE ( Nhru_var_dble(Nhru, NhruOutVars) )
        Nhru_var_dble = 0.0D0
      ENDIF
      IF ( Integer_vars==ACTIVE ) THEN
        ALLOCATE ( Nhru_var_int(Nhru, NhruOutVars) )
        Nhru_var_int = 0
      ENDIF

      Daily_flag = OFF
      IF ( NhruOut_freq==DAILY .OR. NhruOut_freq==DAILY_MONTHLY ) THEN
        Daily_flag = ACTIVE
        ALLOCATE ( Dailyunit(NhruOutVars) )
        Dailyunit = 0
      ENDIF

      Monthly_flag = OFF
      IF ( NhruOut_freq==MONTHLY .OR. NhruOut_freq==DAILY_MONTHLY .OR. NhruOut_freq==MEAN_MONTHLY ) Monthly_flag = ACTIVE

      IF ( NhruOut_freq>MEAN_MONTHLY ) THEN
        Yeardays = 0
        ALLOCATE ( Nhru_var_yearly(Nhru, NhruOutVars), Yearlyunit(NhruOutVars) )
        Nhru_var_yearly = 0.0D0
        Yearlyunit = 0
        IF ( NhruOut_format==1 ) THEN
          WRITE ( Output_fmt3, 9003 ) Nhru
        ELSEIF ( NhruOut_format==2 ) THEN
          WRITE ( Output_fmt3, 9010 ) Nhru
        ELSEIF ( NhruOut_format==3 ) THEN
          WRITE ( Output_fmt3, 9009 ) Nhru
        ELSEIF ( NhruOut_format==4 ) THEN
          WRITE ( Output_fmt3, 9008 ) Nhru
        ELSEIF ( NhruOut_format==5 ) THEN
          WRITE ( Output_fmt3, 9011 ) Nhru
        ENDIF
        WRITE ( Output_date_fmt3, "(A)" ) '(I0)'
        WRITE ( Output_fmt3int, 9019 ) NhruOutNcol - 1
      ELSEIF ( Monthly_flag==ACTIVE ) THEN
        Monthdays = 0.0D0
        ALLOCATE ( Nhru_var_monthly(Nhru, NhruOutVars), Monthlyunit(NhruOutVars) )
        Nhru_var_monthly = 0.0D0
        Monthlyunit = 0
      ENDIF

      IF ( NhruOutON_OFF==2 ) THEN
        IF ( getparam(MODNAME, 'nhm_id', Nhru, 'integer', Nhm_id)/=0 ) CALL read_error(2, 'nhm_id')
      ENDIF
      WRITE ( Output_fmt2, 9002 ) Nhru
      ALLOCATE ( Nhru_var_daily(Nhru, NhruOutVars) )
      Nhru_var_daily = 0.0
      DO jj = 1, NhruOutVars
        IF ( Daily_flag==ACTIVE ) THEN
          fileName = NhruOutBaseFileName(:numchars(NhruOutBaseFileName))//NhruOutVar_names(jj)(:Nc_vars(jj))//'.csv'
          CALL PRMS_open_output_file(Dailyunit(jj), fileName, 'xxx', 0, ios)
          IF ( ios/=0 ) CALL error_stop('in nhru_summary, daily', ERROR_open_out)
          IF ( NhruOutON_OFF==1 ) THEN
            WRITE ( Dailyunit(jj), Output_fmt2 ) (j, j=1,Nhru)
          ELSE
            WRITE ( Dailyunit(jj), Output_fmt2 ) (Nhm_id(j), j=1,Nhru)
          ENDIF
        ENDIF
        IF ( NhruOut_freq>MEAN_MONTHLY ) THEN
          IF ( NhruOut_freq==MEAN_YEARLY ) THEN
            fileName = NhruOutBaseFileName(:numchars(NhruOutBaseFileName))//NhruOutVar_names(jj)(:Nc_vars(jj))//'_meanyearly.csv'
            CALL PRMS_open_output_file(Yearlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nhru_summary, mean yearly', ERROR_open_out)
          ELSE  !IF ( NhruOut_freq==YEARLY ) THEN
            fileName = NhruOutBaseFileName(:numchars(NhruOutBaseFileName))//NhruOutVar_names(jj)(:Nc_vars(jj))//'_yearly.csv'
            CALL PRMS_open_output_file(Yearlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nhru_summary, yearly', ERROR_open_out)
          ENDIF
          IF ( NhruOutON_OFF==1 ) THEN
            WRITE ( Yearlyunit(jj), Output_fmt2 ) (j, j=1,Nhru)
          ELSE
            WRITE ( Yearlyunit(jj), Output_fmt2 ) (Nhm_id(j), j=1,Nhru)
          ENDIF
        ENDIF
        IF ( Monthly_flag==ACTIVE ) THEN
          IF ( NhruOut_freq==MEAN_MONTHLY ) THEN
            fileName = NhruOutBaseFileName(:numchars(NhruOutBaseFileName))//NhruOutVar_names(jj)(:Nc_vars(jj))// &
     &                 '_meanmonthly.csv'
            CALL PRMS_open_output_file(Monthlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nhru_summary, mean monthly', ERROR_open_out)
          ELSE
            fileName = NhruOutBaseFileName(:numchars(NhruOutBaseFileName))//NhruOutVar_names(jj)(:Nc_vars(jj))//'_monthly.csv'
            CALL PRMS_open_output_file(Monthlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nhru_summary, monthly', ERROR_open_out)
          ENDIF
          IF ( NhruOutON_OFF==1 ) THEN
            WRITE ( Monthlyunit(jj), Output_fmt2 ) (j, j=1,Nhru)
          ELSE
            WRITE ( Monthlyunit(jj), Output_fmt2 ) (Nhm_id(j), j=1,Nhru)
          ENDIF
        ENDIF
      ENDDO

 9001 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',ES10.3))')
 9002 FORMAT ('("Date"',I0,'('', ''I0))')
 9003 FORMAT ('(I4,', I0,'('','',ES10.3))')
 9004 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',I0))')
 9005 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.4))')
 9006 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.3))')
 9007 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.2))')
 9008 FORMAT ('(I4,', I0,'('','',F0.4))')
 9009 FORMAT ('(I4,', I0,'('','',F0.3))')
 9010 FORMAT ('(I4,', I0,'('','',F0.2))')
 9011 FORMAT ('(I4,', I0,'('','',F0.5))')
 9012 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.5))')
 9013 FORMAT ('(I4, 2(''-'',I2.2))')
 9014 FORMAT ('(',I0,'(ES10.3,'',''),ES10.3)')
 9015 FORMAT ('(',I0,'(F0.4,'',''),F0.4)')
 9016 FORMAT ('(',I0,'(F0.3,'',''),F0.3)')
 9017 FORMAT ('(',I0,'(F0.2,'',''),F0.2)')
 9018 FORMAT ('(',I0,'(I0,'',''),I0)')
 9019 FORMAT ('(I4,', I0,'('','',I0),I0)')

      END SUBROUTINE nhru_summaryinit

!***********************************************************************
!     Output set of declared variables in CSV format
!***********************************************************************
      SUBROUTINE nhru_summaryrun()
      USE PRMS_NHRU_SUMMARY
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order
      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Nowday, Modays
      IMPLICIT NONE
! FUNCTIONS AND SUBROUTINES
      INTRINSIC :: SNGL, DBLE
      INTEGER, EXTERNAL :: getvar
      EXTERNAL :: read_error
! Local Variables
      INTEGER :: j, i, jj, write_month, last_day, write_date
!***********************************************************************
      IF ( Begin_results==OFF ) THEN
        IF ( Nowyear==Begyr .AND. Nowmonth==Start_month .AND. Nowday==Start_day ) THEN
          Begin_results = ACTIVE
        ELSE
          RETURN
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
! need getvars for each variable (only can have short string)
      DO jj = 1, NhruOutVars
        IF ( Nhru_var_type(jj)==REAL_TYPE ) THEN
          IF ( getvar(MODNAME, NhruOutVar_names(jj)(:Nc_vars(jj)), Nhru, 'real', Nhru_var_daily(1, jj))/=0 ) &
     &         CALL read_error(4, NhruOutVar_names(jj)(:Nc_vars(jj)))
        ELSEIF ( Nhru_var_type(jj)==DBLE_TYPE ) THEN
          IF ( getvar(MODNAME, NhruOutVar_names(jj)(:Nc_vars(jj)), Nhru, 'double', Nhru_var_dble(1, jj))/=0 ) &
     &         CALL read_error(4, NhruOutVar_names(jj)(:Nc_vars(jj)))
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            Nhru_var_daily(i, jj) = SNGL( Nhru_var_dble(i, jj) )
          ENDDO
        ELSEIF ( Nhru_var_type(jj)==INT_TYPE ) THEN
          IF ( getvar(MODNAME, NhruOutVar_names(jj)(:Nc_vars(jj)), Nhru, 'integer', Nhru_var_int(1, jj))/=0 ) &
     &         CALL read_error(4, NhruOutVar_names(jj)(:Nc_vars(jj)))
          IF ( NhruOut_freq>DAILY ) THEN
            DO j = 1, Active_hrus
              i = Hru_route_order(j)
              Nhru_var_daily(i, jj) = FLOAT( Nhru_var_int(i, jj) )
            ENDDO
          ENDIF
        ENDIF
        IF ( Daily_flag==ACTIVE ) THEN
          write_date = 1
          IF ( outputSelectDatesON_OFF==ACTIVE ) THEN
            write_date = 0
            IF ( Nowyear==dates_next_year .AND. Nowmonth==dates_next_month.AND. Nowday==dates_next_day ) write_date = 1
          ENDIF
          IF ( write_date==1 ) THEN
            IF ( Nhru_var_type(jj)/=INT_TYPE ) THEN
              IF ( NhruOutNcol==Nhru) THEN
                WRITE ( Dailyunit(jj), Output_fmt) Nowyear, Nowmonth, Nowday, (Nhru_var_daily(j,jj), j=1,Nhru)
              ELSE
                WRITE ( Dailyunit(jj), Output_date_fmt) Nowyear, Nowmonth, Nowday
                WRITE ( Dailyunit(jj), Output_grid_fmt) (Nhru_var_daily(j,jj), j=1,Nhru)
              ENDIF
            ELSE
              IF ( NhruOutNcol==Nhru) THEN
                WRITE ( Dailyunit(jj), Output_fmtint) Nowyear, Nowmonth, Nowday, (Nhru_var_int(j,jj), j=1,Nhru)
              ELSE
                WRITE ( Dailyunit(jj), Output_date_fmt) Nowyear, Nowmonth, Nowday
                WRITE ( Dailyunit(jj), Output_grid_fmtint) (Nhru_var_int(j,jj), j=1,Nhru)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( outputSelectDatesON_OFF==ACTIVE ) THEN
        IF ( Daily_flag==ACTIVE ) THEN
          IF ( write_date==1 ) CALL read_event_date(selectDates_unit, dates_next_year, dates_next_month, dates_next_day)
        ENDIF
      ENDIF

      write_month = OFF
      IF ( NhruOut_freq>MEAN_MONTHLY ) THEN
        last_day = OFF
        IF ( Nowyear==End_year .AND. Nowmonth==End_month .AND. Nowday==End_day ) last_day = ACTIVE
        IF ( Lastyear/=Nowyear .OR. last_day==ACTIVE ) THEN
          IF ( (Nowmonth==Start_month .AND. Nowday==Start_day) .OR. last_day==ACTIVE ) THEN
            DO jj = 1, NhruOutVars
              IF ( NhruOut_freq==MEAN_YEARLY ) THEN
                DO j = 1, Active_hrus
                  i = Hru_route_order(j)
                  Nhru_var_yearly(i, jj) = Nhru_var_yearly(i, jj)/Yeardays
                ENDDO
              ENDIF
              IF ( Nhru_var_type(jj)/=INT_TYPE ) THEN
                IF ( NhruOutNcol==Nhru) THEN
                  WRITE ( Yearlyunit(jj), Output_fmt3) Lastyear, (Nhru_var_yearly(j,jj), j=1,Nhru)
                ELSE
                  WRITE ( Yearlyunit(jj), Output_date_fmt3) Lastyear
                  WRITE ( Yearlyunit(jj), Output_grid_fmt) (Nhru_var_yearly(j,jj), j=1,Nhru)
                ENDIF
              ELSE
                DO i = 1, Nhru
                  Nhru_var_int(i, jj) = INT( Nhru_var_yearly(i, jj) )
                ENDDO
                IF ( NhruOutNcol==Nhru) THEN
                  WRITE ( Yearlyunit(jj), Output_fmt3int) Lastyear, (Nhru_var_int(j,jj), j=1,Nhru)
                ELSE
                  WRITE ( Yearlyunit(jj), Output_date_fmt3) Lastyear
                  WRITE ( Yearlyunit(jj), Output_grid_fmtint) (Nhru_var_int(j,jj), j=1,Nhru)
                ENDIF
              ENDIF
            ENDDO
            Nhru_var_yearly = 0.0D0
            Yeardays = 0
            Lastyear = Nowyear
          ENDIF
        ENDIF
        Yeardays = Yeardays + 1
      ELSEIF ( Monthly_flag==ACTIVE ) THEN
        ! check for last day of month and simulation
        IF ( Nowday==Modays(Nowmonth) ) THEN
          write_month = ACTIVE
        ELSEIF ( Nowyear==End_year ) THEN
          IF ( Nowmonth==End_month ) THEN
            IF ( Nowday==End_day ) write_month = ACTIVE
          ENDIF
        ENDIF
        Monthdays = Monthdays + 1.0D0
      ENDIF

      IF ( NhruOut_freq>MEAN_MONTHLY ) THEN
        DO jj = 1, NhruOutVars
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            Nhru_var_yearly(i, jj) = Nhru_var_yearly(i, jj) + DBLE( Nhru_var_daily(i, jj) )
          ENDDO
        ENDDO
        RETURN
      ENDIF

      IF ( Monthly_flag==ACTIVE ) THEN
        DO jj = 1, NhruOutVars
          DO j = 1, Active_hrus
            i = Hru_route_order(j)
            Nhru_var_monthly(i, jj) = Nhru_var_monthly(i, jj) + DBLE( Nhru_var_daily(i, jj) )
            IF ( write_month==ACTIVE ) THEN
              IF ( NhruOut_freq==MEAN_MONTHLY ) Nhru_var_monthly(i, jj) = Nhru_var_monthly(i, jj)/Monthdays
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF ( write_month==ACTIVE ) THEN
        DO jj = 1, NhruOutVars
          IF ( Nhru_var_type(jj)/=INT_TYPE ) THEN
            IF ( NhruOutNcol==Nhru) THEN
              WRITE ( Monthlyunit(jj), Output_fmt) Nowyear, Nowmonth, Nowday, (Nhru_var_monthly(j,jj), j=1,Nhru)
            ELSE
              WRITE ( Monthlyunit(jj), Output_date_fmt) Nowyear, Nowmonth, Nowday
              WRITE ( Monthlyunit(jj), Output_grid_fmt) (Nhru_var_monthly(j,jj), j=1,Nhru)
            ENDIF
          ELSE
            DO i = 1, Nhru
              Nhru_var_int(i, jj) = INT( Nhru_var_monthly(i, jj) )
            ENDDO
            IF ( NhruOutNcol==Nhru) THEN
              WRITE ( Monthlyunit(jj), Output_fmtint) Nowyear, Nowmonth, Nowday, (Nhru_var_int(j,jj), j=1,Nhru)
            ELSE
              WRITE ( Monthlyunit(jj), Output_date_fmt) Nowyear, Nowmonth, Nowday
              WRITE ( Monthlyunit(jj), Output_grid_fmtint) (Nhru_var_int(j,jj), j=1,Nhru)
            ENDIF
          ENDIF
        ENDDO
        Monthdays = 0.0D0
        Nhru_var_monthly = 0.0D0
      ENDIF

      END SUBROUTINE nhru_summaryrun

!*****************************
! Read event for a source type
!*****************************
      SUBROUTINE read_event_date(Iunit, Next_yr, Next_mo, Next_day)
      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Nowday
      USE PRMS_CONSTANTS, ONLY: ERROR_water_use, ACTIVE, OFF
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit
      INTEGER, INTENT (INOUT) :: Next_yr, Next_mo, Next_day
! Funcions
      EXTERNAL :: check_event, set_transfers, is_eof
! Local Variables
      INTEGER keep_reading
!*******************************************************************************
      IF ( Next_mo==0 ) RETURN ! already found end of file
      keep_reading = ACTIVE
      DO WHILE ( keep_reading==ACTIVE )
        IF ( Next_yr==Nowyear .AND. Next_mo==Nowmonth .AND. Next_day==Nowday ) THEN
          READ ( Iunit, * ) Next_yr, Next_mo, Next_day
          CALL is_eof(Iunit, Next_yr, Next_mo, Next_day)
          IF ( Next_mo==0 ) keep_reading = OFF
        ELSE
          keep_reading = OFF
        ENDIF
      ENDDO
      END SUBROUTINE read_event_date
