!***********************************************************************
!     Output a set of declared variables by segment in CSV format
!***********************************************************************
      MODULE PRMS_NSEGMENT_SUMMARY
      USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH, ERROR_control, ERROR_open_out, &
     &    DAILY, MONTHLY, DAILY_MONTHLY, MEAN_MONTHLY, MEAN_YEARLY, YEARLY, &
     &    ACTIVE, OFF, REAL_TYPE, DBLE_TYPE, DOCUMENTATION
      USE PRMS_MODULE, ONLY: Model, Nsegment, Inputerror_flag, NsegmentOutON_OFF, &
     &    Start_year, Start_month, Start_day, End_year, End_month, End_day, Prms_warmup
      IMPLICIT NONE
! Module Variables
      character(len=*), parameter :: MODDESC = 'Output Summary'
      character(len=*), parameter :: MODNAME = 'nsegment_summary'
      character(len=*), parameter :: Version_nsegment_summary = '2020-12-02'
      INTEGER, SAVE :: Begin_results, Begyr, Lastyear
      INTEGER, SAVE, ALLOCATABLE :: Dailyunit(:), Nc_vars(:), Nsegment_var_type(:)
      REAL, SAVE, ALLOCATABLE :: Nsegment_var_daily(:, :)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Nsegment_var_dble(:, :)
      CHARACTER(LEN=48), SAVE :: Output_fmt, Output_fmt2, Output_fmt3
      INTEGER, SAVE :: Daily_flag, Double_vars, Yeardays, Monthly_flag
      DOUBLE PRECISION, SAVE :: Monthdays
      INTEGER, SAVE, ALLOCATABLE :: Monthlyunit(:), Yearlyunit(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Nsegment_var_monthly(:, :), Nsegment_var_yearly(:, :)
! Paramters
      INTEGER, SAVE, ALLOCATABLE :: Nhm_seg(:)
! Control Parameters
      INTEGER, SAVE :: NsegmentOutVars, NsegmentOut_freq, NsegmentOut_format
      CHARACTER(LEN=36), SAVE, ALLOCATABLE :: NsegmentOutVar_names(:)
      CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: NsegmentOutBaseFileName
      END MODULE PRMS_NSEGMENT_SUMMARY

!     ******************************************************************
!     nsegment results module
!     ******************************************************************
      SUBROUTINE nsegment_summary()
      USE PRMS_CONSTANTS, ONLY: MEAN_MONTHLY, ACTIVE, RUN, DECL, INIT, CLEAN
      USE PRMS_MODULE, ONLY: Process_flag
      USE PRMS_NSEGMENT_SUMMARY
      IMPLICIT NONE
! Functions
      EXTERNAL :: nsegment_summarydecl, nsegment_summaryinit, nsegment_summaryrun
! Local Variables
      INTEGER :: i
!***********************************************************************
      IF ( Process_flag==RUN ) THEN
        CALL nsegment_summaryrun()
      ELSEIF ( Process_flag==DECL ) THEN
        CALL nsegment_summarydecl()
      ELSEIF ( Process_flag==INIT ) THEN
        CALL nsegment_summaryinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        DO i = 1, NsegmentOutVars
          IF ( Daily_flag==ACTIVE ) THEN
            IF ( Dailyunit(i)>0 ) CLOSE ( Dailyunit(i) )
          ENDIF
          IF ( NsegmentOut_freq>MEAN_MONTHLY ) THEN
            IF ( Yearlyunit(i)>0 ) CLOSE ( Yearlyunit(i) )
          ENDIF
          IF ( Monthly_flag==ACTIVE) THEN
            IF ( Monthlyunit(i)>0 ) CLOSE ( Monthlyunit(i) )
          ENDIF
        ENDDO
      ENDIF

      END SUBROUTINE nsegment_summary

!***********************************************************************
!     declare parameters and variables
!***********************************************************************
      SUBROUTINE nsegment_summarydecl()
      USE PRMS_NSEGMENT_SUMMARY
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: control_string_array, control_integer, control_string, declparam
      EXTERNAL :: read_error, print_module, error_stop
! Local Variables
      INTEGER :: i
!***********************************************************************
      CALL print_module(MODDESC, MODNAME, Version_nsegment_summary)

      IF ( control_integer(NsegmentOutVars, 'nsegmentOutVars')/=0 ) NsegmentOutVars = 0
      ! 1 = daily, 2 = monthly, 3 = both, 4 = mean monthly, 5 = mean yearly, 6 = yearly total
      IF ( control_integer(NsegmentOut_freq, 'nsegmentOut_freq')/=0 ) NsegmentOut_freq = 0
      IF ( NsegmentOut_freq<DAILY .OR. NsegmentOut_freq>YEARLY ) CALL error_stop('invalid nsegmentOut_freq value', ERROR_control)
      ! 1 = ES10.3; 2 = F0.2; 3 = F0.3; 4 = F0.4; 5 = F0.5
      IF ( control_integer(NsegmentOut_format, 'nsegmentOut_format')/=0 ) NsegmentOut_format = 1
      IF ( NsegmentOut_format<1 .OR. NsegmentOut_format>5 ) CALL error_stop('invalid nsegmentOut_format value', ERROR_control)

      IF ( NsegmentOutVars==0 ) THEN
        IF ( Model/=DOCUMENTATION ) CALL error_stop('nsegment_summary requested with nsegmentOutVars equal 0', ERROR_control)
      ELSE
        ALLOCATE ( NsegmentOutVar_names(NsegmentOutVars), Nsegment_var_type(NsegmentOutVars), Nc_vars(NsegmentOutVars) )
        NsegmentOutVar_names = ' '
        DO i = 1, NsegmentOutVars
          IF ( control_string_array(NsegmentOutVar_names(i), 'nsegmentOutVar_names', i)/=0 ) &
     &         CALL read_error(5, 'nsegmentOutVar_names')
        ENDDO
        IF ( control_string(NsegmentOutBaseFileName, 'nsegmentOutBaseFileName')/=0 ) CALL read_error(5, 'nsegmentOutBaseFileName')
      ENDIF

      IF ( NsegmentOutON_OFF==2 .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Nhm_seg(Nsegment) )
        IF (declparam(MODNAME, 'nhm_seg', 'nsegment', 'integer', &
     &       '0', '0', '9999999', &
     &       'National Hydrologic Model segment ID', 'National Hydrologic Model segment ID', &
     &       'none')/=0 ) CALL read_error(1, 'nhm_seg')
      ENDIF

      END SUBROUTINE nsegment_summarydecl

!***********************************************************************
!     Initialize module values
!***********************************************************************
      SUBROUTINE nsegment_summaryinit()
      USE PRMS_NSEGMENT_SUMMARY
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getvartype, numchars, getvarsize, getparam
      EXTERNAL :: read_error, PRMS_open_output_file, error_stop
! Local Variables
      INTEGER :: ios, ierr, size, jj, j
      CHARACTER(LEN=MAXFILE_LENGTH) :: fileName
!***********************************************************************
      Begin_results = ACTIVE
      IF ( Prms_warmup>0 ) Begin_results = OFF
      Begyr = Start_year + Prms_warmup
      Lastyear = Begyr

      IF ( NsegmentOut_format==1 ) THEN
        WRITE ( Output_fmt, 9001 ) Nsegment
      ELSEIF ( NsegmentOut_format==2 ) THEN
        WRITE ( Output_fmt, 9007 ) Nsegment
      ELSEIF ( NsegmentOut_format==3 ) THEN
        WRITE ( Output_fmt, 9006 ) Nsegment
      ELSEIF ( NsegmentOut_format==4 ) THEN
        WRITE ( Output_fmt, 9005 ) Nsegment
      ELSEIF ( NsegmentOut_format==5 ) THEN
        WRITE ( Output_fmt, 9012 ) Nsegment
      ENDIF

      Double_vars = OFF
      ierr = 0
      DO jj = 1, NsegmentOutVars
        Nc_vars(jj) = numchars(NsegmentOutVar_names(jj))
        Nsegment_var_type(jj) = getvartype(NsegmentOutVar_names(jj)(:Nc_vars(jj)) )
        IF ( Nsegment_var_type(jj)==DBLE_TYPE ) Double_vars = ACTIVE
        IF ( Nsegment_var_type(jj)/=REAL_TYPE .AND. Nsegment_var_type(jj)/=DBLE_TYPE ) THEN
          PRINT *, 'ERROR, invalid nsegment_summary variable:', NsegmentOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only real or double variables allowed'
          ierr = 1
        ENDIF
        size = getvarsize(NsegmentOutVar_names(jj)(:Nc_vars(jj)) )
        IF ( size/=Nsegment ) THEN
          PRINT *, 'ERROR, invalid nsegment_summary variable:', NsegmentOutVar_names(jj)(:Nc_vars(jj))
          PRINT *, '       only variables dimensioned by nsegment are allowed'
          ierr = 1
        ENDIF
      ENDDO
      IF ( ierr==1 ) ERROR STOP ERROR_control
      IF ( Double_vars==ACTIVE ) THEN
        ALLOCATE ( Nsegment_var_dble(Nsegment, NsegmentOutVars) )
        Nsegment_var_dble = 0.0D0
      ENDIF

      Daily_flag = OFF
      IF ( NsegmentOut_freq==DAILY .OR. NsegmentOut_freq==DAILY_MONTHLY ) THEN
        Daily_flag = ACTIVE
        ALLOCATE ( Dailyunit(NsegmentOutVars) )
        Dailyunit = 0
      ENDIF

      Monthly_flag = OFF
      IF ( NsegmentOut_freq==MONTHLY .OR. NsegmentOut_freq==DAILY_MONTHLY &
     &     .OR. NsegmentOut_freq==MEAN_MONTHLY ) Monthly_flag = ACTIVE

      IF ( NsegmentOut_freq>MEAN_MONTHLY ) THEN
        Yeardays = 0
        ALLOCATE ( Nsegment_var_yearly(Nsegment, NsegmentOutVars), &
     &             Yearlyunit(NsegmentOutVars) )
        Nsegment_var_yearly = 0.0D0
        Yearlyunit = 0
        IF ( NsegmentOut_format==1 ) THEN
          WRITE ( Output_fmt3, 9003 ) Nsegment
        ELSEIF ( NsegmentOut_format==2 ) THEN
          WRITE ( Output_fmt3, 9010 ) Nsegment
        ELSEIF ( NsegmentOut_format==3 ) THEN
          WRITE ( Output_fmt3, 9009 ) Nsegment
        ELSEIF ( NsegmentOut_format==4 ) THEN
          WRITE ( Output_fmt3, 9008 ) Nsegment
        ELSEIF ( NsegmentOut_format==5 ) THEN
          WRITE ( Output_fmt3, 9011 ) Nsegment
        ENDIF
      ENDIF
      IF ( Monthly_flag==ACTIVE ) THEN
        Monthdays = 0.0D0
        ALLOCATE ( Nsegment_var_monthly(Nsegment, NsegmentOutVars), Monthlyunit(NsegmentOutVars) )
        Nsegment_var_monthly = 0.0D0
        Monthlyunit = 0
      ENDIF

      IF ( NsegmentOutON_OFF==2 ) THEN
        IF ( getparam(MODNAME, 'nhm_seg', Nsegment, 'integer', Nhm_seg)/=0 ) CALL read_error(2, 'nhm_seg')
      ENDIF
      WRITE ( Output_fmt2, 9002 ) Nsegment
      ALLOCATE ( Nsegment_var_daily(Nsegment, NsegmentOutVars) )
      Nsegment_var_daily = 0.0
      DO jj = 1, NsegmentOutVars
        IF ( Daily_flag==ACTIVE ) THEN
          fileName = NsegmentOutBaseFileName(:numchars(NsegmentOutBaseFileName))//NsegmentOutVar_names(jj)(:Nc_vars(jj))//'.csv'
          !print *, fileName
          CALL PRMS_open_output_file(Dailyunit(jj), fileName, 'xxx', 0, ios)
          IF ( ios/=0 ) CALL error_stop('in nsegment_summary, daily', ERROR_open_out)
          IF ( NsegmentOutON_OFF==1 ) THEN
            WRITE ( Dailyunit(jj), Output_fmt2 ) (j, j=1,Nsegment)
          ELSE
            WRITE ( Dailyunit(jj), Output_fmt2 ) (Nhm_seg(j), j=1,Nsegment)
          ENDIF
        ENDIF
        IF ( NsegmentOut_freq>MEAN_MONTHLY ) THEN
          IF ( NsegmentOut_freq==MEAN_YEARLY ) THEN
            fileName = NsegmentOutBaseFileName(:numchars(NsegmentOutBaseFileName))//NsegmentOutVar_names(jj)(:Nc_vars(jj))// &
      &                 '_meanyearly.csv'
            CALL PRMS_open_output_file(Yearlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nsegment_summary, mean yearly', ERROR_open_out)
          ELSE  !IF ( NsegmentOut_freq==YEARLY ) THEN
            fileName = NsegmentOutBaseFileName(:numchars(NsegmentOutBaseFileName))//NsegmentOutVar_names(jj)(:Nc_vars(jj))// &
      &                 '_yearly.csv'
            CALL PRMS_open_output_file(Yearlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nsegment_summary, yearly', ERROR_open_out)
          ENDIF
          IF ( NsegmentOutON_OFF==1 ) THEN
            WRITE ( Yearlyunit(jj), Output_fmt2 ) (j, j=1,Nsegment)
          ELSE
            WRITE ( Yearlyunit(jj), Output_fmt2 ) (Nhm_seg(j), j=1,Nsegment)
          ENDIF
        ENDIF
        IF ( Monthly_flag==ACTIVE ) THEN
          IF ( NsegmentOut_freq==MEAN_MONTHLY ) THEN
            fileName = NsegmentOutBaseFileName(:numchars(NsegmentOutBaseFileName))//NsegmentOutVar_names(jj)(:Nc_vars(jj))// &
     &                 '_meanmonthly.csv'
            CALL PRMS_open_output_file(Monthlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nsegment_summary, mean monthly', ERROR_open_out)
          ELSE
            fileName = NsegmentOutBaseFileName(:numchars(NsegmentOutBaseFileName))//NsegmentOutVar_names(jj)(:Nc_vars(jj))// &
     &                  '_monthly.csv'
            CALL PRMS_open_output_file(Monthlyunit(jj), fileName, 'xxx', 0, ios)
            IF ( ios/=0 ) CALL error_stop('in nsegment_summary, monthly', ERROR_open_out)
          ENDIF
          IF ( NsegmentOutON_OFF==1 ) THEN
            WRITE ( Monthlyunit(jj), Output_fmt2 ) (j, j=1,Nsegment)
          ELSE
            WRITE ( Monthlyunit(jj), Output_fmt2 ) (Nhm_seg(j), j=1,Nsegment)
          ENDIF
        ENDIF
      ENDDO

 9001 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',ES10.3))')
 9002 FORMAT ('("Date"',I0,'('', ''I0))')
 9003 FORMAT ('(I4,', I0,'('','',ES10.3))')
 9005 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.4))')
 9006 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.3))')
 9007 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.2))')
 9008 FORMAT ('(I4,', I0,'('','',F0.4))')
 9009 FORMAT ('(I4,', I0,'('','',F0.3))')
 9010 FORMAT ('(I4,', I0,'('','',F0.2))')
 9011 FORMAT ('(I4,', I0,'('','',F0.5))')
 9012 FORMAT ('(I4, 2(''-'',I2.2),',I0,'('','',F0.5))')

      END SUBROUTINE nsegment_summaryinit

!***********************************************************************
!     Output set of declared variables in CSV format
!***********************************************************************
      SUBROUTINE nsegment_summaryrun()
      USE PRMS_NSEGMENT_SUMMARY
      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Nowday, Modays
      IMPLICIT NONE
! FUNCTIONS AND SUBROUTINES
      INTRINSIC SNGL, DBLE
      INTEGER, EXTERNAL :: getvar
      EXTERNAL :: read_error
! Local Variables
      INTEGER :: j, i, jj, write_month, last_day
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
      DO jj = 1, NsegmentOutVars
        IF ( Nsegment_var_type(jj)==REAL_TYPE ) THEN
          IF ( getvar(MODNAME, NsegmentOutVar_names(jj)(:Nc_vars(jj)), Nsegment, 'real', Nsegment_var_daily(1, jj))/=0 ) &
     &         CALL read_error(4, NsegmentOutVar_names(jj)(:Nc_vars(jj)))
        ELSEIF ( Nsegment_var_type(jj)==DBLE_TYPE ) THEN
          IF ( getvar(MODNAME, NsegmentOutVar_names(jj)(:Nc_vars(jj)), Nsegment, 'double', Nsegment_var_dble(1, jj))/=0 ) &
     &         CALL read_error(4, NsegmentOutVar_names(jj)(:Nc_vars(jj)))
          DO i = 1, Nsegment
            Nsegment_var_daily(i, jj) = SNGL( Nsegment_var_dble(i, jj) )
          ENDDO
        ENDIF
        IF ( Daily_flag==ACTIVE ) WRITE ( Dailyunit(jj), Output_fmt) Nowyear, Nowmonth, Nowday, &
     &                                                           (Nsegment_var_daily(j,jj), j=1,Nsegment)
      ENDDO

      write_month = OFF
      IF ( NsegmentOut_freq>MEAN_MONTHLY ) THEN
        last_day = OFF
        IF ( Nowyear==End_year .AND. Nowmonth==End_month .AND. Nowday==End_day ) last_day = ACTIVE
        IF ( Lastyear/=Nowyear .OR. last_day==ACTIVE ) THEN
          IF ( (Nowmonth==Start_month .AND. Nowday==Start_day) .OR. last_day==ACTIVE ) THEN
            DO jj = 1, NsegmentOutVars
              IF ( NsegmentOut_freq==MEAN_YEARLY ) THEN
                DO i = 1, Nsegment
                  Nsegment_var_yearly(i, jj) = Nsegment_var_yearly(i, jj)/Yeardays
                ENDDO
              ENDIF
              WRITE ( Yearlyunit(jj), Output_fmt3) Lastyear, (Nsegment_var_yearly(j,jj), j=1,Nsegment)
            ENDDO
            Nsegment_var_yearly = 0.0D0
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

      IF ( NsegmentOut_freq>MEAN_MONTHLY ) THEN
        DO jj = 1, NsegmentOutVars
          DO i = 1, Nsegment
            Nsegment_var_yearly(i, jj) = Nsegment_var_yearly(i, jj) + DBLE( Nsegment_var_daily(i, jj) )
          ENDDO
        ENDDO
        RETURN
      ENDIF

      IF ( Monthly_flag==ACTIVE ) THEN
        DO jj = 1, NsegmentOutVars
          DO i = 1, Nsegment
            Nsegment_var_monthly(i, jj) = Nsegment_var_monthly(i, jj) + DBLE( Nsegment_var_daily(i, jj) )
            IF ( write_month==ACTIVE ) THEN
              IF ( NsegmentOut_freq==MEAN_MONTHLY ) Nsegment_var_monthly(i, jj) = Nsegment_var_monthly(i, jj)/Monthdays
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF ( write_month==ACTIVE ) THEN
        DO jj = 1, NsegmentOutVars
          WRITE ( Monthlyunit(jj), Output_fmt) Nowyear, Nowmonth, Nowday, &
     &                                         (Nsegment_var_monthly(j,jj), j=1,Nsegment)
        ENDDO
        Monthdays = 0.0D0
        Nsegment_var_monthly = 0.0D0
      ENDIF

      END SUBROUTINE nsegment_summaryrun
