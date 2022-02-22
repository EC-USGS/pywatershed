      ! utils_prms.f90 2020-07-30
!***********************************************************************
!     Read CBH File to current time
!***********************************************************************
      SUBROUTINE find_current_time(Iunit, Year, Month, Day, Iret, Cbh_binary_flag)
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Iunit, Year, Month, Day, Cbh_binary_flag
      INTEGER, INTENT(OUT) :: Iret
! Local Variables
      INTEGER :: yr, mo, dy
      CHARACTER(LEN=20), PARAMETER :: fmt1 = '(A, I5, 2("/",I2.2))'
!***********************************************************************
      Iret = 0
      DO
        IF ( Cbh_binary_flag==0 ) THEN
          READ ( Iunit, *, IOSTAT=Iret ) yr, mo, dy
        ELSE
          READ ( Iunit, IOSTAT=Iret ) yr, mo, dy
        ENDIF
        IF ( Iret==-1 ) PRINT fmt1, 'ERROR, end-of-file found reading input file for date:', Year, Month, Day
        IF ( Iret/=0 ) RETURN
        IF ( yr==Year .AND. mo==Month .AND. dy==Day ) EXIT
      ENDDO
      BACKSPACE Iunit
      END SUBROUTINE find_current_time

!***********************************************************************
!     Read File dynamic paramter file to current time
!***********************************************************************
      SUBROUTINE find_current_file_time(Iunit, Year, Month, Day, Year_file, Month_file, Day_file)
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Iunit, Year, Month, Day
      INTEGER, INTENT(OUT) :: Year_file, Month_file, Day_file
! Local Variables
      INTEGER :: i, ios
!***********************************************************************
! find first value for simulation time period
      READ ( Iunit, *, IOSTAT=ios ) Year_file, Month_file, Day_file
      IF ( ios/=0 ) THEN
        Year_file = 0
        Month_file = 0
        Day_file = 0
        RETURN
      ENDIF
      IF ( Year_file<Year ) THEN
        i = 0
        DO WHILE ( i==0 )
          READ ( Iunit, *, IOSTAT=ios ) Year_file, Month_file, Day_file
          IF ( ios/=0 ) THEN
            Year_file = 0
            Month_file = 0
            Day_file = 0
            RETURN
          ENDIF
          IF ( Year_file>=Year ) i = 1
        ENDDO
      ENDIF
      IF ( Year_file==Year ) THEN
        IF ( Month_file<Month ) THEN
          i = 0
          DO WHILE ( i==0 )
            READ ( Iunit, *, IOSTAT=ios ) Year_file, Month_file, Day_file
            IF ( ios/=0 ) THEN
              Year_file = 0
              Month_file = 0
              Day_file = 0
              RETURN
            ENDIF
            IF ( Month_file>=Month .OR. Year_file/=Year ) i = 1
          ENDDO
        ENDIF
        IF ( Year_file==Year .AND. Month_file==Month ) THEN
          IF ( Day_file<Day ) THEN
            i = 0
            DO WHILE ( i==0 )
              READ ( Iunit, *, IOSTAT=ios ) Year_file, Month_file, Day_file
              IF ( ios/=0 ) THEN
                Year_file = 0
                Month_file = 0
                Day_file = 0
                RETURN
              ENDIF
              IF ( Day_file>=Day ) i = 1
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      BACKSPACE Iunit
      END SUBROUTINE find_current_file_time

!***********************************************************************
!     Read File to line before data starts in file
!***********************************************************************
      SUBROUTINE find_header_end(Iunit, Fname, Paramname, Iret, Cbh_flag, Cbh_binary_flag)
      USE PRMS_CONSTANTS, ONLY: DEBUG_less
      USE PRMS_MODULE, ONLY: Nhru, Orad_flag, Print_debug
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Cbh_flag, Cbh_binary_flag
      INTEGER, INTENT(OUT) :: Iunit
      INTEGER, INTENT(INOUT) :: Iret
      CHARACTER(LEN=*), INTENT(IN) :: Fname, Paramname
! Functions
      INTRINSIC :: trim
      INTEGER, EXTERNAL :: get_ftnunit
! Local Variables
      INTEGER :: i, ios, dim
      CHARACTER(LEN=4) :: dum
      CHARACTER(LEN=80) :: dum2
!***********************************************************************
      IF ( Iret/=2 ) Iret = 0
      Iunit = get_ftnunit(7777)
      IF ( Cbh_binary_flag==0 ) THEN
        OPEN ( Iunit, FILE=trim(Fname), STATUS='OLD', IOSTAT=ios )
      ELSE
        OPEN ( Iunit, FILE=trim(Fname), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=ios, ACCESS='STREAM' )
      ENDIF
      IF ( ios/=0 ) THEN
        IF ( Iret==2 ) THEN ! this signals climate_hru to ignore the Humidity CBH file, could add other files
          Iret = 0
          IF ( Print_debug>DEBUG_less ) &
     &         WRITE ( *, '(/,A,/,A,/,A)' ) 'WARNING, optional CBH file not found, will use associated parameter values'
        ELSE
          WRITE ( *, '(/,A,/,A,/,A)' ) 'ERROR reading file:', Fname, 'check to be sure the input file exists'
          Iret = 1
        ENDIF
      ELSE
! read to line before data starts in each file
        i = 0
        DO WHILE ( i==0 )
          IF ( Cbh_binary_flag==0 ) THEN
            READ ( Iunit, FMT='(A4)', IOSTAT=ios ) dum
          ELSE
            READ ( Iunit, IOSTAT=ios ) dum2
            READ ( dum2, '(A4)' ) dum
          ENDIF
          IF ( ios/=0 ) THEN
            WRITE ( *, '(/,A,/,A,/,A)' ) 'ERROR reading file:', Fname, 'check to be sure the input file is in correct format'
            Iret = 1
            EXIT
          ELSEIF ( dum=='####' ) THEN
            IF ( Cbh_flag==0 ) EXIT
            BACKSPACE Iunit
            BACKSPACE Iunit
            IF ( Orad_flag==1 .AND. Paramname(:5)=='swrad' ) BACKSPACE Iunit ! backspace again as swrad CBH file contains orad as last column
            IF ( Cbh_binary_flag==0 ) THEN
              READ ( Iunit, *, IOSTAT=ios ) dum, dim
            ELSE
              READ ( Iunit, IOSTAT=ios ) dum2
              READ ( dum2, * ) dum, dim
            ENDIF
            IF ( ios/=0 ) THEN
              WRITE ( *, '(/,A,/,A,/,A)' ) 'ERROR reading file:', Fname, 'check to be sure dimension line is in correct format'
              Iret = 1
              EXIT
            ENDIF
            IF ( dim/=Nhru ) THEN
              PRINT '(/,2(A,I0))', '***CBH file dimension incorrect*** nhru= ', Nhru, ' CBH dimension= ', dim, ' File: '//Fname
              PRINT *, 'ERROR: update Control File with correct CBH files'
              Iret = 1
              EXIT
            ENDIF
            IF ( Cbh_binary_flag==0 ) THEN
              READ ( Iunit, FMT='(A4)', IOSTAT=ios ) dum
            ELSE
              READ ( Iunit, IOSTAT=ios ) dum
            ENDIF
            IF ( ios/=0 ) THEN
              WRITE ( *, '(/,A,/,A,/)' ) 'ERROR reading file:', Fname
              Iret = 1
              EXIT
            ENDIF
            IF ( Orad_flag==1 .AND. Paramname(:5)=='swrad' ) READ ( Iunit, FMT='(A4)' ) dum ! read again as swrad CBH file contains orad as last column
            i = 1
          ENDIF
        ENDDO
      ENDIF

      END SUBROUTINE find_header_end

!**********************
! Check for end of file
!**********************
      SUBROUTINE is_eof(Iunit, Next_yr, Next_mo, Next_day)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit
      INTEGER, INTENT(OUT) :: Next_yr, Next_mo, Next_day
! Local Variables
      INTEGER :: ios, i
      CHARACTER(LEN=80) :: dum
!*******************************************************************************
      Next_yr = 0
      Next_mo = 0
      Next_day = 0
      i = 0
      DO WHILE ( i==0 )
        READ ( Iunit, '(A)', iostat=ios ) dum
        IF ( ios/=0 ) RETURN
        IF ( dum(:4)=='####' ) CYCLE
        IF ( dum(:2)/='//' ) i = 1
      ENDDO
      READ ( dum, *, iostat=ios ) Next_yr, Next_mo, Next_day
      IF ( ios/=0 ) THEN
        Next_yr = 0
        Next_mo = 0
        Next_day = 0
      ELSE
        BACKSPACE Iunit
      ENDIF
      END SUBROUTINE is_eof

!***********************************************************************
!     Determine an unopened FORTRAN File Unit
!***********************************************************************
      INTEGER FUNCTION get_ftnunit(Iunit)
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Iunit
! Local Variables
      INTEGER :: good_unit
      LOGICAL :: opend
!***********************************************************************
      good_unit = Iunit
      opend = .TRUE.
      DO WHILE ( opend )
        good_unit = good_unit + 1
        INQUIRE ( UNIT=good_unit, OPENED=opend )
      ENDDO
      get_ftnunit = good_unit
      END FUNCTION get_ftnunit

!***********************************************************************
! Convert Fahrenheit to Celsius
!***********************************************************************
      REAL FUNCTION f_to_c(Temp)
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Temp
!***********************************************************************
      f_to_c = (Temp-32.0)/1.8
      END FUNCTION f_to_c

!***********************************************************************
! Convert Celsius to Fahrenheit
!***********************************************************************
      REAL FUNCTION c_to_f(Temp)
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Temp
!***********************************************************************
      c_to_f = Temp*1.8 + 32.0
      END FUNCTION c_to_f

!***********************************************************************
      SUBROUTINE write_integer_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit, Dimen
      INTEGER, INTENT(IN) :: Values(Dimen)
      CHARACTER(LEN=*), INTENT(IN) :: Parm_name, Dimen_name
! Local Variables
      INTEGER i
      CHARACTER(LEN=48), PARAMETER :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, /, "1")'
!***********************************************************************
      WRITE ( Iunit, fmt1 ) Parm_name, Dimen_name, Dimen
      DO i = 1, Dimen
        WRITE ( Iunit, '(I0)' ) Values(i)
      ENDDO
      END SUBROUTINE write_integer_param

!***********************************************************************
      SUBROUTINE write_real_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit, Dimen
      REAL, INTENT(IN) :: Values(Dimen)
      CHARACTER(LEN=*), INTENT(IN) :: Parm_name, Dimen_name
! Local Variables
      INTEGER i
      CHARACTER(LEN=48), PARAMETER :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, /, "2")'
!***********************************************************************
      WRITE ( Iunit, fmt1) Parm_name, Dimen_name, Dimen
      DO i = 1, Dimen
        WRITE ( Iunit, * ) Values(i)
      ENDDO
      END SUBROUTINE write_real_param

!***********************************************************************
      SUBROUTINE write_double_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit, Dimen
      DOUBLE PRECISION, INTENT(IN) :: Values(Dimen)
      CHARACTER(LEN=*), INTENT(IN) :: Parm_name, Dimen_name
! Local Variables
      INTEGER i
      CHARACTER(LEN=40), PARAMETER :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, "3")'
!***********************************************************************
      WRITE ( Iunit, fmt1 ) Parm_name, Dimen_name, Dimen
      DO i = 1, Dimen
        WRITE ( Iunit, * ) Values(i)
      ENDDO
      END SUBROUTINE write_double_param

!***********************************************************************
      SUBROUTINE write_2D_double_param(Iunit, Parm_name, Dimen_name1, Dimen1, &
     &                                 Dimen_name2, Dimen2, Values)
!***********************************************************************
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Dimen1, Dimen2, Iunit
      DOUBLE PRECISION, INTENT(IN) :: Values(Dimen1, Dimen2)
      CHARACTER(LEN=*), INTENT(IN) :: Parm_name
      CHARACTER(LEN=*), INTENT(IN) :: Dimen_name1, Dimen_name2
! Local Variables
      INTEGER i, j
      CHARACTER(LEN=46), PARAMETER :: fmt1 = '("####", /, A, /, "2", /, A, /, A, /, I0, "3")'
!***********************************************************************
      WRITE ( Iunit, fmt1 ) Parm_name, Dimen_name1, Dimen_name2, Dimen1*Dimen2
      DO i = 1, Dimen2
        DO j = 1, Dimen1
          WRITE ( Iunit, * ) Values(j, i)
        ENDDO
      ENDDO
      END SUBROUTINE write_2D_double_param

!***********************************************************************
      SUBROUTINE write_2D_double_array_grid(Iunit, Parm_name, Dimen_name1, &
     &                                      Dimen1, Dimen_name2, Dimen2, Values)
!***********************************************************************
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit, Dimen1, Dimen2
      DOUBLE PRECISION, INTENT(IN) :: Values(Dimen1, Dimen2)
      CHARACTER(LEN=*), INTENT(IN) :: Parm_name, Dimen_name1, Dimen_name2
! Local Variables
      INTEGER i, j
      CHARACTER(LEN=12) :: fmt
!***********************************************************************
      WRITE ( Iunit, 9001) Parm_name, Dimen_name1, Dimen_name2, Dimen1*Dimen2
      WRITE ( fmt, 9002 ) Dimen2
      DO i = 1, Dimen2
        WRITE ( Iunit, fmt ) (Values(j, i), j=1,Dimen1)
      ENDDO

 9001 FORMAT ( '####', /, A, /, '2', /, A, /, A, /, I0, /, '3' )
 9002 FORMAT ( '(', I5, 'F10.5)' )
      END SUBROUTINE write_2D_double_array_grid

!**********************************************************************
!     Version Check
!**********************************************************************
      SUBROUTINE version_check(Module_version, Length, Param_version)
      IMPLICIT NONE
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Module_version, Param_version
      INTEGER, INTENT(IN) :: Length
!**********************************************************************
      IF ( Module_version(13:Length+12)/=Param_version(:Length) ) THEN
        PRINT 9001, Module_version(13:Length+12), Param_version(:Length)
        PRINT *, 'Enter return to continue'
        READ (*, *)
      ENDIF
 9001 FORMAT ('WARNING, module versions are not identical', /, &
     &        'Executable version: ', A, /, &
     &        'Parameter File version: ', A, /)
      END SUBROUTINE version_check

!**********************************************************************
!     Parameter or Variable delcare or read error
!**********************************************************************
      SUBROUTINE read_error(Iflag, Name)
      USE PRMS_CONSTANTS, ONLY: ERROR_decl_get
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iflag
      CHARACTER(LEN=*), INTENT(IN) :: Name
!**********************************************************************
      PRINT '(/,A,/)', 'Due to error condition simulation halted'
      IF ( Iflag==1 ) THEN
        PRINT *, 'Declare error for parameter: ', Name
      ELSEIF ( Iflag==2 ) THEN
        PRINT *, 'Get error for parameter: ', Name
      ELSEIF ( Iflag==3 ) THEN
        PRINT *, 'Declare error for variable: ', Name
      ELSEIF ( Iflag==4 ) THEN
        PRINT *, 'Get error for variable: ', Name
      ELSEIF ( Iflag==5 ) THEN
        PRINT *, 'Read error for control parameter: ', Name
      ELSEIF ( Iflag==6 ) THEN
        PRINT *, 'Read error for dimension parameter: ', Name
      ELSEIF ( Iflag==7 ) THEN
        PRINT *, 'Declare error for dimension parameter: ', Name
      ELSEIF ( Iflag==8 ) THEN
        PRINT *, 'Declare error for Data File variable: ', Name
      ELSEIF ( Iflag==9 ) THEN
        PRINT *, 'Read error for Data File variable: ', Name
      ELSEIF ( Iflag==10 ) THEN
        PRINT *, 'Open error of Control File ', Name
      ELSEIF ( Iflag==11 ) THEN
        PRINT *, 'Read error of Parameter File ', Name
      ELSEIF ( Iflag==12 ) THEN
        PRINT *, 'Read error of Control File ', Name
      ELSEIF ( Iflag==13 ) THEN
        PRINT *, 'Read error of Data File ', Name
      ELSEIF ( Iflag==14 ) THEN
        PRINT *, 'Control parameter not found: ', Name
      ELSEIF ( Iflag==15 ) THEN
        PRINT *, 'ERROR, control ', Name, ' expected and is not available in PRMS'
      ELSEIF ( Iflag==16 ) THEN
        PRINT *, 'ERROR, declared parameter ', Name
      ENDIF
      ERROR STOP ERROR_decl_get
      END SUBROUTINE read_error

!**********************************************************************
!     Module error
!**********************************************************************
      SUBROUTINE module_error(Modname, Arg, Retcode)
      USE PRMS_CONSTANTS, ONLY: ERROR_module
      IMPLICIT NONE
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Modname, Arg
      INTEGER, INTENT(IN) :: Retcode
!**********************************************************************
      PRINT 9001, Modname, Arg, Retcode
      ERROR STOP ERROR_module
 9001 FORMAT ('ERROR in ', A, ' module, arg = ', A, /, 'Return val = ', I0)
      END SUBROUTINE module_error

!***********************************************************************
! Compute saturation vapor pressure over water in millibars
! 6th order Polynominal method (Flatau et. all., 1992) valid: -50 to 50C
! 1 kPa = 10 millibars
! Flatau, P.j., Walko, R.L., Cotton, W.R., 1992, Polynomial Fits to
!   saturation vapor pressure: Jornal of Applied Meteorology, v. 31, p. 1507-1513
!***********************************************************************
      REAL FUNCTION sat_vapor_press_poly(Tempc)
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Tempc
!***********************************************************************
      sat_vapor_press_poly = 6.11176750 + 0.443986062*Tempc &
     &                       + 0.0143053301*Tempc**2 &
     &                       + 0.265027242E-03*Tempc**3 &
     &                       + 0.302246994E-05*Tempc**4 &
     &                       + 0.203886313E-07*Tempc**5 &
     &                       + 0.638780966E-10*Tempc**6
! Mastin documentation for potet_dpm
!      sat_vapor_press_poly = 23.38*exp(18.1-5303.3/(Tempc+273.0))
! Mastin documentation for param_leaf-loss.aml
!      sat_vapor_press_poly = 6.1078*EXP(17.269*Tempc/(237.30D0+Tempc))
! Buck Research Manual (1996)
!      sat_vapor_press_poly = 6.1121D0*EXP((18.678D0-Tempc/234.5D0)*Tempc/(257.14+Tempc))
! WMO 2008, CIMO Guide
!      sat_vapor_press_poly = 6.112*EXP(17.62*Tempc/(243.12+Tempc))
! Irmak and others (2012), equation 12
!      sat_vapor_press_poly = 0.6108*EXP(17.27*Tempc/(237.3+Tempc))
      END FUNCTION sat_vapor_press_poly

!***********************************************************************
! Compute saturation vapor pressure over water
! Irmak and others (2012), equation 12
!***********************************************************************
      REAL FUNCTION sat_vapor_press(Tempc)
      IMPLICIT NONE
! Arguments
      REAL, INTENT(IN) :: Tempc
!***********************************************************************
      sat_vapor_press = 6.1078*EXP( (17.26939*Tempc)/(237.3+Tempc) )
      END FUNCTION sat_vapor_press

!***********************************************************************
! leap_day - is the year a leap year: (1=yes; 0=no)
!***********************************************************************
      INTEGER FUNCTION leap_day(Year)
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(IN) :: Year
      ! Functions
      INTRINSIC MOD
!***********************************************************************
      leap_day = 0
      ! Check if leapyear - Start by identifying all years not divisible by 4
      IF ( MOD(Year,4)==0 ) THEN
        leap_day = 1
        IF ( MOD(Year,100)==0 ) THEN
          IF ( MOD(Year,400)/=0 ) leap_day = 0
        ENDIF
      ENDIF
      END FUNCTION leap_day

!***********************************************************************
! write_outfile - print to model output file
!***********************************************************************
      SUBROUTINE write_outfile(String)
      USE PRMS_CONSTANTS, ONLY: DEBUG_minimum
      USE PRMS_MODULE, ONLY: PRMS_output_unit, Print_debug
      IMPLICIT NONE
      ! Functions
      INTRINSIC :: LEN_TRIM
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: String
      ! Local variable
      INTEGER nchars
!***********************************************************************
      IF ( Print_debug==DEBUG_minimum ) RETURN
      nchars = LEN_TRIM(String)
      IF ( nchars>0 ) THEN
        WRITE ( PRMS_output_unit, '(A)' ) String(:nchars)
      ELSE
        WRITE ( PRMS_output_unit, '(/)' )
      ENDIF
      END SUBROUTINE write_outfile

!***********************************************************************
! julian_day
! computes the Julian date given a Gregorian calendar date
! (Year, Month, Day) relative to: calendar (Jan 1),
! solar (12/22 in Northern; 6/21 in Southern) and
! water year (10/1 in Northern; 4/1 in Southern) start dates.
! The Julian day starts at noon of the Gregorian day and
! extends to noon the next Gregorian day.
!***********************************************************************
      INTEGER FUNCTION julian_day(Date_type, Year_type)
      USE PRMS_CONSTANTS, ONLY: NORTHERN, ERROR_time
      USE PRMS_MODULE, ONLY: Starttime, Endtime
      USE PRMS_BASIN, ONLY: Hemisphere
      USE PRMS_SET_TIME, ONLY: Nowtime
      IMPLICIT NONE
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Date_type ! "start", "end", "now"
      CHARACTER(LEN=*), INTENT(IN) :: Year_type ! "calendar", "solar", "water", "absolute"
      ! Functions
      INTEGER, EXTERNAL :: compute_julday
      ! Local Variables
      INTEGER :: reftime_year, reftime_month, reftime_day, time_array(6)
      INTEGER :: year, month, day, absolute_julday, relative_julday, length, found
!***********************************************************************
      IF ( Date_type(:3)=='end' ) THEN
        time_array = Endtime
      ELSEIF ( Date_type(:3)=='now' ) THEN
        time_array = Nowtime
      ELSEIF ( Date_type(:5)=='start' ) THEN
        time_array = Starttime
      ELSE
        PRINT *, 'ERROR, invalid argument to compute Julian Day: ', Date_type
        ERROR STOP ERROR_time
      ENDIF
      year = time_array(1)
      month = time_array(2)
      day = time_array(3)

      found = 0
      length = LEN(Year_type)
      ! set reftime depending on type arg
      IF ( length>4 ) THEN
        IF ( Year_type(:5)=='solar' ) THEN
          found = 1
          IF ( Hemisphere==NORTHERN ) THEN
            IF ( month==12 .AND. day>21 ) THEN
              reftime_year = year
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 12
            reftime_day = 21
          ELSE ! Southern
            IF ( month==6 .AND. day>20 ) THEN
              reftime_year = year;
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 6
            reftime_day = 20
          ENDIF
        ELSEIF ( Year_type(:5)=='water' ) THEN
          found = 1
          IF ( Hemisphere==NORTHERN ) THEN
            IF ( month>9 ) THEN
              reftime_year = year
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 9
            reftime_day = 30
          ELSE ! Southern
            IF ( month>3 ) THEN
              reftime_year = year
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 3
            reftime_day = 31
          ENDIF
        ENDIF
      ENDIF
      IF ( found==0 .AND. length>5 ) THEN
        IF ( Year_type(:6)=='spring' ) THEN
          found = 1
          IF ( Hemisphere==NORTHERN ) THEN
            IF ( month>3 .OR. (month==3 .AND. day>20) ) THEN
              reftime_year = year
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 3
            reftime_day = 20
          ELSE ! Southern
            IF ( month>9 .OR. (month==9 .AND. day>22) ) THEN
              reftime_year = year
            ELSE
              reftime_year = year - 1
            ENDIF
            reftime_month = 9
            reftime_day = 22
          ENDIF
        ENDIF
      ENDIF
      IF ( found==0 .AND. length>7 ) THEN
        IF ( Year_type(:8)=='calendar' ) THEN
          found = 1
          reftime_year = year - 1
          reftime_month = 12
          reftime_day = 31
        ENDIF
      ENDIF
      IF ( found==0 ) THEN
        PRINT *, 'ERROR, invalid year type argument to compute Julian Day: ', Year_type
        ERROR STOP ERROR_time
      ENDIF

      ! set actual Julian Day
      absolute_julday = compute_julday(year, month, day)

      relative_julday = 0
      IF ( length==8 ) THEN
        IF ( Year_type(:8)=='calendar' ) relative_julday = compute_julday(reftime_year, reftime_month, reftime_day)
      ELSE
        relative_julday = compute_julday(reftime_year, reftime_month, reftime_day)
      ENDIF
      julian_day = absolute_julday - relative_julday

      END FUNCTION julian_day

!***********************************************************************
! compute_julday
! computes the Julian Day given a Gregorian calendar date
!***********************************************************************
      INTEGER FUNCTION compute_julday(Year, Month, Day)
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(IN) :: Year, Month, Day
      ! Local Variables
      !INTEGER yr, mo
!!***********************************************************************
!      mo = Month
!      yr = Year
!      IF ( Month < 3 ) THEN
!        mo = mo + 12
!        yr = yr - 1
!      ENDIF
!      compute_julday = Day + (153*mo - 457) / 5 + 365*yr + (yr/4) - (yr/100) + (yr/400) + 1721118.5
      compute_julday = Day - 32075 + 1461*(Year+4800+(Month-14)/12)/4 + 367*(Month-2-(Month-14)/12*12) &
     &    /12-3*((Year+4900+(Month-14)/12)/100)/4

      END FUNCTION compute_julday

!***********************************************************************
! compute_gregorian
! computes the Gregorian calendar date given the Julian Day
!***********************************************************************
      SUBROUTINE compute_gregorian(Julday, Year, Month, Day)
      !USE PRMS_CONSTANTS, ONLY: DAYS_YR
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(OUT) :: Year, Month, Day
      INTEGER, INTENT(IN) :: Julday
      ! Functions
      INTRINSIC FLOOR, NINT
      ! Local Variables
      INTEGER m, n
!***********************************************************************
      m = Julday + 68569
      n = 4*m/146097
      m = m - (146097*n+3)/4
      Year = 4000*(m+1)/1461001
      m = m - 1461*Year/4+31
      Month = 80*m/2447
      Day = m - 2447*Month/80
      m = Month/11
      Month = Month + 2 - 12*m
      Year = 100*(n-49) + Year + m

      !
      !Z = Julday + 0.5
      !W = FLOOR((Z - 1867216.25)/36524.25)
      !X = FLOOR(W/4.0)
      !A = Z + 1.0 + W - X
      !B = A + 1524.0
      !C = FLOOR((B - 122.1)/DAYS_YR)
      !D = FLOOR(DAYS_YR*C)
      !E = FLOOR((B - D)/30.6001)
      !F = FLOOR(30.6001*E)
      !Day = NINT(B - D -F)
      !Month = NINT(E - 1)
      !IF ( Month>12 ) Month = Month - 12
      !IF ( Month<3 ) THEN
      !  Year = NINT(C - 4715)
      !ELSE
      !  Year = NINT(C - 4716)
      !ENDIF
      END SUBROUTINE compute_gregorian

!***********************************************************************
! julday_in_year
! computes the Julian Day of a date
!***********************************************************************
      INTEGER FUNCTION julday_in_year(Year, Month, Day)
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(IN) :: Year, Month, Day
      ! Functions
      INTEGER, EXTERNAL :: leap_day
      ! Local Variables
      INTEGER daypmo(12), i
      DATA daypmo/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!***********************************************************************
      daypmo(2) = 28
      IF ( leap_day(Year)==1 ) daypmo(2) = 29
      julday_in_year = Day
      DO i = 1, Month - 1
        julday_in_year = julday_in_year + daypmo(i)
      ENDDO
      END FUNCTION julday_in_year

!***********************************************************************
!     Open PRMS input File and assign unit number
!***********************************************************************
      SUBROUTINE PRMS_open_input_file(Iunit, Fname, Paramname, Ftype, Iret)
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Ftype
      INTEGER, INTENT(OUT) :: Iunit, Iret
      CHARACTER(LEN=*), INTENT(IN) :: Fname, Paramname
! Functions
      INTEGER, EXTERNAL :: get_ftnunit, numchars
! Local Variables
      INTEGER :: ios, nchars
!***********************************************************************
      Iret = 0
      Iunit = get_ftnunit(7777)
      nchars = numchars(Fname)
      IF ( Ftype==0 ) THEN
        OPEN ( Iunit, FILE=Fname(:nchars), STATUS='OLD', IOSTAT=ios )
      ELSE
        OPEN ( Iunit, FILE=Fname(:nchars), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=ios, ACCESS='STREAM' )
      ENDIF
      IF ( ios/=0 ) THEN
        WRITE ( *, '(/,2A,/,A,/,2A,/)' ) 'ERROR opening input file: ', Fname(:nchars), &
     &                                   'check to be sure the input file exists', &
     &                                   'file specified by control parameter: ', Paramname
        Iret = 1
      ENDIF
      END SUBROUTINE PRMS_open_input_file

!***********************************************************************
!     Open PRMS output file and assign unit number
!***********************************************************************
      SUBROUTINE PRMS_open_output_file(Iunit, Fname, Paramname, Ftype, Iret)
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(IN) :: Ftype ! 0=text; 1=BINARY
      INTEGER, INTENT(OUT) :: Iunit, Iret
      CHARACTER(LEN=*), INTENT(IN) :: Fname, Paramname
! Functions
      INTEGER, EXTERNAL :: get_ftnunit, numchars
! Local Variables
      INTEGER :: ios, nchars
!***********************************************************************
      Iret = 0
      Iunit = get_ftnunit(8888)
      nchars = numchars(Fname)

      IF ( Ftype==0 ) THEN
        OPEN ( Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios )
      ELSE
        OPEN ( Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios, FORM='UNFORMATTED', ACCESS='STREAM' )
      ENDIF

      IF ( ios/=0 ) THEN
        WRITE ( *, '(/,A,/,A,/)' ) 'ERROR opening output file:', Fname(:nchars), &
     &                             'check to be sure the pathname is valid and the file is not open'
        WRITE ( *, '(2A,/)' ) 'file specified by control parameter: ', Paramname
        Iret = 1
      ENDIF

      END SUBROUTINE PRMS_open_output_file

!***********************************************************************
!     Open PRMS module output file and assign unit number
!***********************************************************************
      SUBROUTINE PRMS_open_module_file(Iunit, Fname)
      USE PRMS_CONSTANTS, ONLY: ERROR_open_in
      IMPLICIT NONE
! Argument
      INTEGER, INTENT(OUT) :: Iunit
      CHARACTER(LEN=*), INTENT(IN) :: Fname
! Functions
      INTEGER, EXTERNAL :: get_ftnunit, numchars
! Local Variables
      INTEGER :: ios, nchars
!***********************************************************************
      Iunit = get_ftnunit(8888)
      nchars = numchars(Fname)
      OPEN ( Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios )
      IF ( ios/=0 ) THEN
        WRITE ( *, '(/,A,/,A,/)' ) 'ERROR opening water balance output file:', Fname(:nchars), &
     &                             'check to be sure the pathname is valid and the file is not open'
        ERROR STOP ERROR_open_in
      ENDIF
      END SUBROUTINE PRMS_open_module_file

!***********************************************************************
!     Determine number of characters in a string
!***********************************************************************
      INTEGER FUNCTION numchars(String)
      USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH, ERROR_control
      IMPLICIT NONE
! Argument
      CHARACTER(LEN=*), INTENT(IN) :: String
! Functions
      INTRINSIC INDEX, CHAR, LEN_TRIM
!***********************************************************************
      numchars = INDEX( String, CHAR(0) )
      IF ( numchars==0 ) numchars = INDEX( String, ' ' )
      numchars = numchars - 1
      IF ( numchars==-1 ) numchars = LEN_TRIM( String )
      IF ( numchars>MAXFILE_LENGTH ) THEN
        PRINT *, 'PRMS code error, string longer than:', MAXFILE_LENGTH, ' referenced'
        PRINT *, 'string length:', numchars, ' value: ', String
        PRINT *, 'Contact PRMS program support'
        ERROR STOP ERROR_control
      ENDIF
      END FUNCTION numchars

!***********************************************************************
! print_module
! print module version information to user's screen
!***********************************************************************
      SUBROUTINE print_module(Description, Modname, Versn)
      USE PRMS_CONSTANTS, ONLY: DEBUG_minimum
      USE PRMS_MODULE, ONLY: Print_debug, PRMS_output_unit
      IMPLICIT NONE
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Description, Modname, Versn
      ! Functions
      INTRINSIC TRIM, LEN_TRIM, MAX
      ! Local Variables
      INTEGER nvers, nmod, nblanks, nblanks2
      CHARACTER(LEN=24) :: blanks
      CHARACTER(LEN=68) :: string
!***********************************************************************
      IF ( Print_debug==DEBUG_minimum ) RETURN
      nvers = LEN_TRIM(Description)
      nblanks = MIN(32 - nvers, 32)
      nmod = LEN_TRIM(Modname)
      nblanks2 = MIN(20 - nmod, 20)
      blanks = ' '
      string = Description//blanks(:nblanks)//Modname//blanks(:nblanks2)//Versn
      PRINT '(A)', TRIM( string )
      WRITE ( PRMS_output_unit, '(A)' ) TRIM( string )
      END SUBROUTINE print_module

!***********************************************************************
! check restart file module order
!***********************************************************************
      SUBROUTINE check_restart(Modname, Restart_module)
      USE PRMS_CONSTANTS, ONLY: ERROR_restart
      IMPLICIT NONE
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Modname, Restart_module
!***********************************************************************
      IF ( Restart_module/=Modname ) THEN
        PRINT *, 'ERROR READING RESTART FILE, expecting module: ', Modname, ' found: ', Restart_module
        ERROR STOP ERROR_restart
      ENDIF
      END SUBROUTINE check_restart

!***********************************************************************
! check restart file dimensions order
!***********************************************************************
      SUBROUTINE check_restart_dimen(Dimen, Oldval, Newval, ierr)
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(IN) :: Oldval, Newval
      INTEGER, INTENT(INOUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN) :: Dimen
!***********************************************************************
      IF ( Oldval/=Newval ) THEN
        PRINT *, 'ERROR READING RESTART FILE, for dimension ', Dimen
        PRINT '(A,I0,A,I0,/)', '      restart value= ', Oldval, ' new value= ', Newval
        ierr = 1
      ENDIF
      END SUBROUTINE check_restart_dimen

!***********************************************************************
! Print date
!***********************************************************************
      SUBROUTINE print_date(Flag)
      USE PRMS_MODULE, ONLY: Nowyear, Nowmonth, Nowday
      USE PRMS_SET_TIME, ONLY: Nowhour, Nowminute
      IMPLICIT NONE
      ! Arguments
      INTEGER, INTENT(IN) :: Flag
!***********************************************************************
      IF ( Flag==1 ) THEN
        PRINT 9001, Nowyear, Nowmonth, Nowday, Nowhour, Nowminute
      ELSEIF ( Flag==0 ) THEN
        PRINT 9001, Nowyear, Nowmonth, Nowday
      ELSE
        WRITE ( Flag, 9001 ) Nowyear, Nowmonth, Nowday
      ENDIF
 9001 FORMAT ('    Date: ', I4, 2('/', I2.2), I3.2, ':', I2.2, /)
      END SUBROUTINE print_date

!***********************************************************************
!     Check parameter value limits
!***********************************************************************
      SUBROUTINE check_param_limits(Indx, Param, Param_value, Lower_val, Upper_val, Iret)
! Arguments
      INTEGER, INTENT(IN) :: Indx
      REAL, INTENT(IN) :: Param_value, Lower_val, Upper_val
      CHARACTER(LEN=*), INTENT(IN) :: Param
      INTEGER, INTENT(INOUT) :: Iret
!***********************************************************************
      IF ( Param_value<Lower_val .OR. Param_value>Upper_val ) THEN
        PRINT *, 'ERROR, bad value, parameter: ', Param
        PRINT *, '       value:  ', Param_value, '; array index:', Indx
        PRINT *, '       minimum:', Lower_val, '  ; maximum:', Upper_val
        PRINT *, ' '
        Iret = 1
      ENDIF
      END SUBROUTINE check_param_limits

!***********************************************************************
!     Check parameter value against dimension
!***********************************************************************
      SUBROUTINE checkdim_param_limits(Indx, Param, Dimen, Param_value, Lower_val, Upper_val, Iret)
! Arguments
      INTEGER, INTENT(IN) :: Indx, Param_value, Lower_val, Upper_val
      CHARACTER(LEN=*), INTENT(IN) :: Param, Dimen
      INTEGER, INTENT(INOUT) :: Iret
!***********************************************************************
      IF ( Param_value<Lower_val .OR. Param_value>Upper_val ) THEN
        PRINT *, 'ERROR, out-of-bounds value for bounded parameter: ', Param
        PRINT '(A,I0,A,I0)', '       value:  ', Param_value, '; array index: ', Indx
        PRINT '(A,I0,3A,I0,/)', '       minimum: ', Lower_val, '; maximum is dimension ', Dimen, ' = ', Upper_val
        Iret = 1
      ENDIF
      END SUBROUTINE checkdim_param_limits

!***********************************************************************
!     Check parameter value against bounded dimension
!***********************************************************************
      SUBROUTINE checkdim_bounded_limits(Param, Bound, Param_value, Num_values, Lower_val, Upper_val, Iret)
! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Param, Bound
      INTEGER, INTENT(IN) :: Num_values, Param_value(Num_values), Lower_val, Upper_val
      INTEGER, INTENT(OUT) :: Iret
! Local Variable
      INTEGER :: i
!***********************************************************************
      DO i = 1, Num_values
        IF ( Param_value(i)<Lower_val .OR. Param_value(i)>Upper_val ) THEN
          PRINT *, 'ERROR, out-of-bounds value for bounded parameter: ', Param
          PRINT '(A,I0,A,I0)', '       value:  ', Param_value(i), '; array index: ', i
          PRINT '(A,I0,3A,I0,/)', '       minimum: ', Lower_val, '; maximum is dimension ', Bound, ' = ', Upper_val
          Iret = 1
        ENDIF
      ENDDO
    END SUBROUTINE checkdim_bounded_limits

!***********************************************************************
!     Check parameter value < 0.0
!***********************************************************************
      SUBROUTINE check_param_zero(Indx, Param, Param_value, Iret)
! Arguments
      INTEGER, INTENT(IN) :: Indx
      REAL, INTENT(IN) :: Param_value
      CHARACTER(LEN=*), INTENT(IN) :: Param
      INTEGER, INTENT(INOUT) :: Iret
!***********************************************************************
      IF ( Param_value<0.0 ) THEN
        PRINT *, 'ERROR, value < 0.0 for parameter: ', Param
        PRINT *, '       value:', Param_value, '; HRU:', Indx
        PRINT *, ' '
        Iret = 1
      ENDIF
      END SUBROUTINE check_param_zero

!***********************************************************************
! Print ERROR message and stop
!***********************************************************************
      SUBROUTINE error_stop(Msg, Ierr)
      IMPLICIT NONE
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Msg
      INTEGER, INTENT(IN) :: Ierr
!***********************************************************************
      PRINT '(/,A,I0,2A,/)', 'ERROR ', Ierr, ', ', Msg
      ERROR STOP -1
      END SUBROUTINE error_stop
