! utils_prms.f90 2020-07-30

submodule(prms_utils) sm_utils_prms

contains

!***********************************************************************
!     Read CBH File to current time
!***********************************************************************
  module subroutine find_current_time(Iunit, Year, Month, Day, Iret, Cbh_binary_flag)
    implicit none
    ! Argument
    integer, intent(IN) :: Iunit, Year, Month, Day, Cbh_binary_flag
    integer, intent(OUT) :: Iret
    ! Local Variables
    integer :: yr, mo, dy
    character(LEN=20), parameter :: fmt1 = '(A, I5, 2("/",I2.2))'
    !***********************************************************************
    Iret = 0
    do
      if (Cbh_binary_flag == 0) then
        read (Iunit, *, IOSTAT=Iret) yr, mo, dy
      else
        read (Iunit, IOSTAT=Iret) yr, mo, dy
      end if
      if (Iret == -1) print fmt1, 'ERROR, end-of-file found reading input file for date:', Year, Month, Day
      if (Iret /= 0) return
      if (yr == Year .and. mo == Month .and. dy == Day) exit
    end do
    backspace Iunit
  end subroutine find_current_time

!***********************************************************************
!     Read File dynamic paramter file to current time
!***********************************************************************
  module subroutine find_current_file_time(Iunit, Year, Month, Day, Year_file, Month_file, Day_file)
    implicit none
    ! Argument
    integer, intent(IN) :: Iunit, Year, Month, Day
    integer, intent(OUT) :: Year_file, Month_file, Day_file
    ! Local Variables
    integer :: i, ios
    !***********************************************************************
    ! find first value for simulation time period
    read (Iunit, *, IOSTAT=ios) Year_file, Month_file, Day_file
    if (ios /= 0) then
      Year_file = 0
      Month_file = 0
      Day_file = 0
      return
    end if
    if (Year_file < Year) then
      i = 0
      do while (i == 0)
        read (Iunit, *, IOSTAT=ios) Year_file, Month_file, Day_file
        if (ios /= 0) then
          Year_file = 0
          Month_file = 0
          Day_file = 0
          return
        end if
        if (Year_file >= Year) i = 1
      end do
    end if
    if (Year_file == Year) then
      if (Month_file < Month) then
        i = 0
        do while (i == 0)
          read (Iunit, *, IOSTAT=ios) Year_file, Month_file, Day_file
          if (ios /= 0) then
            Year_file = 0
            Month_file = 0
            Day_file = 0
            return
          end if
          if (Month_file >= Month .or. Year_file /= Year) i = 1
        end do
      end if
      if (Year_file == Year .and. Month_file == Month) then
        if (Day_file < Day) then
          i = 0
          do while (i == 0)
            read (Iunit, *, IOSTAT=ios) Year_file, Month_file, Day_file
            if (ios /= 0) then
              Year_file = 0
              Month_file = 0
              Day_file = 0
              return
            end if
            if (Day_file >= Day) i = 1
          end do
        end if
      end if
    end if
    backspace Iunit
  end subroutine find_current_file_time

!***********************************************************************
!     Read File to line before data starts in file
!***********************************************************************
  module subroutine find_header_end(Iunit, Fname, Paramname, Iret, Cbh_flag, Cbh_binary_flag)
    use PRMS_CONSTANTS, only: DEBUG_less
    use PRMS_MODULE, only: Nhru, Orad_flag, Print_debug
    implicit none
    ! Argument
    integer, intent(IN) :: Cbh_flag, Cbh_binary_flag
    integer, intent(OUT) :: Iunit
    integer, intent(INOUT) :: Iret
    character(LEN=*), intent(IN) :: Fname, Paramname
    ! Functions
    intrinsic :: trim

    ! Local Variables
    integer :: i, ios, dim
    character(LEN=4) :: dum
    character(LEN=80) :: dum2
    !***********************************************************************
    if (Iret /= 2) Iret = 0
    Iunit = get_ftnunit(7777)
    if (Cbh_binary_flag == 0) then
      open (Iunit, FILE=trim(Fname), STATUS='OLD', IOSTAT=ios)
    else
      open (Iunit, FILE=trim(Fname), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=ios)
    end if
    if (ios /= 0) then
      if (Iret == 2) then ! this signals climate_hru to ignore the Humidity CBH file, could add other files
        Iret = 0
        if (Print_debug > DEBUG_less) &
   &         write (*, '(/,A,/,A,/,A)') 'WARNING, optional CBH file not found, will use associated parameter values'
      else
        write (*, '(/,A,/,A,/,A)') 'ERROR reading file:', Fname, 'check to be sure the input file exists'
        Iret = 1
      end if
    else
      ! read to line before data starts in each file
      i = 0
      do while (i == 0)
        if (Cbh_binary_flag == 0) then
          read (Iunit, FMT='(A4)', IOSTAT=ios) dum
        else
          read (Iunit, IOSTAT=ios) dum2
          read (dum2, '(A4)') dum
        end if
        if (ios /= 0) then
          write (*, '(/,A,/,A,/,A)') 'ERROR reading file:', Fname, 'check to be sure the input file is in correct format'
          Iret = 1
          exit
        elseif (dum == '####') then
          if (Cbh_flag == 0) exit
          backspace Iunit
          backspace Iunit
          if (Orad_flag == 1 .and. Paramname(:5) == 'swrad') backspace Iunit ! backspace again as swrad CBH file contains orad as last column
          if (Cbh_binary_flag == 0) then
            read (Iunit, *, IOSTAT=ios) dum, dim
          else
            read (Iunit, IOSTAT=ios) dum2
            read (dum2, *) dum, dim
          end if
          if (ios /= 0) then
            write (*, '(/,A,/,A,/,A)') 'ERROR reading file:', Fname, 'check to be sure dimension line is in correct format'
            Iret = 1
            exit
          end if
          if (dim /= Nhru) then
            print '(/,2(A,I0))', '***CBH file dimension incorrect*** nhru= ', Nhru, ' CBH dimension= ', dim, ' File: '//Fname
            print *, 'ERROR: update Control File with correct CBH files'
            Iret = 1
            exit
          end if
          if (Cbh_binary_flag == 0) then
            read (Iunit, FMT='(A4)', IOSTAT=ios) dum
          else
            read (Iunit, IOSTAT=ios) dum
          end if
          if (ios /= 0) then
            write (*, '(/,A,/,A,/)') 'ERROR reading file:', Fname
            Iret = 1
            exit
          end if
          if (Orad_flag == 1 .and. Paramname(:5) == 'swrad') read (Iunit, FMT='(A4)') dum ! read again as swrad CBH file contains orad as last column
          i = 1
        end if
      end do
    end if

  end subroutine find_header_end

!**********************
! Check for end of file
!**********************
  module subroutine is_eof(Iunit, Next_yr, Next_mo, Next_day)
    implicit none
    ! Arguments
    integer, intent(IN) :: Iunit
    integer, intent(OUT) :: Next_yr, Next_mo, Next_day
    ! Local Variables
    integer :: ios, i
    character(LEN=80) :: dum
    !*******************************************************************************
    Next_yr = 0
    Next_mo = 0
    Next_day = 0
    i = 0
    do while (i == 0)
      read (Iunit, '(A)', iostat=ios) dum
      if (ios /= 0) return
      if (dum(:4) == '####') cycle
      if (dum(:2) /= '//') i = 1
    end do
    read (dum, *, iostat=ios) Next_yr, Next_mo, Next_day
    if (ios /= 0) then
      Next_yr = 0
      Next_mo = 0
      Next_day = 0
    else
      backspace Iunit
    end if
  end subroutine is_eof

!***********************************************************************
!     Determine an unopened FORTRAN File Unit
!***********************************************************************
  integer module function get_ftnunit(Iunit)
    implicit none
    ! Argument
    integer, intent(IN) :: Iunit
    ! Local Variables
    integer :: good_unit
    logical :: opend
    !***********************************************************************
    good_unit = Iunit
    opend = .true.
    do while (opend)
      good_unit = good_unit + 1
      inquire (UNIT=good_unit, OPENED=opend)
    end do
    get_ftnunit = good_unit
  end function get_ftnunit

!***********************************************************************
! Convert Fahrenheit to Celsius
!***********************************************************************
  real module function f_to_c(Temp)
    implicit none
    ! Arguments
    real, intent(IN) :: Temp
    !***********************************************************************
    f_to_c = (Temp - 32.0) / 1.8
  end function f_to_c

!***********************************************************************
! Convert Celsius to Fahrenheit
!***********************************************************************
  real module function c_to_f(Temp)
    implicit none
    ! Arguments
    real, intent(IN) :: Temp
    !***********************************************************************
    c_to_f = Temp * 1.8 + 32.0
  end function c_to_f

!***********************************************************************
  module subroutine write_integer_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
    implicit none
    ! Arguments
    integer, intent(IN) :: Iunit, Dimen
    integer, intent(IN) :: Values(Dimen)
    character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    ! Local Variables
    integer i
    character(LEN=48), parameter :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, /, "1")'
    !***********************************************************************
    write (Iunit, fmt1) Parm_name, Dimen_name, Dimen
    do i = 1, Dimen
      write (Iunit, '(I0)') Values(i)
    end do
  end subroutine write_integer_param

!***********************************************************************
  module subroutine write_real_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
    implicit none
    ! Arguments
    integer, intent(IN) :: Iunit, Dimen
    real, intent(IN) :: Values(Dimen)
    character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    ! Local Variables
    integer i
    character(LEN=48), parameter :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, /, "2")'
    !***********************************************************************
    write (Iunit, fmt1) Parm_name, Dimen_name, Dimen
    do i = 1, Dimen
      write (Iunit, *) Values(i)
    end do
  end subroutine write_real_param

!***********************************************************************
  module subroutine write_double_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
!***********************************************************************
    implicit none
    ! Arguments
    integer, intent(IN) :: Iunit, Dimen
    double precision, intent(IN) :: Values(Dimen)
    character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    ! Local Variables
    integer i
    character(LEN=40), parameter :: fmt1 = '("####", /, A, /, "1", /, A, /, I0, "3")'
    !***********************************************************************
    write (Iunit, fmt1) Parm_name, Dimen_name, Dimen
    do i = 1, Dimen
      write (Iunit, *) Values(i)
    end do
  end subroutine write_double_param

!***********************************************************************
  module subroutine write_2D_double_param(Iunit, Parm_name, Dimen_name1, Dimen1, Dimen_name2, Dimen2, Values)
!***********************************************************************
    implicit none
    ! Arguments
    integer, intent(IN) :: Dimen1, Dimen2, Iunit
    double precision, intent(IN) :: Values(Dimen1, Dimen2)
    character(LEN=*), intent(IN) :: Parm_name
    character(LEN=*), intent(IN) :: Dimen_name1, Dimen_name2
    ! Local Variables
    integer i, j
    character(LEN=46), parameter :: fmt1 = '("####", /, A, /, "2", /, A, /, A, /, I0, "3")'
    !***********************************************************************
    write (Iunit, fmt1) Parm_name, Dimen_name1, Dimen_name2, Dimen1 * Dimen2
    do i = 1, Dimen2
      do j = 1, Dimen1
        write (Iunit, *) Values(j, i)
      end do
    end do
  end subroutine write_2D_double_param

!***********************************************************************
  module subroutine write_2D_double_array_grid(Iunit, Parm_name, Dimen_name1, Dimen1, Dimen_name2, Dimen2, Values)
!***********************************************************************
    implicit none
    ! Arguments
    integer, intent(IN) :: Iunit, Dimen1, Dimen2
    double precision, intent(IN) :: Values(Dimen1, Dimen2)
    character(LEN=*), intent(IN) :: Parm_name, Dimen_name1, Dimen_name2
    ! Local Variables
    integer i, j
    character(LEN=12) :: fmt
    !***********************************************************************
    write (Iunit, 9001) Parm_name, Dimen_name1, Dimen_name2, Dimen1 * Dimen2
    write (fmt, 9002) Dimen2
    do i = 1, Dimen2
      write (Iunit, fmt) (Values(j, i), j=1, Dimen1)
    end do

 9001 format('####', /, A, /, '2', /, A, /, A, /, I0, /, '3')
 9002 format('(', I5, 'F10.5)')
  end subroutine write_2D_double_array_grid

!**********************************************************************
!     Version Check
!**********************************************************************
  module subroutine version_check(Module_version, Length, Param_version)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Module_version, Param_version
    integer, intent(IN) :: Length
    !**********************************************************************
    if (Module_version(13:Length + 12) /= Param_version(:Length)) then
      print 9001, Module_version(13:Length + 12), Param_version(:Length)
      print *, 'Enter return to continue'
      read (*, *)
    end if
9001 format('WARNING, module versions are not identical', /, &
      &        'Executable version: ', A, /, &
      &        'Parameter File version: ', A,/)
  end subroutine version_check

!**********************************************************************
!     Parameter or Variable delcare or read error
!**********************************************************************
  module subroutine read_error(Iflag, Name)
    use PRMS_CONSTANTS, only: ERROR_decl_get
    implicit none
    ! Arguments
    integer, intent(IN) :: Iflag
    character(LEN=*), intent(IN) :: Name
    !**********************************************************************
    print '(/,A,/)', 'Due to error condition simulation halted'
    if (Iflag == 1) then
      print *, 'Declare error for parameter: ', Name
    elseif (Iflag == 2) then
      print *, 'Get error for parameter: ', Name
    elseif (Iflag == 3) then
      print *, 'Declare error for variable: ', Name
    elseif (Iflag == 4) then
      print *, 'Get error for variable: ', Name
    elseif (Iflag == 5) then
      print *, 'Read error for control parameter: ', Name
    elseif (Iflag == 6) then
      print *, 'Read error for dimension parameter: ', Name
    elseif (Iflag == 7) then
      print *, 'Declare error for dimension parameter: ', Name
    elseif (Iflag == 8) then
      print *, 'Declare error for Data File variable: ', Name
    elseif (Iflag == 9) then
      print *, 'Read error for Data File variable: ', Name
    elseif (Iflag == 10) then
      print *, 'Open error of Control File ', Name
    elseif (Iflag == 11) then
      print *, 'Read error of Parameter File ', Name
    elseif (Iflag == 12) then
      print *, 'Read error of Control File ', Name
    elseif (Iflag == 13) then
      print *, 'Read error of Data File ', Name
    elseif (Iflag == 14) then
      print *, 'Control parameter not found: ', Name
    elseif (Iflag == 15) then
      print *, 'ERROR, control ', Name, ' expected and is not available in PRMS'
    elseif (Iflag == 16) then
      print *, 'ERROR, declared parameter ', Name
    end if
    ERROR stop ERROR_decl_get
  end subroutine read_error

!**********************************************************************
!     Module error
!**********************************************************************
  module subroutine module_error(Modname, Arg, Retcode)
    use PRMS_CONSTANTS, only: ERROR_module
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Arg
    integer, intent(IN) :: Retcode
    !**********************************************************************
    print 9001, Modname, Arg, Retcode
    ERROR stop ERROR_module
 9001 format('ERROR in ', A, ' module, arg = ', A, /, 'Return val = ', I0)
  end subroutine module_error

!***********************************************************************
! Compute saturation vapor pressure over water in millibars
! 6th order Polynominal method (Flatau et. all., 1992) valid: -50 to 50C
! 1 kPa = 10 millibars
! Flatau, P.j., Walko, R.L., Cotton, W.R., 1992, Polynomial Fits to
!   saturation vapor pressure: Jornal of Applied Meteorology, v. 31, p. 1507-1513
!***********************************************************************
  real module function sat_vapor_press_poly(Tempc)
    implicit none
    ! Arguments
    real, intent(IN) :: Tempc
    !***********************************************************************
    sat_vapor_press_poly = 6.11176750 + 0.443986062 * Tempc &
    &                       + 0.0143053301 * Tempc**2 &
    &                       + 0.265027242E-03 * Tempc**3 &
    &                       + 0.302246994E-05 * Tempc**4 &
    &                       + 0.203886313E-07 * Tempc**5 &
    &                       + 0.638780966E-10 * Tempc**6
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
  end function sat_vapor_press_poly

!***********************************************************************
! Compute saturation vapor pressure over water
! Irmak and others (2012), equation 12
!***********************************************************************
  real module function sat_vapor_press(Tempc)
    implicit none
    ! Arguments
    real, intent(IN) :: Tempc
    !***********************************************************************
    sat_vapor_press = 6.1078 * exp((17.26939 * Tempc) / (237.3 + Tempc))
  end function sat_vapor_press

!***********************************************************************
! leap_day - is the year a leap year: (1=yes; 0=no)
!***********************************************************************
  integer module function leap_day(Year)
    implicit none
    ! Arguments
    integer, intent(IN) :: Year
    ! Functions
    intrinsic MOD
    !***********************************************************************
    leap_day = 0
    ! Check if leapyear - Start by identifying all years not divisible by 4
    if (mod(Year, 4) == 0) then
      leap_day = 1
      if (mod(Year, 100) == 0) then
        if (mod(Year, 400) /= 0) leap_day = 0
      end if
    end if
  end function leap_day

!***********************************************************************
! write_outfile - print to model output file
!***********************************************************************
  module subroutine write_outfile(String)
    use PRMS_CONSTANTS, only: DEBUG_minimum
    use PRMS_MODULE, only: PRMS_output_unit, Print_debug
    implicit none
    ! Functions
    intrinsic :: LEN_TRIM
    ! Arguments
    character(LEN=*), intent(IN) :: String
    ! Local variable
    integer nchars
    !***********************************************************************
    if (Print_debug == DEBUG_minimum) return
    nchars = len_trim(String)
    if (nchars > 0) then
      write (PRMS_output_unit, '(A)') String(:nchars)
    else
      write (PRMS_output_unit, '(/)')
    end if
  end subroutine write_outfile

!***********************************************************************
! julian_day
! computes the Julian date given a Gregorian calendar date
! (Year, Month, Day) relative to: calendar (Jan 1),
! solar (12/22 in Northern; 6/21 in Southern) and
! water year (10/1 in Northern; 4/1 in Southern) start dates.
! The Julian day starts at noon of the Gregorian day and
! extends to noon the next Gregorian day.
!***********************************************************************
  integer module function julian_day(Date_type, Year_type)
    use PRMS_CONSTANTS, only: NORTHERN, ERROR_time
    use PRMS_MODULE, only: Starttime, Endtime
    use PRMS_BASIN, only: Hemisphere
    use PRMS_SET_TIME, only: Nowtime
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Date_type ! "start", "end", "now"
    character(LEN=*), intent(IN) :: Year_type ! "calendar", "solar", "water", "absolute"

    ! Local Variables
    integer :: reftime_year, reftime_month, reftime_day, time_array(6)
    integer :: year, month, day, absolute_julday, relative_julday, length, found
    !***********************************************************************
    if (Date_type(:3) == 'end') then
      time_array = Endtime
    elseif (Date_type(:3) == 'now') then
      time_array = Nowtime
    elseif (Date_type(:5) == 'start') then
      time_array = Starttime
    else
      print *, 'ERROR, invalid argument to compute Julian Day: ', Date_type
      ERROR stop ERROR_time
    end if
    year = time_array(1)
    month = time_array(2)
    day = time_array(3)

    found = 0
    length = len(Year_type)
    ! set reftime depending on type arg
    if (length > 4) then
      if (Year_type(:5) == 'solar') then
        found = 1
        if (Hemisphere == NORTHERN) then
          if (month == 12 .and. day > 21) then
            reftime_year = year
          else
            reftime_year = year - 1
          end if
          reftime_month = 12
          reftime_day = 21
        else ! Southern
          if (month == 6 .and. day > 20) then
            reftime_year = year;
          else
            reftime_year = year - 1
          end if
          reftime_month = 6
          reftime_day = 20
        end if
      elseif (Year_type(:5) == 'water') then
        found = 1
        if (Hemisphere == NORTHERN) then
          if (month > 9) then
            reftime_year = year
          else
            reftime_year = year - 1
          end if
          reftime_month = 9
          reftime_day = 30
        else ! Southern
          if (month > 3) then
            reftime_year = year
          else
            reftime_year = year - 1
          end if
          reftime_month = 3
          reftime_day = 31
        end if
      end if
    end if
    if (found == 0 .and. length > 5) then
      if (Year_type(:6) == 'spring') then
        found = 1
        if (Hemisphere == NORTHERN) then
          if (month > 3 .or. (month == 3 .and. day > 20)) then
            reftime_year = year
          else
            reftime_year = year - 1
          end if
          reftime_month = 3
          reftime_day = 20
        else ! Southern
          if (month > 9 .or. (month == 9 .and. day > 22)) then
            reftime_year = year
          else
            reftime_year = year - 1
          end if
          reftime_month = 9
          reftime_day = 22
        end if
      end if
    end if
    if (found == 0 .and. length > 7) then
      if (Year_type(:8) == 'calendar') then
        found = 1
        reftime_year = year - 1
        reftime_month = 12
        reftime_day = 31
      end if
    end if
    if (found == 0) then
      print *, 'ERROR, invalid year type argument to compute Julian Day: ', Year_type
      ERROR stop ERROR_time
    end if

    ! set actual Julian Day
    absolute_julday = compute_julday(year, month, day)

    relative_julday = 0
    if (length == 8) then
      if (Year_type(:8) == 'calendar') relative_julday = compute_julday(reftime_year, reftime_month, reftime_day)
    else
      relative_julday = compute_julday(reftime_year, reftime_month, reftime_day)
    end if
    julian_day = absolute_julday - relative_julday

  end function julian_day

!***********************************************************************
! compute_julday
! computes the Julian Day given a Gregorian calendar date
!***********************************************************************
  integer module function compute_julday(Year, Month, Day)
    implicit none
    ! Arguments
    integer, intent(IN) :: Year, Month, Day
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
    compute_julday = Day - 32075 + 1461 * (Year + 4800 + (Month - 14) / 12) / 4 + 367 * (Month - 2 - (Month - 14) / 12 * 12) &
    &    / 12 - 3 * ((Year + 4900 + (Month - 14) / 12) / 100) / 4

  end function compute_julday

!***********************************************************************
! compute_gregorian
! computes the Gregorian calendar date given the Julian Day
!***********************************************************************
  module subroutine compute_gregorian(Julday, Year, Month, Day)
    !USE PRMS_CONSTANTS, ONLY: DAYS_YR
    implicit none
    ! Arguments
    integer, intent(OUT) :: Year, Month, Day
    integer, intent(IN) :: Julday
    ! Functions
    intrinsic FLOOR, NINT
    ! Local Variables
    integer m, n
    !***********************************************************************
    m = Julday + 68569
    n = 4 * m / 146097
    m = m - (146097 * n + 3) / 4
    Year = 4000 * (m + 1) / 1461001
    m = m - 1461 * Year / 4 + 31
    Month = 80 * m / 2447
    Day = m - 2447 * Month / 80
    m = Month / 11
    Month = Month + 2 - 12 * m
    Year = 100 * (n - 49) + Year + m

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
  end subroutine compute_gregorian

!***********************************************************************
! julday_in_year
! computes the Julian Day of a date
!***********************************************************************
  integer module function julday_in_year(Year, Month, Day)
    implicit none
    ! Arguments
    integer, intent(IN) :: Year, Month, Day

    ! Local Variables
    integer daypmo(12), i
    data daypmo/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
    !***********************************************************************
    daypmo(2) = 28
    if (leap_day(Year) == 1) daypmo(2) = 29
    julday_in_year = Day
    do i = 1, Month - 1
      julday_in_year = julday_in_year + daypmo(i)
    end do
  end function julday_in_year

!***********************************************************************
!     Open PRMS input File and assign unit number
!***********************************************************************
  module subroutine PRMS_open_input_file(Iunit, Fname, Paramname, Ftype, Iret)
    implicit none
    ! Argument
    integer, intent(IN) :: Ftype
    integer, intent(OUT) :: Iunit, Iret
    character(LEN=*), intent(IN) :: Fname, Paramname
    ! Local Variables
    integer :: ios, nchars
    !***********************************************************************
    Iret = 0
    Iunit = get_ftnunit(7777)
    nchars = numchars(Fname)
    if (Ftype == 0) then
      open (Iunit, FILE=Fname(:nchars), STATUS='OLD', IOSTAT=ios)
    else
      OPEN ( Iunit, FILE=Fname(:nchars), STATUS='OLD', FORM='UNFORMATTED', IOSTAT=ios, ACCESS='STREAM' )
    end if
    if (ios /= 0) then
      write (*, '(/,2A,/,A,/,2A,/)') 'ERROR opening input file: ', Fname(:nchars), &
     &                               'check to be sure the input file exists', &
     &                               'file specified by control parameter: ', Paramname
      Iret = 1
    end if
  end subroutine PRMS_open_input_file

!***********************************************************************
!     Open PRMS output file and assign unit number
!***********************************************************************
  module subroutine PRMS_open_output_file(Iunit, Fname, Paramname, Ftype, Iret)
    implicit none
    ! Argument
    integer, intent(IN) :: Ftype ! 0=text; 1=BINARY
    integer, intent(OUT) :: Iunit, Iret
    character(LEN=*), intent(IN) :: Fname, Paramname
    ! Local Variables
    integer :: ios, nchars
    !***********************************************************************
    Iret = 0
    Iunit = get_ftnunit(8888)
    nchars = numchars(Fname)

    if (Ftype == 0) then
      open (Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios)
    else
      OPEN ( Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios, FORM='UNFORMATTED', ACCESS='STREAM' )
    end if

    if (ios /= 0) then
      write (*, '(/,A,/,A,/)') 'ERROR opening output file:', Fname(:nchars), &
      &                             'check to be sure the pathname is valid and the file is not open'
      write (*, '(2A,/)') 'file specified by control parameter: ', Paramname
      Iret = 1
    end if

  end subroutine PRMS_open_output_file

!***********************************************************************
!     Open PRMS module output file and assign unit number
!***********************************************************************
  module subroutine PRMS_open_module_file(Iunit, Fname)
    use PRMS_CONSTANTS, only: ERROR_open_in
    implicit none
    ! Argument
    integer, intent(OUT) :: Iunit
    character(LEN=*), intent(IN) :: Fname
    ! Local Variables
    integer :: ios, nchars
    !***********************************************************************
      Iunit = get_ftnunit(8888)
      nchars = numchars(Fname)
    open (Iunit, FILE=Fname(:nchars), STATUS='REPLACE', IOSTAT=ios)
    if (ios /= 0) then
      write (*, '(/,A,/,A,/)') 'ERROR opening water balance output file:', Fname(:nchars), &
     &                         'check to be sure the pathname is valid and the file is not open'
      ERROR stop ERROR_open_in
    end if
  end subroutine PRMS_open_module_file

!***********************************************************************
!     Determine number of characters in a string
!***********************************************************************
  integer module function numchars(String)
    use PRMS_CONSTANTS, only: MAXFILE_LENGTH, ERROR_control
    implicit none
    ! Argument
    character(LEN=*), intent(IN) :: String
    ! Functions
    intrinsic INDEX, CHAR, LEN_TRIM
    !***********************************************************************
    numchars = index(String, char(0))
    if (numchars == 0) numchars = index(String, ' ')
    numchars = numchars - 1
    if (numchars == -1) numchars = len_trim(String)
    if (numchars > MAXFILE_LENGTH) then
      print *, 'PRMS code error, string longer than:', MAXFILE_LENGTH, ' referenced'
      print *, 'string length:', numchars, ' value: ', String
      print *, 'Contact PRMS program support'
      ERROR stop ERROR_control
    end if
  end function numchars

!***********************************************************************
! print_module
! print module version information to user's screen
!***********************************************************************
  module subroutine print_module(Description, Modname, Versn)
    use PRMS_CONSTANTS, only: DEBUG_minimum
    use PRMS_MODULE, only: Print_debug, PRMS_output_unit
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Description, Modname, Versn
    ! Functions
    intrinsic TRIM, LEN_TRIM, MAX
    ! Local Variables
    integer nvers, nmod, nblanks, nblanks2
    character(LEN=24) :: blanks
    character(LEN=68) :: string
    !***********************************************************************
    if (Print_debug == DEBUG_minimum) return
    nvers = len_trim(Description)
    nblanks = min(32 - nvers, 32)
    nmod = len_trim(Modname)
    nblanks2 = min(20 - nmod, 20)
    blanks = ' '
    string = Description//blanks(:nblanks)//Modname//blanks(:nblanks2)//Versn
    print '(A)', trim(string)
    write (PRMS_output_unit, '(A)') trim(string)
  end subroutine print_module

!***********************************************************************
! check restart file module order
!***********************************************************************
  module subroutine check_restart(Modname, Restart_module)
    use PRMS_CONSTANTS, only: ERROR_restart
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Restart_module
    !***********************************************************************
    if (Restart_module /= Modname) then
      print *, 'ERROR READING RESTART FILE, expecting module: ', Modname, ' found: ', Restart_module
      ERROR stop ERROR_restart
    end if
  end subroutine check_restart

!***********************************************************************
! check restart file dimensions order
!***********************************************************************
  module subroutine check_restart_dimen(Dimen, Oldval, Newval, ierr)
    implicit none
    ! Arguments
    integer, intent(IN) :: Oldval, Newval
    integer, intent(INOUT) :: ierr
    character(LEN=*), intent(IN) :: Dimen
    !***********************************************************************
    if (Oldval /= Newval) then
      print *, 'ERROR READING RESTART FILE, for dimension ', Dimen
      print '(A,I0,A,I0,/)', '      restart value= ', Oldval, ' new value= ', Newval
      ierr = 1
    end if
  end subroutine check_restart_dimen

!***********************************************************************
! Print date
!***********************************************************************
  module subroutine print_date(Flag)
    use PRMS_MODULE, only: Nowyear, Nowmonth, Nowday
    use PRMS_SET_TIME, only: Nowhour, Nowminute
    implicit none
    ! Arguments
    integer, intent(IN) :: Flag
    !***********************************************************************
    if (Flag == 1) then
      print 9001, Nowyear, Nowmonth, Nowday, Nowhour, Nowminute
    elseif (Flag == 0) then
      print 9001, Nowyear, Nowmonth, Nowday
    else
      write (Flag, 9001) Nowyear, Nowmonth, Nowday
    end if
9001 format('    Date: ', I4, 2('/', I2.2), I3.2, ':', I2.2,/)
  end subroutine print_date

!***********************************************************************
!     Check parameter value limits
!***********************************************************************
  module subroutine check_param_limits(Indx, Param, Param_value, Lower_val, Upper_val, Iret)
    ! Arguments
    integer, intent(IN) :: Indx
    real, intent(IN) :: Param_value, Lower_val, Upper_val
    character(LEN=*), intent(IN) :: Param
    integer, intent(INOUT) :: Iret
    !***********************************************************************
    if (Param_value < Lower_val .or. Param_value > Upper_val) then
      print *, 'ERROR, bad value, parameter: ', Param
      print *, '       value:  ', Param_value, '; array index:', Indx
      print *, '       minimum:', Lower_val, '  ; maximum:', Upper_val
      print *, ' '
      Iret = 1
    end if
  end subroutine check_param_limits

!***********************************************************************
!     Check parameter value against dimension
!***********************************************************************
  module subroutine checkdim_param_limits(Indx, Param, Dimen, Param_value, Lower_val, Upper_val, Iret)
    ! Arguments
    integer, intent(IN) :: Indx, Param_value, Lower_val, Upper_val
    character(LEN=*), intent(IN) :: Param, Dimen
    integer, intent(INOUT) :: Iret
    !***********************************************************************
    if (Param_value < Lower_val .or. Param_value > Upper_val) then
      print *, 'ERROR, out-of-bounds value for bounded parameter: ', Param
      print '(A,I0,A,I0)', '       value:  ', Param_value, '; array index: ', Indx
      print '(A,I0,3A,I0,/)', '       minimum: ', Lower_val, '; maximum is dimension ', Dimen, ' = ', Upper_val
      Iret = 1
    end if
  end subroutine checkdim_param_limits

!***********************************************************************
!     Check parameter value against bounded dimension
!***********************************************************************
  module subroutine checkdim_bounded_limits(Param, Bound, Param_value, Num_values, Lower_val, Upper_val, Iret)
    ! Arguments
    character(LEN=*), intent(IN) :: Param, Bound
    integer, intent(IN) :: Num_values, Param_value(Num_values), Lower_val, Upper_val
    integer, intent(OUT) :: Iret
    ! Local Variable
    integer :: i
    !***********************************************************************
    do i = 1, Num_values
      if (Param_value(i) < Lower_val .or. Param_value(i) > Upper_val) then
        print *, 'ERROR, out-of-bounds value for bounded parameter: ', Param
        print '(A,I0,A,I0)', '       value:  ', Param_value(i), '; array index: ', i
        print '(A,I0,3A,I0,/)', '       minimum: ', Lower_val, '; maximum is dimension ', Bound, ' = ', Upper_val
        Iret = 1
      end if
    end do
  end subroutine checkdim_bounded_limits

!***********************************************************************
!     Check parameter value < 0.0
!***********************************************************************
  module subroutine check_param_zero(Indx, Param, Param_value, Iret)
    ! Arguments
    integer, intent(IN) :: Indx
    real, intent(IN) :: Param_value
    character(LEN=*), intent(IN) :: Param
    integer, intent(INOUT) :: Iret
    !***********************************************************************
    if (Param_value < 0.0) then
      print *, 'ERROR, value < 0.0 for parameter: ', Param
      print *, '       value:', Param_value, '; HRU:', Indx
      print *, ' '
      Iret = 1
    end if
  end subroutine check_param_zero

!***********************************************************************
! Print ERROR message and stop
!***********************************************************************
  module subroutine error_stop(Msg, Ierr)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Msg
    integer, intent(IN) :: Ierr
    !***********************************************************************
    print '(/,A,I0,2A,/)', 'ERROR ', Ierr, ', ', Msg
    ERROR stop - 1
  end subroutine error_stop
end submodule
