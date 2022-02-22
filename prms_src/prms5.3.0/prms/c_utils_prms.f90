! utils_prms.f90 2020-07-30

module prms_utils
  interface
    module subroutine find_current_time(Iunit, Year, Month, Day, Iret, Cbh_binary_flag)
      integer, intent(IN) :: Iunit, Year, Month, Day, Cbh_binary_flag
      integer, intent(OUT) :: Iret
    end subroutine
  end interface

  interface
    module subroutine find_current_file_time(Iunit, Year, Month, Day, Year_file, Month_file, Day_file)
      integer, intent(IN) :: Iunit, Year, Month, Day
      integer, intent(OUT) :: Year_file, Month_file, Day_file
    end subroutine
  end interface

  interface
    module subroutine find_header_end(Iunit, Fname, Paramname, Iret, Cbh_flag, Cbh_binary_flag)
      integer, intent(IN) :: Cbh_flag, Cbh_binary_flag
      integer, intent(OUT) :: Iunit
      integer, intent(INOUT) :: Iret
      character(LEN=*), intent(IN) :: Fname, Paramname
    end subroutine
  end interface

  interface
    module subroutine is_eof(Iunit, Next_yr, Next_mo, Next_day)
      integer, intent(IN) :: Iunit
      integer, intent(OUT) :: Next_yr, Next_mo, Next_day
    end subroutine
  end interface

  interface
    integer module function get_ftnunit(Iunit)
      integer, intent(IN) :: Iunit
    end function
  end interface

  interface
    real module function f_to_c(Temp)
      real, intent(IN) :: Temp
    end function
  end interface

  interface
    real module function c_to_f(Temp)
      real, intent(IN) :: Temp
    end function
  end interface

  interface
    module subroutine write_integer_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
      integer, intent(IN) :: Iunit, Dimen
      integer, intent(IN) :: Values(Dimen)
      character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    end subroutine
  end interface

  interface
    module subroutine write_double_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
      !***********************************************************************
      implicit none
      ! Arguments
      integer, intent(IN) :: Iunit, Dimen
      double precision, intent(IN) :: Values(Dimen)
      character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    end subroutine
  end interface

  interface
    module subroutine write_real_param(Iunit, Parm_name, Dimen_name, Dimen, Values)
      integer, intent(IN) :: Iunit, Dimen
      real, intent(IN) :: Values(Dimen)
      character(LEN=*), intent(IN) :: Parm_name, Dimen_name
    end subroutine
  end interface

  interface
    module subroutine write_2D_double_param(Iunit, Parm_name, Dimen_name1, Dimen1, &
                                            Dimen_name2, Dimen2, Values)
      integer, intent(IN) :: Dimen1, Dimen2, Iunit
      double precision, intent(IN) :: Values(Dimen1, Dimen2)
      character(LEN=*), intent(IN) :: Parm_name
      character(LEN=*), intent(IN) :: Dimen_name1, Dimen_name2
    end subroutine
  end interface

  interface
    module subroutine write_2D_double_array_grid(Iunit, Parm_name, Dimen_name1, &
                                                 Dimen1, Dimen_name2, Dimen2, Values)
      integer, intent(IN) :: Iunit, Dimen1, Dimen2
      double precision, intent(IN) :: Values(Dimen1, Dimen2)
      character(LEN=*), intent(IN) :: Parm_name, Dimen_name1, Dimen_name2
    end subroutine
  end interface

  interface
    module subroutine version_check(Module_version, Length, Param_version)
      character(LEN=*), intent(IN) :: Module_version, Param_version
      integer, intent(IN) :: Length
    end subroutine
  end interface

  interface
    module subroutine read_error(Iflag, Name)
      integer, intent(IN) :: Iflag
      character(LEN=*), intent(IN) :: Name
    end subroutine
  end interface

  interface
    module subroutine module_error(Modname, Arg, Retcode)
      character(LEN=*), intent(IN) :: Modname, Arg
      integer, intent(IN) :: Retcode
    end subroutine
  end interface

  interface
    real module function sat_vapor_press_poly(Tempc)
      real, intent(IN) :: Tempc
    end function
  end interface

  interface
    real module function sat_vapor_press(Tempc)
      real, intent(IN) :: Tempc
    end function
  end interface

  interface
    integer module function leap_day(Year)
      integer, intent(IN) :: Year
    end function
  end interface

  interface
    module subroutine write_outfile(String)
      character(LEN=*), intent(IN) :: String
    end subroutine
  end interface

  interface
    integer module function julian_day(Date_type, Year_type)
      character(LEN=*), intent(IN) :: Date_type ! "start", "end", "now"
      character(LEN=*), intent(IN) :: Year_type ! "calendar", "solar", "water", "absolute"
    end function
  end interface

  interface
    integer module function compute_julday(Year, Month, Day)
      integer, intent(IN) :: Year, Month, Day
    end function
  end interface

  interface
    module subroutine compute_gregorian(Julday, Year, Month, Day)
      integer, intent(OUT) :: Year, Month, Day
      integer, intent(IN) :: Julday
    end subroutine
  end interface

  interface
    integer module function julday_in_year(Year, Month, Day)
      integer, intent(IN) :: Year, Month, Day
    end function
  end interface

  interface
    module subroutine PRMS_open_input_file(Iunit, Fname, Paramname, Ftype, Iret)
      integer, intent(IN) :: Ftype
      integer, intent(OUT) :: Iunit, Iret
      character(LEN=*), intent(IN) :: Fname, Paramname
    end subroutine
  end interface

  interface
    module subroutine PRMS_open_output_file(Iunit, Fname, Paramname, Ftype, Iret)
      integer, intent(IN) :: Ftype ! 0=text; 1=BINARY
      integer, intent(OUT) :: Iunit, Iret
      character(LEN=*), intent(IN) :: Fname, Paramname
    end subroutine
  end interface

  interface
    module subroutine PRMS_open_module_file(Iunit, Fname)
      integer, intent(OUT) :: Iunit
      character(LEN=*), intent(IN) :: Fname
    end subroutine
  end interface

  interface
    integer module function numchars(String)
      character(LEN=*), intent(IN) :: String
    end function
  end interface

  interface
    module subroutine print_module(Description, Modname, Versn)
      character(LEN=*), intent(IN) :: Description, Modname, Versn
    end subroutine
  end interface

  interface
    module subroutine check_restart(Modname, Restart_module)
      character(LEN=*), intent(IN) :: Modname, Restart_module
    end subroutine
  end interface

  interface
    module subroutine check_restart_dimen(Dimen, Oldval, Newval, ierr)
      integer, intent(IN) :: Oldval, Newval
      integer, intent(INOUT) :: ierr
      character(LEN=*), intent(IN) :: Dimen
    end subroutine
  end interface

  interface
    module subroutine print_date(Flag)
      integer, intent(IN) :: Flag
    end subroutine
  end interface

  interface
    module subroutine check_param_limits(Indx, Param, Param_value, Lower_val, Upper_val, Iret)
      integer, intent(IN) :: Indx
      real, intent(IN) :: Param_value, Lower_val, Upper_val
      character(LEN=*), intent(IN) :: Param
      integer, intent(INOUT) :: Iret
    end subroutine
  end interface

  interface
    module subroutine checkdim_param_limits(Indx, Param, Dimen, Param_value, Lower_val, Upper_val, Iret)
      integer, intent(IN) :: Indx, Param_value, Lower_val, Upper_val
      character(LEN=*), intent(IN) :: Param, Dimen
      integer, intent(INOUT) :: Iret
    end subroutine
  end interface

  interface
    module subroutine checkdim_bounded_limits(Param, Bound, Param_value, Num_values, Lower_val, Upper_val, Iret)
      character(LEN=*), intent(IN) :: Param, Bound
      integer, intent(IN) :: Num_values, Param_value(Num_values), Lower_val, Upper_val
      integer, intent(OUT) :: Iret
    end subroutine
  end interface

  interface
    module subroutine check_param_zero(Indx, Param, Param_value, Iret)
      integer, intent(IN) :: Indx
      real, intent(IN) :: Param_value
      character(LEN=*), intent(IN) :: Param
      integer, intent(INOUT) :: Iret
    end subroutine
  end interface

  interface
    module subroutine error_stop(Msg, Ierr)
      character(LEN=*), intent(IN) :: Msg
      integer, intent(IN) :: Ierr
    end subroutine
  end interface

end module
