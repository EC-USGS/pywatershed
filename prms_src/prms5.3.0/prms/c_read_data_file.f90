!***********************************************************************
! Read PRMS Data File
!***********************************************************************
module PRMS_DATA_FILE
  implicit none

  integer, save :: Num_datafile_types, Num_datafile_columns, Datafile_unit
  character(LEN=16), allocatable, save :: Data_varname(:)
  integer, allocatable, save :: Data_varnum(:)
  real, allocatable, save :: Data_line_values(:)
  character(LEN=1), allocatable :: data_transfer(:)

  interface
    module subroutine read_prms_data_file()
    end subroutine

    module subroutine read_data_line()
    end subroutine

    module subroutine check_data_variables(Varname, Numvalues, Values, Iflag, Iret)
      character(LEN=*), intent(IN) :: Varname
      integer, intent(IN) :: Numvalues, Iflag
      integer, intent(OUT) :: Iret
      real, intent(IN) :: Values(Numvalues)
    end subroutine

    module subroutine read_data_file_line(Iret)
      integer, intent(OUT) :: Iret
    end subroutine
  end interface
end module PRMS_DATA_FILE
