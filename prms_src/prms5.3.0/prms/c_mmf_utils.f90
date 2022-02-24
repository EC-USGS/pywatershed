!***********************************************************************
! DONE:
!  getdim, declfix, declmodule, decldim, declparam, getparam
!  control_string, control_integer, control_string_array, getvartype
!***********************************************************************
module PRMS_MMFAPI
  use PRMS_CONSTANTS, only: MAXCONTROL_LENGTH, MAXFILE_LENGTH, MAXLINE_LENGTH
  implicit none

  integer, parameter :: MAXVARIABLES = 512

  integer, save:: Num_variables

  type PRMS_variable
    character(LEN=MAXCONTROL_LENGTH) :: variable_name
    character(LEN=MAXLINE_LENGTH) :: description
    integer :: numvals, data_flag, decl_flag, get_flag, var_name_nchars, id_num, dim_flag
    character(LEN=MAXCONTROL_LENGTH) :: data_type, dimen_names, module_name, units
    integer, pointer :: values_int_0d ! Scalars
    integer, pointer :: values_int_1d(:)
    integer, pointer :: values_int_2d(:, :)
    real, pointer :: values_real_0d
    real, pointer :: values_real_1d(:)
    real, pointer :: values_real_2d(:, :)
    double precision, pointer :: values_dble_0d
    double precision, pointer :: values_dble_1d(:)
    double precision, pointer :: values_dble_2d(:)
  end type PRMS_variable

  type(PRMS_variable), save, allocatable :: Variable_data(:)

  interface
    module subroutine dattim(String, Datetime)
      character(LEN=*), intent(IN) :: String
      integer, intent(OUT) :: Datetime(6)
    end subroutine

    module subroutine declvar(Modname, Varname, Dimenname, Numvalues, Data_type, Desc, Units)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Data_type, Desc, Units
      integer, intent(IN) :: Numvalues
    end subroutine

    double precision module function deltim()
    end function

    integer module function find_variable(Modname, Varname, Numvalues, Data_type)
      character(LEN=*), intent(IN) :: Modname, Varname, Data_type
      integer, intent(IN) :: Numvalues
    end function

    integer module function getvar(Modname, Varname, Numvalues, Data_type, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Data_type
      integer, intent(IN) :: Numvalues
      ! values could be any data type
      real, intent(OUT) :: Values(:)
    end function

    integer module function getvarnvals(Varname)
      character(LEN=*), intent(IN) :: Varname
    end function

    integer module function getvarsize(Varname)
      character(LEN=*), intent(IN) :: Varname
    end function

    integer module function getvartype(Varname)
      character(LEN=*), intent(IN) :: Varname
    end function

    integer module function getvar_id(Varname)
      character(LEN=*), intent(IN) :: Varname
    end function

    module subroutine set_data_type(Data_type, Type_flag)
      character(LEN=*), intent(IN) :: Data_type
      integer, intent(OUT) :: Type_flag
    end subroutine

    module function var_is_scalar(Varname) result(res)
      logical :: res
      character(len=*), intent(in) :: Varname
    end function
  end interface

  interface getvar_dble
    module subroutine getvar_dble_0d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      double precision, intent(OUT) :: Values
    end subroutine

    module subroutine getvar_dble_1d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      double precision, intent(OUT) :: Values(:)
    end subroutine

    ! module SUBROUTINE getvar_dble_2d(Modname, Varname, Numvalues, Values)
    !   CHARACTER(LEN=*), INTENT(IN) :: Modname, Varname
    !   INTEGER, INTENT(IN) :: Numvalues
    !   DOUBLE PRECISION, INTENT(OUT) :: Values(:, :)
    ! END SUBROUTINE
  end interface

  interface getvar_int
    module subroutine getvar_int_0d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values
    end subroutine

    module subroutine getvar_int_1d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values(:)
    end subroutine

    module subroutine getvar_int_2d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values(:, :)
    end subroutine
  end interface

  interface getvar_real
    module subroutine getvar_real_0d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values
    end subroutine

    module subroutine getvar_real_1d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values(:)
    end subroutine

    module subroutine getvar_real_2d(Modname, Varname, Numvalues, Values)
      character(LEN=*), intent(IN) :: Modname, Varname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values(:, :)
    end subroutine
  end interface

  interface declvar_dble
    module subroutine declvar_dble_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      double precision, target :: Values
    end subroutine

    module subroutine declvar_dble_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      double precision, target :: Values(:)
    end subroutine
  end interface

  interface declvar_int
    module subroutine declvar_int_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      integer, target :: Values
    end subroutine

    module subroutine declvar_int_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      integer, target :: Values(:)
    end subroutine

    module subroutine declvar_int_2d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      integer, target :: Values(:, :)
    end subroutine
  end interface

  interface declvar_real
    module subroutine declvar_real_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      real, target :: Values
    end subroutine

    module subroutine declvar_real_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      real, target :: Values(:)
    end subroutine

    module subroutine declvar_real_2d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
      character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
      integer, intent(IN) :: Numvalues
      real, target :: Values(:, :)
    end subroutine
  end interface

  interface getvalues_dbl
    module subroutine getvalues_dbl_0d(param_id, Values)
      integer, intent(IN) :: param_id
      double precision, intent(OUT) :: Values
    end subroutine

    module subroutine getvalues_dbl_1d(param_id, Values)
      integer, intent(IN) :: param_id
      double precision, intent(OUT) :: Values(:)
    end subroutine
  end interface

  interface getvalues_int
    module subroutine getvalues_int_0d(param_id, Values)
      integer, intent(IN) :: param_id
      integer, intent(OUT) :: Values
    end subroutine

    module subroutine getvalues_int_1d(param_id, Values)
      integer, intent(IN) :: param_id
      integer, intent(OUT) :: Values(:)
    end subroutine

    module subroutine getvalues_int_2d(param_id, Values)
      integer, intent(IN) :: param_id
      integer, intent(OUT) :: Values(:, :)
    end subroutine
  end interface

end module PRMS_MMFAPI
