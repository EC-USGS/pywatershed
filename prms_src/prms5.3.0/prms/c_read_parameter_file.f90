module PRMS_READ_PARAM_FILE
  use PRMS_CONSTANTS, only: MAXCONTROL_LENGTH, MAXFILE_LENGTH, MAXLINE_LENGTH
  implicit none
  ! DANGER, DANGER, hard coded maximum number of paraemters and dimensions, DANGER, DANGER
  integer, parameter :: MAXDIMENSIONS = 48, MAXPARAMETERS = 256
  integer, save :: Num_parameters, Num_dimensions

  type PRMS_parameter
    character(LEN=MAXCONTROL_LENGTH) :: param_name
    character(LEN=512) :: short_description, long_description
    integer :: numvals, data_flag, decl_flag, read_flag, nchars, id_num, scalar_flag
    integer :: default_int, maximum_int, minimum_int, num_dimens, num_dim1, num_dim2
    character(LEN=MAXCONTROL_LENGTH) :: max_value, min_value, def_value, data_type
    character(LEN=MAXCONTROL_LENGTH) :: dimen_names, module_name, units
    real, pointer :: values_real_0d ! Scalars
    real, pointer :: values_real_1d(:)
    real, pointer :: values_real_2d(:, :)
    integer, pointer :: values_int_0d  ! Scalars
    integer, pointer :: values_int_1d(:)
    integer, pointer :: values_int_2d(:, :)
    real :: maximum, minimum, default_real
  end type PRMS_parameter

  type(PRMS_parameter), save, allocatable :: Parameter_data(:)

  type PRMS_dimension
    character(LEN=16) :: name
    integer :: value, default, maximum, length
    character(LEN=MAXFILE_LENGTH) :: description
  end type PRMS_dimension

  type(PRMS_dimension), save, allocatable :: Dimension_data(:)

  ! Local Variables
  character(len=*), parameter :: MODDESC = 'Read values from Parameter File'

  interface
    module subroutine check_parameters_declared(Parmname, Modname, Iret)
      character(LEN=*), intent(IN) :: Parmname, Modname
      integer, intent(OUT) :: Iret
    end subroutine

    integer module function decldim(Dimname, Defval, Maxval, Desc)
      integer, intent(IN) :: Defval, Maxval
      character(LEN=*), intent(IN) :: Dimname, Desc
    end function

    integer module function declfix(Dimname, Defval, Maxval, Desc)
      integer, intent(IN) :: Defval, Maxval
      character(LEN=*), intent(IN) :: Dimname, Desc
    end function

    integer module function declparam(Modname, Paramname, Dimenname, Datatype, &
                                      Defvalue, Minvalue, Maxvalue, Descshort, Desclong, Units)
      character(LEN=*), intent(IN) :: Modname, Paramname, Dimenname, Datatype
      character(LEN=*), intent(IN) :: Defvalue, Minvalue, Maxvalue, Descshort, Desclong, Units
    end function

    integer module function getdim(Dimname)
      character(LEN=*), intent(IN) :: Dimname
    end function

    integer module function getparamstring(Paramname, Numvalues, Data_type, Array_index, String)
      character(LEN=*), intent(IN) :: Paramname, Data_type
      integer, intent(IN) :: Numvalues, Array_index
      character(LEN=*), intent(OUT) :: String
    end function

    module subroutine setdimension(Dimname, Dim)
      integer, intent(IN) :: Dim
      character(LEN=*), intent(IN) :: Dimname
    end subroutine

    module subroutine setparam(Paramname, Numvalues, Data_type, Num_dims, Dim_string, Values, Ivalues)
      integer, intent(IN) :: Numvalues, Data_type, Num_dims, Ivalues(:)
      character(LEN=*), intent(IN) :: Paramname, Dim_string(Num_dims)
      real, intent(IN) :: Values(:)
    end subroutine

    module subroutine setup_params()
    end subroutine

    module subroutine setup_dimens()
    end subroutine

    module subroutine read_parameter_file_dimens()
    end subroutine

    module subroutine read_parameter_file_params()
    end subroutine

    module subroutine check_parameters()
    end subroutine
  end interface

  interface getparam_int
    module function getparam_int_0d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values
    end function

    module function getparam_int_1d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values(:)
    end function

    module function getparam_int_2d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      integer, intent(OUT) :: Values(:, :)
    end function
  end interface

  interface getparam_real
    module function getparam_real_0d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res  ! Function result
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values
    end function

    module function getparam_real_1d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res  ! Function result
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values(:)
    end function

    module function getparam_real_2d(Modname, Paramname, Numvalues, Values) result(res)
      integer :: res  ! Function result
      character(LEN=*), intent(IN) :: Modname, Paramname
      integer, intent(IN) :: Numvalues
      real, intent(OUT) :: Values(:, :)
    end function
  end interface
end module PRMS_READ_PARAM_FILE
