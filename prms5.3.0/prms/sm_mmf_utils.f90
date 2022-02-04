submodule(PRMS_MMFAPI) sm_mmf_utils
  use PRMS_CONSTANTS, only: MAXCONTROL_LENGTH, MAXFILE_LENGTH, MAXLINE_LENGTH, ERROR_var
  implicit none

contains

!***********************************************************************
! Set data type flag
!***********************************************************************
  module subroutine set_data_type(Data_type, Type_flag)
    use PRMS_CONSTANTS, only: ERROR_var
    use prms_utils, only: numchars
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Data_type
    integer, intent(OUT) :: Type_flag
    ! Local Variables
    integer :: string_length
    !***********************************************************************
    string_length = numchars(Data_type)
    if (string_length > 3 .and. Data_type(:4) == 'real') then
      Type_flag = 2
    elseif (string_length > 5 .and. Data_type(:6) == 'double') then
      Type_flag = 3
    elseif (string_length > 5 .and. Data_type(:6) == 'string') then
      Type_flag = 4
    elseif (string_length > 6 .and. Data_type(:7) == 'integer') then
      Type_flag = 1
    else
      print *, 'ERROR, invalid data type: ', Data_type
      print *, '       valid values are real, double, string, integer'
      ERROR stop ERROR_var
    end if
  end subroutine set_data_type

!***********************************************************************
! declvar - set up memory for variables
!***********************************************************************
  module subroutine declvar(Modname, Varname, Dimenname, Numvalues, Data_type, Desc, Units)
    use PRMS_CONSTANTS, only: ERROR_var
    use prms_utils, only: error_stop, numchars
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Data_type, Desc, Units
    integer, intent(IN) :: Numvalues
    ! Local Variables
    integer :: type_flag
    !***********************************************************************
    ! need to declare parameters first, but don't know how many, know how many in Parameter File
    Num_variables = Num_variables + 1
    if (Num_variables > MAXVARIABLES) then
      print '(A,I0)', 'PRMS ERROR, maximum number of declared variables exceeded: ', MAXVARIABLES
      call error_stop('maximum number of declared variables exceeded', ERROR_var)
    end if
    Variable_data(Num_variables)%get_flag = 0
    Variable_data(Num_variables)%decl_flag = 1
    Variable_data(Num_variables)%variable_name = Varname
    Variable_data(Num_variables)%var_name_nchars = numchars(Varname)
    Variable_data(Num_variables)%description = Desc
    Variable_data(Num_variables)%units = Units
    Variable_data(Num_variables)%dimen_names = Dimenname
    Variable_data(Num_variables)%module_name = Modname
    Variable_data(Num_variables)%numvals = Numvalues
    Variable_data(Num_variables)%data_type = Data_type
    Variable_data(Num_variables)%id_num = Num_variables
    Variable_data(Num_variables)%dim_flag = 1
    call set_data_type(Data_type, type_flag)
    if (type_flag < 1 .or. type_flag > 3) then
      print *, 'ERROR, data type not implemented: ', Data_type, ' Variable: ', &
   &           Varname(:Variable_data(Num_variables)%var_name_nchars)
      ERROR stop ERROR_var
    end if
    Variable_data(Num_variables)%data_flag = type_flag
  end subroutine declvar

!***********************************************************************
! declvar_dble_0d - set up memory for scalar double precision variables
!***********************************************************************
  module subroutine declvar_dble_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    double precision, target :: Values
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'double', Desc, Units)
    Variable_data(Num_variables)%values_dble_0d => Values
    Variable_data(Num_variables)%dim_flag = 0
  end subroutine declvar_dble_0d

!***********************************************************************
! declvar_dble_1d - set up memory for 1-dimensional double precision variables
!***********************************************************************
  module subroutine declvar_dble_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
      ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    double precision, target :: Values(:)
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'double', Desc, Units)
    Variable_data(Num_variables)%values_dble_1d => Values(:Numvalues)
    Variable_data(Num_variables)%dim_flag = 1
  end subroutine declvar_dble_1d

!***********************************************************************
! declvar_dble_2d - set up memory for 2-dimensional double precision variables
!***********************************************************************
!  module subroutine declvar_dble_2d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
!      implicit none
!      ! Arguments
!      CHARACTER(LEN=*), INTENT(IN) :: Modname, Varname, Dimenname, Desc, Units
!      INTEGER, INTENT(IN) :: Numvalues
!      DOUBLE PRECISION, TARGET :: Values(:, :)
      ! Local variables
!***********************************************************************
!      call declvar(Modname, Varname, Dimenname, Numvalues, 'double', Desc, Units)
!      Variable_data(Num_variables)%values_dble_2d => Values
!      Variable_data(Num_variables)%dim_flag = 2
!  end subroutine declvar_dble_2d

!***********************************************************************
! declvar_real_0d - set up memory for scalar real variables
!***********************************************************************
  module subroutine declvar_real_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    real, target :: Values
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'real', Desc, Units)
    Variable_data(Num_variables)%values_real_0d => Values
    Variable_data(Num_variables)%dim_flag = 0
  end subroutine declvar_real_0d

!***********************************************************************
! declvar_real_1d - set up memory for 1-dimensional real variables
!***********************************************************************
  module subroutine declvar_real_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    real, target :: Values(:)
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'real', Desc, Units)
    Variable_data(Num_variables)%values_real_1d => Values(:Numvalues)
    Variable_data(Num_variables)%dim_flag = 1
  end subroutine declvar_real_1d

!***********************************************************************
! declvar_real_2d - set up memory for 2-dimensional real variables
!***********************************************************************
  module subroutine declvar_real_2d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    real, target :: Values(:, :)
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'real', Desc, Units)
    Variable_data(Num_variables)%values_real_2d => Values
    Variable_data(Num_variables)%dim_flag = 2
  end subroutine declvar_real_2d

!***********************************************************************
! declvar_int - set up memory for scalar integer variables
!***********************************************************************
  module subroutine declvar_int_0d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    integer, target :: Values
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'integer', Desc, Units)
    Variable_data(Num_variables)%values_int_0d => Values
    Variable_data(Num_variables)%dim_flag = 0
  end subroutine declvar_int_0d

!***********************************************************************
! declvar_int - set up memory for 1-dimensional integer variables
!***********************************************************************
  module subroutine declvar_int_1d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    integer, target :: Values(:)
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'integer', Desc, Units)
    Variable_data(Num_variables)%values_int_1d => Values(:Numvalues)
    Variable_data(Num_variables)%dim_flag = 1
  end subroutine declvar_int_1d

!***********************************************************************
! declvar_int_2d - set up memory for two-dimensional integer variables
!***********************************************************************
  module subroutine declvar_int_2d(Modname, Varname, Dimenname, Numvalues, Desc, Units, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Dimenname, Desc, Units
    integer, intent(IN) :: Numvalues
    integer, target :: Values(:, :)
    !***********************************************************************
    call declvar(Modname, Varname, Dimenname, Numvalues, 'integer', Desc, Units)
    Variable_data(Num_variables)%values_int_2d => Values
    Variable_data(Num_variables)%dim_flag = 2
  end subroutine declvar_int_2d

!***********************************************************************
! getvar_dble_0d - get double precision scalar variable values
!***********************************************************************
  module subroutine getvar_dble_0d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    double precision, intent(OUT) :: Values
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'double')
    Values = Variable_data(var_id)%values_dble_0d
  end subroutine getvar_dble_0d

!***********************************************************************
! getvar_dble - get double precision 1-dimensional variable values
!***********************************************************************
  module subroutine getvar_dble_1d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    double precision, intent(OUT) :: Values(:)
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'double')
    Values = Variable_data(var_id)%values_dble_1d
  end subroutine getvar_dble_1d

!***********************************************************************
! getvar_dble_2d - get double precision 2-dimensional variable values
!***********************************************************************
! module SUBROUTINE getvar_dble_2d(Modname, Varname, Numvalues, Values)
!      IMPLICIT NONE
!      ! Arguments
!      CHARACTER(LEN=*), INTENT(IN) :: Modname, Varname
!      INTEGER, INTENT(IN) :: Numvalues
!      DOUBLE PRECISION, INTENT(OUT) :: Values(:, :)
!      ! Local Variables
!      INTEGER :: var_id
!***********************************************************************
!      var_id = find_variable(Modname, Varname, Numvalues, 'double')
!      Values = Variable_data(var_id)%values_dble_2d
!  end subroutine getvar_dble_2d

!***********************************************************************
! getvar_real_0d - get scalar single precision variable values
!***********************************************************************
  module subroutine getvar_real_0d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    real, intent(OUT) :: Values
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    !print *, modname, varname, numvalues
      var_id = find_variable(Modname, Varname, Numvalues, 'real')
      Values = Variable_data(var_id)%values_real_0d
  end subroutine getvar_real_0d

!***********************************************************************
! getvar_real - get single precision variable values
!***********************************************************************
  module subroutine getvar_real_1d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    real, intent(OUT) :: Values(:)
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'real')
    Values = Variable_data(var_id)%values_real_1d
  end subroutine getvar_real_1d

  module subroutine getvar_real_2d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    real, intent(OUT) :: Values(:, :)
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'real')
    Values = Variable_data(var_id)%values_real_2d
  end subroutine getvar_real_2d

!***********************************************************************
! getvar_int_0d - get scalar integer variable values
!***********************************************************************
  module subroutine getvar_int_0d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    integer, intent(OUT) :: Values
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'integer')
    Values = Variable_data(var_id)%values_int_0d
  end subroutine getvar_int_0d

!***********************************************************************
! getvar_int - get integer 1-dimensional variable values
!***********************************************************************
  module subroutine getvar_int_1d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    integer, intent(OUT) :: Values(:)
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'integer')
    Values = Variable_data(var_id)%values_int_1d
  end subroutine getvar_int_1d

  module subroutine getvar_int_2d(Modname, Varname, Numvalues, Values)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname
    integer, intent(IN) :: Numvalues
    integer, intent(OUT) :: Values(:, :)
    ! Local Variables
    integer :: var_id
    !***********************************************************************
    var_id = find_variable(Modname, Varname, Numvalues, 'integer')
    Values = Variable_data(var_id)%values_int_2d
  end subroutine

!***********************************************************************
! getvar - get variable values
!***********************************************************************
!  integer module function getvar(Modname, Varname, Numvalues, Data_type, Values)
!    implicit none
!    ! Arguments
!    character(LEN=*), intent(IN) :: Modname, Varname, Data_type
!    integer, intent(IN) :: Numvalues
!    ! values could be any data type
!    real, intent(OUT) :: Values(:)

!    ! Local Variables
!    integer :: var_id, var_type
!    integer, allocatable :: itemp(:)
!    real, allocatable :: temp(:)
!    double precision, allocatable :: dtemp(:)
!    !***********************************************************************
!    var_id = find_variable(Modname, Varname, Numvalues, Data_type)
!    var_type = Variable_data(var_id)%data_flag

!    if (var_type == 1) then
!      allocate (itemp(Numvalues))
!      itemp = Variable_data(var_id)%values_int_1d
!      Values = transfer(itemp, Values)
!      deallocate (itemp)
!    elseif (var_type == 2) then
!      allocate (temp(Numvalues))
!      temp = Variable_data(var_id)%values_real_1d
!      Values = transfer(temp, Values)
!      deallocate (temp)
!    elseif (var_type == 3) then
!      allocate (dtemp(Numvalues))
!      dtemp = Variable_data(var_id)%values_dble_1d
!      Values = transfer(dtemp, Values)
!      deallocate (dtemp)
!    end if

    !Values = TRANSFER(Variable_data(var_id)%values_dble,Values)
    !do i = 1, Numvalues
    !    !Values(i) = SNGL(Variable_data(var_id)%values_dble(i))
    !    Values(i) = temp(i)
    !    IF ( values(i)<0.0 ) values(i) = 0.0
    !    !print *, Variable_data(var_id)%values_dble(i)
    !ENDDO
    !values = temp

!    getvar = 0
!  end function getvar

!***********************************************************************
! find_variable - find variable in data structure
!***********************************************************************
  integer module function find_variable(Modname, Varname, Numvalues, Data_type)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Modname, Varname, Data_type
    integer, intent(IN) :: Numvalues
    ! Functions
    intrinsic :: TRIM
    ! Local Variables
    integer :: found, i, ierr
    !***********************************************************************
    ierr = 0
    found = 0
    find_variable = 1
    do i = 1, Num_variables
      if (Varname == trim(Variable_data(i)%variable_name)) then
        found = 1
        if (Variable_data(i)%numvals /= Numvalues) then
          ierr = 1
          print *, 'ERROR in: ', Modname, ', Variable: ', Varname, &
   &               ' number of values in getvar does not match declared number of values'
        end if
        if (trim(Variable_data(i)%data_type) /= Data_type) then
          ierr = 1
          print *, 'ERROR in: ', Modname, ', Variable: ', Varname, ' data type does in getvar not match declared data type'
        end if
        find_variable = i
        exit
      end if
    end do

    if (found == 0) then
      print *, 'ERROR in: ', Modname, ', Variable: ', Varname, ' not declared'
      ierr = 1
    end if
    if (ierr == 1) ERROR stop ERROR_var

  end function find_variable

!***********************************************************************
! getvar_id - get variable index
!***********************************************************************
  integer module function getvar_id(Varname)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Varname
    ! Functions
    intrinsic :: TRIM
    ! Local Variables
    integer :: i
    !***********************************************************************
    getvar_id = 1
    do i = 1, Num_variables
      if (Varname == trim(Variable_data(i)%variable_name)) then
        getvar_id = Variable_data(i)%id_num
        return
      end if
    end do
    print *, 'ERROR variable: ', Varname, ' not available'
    ERROR stop ERROR_var
  end function getvar_id

!***********************************************************************
! getvartype - get variable type
!***********************************************************************
  integer module function getvartype(Varname)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Varname
    ! Functions
    intrinsic TRIM
    ! Local Variables
    integer :: i
    !***********************************************************************
    getvartype = 1
    do i = 1, Num_variables
      if (Varname == trim(Variable_data(i)%variable_name)) then
        getvartype = Variable_data(i)%data_flag
        getvartype = getvartype
        return
      end if
    end do
    print *, 'ERROR variable: ', Varname, ' not available'
    ERROR stop ERROR_var
  end function getvartype

!***********************************************************************
! getvarnvals - get variable number of values
!***********************************************************************
  integer module function getvarnvals(Varname)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Varname
    ! Functions
    intrinsic :: TRIM
    ! Local Variables
    integer :: i
    !***********************************************************************
    getvarnvals = 1
    do i = 1, Num_variables
      if (Varname == trim(Variable_data(i)%variable_name)) then
        getvarnvals = Variable_data(i)%numvals
        return
      end if
    end do
    print *, 'ERROR in: getvarnvals, Variable: ', Varname, ' not declared'
    ERROR stop ERROR_var
  end function getvarnvals

  module function var_is_scalar(Varname) result(res)
    implicit none

    ! Arguments
    logical :: res
    character(len=*), intent(in) :: Varname

    ! Local variables
    integer :: ii

    ! ********************************************************************
    res = .false.
    do ii = 1, Num_variables
      if (Varname == trim(Variable_data(ii)%variable_name)) then
        if (trim(Variable_data(ii)%dimen_names) == 'one') then
          res = .true.
        end if

        return
      end if
    end do
  end function

!***********************************************************************
! getvalues_int - get values from parameter data structure
!***********************************************************************
  module subroutine getvalues_int_0d(param_id, Values)
    implicit none
    ! Arguments
    integer, intent(IN) :: param_id
    integer, intent(OUT) :: Values
    !***********************************************************************
    Values = Variable_data(param_id)%values_int_0d
  end subroutine

  module subroutine getvalues_int_1d(param_id, Values)
    implicit none
    ! Arguments
    integer, intent(IN) :: param_id
    integer, intent(OUT) :: Values(:)
    !***********************************************************************
    Values = Variable_data(param_id)%values_int_1d
  end subroutine

  module subroutine getvalues_int_2d(param_id, Values)
    implicit none
    ! Arguments
    integer, intent(IN) :: param_id
    integer, intent(OUT) :: Values(:, :)
    !***********************************************************************
    Values = Variable_data(param_id)%values_int_2d
  end subroutine

!***********************************************************************
! getvalues_dbl - get values from parameter data structure
!***********************************************************************
  module subroutine getvalues_dbl_0d(param_id, Values)
    implicit none
    ! Arguments
    integer, intent(IN) :: param_id
    double precision, intent(OUT) :: Values
    !***********************************************************************
    values = Variable_data(param_id)%values_dble_0d
  end subroutine

  module subroutine getvalues_dbl_1d(param_id, Values)
    implicit none
    ! Arguments
    integer, intent(IN) :: param_id
    double precision, intent(OUT) :: Values(:)
    !***********************************************************************
    values = Variable_data(param_id)%values_dble_1d
  end subroutine

!***********************************************************************
! timestep_hours - time step increment in hours
!***********************************************************************
  double precision module function deltim()
    implicit none
    !***********************************************************************
    !deltim = lisfunction() ! need to make routine to get time step increment
    deltim = 24.0D0
  end function deltim

!***********************************************************************
! dattim - get start, end, or current date and time
!***********************************************************************
  module subroutine dattim(String, Datetime)
    use PRMS_CONSTANTS, only: ERROR_time
    use PRMS_MODULE, only: Endtime, Starttime, Nowyear, Nowmonth, Nowday
    use PRMS_SET_TIME, only: Julian_day_absolute
    use prms_utils, only: compute_gregorian, error_stop
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: String
    integer, intent(OUT) :: Datetime(6)

    ! Local variable
    integer :: string_length
    !***********************************************************************
    Datetime = 0
    string_length = len(String)
    if (String(:3) == 'end') then
      Datetime = Endtime
    elseif (String(:3) == 'now') then
      call compute_gregorian(Julian_day_absolute, Nowyear, Nowmonth, Nowday)
      Datetime(1) = Nowyear
      Datetime(2) = Nowmonth
      Datetime(3) = Nowday
      ! Datetime = LIS function
    elseif (string_length > 4) then
      if (String(:5) == 'start') then
        Datetime = Starttime
      else
        call error_stop('invalid call to dattim', ERROR_time)
      end if
    end if
  end subroutine dattim

!***********************************************************************
! getvarsize - return the number of values for a parameter
!***********************************************************************
  integer module function getvarsize(Varname)
    implicit none
    ! Arguments
    character(LEN=*), intent(IN) :: Varname
    ! Functions
    intrinsic :: TRIM
    ! Local Variables
    integer :: found, i
    !***********************************************************************
    found = 0
    do i = 1, Num_variables
      if (Varname == trim(Variable_data(i)%variable_name)) then
        found = i
        getvarsize = Variable_data(i)%numvals
        exit
      end if
    end do

    if (found == 0) then
      print *, 'ERROR, Variable: ', Varname, ' not declared'
      ERROR stop ERROR_var
    end if
  end function

end submodule
