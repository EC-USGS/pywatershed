!***********************************************************************
!     Output statvar file
!***********************************************************************
      MODULE PRMS_STATVAR_OUT
        USE PRMS_CONSTANTS, ONLY: MAXCONTROL_LENGTH
        IMPLICIT NONE
! Module Variables
        character(len=*), parameter :: MODDESC = 'Statistics Variables Output'
        character(len=*), parameter :: MODNAME = 'statvar_out'
        character(len=*), parameter :: Version_statvar_out = '2021-09-08'
        character(len=*), parameter :: stamp = '(7(I0,1X),'
        character(len=24) :: Output_fmt
        integer, save :: Statvar_unit
        integer, save, allocatable :: Nc_vars(:), Stat_var_type(:), Statvar_id(:), Statvar_nvals(:)
        real, save, allocatable :: Stat_var_values(:)
! Control Parameters
        integer, save :: nstatVars, statvarOut_format
        !integer, save :: Num_statvar_elements, Num_statvar_names
        character(len=MAXCONTROL_LENGTH), allocatable, save :: statVar_element(:), statVar_names(:)
      END MODULE PRMS_STATVAR_OUT

!     ******************************************************************
!     Statistics Variables module
!     ******************************************************************
      SUBROUTINE statvar_out()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN
      USE PRMS_MODULE, ONLY: Process_flag
      USE PRMS_STATVAR_OUT, ONLY: Statvar_unit
      IMPLICIT NONE
! Functions
      EXTERNAL :: statvar_outinit, statvar_outrun
!***********************************************************************
      IF ( Process_flag==RUN ) THEN
        CALL statvar_outrun()
      ELSEIF ( Process_flag==DECL ) THEN
        CALL statvar_outdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        CALL statvar_outinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        CLOSE ( Statvar_unit )
      ENDIF

      END SUBROUTINE statvar_out

!***********************************************************************
!     declare parameters and variables
!***********************************************************************
      SUBROUTINE statvar_outdecl()
      USE PRMS_CONSTANTS, ONLY: ERROR_control
      use PRMS_CONTROL_FILE, only: control_integer, control_string_array
      USE PRMS_STATVAR_OUT
      use prms_utils, only: error_stop, print_module, read_error
      IMPLICIT NONE
! Local Variables
      INTEGER :: i
!***********************************************************************
      CALL print_module(MODDESC, MODNAME, Version_statvar_out)

      IF ( control_integer(nstatVars, 'nstatVars')/=0 ) nstatVars = 0
      IF ( nstatVars==0 ) CALL error_stop('statvar output requested with nstatVars equal 0', ERROR_control)

      ALLOCATE ( statVar_element(nstatVars), statVar_names(nstatVars) )
      statVar_element = ' '
      statVar_names = ' '
      DO i = 1, nstatVars
        IF ( control_string_array(statVar_element(i), 'statVar_element', i)/=0 ) CALL read_error(5, 'statVar_element')
        IF ( control_string_array(statVar_names(i), 'statVar_names', i)/=0 ) CALL read_error(5, 'statVar_names')
      ENDDO
      !print *, statvar_names

      ! 1 = ES10.3; 2 = F0.2; 3 = F0.3; 4 = F0.4; 5 = F0.5
      IF ( control_integer(statvarOut_format, 'statvarOut_format')/=0 ) statvarOut_format = 1
      IF ( statvarOut_format<1 .OR. statvarOut_format>5 ) CALL error_stop('invalid statvarOut_format value', ERROR_control)

      IF ( statvarOut_format==1 ) THEN
        WRITE ( Output_fmt, 9001 ) stamp, nstatVars
      ELSEIF ( statvarOut_format==2 ) THEN
        WRITE ( Output_fmt, 9002 ) stamp, nstatVars
      ELSEIF ( statvarOut_format==3 ) THEN
        WRITE ( Output_fmt, 9003 ) stamp, nstatVars
      ELSEIF ( statvarOut_format==4 ) THEN
        WRITE ( Output_fmt, 9004 ) stamp, nstatVars
      ELSEIF ( statvarOut_format==5 ) THEN
        WRITE ( Output_fmt, 9005 ) stamp, nstatVars
      ENDIF
      !print *, output_fmt
 9001 FORMAT (A,I0,'(1X,ES10.3))')
 9002 FORMAT (A,I0,'(1X,F0.2))')
 9003 FORMAT (A,I0,'(1X,F0.3))')
 9004 FORMAT (A,I0,'(1X,F0.4))')
 9005 FORMAT (A,I0,'(1X,F0.5))')

    END SUBROUTINE statvar_outdecl

!***********************************************************************
!     statvar_outinit - Initialize statvar_out module
!***********************************************************************
      SUBROUTINE statvar_outinit()
      USE PRMS_CONSTANTS, ONLY: ERROR_control, ERROR_open_out
      use PRMS_MMFAPI, only: getvartype, getvarnvals
      USE PRMS_MODULE, ONLY: Stat_var_file
      USE PRMS_STATVAR_OUT
      use prms_utils, only: error_stop, numchars, PRMS_open_output_file
      IMPLICIT NONE
! Local Variables
      INTEGER :: jj, ios, ierr, itype
!***********************************************************************
      ierr = 0
      !IF ( Num_statvar_names/=nstatVars ) THEN
      !  PRINT *, 'ERROR, number of names specified for statVar_names not equal to nstatVars'
      !  PRINT *, '       nstatVars =', nstatVars, ' number of names = ', Num_statvar_names
      !ENDIF
      !IF ( Num_statvar_elements/=nstatVars ) THEN
      !  PRINT *, 'ERROR, number of elements specified for statVar_elements not equal to nstatVars'
      !  PRINT *, '       nstatVars =', nstatVars, ' number of names = ', Num_statvar_elements
      !ENDIF

      ALLOCATE ( Stat_var_type(nstatVars), Nc_vars(nstatVars), Statvar_id(nstatVars) )
      ALLOCATE ( Statvar_nvals(nstatVars), Stat_var_values(nstatVars) )

      CALL PRMS_open_output_file(Statvar_unit, Stat_var_file, 'statvar.out', 0, ios)
      IF ( ios/=0 ) CALL error_stop('opening statvar file', ERROR_open_out)
      WRITE ( Statvar_unit, '(I0)' ) nstatVars

      DO jj = 1, nstatVars
        Nc_vars(jj) = numchars(statVar_names(jj))
        itype = getvartype(statVar_names(jj)(:Nc_vars(jj)))
        IF ( itype<1 .AND. itype>3 ) THEN
          ierr = 1
          PRINT *, 'ERROR, invalid statvar type: 1, 2, and 3 allowed', itype
          CYCLE
        ENDIF
        Stat_var_type(jj) = itype
        Statvar_nvals(jj) = getvarnvals(statVar_names(jj)(:Nc_vars(jj)))
        READ ( statVar_element(jj), * ) Statvar_id(jj)
        WRITE ( Statvar_unit, '(A,1X,I0)' ) statVar_names(jj)(:Nc_vars(jj)), Statvar_id(jj)
      ENDDO
      IF ( ierr==1 ) CALL error_stop('statvar variable issue', ERROR_control)

      Stat_var_values = 0.0

      END SUBROUTINE statvar_outinit

!***********************************************************************
!     statvar_outrun - Output to a file the statvar variables
!***********************************************************************
      SUBROUTINE statvar_outrun()
      USE PRMS_CONSTANTS, ONLY: ERROR_write
      use PRMS_MMFAPI, only: getvarnvals, getvar_dble, getvar_int, getvar_real, var_is_scalar
      USE PRMS_STATVAR_OUT
      USE PRMS_MODULE, ONLY: Timestep
      USE PRMS_SET_TIME, ONLY: Nowtime
      use prms_utils, only: error_stop, read_error
      IMPLICIT NONE
! FUNCTIONS AND SUBROUTINES
      INTRINSIC SNGL, FLOAT, TRIM
! Local Variables
      INTEGER :: jj, nvals, nc, nvalues
      INTEGER, ALLOCATABLE, TARGET :: values_int(:)
      REAL, ALLOCATABLE, TARGET :: values_real(:)
      DOUBLE PRECISION, ALLOCATABLE, TARGET :: values_dble(:)
!***********************************************************************
      DO jj = 1, nstatVars
        nc = Nc_vars(jj)
        nvals = Statvar_nvals(jj)
        nvalues = getvarnvals(statVar_names(jj)(:nc))
        IF ( Statvar_id(jj)>nvalues ) THEN
          PRINT *, 'ERROR, element id number exceeds dimension for variable: ', statVar_names(jj)(:nc)
          STOP
        ENDIF

        IF ( Stat_var_type(jj)==3 ) THEN
          ALLOCATE ( values_dble(nvals) )
          IF ( var_is_scalar(TRIM(statVar_names(jj))) ) THEN
            CALL getvar_dble(MODNAME, statVar_names(jj)(:nc), nvals, values_dble(1))
          ELSE
            CALL getvar_dble(MODNAME, statVar_names(jj)(:nc), nvals, values_dble)
          ENDIF

          Stat_var_values(jj) = SNGL(values_dble(Statvar_id(jj)))
          !print *, statVar_names(jj)(:nc), nvals, Stat_var_values(jj), 'double'
          DEALLOCATE ( values_dble )
        ELSEIF ( Stat_var_type(jj)==2 ) THEN
          ALLOCATE ( values_real(nvals) )
          IF ( var_is_scalar(TRIM(statVar_names(jj))) ) THEN
            CALL getvar_real(MODNAME, statVar_names(jj)(:nc), nvals, values_real(1))
          ELSE
            CALL getvar_real(MODNAME, statVar_names(jj)(:nc), nvals, values_real)
          ENDIF

          Stat_var_values(jj) = values_real(Statvar_id(jj))
          !print *, statVar_names(jj)(:nc), nvals,Stat_var_values(jj), 'real'
          DEALLOCATE ( values_real )
        ELSEIF ( Stat_var_type(jj)==1 ) THEN
          ALLOCATE ( values_int(nvals) )
          IF ( var_is_scalar(TRIM(statVar_names(jj))) ) THEN
            CALL getvar_int(MODNAME, statVar_names(jj)(:nc), nvals, values_int(1))
          ELSE
            CALL getvar_int(MODNAME, statVar_names(jj)(:nc), nvals, values_int)
          ENDIF

          Stat_var_values(jj) = FLOAT(values_int(Statvar_id(jj)))
          !print *, statVar_names(jj)(:nc), nvals,Stat_var_values(jj), 'integer'
          DEALLOCATE ( values_int )
        ELSE
          CALL error_stop('statvar_out write issue', ERROR_write)
        ENDIF
      ENDDO

      WRITE ( Statvar_unit, Output_fmt ) Timestep, Nowtime, (Stat_var_values(jj), jj=1, nstatVars)

      END SUBROUTINE statvar_outrun
