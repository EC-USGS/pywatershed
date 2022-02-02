      PROGRAM PRMS_FORTRAN
!***********************************************************************
! PRMS main that controls time loop
!***********************************************************************
      USE PRMS_MODULE, ONLY: Number_timesteps
      IMPLICIT NONE
! Functions
      EXTERNAL :: call_modules
! Local Variables
      INTEGER :: i
!***********************************************************************
      CALL call_modules('setdims')

      CALL call_modules('decl')

      CALL call_modules('init')

      DO i = 1, Number_timesteps
        CALL call_modules('run')
      ENDDO
      CALL call_modules('clean')

      END PROGRAM PRMS_FORTRAN
