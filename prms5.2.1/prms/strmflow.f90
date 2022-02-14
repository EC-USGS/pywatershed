!***********************************************************************
! Computes daily streamflow as the sum of surface runoff,
! shallow-subsurface flow (interflow), and ground-water flow
!***********************************************************************
      INTEGER FUNCTION strmflow()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, CFS2CMS_CONV
      USE PRMS_MODULE, ONLY: Process_flag
      USE PRMS_BASIN, ONLY: Active_area
      USE PRMS_GWFLOW, ONLY: Basin_gwflow
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_cfs, Basin_cms, Basin_stflow_in, &
     &    Basin_sroff_cfs, Basin_ssflow_cfs, Basin_gwflow_cfs, Basin_stflow_out, Basin_sroff
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      IMPLICIT NONE
! Functions
      EXTERNAL :: print_module
! Local Variables
      character(len=*), parameter :: MODDESC = 'Streamflow Routing'
      character(len=*), parameter :: MODNAME = 'strmflow'
      character(len=*), parameter :: Version_strmflow = '2020-08-03'
      double precision :: area_fac
!***********************************************************************
      strmflow = 0

      IF ( Process_flag==RUN ) THEN
!   Compute daily flow.
        area_fac = Cfs_conv*Active_area
        Basin_stflow_in = Basin_sroff + Basin_gwflow + Basin_ssflow
        Basin_stflow_out = Basin_stflow_in
        Basin_cfs = Basin_stflow_in*area_fac
        Basin_cms = Basin_cfs*CFS2CMS_CONV
        Basin_sroff_cfs = Basin_sroff*area_fac
        Basin_ssflow_cfs = Basin_ssflow*area_fac
        Basin_gwflow_cfs = Basin_gwflow*area_fac

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_strmflow)
      ENDIF

      END FUNCTION strmflow
