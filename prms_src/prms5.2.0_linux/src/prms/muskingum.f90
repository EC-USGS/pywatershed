!***********************************************************************
! Routes water between segments in the system using Muskingum routing
!
!   The Muskingum equation is described in 'Hydrology for Engineers', 3rd ed.
!   by Linsley, R.K, Kohler, M.A., and Paulhus, J.L.H., 1982 p. 275 and in
!   'Water in Environmental Planning' by Dunne, T., and Leopold, L.B. 1978
!   p. 357.
!
!   Note that the Muskingum equation assumes a linear relation of storage
!   to the inflow/outflow relation and therefore the relation is the same
!   throughout the range of the hydrograph.  The route_time parameter in
!   the fixroute module is replaced by two new parameters, K_coef and
!   x_coef, which are described below:
!
!   The Muskingum method is based on the equation: S = K[xI + (1 - x)O]
!     where S is storage, K is the storage coefficient, x is a coefficient
!     between 0 and .5, I is inflow, and O is outflow.
!
!   Solving for the outflow at day 2,O2; and knowing the inflow at day 1,
!   I1; the inflow at day 2,I2; and the outflow at day 1, O1; the storage
!   equation can be written as follows:
!
!        O2 = czero*I2 + cone*I1 + ctwo*O1
!
!     where czero = -((Kx - 12)    / (K - Kx + 12))
!           cone  =  (Kx + 12)     / (K - Kx + 12)
!           ctwo  =  (K - Kx - 12) / (K - Kx + 12)
!
!     assuming a time step of one day and K is in units of hours
!
!   This module is based on the "musroute.f" module. It differs in three
!   basic ways:
!
!   1. This module uses an internal routing time step of one hour.
!      The old muskingum module ran on the same daily time step as
!      the rest of PRMS. The problem with this is that there is no
!      ability to distinguish where the flood wave (front of the flow
!      change) within the segment. For example, if there is a series
!      of 4 1-day long segments, a flood wave will make it to the bottom
!      of these in 1 day. If the same system is modeled as 1 4-day long
!      segment, it will take 4 days.
!
!   2. The X parameter has been removed as a specified input and is now computed. To
!      my knowledge, no modeler had ever set this to anything other than the default
!      value (0.2) anyway. Always using the default value can lead to problems
!      with the C coffecients which can result in mass balance problems or negative
!      flow values.
!
!      To solve this problem, I assume that the C coefficients must
!      always be between 0 and 1. By setting the C coefficients equal to 0 and 1,
!      various limits on the time step (ts), X, and K can be determined. There are
!      two of these limits which are of interest:
!
!      When C0 = 0:
!             ts
!        K = -----
!             2X
!
!      When C2 = 0:
!            ts
!       K = -----
!           2(1-X)
!
!      Determining a value of K half way between these two limits (by averaging)
!      and solving for X using the quadratic formula results in:
!
!            1-sqrt(1-(ts/K))
!       X = ------------------
!                  2
!
!       So when ts is fixed at one hour and K is fixed as the average (or expected)
!       travel time corresponding to the segment (for each segment in the stream
!       network), a value of X can be computed (for each segment in the stream
!       network) which will result in both conservation of mass and non-negative
!       flows. Another benefit is that only one input parameter (K) needs to be
!       input to the module.
!
!   3. If the travel time of a segment is less than or equal to the routing
!      time step (one hour), then the outflow of the segment is set to the
!      value of the inflow.
!
!***********************************************************************
      MODULE PRMS_MUSKINGUM
      USE PRMS_CONSTANTS, ONLY: ACTIVE, NEARZERO, CFS2CMS_CONV, &
     &    OUTFLOW_SEGMENT, ERROR_streamflow, strmflow_muskingum_module
      USE PRMS_MODULE, ONLY: Nsegment, Init_vars_from_file, Strmflow_flag, Glacier_flag
      IMPLICIT NONE
      character(len=*), parameter :: MODDESC = 'Streamflow Routing'
      character(len=14), parameter :: MODNAME = 'muskingum_mann'
      character(len=*), parameter :: Version_muskingum = '2020-12-02'
!   Local Variables
      DOUBLE PRECISION, PARAMETER :: ONE_24TH = 1.0D0 / 24.0D0
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Currinsum(:), Pastin(:), Pastout(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Outflow_ts(:), Inflow_ts(:)
      END MODULE PRMS_MUSKINGUM

!***********************************************************************
!     Main muskingum routine
!***********************************************************************
      INTEGER FUNCTION muskingum()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: muskingum_decl, muskingum_init, muskingum_run
      EXTERNAL :: muskingum_restart
!***********************************************************************
      muskingum = 0

      IF ( Process_flag==RUN ) THEN
        muskingum = muskingum_run()
      ELSEIF ( Process_flag==DECL ) THEN
        muskingum = muskingum_decl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL muskingum_restart(READ_INIT)
        muskingum = muskingum_init()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL muskingum_restart(SAVE_INIT)
      ENDIF

      END FUNCTION muskingum

!***********************************************************************
!     muskingum_decl - Declare parameters and variables and allocate arrays
!   Declared Parameters
!     tosegment, hru_segment, obsin_segment, K_coef, x_coef
!***********************************************************************
      INTEGER FUNCTION muskingum_decl()
      USE PRMS_MUSKINGUM
      IMPLICIT NONE
! Functions
      EXTERNAL :: print_module
!***********************************************************************
      muskingum_decl = 0

      IF ( Strmflow_flag==strmflow_muskingum_module ) THEN
        CALL print_module(MODDESC, MODNAME(:9), Version_muskingum)
      ELSE ! muskingum_mann
        CALL print_module(MODDESC, MODNAME, Version_muskingum)
      ENDIF

      ALLOCATE ( Currinsum(Nsegment) )
      ALLOCATE ( Pastin(Nsegment), Pastout(Nsegment) )
      ALLOCATE ( Outflow_ts(Nsegment), Inflow_ts(Nsegment) )

      END FUNCTION muskingum_decl

!***********************************************************************
!    muskingum_init - Get and check parameter values and initialize variables
!***********************************************************************
      INTEGER FUNCTION muskingum_init()
      USE PRMS_MUSKINGUM
      USE PRMS_BASIN, ONLY: Basin_area_inv
      USE PRMS_FLOWVARS, ONLY: Seg_outflow
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_ROUTING, ONLY: Basin_segment_storage
      IMPLICIT NONE
! Local Variables
      INTEGER :: i
!***********************************************************************
      muskingum_init = 0

      !Seg_outflow will have been initialized to Segment_flow_init in PRMS_ROUTING
      IF ( Init_vars_from_file==0 ) Outflow_ts = 0.0D0

      Basin_segment_storage = 0.0D0
      DO i = 1, Nsegment
        Basin_segment_storage = Basin_segment_storage + Seg_outflow(i)
      ENDDO
      Basin_segment_storage = Basin_segment_storage*Basin_area_inv/Cfs_conv

      END FUNCTION muskingum_init

!***********************************************************************
!     muskingum_run - Compute routing summary values
!***********************************************************************
      INTEGER FUNCTION muskingum_run()
      USE PRMS_MUSKINGUM
      USE PRMS_BASIN, ONLY: Basin_area_inv, Basin_gl_cfs, Basin_gl_ice_cfs
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_cms, Basin_gwflow_cfs, Basin_ssflow_cfs, &
     &    Basin_stflow_out, Basin_cfs, Basin_stflow_in, Basin_sroff_cfs, Seg_inflow, Seg_outflow, &
     &    Seg_upstream_inflow, Seg_lateral_inflow, Flow_out, Basin_sroff
      USE PRMS_OBS, ONLY: Streamflow_cfs
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_ROUTING, ONLY: Use_transfer_segment, Segment_delta_flow, Basin_segment_storage, &
     &    Obsin_segment, Segment_order, Tosegment, C0, C1, C2, Ts, Ts_i, Obsout_segment, &
     &    Flow_to_ocean, Flow_to_great_lakes, Flow_out_region, Flow_out_NHM, Segment_type, Flow_terminus, &
     &    Flow_to_lakes, Flow_replacement, Flow_in_region, Flow_in_nation, Flow_headwater, Flow_in_great_lakes
      USE PRMS_GLACR, ONLY: Basin_gl_top_melt, Basin_gl_ice_melt
      USE PRMS_GWFLOW, ONLY: Basin_gwflow
      IMPLICIT NONE
! Functions
      INTRINSIC :: MOD
      EXTERNAL :: error_stop
! Local Variables
      INTEGER :: i, j, iorder, toseg, imod, tspd, segtype
      DOUBLE PRECISION :: area_fac, segout, currin
!***********************************************************************
      muskingum_run = 0

!     SET yesterdays inflows and outflows into temp (past arrays)
!     values may be 0.0 as intial, > 0.0 for runtime and dynamic
!     initial condtions. Then set outlfow and inflow for this time
!     step to 0.0
!
!     upstream_inflow and outflow will vary by hour
!     lateral_inflow and everything else will vary by day
!
!     Compute surface runoff, ssflow, and gwflow going to each segment
!     This is todays "seg_inflow" before additional water is routed to
!     a new (if any is routed)
!
!     For each HRU if the lateral flow for this HRU goes to the
!     segment being evaluated (segment i) then sum flows
!
!     Do these calculations once for the current day, before the hourly
!     routing starts.
!
!       Out2   =      In2*C0    +        In1*C1    +          Out1*C2
!     Seg_outflow = Seg_inflow*Czero + Pastinflow*Cone + Pastoutflow*Ctwo
!       C0, C1, and C2: initialized in the "init" part of this module
!
      Pastin = Seg_inflow
      Pastout = Seg_outflow
      Seg_inflow = 0.0D0
      Seg_outflow = 0.0D0
      Inflow_ts = 0.0D0
      Currinsum = 0.0D0

! 24 hourly timesteps per day
      DO j = 1, 24

        Seg_upstream_inflow = 0.0D0
        DO i = 1, Nsegment
          iorder = Segment_order(i)

! current inflow to the segment is the time weighted average of the outflow
! of the upstream segments plus the lateral HRU inflow plus any gains.
          currin = Seg_lateral_inflow(iorder) !note, this routes to inlet
          IF ( Obsin_segment(iorder)>0 ) Seg_upstream_inflow(iorder) = Streamflow_cfs(Obsin_segment(iorder))
          currin = currin + Seg_upstream_inflow(iorder)
          Seg_inflow(iorder) = Seg_inflow(iorder) + currin
          Inflow_ts(iorder) = Inflow_ts(iorder) + currin
          Currinsum(iorder) = Currinsum(iorder) + Seg_upstream_inflow(iorder)

          ! Check to see if this segment is to be routed on this time step
          tspd = Ts_i(iorder)
          imod = MOD( j, tspd )
          IF ( imod==0 ) THEN
            Inflow_ts(iorder) = (Inflow_ts(iorder) / Ts(iorder))
! Compute routed streamflow
            IF ( Ts_i(iorder)>0 ) THEN
! Muskingum routing equation
              Outflow_ts(iorder) = Inflow_ts(iorder)*C0(iorder) + Pastin(iorder)*C1(iorder) + Outflow_ts(iorder)*C2(iorder)
            ELSE
! If travel time (K_coef paremter) is less than or equal to
! time step (one hour), then the outflow is equal to the inflow
! Outflow_ts is the value from last hour
              Outflow_ts(iorder) = Inflow_ts(iorder)
            ENDIF

            ! pastin is equal to the Inflow_ts on the previous routed timestep
            Pastin(iorder) = Inflow_ts(iorder)

! because the upstream inflow from streams is used, reset it to zero so new average
! can be computed next routing timestep.
            Inflow_ts(iorder) = 0.0D0
          ENDIF

          IF ( Obsout_segment(iorder)>0 ) Outflow_ts(iorder) = Streamflow_cfs(Obsout_segment(iorder))

          ! water-use removed/added in routing module
          ! check for negative flow
          IF ( Outflow_ts(iorder)<0.0 ) THEN
            IF ( Use_transfer_segment==1 ) THEN
              PRINT *, 'ERROR, transfer(s) from stream segment:', iorder, ' causes outflow to be negative'
              PRINT *, '       outflow =', Outflow_ts(iorder), ' must fix water-use stream segment transfer file'
            ELSE
              PRINT *, 'ERROR, outflow from segment:', iorder, ' is negative:', Outflow_ts(iorder)
              PRINT *, '       routing parameters may be invalid'
            ENDIF
            CALL error_stop('negative streamflow in muskingum', ERROR_streamflow)
          ENDIF

          ! Seg_outflow (the mean daily flow rate for each segment) will be the average of the hourly values.
          Seg_outflow(iorder) = Seg_outflow(iorder) + Outflow_ts(iorder)
          ! pastout is equal to the Inflow_ts on the previous routed timestep
          Pastout(iorder) = Outflow_ts(iorder)

! Add current timestep's flow rate to sum the upstream flow rates.
! This can be thought of as a volume because it is a volumetric rate
! (cubic feet per second) over a time step of an hour. Down below when
! this value is used, it will be divided by the number of hours in the
! segment's simulation time step, giving the mean flow rate over that
! period of time.
          toseg = Tosegment(iorder)
          IF ( toseg>0 ) Seg_upstream_inflow(toseg) = Seg_upstream_inflow(toseg) + Outflow_ts(iorder)

        ENDDO ! segment

      ENDDO  ! timestep

      Basin_segment_storage = 0.0D0
      Flow_out = 0.0D0
      Flow_to_lakes = 0.0D0
      Flow_to_ocean = 0.0D0
      Flow_to_great_lakes = 0.0D0
      Flow_out_region = 0.0D0
      Flow_out_NHM = 0.0D0
      Flow_in_region = 0.0D0
      Flow_terminus = 0.0D0
      Flow_in_nation = 0.0D0
      Flow_headwater = 0.0D0
      Flow_in_great_lakes = 0.0D0
      Flow_replacement = 0.0D0
      DO i = 1, Nsegment
        Seg_outflow(i) = Seg_outflow(i) * ONE_24TH
        segout = Seg_outflow(i)
        segtype = Segment_type(i)
        Seg_inflow(i) = Seg_inflow(i) * ONE_24TH
        Seg_upstream_inflow(i) = Currinsum(i) * ONE_24TH
! Flow_out is the total flow out of the basin, which allows for multiple outlets
! includes closed basins (tosegment=0)
        IF ( segtype==1 ) THEN
          Flow_headwater = Flow_headwater + segout
        ELSEIF ( segtype==2 ) THEN
          Flow_to_lakes = Flow_to_lakes + segout
        ELSEIF ( segtype==3 ) THEN
          Flow_replacement = Flow_replacement + segout
        ELSEIF ( segtype==4 ) THEN
          Flow_in_nation = Flow_in_nation + segout
        ELSEIF ( segtype==5 ) THEN
          Flow_out_NHM = Flow_out_NHM + segout
        ELSEIF ( segtype==6 ) THEN
          Flow_in_region = Flow_in_region + segout
        ELSEIF ( segtype==7 ) THEN
          Flow_out_region = Flow_out_region + segout
        ELSEIF ( segtype==8 ) THEN
          Flow_to_ocean = Flow_to_ocean + segout
        ELSEIF ( segtype==9 ) THEN
          Flow_terminus = Flow_terminus + segout
        ELSEIF ( segtype==10 ) THEN
          Flow_in_great_lakes = Flow_in_great_lakes + segout
        ELSEIF ( segtype==11 ) THEN
          Flow_to_great_lakes = Flow_to_great_lakes + segout
        ENDIF
        IF ( Tosegment(i)==OUTFLOW_SEGMENT ) Flow_out = Flow_out + segout
        Segment_delta_flow(i) = Segment_delta_flow(i) + Seg_inflow(i) - segout
!        IF ( Segment_delta_flow(i) < 0.0D0 ) PRINT *, 'negative delta flow', Segment_delta_flow(i)
        Basin_segment_storage = Basin_segment_storage + Segment_delta_flow(i)
      ENDDO

      area_fac = Cfs_conv/Basin_area_inv
      Basin_stflow_in = Basin_sroff + Basin_gwflow + Basin_ssflow ! not equal to basin_stflow_out if replacement flows
      Basin_cfs = Flow_out
      Basin_stflow_out = Basin_cfs / area_fac
      Basin_cms = Basin_cfs*CFS2CMS_CONV
      IF ( Glacier_flag==ACTIVE ) THEN
        Basin_stflow_in = Basin_stflow_in + Basin_gl_top_melt
        Basin_gl_ice_cfs = Basin_gl_ice_melt*area_fac
        Basin_gl_cfs = Basin_gl_top_melt*area_fac
      ENDIF
      Basin_sroff_cfs = Basin_sroff*area_fac
      Basin_ssflow_cfs = Basin_ssflow*area_fac
      Basin_gwflow_cfs = Basin_gwflow*area_fac
      Basin_segment_storage = Basin_segment_storage/area_fac

      END FUNCTION muskingum_run

!***********************************************************************
!     muskingum_restart - write or read restart file
!***********************************************************************
      SUBROUTINE muskingum_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_MUSKINGUM
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Function
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=14) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Outflow_ts
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Outflow_ts
      ENDIF
      END SUBROUTINE muskingum_restart
