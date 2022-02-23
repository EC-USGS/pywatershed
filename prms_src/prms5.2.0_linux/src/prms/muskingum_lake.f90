! need to add lake shape and changing area in same manner as depression storage

!***********************************************************************
! Routes water between segments and lakes in the stream network;
! using Muskingum routing for stream segments and 5 lake routing methods
!
! gwflow goes to GWR instead of to the lake unless specified as
! going to stream segment associated with the lake, which would be a
! problem
!
! nlake_hrus set to nlake for version 5.0.0, nlake_hrus to be added in 5.0.1
! in future this module may be used for muskingum only, so would need to
! check lake_route_flag = 1 (ACTIVE) in a bunch of places
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
      MODULE PRMS_MUSKINGUM_LAKE
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, NEARZERO, CFS2CMS_CONV, &
     &    OUTFLOW_SEGMENT, DOCUMENTATION, NEARZERO, DNEARZERO, CFS2CMS_CONV, LAKE, ERROR_water_use, &
     &    ERROR_streamflow, ERROR_dim, CASCADE_OFF, CASCADE_HRU_SEGMENT
      USE PRMS_MODULE, ONLY: Model, Nsegment, Nhru, Nratetbl, Nlake, &
     &    Print_debug, Init_vars_from_file, Strmflow_flag, Cascade_flag, Inputerror_flag, &
     &    Glacier_flag, Lake_transfer_water_use, Lake_add_water_use
      IMPLICIT NONE
!   Local Variables
      DOUBLE PRECISION, PARAMETER :: ONE_24TH = 1.0D0 / 24.0D0
      character(len=*), parameter :: MODDESC = 'Streamflow & Lake Routing'
      character(len=14), parameter :: MODNAME = 'muskingum_lake'
      character(len=*), parameter :: Version_muskingum_lake = '2020-12-02'
      INTEGER, SAVE :: Obs_flag, Linear_flag, Weir_flag, Gate_flag, Puls_flag
      INTEGER, SAVE :: Secondoutflow_flag
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Currinsum(:), Pastin(:), Pastout(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Outflow_ts(:), Inflow_ts(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: C24(:, :), S24(:, :), Wvd(:, :)
!   Dimensions
      INTEGER, SAVE :: Mxnsos, Ngate, Nstage, Ngate2, Nstage2, Ngate3, Nstage3, Ngate4, Nstage4
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_2ndstflow
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Din1(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_outcfs(:), Lake_outcms(:), Lake_outvol(:), Lake_invol(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_sto(:), Lake_inflow(:), Lake_outflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_stream_in(:), Lake_lateral_inflow(:), Lake_gwflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_precip(:), Lake_sroff(:), Lake_interflow(:), Lake_outvol_ts(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_seep_in(:), Lake_evap(:), Lake_2gw(:), Lake_outq2(:)
!   Declared Parameters
      ! lake_segment_id only required if cascades are active, otherwise use hru_segment
      INTEGER, SAVE, ALLOCATABLE :: Obsout_lake(:), Lake_out2(:), Nsos(:), Ratetbl_lake(:), Lake_segment_id(:)
      REAL, SAVE, ALLOCATABLE :: Lake_qro(:), Lake_coef(:), Elev_outflow(:), Weir_coef(:), Weir_len(:)
      REAL, SAVE, ALLOCATABLE :: Lake_out2_a(:), Lake_out2_b(:), O2(:, :), S2(:, :)
      REAL, SAVE, ALLOCATABLE :: Lake_din1(:), Lake_init(:), Lake_vol_init(:)
      REAL, SAVE, ALLOCATABLE :: Rate_table(:, :), Rate_table2(:, :), Rate_table3(:, :), Rate_table4(:, :)
      REAL, SAVE, ALLOCATABLE :: Tbl_stage(:), Tbl_gate(:), Tbl_stage2(:), Tbl_gate2(:)
      REAL, SAVE, ALLOCATABLE :: Tbl_stage3(:), Tbl_gate3(:), Tbl_stage4(:), Tbl_gate4(:)
      END MODULE PRMS_MUSKINGUM_LAKE

!***********************************************************************
!     Main muskingum_lake routine
!***********************************************************************
      INTEGER FUNCTION muskingum_lake()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, SETDIMENS, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: muskingum_lake_decl, muskingum_lake_init, muskingum_lake_run, muskingum_lake_setdims
      EXTERNAL :: muskingum_lake_restart
!***********************************************************************
      muskingum_lake = 0

      IF ( Process_flag==RUN ) THEN
        muskingum_lake  = muskingum_lake_run()
      ELSEIF ( Process_flag==SETDIMENS ) THEN
        muskingum_lake = muskingum_lake_setdims()
      ELSEIF ( Process_flag==DECL ) THEN
        muskingum_lake  = muskingum_lake_decl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL muskingum_lake_restart(READ_INIT)
        muskingum_lake = muskingum_lake_init()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL muskingum_lake_restart(SAVE_INIT)
      ENDIF

      END FUNCTION muskingum_lake

!***********************************************************************
!     declares Lake routing specific dimensions
!***********************************************************************
      INTEGER FUNCTION muskingum_lake_setdims()
      USE PRMS_CONSTANTS, ONLY: MAXDIM
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: decldim
      EXTERNAL :: read_error
!***********************************************************************
      muskingum_lake_setdims = 0

      IF ( decldim('ngate', 0, MAXDIM, &
     &     'Maximum number of reservoir gate-opening values (columns) for lake rating table 1')/=0 ) &
     &     CALL read_error(7, 'ngate')
      IF ( decldim('nstage', 0, MAXDIM, &
     &     'Maximum number of lake elevations values (rows) for lake rating table 1')/=0 ) CALL read_error(7, 'nstage')

      IF ( decldim('ngate2', 0, MAXDIM, &
     &     'Maximum number of reservoir gate-opening values (columns) for lake rating table 2')/=0 ) &
     &     CALL read_error(7, 'ngate2')
      IF ( decldim('nstage2', 0, MAXDIM, &
     &     'Maximum number of lake elevations values (rows) for lake rating table 2')/=0 ) CALL read_error(7, 'nstage2')

      IF ( decldim('ngate3', 0, MAXDIM, &
     &     'Maximum number of reservoir gate-opening values (columns) for lake rating table 3')/=0 ) &
     &     CALL read_error(7, 'ngate3')
      IF ( decldim('nstage3', 0, MAXDIM, &
     &     'Maximum number of lake elevations values (rows) for lake rating table 3')/=0 ) CALL read_error(7, 'nstage3')

      IF ( decldim('ngate4', 0, MAXDIM, &
     &     'Maximum number of reservoir gate-opening values (columns) for lake rating table 4')/=0 ) &
     &     CALL read_error(7, 'ngate4')
      IF ( decldim('nstage4', 0, MAXDIM, &
     &     'Maximum number of lake elevations values (rows) for lake rating table 4')/=0 ) CALL read_error(7, 'nstage4')

      IF ( decldim('mxnsos', 0, MAXDIM, &
     &     'Maximum number of storage/outflow table values for storage-detention reservoirs and lakes connected to'// &
     &     ' the stream network using Puls routing')/=0 ) CALL read_error(7, 'mxnsos')

      END FUNCTION muskingum_lake_setdims

!***********************************************************************
!     muskingum_lake_decl - Declare parameters and variables and allocate arrays
!   Declared Parameters
!     tosegment, hru_segment, obsin_segment, K_coef, x_coef, segment_type
!     lake_type, lake_init, lake_qro, lake_din1, lake_coef, o2, s2, nsos, hru_area, lake_segment_id
!     tbl_stage, tbl_gate, lake_vol_init, rate_table, weir_coef, weir_len, elev_outflow, elevlake_init
!     lake_out2, lake_out2_a, lake_out2_b
!***********************************************************************
      INTEGER FUNCTION muskingum_lake_decl()
      USE PRMS_MUSKINGUM_LAKE
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar, getdim
      EXTERNAL :: read_error, print_module, error_stop
!***********************************************************************
      muskingum_lake_decl = 0

      CALL print_module(MODDESC, MODNAME, Version_muskingum_lake)

      ! Dimension for Puls routing
      Mxnsos = getdim('mxnsos')
      IF ( Mxnsos==-1 ) CALL read_error(1, 'mxnsos')
      IF ( Model==DOCUMENTATION .AND. Mxnsos<1 ) Mxnsos = 1

! Nlake_hrus set to nlake in call_modules for 5.0.0, dimension nlake_hrus will be in 5.0.1
! 5.0.0 assumes lakes are 1 HRU
      IF ( Mxnsos>0 ) ALLOCATE ( Wvd(Mxnsos, Nlake), S24(Mxnsos, Nlake), C24(Mxnsos, Nlake) )

      Ngate = 0
      Nstage = 0
      Ngate2 = 0
      Nstage2 = 0
      Ngate3 = 0
      Nstage3 = 0
      Ngate4 = 0
      Nstage4 = 0
      IF ( Model==DOCUMENTATION ) Nratetbl = 4
      IF ( Nratetbl>4 ) THEN
        PRINT *, 'dimension nratetbl specified as:', Nratetbl
        CALL error_stop('lake routing allows a maximum of 4 rating tables', ERROR_dim)
      ENDIF
      IF ( Nratetbl>0 ) THEN
        Ngate = getdim('ngate')
        IF ( Ngate==-1 ) CALL read_error(6, 'ngate')
        Nstage = getdim('nstage')
        IF ( Nstage==-1 ) CALL read_error(6, 'nstage')
        IF ( Nratetbl>1 ) THEN
          Ngate2 = getdim('ngate2')
          IF ( Ngate2==-1 ) CALL read_error(6, 'ngate2')
          Nstage2 = getdim('nstage2')
          IF ( Nstage2==-1 ) CALL read_error(6, 'nstage2')
          IF ( Nratetbl>2 ) THEN
            Ngate3 = getdim('ngate3')
            IF ( Ngate3==-1 ) CALL read_error(6, 'ngate3')
            Nstage3 = getdim('nstage3')
            IF ( Nstage3==-1 ) CALL read_error(6, 'nstage3')
            IF ( Nratetbl==4 ) THEN
              Ngate4 = getdim('ngate4')
              IF ( Ngate4==-1 ) CALL read_error(6, 'ngate4')
              Nstage4 = getdim('nstage4')
              IF ( Nstage4==-1 ) CALL read_error(6, 'nstage4')
            ENDIF
          ENDIF
        ENDIF
        IF ( Model==DOCUMENTATION ) THEN
          IF ( Nstage==0 ) Nstage = 1
          IF ( Ngate==0 ) Ngate = 1
          IF ( Nstage2==0 ) Nstage2 = 1
          IF ( Ngate2==0 ) Ngate2 = 1
          IF ( Nstage3==0 ) Nstage3 = 1
          IF ( Ngate3==0 ) Ngate3 = 1
          IF ( Nstage4==0 ) Nstage4 = 1
          IF ( Ngate4==0 ) Ngate4 = 1
        ELSE
          IF ( Nstage<1 .OR. Ngate<1 ) CALL error_stop('nratetbl>0 and nstage or ngate = 0', ERROR_dim)
        ENDIF
        IF ( Nratetbl>1 ) THEN
          IF ( Nstage2<1.OR.Ngate2<1 ) CALL error_stop('nratetbl>1 and nstage2 or ngate2 = 0', ERROR_dim)
        ENDIF
        IF ( Nratetbl>2 ) THEN
          IF ( Nstage3<1 .OR. Ngate3<1 ) CALL error_stop('nratetbl>2 and nstage3 or ngate3 = 0', ERROR_dim)
        ENDIF
        IF ( Nratetbl>3 ) THEN
          IF ( Nstage4<1 .OR. Ngate4<1 )  CALL error_stop('nratetbl>3 and nstage4 or ngate4 = 0', ERROR_dim)
        ENDIF
      ENDIF

      ALLOCATE ( Currinsum(Nsegment) )
      ALLOCATE ( Pastin(Nsegment), Pastout(Nsegment) )
      ALLOCATE ( Outflow_ts(Nsegment), Inflow_ts(Nsegment) )

      ! Lake declared variables
      ALLOCATE ( Lake_inflow(Nlake) )
      IF ( declvar(MODNAME, 'lake_inflow', 'nlake', Nlake, 'double', &
     &     'Total inflow to each lake', &
     &     'cfs', Lake_inflow)/=0 ) CALL read_error(3, 'lake_inflow')

      ALLOCATE ( Lake_outflow(Nlake) )
      IF ( declvar(MODNAME, 'lake_outflow', 'nlake', Nlake, 'double', &
     &     'Evaporation and seepage from each lake', &
     &     'cfs', Lake_outflow)/=0 ) CALL read_error(3, 'lake_outflow')

      ALLOCATE ( Lake_outcfs(Nlake) )
      IF ( declvar(MODNAME, 'lake_outcfs', 'nlake', Nlake, 'double', &
     &     'Streamflow leaving each lake, includes any second outlet flow', &
     &     'cfs', Lake_outcfs)/=0 ) CALL read_error(3, 'lake_outcfs')

      ALLOCATE ( Lake_outcms(Nlake) )
      IF ( declvar(MODNAME, 'lake_outcms', 'nlake', Nlake, 'double', &
     &     'Streamflow leaving each lake, includes any second outlet flow', &
     &     'cms', Lake_outcms)/=0 ) CALL read_error(3, 'lake_outcms')

! Declared Variables for Puls or linear routing
      ALLOCATE ( Lake_sto(Nlake) )
      IF ( declvar(MODNAME, 'lake_sto', 'nlake', Nlake, 'double', &
     &     'Storage in each lake using Puls or linear storage routing', &
     &     'cfs-days', Lake_sto)/=0 ) CALL read_error(3, 'lake_sto')

      ALLOCATE ( Din1(Nlake) )
      IF ( declvar(MODNAME, 'din1', 'nlake', Nlake, 'double', &
     &     'Inflow from the previous time step to each lake using Puls or linear storage routing', &
     &     'cfs', Din1)/=0 ) CALL read_error(3, 'din1')

      ALLOCATE ( Lake_stream_in(Nlake) )
      IF ( declvar(MODNAME, 'lake_stream_in', 'nlake', Nlake, 'double', &
     &     'Total streamflow into each lake', &
     &     'cfs', Lake_stream_in)/=0 ) CALL read_error(3, 'lake_stream_in')

      ALLOCATE ( Lake_precip(Nlake) )
      IF ( declvar(MODNAME, 'lake_precip', 'nlake', Nlake, 'double', &
     &     'Total precipitation into each lake', &
     &     'cfs', Lake_precip)/=0 ) CALL read_error(3, 'lake_precip')

      IF ( Cascade_flag>CASCADE_OFF .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Lake_lateral_inflow(Nlake) )
        IF ( declvar(MODNAME, 'lake_lateral_inflow', 'nlake', Nlake, 'double', &
     &       'Lateral inflow to each lake', &
     &       'cfs', Lake_lateral_inflow)/=0 ) CALL read_error(3, 'lake_lateral_inflow')
        ALLOCATE ( Lake_sroff(Nlake) )
        IF ( declvar(MODNAME, 'lake_sroff', 'nlake', Nlake, 'double', &
     &       'Total surface runoff into each lake', &
     &       'cfs', Lake_sroff)/=0 ) CALL read_error(3, 'lake_sroff')
        ALLOCATE ( Lake_interflow(Nlake) )
        IF ( declvar(MODNAME, 'lake_interflow', 'nlake', Nlake,'double', &
     &       'Total interflow into each lake', &
     &       'cfs', Lake_interflow)/=0 ) CALL read_error(3, 'lake_interflow')
        ALLOCATE ( Lake_gwflow(Nlake) )
        IF ( declvar(MODNAME, 'lake_gwflow', 'nlake', Nlake,'double', &
     &       'Total groundwater flow into each lake', &
     &       'cfs', Lake_gwflow)/=0 ) CALL read_error(3, 'lake_gwflow')
      ENDIF

      ALLOCATE ( Lake_evap(Nlake) )
      IF ( declvar(MODNAME, 'lake_evap', 'nlake', Nlake, 'double', &
     &     'Total evaporation from each lake', &
     &     'cfs', Lake_evap)/=0 ) CALL read_error(3, 'lake_evap')

! Declared Variables for broad-crested weir or gate opening routing
      ALLOCATE ( Lake_2gw(Nlake) )
      IF ( declvar(MODNAME, 'lake_2gw', 'nlake', Nlake, 'double', &
     &     'Total seepage from each lake using broad-crested weir or gate opening routing', &
     &     'cfs', Lake_2gw)/=0 ) CALL read_error(3, 'lake_2gw')

      ALLOCATE ( Lake_seep_in(Nlake) )
      IF ( declvar(MODNAME, 'lake_seep_in', 'nlake', Nlake, 'double', &
     &     'Total seepage into each lake using broad-crested weir or gate opening routing', &
     &     'cfs', Lake_seep_in)/=0 ) CALL read_error(3, 'lake_seep_in')

      ALLOCATE ( Lake_invol(Nlake) )
      IF ( declvar(MODNAME, 'lake_invol', 'nlake', Nlake, 'double', &
     &     'Inflow to each lake using broad-crested weir or gate opening routing', &
     &     'acre-feet', Lake_invol)/=0 ) CALL read_error(3, 'lake_invol')

! Declared Variables for gate opening routing
      ALLOCATE ( Lake_outvol(Nlake) )
      IF ( declvar(MODNAME, 'lake_outvol', 'nlake', Nlake, 'double', &
     &     'Outflow from each lake using broad-crested weir or gate opening routing', &
     &     'acre-inches', Lake_outvol)/=0 ) CALL read_error(3, 'lake_outvol')

      ALLOCATE ( Lake_outvol_ts(Nlake) )
      IF ( declvar(MODNAME, 'lake_outvol_ts', 'nlake', Nlake, 'double', &
     &     'Outflow from each lake using broad-crested weir or gate opening routing for the time step', &
     &     'acre-inches', Lake_outvol_ts)/=0 ) CALL read_error(3, 'lake_outvol_ts')

! Declared Variables for lakes with a second outlet and gate opening routing
      IF ( Nratetbl>0 .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'basin_2ndstflow', 'one', 1, 'double', &
     &       'Basin volume-weighted average streamflow from each lake with a second outlet', &
     &       'inches', Basin_2ndstflow)/=0 ) CALL read_error(3, 'basin_2ndstflow')
        ALLOCATE ( Lake_outq2(Nlake) )
        IF ( declvar(MODNAME, 'lake_outq2', 'nlake', Nlake, 'double', &
     &       'Streamflow from second outlet for each lake with a second outlet', &
     &       'cfs', Lake_outq2)/=0 ) CALL read_error(3, 'lake_outq2')
      ENDIF

! Declared Parameters
      ALLOCATE ( Lake_segment_id(Nsegment) )
      IF ( Cascade_flag>CASCADE_OFF ) THEN
        IF ( declparam(MODNAME, 'lake_segment_id', 'nsegment', 'integer', &
     &       '0', 'bounded', 'nlake', &
     &       'Index of lake associated with a segment', &
     &       'Index of lake associated with a segment', &
     &       'none')/=0 ) CALL read_error(1, 'lake_segment_id')
      ENDIF

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
        ALLOCATE ( Lake_qro(Nlake) )
        IF ( declparam(MODNAME, 'lake_qro', 'nlake', 'real', &
     &       '0.1', '0.0', '1.0E7', &
     &       'Initial daily mean outflow from each lake', &
     &       'Initial daily mean outflow from each lake', &
     &       'cfs')/=0 ) CALL read_error(1, 'lake_qro')

! Declared Parameters for Puls or linear routing
        ALLOCATE ( Lake_init(Nlake) )
        IF ( declparam(MODNAME, 'lake_init', 'nlake', 'real', &
     &       '0.0', '0.0', '1.0E7', &
     &       'Initial storage in each lake', &
     &       'Initial storage in each lake using Puls or linear storage routing', &
     &       'cfs-days')/=0 ) CALL read_error(1, 'lake_init')

        ALLOCATE ( Lake_din1(Nlake) )
        IF ( declparam(MODNAME, 'lake_din1', 'nlake', 'real', &
     &       '0.1', '0.0', '1.0E7', &
     &       'Initial inflow to each lake', &
     &       'Initial inflow to each lake using Puls or linear storage routing', &
     &       'cfs')/=0 ) CALL read_error(1, 'lake_din1')
      ENDIF

! Declared Parameters for linear routing
      ALLOCATE ( Lake_coef(Nlake) )
      IF ( declparam(MODNAME, 'lake_coef', 'nlake', 'real', &
     &     '0.1', '0.0001', '1.0', &
     &     'Linear lake routing coefficient', &
     &     'Coefficient in equation to route storage to streamflow for each lake using linear routing', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'lake_coef')

! Declared Parameters for Puls routing
      IF ( Mxnsos>0 ) THEN
        ALLOCATE ( O2(Mxnsos, Nlake) )
        IF ( declparam(MODNAME, 'o2', 'mxnsos,nlake', 'real', &
     &       '0.0', '0.0', '1.0E7', &
     &       'Outflow values in outflow/storage tables for Puls routing', &
     &       'Outflow values in outflow/storage tables for each lake using Puls routing', &
     &       'cfs')/=0 ) CALL read_error(1, 'o2')

        ALLOCATE ( S2(Mxnsos, Nlake) )
        IF ( declparam(MODNAME, 's2', 'mxnsos,nlake', 'real', &
     &       '0.0', '0.0', '1.0E7', &
     &       'Storage values in outflow/storage tables for Puls routing', &
     &       'Storage values in outflow/storage table for each lake using Puls routing', &
     &       'cfs-days')/=0 ) CALL read_error(1, 's2')

        ALLOCATE ( Nsos(Nlake) )
        IF ( declparam(MODNAME, 'nsos', 'nlake', 'integer', &
     &       '0', 'bounded', 'mxnsos', &
     &       'Number of storage/outflow values in table for Puls routing', &
     &       'Number of storage/outflow values in table for each lake using Puls routing', &
     &       'none')/=0 ) CALL read_error(1, 'nsos')
      ENDIF

! Declared Parameters for broad-crested weir or gate opening routing
      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
        ALLOCATE ( Lake_vol_init(Nlake) )
        IF ( declparam(MODNAME, 'lake_vol_init', 'nlake', 'real', &
     &       '0.0', '0.0', '1.0E7', &
     &       'Initial lake volume', &
     &       'Initial lake volume for each lake using broad-crested weir or gate opening routing', &
     &       'acre-feet')/=0 ) CALL read_error(1, 'lake_vol_init')
      ENDIF

! Declared Parameters for broad-crested weir routing
      ALLOCATE ( Weir_coef(Nlake) )
      IF ( declparam(MODNAME, 'weir_coef', 'nlake', 'real', &
     &     '2.7', '2.0', '3.0', &
     &     'Broad-crested weir coefficent', &
     &     'Coefficient for lakes using broad-crested weir routing', &
     &     'none')/=0 ) CALL read_error(1, 'weir_coef')

      ALLOCATE ( Weir_len(Nlake) )
      IF ( declparam(MODNAME, 'weir_len', 'nlake', 'real', &
     &     '5.0', '1.0', '1000.0', &
     &     'Broad-crested weir length', &
     &     'Weir length for lakes using broad-crested weir routing', &
     &     'feet')/=0 ) CALL read_error(1, 'weir_len')

      ALLOCATE ( Elev_outflow(Nlake) )
      IF ( declparam(MODNAME, 'elev_outflow', 'nlake', 'real', &
     &     '0.0', '-300.0', '10000.0', &
     &     'Elevation of the main outflow point', &
     &     'Elevation of the main outflow point for each lake using broad-crested weir routing', &
     &     'feet')/=0 ) CALL read_error(1, 'elev_outflow')

! Declared Parameters for gate opening routing
      IF ( Nratetbl>0 ) THEN
        ALLOCATE ( Ratetbl_lake(Nratetbl), Rate_table(Nstage,Ngate), Tbl_stage(Nstage), Tbl_gate(Ngate) )
        IF ( declparam(MODNAME, 'ratetbl_lake', 'nratetbl', 'integer', &
     &       '0', 'bounded', 'nlake', &
     &       'Index of lake associated with each rating table', &
     &       'Index of lake associated with each rating table for'// &
     &       ' each lake using gate opening routing', &
     &       'none')/=0 ) CALL read_error(1, 'ratetbl_lake')
        IF ( declparam(MODNAME, 'rate_table', 'nstage,ngate', 'real', &
     &       '5.0', '-100.0', '1000.0', &
     &       'Rating table 1 with stage (rows) and gate opening (cols)', &
     &       'Rating table with stage (rows) and gate opening'// &
     &       ' (cols) for rating table 1 for lakes using gate opening routing and nratetbl>0', &
     &       'cfs')/=0 ) CALL read_error(1, 'rate_table')
        IF ( declparam(MODNAME, 'tbl_stage', 'nstage', 'real', &
     &       '5.0', '-100.0', '1000.0', &
     &       'Stage values for each row of rating table 1', &
     &       'Stage values for each row for rating table 1 for lakes using gate opening routing and nratetbl>0', &
     &       'feet')/=0 ) CALL read_error(1, 'tbl_stage')
        IF ( declparam(MODNAME, 'tbl_gate', 'ngate', 'real', &
     &       '0.0', '0.0', '20.0', &
     &       'Gate openings for each column of rating table 1', &
     &       'Gate openings for each column for rating table 1 for lakes using gate opening routing and nratetbl>0', &
     &       'inches')/=0 ) CALL read_error(1, 'tbl_gate')

        IF ( Nratetbl>1 ) THEN
          ALLOCATE ( Rate_table2(Nstage2,Ngate2), Tbl_stage2(Nstage2), Tbl_gate2(Ngate2) )
          IF ( declparam(MODNAME, 'rate_table2', 'nstage2,ngate2', 'real', &
     &         '5.0', '-100.0', '1000.0', &
     &         'Rating table 2 with stage (rows) and gate opening (cols)', &
     &         'Rating table with stage (rows) and gate opening'// &
     &         ' (cols) for rating table 2 for lakes using gate opening routing and nratetbl>1', &
     &         'cfs')/=0 ) CALL read_error(1, 'rate_table2')
          IF ( declparam(MODNAME, 'tbl_stage2', 'nstage2', 'real', &
     &         '5.0', '-100.0', '1000.0', &
     &         'Stage values for each row of rating table 2', &
     &         'Stage values for each row for rating table 2 for lakes using gate opening routing and nratetbl>1', &
     &         'feet')/=0 ) CALL read_error(1, 'tbl_stage2')
          IF ( declparam(MODNAME, 'tbl_gate2', 'ngate2', 'real', &
     &         '0.0', '0.0', '20.0', &
     &         'Gate openings for each column of rating table 2', &
     &         'Gate openings for each column for rating table 2 for lakes using gate opening routing and nratetbl>1', &
     &         'inches')/=0 ) CALL read_error(1, 'tbl_gate2')

          IF ( Nratetbl>2 ) THEN
            ALLOCATE ( Rate_table3(Nstage3,Ngate3), Tbl_stage3(Nstage3), Tbl_gate3(Ngate3) )
            IF ( declparam(MODNAME, 'rate_table3', 'nstage3,ngate3', 'real', &
     &           '5.0', '-100.0', '1000.0', &
     &           'Rating table 3 with stage (rows) and gate opening (cols)', &
     &           'Rating table with stage (rows) and gate opening'// &
     &           ' (cols) for rating table 3 for lakes using gate opening routing and nratetbl>2', &
     &           'cfs')/=0 ) CALL read_error(1, 'rate_table3')
            IF ( declparam(MODNAME, 'tbl_stage3', 'nstage3', 'real', &
     &           '5.0', '-100.0', '1000.0', &
     &           'Stage values for each row of rating table 3', &
     &           'Stage values for each row for rating table 3 for lakes using gate opening routing and nratetbl>2', &
     &           'feet')/=0 ) CALL read_error(1, 'tbl_stage3')
            IF ( declparam(MODNAME, 'tbl_gate3', 'ngate3', 'real', &
     &           '0.0', '0.0', '20.0', &
     &           'Gate openings for each column of rating table 3', &
     &           'Gate openings for each column for rating table 3 for lakes using gate opening routing and nratetbl>2', &
     &           'inches')/=0 ) CALL read_error(1, 'tbl_gate3')

            IF ( Nratetbl>3 ) THEN
              ALLOCATE ( Rate_table4(Nstage4,Ngate4), Tbl_stage4(Nstage4), Tbl_gate4(Ngate4) )
              IF ( declparam(MODNAME, 'rate_table4', 'nstage4,ngate4', 'real', &
     &             '5.0', '-100.0', '1000.0', &
     &             'Rating table 4 with stage (rows) and gate opening (cols)', &
     &             'Rating table with stage (rows) and gate opening'// &
     &             ' (cols) for rating table 4 for lakes using gate opening routing and nratetbl>3', &
     &             'cfs')/=0 ) CALL read_error(1, 'rate_table4')
              IF ( declparam(MODNAME, 'tbl_stage4', 'nstage4', 'real', &
     &             '5.0', '-100.0', '1000.0', &
     &             'Stage values for each row of rating table 4', &
     &             'Stage values for each row for rating table 4 for lakes using gate opening routing and nratetbl>3', &
     &             'feet')/=0 ) CALL read_error(1, 'tbl_stage4')
              IF ( declparam(MODNAME, 'tbl_gate4', 'ngate4', 'real', &
     &             '0.0', '0.0', '20.0', &
     &             'Gate openings for each column of rating table 4', &
     &             'Gate openings for each column for rating table 4 for lakes using gate opening routing and nratetbl>3', &
     &             'inches')/=0 ) CALL read_error(1, 'tbl_gate4')
            ENDIF
          ENDIF
        ENDIF
      ENDIF

! Declared Parameters for lakes with lake outflow set to measured streamflow
      ALLOCATE ( Obsout_lake(Nlake) )
      IF ( declparam(MODNAME, 'obsout_lake', 'nlake', 'integer', &
     &     '0', 'bounded', 'nobs', &
     &     'Index of streamflow measurement station that specifies outflow from a lake', &
     &     'Index of streamflow measurement station that specifies outflow from each lake using measured flow replacement', &
     &     'none')/=0 ) CALL read_error(1, 'obsout_lake')

      IF ( Nratetbl>0 ) THEN
! Declared Parameters for lakes with a second outlet and gate opening routing
        ALLOCATE ( Lake_out2(Nlake) )
        IF ( declparam(MODNAME, 'lake_out2', 'nlake', 'integer', &
     &       '0', '0', '1', &
     &       'Switch to specify a second outlet from a lake', &
     &       'Switch to specify a second outlet from each lake using gate opening routing (0=no; 1=yes)', &
     &       'none')/=0 ) CALL read_error(1, 'lake_out2')

        ALLOCATE ( Lake_out2_a(Nlake) )
        IF ( declparam(MODNAME, 'lake_out2_a', 'nlake', 'real', &
     &       '1.0', '0.0', '10000.0', &
     &       'Outflow coefficient A for each lake with second outlet', &
     &       'Coefficient A in outflow equation for each lake with a second outlet using gate opening routing', &
     &       'cfs/ft')/=0 ) CALL read_error(1, 'lake_out2_a')

        ALLOCATE ( Lake_out2_b(Nlake) )
        IF ( declparam(MODNAME, 'lake_out2_b', 'nlake', 'real', &
     &       '100.0', '0.0', '10000.0', &
     &       'Outflow coefficient A for each lake with second outlet', &
     &       'Coefficient B in outflow equation for each lake with a second outlet using gate opening routing', &
     &       'cfs')/=0 ) CALL read_error(1, 'lake_out2_b')
      ENDIF

      END FUNCTION muskingum_lake_decl

!***********************************************************************
!    muskingum_lake_init - Get and check parameter values and initialize variables
!***********************************************************************
      INTEGER FUNCTION muskingum_lake_init()
      USE PRMS_MUSKINGUM_LAKE
      USE PRMS_BASIN, ONLY: Basin_area_inv, Active_hrus, Hru_route_order, Gwr_type, &
     &    Lake_hru_id, Weir_gate_flag, Lake_type, Puls_lin_flag
      USE PRMS_FLOWVARS, ONLY: Seg_outflow, Basin_lake_stor, Lake_vol
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_ROUTING, ONLY: Basin_segment_storage, Segment_type, Hru_segment
      IMPLICIT NONE
! Functions
      INTRINSIC :: ABS, NINT, DBLE, DABS
      EXTERNAL :: read_error, error_stop
      INTEGER, EXTERNAL :: getparam
! Local Variables
      INTEGER :: i, ierr, j, jj, kk, ii, jjj
      DOUBLE PRECISION :: tmp
!***********************************************************************
      muskingum_lake_init = 0

      IF ( Init_vars_from_file==0 ) Outflow_ts = 0.0D0

      Basin_segment_storage = 0.0D0
      DO i = 1, Nsegment
        Basin_segment_storage = Basin_segment_storage + Seg_outflow(i)
      ENDDO
      Basin_segment_storage = Basin_segment_storage*Basin_area_inv/Cfs_conv

! Weir_gate_flag and Puls_lin_flag set in basin as needed for restart and gwflow module
      Obs_flag = OFF
      Linear_flag = OFF
      Weir_flag = OFF
      Gate_flag = OFF
      Puls_flag = OFF
      DO i = 1, Nlake
        IF ( Lake_type(i)==1 ) THEN
           Puls_flag = ACTIVE
        ELSEIF ( Lake_type(i)==2 ) THEN
          Linear_flag = ACTIVE
        ELSEIF ( Lake_type(i)==4 ) THEN
          Weir_flag = ACTIVE
        ELSEIF ( Lake_type(i)==5 ) THEN
          Gate_flag = ACTIVE
        ELSEIF ( Lake_type(i)==6 ) THEN
          Obs_flag = ACTIVE
        ELSEIF ( Lake_type(i)/=3 ) THEN
          PRINT *, 'ERROR, invalid lake_type for lake:', i, Lake_type(i)
          Inputerror_flag = 1
        ENDIF
      ENDDO

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
        IF ( getparam(MODNAME, 'lake_qro', Nlake, 'real', Lake_qro)/=0 ) CALL read_error(2, 'lake_qro')
        DO j = 1, Nlake
          Lake_outcfs(j) = Lake_qro(j)
          Lake_outcms(j) = Lake_qro(j)*CFS2CMS_CONV
        ENDDO
      ENDIF

      Lake_outvol = 0.0D0
      Lake_outvol_ts = 0.0D0
      Lake_invol = 0.0D0
      Lake_precip = 0.0D0
      Lake_seep_in = 0.0D0
      Lake_evap = 0.0D0
      Lake_2gw = 0.0D0
      Lake_inflow = 0.0D0
      Lake_outflow = 0.0D0
      IF ( Gate_flag==ACTIVE ) Lake_outq2 = 0.0D0
      Basin_2ndstflow = 0.0D0
      Lake_stream_in = 0.0D0
      Basin_lake_stor = 0.0D0
      IF ( Cascade_flag>CASCADE_OFF ) THEN
        Lake_lateral_inflow = 0.0D0
        Lake_sroff = 0.0D0
        Lake_interflow = 0.0D0
        Lake_gwflow = 0.0D0
      ENDIF

      IF ( Cascade_flag==CASCADE_OFF .OR. Cascade_flag==CASCADE_HRU_SEGMENT ) THEN
        Lake_segment_id = 0
        DO jjj = 1, Active_hrus
          j = Hru_route_order(jjj)
          i = Hru_segment(j)
          IF ( i>0 ) THEN
            IF ( Segment_type(i)==2 ) Lake_segment_id(i) = Lake_hru_id(j)
          ENDIF
        ENDDO
      ELSE
        ! if cascades are active, must input new parameter lake_segment_id
        IF ( getparam(MODNAME, 'lake_segment_id', Nsegment, 'integer', Lake_segment_id)/=0 ) CALL read_error(2, 'lake_segment_id')
      ENDIF
      DO j = 1, Nsegment
        IF ( Lake_segment_id(j)>0 .AND. Segment_type(j)/=2 ) THEN
          PRINT *, 'ERROR, segment_type not equal to 2 when the segment is associated with a lake'
          PRINT *, '       segment:', j, ' lake:', Lake_segment_id(j)
          Inputerror_flag = 1
        ENDIF
        IF ( Lake_segment_id(j)==0 .AND. Segment_type(j)==2 ) THEN
          PRINT *, 'ERROR, segment_type equals 2 when the segment is not associated with a lake'
          PRINT *, '       segment:', j, ' lake:', Lake_segment_id(j)
          Inputerror_flag = 1
        ENDIF
      ENDDO

      Secondoutflow_flag = OFF
      IF ( Gate_flag==ACTIVE ) THEN
        IF ( Nratetbl<1 ) CALL error_stop('nratetbl = 0 and gate opening routing requested', ERROR_dim)
        IF ( getparam(MODNAME, 'rate_table', Nstage*Ngate, 'real', Rate_table)/=0 ) CALL read_error(2, 'rate_table')
        IF ( getparam(MODNAME, 'tbl_stage', Nstage, 'real', Tbl_stage)/=0 ) CALL read_error(2, 'tbl_stage')
        IF ( getparam(MODNAME, 'tbl_gate', Ngate, 'real', Tbl_gate)/=0 ) CALL read_error(2, 'tbl_gate')
        IF ( getparam(MODNAME, 'ratetbl_lake', Nratetbl, 'integer', Ratetbl_lake)/=0 ) CALL read_error(2, 'ratetbl_lake')
        IF ( Gate_flag==1 ) THEN
          IF ( getparam(MODNAME, 'lake_out2', Nlake, 'integer', Lake_out2)/=0  ) CALL read_error(2, 'lake_out2')
          DO j = 1, Nlake
            IF ( Lake_out2(j)==1 ) Secondoutflow_flag = ACTIVE
          ENDDO
          IF ( Secondoutflow_flag==ACTIVE ) THEN
            IF ( getparam(MODNAME, 'lake_out2_a', Nlake, 'real', Lake_out2_a)/=0 ) CALL read_error(2, 'lake_out2_a')
            IF ( getparam(MODNAME, 'lake_out2_b', Nlake, 'real', Lake_out2_b)/=0 ) CALL read_error(2, 'lake_out2_b')
          ENDIF
        ENDIF

        IF ( Nratetbl>1 ) THEN
          IF ( getparam(MODNAME, 'rate_table2', Nstage2*Ngate2, 'real', Rate_table2)/=0 ) CALL read_error(2, 'rate_table2')
          IF ( getparam(MODNAME, 'tbl_stage2', Nstage2, 'real', Tbl_stage2)/=0 ) CALL read_error(2, 'tbl_stage2')
          IF ( getparam(MODNAME, 'tbl_gate2', Ngate2, 'real', Tbl_gate2)/=0 ) CALL read_error(2, 'tbl_gate2')

          IF ( Nratetbl>2 ) THEN
            IF ( getparam(MODNAME, 'rate_table3', Nstage3*Ngate3, 'real', Rate_table3)/=0 ) &
     &           CALL read_error(2, 'rate_table3')
            IF ( getparam(MODNAME, 'tbl_stage3', Nstage3, 'real', Tbl_stage3)/=0 ) CALL read_error(2, 'tbl_stage3')
            IF ( getparam(MODNAME, 'tbl_gate3', Ngate3, 'real', Tbl_gate3)/=0 ) CALL read_error(2, 'tbl_gate3')

            IF ( Nratetbl>3 ) THEN
              IF ( getparam(MODNAME, 'rate_table4', Nstage4*Ngate4, 'real', Rate_table4)/=0 ) &
     &             CALL read_error(2, 'rate_table4')
              IF ( getparam(MODNAME, 'tbl_stage4', Nstage4, 'real', Tbl_stage4)/=0 ) CALL read_error(2, 'tbl_stage4')
              IF ( getparam(MODNAME, 'tbl_gate4', Ngate4, 'real', Tbl_gate4)/=0 ) CALL read_error(2, 'tbl_gate4')
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      IF ( Puls_lin_flag==ACTIVE ) THEN
        IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
          IF ( getparam(MODNAME, 'lake_init', Nlake, 'real', Lake_init)/=0 ) CALL read_error(2, 'lake_init')
          IF ( getparam(MODNAME, 'lake_din1', Nlake, 'real', Lake_din1)/=0 ) CALL read_error(2, 'lake_din1')
          DO i = 1, Nlake
            Lake_sto(i) = DBLE( Lake_init(i) )
            Din1(i) = DBLE( Lake_din1(i) )
          ENDDO
        ENDIF
        DO i = 1, Nlake
          IF ( Lake_type(i)==1 ) THEN
            kk = Nsos(i)
            IF ( kk<1 ) THEN
              PRINT *, 'ERROR, lake_type = 1, but, nsos<1, lake:', i, ' nsos:', kk, ' mxnsos:', Mxnsos
              Inputerror_flag = 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF ( Puls_flag==ACTIVE ) THEN
        IF ( Mxnsos==0 ) CALL error_stop('dimension mxnsos = 0 and Puls routing requested', ERROR_dim)
        IF ( getparam(MODNAME, 'o2', Mxnsos*Nlake, 'real', O2)/=0 ) CALL read_error(2, 'o2')
        IF ( getparam(MODNAME, 's2', Mxnsos*Nlake, 'real', S2)/=0 ) CALL read_error(2, 's2')
        IF ( getparam(MODNAME, 'nsos', Nlake, 'integer', Nsos)/=0 ) CALL read_error(2, 'nsos')
      ENDIF

      IF ( Linear_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'lake_coef', Nlake, 'real', Lake_coef)/=0 ) CALL read_error(2, 'lake_coef')
      ENDIF

      IF ( Weir_gate_flag==ACTIVE ) THEN
        IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
          IF ( getparam(MODNAME, 'lake_vol_init', Nlake, 'real', Lake_vol_init)/=0 ) CALL read_error(2, 'lake_vol_init')
          DO i = 1, Nlake
            Lake_vol(i) = DBLE( Lake_vol_init(i) )
          ENDDO
        ENDIF
      ENDIF
      DO j = 1, Nhru
        IF ( Gwr_type(j)==LAKE ) THEN
          jjj = Lake_hru_id(j)
          IF ( jjj==0 ) THEN
            PRINT *, 'ERROR, GWR specified as a lake but lake_hru_id value = 0, GWR:', j
            Inputerror_flag = 1
          ENDIF
        ENDIF
        IF ( Lake_hru_id(j)>0 .AND. Gwr_type(j)/=LAKE ) THEN
          PRINT *, 'ERROR, GWR specified as associated with a lake but gwr_type = 0, GWR:', j
          Inputerror_flag = 1
        ENDIF
      ENDDO

      IF ( Weir_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'weir_coef', Nlake, 'real', Weir_coef)/=0 ) CALL read_error(2, 'weir_coef')
        IF ( getparam(MODNAME, 'weir_len', Nlake, 'real', Weir_len)/=0 ) CALL read_error(2, 'weir_len')
        IF ( getparam(MODNAME, 'elev_outflow', Nlake, 'real', Elev_outflow)/=0 ) CALL read_error(2, 'elev_outflow')
      ENDIF

      IF ( Obs_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'obsout_lake', Nlake, 'integer', Obsout_lake)/=0 ) CALL read_error(2, 'obsout_lake')
      ELSE
        Obsout_lake = 1
      ENDIF

      DO j = 1, Nlake
        ierr = 0
        IF ( Lake_type(j)==1 .OR. Lake_type(j)==2 ) THEN
          IF ( Lake_type(j)==1 ) THEN
            kk = Nsos(j)
            IF ( kk<1 ) THEN
              PRINT *, 'ERROR, lake_type = 1, but, nsos<1, lake:', j, ' nsos:', kk, ' mxnsos:', Mxnsos
              ierr = 1
            ENDIF
          ENDIF
!        ELSEIF ( Weir_gate_flag==ACTIVE ) THEN
!          IF ( Lake_type(j)==4 ) THEN
!            IF ( Elev_outflow(j)<0.0 ) THEN
!              PRINT *, 'ERROR, elev_outflow < 0.0 for lake:', j, Elev_outflow(j)
!              ierr = 1
!            ENDIF
!          ENDIF
        ELSEIF ( Lake_type(j)==6 ) THEN
          IF ( Obsout_lake(j)==0 ) THEN
            PRINT *, 'ERROR, obsout_lake value = 0 for lake:', j, Obsout_lake(j)
            ierr = 1
          ENDIF
        ENDIF
        IF ( ierr==1 ) THEN
          Inputerror_flag = 1
          CYCLE
        ENDIF
        IF ( Lake_type(j)==1 ) THEN
          kk = Nsos(j)
          DO ii = 1, kk
            Wvd(ii, j) = DBLE( S2(ii, j) + O2(ii, j)*0.5 )
          ENDDO
          DO jj = 2, kk
            tmp = Wvd(jj, j) - Wvd(jj-1, j)
            IF ( DABS(tmp)<DNEARZERO ) tmp = 0.00001D0
            S24(jj, j) = DBLE( (O2(jj, j)-O2(jj-1, j)) )/tmp
            C24(jj, j) = DBLE( O2(jj, j) ) - S24(jj, j)*Wvd(jj, j)
          ENDDO
        ELSEIF ( Weir_gate_flag==ACTIVE ) THEN
          Basin_lake_stor = Basin_lake_stor + Lake_vol(j)*12.0D0
        ENDIF
      ENDDO
      Basin_lake_stor = Basin_lake_stor*Basin_area_inv

      !DEALLOCATE ( Lake_init, Lake_vol_init )
      !DEALLOCATE ( Lake_din1, Lake_qro )
      IF ( Mxnsos>0 ) DEALLOCATE ( O2, S2 )

      END FUNCTION muskingum_lake_init

!***********************************************************************
!     muskingum_lake_run - Compute routing summary values
!***********************************************************************
      INTEGER FUNCTION muskingum_lake_run()
      USE PRMS_MUSKINGUM_LAKE
      USE PRMS_BASIN, ONLY: Basin_area_inv, Hru_route_order, Active_hrus, &
     &    Lake_area, Lake_type, Hru_area_dble, Lake_hru_id, Hru_type, Weir_gate_flag, &
     &    Basin_gl_cfs, Basin_gl_ice_cfs
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_cms, Basin_gwflow_cfs, Basin_ssflow_cfs, &
     &    Basin_stflow_out, Basin_cfs, Basin_stflow_in, Basin_sroff_cfs, Seg_inflow, Seg_outflow, &
     &    Seg_upstream_inflow, Seg_lateral_inflow, Flow_out, Basin_lake_stor, Hru_actet, Lake_vol, Basin_sroff
      USE PRMS_OBS, ONLY: Streamflow_cfs
      USE PRMS_SET_TIME, ONLY: Cfs_conv, Nowyear, Nowmonth, Nowday
      USE PRMS_WATER_USE, ONLY: Lake_transfer, Lake_gain
      USE PRMS_ROUTING, ONLY: Use_transfer_segment, Segment_delta_flow, Basin_segment_storage, &
     &    Obsin_segment, Segment_order, Tosegment, C0, C1, C2, Ts, Ts_i, Obsout_segment, &
     &    Flow_to_ocean, Flow_to_great_lakes, Flow_out_region, Flow_out_NHM, Segment_type, Flow_terminus, &
     &    Flow_to_lakes, Flow_replacement, Flow_in_region, Flow_in_nation, Flow_headwater, Flow_in_great_lakes
      USE PRMS_SRUNOFF, ONLY: Hortonian_lakes
      USE PRMS_GLACR, ONLY: Basin_gl_top_melt, Basin_gl_ice_melt
      USE PRMS_SOILZONE, ONLY: Upslope_dunnianflow, Upslope_interflow
      USE PRMS_GWFLOW, ONLY: Basin_gwflow, Lake_seepage, Gw_seep_lakein, Gw_upslope
      IMPLICIT NONE
! Functions
      INTRINSIC :: MOD, DBLE
      EXTERNAL :: route_lake, error_stop
! Local Variables
      INTEGER :: i, j, iorder, toseg, imod, tspd, segtype, lakeid, k, jj
      DOUBLE PRECISION :: area_fac, segout, currin, tocfs, lake_in_ts
!***********************************************************************
      muskingum_lake_run = 0

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

      IF ( Secondoutflow_flag==ACTIVE ) THEN
        Basin_2ndstflow = 0.0D0
        Lake_outq2 = 0.0D0
      ENDIF
      Basin_lake_stor = 0.0D0
      Lake_inflow = 0.0D0
      Lake_outflow = 0.0D0
      Lake_stream_in = 0.0D0
      Lake_precip = 0.0D0
      IF ( Cascade_flag>CASCADE_OFF ) THEN
        Lake_lateral_inflow = 0.0D0
        Lake_sroff = 0.0D0
        Lake_interflow = 0.0D0
        Lake_gwflow = 0.0D0
      ENDIF
      IF ( Weir_gate_flag==ACTIVE ) THEN
        Lake_seep_in = 0.0D0
        Lake_2gw = 0.0D0
      ENDIF
      Lake_evap = 0.0D0
      ! shouldn't have snowpack, all precipitation should be added directly to lake
      ! units of lake_inflow = cfs
      DO jj = 1, Active_hrus
        k = Hru_route_order(jj)
        IF ( Hru_type(k)/=LAKE ) CYCLE
        tocfs = Hru_area_dble(k)*Cfs_conv
        lakeid = Lake_hru_id(k)
        Lake_precip(lakeid) = Lake_precip(lakeid) + tocfs*DBLE(Hru_ppt(k))
        IF ( Cascade_flag>CASCADE_OFF ) THEN
          Lake_sroff(lakeid) = Lake_sroff(lakeid) + tocfs*(Hortonian_lakes(k)+Upslope_dunnianflow(k))
          Lake_interflow(lakeid) = Lake_interflow(lakeid) + tocfs*Upslope_interflow(k)
          Lake_gwflow(lakeid) = Lake_gwflow(lakeid) + tocfs*Gw_upslope(k)
        ENDIF
        Lake_evap(lakeid) = Lake_evap(lakeid) + tocfs*Hru_actet(k)
      ENDDO
      DO lakeid = 1, Nlake
        Lake_inflow(lakeid) = Lake_precip(lakeid)
        IF ( Cascade_flag>CASCADE_OFF ) THEN
          Lake_lateral_inflow(lakeid) = Lake_sroff(lakeid) + Lake_interflow(lakeid) + Lake_gwflow(lakeid)
          Lake_inflow(lakeid) = Lake_inflow(lakeid) + Lake_lateral_inflow(lakeid)
        ENDIF
        IF ( Lake_add_water_use==ACTIVE ) THEN
          IF ( Lake_gain(lakeid)>0.0 ) Lake_inflow(lakeid) = Lake_inflow(lakeid) + Lake_gain(lakeid)
        ENDIF
        Lake_outflow(lakeid) = Lake_evap(lakeid)
        IF ( Lake_transfer_water_use==ACTIVE ) THEN
          IF ( Lake_transfer(lakeid)>0.0 ) Lake_outflow(lakeid) = Lake_outflow(lakeid) - Lake_transfer(lakeid)
        ENDIF
        IF ( Weir_gate_flag==ACTIVE ) THEN
          tocfs = Lake_area(lakeid)*Cfs_conv
          Lake_seep_in(lakeid) = tocfs*Gw_seep_lakein(lakeid)
          Lake_2gw(lakeid) = tocfs*Lake_seepage(lakeid)
          Lake_inflow(lakeid) = Lake_inflow(lakeid) + Lake_seep_in(lakeid)
          Lake_outflow(lakeid) = Lake_outflow(lakeid) + Lake_2gw(lakeid)
        ENDIF
      ENDDO

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
            IF ( Segment_type(iorder)==LAKE ) THEN ! TS must equal 24
              lakeid = Lake_segment_id(iorder)
              lake_in_ts = Lake_inflow(lakeid) + currin
! what about water use?
              CALL route_lake(lakeid, Lake_type(lakeid), Lake_area(lakeid), lake_in_ts, &
     &                        Lake_outvol_ts(lakeid), Lake_vol(lakeid))

              IF ( Lake_transfer_water_use==ACTIVE ) THEN
                IF ( Lake_transfer(lakeid)>0.0 ) THEN
                  IF ( Lake_vol(lakeid)<0.0 ) THEN
                    PRINT *, 'ERROR, not enough storage for transfer from Lake storage:', &
     &                       lakeid, ' Date:', Nowyear, Nowmonth, Nowday
                    PRINT *, '       storage: ', Lake_vol(lakeid), '; transfer: ', Lake_transfer(lakeid)/Cfs_conv
                    ERROR STOP ERROR_water_use
                  ENDIF
                ENDIF
              ENDIF

              Outflow_ts(iorder) = Lake_outcfs(lakeid)
              Seg_outflow(iorder) = Lake_outcfs(lakeid)
              Lake_stream_in(lakeid) = Seg_upstream_inflow(iorder)
            ELSE
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
            IF ( Use_transfer_segment==ACTIVE ) THEN
              PRINT *, 'ERROR, transfer(s) from stream segment:', iorder, ' causes outflow to be negative'
              PRINT *, '       outflow =', Outflow_ts(iorder), ' must fix water-use stream segment transfer file'
            ELSE
              PRINT *, 'ERROR, outflow from segment:', iorder, ' is negative:', Outflow_ts(iorder)
              PRINT *, '       routing parameters may be invalid'
            ENDIF
            CALL error_stop('negative streamflow in muskingum_lake', ERROR_streamflow)
          ENDIF

          IF ( Segment_type(iorder)/=2 ) THEN
            ! Seg_outflow (the mean daily flow rate for each segment) will be the average of the hourly values.
            Seg_outflow(iorder) = Seg_outflow(iorder) + Outflow_ts(iorder)
            ! pastout is equal to the Inflow_ts on the previous routed timestep
            Pastout(iorder) = Outflow_ts(iorder)
          ENDIF

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
      ! add water balance check
      DO i = 1, Nsegment
        segtype = Segment_type(i)
        IF ( segtype/=2 ) Seg_outflow(i) = Seg_outflow(i) * ONE_24TH ! lake values set above
        segout = Seg_outflow(i)
        Seg_inflow(i) = Seg_inflow(i) * ONE_24TH
        Seg_upstream_inflow(i) = Currinsum(i) * ONE_24TH
! Flow_out is the total flow out of the basin, which allows for multiple outlets
! includes closed basins (tosegment=0)
        IF ( segtype==1 ) THEN
          Flow_headwater = Flow_headwater + segout
        ELSEIF ( segtype==2 ) THEN
          Flow_to_lakes = Flow_to_lakes + segout
          lakeid = Lake_segment_id(i)
          Lake_outcms(lakeid) = Lake_outcfs(lakeid)*CFS2CMS_CONV
          Basin_lake_stor = Basin_lake_stor + Lake_vol(Lakeid)*12.0D0
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
      Basin_2ndstflow = Basin_2ndstflow*Basin_area_inv
      Basin_lake_stor = Basin_lake_stor*Basin_area_inv

!      write(77,'(10f11.3)') lake_vol, lake_outflow, lake_stream_in, lake_outcfs, elevlake
      END FUNCTION muskingum_lake_run

!     ***********************************
!     * Route Lake
!     ***********************************
      SUBROUTINE route_lake(Lakeid, Laketype, Lake_area, Lake_in_ts, Lake_outvol_ts, Lake_vol)
      USE PRMS_CONSTANTS, ONLY: ERROR_lake, ACTIVE, DEBUG_less
      USE PRMS_MODULE, ONLY: Nratetbl, Print_debug
      USE PRMS_MUSKINGUM_LAKE, ONLY: Lake_sto, Lake_outcfs, Lake_inflow, Din1, Lake_outflow, &
     &    Wvd, Nsos, S24, C24, Lake_coef, Lake_invol, Weir_coef, Weir_len, Elev_outflow, &
     &    Tbl_gate, Tbl_stage, Rate_table, Tbl_gate2, Tbl_stage2, Rate_table2, &
     &    Tbl_gate3, Tbl_stage3, Rate_table3, Tbl_gate4, Tbl_stage4, Rate_table4, &
     &    Obsout_lake, Ratetbl_lake, Ngate, Nstage, Ngate2, Nstage2, Ngate3, Nstage3, Ngate4, Nstage4, &
     &    Secondoutflow_flag, Lake_out2, Lake_out2_a, Lake_out2_b, Lake_outq2, Basin_2ndstflow
      USE PRMS_OBS, ONLY: Gate_ht, Streamflow_cfs
      USE PRMS_ROUTING, ONLY: Cfs2acft
      USE PRMS_GWFLOW, ONLY: Elevlake
      IMPLICIT NONE
! Functions
      INTRINSIC :: EXP, DBLE, SNGL, DABS
      EXTERNAL :: table_comp
! Arguments
      INTEGER, INTENT(IN) :: Lakeid, Laketype
      DOUBLE PRECISION, INTENT(IN) :: Lake_area, Lake_in_ts
      DOUBLE PRECISION, INTENT(INOUT) :: Lake_outvol_ts, Lake_vol ! acft
! Local Variables
      INTEGER :: n, jjj, i
      REAL :: q1, q3, elevold, head, new_elevlake, head2, scnd_cfs1, scnd_cfs2
      DOUBLE PRECISION :: avin, s2o2, q2, lake_out, diff_vol, lake_out1
      DOUBLE PRECISION :: xkt, coef2, lake_storage
!***********************************************************************
      !!! ? adjust storage at end of time step ??
      ! q2 = lake out in cfs
      q2 = 0.0D0
      lake_storage = 0.0D0
      IF ( Laketype==1 .OR. Laketype==2 ) THEN
        lake_storage = Lake_sto(Lakeid)
      ELSE
        Lake_outcfs(Lakeid) = 0.0D0
      ENDIF
!   Compute outflow using Puls routing method
      IF ( Laketype==1 ) THEN
        !rsr, why half of current in and last in???
        avin = (Lake_inflow(Lakeid)+Din1(Lakeid))*0.5D0
        s2o2 = lake_storage - (Lake_outflow(Lakeid)+Lake_outcfs(Lakeid))*0.5D0
        s2o2 = s2o2 + avin
        Din1(Lakeid) = Lake_inflow(Lakeid)
        n = Nsos(Lakeid)
        DO jjj = 2, n
          IF ( s2o2<Wvd(jjj, Lakeid) ) THEN
            n = jjj
            EXIT
          ENDIF
        ENDDO
        q2 = S24(n, Lakeid)*s2o2 + C24(n, Lakeid)
        Lake_outcfs(Lakeid) = 0.0D0 ! set to q2 below
        lake_storage = s2o2 - q2*0.5D0
        Lake_sto(Lakeid) = lake_storage

!   Compute outflow using linear reservoir method, ?? need to fix for hourly loop
      ELSEIF ( Laketype==2 ) THEN
        !rsr, why half of current in and last in???
        avin = (Lake_inflow(Lakeid)+Din1(Lakeid)-Lake_outflow(Lakeid))*0.5D0
        xkt = DBLE( Lake_coef(Lakeid) )
        coef2 = 1.0D0 - EXP(-xkt)
        q2 = (avin*(1.0D0-(coef2/xkt))) + lake_storage*coef2
        Din1(Lakeid) = Lake_inflow(Lakeid)
        lake_storage = lake_storage + avin - q2
        Lake_sto(Lakeid) = lake_storage

!   Compute using flow through reservoir
      ELSEIF ( Laketype==3 ) THEN
        q2 = Lake_inflow(Lakeid) - Lake_outflow(Lakeid)

!   Set outflow as a measured flow
      ELSEIF ( Laketype==6 ) THEN
        q2 = Streamflow_cfs(Obsout_lake(Lakeid))
        IF ( q2<0.0D0 ) THEN
          IF ( Print_debug>DEBUG_less ) THEN
            PRINT *, 'WARNING: specified observed runoff value for outflow from lake < 0:', Lakeid, ' value:', q2
            PRINT *, 'runoff id:', Obsout_lake(Lakeid), ' outflow set to 0.0'
          ENDIF
          q2 = 0.0D0
        ENDIF

      ELSE ! 4 or 5; broad-crested weir or gate opening
        elevold = Elevlake(Lakeid)

        ! units lake_invol = acft
        Lake_invol(Lakeid) = Lake_in_ts * Cfs2acft

        ! units lake_out = acft
        lake_out = Lake_outflow(Lakeid) * Cfs2acft ! evap and seepage
        diff_vol = Lake_invol(Lakeid) - lake_out
        q1 = 0.0
        q3 = 0.0

!   Compute using lake surface elevation and broad crested weir
        IF ( Laketype==4 ) THEN
          head = elevold - Elev_outflow(Lakeid)
          IF ( head<0.0 ) head = 0.0
          q1 = (head**1.5) * Weir_coef(Lakeid) * Weir_len(Lakeid)
          lake_out1 = DBLE(q1)*Cfs2acft

          ! new_elevlake has units of feet
          new_elevlake = elevold + SNGL( (diff_vol-lake_out1)/Lake_area )

          head2 = (new_elevlake+elevold)*0.5 - Elev_outflow(Lakeid)
          IF ( head2<0.0 ) head2 = 0.0
          q3 = (head2**1.5) * Weir_coef(Lakeid) * Weir_len(Lakeid)

!  Compute using a rating table of lake surface elevation & gate opening
        ELSE ! type = 5
          DO i = 1, Nratetbl
            IF ( Lakeid==Ratetbl_lake(i) ) THEN
              IF ( i==1 ) THEN
                CALL table_comp(Ngate, Nstage, Tbl_gate, Tbl_stage, &
     &                          Rate_table, elevold, Gate_ht(i), q1, Lake_area)
              ELSEIF ( i==2 ) THEN
                CALL table_comp(Ngate2, Nstage2, Tbl_gate2, Tbl_stage2, &
     &                          Rate_table2, elevold, Gate_ht(i), q1, Lake_area)
              ELSEIF ( i==3 ) THEN
                CALL table_comp(Ngate3, Nstage3, Tbl_gate3, Tbl_stage3, &
     &                          Rate_table3, elevold, Gate_ht(i), q1, Lake_area)
              ELSEIF ( i==4 ) THEN
                CALL table_comp(Ngate4, Nstage4, Tbl_gate4, Tbl_stage4, &
     &                          Rate_table4, elevold, Gate_ht(i), q1, Lake_area)
              ENDIF
            ENDIF
          ENDDO
          scnd_cfs1 = 0.0D0
          IF ( Secondoutflow_flag==ACTIVE ) THEN
!  if lake has a second outlet then outflow in cfs is computed by
!           Q = Lake_out2_a*Elevlake - Lake_out2_b
!               (as per Rob Dudley email 7 Sep 2006)
            IF ( Lake_out2(Lakeid)==1 ) scnd_cfs1 = (Lake_out2_a(Lakeid)*elevold) - Lake_out2_b(Lakeid)
          ENDIF

          lake_out1 = DBLE(q1+scnd_cfs1)*Cfs2acft

          ! new_elevlake has units of feet
          new_elevlake = elevold + SNGL( (diff_vol-lake_out1)/Lake_area )

          DO i = 1, Nratetbl
            IF ( Lakeid==Ratetbl_lake(i) ) THEN
              IF ( i==1 ) THEN
                CALL table_comp(Ngate, Nstage, Tbl_gate, Tbl_stage, Rate_table, &
     &                          new_elevlake, Gate_ht(i), q3, Lake_area)
              ELSEIF ( i==2 ) THEN
                CALL table_comp(Ngate2, Nstage2, Tbl_gate2, Tbl_stage2, Rate_table2, &
     &                          new_elevlake, Gate_ht(i), q3, Lake_area)
              ELSEIF ( i==3 ) THEN
                CALL table_comp(Ngate3, Nstage3, Tbl_gate3, Tbl_stage3, Rate_table3, &
     &                          new_elevlake, Gate_ht(i), q3, Lake_area)
              ELSEIF ( i==4 ) THEN
                CALL table_comp(Ngate4, Nstage4, Tbl_gate4, Tbl_stage4, Rate_table4, &
     &                          new_elevlake, Gate_ht(i), q3, Lake_area)
              ENDIF
            ENDIF
          ENDDO

          IF ( Secondoutflow_flag==ACTIVE ) THEN
            IF ( Lake_out2(lakeid)==1 ) THEN
              scnd_cfs2 = (Lake_out2_a(Lakeid)*new_elevlake) - Lake_out2_b(Lakeid)
            ELSE
              scnd_cfs2 = 0.0D0
            ENDIF
            Lake_outq2(Lakeid) = (scnd_cfs1+scnd_cfs2)*0.5D0
            Basin_2ndstflow = Basin_2ndstflow + Lake_outq2(Lakeid)*Cfs2acft*12.0D0
          ENDIF
        ENDIF

        q2 = DBLE( (q1+q3)*0.5 )
!       !sanity check, rsr
        IF ( q2<0.0D0 ) PRINT *, 'q2<0', q2, ' lake:', Lakeid
        IF ( Secondoutflow_flag==ACTIVE ) q2 = q2 + Lake_outq2(Lakeid)

        Lake_outvol_ts = q2*Cfs2acft + lake_out
        Lake_vol = Lake_vol + Lake_invol(Lakeid) - Lake_outvol_ts
        IF ( Lake_vol<0.0D0 ) THEN
          Lake_outvol_ts = DABS(Lake_vol)
          Lake_vol = 0.0D0
        ELSE
        ! adjust lake elevation with stream and lateral inflows
        ! and streamflow, any second outlet, GWR, and evaporation outflows
          Elevlake(Lakeid) = Elevlake(Lakeid) + SNGL( (Lake_invol(Lakeid)-Lake_outvol_ts)/Lake_area )
        ENDIF
      ENDIF
      IF ( lake_storage<0.0D0 ) THEN
        PRINT *, 'ERROR: lake storage < 0 lake:', Lakeid, '; storage:', lake_storage
        CALL error_stop('in muskingum_lake', ERROR_lake)
      ENDIF

      Lake_outcfs(Lakeid) = q2

      END SUBROUTINE route_lake

!=====================================================================
!    Rating table computation
!=====================================================================
      SUBROUTINE table_comp(Ngate, Nstage, Tbl_gate, Tbl_stage, Rate_table, Elevlake, Gate_ht, Q2, Lake_area)
      USE PRMS_CONSTANTS, ONLY: DEBUG_less
      USE PRMS_MODULE, ONLY: Print_debug
      USE PRMS_ROUTING, ONLY: Cfs2acft
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Ngate, Nstage
      REAL, INTENT(IN) :: Tbl_gate(Ngate), Tbl_stage(Nstage), Rate_table(Nstage, Ngate), Gate_ht, Elevlake
      DOUBLE PRECISION, INTENT(IN) :: Lake_area
      REAL, INTENT(OUT) :: Q2
! Functions
      INTRINSIC :: SNGL
! Local Variables
      INTEGER m, mm, stg1, stg2, gate1, gate2
      REAL :: diff_q_stg1, diff_q_stg2, ratiog, ratios, q_stg1, q_stg2, diffq
!***********************************************************************
      IF ( Elevlake<Tbl_stage(Nstage) ) THEN
        Q2 = 0.0

      ELSEIF ( Elevlake>Tbl_stage(1) ) THEN
        ! lake elevation is > maximum stage, spill all water
        Q2 = (Elevlake-Tbl_stage(1))*SNGL(Lake_area/Cfs2acft)
        IF ( Print_debug>DEBUG_less ) THEN
          PRINT *, 'WARNING, lake elevation > maximum stage in rating table all water above rating table spills'
          PRINT *, 'Lake elevation:', Elevlake, ' Rating table stage:', Tbl_stage(1), ' discharge to stream:', Q2
        ENDIF
      ELSE
        stg2 = 1
        stg1 = 0
        DO m = 1, Nstage
          IF ( Elevlake>Tbl_stage(m) ) THEN
            IF ( m==1 ) THEN
              stg2 = 1
              stg1 = 1
            ELSE
              stg2 = m
              stg1 = m - 1
            ENDIF
            EXIT
          ENDIF
        ENDDO

        gate2 = Ngate
        gate1 = Ngate - 1
        IF ( Gate_ht<=Tbl_gate(Ngate) ) THEN
          DO mm = 1, Ngate
            IF ( Tbl_gate(mm)>Gate_ht ) THEN
              IF ( mm==1 ) THEN
                gate2 = 1
                gate1 = 1
              ELSE
                gate2 = mm
                gate1 = mm - 1
              ENDIF
              EXIT
            ENDIF
          ENDDO
        ENDIF

        IF ( stg1==0 ) THEN
          Q2 = Rate_table(1, gate2)

        ELSE
          diff_q_stg2 = Rate_table(stg2, gate2) - Rate_table(stg2, gate1)
          diff_q_stg1 = Rate_table(stg1, gate2) - Rate_table(stg1, gate1)

 !rsr, possible divide by 0.0???
          ratiog = (Gate_ht-Tbl_gate(gate1))/(Tbl_gate(gate2)-Tbl_gate(gate1))
          q_stg2 = (ratiog*diff_q_stg2) + Rate_table(stg2, gate1)
          q_stg1 = (ratiog*diff_q_stg1) + Rate_table(stg1, gate1)

 !rsr, possible divide by 0.0???
          ratios = (Elevlake-Tbl_stage(stg2))/(Tbl_stage(stg1)-Tbl_stage(stg2))
          diffq = q_stg1 - q_stg2
          Q2 = q_stg2 + (ratios*diffq)
        ENDIF
      ENDIF

      END SUBROUTINE table_comp

!***********************************************************************
!     muskingum_lake_restart - write or read restart file
!***********************************************************************
      SUBROUTINE muskingum_lake_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: ACTIVE, SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_BASIN, ONLY: Puls_lin_flag
      USE PRMS_MUSKINGUM_LAKE
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
        IF ( Puls_lin_flag==ACTIVE ) THEN
          WRITE ( Restart_outunit ) Din1
          WRITE ( Restart_outunit ) Lake_sto
        ENDIF
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Outflow_ts
        IF ( Puls_lin_flag==ACTIVE ) THEN
          READ ( Restart_inunit ) Din1
          READ ( Restart_inunit ) Lake_sto
        ENDIF
      ENDIF
      END SUBROUTINE muskingum_lake_restart
