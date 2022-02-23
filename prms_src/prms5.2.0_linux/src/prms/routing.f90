!***********************************************************************
! Defines stream and lake routing parameters and variables
!***********************************************************************
      MODULE PRMS_ROUTING
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ACTIVE, OFF, FT2_PER_ACRE, &
     &    NEARZERO, DNEARZERO, OUTFLOW_SEGMENT, ERROR_param, &
     &    strmflow_muskingum_mann_module, strmflow_muskingum_lake_module, &
     &    strmflow_muskingum_module, strmflow_in_out_module, CASCADE_OFF, CASCADE_HRU_SEGMENT
      USE PRMS_MODULE, ONLY: Nhru, Nsegment, Model, Init_vars_from_file, &
     &    Strmflow_flag, Cascade_flag, Print_debug, Glacier_flag, &
     &    Water_use_flag, Segment_transferON_OFF, Inputerror_flag, Parameter_check_flag
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Streamflow Routing Init'
      character(len=7), parameter :: MODNAME = 'routing'
      character(len=*), parameter :: Version_routing = '2020-12-02'
      DOUBLE PRECISION, SAVE :: Cfs2acft
      DOUBLE PRECISION, SAVE :: Segment_area
      INTEGER, SAVE :: Use_transfer_segment, Noarea_flag, Hru_seg_cascades
      INTEGER, SAVE, ALLOCATABLE :: Segment_order(:), Segment_up(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Segment_hruarea(:)
      !CHARACTER(LEN=32), SAVE :: Outfmt
      INTEGER, SAVE, ALLOCATABLE :: Ts_i(:)
      REAL, SAVE, ALLOCATABLE :: Ts(:), C0(:), C1(:), C2(:)
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_segment_storage
      DOUBLE PRECISION, SAVE :: Flow_to_lakes, Flow_to_ocean, Flow_to_great_lakes, Flow_out_region
      DOUBLE PRECISION, SAVE :: Flow_in_region, Flow_in_nation, Flow_headwater, Flow_out_NHM
      DOUBLE PRECISION, SAVE :: Flow_in_great_lakes, Flow_replacement, Flow_terminus
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Seginc_ssflow(:), Seginc_sroff(:), Segment_delta_flow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Seginc_gwflow(:), Seginc_swrad(:), Seginc_potet(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_outflow(:), Seg_ssflow(:), Seg_sroff(:), Seg_gwflow(:)
!   Declared Parameters
      INTEGER, SAVE, ALLOCATABLE :: Segment_type(:), Tosegment(:), Hru_segment(:), Obsin_segment(:), Obsout_segment(:)
      REAL, SAVE, ALLOCATABLE :: K_coef(:), X_coef(:)
      REAL, SAVE, ALLOCATABLE :: Seg_depth(:), Mann_n(:), Segment_flow_init(:)
      REAL, SAVE, ALLOCATABLE :: Seg_length(:), Seg_slope(:)
      END MODULE PRMS_ROUTING

!***********************************************************************
!     Main routing routine
!***********************************************************************
      INTEGER FUNCTION routing()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Init_vars_from_file, Save_vars_to_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: routingdecl, routinginit, route_run
      EXTERNAL :: routing_restart
!***********************************************************************
      routing = 0

      IF ( Process_flag==RUN ) THEN
        routing = route_run()
      ELSEIF ( Process_flag==DECL ) THEN
        routing = routingdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL routing_restart(READ_INIT)
        routing = routinginit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL routing_restart(SAVE_INIT)
      ENDIF

      END FUNCTION routing

!***********************************************************************
!     routingdecl - set up parameters
!***********************************************************************
      INTEGER FUNCTION routingdecl()
      USE PRMS_ROUTING
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      routingdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_routing)

! Declared Variables
      ALLOCATE ( Hru_outflow(Nhru) )
      IF ( declvar(MODNAME, 'hru_outflow', 'nhru', Nhru, 'double', &
     &     'Total flow leaving each HRU', &
     &     'cfs', Hru_outflow)/=0 ) CALL read_error(3, 'hru_outflow')

      IF ( declvar(MODNAME, 'flow_to_lakes', 'one', 1, 'double', &
     &     'Total flow to lakes (segment_type=2)', &
     &     'cfs', Flow_to_lakes)/=0 ) CALL read_error(3, 'flow_to_lakes')

      IF ( declvar(MODNAME, 'flow_terminus', 'one', 1, 'double', &
     &     'Total flow to terminus segments (segment_type=9)', &
     &     'cfs', Flow_terminus)/=0 ) CALL read_error(3, 'flow_terminus')

      IF ( declvar(MODNAME, 'flow_to_ocean', 'one', 1, 'double', &
     &     'Total flow to oceans (segment_type=8)', &
     &     'cfs', Flow_to_ocean)/=0 ) CALL read_error(3, 'flow_to_ocean')

      IF ( declvar(MODNAME, 'flow_to_great_lakes', 'one', 1, 'double', &
     &     'Total flow to Great Lakes (segment_type=11)', &
     &     'cfs', Flow_to_great_lakes)/=0 ) CALL read_error(3, 'Flow_to_great_lakes')

      IF ( declvar(MODNAME, 'flow_out_region', 'one', 1, 'double', &
     &     'Total flow out of region (segment_type=7)', &
     &     'cfs', Flow_out_region)/=0 ) CALL read_error(3, 'flow_out_region')

      IF ( declvar(MODNAME, 'flow_out_NHM', 'one', 1, 'double', &
     &     'Total flow out of model domain to Mexico or Canada (segment_type=5)', &
     &     'cfs', Flow_out_NHM)/=0 ) CALL read_error(3, 'flow_out_NHM')

      IF ( declvar(MODNAME, 'flow_in_region', 'one', 1, 'double', &
     &     'Total flow into region (segment_type=6)', &
     &     'cfs', Flow_in_region)/=0 ) CALL read_error(3, 'flow_in_region')

      IF ( declvar(MODNAME, 'flow_in_nation', 'one', 1, 'double', &
     &     'Total flow into model domain from Mexico or Canada (segment_type=4)', &
     &     'cfs', Flow_in_nation)/=0 ) CALL read_error(3, 'flow_in_nation')

      IF ( declvar(MODNAME, 'flow_headwater', 'one', 1, 'double', &
     &     'Total flow out of headwater segments (segment_type=1)', &
     &     'cfs', Flow_headwater)/=0 ) CALL read_error(3, 'flow_headwater')

      IF ( declvar(MODNAME, 'flow_in_great_lakes', 'one', 1, 'double', &
     &     'Total flow into model domain from Great Lakes (segment_type=10)', &
     &     'cfs', Flow_in_great_lakes)/=0 ) CALL read_error(3, 'flow_in_great_lakes')

      IF ( declvar(MODNAME, 'flow_replacement', 'one', 1, 'double', &
     &     'Total flow out from replacement flow (segment_type=3)', &
     &     'cfs', Flow_replacement)/=0 ) CALL read_error(3, 'flow_replacement')

      ! 0 = normal; 1 = headwater; 2 = lake; 3 = replacement flow; 4 = inbound to nation;
      ! 5 = outbound from nation; 6 = inbound to region; 7 = outbound from region;
      ! 8 = drains to ocean; 9 = sink (terminus to soil); 10 = inbound from Great Lakes;
      ! 11 = outbound to Great Lakes; 12 = ephemeral; + 100 user updated; 1000 user virtual segment
      ! 100 = user normal; 101 - 108 = not used; 109 sink (tosegment used by Lumen)

      IF ( Strmflow_flag==strmflow_muskingum_mann_module .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Mann_n(Nsegment) )
        IF ( declparam( MODNAME, 'mann_n', 'nsegment', 'real', &
     &       '0.04', '0.001', '0.15', &
     &       'Mannings roughness coefficient', &
     &       'Mannings roughness coefficient for each segment', &
     &       'dimensionless')/=0 ) CALL read_error(1, 'mann_n')

        ALLOCATE ( Seg_length(Nsegment) )
        IF ( declparam( MODNAME, 'seg_length', 'nsegment', 'real', &
     &       '1000.0', '1.0', '100000.0', &
     &       'Length of each segment', &
     &       'Length of each segment', &
     &       'meters')/=0 ) CALL read_error(1, 'seg_length')

        ALLOCATE ( Seg_depth(Nsegment) )
        IF ( declparam(MODNAME, 'seg_depth', 'nsegment', 'real', &
     &       '1.0', '0.03', '250.0', &
     &       'Segment river depth', &
     &       'Segment river depth at bankfull; shallowest depth from Blackburn-Lynch (2017);'//&
     &       'Congo is deepest at 250 m but in the US it is probably the Hudson at 66 m', &
     &       'meters')/=0 ) CALL read_error(1, 'seg_depth')
      ENDIF

      IF ( Strmflow_flag==strmflow_muskingum_mann_module .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Seg_slope(Nsegment) )
        IF ( declparam( MODNAME, 'seg_slope', 'nsegment', 'real', &
     &       '0.0001', '0.0000001', '2.0', &
     &       'Surface slope of each segment', &
     &       'Surface slope of each segment as approximation for bed slope of stream', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'seg_slope')
      ENDIF

      ALLOCATE ( Segment_type(Nsegment) )
      IF ( declparam(MODNAME, 'segment_type', 'nsegment', 'integer', &
     &     '0', '0', '111', &
     &     'Segment type', &
     &     'Segment type (0=segment; 1=headwater; 2=lake; 3=replace inflow; 4=inbound to NHM;'// &
     &     ' 5=outbound from NHM; 6=inbound to region; 7=outbound from region; 8=drains to ocean;'// &
     &     ' 9=sink; 10=inbound from Great Lakes; 11=outbound to Great Lakes)', &
     &     'none')/=0 ) CALL read_error(1, 'segment_type')

      ! user updated values if different than tosegment_orig
      ! -5 = outbound from NHM; -6 = inbound from region; -7 = outbound from region;
      ! -8 = drains to ocean; -11 = drains to Great Lake
      ALLOCATE ( Tosegment(Nsegment) )
      IF ( declparam(MODNAME, 'tosegment', 'nsegment', 'integer', &
     &     '0', '0', '9999999', &
     &     'The index of the downstream segment', &
     &     'Index of downstream segment to which the segment'// &
     &     ' streamflow flows, for segments that do not flow to another segment enter 0', &
     &     'none')/=0 ) CALL read_error(1, 'tosegment')

      IF ( Cascade_flag==CASCADE_OFF .OR. Cascade_flag==CASCADE_HRU_SEGMENT .OR. Model==DOCUMENTATION ) THEN
        Hru_seg_cascades = ACTIVE
        ALLOCATE ( Hru_segment(Nhru) )
        IF ( declparam(MODNAME, 'hru_segment', 'nhru', 'integer', &
     &       '0', 'bounded', 'nsegment', &
     &       'Segment index for HRU lateral inflows', &
     &       'Segment index to which an HRU contributes lateral flows'// &
     &       ' (surface runoff, interflow, and groundwater discharge)', &
     &       'none')/=0 ) CALL read_error(1, 'hru_segment')
      ELSE
        Hru_seg_cascades = OFF
      ENDIF

      ALLOCATE ( Obsin_segment(Nsegment) )
      IF ( declparam(MODNAME, 'obsin_segment', 'nsegment', 'integer', &
     &     '0', 'bounded', 'nobs', &
     &     'Index of measured streamflow station that replaces inflow to a segment', &
     &     'Index of measured streamflow station that replaces inflow to a segment', &
     &     'none')/=0 ) CALL read_error(1, 'obsin_segment')

      ALLOCATE ( Obsout_segment(Nsegment) )
      IF ( declparam(MODNAME, 'obsout_segment', 'nsegment', 'integer', &
     &     '0', 'bounded', 'nobs', &
     &     'Index of measured streamflow station that replaces outflow from a segment', &
     &     'Index of measured streamflow station that replaces outflow from a segment', &
     &     'none')/=0 ) CALL read_error(1, 'obsout_segment')

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 ) THEN
        ALLOCATE ( Segment_flow_init(Nsegment) )
        IF ( declparam(MODNAME, 'segment_flow_init', 'nsegment', 'real', &
     &       '0.0', '0.0', '1.0E7', &
     &       'Initial flow in each stream segment', &
     &       'Initial flow in each stream segment', &
     &       'cfs')/=0 ) CALL read_error(1, 'segment_flow_init')
      ENDIF

      IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &     Strmflow_flag==strmflow_muskingum_mann_module ) ALLOCATE ( K_coef(Nsegment) )
      IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &     Model==DOCUMENTATION ) THEN
        IF ( declparam(MODNAME, 'K_coef', 'nsegment', 'real', &
     &       '1.0', '0.01', '24.0', &
     &       'Muskingum storage coefficient', &
     &       'Travel time of flood wave from one segment to the next downstream segment,'// &
     &       ' called the Muskingum storage coefficient; enter 1.0 for reservoirs,'// &
     &       ' diversions, and segment(s) flowing out of the basin', &
     &       'hours')/=0 ) CALL read_error(1, 'K_coef')
      ENDIF

      IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &     Strmflow_flag==strmflow_muskingum_mann_module .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( X_coef(Nsegment) )
        IF ( declparam(MODNAME, 'x_coef', 'nsegment', 'real', &
     &       '0.2', '0.0', '0.5', &
     &       'Routing weighting factor', &
     &       'The amount of attenuation of the flow wave, called the'// &
     &       ' Muskingum routing weighting factor; enter 0.0 for'// &
     &       ' reservoirs, diversions, and segment(s) flowing out of the basin', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'x_coef')
        IF ( declvar(MODNAME, 'basin_segment_storage', 'one', 1, 'double', &
     &       'Basin area-weighted average storage in the stream network', &
     &       'inches', Basin_segment_storage)/=0 ) CALL read_error(3, 'basin_segment_storage')
        ALLOCATE ( Segment_delta_flow(Nsegment) )
        IF ( declvar(MODNAME, 'segment_delta_flow', 'nsegment', Nsegment, 'double', &
     &       'Cumulative flow in minus flow out for each stream segment', &
     &       'cfs', Segment_delta_flow)/=0 ) CALL read_error(3, 'segment_delta_flow')
      ENDIF

      IF ( Hru_seg_cascades==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Seginc_potet(Nsegment) )
        IF ( declvar(MODNAME, 'seginc_potet', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average potential ET for each segment'// &
     &       ' from HRUs contributing flow to the segment', &
     &       'inches', Seginc_potet)/=0 ) CALL read_error(3, 'seginc_potet')

        ALLOCATE ( Seginc_swrad(Nsegment) )
        IF ( declvar(MODNAME, 'seginc_swrad', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average solar radiation for each segment'// &
     &       ' from HRUs contributing flow to the segment', &
     &       'Langleys', Seginc_swrad)/=0 ) CALL read_error(3, 'seginc_swrad')

        ALLOCATE ( Seginc_ssflow(Nsegment) )
        IF ( declvar(MODNAME, 'seginc_ssflow', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average interflow for each segment from'// &
     &       ' HRUs contributing flow to the segment', &
     &       'cfs', Seginc_ssflow)/=0 ) CALL read_error(3, 'seginc_ssflow')

        ALLOCATE ( Seginc_gwflow(Nsegment) )
        IF ( declvar(MODNAME, 'seginc_gwflow', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average groundwater discharge for each'// &
     &       ' segment from HRUs contributing flow to the segment', &
     &       'cfs', Seginc_gwflow)/=0 ) CALL read_error(3, 'seginc_gwflow')

        ALLOCATE ( Seginc_sroff(Nsegment) )
        IF ( declvar(MODNAME, 'seginc_sroff', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average surface runoff for each'// &
     &       ' segment from HRUs contributing flow to the segment', &
     &       'cfs', Seginc_sroff)/=0 ) CALL read_error(3, 'seginc_sroff')

        ALLOCATE ( Seg_ssflow(Nsegment) )
        IF ( declvar(MODNAME, 'seg_ssflow', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average interflow for each segment from'// &
     &       ' HRUs contributing flow to the segment and upstream HRUs', &
     &       'inches', Seg_ssflow)/=0 ) CALL read_error(3, 'seg_ssflow')

        ALLOCATE ( Seg_gwflow(Nsegment) )
        IF ( declvar(MODNAME, 'seg_gwflow', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average groundwater discharge for each segment from'// &
     &       ' HRUs contributing flow to the segment and upstream HRUs', &
     &       'inches', Seg_gwflow)/=0 ) CALL read_error(3, 'seg_gwflow')

        ALLOCATE ( Seg_sroff(Nsegment) )
        IF ( declvar(MODNAME, 'seg_sroff', 'nsegment', Nsegment, 'double', &
     &       'Area-weighted average surface runoff for each segment from'// &
     &       ' HRUs contributing flow to the segment and upstream HRUs', &
     &       'inches', Seg_sroff)/=0 ) CALL read_error(3, 'seg_sroff')
      ENDIF

      ! local arrays
      ALLOCATE ( Segment_order(Nsegment), Segment_up(Nsegment), Segment_hruarea(Nsegment) )

      END FUNCTION routingdecl

!**********************************************************************
!     routinginit - check for validity of parameters
!**********************************************************************
      INTEGER FUNCTION routinginit()
      USE PRMS_ROUTING
      USE PRMS_SET_TIME, ONLY: Timestep_seconds
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_area_dble !, Active_area
      USE PRMS_FLOWVARS, ONLY: Seg_outflow, Seg_inflow
      IMPLICIT NONE
! Functions
      INTRINSIC :: MOD
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: read_error
! Local Variables
      INTEGER :: i, j, test, lval, toseg, iseg, isegerr, ierr, eseg
      REAL :: k, x, d, x_max, velocity
      INTEGER, ALLOCATABLE :: x_off(:)
      CHARACTER(LEN=10) :: buffer
!**********************************************************************
      routinginit = 0

      Use_transfer_segment = OFF
      IF ( Water_use_flag==ACTIVE .AND. Segment_transferON_OFF==ACTIVE ) Use_transfer_segment = ACTIVE

      IF ( Init_vars_from_file==0 ) THEN
        Basin_segment_storage = 0.0D0
        IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &       Strmflow_flag==strmflow_muskingum_mann_module ) Segment_delta_flow = 0.0D0
      ENDIF

      IF ( Hru_seg_cascades==ACTIVE ) THEN
        Seginc_potet = 0.0D0
        Seginc_gwflow = 0.0D0
        Seginc_ssflow = 0.0D0
        Seginc_sroff = 0.0D0
        Seginc_swrad = 0.0D0
        Seg_gwflow = 0.0D0
        Seg_ssflow = 0.0D0
        Seg_sroff = 0.0D0
      ENDIF
      Hru_outflow = 0.0D0
      Flow_to_ocean = 0.0D0
      Flow_to_great_lakes = 0.0D0
      Flow_out_region = 0.0D0
      Flow_out_NHM = 0.0D0
      Flow_terminus = 0.0D0
      Flow_to_lakes = 0.0D0
      Flow_in_nation = 0.0D0
      Flow_in_region = 0.0D0
      Flow_headwater = 0.0D0
      Flow_in_great_lakes = 0.0D0
      Flow_replacement = 0.0D0

      Cfs2acft = Timestep_seconds/FT2_PER_ACRE

      IF ( getparam(MODNAME, 'segment_type', Nsegment, 'integer', Segment_type)/=0 ) CALL read_error(2, 'segment_type')
      DO i = 1, Nsegment
        Segment_type(i) = MOD( Segment_type(i), 100 )
      ENDDO

      IF ( Strmflow_flag==strmflow_muskingum_mann_module ) THEN
        IF ( getparam(MODNAME, 'mann_n', Nsegment, 'real', Mann_n)/=0 ) CALL read_error(2, 'mann_n')
        IF ( getparam( MODNAME, 'seg_length', Nsegment, 'real', Seg_length)/=0 ) CALL read_error(2, 'seg_length')
        IF ( getparam( MODNAME, 'seg_slope', Nsegment, 'real', Seg_slope)/=0 ) CALL read_error(2, 'seg_slope')
! find segments that are too short and print them out as they are found
        ierr = 0
        DO i = 1, Nsegment
           IF ( Seg_length(i)<NEARZERO ) THEN
              PRINT *, 'ERROR, seg_length too small for segment:', i, ', value:', Seg_length(i)
              ierr = 1
           ENDIF
        ENDDO
! exit if there are any segments that are too short
        IF ( ierr==1 ) THEN
           Inputerror_flag = ierr
           RETURN
        ENDIF
        IF ( getparam(MODNAME, 'seg_depth', Nsegment, 'real', seg_depth)/=0 ) CALL read_error(2, 'seg_depth')
      ENDIF

      IF ( getparam(MODNAME, 'tosegment', Nsegment, 'integer', Tosegment)/=0 ) CALL read_error(2, 'tosegment')
      IF ( getparam(MODNAME, 'obsin_segment', Nsegment, 'integer', Obsin_segment)/=0 ) CALL read_error(2, 'obsin_segment')
      IF ( getparam(MODNAME, 'obsout_segment', Nsegment, 'integer', Obsout_segment)/=0 ) CALL read_error(2, 'obsout_segment')

      IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &     Strmflow_flag==strmflow_muskingum_mann_module ) THEN
        IF ( getparam(MODNAME, 'x_coef', Nsegment, 'real', X_coef)/=0 ) CALL read_error(2, 'x_coef')
        ALLOCATE ( C1(Nsegment), C2(Nsegment), C0(Nsegment), Ts(Nsegment), Ts_i(Nsegment) )
        IF ( Init_vars_from_file==0 ) Segment_delta_flow = 0.0D0
        IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module ) THEN
          IF ( getparam(MODNAME, 'K_coef', Nsegment, 'real', K_coef)/=0 ) CALL read_error(2, 'K_coef')
        ENDIF
      ENDIF

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 ) THEN
        IF ( getparam(MODNAME, 'segment_flow_init',  Nsegment, 'real', Segment_flow_init)/=0 ) &
     &       CALL read_error(2,'segment_flow_init')
        DO i = 1, Nsegment
          Seg_outflow(i) = Segment_flow_init(i)
          IF ( tosegment(i)>0 ) Seg_inflow(tosegment(i)) = Seg_outflow(i)
        ENDDO
        DEALLOCATE ( Segment_flow_init )
      ENDIF

! if cascades are active then ignore hru_segment
      Noarea_flag = OFF
      IF ( Hru_seg_cascades==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'hru_segment', Nhru, 'integer', Hru_segment)/=0 ) CALL read_error(2, 'hru_segment')
        Segment_hruarea = 0.0D0
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          iseg = Hru_segment(i)
          IF ( iseg>0 ) Segment_hruarea(iseg) = Segment_hruarea(iseg) + Hru_area_dble(i)
        ENDDO
        Segment_area = 0.0D0
        DO j = 1, Nsegment
          Segment_area = Segment_area + Segment_hruarea(j)
          IF ( Segment_hruarea(j)<DNEARZERO ) THEN
            Noarea_flag = ACTIVE
            IF ( Parameter_check_flag>0 ) THEN
              WRITE ( buffer, '(I10)' ) j
              CALL write_outfile('WARNING, No HRUs are associated with segment:'//buffer)
              IF ( Tosegment(j)==0 ) PRINT *, 'WARNING, No HRUs and tosegment=0 for segment:', j
            ENDIF
          ENDIF
        ENDDO
!        IF ( Active_area/=Segment_area ) PRINT *, 'Not all area in model domain included with segments, basin area =', &
!     &                                            Active_area, ' segment area = ', Segment_area
      ENDIF

      isegerr = 0
      Segment_up = 0
      ! Begin the loops for ordering segments
      DO j = 1, Nsegment
        iseg = Obsin_segment(j)
        toseg = Tosegment(j)
        IF ( toseg==j ) THEN
          PRINT *, 'ERROR, tosegment value (', toseg, ') equals itself for segment:', j
          isegerr = 1
        ELSEIF ( toseg>OUTFLOW_SEGMENT ) THEN
          IF ( Tosegment(toseg)==j ) THEN
            PRINT *, 'ERROR, circle found, segment:', j, ' sends flow to segment:', toseg, ' that sends it flow'
            isegerr = 1
          ELSE
            ! load segment_up with last stream segment that flows into a segment
            Segment_up(toseg) = j
          ENDIF
        ENDIF
      ENDDO

      IF ( Parameter_check_flag>0 ) THEN
        DO i = 1, Nsegment
          IF ( Segment_up(i)==0 .AND. Tosegment(i)==OUTFLOW_SEGMENT ) &
     &         PRINT *, 'WARNING, no other segment flows into segment:',  i, ' and tosegment=0'
        ENDDO
      ENDIF

      IF ( isegerr==1 ) THEN
        Inputerror_flag = 1
        RETURN
      ENDIF

      ! Begin the loops for ordering segments
      ALLOCATE ( x_off(Nsegment) )
      x_off = 0
      Segment_order = 0
      lval = 0
      iseg = 0
      eseg = 0
      DO WHILE ( lval<Nsegment )
        ierr = 1
        DO i = 1, Nsegment
          ! If segment "i" has not been crossed out consider it, else continue
          IF ( x_off(i)==1 ) CYCLE
          iseg = i
          ! Test to see if segment "i" is the tosegment from other segments
          test = 1
          DO j = 1, Nsegment
            IF ( Tosegment(j)==i ) THEN
              ! If segment "i" is a to segment, test to see if the originating
              ! segment has been crossed off the list.  if all have been, then
              ! put the segment in as an ordered segment
              IF ( x_off(j)==0 ) THEN
                test = 0
                eseg = j
                EXIT
              ENDIF
            ENDIF
          ENDDO
          IF ( test==1 ) THEN
            lval = lval + 1
            Segment_order(lval) = i
            x_off(i) = 1
            ierr = 0
          ENDIF
        ENDDO
        IF ( ierr==1 ) THEN
          PRINT *, 'ERROR, circular segments involving', iseg, 'and', eseg
          ERROR STOP ERROR_param
        ENDIF
      ENDDO
!      IF ( Print_debug==20 ) THEN
!        PRINT *, 'Stream Network Routing Order:'
!        PRINT '(10I5)', Segment_order
!        PRINT *, 'tosegment:'
!        PRINT '(10I5)', Tosegment
!      ENDIF
      DEALLOCATE ( x_off )

      IF ( Strmflow_flag==strmflow_in_out_module ) RETURN
!
!      Compute the three constants in the Muskingum routing equation based
!      on the values of K_coef and a routing period of 1 hour. See the notes
!      at the top of this file.
!
      C0 = 0.0
      C1 = 0.0
      C2 = 0.0
!make sure K>0
      Ts = 1.0
      ierr = 0
      DO i = 1, Nsegment
        IF ( Strmflow_flag==strmflow_muskingum_mann_module ) THEN
          velocity = (1./Mann_n(i))*SQRT(Seg_slope(i))*Seg_depth(i)**(2./3.) ! simplify if say width>>depth
          K_coef(i) = Seg_length(i)/(velocity*60.*60.) !want in hours, length should include sloped length
          !K_coef(i) = Seg_length(i)*sqrt(1+ Seg_slope(i)**2)/(velocity*60.*60.) !want in hours
        ENDIF

        IF ( Segment_type(i)==2 .AND. K_coef(i)<24.0 ) THEN
          IF ( Parameter_check_flag>0 ) PRINT *, 'WARNING, K_coef must be specified = 24.0 for lake segments'
          K_coef(i) = 24.0
        ENDIF
        IF ( K_coef(i)<0.01 ) K_coef(i) = 0.01 !make compliant with old version of K_coef
        IF ( K_coef(i)>24.0 ) K_coef(i) = 24.0
        k = K_coef(i)
        x = X_coef(i)

! check the values of k and x to make sure that Muskingum routing is stable

        IF ( k<1.0 ) THEN
!          IF ( Parameter_check_flag>0 ) THEN
!            PRINT '(/,A,I6,A,F6.2,/,9X,A,/)', 'WARNING, segment ', i, ' has K_coef < 1.0,', k, &
!     &              'this may produce unstable results'
!              ierr = 1
!          ENDIF
!          Ts(i) = 0.0 ! not sure why this was set to zero, causes divide by 0 if K_coef < 1, BUG FIX 10/18/2016 RSR
          Ts_i(i) = -1

        ELSEIF ( k<2.0 ) THEN
          Ts(i) = 1.0
          Ts_i(i) = 1

        ELSEIF ( k<3.0 ) THEN
          Ts(i) = 2.0
          Ts_i(i) = 2

        ELSEIF ( k<4.0 ) THEN
          Ts(i) = 3.0
          Ts_i(i) = 3

        ELSEIF ( k<6.0 ) THEN
          Ts(i) = 4.0
          Ts_i(i) = 4

        ELSEIF ( k<8.0 ) THEN
          Ts(i) = 6.0
          Ts_i(i) = 6

        ELSEIF ( k<12.0 ) THEN
          Ts(i) = 8.0
          Ts_i(i) = 8

        ELSEIF ( k<24.0 ) THEN
          Ts(i) = 12.0
          Ts_i(i) = 12

        ELSE
          Ts(i) = 24.0
          Ts_i(i) = 24

        ENDIF

!  x must be <= t/(2K) the C coefficents will be negative. Check for this for all segments
!  with Ts >= minimum Ts (1 hour).
        IF ( Ts(i)>0.1 ) THEN
          x_max = Ts(i) / (2.0 * k)
          IF ( x>x_max ) THEN
            PRINT *, 'ERROR, x_coef value is too large for stable routing for segment:', i, ' x_coef:', x
            PRINT *, '       a maximum value of:', x_max, ' is suggested'
            Inputerror_flag = 1
            CYCLE
          ENDIF
        ENDIF

        d = k - (k * x) + (0.5 * Ts(i))
        IF ( ABS(d)<NEARZERO ) THEN
          IF ( Parameter_check_flag>0 ) PRINT *, 'WARNING, segment ', i, ' computed value d <', NEARZERO, ', set to 0.0001'
          d = 0.0001
        ENDIF
        C0(i) = (-(k * x) + (0.5 * Ts(i))) / d
        C1(i) = ((k * x) + (0.5 * Ts(i))) / d
        C2(i) = (k - (k * x) - (0.5 * Ts(i))) / d

        ! the following code was in the original musroute, but, not in Linsley and others
        ! rsr, 3/1/2016 - having < 0 coefficient can cause negative flows as found by Jacob in GCPO headwater
!  if c2 is <= 0.0 then  short travel time though reach (less daily
!  flows), thus outflow is mainly = inflow w/ small influence of previous
!  inflow. Therefore, keep c0 as is, and lower c1 by c2, set c2=0

!  if c0 is <= 0.0 then long travel time through reach (greater than daily
!  flows), thus mainly dependent on yesterdays flows.  Therefore, keep
!  c2 as is, reduce c1 by c0 and set c0=0
! SHORT travel time
        IF ( C2(i)<0.0 ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT '(/,A)', 'WARNING, c2 < 0, set to 0, c1 set to c1 + c2'
            PRINT *, '        old c2:', C2(i), '; old c1:', C1(i), '; new c1:', C1(i) + C2(i)
            PRINT *, '        K_coef:', K_coef(i), '; x_coef:', x_coef(i)
          ENDIF
          C1(i) = C1(i) + C2(i)
          C2(i) = 0.0
        ENDIF

! LONG travel time
        IF ( C0(i)<0.0 ) THEN
          IF ( Parameter_check_flag>0 ) THEN
            PRINT '(/,A)', 'WARNING, c0 < 0, set to 0, c0 set to c1 + c0'
            PRINT *, '      old c0:', C0(i), 'old c1:', C1(i), 'new c1:', C1(i) + C0(i)
            PRINT *, '        K_coef:', K_coef(i), '; x_coef:', x_coef(i)
          ENDIF
          C1(i) = C1(i) + C0(i)
          C0(i) = 0.0
        ENDIF

      ENDDO
      IF ( ierr==1 ) PRINT '(/,A,/)', '***Recommend that the Muskingum parameters be adjusted in the Parameter File'
      DEALLOCATE ( K_coef, X_coef)

      END FUNCTION routinginit

!***********************************************************************
!     route_run - Computes segment flow states and fluxes
!***********************************************************************
      INTEGER FUNCTION route_run()
      USE PRMS_ROUTING
      USE PRMS_BASIN, ONLY: Hru_area, Hru_route_order, Active_hrus
      USE PRMS_CLIMATEVARS, ONLY: Swrad, Potet
      USE PRMS_SET_TIME, ONLY: Timestep_seconds, Cfs_conv
      USE PRMS_FLOWVARS, ONLY: Ssres_flow, Sroff, Seg_lateral_inflow !, Seg_outflow
      USE PRMS_WATER_USE, ONLY: Segment_transfer, Segment_gain
      USE PRMS_GWFLOW, ONLY: Gwres_flow
      USE PRMS_SRUNOFF, ONLY: Strm_seg_in
      USE PRMS_GLACR, ONLY: Glacr_flow
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
! Local Variables
      INTEGER :: i, j, jj, this_seg
      DOUBLE PRECISION :: tocfs
      LOGICAL :: found
!***********************************************************************
      route_run = 0

      Cfs2acft = Timestep_seconds/FT2_PER_ACRE

! seg variables are not computed if cascades are active as hru_segment is ignored
      IF ( Hru_seg_cascades==ACTIVE ) THEN
        ! add hru_ppt, hru_actet
        Seginc_gwflow = 0.0D0
        Seginc_ssflow = 0.0D0
        Seginc_sroff = 0.0D0
        Seginc_swrad = 0.0D0
        Seginc_potet = 0.0D0
        Seg_gwflow = 0.0D0
        Seg_sroff = 0.0D0
        Seg_ssflow = 0.0D0
      ENDIF
      IF ( Cascade_flag==CASCADE_OFF ) THEN
        Seg_lateral_inflow = 0.0D0
      ELSE ! use strm_seg_in for cascade_flag = 1 (CASCADE_NORMAL) or 2 (CASCADE_HRU_SEGMENT)
        Seg_lateral_inflow = Strm_seg_in
      ENDIF

      DO jj = 1, Active_hrus
        j = Hru_route_order(jj)
        tocfs = DBLE( Hru_area(j) )*Cfs_conv
        Hru_outflow(j) = DBLE( (Sroff(j) + Ssres_flow(j) + Gwres_flow(j)) )*tocfs
        ! Note: glacr_flow (from glacier or snowfield) is added as a gain, outside stream network addition
        ! glacr_flow in inch^3, 1728=12^3
        IF ( Glacier_flag==ACTIVE ) Hru_outflow(j) = Hru_outflow(j) + Glacr_flow(j)/1728.0/Timestep_seconds
        IF ( Hru_seg_cascades==ACTIVE ) THEN
          i = Hru_segment(j)
          IF ( i>0 ) THEN
            Seg_gwflow(i) = Seg_gwflow(i) + Gwres_flow(j)
            Seg_sroff(i) = Seg_sroff(i) + Sroff(j)
            Seg_ssflow(i) = Seg_ssflow(i) + Ssres_flow(j)
            ! if cascade_flag = CASCADE_HRU_SEGMENT, seg_lateral_inflow set with strm_seg_in
            IF ( Cascade_flag==CASCADE_OFF ) Seg_lateral_inflow(i) = Seg_lateral_inflow(i) + Hru_outflow(j)
            Seginc_sroff(i) = Seginc_sroff(i) + DBLE( Sroff(j) )*tocfs
            Seginc_ssflow(i) = Seginc_ssflow(i) + DBLE( Ssres_flow(j) )*tocfs
            Seginc_gwflow(i) = Seginc_gwflow(i) + DBLE( Gwres_flow(j) )*tocfs
            Seginc_swrad(i) = Seginc_swrad(i) + DBLE( Swrad(j)*Hru_area(j) )
            Seginc_potet(i) = Seginc_potet(i) + DBLE( Potet(j)*Hru_area(j) )
          ENDIF
        ENDIF
      ENDDO

      IF ( Use_transfer_segment==ACTIVE ) THEN
        DO i = 1, Nsegment
          Seg_lateral_inflow(i) = Seg_lateral_inflow(i) + DBLE( Segment_gain(i) - Segment_transfer(i) )
        ENDDO
      ENDIF

      IF ( Cascade_flag>CASCADE_OFF ) RETURN

! Divide solar radiation and PET by sum of HRU area to get avarage
      IF ( Noarea_flag==OFF ) THEN
        DO i = 1, Nsegment
          Seginc_swrad(i) = Seginc_swrad(i)/Segment_hruarea(i)
          Seginc_potet(i) = Seginc_potet(i)/Segment_hruarea(i)
        ENDDO

! If there are no HRUs associated with a segment, then figure out some
! other way to get the solar radiation, the following is not great
      ELSE !     IF ( Noarea_flag==ACTIVE ) THEN
        DO i = 1, Nsegment
! This reworked by markstrom
          IF ( Segment_hruarea(i)>NEARZERO ) THEN
            Seginc_swrad(i) = Seginc_swrad(i)/Segment_hruarea(i)
            Seginc_potet(i) = Seginc_potet(i)/Segment_hruarea(i)
          ELSE

! Segment does not have any HRUs, check upstream segments.
            this_seg = i
            found = .false.
            do
              if (Segment_hruarea(this_seg) <= NEARZERO) then

                 ! Hit the headwater segment without finding any HRUs (i.e. sources of streamflow)
                 if (segment_up(this_seg) == 0) then
                     found = .false.
                     exit
                 endif

                 ! There is an upstream segment, check that segment for HRUs
                 this_seg = segment_up(this_seg)
              else
                  ! This segment has HRUs so there will be swrad and potet
                  Seginc_swrad(i) = Seginc_swrad(this_seg)/Segment_hruarea(this_seg)
                  Seginc_potet(i) = Seginc_potet(this_seg)/Segment_hruarea(this_seg)
                  found = .true.
                  exit
              endif
            enddo

            if (.not. found) then
! Segment does not have any upstream segments with HRUs, check downstream segments.

              this_seg = i
              found = .false.
              do
                if (Segment_hruarea(this_seg) <= NEARZERO) then

                   ! Hit the terminal segment without finding any HRUs (i.e. sources of streamflow)
                   if (tosegment(this_seg) == OUTFLOW_SEGMENT) then
                     found = .false.
                     exit
                   endif

                   ! There is a downstream segment, check that segment for HRUs
                   this_seg = tosegment(this_seg)
                else
                    ! This segment has HRUs so there will be swrad and potet
                    Seginc_swrad(i) = Seginc_swrad(this_seg)/Segment_hruarea(this_seg)
                    Seginc_potet(i) = Seginc_potet(this_seg)/Segment_hruarea(this_seg)
                    found = .true.
                    exit
                endif
              enddo

              if (.not. found) then
!                write(*,*) "route_run: no upstream or downstream HRU found for segment ", i
!                write(*,*) "    no values for seginc_swrad and seginc_potet"
                Seginc_swrad(i) = -99.9
                Seginc_potet(i) = -99.9
              endif
            endif
          ENDIF
        ENDDO
      ENDIF

      END FUNCTION route_run

!***********************************************************************
!     routing_restart - write or read restart file
!***********************************************************************
      SUBROUTINE routing_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_ROUTING
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variables
      CHARACTER(LEN=7) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_segment_storage
        IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &       Strmflow_flag==strmflow_muskingum_mann_module ) WRITE ( Restart_outunit ) Segment_delta_flow
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_segment_storage
        IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module .OR. &
     &     Strmflow_flag==strmflow_muskingum_mann_module ) READ ( Restart_inunit ) Segment_delta_flow
      ENDIF
      END SUBROUTINE routing_restart
