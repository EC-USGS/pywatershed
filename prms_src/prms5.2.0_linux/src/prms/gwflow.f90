!***********************************************************************
! Sums inflow to and outflow from PRMS ground-water reservoirs; outflow
! can be routed to downslope ground-water reservoirs and stream
! segments
!
! Can be used for depression storage
!***********************************************************************
! Modified 7/1997 J. Vaccaro to set a minimum value for groundwater flow
! by reading in a minimum ground-water storage value for each groundwater
! reservoir, if this value is set=0, then standard PRMS routine module.
! A minimum may represent an injection well, intrabasin transfer,
! contribution from larger regional gw system, or past residual storage
! modified 10/1/2008 rsregan to include Vaccaro code
!***********************************************************************
      MODULE PRMS_GWFLOW
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, LAKE, SWALE, DEBUG_less, DNEARZERO, &
     &    DOCUMENTATION, CASCADEGW_OFF, ERROR_water_use
      USE PRMS_MODULE, ONLY: Model, Nhru, Ngw, Nlake, Print_debug, Init_vars_from_file, &
     &    Dprst_flag, Cascadegw_flag, Lake_route_flag, Inputerror_flag, Gwr_swale_flag, &
     &    Gwr_add_water_use, Gwr_transfer_water_use
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Groundwater'
      character(len=6), parameter :: MODNAME = 'gwflow'
      character(len=*), parameter :: Version_gwflow = '2020-12-02'
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gwstor_minarea(:), Gwin_dprst(:)
      DOUBLE PRECISION, SAVE :: Basin_gw_upslope
      INTEGER, SAVE :: Gwminarea_flag
      DOUBLE PRECISION, SAVE :: Basin_dnflow
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_seepage_max(:)
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_gwstor, Basin_gwflow, Basin_gwsink
      DOUBLE PRECISION, SAVE :: Basin_gwin, Basin_lake_seep
      DOUBLE PRECISION, SAVE :: Basin_gwstor_minarea_wb
      REAL, SAVE, ALLOCATABLE :: Gwres_flow(:), Gwres_sink(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gw_upslope(:), Gwres_in(:)
      REAL, SAVE, ALLOCATABLE :: Hru_gw_cascadeflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gw_in_soil(:), Gw_in_ssr(:), Hru_storage(:), Hru_lateral_flow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gwstor_minarea_wb(:), Hru_streamflow_out(:), Lakein_gwflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Lake_seepage(:), Gw_seep_lakein(:), Lake_seepage_gwr(:)
      REAL, SAVE, ALLOCATABLE :: Elevlake(:)
!   Declared Parameters
      REAL, SAVE, ALLOCATABLE :: Gwflow_coef(:), Gwsink_coef(:)
      REAL, SAVE, ALLOCATABLE :: Gwstor_init(:), Gwstor_min(:)
      REAL, SAVE, ALLOCATABLE :: Lake_seep_elev(:), Elevlake_init(:), Gw_seep_coef(:)
      END MODULE PRMS_GWFLOW

!***********************************************************************
!     Main gwflow routine
!***********************************************************************
      INTEGER FUNCTION gwflow()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, READ_INIT, SAVE_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Init_vars_from_file, Save_vars_to_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: gwflowdecl, gwflowinit, gwflowrun
      EXTERNAL gwflow_restart
!***********************************************************************
      gwflow = 0

      IF ( Process_flag==RUN ) THEN
        gwflow = gwflowrun()
      ELSEIF ( Process_flag==DECL ) THEN
        gwflow = gwflowdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL gwflow_restart(READ_INIT)
        gwflow = gwflowinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL gwflow_restart(SAVE_INIT)
      ENDIF

      END FUNCTION gwflow

!***********************************************************************
!     gwflowdecl - set up parameters for groundwater computations
!   Declared Parameters
!     gwstor_init, gwflow_coef, gwsink_coef
!     lake_seep_elev, elevlake_init, gw_seep_coef
!***********************************************************************
      INTEGER FUNCTION gwflowdecl()
      USE PRMS_GWFLOW
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      gwflowdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_gwflow)

! cascading variables and parameters
      IF ( Cascadegw_flag>CASCADEGW_OFF .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Gw_upslope(Ngw) )
        IF ( declvar(MODNAME, 'gw_upslope', 'ngw', Ngw, 'double', &
     &       'Groundwater flow received from upslope GWRs for each GWR', &
     &       'acre-inches', Gw_upslope)/=0 ) CALL read_error(3, 'gw_upslope')

        ALLOCATE ( Hru_gw_cascadeflow(Ngw) )
        IF ( declvar(MODNAME, 'hru_gw_cascadeflow', 'ngw', Ngw, 'real', &
     &       'Cascading groundwater flow from each GWR', &
     &       'inches', Hru_gw_cascadeflow)/=0 ) CALL read_error(3, 'hru_gw_cascadeflow')

        IF ( (Nlake>0.AND.Cascadegw_flag>CASCADEGW_OFF) .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Lakein_gwflow(Nlake) )
          IF ( declvar(MODNAME, 'lakein_gwflow', 'nlake', Nlake, 'double', &
     &         'Groundwater flow received from upslope GWRs for each Lake GWR', &
     &         'acre-inches', Lakein_gwflow)/=0 ) CALL read_error(3, 'lakein_gwflow')
        ENDIF
      ENDIF

      ALLOCATE ( Gwres_flow(Ngw) )
      IF ( declvar(MODNAME, 'gwres_flow', 'ngw', Ngw, 'real', &
     &     'Groundwater discharge from each GWR to the stream network', &
     &     'inches', Gwres_flow)/=0 ) CALL read_error(3, 'gwres_flow')

      ALLOCATE ( Gwres_in(Ngw) )
      IF ( declvar(MODNAME, 'gwres_in', 'ngw', Ngw, 'double', &
     &     'Total inflow to each GWR from associated capillary and gravity reservoirs', &
     &     'acre-inches', Gwres_in)/=0 ) CALL read_error(3, 'gwres_in')

      ALLOCATE ( Gwres_sink(Ngw) )
      IF ( declvar(MODNAME, 'gwres_sink', 'ngw', Ngw, 'real', &
     &     'Outflow from GWRs to the groundwater sink; water is'// &
     &     ' considered underflow or flow to deep aquifers and does'// &
     &     ' not flow to the stream network', &
     &     'inches', Gwres_sink)/=0 ) CALL read_error(3, 'gwres_sink')

      ALLOCATE ( Gw_in_soil(Ngw) )
      IF ( declvar(MODNAME, 'gw_in_soil', 'ngw', Ngw, 'double', &
     &     'Drainage from capillary reservoir excess water for each GWR', &
     &     'acre-inches', Gw_in_soil)/=0 ) CALL read_error(3, 'gw_in_soil')

      ALLOCATE ( Gw_in_ssr(Ngw) )
      IF ( declvar(MODNAME, 'gw_in_ssr', 'ngw', Ngw, 'double', &
     &     'Drainage from gravity reservoir excess water for each GWR', &
     &     'acre-inches', Gw_in_ssr)/=0 ) CALL read_error(3, 'gw_in_ssr')

      IF ( declvar(MODNAME, 'basin_gwstor', 'one', 1, 'double', &
     &     'Basin area-weighted average of storage in GWRs', &
     &     'inches', Basin_gwstor)/=0 ) CALL read_error(3, 'basin_gwstor')

      IF ( declvar(MODNAME, 'basin_gwin', 'one', 1, 'double', &
     &     'Basin area-weighted average of inflow to GWRs', &
     &     'inches', Basin_gwin)/=0 ) CALL read_error(3, 'basin_gwin')

      IF ( declvar(MODNAME, 'basin_gwflow', 'one', 1, 'double', &
     &     'Basin area-weighted average of groundwater flow to the stream network', &
     &     'inches', Basin_gwflow)/=0 ) CALL read_error(3, 'basin_gwflow')

      IF ( declvar(MODNAME, 'basin_gwsink', 'one', 1, 'double', &
     &     'Basin area-weighted average of GWR outflow to the groundwater sink', &
     &     'inches', Basin_gwsink)/=0 ) CALL read_error(3, 'basin_gwsink')

      ALLOCATE ( Hru_streamflow_out(Nhru) )
      IF ( declvar(MODNAME, 'hru_streamflow_out', 'nhru', Nhru, 'double', &
     &     'Total flow to stream network from each HRU', &
     &     'cfs', Hru_streamflow_out)/=0 ) CALL read_error(3, 'Hru_streamflow_out')

      ALLOCATE ( Hru_lateral_flow(Nhru) )
      IF ( declvar(MODNAME, 'hru_lateral_flow', 'nhru', Nhru, 'double', &
     &     'Lateral flow to stream network from each HRU', &
     &     'inches', Hru_lateral_flow)/=0 ) CALL read_error(3, 'Hru_lateral_flow')

      ALLOCATE ( Hru_storage(Nhru) )
      IF ( declvar(MODNAME, 'hru_storage', 'nhru', Nhru, 'double', &
     &     'Storage for each HRU', &
     &     'inches', Hru_storage)/=0 ) CALL read_error(3, 'hru_storage')

      ALLOCATE ( Gwstor_minarea(Ngw) )
      IF ( Dprst_flag==1 ) ALLOCATE ( Gwin_dprst(Ngw) )

      IF ( Lake_route_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'basin_lake_seep', 'one', 1, 'double', &
     &       'Basin area-weighted average of lake-bed seepage to GWRs', &
     &       'acre-feet', Basin_lake_seep)/=0 ) CALL read_error(3, 'basin_lake_seep')

        ALLOCATE ( Lake_seepage(Nlake), Lake_seepage_max(Nlake) )
        IF ( declvar(MODNAME, 'lake_seepage', 'nlake', Nlake, 'double', &
     &       'Lake-bed seepage from each lake to associated GWRs', &
     &       'acre-feet', Lake_seepage)/=0 ) CALL read_error(3, 'lake_seepage')

        ALLOCATE ( Gw_seep_lakein(Nlake) )
        IF ( declvar(MODNAME, 'gw_seep_lakein', 'nlake', Nlake, 'double', &
     &       'Groundwater discharge to any associated lake for each GWR', &
     &       'acre-feet', Gw_seep_lakein)/=0 ) CALL read_error(3, 'gw_seep_lakein')

        ALLOCATE ( Lake_seepage_gwr(Ngw) )
        IF ( declvar(MODNAME, 'lake_seepage_gwr', 'ngw', Ngw, 'double', &
     &       'Net lake-bed seepage to associated GWRs', &
     &       'inches', Lake_seepage_gwr)/=0 ) CALL read_error(3, 'lake_seepage_gwr')

        ALLOCATE ( Elevlake(Nlake) )
        IF ( declvar(MODNAME, 'elevlake', 'nlake', Nlake, 'real', &
     &       'Surface elevation of each lake', &
     &       'feet', Elevlake)/=0 ) CALL read_error(3, 'elevlake')
      ENDIF

      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==6 ) THEN
        ALLOCATE ( Gwstor_init(Ngw) )
        IF ( declparam(MODNAME, 'gwstor_init', 'ngw', 'real', &
     &       '2.0', '0.0', '50.0', &
     &       'Initial storage in each GWR', &
     &       'Storage in each GWR at the beginning of a simulation', &
     &       'inches')/=0 ) CALL read_error(1, 'gwstor_init')
      ENDIF

      ALLOCATE ( Gwflow_coef(Ngw) )
      IF ( declparam(MODNAME, 'gwflow_coef', 'ngw', 'real', &
     &     '0.015', '0.0', '0.5', &
     &     'Groundwater routing coefficient', &
     &     'Linear coefficient in the equation to compute groundwater discharge for each GWR', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'gwflow_coef')

      ALLOCATE ( Gwsink_coef(Ngw) )
      IF ( declparam(MODNAME, 'gwsink_coef', 'ngw', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Groundwater sink coefficient', &
     &     'Linear coefficient in the equation to compute outflow'// &
     &     ' to the groundwater sink for each GWR', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'gwsink_coef')

      IF ( Lake_route_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Lake_seep_elev(Nlake) )
        IF ( declparam(MODNAME, 'lake_seep_elev', 'nlake', 'real', &
     &       '1.0', '-300.0', '10000.0', &
     &       'Elevation over which lakebed seepage to the GWR occurs', &
     &       'Elevation over which lakebed seepage to the GWR occurs for lake HRUs using broad-crested weir'// &
     &       ' routing or gate opening routing', &
     &       'feet')/=0 ) CALL read_error(1, 'lake_seep_elev')

        IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
          ALLOCATE ( Elevlake_init(Nlake) )
          IF ( declparam(MODNAME, 'elevlake_init', 'nlake', 'real', &
     &         '1.0', '-300.0', '10000.0', &
     &         'Initial lake surface elevation', &
     &         'Initial lake surface elevation for each lake using broad-crested weir or gate opening routing', &
     &         'feet')/=0 ) CALL read_error(1, 'elevlake_init')
        ENDIF

        ALLOCATE ( Gw_seep_coef(Ngw) )
        IF ( declparam(MODNAME, 'gw_seep_coef', 'ngw', 'real', &
     &       '0.015', '0.001', '0.05', &
     &       'Linear coefficient to compute seepage and groundwater'// &
     &       ' discharge to and from associated lake HRUs', &
     &       'Linear coefficient in equation to compute lakebed'// &
     &       ' seepage to the GWR and groundwater discharge to each lake using broad-crested weir'// &
     &       ' routing or gate opening routing', &
     &       'fraction/day')/=0 ) CALL read_error(1, 'gw_seep_coef')
      ENDIF

      ALLOCATE ( Gwstor_min(Ngw) )
      IF ( declparam(MODNAME, 'gwstor_min', 'ngw', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Minimum storage in each GWR', &
     &     'Minimum storage in each GWR to ensure storage is greater'// &
     &     ' than specified value to account for inflow from deep'// &
     &     ' aquifers or injection wells with the water source outside the basin', &
     &     'inches')/=0 ) CALL read_error(1, 'gwstor_min')

      ALLOCATE ( Gwstor_minarea_wb(Ngw) )
      IF ( declvar(MODNAME, 'gwstor_minarea_wb', 'ngw', Ngw, 'double', &
     &     'Storage added to each GWR when storage is less than gwstor_min', &
     &     'inches', Gwstor_minarea_wb)/=0 ) CALL read_error(3, 'gwstor_minarea_wb')

      IF ( declvar(MODNAME, 'basin_gwstor_minarea_wb', 'one', 1, 'double', &
     &     'Basin area-weighted average storage added to each GWR when storage is less than gwstor_min', &
     &     'inches', Basin_gwstor_minarea_wb)/=0 ) CALL read_error(3, 'basin_gwstor_minarea_wb')

      END FUNCTION gwflowdecl

!***********************************************************************
!     gwflowinit - Initialize gwflow module - get parameter values,
!                  compute initial values.
!***********************************************************************
      INTEGER FUNCTION gwflowinit()
      USE PRMS_GWFLOW
      USE PRMS_BASIN, ONLY: Gwr_type, Hru_area, Basin_area_inv, Active_gwrs, Gwr_route_order, &
     &                      Lake_hru_id, Weir_gate_flag
      USE PRMS_FLOWVARS, ONLY: Gwres_stor, Pkwater_equiv
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SRUNOFF, ONLY: Hru_impervstor, Dprst_stor_hru
      USE PRMS_SOILZONE, ONLY: Soil_moist_tot
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: read_error
      INTRINSIC :: DBLE
! Local Variables
      INTEGER :: i, j, jjj
!***********************************************************************
      gwflowinit = 0

      IF ( getparam(MODNAME, 'gwflow_coef', Ngw, 'real', Gwflow_coef)/=0 ) CALL read_error(2, 'gwflow_coef')
      IF ( getparam(MODNAME, 'gwsink_coef', Ngw, 'real', Gwsink_coef)/=0 ) CALL read_error(2, 'gwsink_coef')
      IF ( getparam(MODNAME, 'gwstor_min', Ngw, 'real', Gwstor_min)/=0 ) CALL read_error(2, 'gwstor_min')

      Gwminarea_flag = OFF
      Gwstor_minarea = 0.0D0
      Gwstor_minarea_wb = 0.0D0
      Basin_gwstor_minarea_wb = 0.0D0
      IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==6 ) THEN
        IF ( getparam(MODNAME, 'gwstor_init', Ngw, 'real', Gwstor_init)/=0 ) CALL read_error(2, 'gwstor_init')
        Gwres_stor = DBLE( Gwstor_init )
        DEALLOCATE ( Gwstor_init )
      ENDIF
      Hru_storage = 0.0D0
      Basin_gwstor = 0.0D0
      DO j = 1, Active_gwrs
        i = Gwr_route_order(j)
        Basin_gwstor = Basin_gwstor + Gwres_stor(i)*DBLE(Hru_area(i))
        IF ( Gwstor_min(i)>0.0 ) THEN
          Gwminarea_flag = 1
          Gwstor_minarea(i) = DBLE( Gwstor_min(i)*Hru_area(i) )
        ENDIF
        IF ( Gwflow_coef(i)>1.0 ) THEN
          IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, gwflow_coef value > 1.0 for GWR:', i, Gwflow_coef(i)
        ENDIF

        ! GWR's cannot be swales unless gwr_swale_flag > 0
        IF ( Gwr_type(i)==SWALE ) THEN ! rsr, may need to add gwr_type and gwr_segment
          IF ( Gwr_swale_flag==OFF ) THEN
            PRINT *, 'ERROR, GWRs cannot be swales when gwr_swale_flag = 0, GWR:', i
            Inputerror_flag = 1
          ELSEIF ( Gwr_swale_flag==1 ) THEN
            IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, GWR:', i, ' is treated as a swale, flow sent to sink'
          ELSEIF ( Gwr_swale_flag==2 ) THEN
            IF ( Print_debug>DEBUG_less ) PRINT *, 'WARNING, GWR:', i, &
       &         ' is treated as a swale, flow sent to basin_cfs and hru_segment if > 0'
          ELSE
! maybe gwr_swale_flag = 3 abs(hru_segment) so hru_segment could be changed from 0 to allow HRU swales
            PRINT *, 'ERROR, invalid gwr_swale_flag value, specified as:', gwr_swale_flag
            Inputerror_flag = 1
          ENDIF
        ENDIF

        Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Gwres_stor(i) + Pkwater_equiv(i)
        IF ( Dprst_flag==ACTIVE ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)
      ENDDO
      IF ( Gwminarea_flag==OFF ) DEALLOCATE ( Gwstor_minarea )
      Basin_gwstor = Basin_gwstor*Basin_area_inv

      IF ( Dprst_flag==ACTIVE ) Gwin_dprst = 0.0D0

      IF ( Weir_gate_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'gw_seep_coef', Ngw, 'real', Gw_seep_coef)/=0 ) CALL read_error(2, 'gw_seep_coef')
        IF ( getparam(MODNAME, 'lake_seep_elev', Nlake, 'real', Lake_seep_elev)/=0 ) CALL read_error(2, 'lake_seep_elev')
        IF ( Init_vars_from_file==0 .OR. Init_vars_from_file==2 .OR. Init_vars_from_file==4 ) THEN
          IF ( getparam(MODNAME, 'elevlake_init', Nlake, 'real', Elevlake_init)/=0 ) CALL read_error(2, 'elevlake_init')
          Elevlake = Elevlake_init
          DEALLOCATE ( Elevlake_init )
        ENDIF
        Lake_seepage_max = 0.0D0
        Lake_seepage_gwr = 0.0D0
        Lake_seepage = 0.0D0
        Gw_seep_lakein = 0.0D0
        DO i = 1, Active_gwrs
          j = Gwr_route_order(i)
          IF ( Gwr_type(j)==LAKE ) THEN
            jjj = Lake_hru_id(j)
            IF ( jjj==0 ) THEN
              PRINT *, 'ERROR, GWR specified as a lake but lake_hru_id value = 0, GWR:', j
              Inputerror_flag = 1
              CYCLE
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      Basin_gwflow = 0.0D0
      Basin_gwsink = 0.0D0
      Basin_gwin = 0.0D0
      Basin_gw_upslope = 0.0D0
      Basin_dnflow = 0.0D0
      Basin_lake_seep = 0.0D0
! do only once, so restart uses saved values
      IF ( Cascadegw_flag>CASCADEGW_OFF ) THEN
        Gw_upslope = 0.0D0
        Hru_gw_cascadeflow = 0.0
        IF ( Nlake>0 ) Lakein_gwflow = 0.0D0
      ENDIF
      Gwres_flow = 0.0
      Gwres_in = 0.0
      Gwres_sink = 0.0
      Gw_in_ssr = 0.0D0
      Gw_in_soil = 0.0D0
      Hru_streamflow_out = 0.0D0
      Hru_lateral_flow = 0.0D0

      END FUNCTION gwflowinit

!***********************************************************************
!     gwflowrun - Computes groundwater flow to streamflow and to
!                 groundwater sink
!***********************************************************************
      INTEGER FUNCTION gwflowrun()
      USE PRMS_GWFLOW
      USE PRMS_BASIN, ONLY: Active_gwrs, Gwr_route_order, Lake_type, &
     &    Basin_area_inv, Hru_area, Gwr_type, Lake_hru_id, Weir_gate_flag, Hru_area_dble
      USE PRMS_FLOWVARS, ONLY: Soil_to_gw, Ssr_to_gw, Sroff, Ssres_flow, Gwres_stor, Pkwater_equiv, Lake_vol
      USE PRMS_CASCADE, ONLY: Ncascade_gwr
      USE PRMS_SET_TIME, ONLY: Cfs_conv, Nowyear, Nowmonth, Nowday
      USE PRMS_SRUNOFF, ONLY: Dprst_seep_hru, Hru_impervstor, Dprst_stor_hru
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SOILZONE, ONLY: Soil_moist_tot
      USE PRMS_WATER_USE, ONLY: Gwr_transfer, Gwr_gain
      IMPLICIT NONE
! Functions
      EXTERNAL :: rungw_cascade, print_date
      INTRINSIC :: DBLE, DABS, SNGL, MIN
! Local Variables
      INTEGER :: i, j, jj, jjj
      REAL :: dnflow
      DOUBLE PRECISION :: gwin, gwstor, gwsink, gwflow, gwstor_last, seepage, gwarea, inch2acre_feet
!***********************************************************************
      gwflowrun = 0

      IF ( Cascadegw_flag>CASCADEGW_OFF ) THEN
        Gw_upslope = 0.0D0
        Basin_dnflow = 0.0D0
        Basin_gw_upslope = 0.0D0
        IF ( Nlake>0 ) Lakein_gwflow = 0.0D0
      ENDIF

      IF ( Weir_gate_flag==ACTIVE ) THEN
        ! elevlake from last timestep
        Lake_seepage = 0.0D0
        Gw_seep_lakein = 0.0D0
        Basin_lake_seep = 0.0D0
        DO jj = 1, Active_gwrs
          j = Gwr_route_order(jj)
          IF ( Gwr_type(j)/=LAKE ) CYCLE
          jjj = Lake_hru_id(j)
          ! what happens when lake goes dry? need lake bottom elevation ?
          IF ( Lake_type(jjj)==4 .OR. Lake_type(jjj)==5 ) &
     &         Lake_seepage_max(jjj) = DBLE( (Elevlake(jjj)-Lake_seep_elev(jjj))*12.0*Gw_seep_coef(j) ) ! units = inches
        ENDDO
        DO jj = 1, Active_gwrs
          j = Gwr_route_order(jj)
          IF ( Gwr_type(j)==LAKE ) THEN ! only if a weir gate lake
            jjj = Lake_hru_id(j) !! jjj must be > zero due to check above
            IF ( Lake_type(jjj)==4 .OR. Lake_type(jjj)==5 ) THEN
              inch2acre_feet = Hru_area_dble(j)/12.0D0
              seepage = Lake_seepage_max(jjj)
              ! seepage added to GWR
              !rsr, need seepage variable for WB
              IF ( seepage<0.0D0 ) THEN
! water to lake from GWR, negative value of seepage
                IF ( DABS(seepage)>Gwres_stor(j) ) THEN
                  IF ( Print_debug>DEBUG_less ) THEN
                    PRINT *, 'WARNING, GWR storage insufficient for discharge to lake:', jjj, ' GWR:', j
                    CALL print_date(1)
                    PRINT *, 'GWR storage:', Gwres_stor(j), ', computed discharge:', seepage
                    PRINT *, 'Discharge set to available GWR storage'
!?? adjust lake storage and elevation
                    PRINT *, 'Lake elevation, storage, and water balance not adjusted'
                  ENDIF
                  seepage = -Gwres_stor(j)
                ENDIF
                Gw_seep_lakein(jjj) = Gw_seep_lakein(jjj) - seepage*inch2acre_feet ! units, acre-feet
              ELSE
! water from lake to GWR, positive value of seepage
                ! needed because lakes could go dry
! DANGER, needed for multiple HRU lakes, if lake goes dry some GWRs won't receive seepage, withdrawn from lake based on route order
                IF ( Lake_vol(jjj)<seepage*inch2acre_feet ) seepage = Lake_vol(jjj)/inch2acre_feet
                Lake_seepage(jjj) = Lake_seepage(jjj) + seepage*inch2acre_feet ! units, acre-feet
                Basin_lake_seep = Basin_lake_seep + seepage*Hru_area_dble(j) ! units, acre-inches
              ENDIF
              Gwres_stor(j) = Gwres_stor(j) + seepage
              Lake_seepage_gwr(j) = seepage ! can be positive or negative value
              Lake_vol(jjj) = Lake_vol(jjj) - seepage*inch2acre_feet
            ENDIF
          ENDIF
        ENDDO
        Basin_lake_seep = Basin_lake_seep*Basin_area_inv
      ENDIF

      Basin_gwstor_minarea_wb = 0.0D0
      Basin_gwflow = 0.0D0
      Basin_gwstor = 0.0D0
      Basin_gwsink = 0.0D0
      Basin_gwin = 0.0D0
      DO j = 1, Active_gwrs
        i = Gwr_route_order(j)
        gwarea = Hru_area_dble(i)
        gwstor = Gwres_stor(i)*gwarea ! acre-inches
        ! soil_to_gw is for whole HRU, not just perv
        Gw_in_soil(i) = Soil_to_gw(i)*Hru_area(i)
        Gw_in_ssr(i) = Ssr_to_gw(i)*Hru_area(i)
        gwin = Gw_in_soil(i) + Gw_in_ssr(i)
        IF ( Cascadegw_flag>CASCADEGW_OFF ) THEN
          gwin = gwin + Gw_upslope(i)
          Basin_gw_upslope = Basin_gw_upslope + Gw_upslope(i)
        ENDIF
        IF ( Dprst_flag==ACTIVE ) THEN
          !rsr, need basin variable for WB
          Gwin_dprst(i) = Dprst_seep_hru(i)*gwarea
          gwin = gwin + Gwin_dprst(i)
        ENDIF
        IF ( Gwr_add_water_use==ACTIVE ) THEN
          IF ( Gwr_gain(i)>0.0 ) gwin = gwin + DBLE(Gwr_gain(i))/Cfs_conv
        ENDIF
        gwstor = gwstor + gwin
        Basin_gwin = Basin_gwin + gwin
        IF ( Gwminarea_flag==ACTIVE ) THEN
          ! check to be sure gwres_stor >= gwstor_minarea before computing outflows
          IF ( gwstor<Gwstor_minarea(i) ) THEN
            IF ( gwstor<0.0D0 ) THEN
              IF ( Print_debug>DEBUG_less ) PRINT *, 'Warning, groundwater reservoir for HRU:', i, &
     &                                       ' is < 0.0 with gwstor_min active', gwstor
!              ERROR STOP ERROR_var
            ENDIF
            gwstor_last = gwstor
            gwstor = Gwstor_minarea(i)
            !rsr, keep track of change in storage for WB
            Gwstor_minarea_wb(i) = gwstor - gwstor_last
            Basin_gwstor_minarea_wb = Basin_gwstor_minarea_wb + Gwstor_minarea_wb(i)
            Gwstor_minarea_wb(i) = Gwstor_minarea_wb(i)/gwarea
            IF ( Print_debug>DEBUG_less ) PRINT *, 'Added to gwres_stor as storage < gwstor_min to GWR:', i, &
     &                                             ' amount:', Gwstor_minarea_wb(i)
          ELSE
            Gwstor_minarea_wb(i) = 0.0D0
          ENDIF
        ENDIF

        IF ( Gwr_transfer_water_use==ACTIVE ) THEN
          IF ( Gwr_transfer(i)>0.0 ) THEN
            IF ( SNGL(gwstor*Cfs_conv)<Gwr_transfer(i) ) THEN
              PRINT *, 'ERROR, not enough storage for transfer from groundwater reservoir storage:', &
     &                  i, ' Date:', Nowyear, Nowmonth, Nowday
              PRINT *, '       storage: ', gwstor, '; transfer: ', Gwr_transfer(i)/Cfs_conv
              ERROR STOP ERROR_water_use
            ENDIF
            gwstor = gwstor - DBLE( Gwr_transfer(i) ) / Cfs_conv 
          ENDIF
        ENDIF
 
        gwsink = 0.0D0
        IF ( gwstor<0.0D0 ) THEN ! could happen with water use
          IF ( Print_debug>DEBUG_less ) PRINT *, 'Warning, groundwater reservoir for HRU:', i, ' is < 0.0', gwstor
          gwflow = 0.0D0
          Gwres_sink(i) = 0.0
        ELSE

! Compute groundwater discharge
          gwflow = gwstor*DBLE( Gwflow_coef(i) )

! Reduce storage by outflow
          gwstor = gwstor - gwflow

          IF ( Gwsink_coef(i)>0.0 ) THEN
            gwsink = MIN( gwstor*DBLE( Gwsink_coef(i) ), gwstor ) ! if gwsink_coef > 1, could have had negative gwstor
            gwstor = gwstor - gwsink
          ENDIF
! if gwr_swale_flag = 1 swale GWR flow goes to sink, 2 included in stream network and cascades
! maybe gwr_swale_flag = 3 abs(hru_segment) so hru_segment could be changed from 0 to allow HRU swales
          IF ( Gwr_swale_flag==ACTIVE ) THEN
            IF ( Gwr_type(i)==SWALE ) THEN
              gwsink = gwsink + gwflow
              gwflow = 0.0D0
            ENDIF
          ENDIF
          Gwres_sink(i) = SNGL( gwsink/gwarea )
          Basin_gwsink = Basin_gwsink + gwsink
        ENDIF
        Basin_gwstor = Basin_gwstor + gwstor

        Gwres_flow(i) = SNGL( gwflow/gwarea )
        IF ( Cascadegw_flag>CASCADEGW_OFF ) THEN
          IF ( Ncascade_gwr(i)>0 ) THEN
            CALL rungw_cascade(i, Ncascade_gwr(i), Gwres_flow(i), dnflow)
            Hru_gw_cascadeflow(i) = dnflow
            Basin_dnflow = Basin_dnflow + dnflow*gwarea
          ELSEIF ( Gwr_type(i)==LAKE ) THEN
            Lakein_gwflow(Lake_hru_id(i)) = Lakein_gwflow(Lake_hru_id(i)) + Gwres_flow(i)
          ENDIF
        ENDIF
        Basin_gwflow = Basin_gwflow + DBLE(Gwres_flow(i))*gwarea

        ! leave gwin in inch-acres
        Gwres_in(i) = gwin
        Gwres_stor(i) = gwstor/gwarea
        Hru_lateral_flow(i) = DBLE( Gwres_flow(i) + Sroff(i) + Ssres_flow(i) )
        ! Cfs_conv converts acre-inches per timestep to cfs
        Hru_streamflow_out(i) = gwarea*Cfs_conv*Hru_lateral_flow(i)
        Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Gwres_stor(i) &
     &                   + Pkwater_equiv(i)
        IF ( Dprst_flag==ACTIVE ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)
      ENDDO

      Basin_gwflow = Basin_gwflow*Basin_area_inv
      Basin_gwstor = Basin_gwstor*Basin_area_inv
      Basin_gwsink = Basin_gwsink*Basin_area_inv
      Basin_gwin = Basin_gwin*Basin_area_inv
      Basin_gw_upslope = Basin_gw_upslope*Basin_area_inv
      Basin_gwstor_minarea_wb = Basin_gwstor_minarea_wb*Basin_area_inv
      Basin_dnflow = Basin_dnflow*Basin_area_inv

      END FUNCTION gwflowrun

!***********************************************************************
!     Compute cascading GW flow
!***********************************************************************
      SUBROUTINE rungw_cascade(Igwr, Ncascade_gwr, Gwres_flow, Dnflow)
      USE PRMS_SRUNOFF, ONLY: Strm_seg_in
      USE PRMS_GWFLOW, ONLY: Gw_upslope
      USE PRMS_CASCADE, ONLY: Gwr_down, Gwr_down_frac, Cascade_gwr_area
      ! Cfs_conv converts acre-inches per timestep to cfs
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      IMPLICIT NONE
! Functions
      INTRINSIC :: IABS, DBLE
! Arguments
      INTEGER, INTENT(IN) :: Igwr, Ncascade_gwr
      REAL, INTENT(INOUT) :: Gwres_flow
      REAL, INTENT(OUT) :: Dnflow
! Local variables
      INTEGER :: j, k
!***********************************************************************
      Dnflow = 0.0
      DO k = 1, Ncascade_gwr
        j = Gwr_down(k, Igwr)
        ! Gwres_flow is in inches
! if gwr_down(k, Igwr) > 0, cascade contributes to a downslope GWR
        IF ( j>0 ) THEN
          Gw_upslope(j) = Gw_upslope(j) + Gwres_flow*Cascade_gwr_area(k, Igwr)
          Dnflow = Dnflow + Gwres_flow*Gwr_down_frac(k, Igwr)
! if gwr_down(k, Igwr) < 0, cascade contributes to a stream
        ELSEIF ( j<0 ) THEN
          j = IABS( j )
          Strm_seg_in(j) = Strm_seg_in(j) + DBLE( Gwres_flow*Cascade_gwr_area(k, Igwr) )*Cfs_conv
        ENDIF
      ENDDO

      ! gwres_flow reduced by cascading flow to HRUs
      Gwres_flow = Gwres_flow - Dnflow
      IF ( Gwres_flow<0.0 ) Gwres_flow = 0.0

      END SUBROUTINE rungw_cascade

!***********************************************************************
!     gwflow_restart - write or read gwflow restart file
!***********************************************************************
      SUBROUTINE gwflow_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_BASIN, ONLY: Weir_gate_flag
      USE PRMS_GWFLOW
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=6) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        IF ( Weir_gate_flag==ACTIVE ) WRITE ( Restart_outunit ) Elevlake
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        IF ( Weir_gate_flag==ACTIVE ) READ ( Restart_inunit ) Elevlake ! could be error if someone turns off weirs for restart
      ENDIF
      END SUBROUTINE gwflow_restart
