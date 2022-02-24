!***********************************************************************
! Computes streamflow and other variables in subbasins
!
!     version: 3.1 (rsregan, April 2008 ) changed code to assume
!                  nhru=nssr=ngw
!     version: 3.0 (rsregan, April 2008 ) removed all but subbasin code
!     version: 2.1 (rsregan, June 2007)
!     version: 2 (lhay, May 2007)
!     version: 1.1 (rsregan, March 2007)
!     version: 1 (lhay, March 2007)
!***********************************************************************

      MODULE PRMS_SUBBASIN
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Output Summary'
      character(len=*), parameter :: MODNAME = 'subbasin'
      character(len=*), parameter :: Version_subbasin = '2022-01-25'
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Qsub(:), Sub_area(:), Laststor(:)
      INTEGER, SAVE, ALLOCATABLE :: Tree(:, :)
!   Declared Variables
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Sub_inq(:), Sub_cfs(:), Sub_cms(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Sub_interflow(:), Sub_gwflow(:), Sub_sroff(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_interflow(:), Subinc_gwflow(:), Subinc_sroff(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_precip(:), Subinc_snowmelt(:), Subinc_pkweqv(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_actet(:), Subinc_potet(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_swrad(:), Subinc_snowcov(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_tmaxc(:), Subinc_tminc(:), Subinc_tavgc(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_deltastor(:), Subinc_wb(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_rain(:), Subinc_snow(:), Subinc_stor(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Subinc_recharge(:), Subinc_szstor_frac(:), Subinc_capstor_frac(:)
!   Declared Parameters
      INTEGER, SAVE, ALLOCATABLE :: Subbasin_down(:), Hru_subbasin(:)
      END MODULE PRMS_SUBBASIN

!***********************************************************************
!     Main daily stream flow routine
!***********************************************************************
      INTEGER FUNCTION subbasin()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT
      USE PRMS_MODULE, ONLY: Process_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: subdecl, subinit, subrun
!***********************************************************************
      subbasin = 0

      IF ( Process_flag==RUN ) THEN
        subbasin = subrun()
      ELSEIF ( Process_flag==DECL ) THEN
        subbasin = subdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        subbasin = subinit()
      ENDIF

      END FUNCTION subbasin

!***********************************************************************
!     subdecl - set up parameters for streamflow, lake computations,
!               and subbasin flow
!   Declared Parameters
!     hru_area, subbasin_down, hru_subbasin
!***********************************************************************
      INTEGER FUNCTION subdecl()
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ERROR_dim
      use PRMS_MMFAPI, only: declvar_dble
      use PRMS_READ_PARAM_FILE, only: declparam
      USE PRMS_MODULE, ONLY: Model, Nsub, Nhru
      USE PRMS_SUBBASIN
      use prms_utils, only: error_stop, print_module, read_error
      IMPLICIT NONE
!***********************************************************************
      subdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_subbasin)

      IF ( Nsub==0 ) THEN
        IF ( Model/=DOCUMENTATION ) CALL error_stop('nsub=0 when subbasin module called', ERROR_dim)
        Nsub = 1
      ENDIF

! Declared Variables
      ALLOCATE ( Sub_interflow(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_interflow', 'nsub', Nsub, &
     &     'Area-weighted average interflow to each subbasin from associated HRUs and from upstream subbasins', &
     &     'cfs', Sub_interflow)

      ALLOCATE ( Sub_gwflow(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_gwflow', 'nsub', Nsub, &
     &     'Area-weighted average groundwater discharge from associated GWRs to each subbasin and from upstream subbasins', &
     &     'cfs', Sub_gwflow)
      ALLOCATE ( Subinc_gwflow(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_gwflow', 'nsub', Nsub, &
     &     'Area-weighted average groundwater discharge from associated GWRs to each subbasin', &
     &     'cfs', Subinc_gwflow)

      ALLOCATE ( Sub_sroff(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_sroff', 'nsub', Nsub, &
     &     'Area-weighted average surface runoff from associated HRUs to each subbasin and from upstream subbasins', &
     &     'cfs', Sub_sroff)

      ALLOCATE ( Subinc_snowcov(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_snowcov', 'nsub', Nsub, &
     &     'Area-weighted average snow-covered area for associated HRUs for each subbasin', &
     &     'decimal fraction', Subinc_snowcov)

      ALLOCATE ( Subinc_interflow(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_interflow', 'nsub', Nsub, &
     &     'Area-weighted average interflow from associated HRUs to each subbasin', &
     &     'cfs', Subinc_interflow)

      ALLOCATE ( Subinc_sroff(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_sroff', 'nsub', Nsub, &
     &    'Area-weighted average surface runoff from associated HRUs to each subbasin', &
     &    'cfs', Subinc_sroff)

      ALLOCATE ( Subinc_precip(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_precip', 'nsub', Nsub, &
     &     'Area-weighted average precipitation from associated HRUs to each subbasin', &
     &     'inches', Subinc_precip)

      ALLOCATE ( Subinc_rain(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_rain', 'nsub', Nsub, &
     &     'Area-weighted average rain on associated HRUs to each subbasin', &
     &     'inches', Subinc_rain)

      ALLOCATE ( Subinc_snow(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_snow', 'nsub', Nsub, &
     &     'Area-weighted average snow on associated HRUs to each subbasin', &
     &     'inches', Subinc_snow)

      ALLOCATE ( Subinc_actet(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_actet', 'nsub', Nsub, &
     &     'Area-weighted average actual ET from associated HRUs to each subbasin', &
     &     'inches', Subinc_actet)

      ALLOCATE ( Subinc_potet(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_potet', 'nsub', Nsub, &
     &     'Area-weighted average potential ET from associated HRUs to each subbasin', &
     &     'inches', Subinc_potet)

      ALLOCATE ( Subinc_swrad(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_swrad', 'nsub', Nsub, &
     &     'Area-weighted average shortwave radiation distributed to associated HRUs of each subbasin', &
     &     'Langleys', Subinc_swrad)

      ALLOCATE ( Subinc_tminc(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_tminc', 'nsub', Nsub, &
     &     'Area-weighted average minimum air temperature for associated HRUs to each subbasin', &
     &     'degrees Celsius', Subinc_tminc)

      ALLOCATE ( Subinc_tmaxc(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_tmaxc', 'nsub', Nsub, &
     &     'Area-weighted average maximum air temperature for associated HRUs to each subbasin', &
     &     'degrees Celsius', Subinc_tmaxc)

      ALLOCATE ( Subinc_tavgc(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_tavgc', 'nsub', Nsub, &
     &     'Area-weighted average air temperature for associated HRUs to each subbasin', &
     &     'degrees Celsius', Subinc_tavgc)

      ALLOCATE ( Subinc_wb(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_wb', 'nsub', Nsub, &
     &     'Water balance for each subbasin', &
     &     'inches', Subinc_wb)

      ALLOCATE ( Subinc_deltastor(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_deltastor', 'nsub', Nsub, &
     &     'Change in storage for each subbasin', &
     &     'inches', Subinc_deltastor)

      ALLOCATE ( Subinc_snowmelt(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_snowmelt', 'nsub', Nsub, &
     &     'Area-weighted average snowmelt from associated HRUs of each subbasin', &
     &     'inches', Subinc_snowmelt)

      ALLOCATE ( Subinc_pkweqv(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_pkweqv', 'nsub', Nsub, &
     &     'Area-weighted average snowpack water equivalent from associated HRUs of each subbasin', &
     &     'inches', Subinc_pkweqv)

      ALLOCATE ( Subinc_recharge(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_recharge', 'nsub', Nsub, &
     &     'Area-weighted average recharge from associated HRUs of each subbasin', &
     &     'inches', Subinc_recharge)

      ALLOCATE ( Subinc_szstor_frac(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_szstor_frac', 'nsub', Nsub, &
     &     'Area-weighted average fraction of soil-zone water content storage for associated HRUs of each subbasin', &
     &     'decimal fraction', Subinc_szstor_frac)

      ALLOCATE ( Subinc_capstor_frac(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_capstor_frac', 'nsub', Nsub, &
     &     'Area-weighted average fraction of capillary reservoir water content storage for associated HRUs of each subbasin', &
     &     'decimal fraction', Subinc_capstor_frac)

      ALLOCATE ( Subinc_stor(Nsub) )
      CALL declvar_dble(MODNAME, 'subinc_stor', 'nsub', Nsub, &
     &     'Area-weighted average total water content in storage reservoirs for associated HRUs of each subbasin', &
     &     'inches', subinc_stor)

      ALLOCATE ( Qsub(Nsub), Tree(Nsub, Nsub), Sub_inq(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_inq', 'nsub', Nsub, &
     &     'Sum of streamflow from upstream subbasins to each subbasin', &
     &     'cfs', Sub_inq)

      ALLOCATE ( Sub_cfs(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_cfs', 'nsub', Nsub, &
     &     'Total streamflow leaving each subbasin', &
     &     'cfs', Sub_cfs)

      ALLOCATE ( Sub_cms(Nsub) )
      CALL declvar_dble(MODNAME, 'sub_cms', 'nsub', Nsub, &
     &     'Total streamflow leaving each subbasin', &
     &     'cms', Sub_cms)

      ALLOCATE ( Subbasin_down(Nsub) )
      IF ( declparam(MODNAME, 'subbasin_down', 'nsub', 'integer', &
     &     '0', 'bounded', 'nsub', &
     &     'Downstream subbasin for each subbasin', &
     &     'Index number for the downstream subbasin whose inflow is outflow from this subbasin', &
     &     'none')/=0 ) CALL read_error(1, 'subbasin_down')

      ALLOCATE ( Hru_subbasin(Nhru) )
      IF ( declparam(MODNAME, 'hru_subbasin', 'nhru', 'integer', &
     &     '0', 'bounded', 'nsub', &
     &     'Index of subbasin assigned to each HRU', &
     &     'Index of subbasin assigned to each HRU', &
     &     'none')/=0 ) CALL read_error(1, 'hru_subbasin')

! Allocate arrays for variables
      ALLOCATE ( Sub_area(Nsub), Laststor(Nsub) )

      END FUNCTION subdecl

!***********************************************************************
!     subinit - Initialize subbasin module - get parameter values,
!               compute initial values
!***********************************************************************
      INTEGER FUNCTION subinit()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, CFS2CMS_CONV, LAKE, DNEARZERO
      use PRMS_READ_PARAM_FILE, only: getparam_int
      USE PRMS_MODULE, ONLY: Nsub, Nhru, Print_debug, &
     &    Inputerror_flag, Dprst_flag, Lake_route_flag, Cascade_flag, Hru_type
      USE PRMS_SUBBASIN
      USE PRMS_BASIN, ONLY: Hru_area_dble, Active_hrus, Hru_route_order, Hru_frac_perv, Lake_hru_id
      USE PRMS_FLOWVARS, ONLY: Ssres_stor, Soil_moist, Pkwater_equiv, Gwres_stor, Sroff, Ssres_flow, Lake_vol
      USE PRMS_SET_TIME, ONLY: Cfs_conv, Cfs2inches
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SRUNOFF, ONLY: Hru_impervstor, Hortonian_lakes, Dprst_stor_hru
      USE PRMS_SOILZONE, ONLY: Lakein_sz
      USE PRMS_GWFLOW, ONLY: Gwres_flow
      USE PRMS_MUSKINGUM_LAKE, ONLY: Lake_outcfs
      use prms_utils, only: PRMS_open_module_file, read_error
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
! Local Variables
      INTEGER :: i, j, k, kk, TREEUNIT
      DOUBLE PRECISION :: harea, gwstor, soilstor, snowstor, landstor, srq, ssq, gwq
!***********************************************************************
      subinit = 0

      IF ( getparam_int(MODNAME, 'hru_subbasin', Nhru, Hru_subbasin)/=0 ) CALL read_error(2, 'hru_subbasin')
      IF ( getparam_int(MODNAME, 'subbasin_down', Nsub, Subbasin_down)/=0 ) CALL read_error(2, 'subbasin_down')

! Determine the tree structure for the internal nodes
      Tree = 0

      DO j = 1, Nsub
        kk = Subbasin_down(j)
        IF ( kk/=0 ) Tree(kk, j) = 1
      ENDDO

      IF ( Print_debug==14 ) THEN
        CALL PRMS_open_module_file(TREEUNIT, 'tree_structure')
        WRITE ( TREEUNIT, 9001 ) (j, j=1, Nsub)
        DO j = 1, Nsub
          WRITE ( TREEUNIT, 9002 ) j, (Tree(j,k), k=1, Nsub)
        ENDDO
      ENDIF

      DO j = 1, Nsub
        DO k = 1, Nsub
          IF ( Tree(j,k)>0 ) THEN
            DO i = 1, Nsub
              Tree(j, i) = Tree(j, i) + Tree(k, i)
              IF ( Tree(j,i)>1 ) THEN
                Tree(j, i) = 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      DO j = Nsub, 1, -1
        DO k = 1, Nsub
          IF ( Tree(j,k)>0 ) THEN
            DO i = 1, Nsub
              Tree(j, i) = Tree(j, i) + Tree(k, i)
              IF ( Tree(j,i)>1 ) Tree(j,i)=1
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      IF ( Print_debug==14 ) THEN
        WRITE ( TREEUNIT, 9003 ) (j, j=1, Nsub)
        DO j = 1, Nsub
          WRITE ( TREEUNIT, 9002 ) j, (Tree(j,k), k=1, Nsub)
        ENDDO
        CLOSE ( TREEUNIT )
      ENDIF

! added some code to allow for restart, but not climate states and fluxes and subinc_deltastor

      Subinc_interflow = 0.0D0
      Subinc_gwflow = 0.0D0
      Subinc_sroff = 0.0D0
      Subinc_stor = 0.0D0
      Sub_area = 0.0D0
      gwstor = 0.0D0
      gwq = 0.0D0
      Qsub = 0.0D0
      DO i = 1, Active_hrus
        j = Hru_route_order(i)
        ! k indicates which HRU is in which subbasin
        k = Hru_subbasin(j)
        IF ( k>0 ) THEN
          harea = Hru_area_dble(j)
          gwstor = Gwres_stor(j)*harea
          gwq = DBLE(Gwres_flow(j))*harea
          IF ( Hru_type(j)/=LAKE ) THEN
            srq = DBLE(Sroff(j))*harea
            ssq = DBLE(Ssres_flow(j))*harea
            soilstor = DBLE(Soil_moist(j)*Hru_frac_perv(j) + Ssres_stor(j))*harea
            snowstor = Pkwater_equiv(j)*harea
            landstor = DBLE(Hru_intcpstor(j)+Hru_impervstor(j))*harea
            IF ( Dprst_flag==ACTIVE ) landstor = landstor + Dprst_stor_hru(j)*harea
          ELSE
            soilstor = 0.0D0
            snowstor = 0.0D0
            landstor = 0.0D0
            ! wrong if multiple HRUs for any lake
            IF ( Lake_route_flag==ACTIVE ) THEN
              landstor = Lake_vol(Lake_hru_id(j))*12.0D0
              srq = Lake_outcfs(Lake_hru_id(j))*Cfs2inches
              ssq = 0.0D0
            ELSEIF ( Cascade_flag>0 ) THEN
              srq = Hortonian_lakes(j)*harea
              ssq = Lakein_sz(j)*harea
            ELSE
              srq = 0.0D0
              ssq = 0.0D0
            ENDIF
          ENDIF
          Qsub(k) = Qsub(k) + srq + ssq + gwq
          Subinc_interflow(k) = Subinc_interflow(k) + ssq
          Subinc_sroff(k) = Subinc_sroff(k) + srq
          Subinc_gwflow(k) = Subinc_gwflow(k) + gwq
          Subinc_stor(k) = Subinc_stor(k) + soilstor + gwstor + snowstor + landstor
          Sub_area(k) = Sub_area(k) + harea
        ENDIF
      ENDDO

      DO i = 1, Nsub
        IF ( Sub_area(i)<DNEARZERO ) THEN
          PRINT *, 'ERROR, subbasin:', i, ' does not include any HRUs'
          Inputerror_flag = 1
        ELSE
          Subinc_stor(i) = Subinc_stor(i)/Sub_area(i)
        ENDIF
        Sub_inq(i) = Qsub(i)*Cfs_conv
        Subinc_interflow(i) = Subinc_interflow(i)*Cfs_conv
        Subinc_gwflow(i) = Subinc_gwflow(i)*Cfs_conv
        Subinc_sroff(i) = Subinc_sroff(i)*Cfs_conv
        ! water balance off if lake or muskingum routing
      ENDDO

! allow for possible restart
      !get cumulative subbasin flows
      DO j = 1, Nsub
        Sub_gwflow(j) = Subinc_gwflow(j)
        Sub_sroff(j) = Subinc_sroff(j)
        Sub_interflow(j) = Subinc_interflow(j)
        Sub_cfs(j) = Sub_inq(j)
        DO k = 1, Nsub
          IF ( Tree(j,k)/=0 ) THEN
            Sub_gwflow(j) = Sub_gwflow(j) + Subinc_gwflow(k)
            Sub_sroff(j) = Sub_sroff(j) + Subinc_sroff(k)
            Sub_interflow(j) = Sub_interflow(j) + Subinc_interflow(k)
            Sub_cfs(j) = Sub_cfs(j) + Sub_inq(k)
          ENDIF
        ENDDO
        Sub_cms(j) = Sub_cfs(j)*CFS2CMS_CONV
      ENDDO

 9001 FORMAT ('Initial Tree Structure for Internal Subbasins', /, &
     &        ' 1 indicates a node that flows directly into subbasin', /, 5X, 500I3)
 9002 FORMAT (I3, 2X, 500I3)
 9003 FORMAT (/, 'Tree Structure for Internal Subbasins', /, &
     &        ' 1 indicates a node that eventually flows into subbasin', /, 5X, 500I3)

      END FUNCTION subinit

!***********************************************************************
!     subrun - Computes basin streamflow and on-channel reservoir
!                  storage and outflows
!***********************************************************************
      INTEGER FUNCTION subrun()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, CFS2CMS_CONV, LAKE, CASCADE_OFF
      USE PRMS_MODULE, ONLY: Nsub, Dprst_flag, Lake_route_flag, Cascade_flag, Hru_type
      USE PRMS_SUBBASIN
      USE PRMS_BASIN, ONLY: Hru_area_dble, Active_hrus, Hru_route_order, Hru_frac_perv, Lake_hru_id
      USE PRMS_SET_TIME, ONLY: Cfs_conv, Cfs2inches
      USE PRMS_CLIMATEVARS, ONLY: Hru_ppt, Swrad, Potet, Tminc, Tmaxc, Tavgc, Hru_rain, Hru_snow
      USE PRMS_FLOWVARS, ONLY: Hru_actet, Ssres_flow, Sroff, Snowcov_area, Snowmelt, &
     &    Ssres_stor, Soil_moist, Pkwater_equiv, Gwres_stor, Lake_vol, Soil_moist, Soil_moist_max
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SRUNOFF, ONLY: Hru_impervstor, Hortonian_lakes, Dprst_stor_hru
      USE PRMS_SOILZONE, ONLY: Lakein_sz, Recharge, Soil_moist_tot, Soil_zone_max
      USE PRMS_GWFLOW, ONLY: Gwres_flow
      USE PRMS_MUSKINGUM_LAKE, ONLY: Lake_outcfs
      IMPLICIT NONE
! Functions
      INTRINSIC :: DBLE
! Local Variables
      INTEGER :: j, jj, k
      DOUBLE PRECISION :: harea, srq, ssq, gwq, dmy, dmy1, subarea
      DOUBLE PRECISION :: soilstor, snowstor, landstor, dmy2
!***********************************************************************
      subrun = 0

!   Daily time step.
!   Compute reservoir routing and basin outflow

      Qsub = 0.0D0
      Subinc_interflow = 0.0D0
      Subinc_sroff = 0.0D0
      Subinc_precip = 0.0D0
      Subinc_rain = 0.0D0
      Subinc_snow = 0.0D0
      Subinc_snowmelt = 0.0D0
      Subinc_pkweqv = 0.0D0
      Subinc_actet = 0.0D0
      Subinc_snowcov = 0.0D0
      Subinc_swrad = 0.0D0
      Subinc_tminc = 0.0D0
      Subinc_tmaxc = 0.0D0
      Subinc_tavgc = 0.0D0
      Subinc_potet = 0.0D0
      Subinc_wb = 0.0D0
      Subinc_recharge = 0.0D0
      Subinc_deltastor = 0.0D0
      Subinc_szstor_frac = 0.0D0
      Subinc_capstor_frac = 0.0D0
      Subinc_gwflow = 0.0D0

      Laststor = Subinc_stor
      Subinc_stor = 0.0D0

      DO jj = 1, Active_hrus
        j = Hru_route_order(jj)
        ! k indicates which HRU is in which subbasin
        k = Hru_subbasin(j)
        IF ( k>0 ) THEN
          harea = Hru_area_dble(j)
          IF ( Hru_type(j)/=LAKE ) THEN
            srq = DBLE(Sroff(j))*harea
            ssq = DBLE(Ssres_flow(j))*harea
            soilstor = DBLE(Soil_moist(j)*Hru_frac_perv(j) + Ssres_stor(j))*harea
            snowstor = Pkwater_equiv(j)*harea
            landstor = DBLE(Hru_intcpstor(j)+Hru_impervstor(j))*harea
            IF ( Dprst_flag==ACTIVE ) landstor = landstor + Dprst_stor_hru(j)*harea
          ELSE
            soilstor = 0.0D0
            snowstor = 0.0D0
            landstor = 0.0D0
            ! wrong if multiple HRUs for any lake
            IF ( Lake_route_flag==ACTIVE ) THEN
              landstor = Lake_vol(Lake_hru_id(j))*12.0D0
              srq = Lake_outcfs(Lake_hru_id(j))*Cfs2inches
              ssq = 0.0D0
            ELSEIF ( Cascade_flag>CASCADE_OFF ) THEN
              srq = Hortonian_lakes(j)*harea
              ssq = Lakein_sz(j)*harea
            ELSE
              srq = 0.0D0
              ssq = 0.0D0
            ENDIF
          ENDIF
          Qsub(k) = Qsub(k) + srq + ssq
          Subinc_interflow(k) = Subinc_interflow(k) + ssq
          Subinc_sroff(k) = Subinc_sroff(k) + srq
          Subinc_precip(k) = Subinc_precip(k) + DBLE(Hru_ppt(j))*harea
          Subinc_rain(k) = Subinc_rain(k) + DBLE(Hru_rain(j))*harea
          Subinc_snow(k) = Subinc_snow(k) + DBLE(Hru_snow(j))*harea
          Subinc_actet(k) = Subinc_actet(k) + DBLE(Hru_actet(j))*harea
          Subinc_snowmelt(k) = Subinc_snowmelt(k) + DBLE(Snowmelt(j))*harea
          Subinc_pkweqv(k) = Subinc_pkweqv(k) + Pkwater_equiv(j)*harea
          Subinc_snowcov(k) = Subinc_snowcov(k) + DBLE(Snowcov_area(j))*harea
          Subinc_potet(k) = Subinc_potet(k) + DBLE(Potet(j))*harea
          Subinc_swrad(k) = Subinc_swrad(k) + DBLE(Swrad(j))*harea
          Subinc_tminc(k) = Subinc_tminc(k) + DBLE(Tminc(j))*harea
          Subinc_tmaxc(k) = Subinc_tmaxc(k) + DBLE(Tmaxc(j))*harea
          Subinc_tavgc(k) = Subinc_tavgc(k) + DBLE(Tavgc(j))*harea
          Subinc_recharge(k) = Subinc_recharge(k) + Recharge(j)*harea
          IF ( Soil_zone_max(j)>0.0 ) Subinc_szstor_frac(k) = Subinc_szstor_frac(k) + Soil_moist_tot(j)/Soil_zone_max(j)*harea
          IF ( Soil_moist_max(j)>0.0 ) Subinc_capstor_frac(k) = Subinc_capstor_frac(k) + Soil_moist(j)/Soil_moist_max(j)*harea
          Subinc_stor(k) = Subinc_stor(k) + soilstor + snowstor + landstor
          gwq = DBLE(Gwres_flow(j))*harea
          Qsub(k) = Qsub(k) + gwq
          Subinc_gwflow(k) = Subinc_gwflow(k) + gwq
          Subinc_stor(k) = Subinc_stor(k) + Gwres_stor(j)*harea
        ENDIF
      ENDDO

      !convert first as subbasins don't have to be in order
      DO j = 1, Nsub
        subarea = Sub_area(j)
        Sub_inq(j) = Qsub(j)*Cfs_conv
        dmy = Subinc_interflow(j)/subarea
        dmy2 = Subinc_sroff(j)/subarea
        Subinc_interflow(j) = Subinc_interflow(j)*Cfs_conv
        Subinc_sroff(j) = Subinc_sroff(j)*Cfs_conv
        Subinc_precip(j) = Subinc_precip(j)/subarea
        Subinc_rain(j) = Subinc_rain(j)/subarea
        Subinc_snow(j) = Subinc_snow(j)/subarea
        Subinc_actet(j) = Subinc_actet(j)/subarea
        Subinc_snowmelt(j) = Subinc_snowmelt(j)/subarea
        Subinc_pkweqv(j) = Subinc_pkweqv(j)/subarea
        Subinc_snowcov(j) = Subinc_snowcov(j)/subarea
        Subinc_potet(j) = Subinc_potet(j)/subarea
        Subinc_swrad(j) = Subinc_swrad(j)/subarea
        Subinc_tminc(j) = Subinc_tminc(j)/subarea
        Subinc_tmaxc(j) = Subinc_tmaxc(j)/subarea
        Subinc_tavgc(j) = Subinc_tavgc(j)/subarea
        Subinc_stor(j) = Subinc_stor(j)/subarea
        Subinc_recharge(j) = Subinc_recharge(j)/subarea
        Subinc_szstor_frac(j) = Subinc_szstor_frac(j)/subarea
        Subinc_capstor_frac(j) = Subinc_capstor_frac(j)/subarea
        Subinc_deltastor(j) = Laststor(j) - Subinc_stor(j)
        dmy1 = Subinc_gwflow(j)/subarea
        Subinc_gwflow(j) = Subinc_gwflow(j)*Cfs_conv
        ! water balance off if lake or muskingum routing
        Subinc_wb(j) = Subinc_precip(j) - Subinc_actet(j) - &
     &                 dmy - dmy1 - dmy2 + Subinc_deltastor(j)
      ENDDO

      !get cumulative subbasin flows
      Sub_gwflow = Subinc_gwflow
      DO j = 1, Nsub
        Sub_sroff(j) = Subinc_sroff(j)
        Sub_interflow(j) = Subinc_interflow(j)
        DO k = 1, Nsub
          IF ( Tree(j,k)/=0 ) THEN
            Sub_gwflow(j) = Sub_gwflow(j) + Subinc_gwflow(k)
            Sub_sroff(j) = Sub_sroff(j) + Subinc_sroff(k)
            Sub_interflow(j) = Sub_interflow(j) + Subinc_interflow(k)
          ENDIF
        ENDDO
      ENDDO

      DO j = 1, Nsub
        Sub_cfs(j) = Sub_inq(j)
        DO k = 1, Nsub
          IF ( Tree(j,k)/=0 ) Sub_cfs(j) = Sub_cfs(j) + Sub_inq(k)
        ENDDO
        Sub_cms(j) = Sub_cfs(j)*CFS2CMS_CONV
      ENDDO

      END FUNCTION subrun
