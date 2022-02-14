!***********************************************************************
! Defines shared watershed and HRU physical parameters and variables
!***********************************************************************
      MODULE PRMS_BASIN
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Basin Definition'
      character(len=*), parameter :: MODNAME = 'basin'
      character(len=*), parameter :: Version_basin = '2021-08-18'
      INTEGER, SAVE :: Numlake_hrus, Active_hrus, Active_gwrs, Numlakes_check
      INTEGER, SAVE :: Hemisphere, Dprst_clos_flag, Dprst_open_flag
      DOUBLE PRECISION, SAVE :: Land_area, Water_area
      DOUBLE PRECISION, SAVE :: Basin_area_inv, Basin_lat, Totarea, Active_area
      REAL, SAVE, ALLOCATABLE :: Hru_elev_meters(:) !, Hru_elev_feet(:)
      REAL, SAVE, ALLOCATABLE :: Dprst_frac_clos(:)
      INTEGER, SAVE, ALLOCATABLE :: Gwr_type(:), Hru_route_order(:), Gwr_route_order(:)
      INTEGER, SAVE :: Weir_gate_flag, Puls_lin_flag
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_area_dble(:), Lake_area(:)
!   Declared Variables
      REAL, SAVE, ALLOCATABLE :: Hru_frac_perv(:)
      REAL, SAVE, ALLOCATABLE :: Dprst_area_max(:)
      REAL, SAVE, ALLOCATABLE :: Hru_perv(:), Hru_imperv(:)
      REAL, SAVE, ALLOCATABLE :: Dprst_area_open_max(:), Dprst_area_clos_max(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_storage(:)
      REAL, SAVE, ALLOCATABLE :: Hru_elev_ts(:)
      DOUBLE PRECISION, SAVE :: Basin_gl_cfs, Basin_gl_ice_cfs
!   Declared Parameters
      INTEGER, SAVE :: Elev_units
      INTEGER, SAVE, ALLOCATABLE :: Hru_type(:), Cov_type(:)
      INTEGER, SAVE, ALLOCATABLE :: Lake_hru_id(:), Lake_type(:) !not needed if no lakes
      REAL, SAVE, ALLOCATABLE :: Hru_area(:), Hru_percent_imperv(:), Hru_elev(:), Hru_lat(:)
      REAL, SAVE, ALLOCATABLE :: Covden_sum(:), Covden_win(:)
      REAL, SAVE, ALLOCATABLE :: Dprst_frac_open(:), Dprst_area(:), Dprst_frac(:)
      END MODULE PRMS_BASIN

!***********************************************************************
!     Main basin routine
!***********************************************************************
      INTEGER FUNCTION basin()
      USE PRMS_CONSTANTS, ONLY: DECL, INIT
      USE PRMS_MODULE, ONLY: Process_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: basdecl, basinit
!***********************************************************************
      basin = 0

      IF ( Process_flag==DECL ) THEN
        basin = basdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        basin = basinit()
      ENDIF

      END FUNCTION basin

!***********************************************************************
!     basdecl - set up parameters
!   Declared Parameters
!     print_debug, hru_area, hru_percent_imperv, hru_type, hru_elev,
!     cov_type, hru_lat, dprst_frac_open, dprst_frac, basin_area
!     lake_hru_id
!***********************************************************************
      INTEGER FUNCTION basdecl()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, DOCUMENTATION
      USE PRMS_MODULE, ONLY: Nhru, Nlake, Model, Dprst_flag, Lake_route_flag, &
     &    PRMS4_flag, GSFLOW_flag, Glacier_flag
      USE PRMS_BASIN
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      basdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_basin)

! Declared Variables
      ALLOCATE ( Hru_elev_ts(Nhru) )
      IF ( Glacier_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'hru_elev_ts', 'nhru', Nhru, 'real', &
     &       'HRU elevation for timestep, which can change for glaciers', &
     &       'elev_units', Hru_elev_ts)/=0 ) CALL read_error(3, 'hru_elev_ts')

        IF ( declvar(MODNAME, 'basin_gl_ice_cfs', 'one', 1, 'double', &
     &       'Basin glacier ice (firn) melt leaving the basin through the stream network', &
     &       'cfs', Basin_gl_ice_cfs)/=0 ) CALL read_error(3, 'basin_gl_ice_cfs')

        IF ( declvar(MODNAME, 'basin_gl_cfs', 'one', 1, 'double', &
     &       'Basin glacier surface melt (rain, snow, ice) leaving the basin through the stream network', &
     &       'cfs', Basin_gl_cfs)/=0 ) CALL read_error(3, 'basin_gl_cfs')
      ENDIF

      ALLOCATE ( Hru_imperv(Nhru) )
      IF ( declvar(MODNAME, 'hru_imperv', 'nhru', Nhru, 'real', &
     &     'Area of HRU that is impervious', &
     &     'acres', Hru_imperv)/=0 ) CALL read_error(3, 'hru_imperv')

      ALLOCATE ( Hru_perv(Nhru) )
      IF ( declvar(MODNAME, 'hru_perv', 'nhru', Nhru, 'real', &
     &     'Area of HRU that is pervious', &
     &     'acres', Hru_perv)/=0 ) CALL read_error(3, 'hru_perv')

      ALLOCATE ( Hru_frac_perv(Nhru) )
     ! IF ( declvar(MODNAME, 'hru_frac_perv', 'nhru', Nhru, 'real', &
     !&     'Fraction of HRU that is pervious', &
     !&     'decimal fraction', Hru_frac_perv)/=0 ) CALL read_error(3, 'hru_frac_perv')

      ALLOCATE ( Hru_storage(Nhru) )
      IF ( declvar(MODNAME, 'hru_storage', 'nhru', Nhru, 'double', &
     &     'Storage for each HRU', &
     &     'inches', Hru_storage)/=0 ) CALL read_error(3, 'hru_storage')

      IF ( Dprst_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        ALLOCATE ( Dprst_area_max(Nhru) )
        IF ( declvar(MODNAME, 'dprst_area_max', 'nhru', Nhru, 'real', &
     &       'Aggregate sum of surface-depression storage areas of each HRU', &
     &       'acres', Dprst_area_max)/=0 ) CALL read_error(1, 'dprst_area_max')

        ALLOCATE ( Dprst_area_open_max(Nhru) )
        IF ( declvar(MODNAME, 'dprst_area_open_max', 'nhru', Nhru, 'real', &
     &       'Aggregate sum of open surface-depression storage areas of each HRU', &
     &       'acres', Dprst_area_open_max)/=0 ) CALL read_error(1, 'dprst_area_open_max')

        ALLOCATE ( Dprst_area_clos_max(Nhru) )
        IF ( declvar(MODNAME, 'dprst_area_clos_max', 'nhru', Nhru, 'real', &
     &       'Aggregate sum of closed surface-depression storage areas of each HRU', &
     &       'acres', Dprst_area_clos_max)/=0 ) CALL read_error(1, 'dprst_area_clos_max')

        ALLOCATE ( Dprst_frac(Nhru) )
        IF ( PRMS4_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Dprst_area(Nhru) )
          IF ( declparam(MODNAME, 'dprst_area', 'nhru', 'real', &
     &         '0.0', '0.0', '1.0E9', &
     &         'Aggregate sum of surface-depression storage areas of each HRU', &
     &         'Aggregate sum of surface-depression storage areas of each HRU', &
     &         'acres')/=0 ) CALL read_error(1, 'dprst_area')
          IF ( declparam(MODNAME, 'dprst_frac_hru', 'nhru', 'real', &
     &         '-1.0', '-1.0', '0.999', &
     &         'Fraction of each HRU area that has surface depressions', &
     &         'Fraction of each HRU area that has surface depressions', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'dprst_frac_hru')
        ENDIF
        IF ( PRMS4_flag==OFF .OR. Model==DOCUMENTATION ) THEN
          IF ( declparam(MODNAME, 'dprst_frac', 'nhru', 'real', &
     &         '0.0', '0.0', '0.999', &
     &         'Fraction of each HRU area that has surface depressions', &
     &         'Fraction of each HRU area that has surface depressions', &
     &         'decimal fraction')/=0 ) CALL read_error(1, 'dprst_frac')
        ENDIF

        ALLOCATE ( Dprst_frac_open(Nhru), Dprst_frac_clos(Nhru) )
        IF ( declparam(MODNAME, 'dprst_frac_open', 'nhru', 'real', &
     &       '1.0', '0.0', '1.0', &
     &       'Fraction of open surface-depression storage area within'// &
     &       ' an HRU that can generate surface runoff as a function of storage volume', &
     &       'Fraction of open surface-depression storage area within'// &
     &       ' an HRU that can generate surface runoff as a function of storage volume', &
     &       'decimal fraction')/=0 ) CALL read_error(1, 'dprst_frac_open')
      ENDIF

      ! local arrays
      ALLOCATE ( Hru_route_order(Nhru) )
! gwflow inactive for GSFLOW mode so arrays not allocated
! when GSFLOW can run in multi-mode will need these arrays
      IF ( GSFLOW_flag==OFF ) ALLOCATE ( Gwr_route_order(Nhru), Gwr_type(Nhru) )
!      ALLOCATE ( Hru_elev_feet(Nhru) )
      ALLOCATE ( Hru_elev_meters(Nhru) )

      ! Declared Parameters
      ALLOCATE ( Hru_area(Nhru), Hru_area_dble(Nhru) )
      IF ( declparam(MODNAME, 'hru_area', 'nhru', 'real', &
     &     '1.0', '0.0001', '1.0E9', &
     &     'HRU area', 'Area of each HRU', &
     &     'acres')/=0 ) CALL read_error(1, 'hru_area')

      IF ( declparam(MODNAME, 'elev_units', 'one', 'integer', &
     &     '0', '0', '1', &
     &     'Elevation units flag', &
     &     'Flag to indicate the units of the elevation values (0=feet; 1=meters)', &
     &     'none')/=0 ) CALL read_error(1, 'elev_units')

      ALLOCATE ( Hru_elev(Nhru) )
      IF ( declparam(MODNAME, 'hru_elev', 'nhru', 'real', &
     &     '0.0', '-1000.0', '30000.0', &
     &     'HRU mean elevation', 'Mean elevation for each HRU', &
     &     'elev_units')/=0 ) CALL read_error(1, 'hru_elev')

      ALLOCATE ( Hru_lat(Nhru) )
      IF ( declparam(MODNAME, 'hru_lat', 'nhru', 'real', &
     &     '40.0', '-90.0', '90.0', &
     &     'HRU latitude', 'Latitude of each HRU', &
     &     'degrees North')/=0 ) CALL read_error(1, 'hru_lat')

      ALLOCATE ( Hru_percent_imperv(Nhru) )
      IF ( declparam(MODNAME, 'hru_percent_imperv', 'nhru', 'real', &
     &     '0.0', '0.0', '0.999', &
     &     'HRU percent impervious', 'Fraction of each HRU area that is impervious', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'hru_percent_imperv')

      ALLOCATE ( Hru_type(Nhru) )
      IF ( declparam(MODNAME, 'hru_type', 'nhru', 'integer', &
     &     '1', '0', '4', &
     &     'HRU type', 'Type of each HRU (0=inactive; 1=land; 2=lake; 3=swale; 4=glacier)', &
     &     'none')/=0 ) CALL read_error(1, 'hru_type')

      ALLOCATE ( Cov_type(Nhru) )
      IF ( declparam(MODNAME, 'cov_type', 'nhru', 'integer', &
     &     '3', '0', '4', &
     &     'Cover type designation for each HRU', &
     &     'Vegetation cover type for each HRU (0=bare soil;'// &
     &     ' 1=grasses; 2=shrubs; 3=trees; 4=coniferous)', &
     &     'none')/=0 ) CALL read_error(1, 'cov_type')

      ALLOCATE ( Covden_sum(Nhru) )
      IF ( declparam(MODNAME, 'covden_sum', 'nhru', 'real', &
     &     '0.5', '0.0', '1.0', &
     &     'Summer vegetation cover density for major vegetation type', &
     &     'Summer vegetation cover density for the major vegetation type in each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'covden_sum')

      ALLOCATE ( Covden_win(Nhru) )
      IF ( declparam(MODNAME, 'covden_win', 'nhru', 'real', &
     &     '0.5', '0.0', '1.0', &
     &     'Winter vegetation cover density for major vegetation type', &
     &     'Winter vegetation cover density for the major vegetation type in each HRU', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'covden_win')

      ALLOCATE ( Lake_hru_id(Nhru) )
      IF ( Nlake>0 .OR. Model==DOCUMENTATION ) THEN
        ! Local array
        ALLOCATE ( Lake_area(Nlake) )
        ! parameters
        IF ( declparam(MODNAME, 'lake_hru_id', 'nhru', 'integer', &
     &       '0', 'bounded', 'nlake', &
     &       'Identification number of the lake associated with an HRU', &
     &       'Identification number of the lake associated with an HRU;'// &
     &       ' more than one HRU can be associated with each lake', &
     &       'none')/=0 ) CALL read_error(1, 'lake_hru_id')
        IF ( (Lake_route_flag==ACTIVE .AND. GSFLOW_flag==OFF ) .OR. Model==DOCUMENTATION ) THEN
          ALLOCATE ( Lake_type(Nlake) )
          IF ( declparam(MODNAME, 'lake_type', 'nlake', 'integer', &
     &         '1', '1', '6', &
     &         'Type of lake routing method', &
     &         'Type of lake routing method (1=Puls routing; 2=linear routing; 3=flow through;'// &
     &         ' 4=broad crested weir; 5=gate opening; 6=measured flow)', &
     &         'none')/=0 ) CALL read_error(1, 'lake_type')
        ENDIF
      ELSE
        ALLOCATE ( Lake_area(1) )
      ENDIF

      END FUNCTION basdecl

!**********************************************************************
!     basinit - check for validity of basin parameters
!               and compute reservoir areas
!**********************************************************************
      INTEGER FUNCTION basinit()
      USE PRMS_CONSTANTS, ONLY: DEBUG_less, ACTIVE, OFF, &
     &    INACTIVE, LAKE, SWALE, FEET, ERROR_basin, DEBUG_minimum, &
     &    NORTHERN, SOUTHERN, FEET2METERS, DNEARZERO
      USE PRMS_MODULE, ONLY: Nhru, Nlake, Print_debug, &
     &    Dprst_flag, Lake_route_flag, PRMS4_flag, GSFLOW_flag, Frozen_flag, PRMS_VERSION, &
     &    Starttime, Endtime, Parameter_check_flag
      USE PRMS_BASIN
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: write_outfile, checkdim_bounded_limits
      INTRINSIC :: DBLE
! Local Variables
      CHARACTER(LEN=69) :: buffer
      INTEGER :: i, j, dprst_frac_flag, lakeid
      REAL :: harea, perv_area
      DOUBLE PRECISION :: basin_imperv, basin_perv, basin_dprst, harea_dble
!**********************************************************************
      basinit = 0

      IF ( getparam(MODNAME, 'hru_area', Nhru, 'real', Hru_area)/=0 ) CALL read_error(2, 'hru_area')
      IF ( getparam(MODNAME, 'hru_elev', Nhru, 'real', Hru_elev)/=0 ) CALL read_error(2, 'hru_elev')
      Hru_elev_ts = Hru_elev
      IF ( getparam(MODNAME, 'hru_lat', Nhru, 'real', Hru_lat)/=0 ) CALL read_error(2, 'hru_lat')
      IF ( getparam(MODNAME, 'hru_type', Nhru, 'integer', Hru_type)/=0 ) CALL read_error(2, 'hru_type')
      IF ( getparam(MODNAME, 'cov_type', Nhru, 'integer', Cov_type)/=0 ) CALL read_error(2, 'cov_type')
      IF ( getparam(MODNAME, 'covden_sum', Nhru, 'real', Covden_sum)/=0 ) CALL read_error(2, 'covden_sum')
      IF ( getparam(MODNAME, 'covden_win', Nhru, 'real', Covden_win)/=0 ) CALL read_error(2, 'covden_win')
      IF ( getparam(MODNAME, 'elev_units', 1, 'integer', Elev_units)/=0 ) CALL read_error(2, 'elev_units')
      IF ( getparam(MODNAME, 'hru_percent_imperv', Nhru, 'real', Hru_percent_imperv)/=0 ) CALL read_error(2, 'hru_percent_imperv')

      dprst_frac_flag = 0
      IF ( Dprst_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'dprst_frac_open', Nhru, 'real', Dprst_frac_open)/=0 ) CALL read_error(2, 'dprst_frac_open')
        IF ( PRMS4_flag==ACTIVE ) THEN
          IF ( getparam(MODNAME, 'dprst_frac_hru', Nhru, 'real', Dprst_frac)/=0 ) CALL read_error(2, 'dprst_frac_hru')
          IF ( Dprst_frac(1)>-1.0 ) THEN
            IF ( Print_debug>DEBUG_less ) PRINT *, 'Using dprst_frac_hru instead of dprst_area'
            dprst_frac_flag = 1
          ELSE
            IF ( getparam(MODNAME, 'dprst_area', Nhru, 'real', Dprst_area)/=0 ) CALL read_error(2, 'dprst_area')
          ENDIF
        ELSE
          IF ( getparam(MODNAME, 'dprst_frac', Nhru, 'real', Dprst_frac)/=0 ) CALL read_error(2, 'Dprst_frac')
        ENDIF
      ENDIF

      Weir_gate_flag = OFF
      Puls_lin_flag = OFF
      Water_area = 0.0D0
      Lake_area = 0.0D0
      Numlakes_check = 0
      Numlake_hrus = 0
      IF ( Nlake>0 ) THEN
        IF ( getparam(MODNAME, 'lake_hru_id', Nhru, 'integer', Lake_hru_id)/=0 ) CALL read_error(1, 'lake_hru_id')
        IF ( Parameter_check_flag>0 ) CALL checkdim_bounded_limits('lake_hru_id', 'nlake', Lake_hru_id, Nhru, 0, Nlake, basinit)
        IF ( Lake_route_flag==ACTIVE ) THEN ! Lake_route_flag set to 0 for GSFLOW mode and if muskingum_lake and nlake = 1
          IF ( getparam(MODNAME, 'lake_type', Nlake, 'integer', Lake_type)/=0 ) CALL read_error(2, 'lake_type')
          DO i = 1, Nlake
            IF ( Lake_type(i)==4 .OR. Lake_type(i)==5 ) THEN
              Weir_gate_flag = ACTIVE
            ELSEIF ( Lake_type(i)==1 .OR. Lake_type(i)==2 ) THEN
              Puls_lin_flag = ACTIVE
            ENDIF
          ENDDO
        ENDIF
      ELSE
        Lake_hru_id = 0
      ENDIF

      Basin_gl_cfs = 0.0D0
      Basin_gl_ice_cfs = 0.0D0
      IF ( Dprst_flag==ACTIVE ) THEN
        Dprst_frac_clos = 0.0
        Dprst_area_open_max = 0.0
        Dprst_area_clos_max = 0.0
        Dprst_area_max = 0.0
      ENDIF
      Dprst_clos_flag = OFF
      Dprst_open_flag = OFF
      basin_perv = 0.0D0
      basin_imperv = 0.0D0
      basin_dprst = 0.0D0
      Totarea = 0.0D0
      Land_area = 0.0D0
      Active_area = 0.0D0
      Basin_lat = 0.0D0
      Hru_frac_perv = 0.0
      Hru_imperv = 0.0
      Hru_perv = 0.0
      Hru_storage = 0.0D0
      Hru_route_order = 0
      j = 0
      DO i = 1, Nhru
        harea = Hru_area(i)
        harea_dble = DBLE( harea )
        Hru_area_dble(i) = harea_dble
        Totarea = Totarea + harea_dble
        perv_area = harea

        IF ( Hru_type(i)==INACTIVE ) CYCLE
! ????????? need to fix for lakes with multiple HRUs and PRMS lake routing ????????
        IF ( Hru_type(i)==LAKE ) THEN
          Numlake_hrus = Numlake_hrus + 1
          Water_area = Water_area + harea_dble
          lakeid = Lake_hru_id(i)
          IF ( lakeid>0 ) THEN
            Lake_area(lakeid) = Lake_area(lakeid) + harea_dble
            IF ( lakeid>Numlakes_check ) Numlakes_check = lakeid
          ELSE
            PRINT *, 'ERROR, hru_type = 2 for HRU:', i, ' and lake_hru_id = 0'
            basinit = 1
            CYCLE
          ENDIF
          IF ( Nlake==0 ) THEN
            PRINT *, 'ERROR, hru_type = 2 for HRU:', i, ' and dimension nlake = 0'
            basinit = 1
            CYCLE
          ENDIF
        ELSE
          Land_area = Land_area + harea_dble ! swale or land or glacier
          IF ( Lake_hru_id(i)>0 ) THEN
            PRINT *, 'ERROR, HRU:', i, ' specifed to be a lake by lake_hru_id but hru_type not equal 2'
            basinit = 1
            CYCLE
          ENDIF
          IF ( Frozen_flag==ACTIVE ) THEN
            IF ( Hru_type(i)==SWALE ) THEN
              PRINT *, 'ERROR, a swale HRU cannot be frozen for CFGI, HRU:', i
              basinit = 1
              CYCLE
            ENDIF
          ENDIF

        ENDIF

        Basin_lat = Basin_lat + DBLE( Hru_lat(i)*harea )
        IF ( Elev_units==FEET ) THEN
!          Hru_elev_feet(i) = Hru_elev(i)
          Hru_elev_meters(i) = Hru_elev(i)*FEET2METERS
        ELSE
!          Hru_elev_feet(i) = Hru_elev(i)*METERS2FEET
          Hru_elev_meters(i) = Hru_elev(i)
        ENDIF
        j = j + 1
        Hru_route_order(j) = i

        IF ( Hru_type(i)==LAKE ) CYCLE

        IF ( Hru_percent_imperv(i)>0.0 ) THEN
          Hru_imperv(i) = Hru_percent_imperv(i)*harea
          basin_imperv = basin_imperv + DBLE( Hru_imperv(i) )
          perv_area = perv_area - Hru_imperv(i)
        ENDIF

        IF ( Dprst_flag==ACTIVE ) THEN
          IF ( dprst_frac_flag==1 .OR. PRMS4_flag==OFF ) THEN
            Dprst_area_max(i) = Dprst_frac(i)*harea
          ELSE
            Dprst_area_max(i) = Dprst_area(i)
            Dprst_frac(i) = Dprst_area_max(i)/harea
          ENDIF
          IF ( Dprst_area_max(i)>0.0 ) THEN
            Dprst_area_open_max(i) = Dprst_area_max(i)*Dprst_frac_open(i)
            Dprst_frac_clos(i) = 1.0 - Dprst_frac_open(i)
            Dprst_area_clos_max(i) = Dprst_area_max(i) - Dprst_area_open_max(i)
            basin_dprst = basin_dprst + DBLE( Dprst_area_max(i) )
            IF ( Dprst_area_clos_max(i)>0.0 ) Dprst_clos_flag = ACTIVE
            IF ( Dprst_area_open_max(i)>0.0 ) Dprst_open_flag = ACTIVE
            perv_area = perv_area - Dprst_area_max(i)
          ENDIF
        ENDIF

        Hru_perv(i) = perv_area
        Hru_frac_perv(i) = perv_area/harea
        IF ( Hru_frac_perv(i)<0.00099 ) THEN
          PRINT *, 'ERROR, pervious fraction must be >= 0.001 for HRU:', i
          PRINT *, '       pervious portion is HRU fraction - impervious fraction - depression fraction'
          PRINT *, '       pervious fraction:', Hru_frac_perv(i)
          PRINT *, '       impervious fraction:', Hru_percent_imperv(i)
          IF ( Dprst_flag==ACTIVE ) PRINT *, '       depression storage fraction:', Dprst_frac(i)
          basinit = 1
        ENDIF
        basin_perv = basin_perv + DBLE( Hru_perv(i) )
      ENDDO
      IF ( Dprst_flag==ACTIVE .AND. PRMS4_flag==ACTIVE ) DEALLOCATE ( Dprst_area )

      IF ( Nlake>0 ) THEN
!        IF ( Numlake_hrus/=Nlake_hrus ) THEN
!          PRINT *, 'ERROR, number of lake HRUs specified in hru_type'
!          PRINT *, 'does not equal dimension nlake:', Nlake, ', number of lake HRUs:', Numlake_hrus
!          PRINT *, 'For PRMS lake routing each lake must be a single HRU'
!          basinit = 1
!        ENDIF
        IF ( Numlakes_check/=Nlake ) THEN
          PRINT *, 'ERROR, number of lakes specified in lake_hru_id'
          PRINT *, 'does not equal dimension nlake:', Nlake, ', number of lakes:', Numlakes_check
          PRINT *, 'For PRMS lake routing each lake must be a single HRU'
          basinit = 1
        ENDIF
        DO i = 1, Numlakes_check
          IF ( Lake_area(i)<DNEARZERO ) THEN
            PRINT *, 'ERROR, Lake:', i, ' has 0 area, thus no value of lake_hru_id is associated with the lake'
            basinit = 1
          ENDIF
        ENDDO
      ENDIF

      IF ( basinit==1 ) ERROR STOP ERROR_basin

      Active_hrus = j
      Active_area = Land_area + Water_area

      Active_gwrs = Active_hrus
      IF ( GSFLOW_flag==OFF ) THEN
        Gwr_type = Hru_type
        Gwr_route_order = Hru_route_order
      ENDIF

      Basin_area_inv = 1.0D0/Active_area
      Basin_lat = Basin_lat*Basin_area_inv
      ! used in solrad modules to winter/summer radiation adjustment
      IF ( Basin_lat>0.0D0 ) THEN
        Hemisphere = NORTHERN
      ELSE
        Hemisphere = SOUTHERN
      ENDIF

      basin_perv = basin_perv*Basin_area_inv
      basin_imperv = basin_imperv*Basin_area_inv

      IF ( Print_debug==2 ) THEN
        PRINT *, ' HRU     Area'
        PRINT ('(I7, F14.5)'), (i, Hru_area(i), i=1, Nhru)
        PRINT *,                           'Model domain area           = ', Totarea
        PRINT *,                           'Active basin area           = ', Active_area
        IF ( Water_area>0.0D0 ) PRINT *,   'Lake area                   = ', Water_area
        PRINT *,                           'Fraction impervious         = ', basin_imperv
        PRINT *,                           'Fraction pervious           = ', basin_perv
        IF ( Water_area>0.0D0 ) PRINT *,   'Fraction lakes              = ', Water_area*Basin_area_inv
        IF ( Dprst_flag==ACTIVE ) PRINT *, 'Fraction storage area       = ', basin_dprst*Basin_area_inv
        PRINT *, ' '
      ENDIF

!     print out start and end times
      IF ( Print_debug>DEBUG_minimum ) THEN
        !CALL write_outfile(' Surface Water and Energy Budgets Simulated by '//PRMS_VERSION)
        CALL write_outfile(' ')
        WRITE (buffer, 9002) 'Start time: ', Starttime
        CALL write_outfile(buffer(:31))
        WRITE (buffer, 9002) 'End time:   ', Endtime
        CALL write_outfile(buffer(:31))
        CALL write_outfile(' ')
        WRITE (buffer, 9003) 'Model domain area:   ', Totarea, '    Active basin area:', Active_area
        CALL write_outfile(buffer)
        WRITE (buffer, 9004)   'Fraction impervious:', basin_imperv, '    Fraction pervious: ', basin_perv
        CALL write_outfile(buffer)
        IF ( Water_area>0.0D0 ) THEN
          WRITE (buffer, 9004) 'Lake area:          ', Water_area,   '    Fraction lakes:    ', Water_area*Basin_area_inv
          CALL write_outfile(buffer)
        ENDIF
        IF ( Dprst_flag==ACTIVE ) THEN
          WRITE (buffer, 9005) 'DPRST area:          ', basin_dprst, '    Fraction DPRST:   ', basin_dprst*Basin_area_inv
          CALL write_outfile(buffer)
        ENDIF
        CALL write_outfile(' ')
      ENDIF

 9002 FORMAT (A, I4.2, 2('/', I2.2), I3.2, 2(':', I2.2))
 9003 FORMAT (2(A,F13.2))
 9004 FORMAT (A, F14.4, A, F12.4)
 9005 FORMAT (A, F13.2, A, F13.4)

      END FUNCTION basinit
