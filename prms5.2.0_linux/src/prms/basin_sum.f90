!***********************************************************************
! Computes daily, monthly, yearly, and total flow summaries of volumes
! and flows for all HRUs
!***********************************************************************
      MODULE PRMS_BASINSUM
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, DOCUMENTATION, &
     &    strmflow_muskingum_module, strmflow_muskingum_lake_module, strmflow_muskingum_mann_module
      USE PRMS_MODULE, ONLY: Nhru, Nobs, Model, Init_vars_from_file, &
     &    Print_debug, End_year, Strmflow_flag, Glacier_flag
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Output Summary'
      character(len=9), parameter :: MODNAME = 'basin_sum'
      character(len=*), parameter :: Version_basin_sum = '2020-12-02'

      INTEGER, SAVE :: BALUNT, Totdays
      INTEGER, SAVE :: Header_prt, Endjday
      CHARACTER(LEN=32) :: Buffer32
      CHARACTER(LEN=40) :: Buffer40
      CHARACTER(LEN=48) :: Buffer48
      CHARACTER(LEN=80) :: Buffer80
      CHARACTER(LEN=120) :: Buffer120
      CHARACTER(LEN=151) :: Buffer151
      CHARACTER(LEN=151), PARAMETER :: DASHS = ' --------------------------------------------------------'// &
     &  '----------------------------------------------------------------------------------------------'
      CHARACTER(LEN=151), PARAMETER :: STARS = ' ********************************************************'// &
     &  '**********************************************************************************************'
      CHARACTER(LEN=151), PARAMETER :: EQULS = ' ========================================================'// &
     &  '=============================================================================================='
      LOGICAL, SAVE :: Dprt, Mprt, Yprt, Tprt
      DOUBLE PRECISION, SAVE :: Basin_swrad_yr, Basin_swrad_tot, Basin_swrad_mo
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Obs_runoff_mo, Obs_runoff_yr, Obs_runoff_tot
      DOUBLE PRECISION, SAVE :: Basin_cfs_mo, Basin_cfs_yr, Basin_cfs_tot
      DOUBLE PRECISION, SAVE :: Basin_net_ppt_yr, Basin_net_ppt_tot, Watbal_sum
      DOUBLE PRECISION, SAVE :: Basin_max_temp_yr, Basin_max_temp_tot
      DOUBLE PRECISION, SAVE :: Basin_min_temp_yr, Basin_min_temp_tot
      DOUBLE PRECISION, SAVE :: Basin_potet_yr, Basin_potet_tot
      DOUBLE PRECISION, SAVE :: Basin_actet_yr, Basin_actet_tot
      DOUBLE PRECISION, SAVE :: Basin_snowmelt_yr, Basin_snowmelt_tot
      DOUBLE PRECISION, SAVE :: Basin_gwflow_yr, Basin_gwflow_tot
      DOUBLE PRECISION, SAVE :: Basin_ssflow_yr, Basin_ssflow_tot
      DOUBLE PRECISION, SAVE :: Basin_sroff_yr, Basin_sroff_tot
      DOUBLE PRECISION, SAVE :: Basin_stflow_yr, Basin_stflow_tot
      DOUBLE PRECISION, SAVE :: Basin_ppt_yr, Basin_ppt_tot, Last_basin_stor
      DOUBLE PRECISION, SAVE :: Basin_intcp_evap_yr, Basin_intcp_evap_tot, Basin_lakeevap_yr
      DOUBLE PRECISION, SAVE :: Obsq_inches_yr, Obsq_inches_tot
      DOUBLE PRECISION, SAVE :: Basin_net_ppt_mo, Obsq_inches_mo
      DOUBLE PRECISION, SAVE :: Basin_max_temp_mo, Basin_min_temp_mo
      DOUBLE PRECISION, SAVE :: Basin_actet_mo
      DOUBLE PRECISION, SAVE :: Basin_snowmelt_mo, Basin_gwflow_mo
      DOUBLE PRECISION, SAVE :: Basin_sroff_mo, Basin_stflow_mo
      DOUBLE PRECISION, SAVE :: Basin_intcp_evap_mo, Basin_storage
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_et_yr(:)
      DOUBLE PRECISION, SAVE :: Basin_storvol, Basin_potet_mo
      DOUBLE PRECISION, SAVE :: Basin_ssflow_mo, Basin_ppt_mo
      DOUBLE PRECISION, SAVE :: Obsq_inches
      DOUBLE PRECISION, SAVE :: Basin_runoff_ratio, Basin_runoff_ratio_mo
      DOUBLE PRECISION, SAVE :: Basin_lakeevap_mo
!   Declared Parameters
      INTEGER, SAVE :: Print_type, Print_freq, Outlet_sta
      END MODULE PRMS_BASINSUM

!***********************************************************************
!     Main basin_sum routine
!***********************************************************************
      INTEGER FUNCTION basin_sum()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, INIT, CLEAN, ACTIVE, OFF, SAVE_INIT, READ_INIT
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: sumbdecl, sumbinit, sumbrun
      EXTERNAL :: basin_sum_restart
!***********************************************************************
      basin_sum = 0

      IF ( Process_flag==RUN ) THEN
        basin_sum = sumbrun()
      ELSEIF ( Process_flag==DECL ) THEN
        basin_sum = sumbdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        IF ( Init_vars_from_file>OFF ) CALL basin_sum_restart(READ_INIT)
        basin_sum = sumbinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL basin_sum_restart(SAVE_INIT)
      ENDIF

      END FUNCTION basin_sum

!***********************************************************************
!     sumbdecl - set up basin summary parameters
!   Declared Parameters
!     print_type, print_freq, outlet_sta
!***********************************************************************
      INTEGER FUNCTION sumbdecl()
      USE PRMS_BASINSUM
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      sumbdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_basin_sum)

      IF ( declvar(MODNAME, 'last_basin_stor', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in all water storage reservoirs from previous time step', &
     &     'inches', Last_basin_stor)/=0 ) CALL read_error(3, 'last_basin_stor')
      IF ( declvar(MODNAME, 'watbal_sum', 'one', 1, 'double', &
     &     'Water balance aggregate', &
     &     'inches', Watbal_sum)/=0 ) CALL read_error(3, 'watbal_sum')
      IF ( declvar(MODNAME, 'obs_runoff_mo', 'one', 1, 'double', &
     &     'Monthly measured streamflow at basin outlet', &
     &     'cfs', Obs_runoff_mo)/=0 ) CALL read_error(3, 'obs_runoff_mo')
      IF ( declvar(MODNAME, 'basin_cfs_mo', 'one', 1, 'double', &
     &     'Monthly total streamflow to stream network', &
     &     'cfs', Basin_cfs_mo)/=0 ) CALL read_error(3, 'basin_cfs_mo')
      IF ( declvar(MODNAME, 'obs_runoff_yr', 'one', 1, 'double', &
     &     'Yearly measured streamflow at basin outlet', &
     &     'cfs', Obs_runoff_yr)/=0 ) CALL read_error(3, 'obs_runoff_yr')
      IF ( declvar(MODNAME, 'basin_cfs_yr', 'one', 1, 'double', &
     &     'Yearly total streamflow to stream network', &
     &     'cfs', Basin_cfs_yr)/=0 ) CALL read_error(3, 'basin_cfs_yr')
      IF ( declvar(MODNAME, 'basin_net_ppt_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average net precipitation', &
     &     'inches', Basin_net_ppt_yr)/=0 ) CALL read_error(3, 'basin_net_ppt_yr')
      IF ( declvar(MODNAME, 'basin_max_temp_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average maximum temperature', &
     &     'temp_units', Basin_max_temp_yr)/=0) CALL read_error(3,'basin_max_temp_yr')
      IF ( declvar(MODNAME, 'basin_min_temp_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average minimum temperature', &
     &     'temp_units', Basin_min_temp_yr)/=0) CALL read_error(3,'basin_min_temp_yr')
      IF ( declvar(MODNAME, 'basin_potet_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average potential ET', &
     &     'temp_units', Basin_potet_yr)/=0 ) CALL read_error(3, 'basin_potet_yr')
      IF ( declvar(MODNAME, 'basin_actet_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average actual ET', &
     &     'inches', Basin_actet_yr)/=0 ) CALL read_error(3, 'basin_actet_yr')
      IF ( declvar(MODNAME, 'basin_snowmelt_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average snowmelt', &
     &     'inches', Basin_snowmelt_yr)/=0) CALL read_error(3,'basin_snowmelt_yr')
      IF ( declvar(MODNAME, 'basin_gwflow_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average groundwater discharge', &
     &     'inches', Basin_gwflow_yr)/=0 ) CALL read_error(3, 'basin_gwflow_yr')
      IF ( declvar(MODNAME, 'basin_ssflow_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average interflow', &
     &     'inches', Basin_ssflow_yr)/=0 ) CALL read_error(3, 'basin_ssflow_yr')
      IF ( declvar(MODNAME, 'basin_sroff_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average overland runoff', &
     &     'inches', Basin_sroff_yr)/=0 ) CALL read_error(3, 'basin_sroff_yr')
      IF ( declvar(MODNAME, 'basin_ppt_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average precipitation', &
     &     'inches', Basin_ppt_yr)/=0 ) CALL read_error(3, 'basin_ppt_yr')
      IF ( declvar(MODNAME, 'basin_stflow_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average streamflow', &
     &     'inches', Basin_stflow_yr)/=0 ) CALL read_error(3, 'basin_stflow_yr')
      IF ( declvar(MODNAME, 'obsq_inches_yr', 'one', 1, 'double', &
     &     'Yearly measured streamflow at specified outlet station', &
     &     'inches', Obsq_inches_yr)/=0 ) CALL read_error(3, 'obsq_inches_yr')
      IF ( declvar(MODNAME, 'basin_intcp_evap_yr', 'one', 1, 'double', &
     &     'Yearly basin area-weighted average canopy evaporation', &
     &     'inches', Basin_intcp_evap_yr)/=0 ) CALL read_error(3, 'basin_intcp_evap_yr')
      IF ( declvar(MODNAME, 'obs_runoff_tot', 'one', 1, 'double', &
     &     'Total simulation measured streamflow at basin outlet', &
     &     'cfs', Obs_runoff_tot)/=0 ) CALL read_error(3, 'obs_runoff_tot')
      IF ( declvar(MODNAME, 'basin_cfs_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average streamflow', &
     &     'inches', Basin_cfs_tot)/=0 ) CALL read_error(3, 'basin_cfs_tot')
      IF ( declvar(MODNAME, 'basin_ppt_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average precipitation', &
     &     'inches', Basin_ppt_tot)/=0 ) CALL read_error(3, 'basin_ppt_tot')
      IF ( declvar(MODNAME, 'basin_max_temp_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average maximum temperature', &
     &     'temp_units', Basin_max_temp_tot)/=0 ) CALL read_error(3, 'basin_max_temp_tot')
      IF ( declvar(MODNAME, 'basin_min_temp_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average minimum temperature', &
     &     'temp_units', Basin_min_temp_tot)/=0 ) CALL read_error(3, 'basin_min_temp_tot')
      IF ( declvar(MODNAME, 'basin_potet_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average potential ET', &
     &     'inches', Basin_potet_tot)/=0 ) CALL read_error(3, 'basin_potet_tot')
      IF ( declvar(MODNAME, 'basin_actet_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average actual ET', &
     &     'inches', Basin_actet_tot)/=0 ) CALL read_error(3, 'basin_actet_tot')
      IF ( declvar(MODNAME, 'basin_snowmelt_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average snowmelt', &
     &     'inches', Basin_snowmelt_tot)/=0 ) CALL read_error(3, 'basin_snowmelt_tot')
      IF ( declvar(MODNAME, 'basin_gwflow_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average groundwater discharge' , &
     &     'inches', Basin_gwflow_tot)/=0 ) CALL read_error(3, 'basin_gwflow_tot')
      IF ( declvar(MODNAME, 'basin_ssflow_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average interflow', &
     &     'inches', Basin_ssflow_tot)/=0 ) CALL read_error(3, 'basin_ssflow_tot')
      IF ( declvar(MODNAME, 'basin_sroff_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average overland flow', &
     &     'inches', Basin_sroff_tot)/=0 ) CALL read_error(3, 'basin_sroff_tot')
      IF ( declvar(MODNAME, 'basin_stflow_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average streamflow', &
     &     'inches', Basin_stflow_tot)/=0 ) CALL read_error(3, 'basin_stflow_tot')
      IF ( declvar(MODNAME, 'obsq_inches_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average measured streamflow at specified outlet station', &
     &     'inches', Obsq_inches_tot)/=0 ) CALL read_error(3, 'obsq_inches_tot')
      IF ( declvar(MODNAME, 'basin_intcp_evap_tot', 'one', 1, 'double', &
     &     'Total simulation basin area-weighted average canopy evaporation', &
     &     'inches', Basin_intcp_evap_tot)/=0 ) CALL read_error(3, 'basin_intcp_evap_tot')

! declare parameters
      IF ( Nobs>0 .OR. Model==DOCUMENTATION ) THEN
        IF ( declparam(MODNAME, 'outlet_sta', 'one', 'integer', &
     &       '0', 'bounded', 'nobs', &
     &       'Index of measurement station to use for basin outlet', &
     &       'Index of measured streamflow station corresponding to the basin outlet', &
     &       'none')/=0 ) CALL read_error(1, 'outlet_sta')
      ENDIF

      IF ( declparam(MODNAME, 'print_type', 'one', 'integer', &
     &     '1', '0', '2', &
     &     'Type of output written to output file', &
     &     'Flag to select the type of results written to the output'// &
     &     ' file (0=measured and simulated flow only;'// &
     &     ' 1=water balance table; 2=detailed output)', &
     &     'none')/=0 ) CALL read_error(1, 'print_type')

      IF ( declparam(MODNAME, 'print_freq', 'one', 'integer', &
     &     '3', '0', '15', &
     &     'Frequency for the output frequency', &
     &     'Flag to select the output frequency; for combinations,'// &
     &     ' add index numbers, e.g., daily plus yearly = 10;'// &
     &     ' yearly plus total = 3 (0=none; 1=run totals; 2=yearly;'// &
     &     ' 4=monthly; 8=daily; or additive combinations)', &
     &     'none')/=0 ) CALL read_error(1, 'print_freq')

      IF ( declvar(MODNAME, 'basin_intcp_evap_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average interception evaporation', &
     &     'inches', Basin_intcp_evap_mo)/=0 ) CALL read_error(3, 'basin_intcp_evap_mo')

      IF ( declvar(MODNAME, 'basin_storage', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in all water storage reservoirs', &
     &     'inches', Basin_storage)/=0 ) CALL read_error(3, 'basin_storage')

!******************basin_storage volume:
      IF ( declvar(MODNAME, 'basin_storvol', 'one', 1, 'double', &
     &     'Basin area-weighted average storage volume in all water storage reservoirs', &
     &     'acre-inches', Basin_storvol)/=0 ) CALL read_error(3, 'basin_storvol')

      IF ( declvar(MODNAME, 'obsq_inches', 'one', 1, 'double', &
     &     'Measured streamflow at specified outlet station', &
     &     'inches', Obsq_inches)/=0 ) CALL read_error(3, 'obsq_inches')

      IF ( declvar(MODNAME, 'basin_ppt_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average precipitation', &
     &     'inches', Basin_ppt_mo)/=0 ) CALL read_error(3, 'basin_ppt_mo')

      IF ( declvar(MODNAME, 'basin_net_ppt_mo', 'one', 1, 'double', &
     &     'Monthly area-weighted average net precipitation', &
     &     'inches', Basin_net_ppt_mo)/=0 ) CALL read_error(3, 'basin_net_ppt_mo')

      IF ( declvar(MODNAME, 'basin_max_temp_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average maximum air temperature', &
     &     'temp_units', Basin_max_temp_mo)/=0 ) CALL read_error(3, 'basin_max_temp_mo')

      IF ( declvar(MODNAME, 'basin_min_temp_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average minimum air temperature', &
     &     'temp_units', Basin_min_temp_mo)/=0 ) CALL read_error(3, 'basin_min_temp_mo')

      IF ( declvar(MODNAME, 'basin_potet_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average potential ET', &
     &     'inches', Basin_potet_mo)/=0 ) CALL read_error(3, 'basin_potet_mo')

      IF ( declvar(MODNAME, 'basin_actet_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average actual ET', &
     &     'inches', Basin_actet_mo)/=0 ) CALL read_error(3, 'basin_actet_mo')

      IF ( declvar(MODNAME, 'basin_snowmelt_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average snowmelt', &
     &     'inches', Basin_snowmelt_mo)/=0 ) CALL read_error(3, 'basin_snowmelt_mo')

      IF ( declvar(MODNAME, 'basin_gwflow_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average groundwater discharge', &
     &     'inches', Basin_gwflow_mo)/=0 ) CALL read_error(3, 'basin_gwflow_mo')

      IF ( declvar(MODNAME, 'basin_ssflow_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average interflow', &
     &     'inches', Basin_ssflow_mo)/=0 ) CALL read_error(3, 'basin_ssflow_mo')

      IF ( declvar(MODNAME, 'basin_sroff_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average surface runoff', &
     &     'inches', Basin_sroff_mo)/=0 ) CALL read_error(3, 'basin_sroff_mo')

      IF ( declvar(MODNAME, 'basin_stflow_mo', 'one', 1, 'double', &
     &     'Monthly basin area-weighted average simulated streamflow', &
     &     'inches', Basin_stflow_mo)/=0 ) CALL read_error(3, 'basin_stflow_mo')

      IF ( declvar(MODNAME, 'obsq_inches_mo', 'one', 1, 'double', &
     &     'Monthly measured streamflow at specified outlet station', &
     &     'inches', Obsq_inches_mo)/=0 ) CALL read_error(3, 'obsq_inches_mo')

      ALLOCATE ( Hru_et_yr(Nhru) )
      IF ( declvar(MODNAME, 'hru_et_yr', 'nhru', Nhru, 'double', &
     &     'Yearly area-weighted average actual ET for each HRU', &
     &     'inches', Hru_et_yr)/=0 ) CALL read_error(3, 'hru_et_yr')

      IF ( declvar(MODNAME, 'basin_runoff_ratio_mo', 'one', 1, 'double', &
     &     'Monthly area-weighted average discharge/precipitation ratio', &
     &     'decimal fraction', Basin_runoff_ratio_mo)/=0 ) CALL read_error(3, 'basin_runoff_ratio_mo')

      IF ( declvar(MODNAME, 'basin_runoff_ratio', 'one', 1, 'double', &
     &     'Basin area-weighted average discharge/precipitation ratio', &
     &     'decimal fraction', Basin_runoff_ratio)/=0 ) CALL read_error(3, 'basin_runoff_ratio')

      END FUNCTION sumbdecl

!***********************************************************************
!     sumbinit - Initialize basinsum module - get parameter values
!***********************************************************************
      INTEGER FUNCTION sumbinit()
      USE PRMS_BASINSUM
      USE PRMS_FLOWVARS, ONLY: Basin_soil_moist, Basin_ssstor, Basin_lake_stor
      USE PRMS_INTCP, ONLY: Basin_intcp_stor
      USE PRMS_SNOW, ONLY: Basin_pweqv
      USE PRMS_SRUNOFF, ONLY: Basin_imperv_stor, Basin_dprst_volcl, Basin_dprst_volop
      USE PRMS_GWFLOW, ONLY: Basin_gwstor
      IMPLICIT NONE
! Functions
      INTRINSIC :: MAX, MOD
      INTEGER, EXTERNAL :: getparam, julian_day
      EXTERNAL :: header_print, read_error, write_outfile, PRMS_open_module_file
! Local Variables
      INTEGER :: pftemp
!***********************************************************************
      sumbinit = 0

      IF ( Nobs>0 ) THEN
        IF ( getparam(MODNAME, 'outlet_sta', 1, 'integer', Outlet_sta) &
     &       /=0 ) CALL read_error(2, 'outlet_sta')
        IF ( Outlet_sta==0 ) Outlet_sta = 1
      ENDIF

      IF ( getparam(MODNAME, 'print_type', 1, 'integer', Print_type) &
     &     /=0 ) CALL read_error(2, 'print_type')

      IF ( getparam(MODNAME, 'print_freq', 1, 'integer', Print_freq) &
     &     /=0 ) CALL read_error(2, 'print_freq')

      IF ( Init_vars_from_file==OFF ) THEN
!  Zero stuff out when Timestep = 0
        Watbal_sum = 0.0D0

        Obs_runoff_mo = 0.0D0
        Basin_cfs_mo = 0.0D0
        Basin_ppt_mo = 0.0D0
        Basin_net_ppt_mo = 0.0D0
        Basin_swrad_mo = 0.0D0
        Basin_max_temp_mo = 0.0D0
        Basin_min_temp_mo = 0.0D0
        Basin_intcp_evap_mo = 0.0D0
        Basin_potet_mo = 0.0D0
        Basin_actet_mo = 0.0D0
        Basin_snowmelt_mo = 0.0D0
        Basin_gwflow_mo = 0.0D0
        Basin_ssflow_mo = 0.0D0
        Basin_sroff_mo = 0.0D0
        Basin_stflow_mo = 0.0D0
        Obsq_inches_mo = 0.0D0
        Basin_runoff_ratio_mo = 0.0D0
        Basin_lakeevap_mo = 0.0D0

        Obs_runoff_yr = 0.0D0
        Basin_cfs_yr = 0.0D0
        Basin_ppt_yr = 0.0D0
        Basin_net_ppt_yr = 0.0D0
        Basin_swrad_yr = 0.0D0
        Basin_max_temp_yr = 0.0D0
        Basin_min_temp_yr = 0.0D0
        Basin_intcp_evap_yr = 0.0D0
        Basin_lakeevap_yr = 0.0D0
        Basin_potet_yr = 0.0D0
        Basin_actet_yr = 0.0D0
        Basin_snowmelt_yr = 0.0D0
        Basin_gwflow_yr = 0.0D0
        Basin_ssflow_yr = 0.0D0
        Basin_sroff_yr = 0.0D0
        Basin_stflow_yr = 0.0D0
        Obsq_inches_yr = 0.0D0

        Obs_runoff_tot = 0.0D0
        Basin_cfs_tot = 0.0D0
        Basin_ppt_tot = 0.0D0
        Basin_net_ppt_tot = 0.0D0
        Basin_swrad_tot = 0.0D0
        Basin_max_temp_tot = 0.0D0
        Basin_min_temp_tot = 0.0D0
        Basin_intcp_evap_tot = 0.0D0
        Basin_potet_tot = 0.0D0
        Basin_actet_tot = 0.0D0
        Basin_snowmelt_tot = 0.0D0
        Basin_gwflow_tot = 0.0D0
        Basin_ssflow_tot = 0.0D0
        Basin_sroff_tot = 0.0D0
        Basin_stflow_tot = 0.0D0
        Obsq_inches_tot = 0.0D0
        Hru_et_yr = 0.0D0
        Totdays = 0
        Basin_storage = 0.0D0
        Basin_storvol = 0.0D0
      ENDIF
      Obsq_inches = 0.0D0
      Basin_runoff_ratio = 0.0D0

!******Set daily print switch
      IF ( Print_freq>7 ) THEN
        Dprt = .TRUE.
      ELSE
        Dprt = .FALSE.
      ENDIF

!******Set monthly print switch
      IF ( Print_freq>3 .AND. (Print_freq<8 .OR. Print_freq>11) ) THEN
        Mprt = .TRUE.
      ELSE
        Mprt = .FALSE.
      ENDIF

!******Set yearly print switch
      IF ( Print_freq==2 .OR. Print_freq==3 .OR. Print_freq==6 .OR. Print_freq==7 .OR. &
     &     Print_freq==10.OR.Print_freq==11.OR.Print_freq==14.OR.Print_freq==15 ) THEN
        Yprt = .TRUE.
      ELSE
        Yprt = .FALSE.
      ENDIF

!******Set total print switch
      pftemp = MOD( Print_freq, 2 )
      IF ( pftemp==1 ) THEN
        Tprt = .TRUE.
      ELSE
        Tprt = .FALSE.
      ENDIF

!******Set header print switch (1 prints a new header after every month
!****** summary, 2 prints a new header after every year summary)
      Header_prt = 0
      IF ( Print_freq==6 .OR. Print_freq==7 .OR. Print_freq==10 .OR. Print_freq==11 ) Header_prt = 1
      IF ( Print_freq>=12 ) Header_prt = 2
      IF ( .NOT.Dprt .AND. Print_type==1 ) Header_prt = 3

      Basin_storage = Basin_soil_moist + Basin_intcp_stor + &
     &                Basin_gwstor + Basin_ssstor + Basin_pweqv + &
     &                Basin_imperv_stor + Basin_lake_stor + &
     &                Basin_dprst_volop + Basin_dprst_volcl
!glacier storage not known at start

      IF ( Print_freq/=0 ) THEN
        CALL header_print(Print_type)
! Print span dashes and initial storage
        IF ( Print_type==1 ) THEN
          WRITE (Buffer48, "(' initial', 27X,F9.3)") Basin_storage
          CALL write_outfile(Buffer48(:44))

        ELSEIF ( Print_type==2 ) THEN
          WRITE (Buffer120, 9001) Basin_intcp_stor, &
     &                            Basin_soil_moist, Basin_pweqv, Basin_gwstor, Basin_ssstor
          CALL write_outfile(Buffer120(:98))
        ENDIF
      ENDIF

      Endjday = julian_day('end', 'calendar')

      IF ( Print_debug==4 ) CALL PRMS_open_module_file(BALUNT, 'basin_sum.dbg')

 9001 FORMAT (' initial', 33X, F6.2, 20X, 2F6.2, F13.2, F6.2)

      END FUNCTION sumbinit

!***********************************************************************
!     sumbrun - Computes summary values
!***********************************************************************
      INTEGER FUNCTION sumbrun()
      USE PRMS_BASINSUM
      USE PRMS_BASIN, ONLY: Active_area, Active_hrus, Hru_route_order
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_lakeevap, &
     &    Basin_actet, Basin_perv_et, Basin_swale_et, Hru_actet, Basin_sroff, &
     &    Basin_ssstor, Basin_soil_moist, Basin_cfs, Basin_stflow_out, Basin_lake_stor
      USE PRMS_CLIMATEVARS, ONLY: Basin_swrad, Basin_ppt, Basin_potet, Basin_tmax, Basin_tmin
      USE PRMS_SET_TIME, ONLY: Jday, Modays, Yrdays, Julwater, Nowyear, Nowmonth, Nowday, Cfs2inches
      USE PRMS_OBS, ONLY: Streamflow_cfs
      USE PRMS_GWFLOW, ONLY: Basin_gwflow, Basin_gwstor, Basin_gwsink, Basin_gwstor_minarea_wb
      USE PRMS_INTCP, ONLY: Basin_intcp_evap, Basin_intcp_stor, Basin_net_ppt
      USE PRMS_SNOW, ONLY: Basin_snowmelt, Basin_pweqv, Basin_snowevap
      USE PRMS_GLACR, ONLY: Basin_gl_storage
      USE PRMS_SRUNOFF, ONLY: Basin_imperv_stor, Basin_imperv_evap, &
     &    Basin_dprst_evap, Basin_dprst_volcl, Basin_dprst_volop
      USE PRMS_ROUTING, ONLY: Basin_segment_storage
      IMPLICIT NONE
! Functions
      INTRINSIC :: SNGL, ABS, ALOG, DBLE
      EXTERNAL :: header_print, write_outfile
! Local variables
      INTEGER :: i, j, wyday, endrun, monthdays
      DOUBLE PRECISION :: wat_bal, obsrunoff
!***********************************************************************
      sumbrun = 0

      wyday = Julwater

      IF ( Nowyear==End_year .AND. Jday==Endjday ) THEN
        endrun = 1
      ELSE
        endrun = 0
      ENDIF

!*****Compute aggregated values

      Last_basin_stor = Basin_storage
      Basin_storage = Basin_soil_moist + Basin_intcp_stor + &
     &                Basin_gwstor + Basin_ssstor + Basin_pweqv + &
     &                Basin_imperv_stor + Basin_lake_stor + Basin_dprst_volop + Basin_dprst_volcl
! Basin_storage doesn't include any processes on glacier
! In glacier module, Basin_gl_storstart is an estimate for starting glacier volume, but only
!   includes glaciers that have depth estimates and these are known to be iffy
      IF ( Glacier_flag==ACTIVE ) Basin_storage = Basin_storage + Basin_gl_storage
      IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Strmflow_flag==strmflow_muskingum_module &
     &     .OR. Strmflow_flag==strmflow_muskingum_mann_module) Basin_storage = Basin_storage + Basin_segment_storage

! volume calculation for storage
      Basin_storvol = Basin_storage*Active_area

      obsrunoff = 0.0D0
      IF ( Nobs>0 ) obsrunoff = Streamflow_cfs(Outlet_sta)
      Obsq_inches = obsrunoff*Cfs2inches

      wat_bal = Last_basin_stor - Basin_storage + Basin_ppt + Basin_gwstor_minarea_wb &
     &          - Basin_actet - Basin_stflow_out - Basin_gwsink

      IF ( Basin_stflow_out>0.0 ) THEN
        Basin_runoff_ratio = Basin_ppt/Basin_stflow_out
      ELSE
        Basin_runoff_ratio = 0.0
      ENDIF

      IF ( Print_debug==4 ) THEN
        WRITE ( BALUNT, "(A,2I4,7F8.4)" ) ' bsto-sm-in-gw-ss-sn-iv ', &
     &          Nowmonth, Nowday, Basin_storage, Basin_soil_moist, &
     &          Basin_intcp_stor, Basin_gwstor, Basin_ssstor, &
     &          Basin_pweqv, Basin_imperv_stor

        WRITE ( BALUNT, "(A,I6,8F8.4)" )' bet-pv-iv-in-sn-lk-sw-dp', Nowday, &
     &          Basin_actet, Basin_perv_et, Basin_imperv_evap, &
     &          Basin_intcp_evap, Basin_snowevap, Basin_lakeevap, &
     &          Basin_swale_et, Basin_dprst_evap

        WRITE ( BALUNT, "(A,I6,7F8.4,/)" ) ' bal-pp-et-st-ls-bs-gs ', &
     &          Nowday, wat_bal, Basin_ppt, Basin_actet, Basin_stflow_out, &
     &          Last_basin_stor, Basin_storage, Basin_gwsink
      ENDIF

      Watbal_sum = Watbal_sum + wat_bal

!******Check for daily print

      IF ( Dprt ) THEN
        IF ( Print_type==0 ) THEN
          WRITE ( Buffer40, "(I7,2I5,F11.2,F12.2)" ) Nowyear, &
     &            Nowmonth, Nowday, obsrunoff, Basin_cfs
          CALL write_outfile(Buffer40)

        ELSEIF ( Print_type==1 ) THEN
          WRITE ( Buffer80, "(I7,2I5,7F9.3)" ) Nowyear, &
     &            Nowmonth, Nowday, Basin_ppt, Basin_actet, Basin_storage, &
     &            Basin_stflow_out, Obsq_inches, wat_bal, Watbal_sum
          CALL write_outfile(Buffer80)

        ELSEIF ( Print_type==2 ) THEN
          WRITE ( Buffer151, 9001 ) Nowyear, Nowmonth, Nowday, Basin_swrad, &
     &            Basin_tmax, Basin_tmin, Basin_ppt, Basin_net_ppt, &
     &            Basin_intcp_stor, Basin_intcp_evap, Basin_potet, &
     &            Basin_actet, Basin_soil_moist, Basin_pweqv, &
     &            Basin_snowmelt, Basin_gwstor, Basin_ssstor, &
     &            Basin_gwflow, Basin_ssflow, Basin_sroff, &
     &            Basin_stflow_out, Basin_cfs, obsrunoff, Basin_lakeevap
          CALL write_outfile(Buffer151)

        ENDIF
      ENDIF
      IF ( Print_debug==4 ) WRITE ( BALUNT, * ) 'wat_bal =', wat_bal, &
     &                                          ' watbal_sum=', Watbal_sum

!******Compute monthly values
      IF ( Nowday==1 ) THEN
        Obs_runoff_mo = 0.0D0
        Basin_cfs_mo = 0.0D0
        Basin_ppt_mo = 0.0D0
        Basin_net_ppt_mo = 0.0D0
        Basin_swrad_mo = 0.0D0
        Basin_max_temp_mo = 0.0D0
        Basin_min_temp_mo = 0.0D0
        Basin_intcp_evap_mo = 0.0D0
        Basin_potet_mo = 0.0D0
        Basin_actet_mo = 0.0D0
        Basin_snowmelt_mo = 0.0D0
        Basin_gwflow_mo = 0.0D0
        Basin_ssflow_mo = 0.0D0
        Basin_sroff_mo = 0.0D0
        Basin_stflow_mo = 0.0D0
        Obsq_inches_mo = 0.0D0
        Basin_lakeevap_mo = 0.0D0
      ENDIF

      Obs_runoff_mo = Obs_runoff_mo + obsrunoff
      Obsq_inches_mo = Obsq_inches_mo + obsrunoff*Cfs2inches
      Basin_cfs_mo = Basin_cfs_mo + Basin_cfs
      Basin_ppt_mo = Basin_ppt_mo + Basin_ppt
      Basin_net_ppt_mo = Basin_net_ppt_mo + Basin_net_ppt
      Basin_swrad_mo = Basin_swrad_mo + Basin_swrad
      Basin_max_temp_mo = Basin_max_temp_mo + Basin_tmax
      Basin_min_temp_mo = Basin_min_temp_mo + Basin_tmin
      Basin_intcp_evap_mo = Basin_intcp_evap_mo + Basin_intcp_evap
      Basin_potet_mo = Basin_potet_mo + Basin_potet
      Basin_actet_mo = Basin_actet_mo + Basin_actet
      Basin_snowmelt_mo = Basin_snowmelt_mo + Basin_snowmelt
      Basin_gwflow_mo = Basin_gwflow_mo + Basin_gwflow
      Basin_ssflow_mo = Basin_ssflow_mo + Basin_ssflow
      Basin_sroff_mo = Basin_sroff_mo + Basin_sroff
      Basin_stflow_mo = Basin_stflow_mo + Basin_stflow_out
      Basin_lakeevap_mo = Basin_lakeevap_mo + Basin_lakeevap

      IF ( Nowday==Modays(Nowmonth) ) THEN
        monthdays = Modays(Nowmonth)
        Basin_swrad_mo = Basin_swrad_mo/monthdays
        Basin_max_temp_mo = Basin_max_temp_mo/monthdays
        Basin_min_temp_mo = Basin_min_temp_mo/monthdays
        Obs_runoff_mo = Obs_runoff_mo/monthdays
        Basin_cfs_mo = Basin_cfs_mo/monthdays
        Basin_runoff_ratio_mo = Basin_ppt_mo/monthdays/Basin_stflow_mo
        Basin_lakeevap_mo = Basin_lakeevap_mo/monthdays

        IF ( Mprt ) THEN
          IF ( Print_type==0 ) THEN
            IF ( Dprt ) CALL write_outfile(DASHS(:40))
            WRITE ( Buffer40, "(I7,I5,F16.2,F12.2)" ) Nowyear, Nowmonth, Obs_runoff_mo, Basin_cfs_mo
            CALL write_outfile(Buffer40)
            IF ( Dprt ) CALL write_outfile(DASHS(:40))

          ELSEIF ( Print_type==1 ) THEN
            IF ( Dprt ) CALL write_outfile(DASHS(:62))
            WRITE ( Buffer80, "(I7,I5,5X,5F9.3)" ) Nowyear, &
     &              Nowmonth, Basin_ppt_mo, Basin_actet_mo, Basin_storage, &
     &              Basin_stflow_mo, Obsq_inches_mo
            CALL write_outfile(Buffer80(:62))
            IF ( Dprt ) CALL write_outfile(DASHS(:62))

          ELSEIF ( Print_type==2 ) THEN
            IF ( Dprt ) CALL write_outfile(DASHS)
            WRITE ( Buffer151, 9006 ) Nowyear, Nowmonth, Basin_swrad_mo, Basin_max_temp_mo, &
     &              Basin_min_temp_mo, Basin_ppt_mo, Basin_net_ppt_mo, &
     &              Basin_intcp_evap_mo, Basin_potet_mo, Basin_actet_mo, &
     &              Basin_soil_moist, Basin_pweqv, Basin_snowmelt_mo, &
     &              Basin_gwstor, Basin_ssstor, Basin_gwflow_mo, &
     &              Basin_ssflow_mo, Basin_sroff_mo, Basin_stflow_mo, &
     &              Basin_cfs_mo, Obs_runoff_mo, Basin_lakeevap_mo
            CALL write_outfile(Buffer151)
            IF ( Dprt ) CALL write_outfile(DASHS)
          ENDIF

        ENDIF
      ENDIF

!******Check for year print

      IF ( Yprt ) THEN
        Obs_runoff_yr = Obs_runoff_yr + obsrunoff
        Obsq_inches_yr = Obsq_inches_yr + obsrunoff*Cfs2inches
        Basin_cfs_yr = Basin_cfs_yr + Basin_cfs
        Basin_ppt_yr = Basin_ppt_yr + Basin_ppt
        Basin_net_ppt_yr = Basin_net_ppt_yr + Basin_net_ppt
        Basin_swrad_yr = Basin_swrad_yr + Basin_swrad
        Basin_max_temp_yr = Basin_max_temp_yr + Basin_tmax
        Basin_min_temp_yr = Basin_min_temp_yr + Basin_tmin
        Basin_intcp_evap_yr = Basin_intcp_evap_yr + Basin_intcp_evap
        Basin_lakeevap_yr = Basin_lakeevap_yr + Basin_lakeevap
        Basin_potet_yr = Basin_potet_yr + Basin_potet
        Basin_actet_yr = Basin_actet_yr + Basin_actet
        Basin_snowmelt_yr = Basin_snowmelt_yr + Basin_snowmelt
        Basin_gwflow_yr = Basin_gwflow_yr + Basin_gwflow
        Basin_ssflow_yr = Basin_ssflow_yr + Basin_ssflow
        Basin_sroff_yr = Basin_sroff_yr + Basin_sroff
        Basin_stflow_yr = Basin_stflow_yr + Basin_stflow_out
        DO j = 1, Active_hrus
          i = Hru_route_order(j)
          Hru_et_yr(i) = Hru_et_yr(i) + DBLE( Hru_actet(i) )
        ENDDO

        IF ( wyday==Yrdays ) THEN
          IF ( Print_type==0 ) THEN

            Obs_runoff_yr = Obs_runoff_yr/Yrdays
            Basin_cfs_yr = Basin_cfs_yr/Yrdays
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS(:40))
            WRITE ( Buffer40, "(I7,F21.2,F12.2)" ) Nowyear, Obs_runoff_yr, Basin_cfs_yr
            CALL write_outfile(Buffer40)
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS(:40))

! ****annual summary here
          ELSEIF ( Print_type==1 ) THEN
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS(:62))
            WRITE ( Buffer80, "(I7,10X,5F9.3)" ) Nowyear, Basin_ppt_yr, &
     &              Basin_actet_yr, Basin_storage, Basin_stflow_yr, Obsq_inches_yr
            CALL write_outfile(Buffer80(:62))
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS(:62))

          ELSEIF ( Print_type==2 ) THEN
            Basin_swrad_yr = Basin_swrad_yr/Yrdays
            Basin_max_temp_yr = Basin_max_temp_yr/Yrdays
            Basin_min_temp_yr = Basin_min_temp_yr/Yrdays
            Obs_runoff_yr = Obs_runoff_yr/Yrdays
            Basin_cfs_yr = Basin_cfs_yr/Yrdays
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS)
            WRITE ( Buffer151, 9007 ) Nowyear, Basin_swrad_yr, Basin_max_temp_yr, &
     &              Basin_min_temp_yr, Basin_ppt_yr, Basin_net_ppt_yr, &
     &              Basin_intcp_stor, Basin_intcp_evap_yr, Basin_potet_yr, Basin_actet_yr, &
     &              Basin_soil_moist, Basin_pweqv, Basin_snowmelt_yr, &
     &              Basin_gwstor, Basin_ssstor, Basin_gwflow_yr, &
     &              Basin_ssflow_yr, Basin_sroff_yr, Basin_stflow_yr, &
     &              Basin_cfs_yr, Obs_runoff_yr, Basin_lakeevap_yr
            CALL write_outfile(Buffer151)
            IF ( Mprt .OR. Dprt ) CALL write_outfile(EQULS)
          ENDIF

          Obs_runoff_yr = 0.0D0
          Basin_cfs_yr = 0.0D0
          Basin_ppt_yr = 0.0D0
          Basin_net_ppt_yr = 0.0D0
          Basin_swrad_yr = 0.0D0
          Basin_max_temp_yr = 0.0D0
          Basin_min_temp_yr = 0.0D0
          Basin_intcp_evap_yr = 0.0D0
          Basin_lakeevap_yr = 0.0D0
          Basin_potet_yr = 0.0D0
          Basin_actet_yr = 0.0D0
          Basin_snowmelt_yr = 0.0D0
          Basin_gwflow_yr = 0.0D0
          Basin_ssflow_yr = 0.0D0
          Basin_sroff_yr = 0.0D0
          Basin_stflow_yr = 0.0D0
          Obsq_inches_yr = 0.0D0
          Hru_et_yr = 0.0D0

        ENDIF
      ENDIF

!******Check for total print

      IF ( Tprt ) THEN
!******Print heading if needed
!        IF ( endrun==1 ) THEN
!            CALL write_outfile(' ')
!            IF ( .NOT.Dprt .OR. .NOT.Mprt .OR. .NOT.Yprt ) CALL header_print(Print_type)
!        ENDIF

        Totdays = Totdays + 1
        Obs_runoff_tot = Obs_runoff_tot + obsrunoff
        Obsq_inches_tot = Obsq_inches_tot + obsrunoff*Cfs2inches
        Basin_cfs_tot = Basin_cfs_tot + Basin_cfs
        Basin_ppt_tot = Basin_ppt_tot + Basin_ppt
        Basin_net_ppt_tot = Basin_net_ppt_tot + Basin_net_ppt
        Basin_swrad_tot = Basin_swrad_tot + Basin_swrad
        Basin_max_temp_tot = Basin_max_temp_tot + Basin_tmax
        Basin_min_temp_tot = Basin_min_temp_tot + Basin_tmin
        Basin_intcp_evap_tot = Basin_intcp_evap_tot + Basin_intcp_evap
        Basin_potet_tot = Basin_potet_tot + Basin_potet
        Basin_actet_tot = Basin_actet_tot + Basin_actet
        Basin_snowmelt_tot = Basin_snowmelt_tot + Basin_snowmelt
        Basin_gwflow_tot = Basin_gwflow_tot + Basin_gwflow
        Basin_ssflow_tot = Basin_ssflow_tot + Basin_ssflow
        Basin_sroff_tot = Basin_sroff_tot + Basin_sroff
        Basin_stflow_tot = Basin_stflow_tot + Basin_stflow_out

        IF ( endrun==1 ) THEN

          IF ( Print_type==0 ) THEN
            Obs_runoff_tot = Obs_runoff_tot/Totdays
            Basin_cfs_tot = Basin_cfs_tot/Totdays
            CALL write_outfile(STARS(:40))
            WRITE ( Buffer48, "(A,F14.2,F12.2)" ) ' Total for run', Obs_runoff_tot, Basin_cfs_tot
            CALL write_outfile(Buffer48(:40))
            CALL write_outfile(STARS(:40))

          ELSEIF ( Print_type==1 ) THEN
            CALL write_outfile(STARS(:62))
            WRITE ( Buffer80, 9005 ) ' Total for run', Basin_ppt_tot, &
     &              Basin_actet_tot, Basin_storage, Basin_stflow_tot, Obsq_inches_tot
            CALL write_outfile(Buffer80(:62))
            CALL write_outfile(STARS(:62))

          ELSEIF ( Print_type==2 ) THEN
            Obs_runoff_tot = Obs_runoff_tot/Totdays
            Basin_cfs_tot = Basin_cfs_tot/Totdays
            CALL write_outfile(STARS)
            WRITE ( Buffer151, 9004 ) ' Total for run', Basin_ppt_tot, &
     &              Basin_net_ppt_tot, Basin_intcp_evap_tot, &
     &              Basin_potet_tot, Basin_actet_tot, Basin_soil_moist, &
     &              Basin_pweqv, Basin_snowmelt_tot, Basin_gwstor, &
     &              Basin_ssstor, Basin_gwflow_tot, Basin_ssflow_tot, &
     &              Basin_sroff_tot, Basin_stflow_tot, Basin_cfs_tot, Obs_runoff_tot, Basin_lakeevap_yr
            CALL write_outfile(Buffer151)
            CALL write_outfile(STARS)
          ENDIF
        ENDIF
      ENDIF

 9001 FORMAT (I6, 2I3, F5.0, 2F5.1, 2F7.2, 2F6.2, 2F7.2, F6.2, F6.3, F7.3, 2F6.3, 3F7.2, F7.4, F9.1, F9.2, F7.2)
 9004 FORMAT (A, 13X, 2F7.2, F12.1, 2F7.2, 2F6.2, F7.2, 2F6.2, 4F7.2, F9.1, F9.2, F7.2)
 9005 FORMAT (A, 3X, 6F9.3)
 9006 FORMAT (I6, I3, 3X, 3F5.1, 2F7.2, F12.1, 2F7.2, 2F6.2, F7.2, 2F6.2, 3F7.2, F9.1, F9.2, 2F7.2)
 9007 FORMAT (I6, 6X, 3F5.1, 2F7.2, 2F6.2, 2F7.2, 2F6.2, F7.2, 2F6.2, 3F7.2, F9.2, F9.2, 2F7.2)

      END FUNCTION sumbrun

!***********************************************************************
! Print headers for tables
! This writes the measured and simulated table header.
!***********************************************************************
      SUBROUTINE header_print(Print_type)
      USE PRMS_BASINSUM, ONLY: DASHS, Buffer80, Print_freq, Header_prt
      IMPLICIT NONE
      EXTERNAL :: write_outfile
! Arguments
      INTEGER, INTENT(IN) :: Print_type
!***********************************************************************
      CALL write_outfile(' ')
!  This writes the water balance table header.
      IF ( Header_prt==3 ) THEN
        CALL write_outfile('   Year Month Day   Precip     ET    Storage S-Runoff M-Runoff')
        WRITE (Buffer80, 9002)
        CALL write_outfile(Buffer80(:62))
        CALL write_outfile(DASHS(:62))

      ELSEIF ( Print_type==0 ) THEN
        IF ( Print_freq==1 ) THEN
          CALL write_outfile('                    Measured   Simulated')
        ELSE
          CALL write_outfile('   Year Month Day   Measured   Simulated')
        ENDIF
        CALL write_outfile('                      (cfs)      (cfs)')
        CALL write_outfile(DASHS(:40))

!  This writes the water balance table header.
      ELSEIF ( Print_type==1 ) THEN
        CALL write_outfile('   Year Month Day   Precip     ET    Storage S-Runoff M-Runoff   Watbal  WBalSum')
        WRITE (Buffer80, 9001)
        CALL write_outfile(Buffer80)
        CALL write_outfile(DASHS(:80))

!  This writes the detailed table header.
      ELSEIF ( Print_type==2 ) THEN
        CALL write_outfile(' Year mo day srad  tmx  tmn    ppt  n-ppt  ints  intl  potet'// &
     &    '  actet  smav pweqv   melt gwsto sssto gwflow ssflow  sroff tot-fl     sim     meas  lkevap')
        CALL write_outfile('             (ly) (F/C)(F/C)  (in)   (in)  (in)  (in)   (in)'// &
     &    '   (in)  (in)  (in)   (in)  (in)  (in)   (in)   (in)   (in)   (in)    (cfs)    (cfs)   (in)')
        CALL write_outfile(DASHS)

      ENDIF

 9001 FORMAT (17X, 7(' (inches)'))
 9002 FORMAT (17X, 5(' (inches)'))

      END SUBROUTINE header_print

!***********************************************************************
!     Write or read restart file
!***********************************************************************
      SUBROUTINE basin_sum_restart(In_out)
      USE PRMS_CONSTANTS, ONLY: SAVE_INIT
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_BASINSUM
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      ! Functions
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=9) :: module_name
!***********************************************************************
      IF ( In_out==SAVE_INIT ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Totdays, Obs_runoff_mo, Obs_runoff_yr, Obs_runoff_tot, Watbal_sum
        WRITE ( Restart_outunit ) Basin_cfs_mo, Basin_cfs_yr, Basin_cfs_tot, Basin_net_ppt_yr, Basin_net_ppt_tot
        WRITE ( Restart_outunit ) Basin_max_temp_yr, Basin_max_temp_tot, Basin_min_temp_yr, Basin_min_temp_tot
        WRITE ( Restart_outunit ) Basin_potet_yr, Basin_potet_tot, Basin_actet_yr, Basin_actet_tot
        WRITE ( Restart_outunit ) Last_basin_stor, Basin_snowmelt_yr, Basin_snowmelt_tot, Basin_gwflow_yr, Basin_gwflow_tot
        WRITE ( Restart_outunit ) Basin_ssflow_yr, Basin_ssflow_tot, Basin_sroff_yr, Basin_sroff_tot
        WRITE ( Restart_outunit ) Basin_stflow_yr, Basin_stflow_tot, Basin_ppt_yr, Basin_ppt_tot
        WRITE ( Restart_outunit ) Basin_intcp_evap_yr, Basin_intcp_evap_tot, Obsq_inches_yr, Obsq_inches_tot, Basin_lakeevap_yr
        WRITE ( Restart_outunit ) Basin_net_ppt_mo, Obsq_inches_mo, Basin_max_temp_mo, Basin_min_temp_mo, Basin_actet_mo
        WRITE ( Restart_outunit ) Basin_snowmelt_mo, Basin_gwflow_mo, Basin_sroff_mo, Basin_stflow_mo
        WRITE ( Restart_outunit ) Basin_intcp_evap_mo, Basin_storage, Basin_storvol, Basin_potet_mo
        WRITE ( Restart_outunit ) Basin_ssflow_mo, Basin_ppt_mo, Basin_runoff_ratio_mo
        WRITE ( Restart_outunit ) Hru_et_yr
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Totdays, Obs_runoff_mo, Obs_runoff_yr, Obs_runoff_tot, Watbal_sum
        READ ( Restart_inunit ) Basin_cfs_mo, Basin_cfs_yr, Basin_cfs_tot, Basin_net_ppt_yr, Basin_net_ppt_tot
        READ ( Restart_inunit ) Basin_max_temp_yr, Basin_max_temp_tot, Basin_min_temp_yr, Basin_min_temp_tot
        READ ( Restart_inunit ) Basin_potet_yr, Basin_potet_tot, Basin_actet_yr, Basin_actet_tot
        READ ( Restart_inunit ) Last_basin_stor, Basin_snowmelt_yr, Basin_snowmelt_tot, Basin_gwflow_yr, Basin_gwflow_tot
        READ ( Restart_inunit ) Basin_ssflow_yr, Basin_ssflow_tot, Basin_sroff_yr, Basin_sroff_tot
        READ ( Restart_inunit ) Basin_stflow_yr, Basin_stflow_tot, Basin_ppt_yr, Basin_ppt_tot
        READ ( Restart_inunit ) Basin_intcp_evap_yr, Basin_intcp_evap_tot, Obsq_inches_yr, Obsq_inches_tot, Basin_lakeevap_yr
        READ ( Restart_inunit ) Basin_net_ppt_mo, Obsq_inches_mo, Basin_max_temp_mo, Basin_min_temp_mo, Basin_actet_mo
        READ ( Restart_inunit ) Basin_snowmelt_mo, Basin_gwflow_mo, Basin_sroff_mo, Basin_stflow_mo
        READ ( Restart_inunit ) Basin_intcp_evap_mo, Basin_storage, Basin_storvol, Basin_potet_mo
        READ ( Restart_inunit ) Basin_ssflow_mo, Basin_ppt_mo, Basin_runoff_ratio_mo
        READ ( Restart_inunit ) Hru_et_yr
      ENDIF
      END SUBROUTINE basin_sum_restart
