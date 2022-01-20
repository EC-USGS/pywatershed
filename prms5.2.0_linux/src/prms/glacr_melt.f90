!**********************************************************************
! Computes glacier runoff for each glacier using three
! linear reservoirs (snow, firn, ice) with time lapses and ability
! to advance or retreat according to Bahr(1997) volume-area scaling.
! This theory has been advanced according to Arendt and others(2006) for
! the scaling constants and Luthi(2009).
!
! ELAs are computed yearly, as well as mass balance, maximum/winter and
! minimum/summer. These can be used for calibration. Note that the calculation
! of these is only appropriate for the northern hemisphere and the code
! needs to be modified if applied in the southern hemisphere. BUT, snow model would
! also be wrong with reset dates for WY Oct 1 resets.
!
! Permafrost/frozen ground calculations should be turned on with glaciers. But they
! can also be turned on without glaciers.
!
! Hru_elev_ts, Hru_elev_feet, Hru_elev_meters, and Hru_slope for glaciers are
! adjusted each timestep. Hru_aspect is not very realistic, since the glacier HRUs
! are not 2-plane and thus have widely different facing areas (they are convex when
! glacierized and concave when deglacierized). An improvement might be to make
! the glaciers 2-planed, but the code would have to be changed vastly to
! reflect this, since there would be two termini HRUs and would want to
! change the glacier_frac value the same with the two termini. For now,
! it makes sense to keep Hru_aspect constant through time as the way the
! glacier predominately faces.
!
! HRUs with glaciers must have parameter glacier_frac(i)=1, unless they
! are at the terminus of the glacier (in which case they can have
! glacier_frac(i)<1). Hru numbering goes from largest HRU ID at top of glacier to
! smallest at ID at bottom (the way Weasel delineation was designed). The parameter
! Glac_HRUnum_down = 1 then in the init function. If the opposite direction,
! then set Glac_HRUnum_down = 0. IDs need to be stacked.
!
! HRUs containing insubstantial (relative to basin) glaciers have their glaciated
! fraction as glrette_frac(i)>0 (but <1)
!
! NOTE: Multiple branches are possible in the melt generation, but basal topography
! calculations will be mathematically unsound as each branch will be considered a
! different glacier.
!
! modified June 2012 by Steve Regan
! modified Jan 2017 by AE Van Beusekom
! dedicated calibration variables removed 2019
!
!***********************************************************************

      MODULE PRMS_GLACR
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, DOCUMENTATION, MONTHS_PER_YEAR, GLACIER, LAND, &
     &    FEET2METERS, METERS2FEET, DNEARZERO, NEARZERO, FEET, METERS, MAX_DAYS_PER_YEAR
      USE PRMS_MODULE, ONLY: Nhru, Model, Init_vars_from_file
      IMPLICIT NONE
      !****************************************************************
      !   Local Variables
      character(len=*), parameter :: MODDESC = 'Glacier Dynamics'
      character(len=10), parameter :: MODNAME = 'glacr_melt'
      character(len=*), parameter :: Version_glacr = '2020-12-02'
      ! Ngl - Number of glaciers counted by termini
      ! Ntp - Number of tops of glaciers, so max glaciers that could ever split in two
      ! Nhrugl - Number of at least partially glacierized hrus at initiation
!#of cells=Nhrugl,#of streams=Ntp,#of cells/stream<=Ntp, #of glaciers<=Nhru
      INTEGER, SAVE :: Nglres, Ngl, Ntp, Nhrugl, MbInit_flag, Output_unit, Fraw_unit, All_unit
      INTEGER, SAVE :: Seven, Four, Glac_HRUnum_down
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_area_inch2(:)
      REAL, PARAMETER :: Gravity = 9.8 ! m/s2
      REAL, PARAMETER :: Aflow = 1.e-25 ! Pa^-3/s, Farinotti 2009 could be 2.4e-24, could be 1e-26 see Patterson 2010
      REAL, PARAMETER :: Density = 917.0 ! kg/m3
      DOUBLE PRECISION, PARAMETER :: Acre_inch2 = 43560.0D0*12.0D0*12.0D0

      !****************************************************************
      !   Declared Variables

      REAL, SAVE, ALLOCATABLE :: Hru_glres_melt(:), Glrette_melt(:), Gl_ice_melt(:), Glacr_elev_init(:)
      REAL, SAVE, ALLOCATABLE :: Basal_elev(:), Basal_slope(:), Keep_gl(:,:), Prev_outi(:, :), Prev_out(:, :)
      REAL, SAVE, ALLOCATABLE :: Ode_glacrva_coef(:), Av_basal_slope(:), Av_fgrad(:), Hru_slope_ts(:)
      REAL, SAVE, ALLOCATABLE :: Hru_mb_yrend(:), Glacr_flow(:), Glacr_slope_init(:), Gl_top_melt(:)
      INTEGER, SAVE, ALLOCATABLE :: Top(:), Term(:), Top_tag(:), Ela(:), Order_flowline(:)
      INTEGER, SAVE, ALLOCATABLE :: Glacr_tag(:), Ikeep_gl(:,:), Tohru(:)
      DOUBLE PRECISION, SAVE :: Basin_gl_ice_melt, Basin_gl_area, Basin_gl_top_melt
      DOUBLE PRECISION, SAVE :: Basin_gl_top_gain, Basin_gl_storvol, Basin_gl_storage
      DOUBLE PRECISION, SAVE :: Basin_gl_storstart
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Hru_mb_yrcumul(:), Delta_volyr(:), Prev_vol(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Prev_area(:), Gl_mb_yrcumul(:), Gl_area(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gl_mb_cumul(:), Glnet_ar_delta(:), Gl_mbc_yrend(:)

      !****************************************************************
      !   Declared Parameters

      REAL, SAVE :: Max_gldepth
      REAL, SAVE, ALLOCATABLE :: Glacrva_coef(:), Glacrva_exp(:), Hru_length(:), Hru_width(:)
      REAL, SAVE, ALLOCATABLE :: Stor_ice(:,:), Stor_snow(:,:), Stor_firn(:,:)
      REAL, SAVE, ALLOCATABLE :: Hru_slope(:), Abl_elev_range(:)

      END MODULE PRMS_GLACR

!***********************************************************************
!     Main glacr routine
!***********************************************************************
      INTEGER FUNCTION glacr()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, RUN, SETDIMENS, DECL, INIT, CLEAN
      USE PRMS_MODULE, ONLY: Process_flag, Save_vars_to_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: glacrdecl, glacrinit, glacrrun, glacrsetdims
      EXTERNAL :: glacr_restart
!***********************************************************************
      glacr = 0

      IF ( Process_flag==RUN ) THEN
        glacr = glacrrun()
      ELSEIF ( Process_flag==SETDIMENS ) THEN
        glacr = glacrsetdims()
      ELSEIF ( Process_flag==DECL ) THEN
        glacr = glacrdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        glacr = glacrinit()
      ELSEIF ( Process_flag==CLEAN ) THEN
        IF ( Save_vars_to_file==ACTIVE ) CALL glacr_restart(0)
      ENDIF

      END FUNCTION GLACR

!***********************************************************************
!     glacrsetdims - declares glacier module specific dimensions
!***********************************************************************
      INTEGER FUNCTION glacrsetdims()
      USE PRMS_GLACR, ONLY: Nglres, Seven, Four, MbInit_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declfix, control_integer
      EXTERNAL :: read_error
!***********************************************************************
      glacrsetdims = 0

      IF ( declfix('nglres', 3, 3, 'Number of reservoirs in a glacier')/=0 ) CALL read_error(7, 'nglres')
      Nglres = 3
      IF ( declfix('seven', 7, 7, 'Need for keeping glacier variable real array ')/=0 ) CALL read_error(7, 'seven')
      Seven = 7
      IF ( declfix('four',4, 4, 'Need for keeping glacier variable integer array')/=0 ) CALL read_error(7, 'four')
      Four = 4

      IF ( control_integer(MbInit_flag, 'mbInit_flag')/=0 ) MbInit_flag = 0

      END FUNCTION glacrsetdims

!***********************************************************************
!     glacrdecl - declare parameters and variables for glacier runoff
!***********************************************************************
      INTEGER FUNCTION glacrdecl()
      USE PRMS_GLACR
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL :: read_error, print_module
!***********************************************************************
      glacrdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_glacr)

      IF ( declvar(MODNAME, 'nhrugl', 'one', 1, 'integer',                &
           'Number of at least partially glacierized HRUs at initiation', &
           'none', Nhrugl)/=0 ) CALL read_error(3, 'nhrugl')

! declare variables
      ALLOCATE ( Hru_slope_ts(Nhru) )
      IF ( declvar(MODNAME, 'hru_slope_ts', 'nhru', Nhru, 'real',   &
     &     'HRU slope for timestep, which can change for glaciers', &
     &     'decimal fraction', Hru_slope_ts)/=0 ) CALL read_error(3, 'hru_slope_ts')

      IF ( declvar(MODNAME, 'basin_gl_top_melt', 'one', 1, 'double',     &
     &     'Basin area-weighted glacier surface melt (snow, ice and rain) coming out of termini of all glaciers and glacierettes', &
     &     'inches', Basin_gl_top_melt)/=0 ) CALL read_error(3, 'basin_gl_top_melt')

      IF ( declvar(MODNAME, 'basin_gl_top_gain', 'one', 1, 'double', &
     &     'Basin area-weighted glacier surface gain (snow and rain minus evaporation) for all glaciers and glacierettes', &
     &     'inches', Basin_gl_top_gain)/=0 ) CALL read_error(3, 'basin_gl_top_gain')

      IF ( declvar(MODNAME, 'basin_gl_ice_melt', 'one', 1, 'double',    &
     &     'Basin area-weighted glacier ice (firn) melt coming out of termini of all glaciers and glacierettes', &
     &     'inches', Basin_gl_ice_melt)/=0 ) CALL read_error(3, 'basin_gl_ice_melt')

      ALLOCATE ( Gl_mb_yrcumul(Nhru) )
      IF ( declvar(MODNAME, 'gl_mb_yrcumul', 'nhru', Nhru, 'double',     &
     &     'Yearly mass balance for each glacier, indexed by Glacr_tag', &
     &     'inches', Gl_mb_yrcumul)/=0 ) CALL read_error(3, 'gl_mb_yrcumul')

      ALLOCATE ( Gl_mb_cumul(Nhru) )
      IF ( declvar(MODNAME, 'gl_mb_cumul', 'nhru', Nhru, 'double', &
     &     'Cumulative mass balance for each glacier since start day, indexed by Glacr_tag', &
     &     'inches', Gl_mb_cumul)/=0 ) CALL read_error(3, 'gl_mb_cumul')

      IF ( declvar(MODNAME, 'basin_gl_area', 'one', 1, 'double', &
     &     'Basin area-weighted average glacier-covered area',   &
     &     'decimal fraction', Basin_gl_area)/=0 ) CALL read_error(3, 'basin_gl_area')

      ALLOCATE ( Gl_area(Nhru) )
      IF ( declvar(MODNAME, 'gl_area', 'nhru', Nhru, 'double', &
     &     'Area of each glacier, indexed by Glacr_tag',       &
     &     'acres', Gl_area)/=0 ) CALL read_error(3, 'gl_area')

      ALLOCATE ( Glnet_ar_delta(Nhru) )
      IF ( declvar(MODNAME, 'glnet_ar_delta', 'nhru', Nhru, 'double',           &
     &     'Sum of area change of each glacier since start year, indexed by Glacr_tag', &
     &     'acres', Glnet_ar_delta)/=0 ) CALL read_error(3, 'glnet_ar_delta')

      ALLOCATE ( Glacr_flow(Nhru) )
      IF ( declvar(MODNAME, 'glacr_flow', 'nhru', Nhru, 'real',         &
     &     'Glacier melt and rain from HRU to stream network, only nonzero at termini HRUs and snowfield HRUs',  &
     &     'inches cubed', Glacr_flow)/=0 ) CALL read_error(3, 'glacr_flow')

      ALLOCATE ( Delta_volyr(Nhru) )
      IF ( declvar(MODNAME, 'delta_volyr', 'nhru', Nhru, 'double',      &
     &     'Year total volume change for each glacier, indexed by Glacr_tag',          &
     &     'inches cubed', Delta_volyr)/=0 ) CALL read_error(3, 'delta_volyr')

      ALLOCATE ( Hru_mb_yrcumul(Nhru) )
      IF ( declvar(MODNAME, 'hru_mb_yrcumul', 'nhru', Nhru, 'double',        &
     &     'Mass balance for a glacier HRU, cumulative for year',       &
     &     'inches', Hru_mb_yrcumul)/=0 ) CALL read_error(3, 'hru_mb_yrcumul')

      ALLOCATE ( Top_tag(Nhru) )
      IF ( declvar(MODNAME, 'top_tag', 'nhru', Nhru, 'integer',         &
     &     'Identifies which glacier top each HRU is fed by. If =-1, then has multiple feeders', &
     &     'none', Top_tag)/=0 ) CALL read_error(3, 'top_tag')

      ALLOCATE ( Glacr_tag(Nhru) )
      IF ( declvar(MODNAME, 'glacr_tag', 'nhru', Nhru, 'integer',       &
     &     'Identifies which glacier each HRU belongs to',              &
     &     'none', Glacr_tag)/=0 ) CALL read_error(3, 'glacr_tag')

      ALLOCATE ( Prev_area(Nhru) )
      IF ( declvar(MODNAME, 'prev_area', 'nhru', Nhru, 'double',        &
     &     'Previous year glacier-covered area above each HRU where all branches of the glacier are included', &
     &     'inches squared', Prev_area)/=0 ) CALL read_error(3, 'prev_area')

      ALLOCATE ( Prev_vol(Nhru) )
      IF ( declvar(MODNAME, 'prev_vol', 'nhru', Nhru, 'double',         &
     &     'Previous volume of each glacier, indexed by Glacr_tag', &
     &     'inches cubed', Prev_vol)/=0 ) CALL read_error(3, 'prev_vol')

      ALLOCATE ( Prev_out(Nhru,Nglres) )
      IF ( declvar(MODNAME, 'prev_out', 'nhru,nglres', Nhru*Nglres, 'real',  &
     &     'Antecedent outflow of the 3 reservoirs in each glacier, indexed by Glacr_tag',&
     &     'inches cubed', Prev_out)/=0 ) CALL read_error(3, 'prev_out')

      ALLOCATE ( Prev_outi(Nhru,Nglres) )
      IF ( declvar(MODNAME, 'prev_outi', 'nhru,nglres', Nhru*Nglres, 'real',  &
     &     'Antecedent outflow of the 3 reservoirs in each glacier for only ice (firn) melt, indexed by Glacr_tag',&
     &     'inches cubed', Prev_outi)/=0 ) CALL read_error(3, 'prev_outi')

      ALLOCATE ( Order_flowline(Nhru) )
      IF ( declvar(MODNAME, 'order_flowline', 'nhru', Nhru, 'integer',  &
     &     'Order of flowlines that belong together as glaciers, Ntp of these', &
     &     'none', Order_flowline)/=0 ) CALL read_error(3, 'order_flowline')

      ALLOCATE ( Ode_glacrva_coef(Nhru) )
      IF ( declvar(MODNAME, 'ode_glacrva_coef', 'nhru', Nhru, 'real',  &
     &     'Estimate of glacrva_coef from ODE basal topography of each glacier, indexed by Glacr_tag', &
     &     'm**(3-2*glacrva_exp)', Ode_glacrva_coef)/=0 ) CALL read_error(3, 'ode_glacrva_coef')

      ALLOCATE ( Ela(Nhru) )
      IF ( declvar(MODNAME, 'ela', 'nhru', Nhru, 'integer',             &
     &     'HRU number at ELA corresponding to each top in each glacier, Ntp of these', &
     &     'none', Ela)/=0 ) CALL read_error(3, 'ela')

      ALLOCATE ( Top(Nhru) )
      IF ( declvar(MODNAME, 'top', 'nhru', Nhru, 'integer',             &
     &     'HRU number at tops of each glacier, Ntp of these',          &
     &     'none', Top)/=0 ) CALL read_error(3, 'top')

      ALLOCATE ( Term(Nhru) )
      IF ( declvar(MODNAME, 'term', 'nhru', Nhru, 'integer',            &
     &     'HRU number at terminus of each glacier, Ngl of these',&
     &     'none', Term)/=0 ) CALL read_error(3, 'term')

      ALLOCATE ( Hru_mb_yrend(Nhru) )
      IF ( declvar(MODNAME, 'hru_mb_yrend', 'nhru', Nhru, 'real',       &
     &     'Glacier HRU mass balance at end of previous hydrological year', &
     &      'inches', Hru_mb_yrend)/=0 ) CALL read_error(3, 'hru_mb_yrend')

      ALLOCATE ( Av_fgrad(Nhru) )
      IF ( declvar(MODNAME, 'av_fgrad', 'nhru', Nhru, 'real',           &
     &     'Glacier average HRU mass balance gradient with elevation at flowline at end of each hydrological year, Ngl of these',&
     &     'decimal fraction', Av_fgrad)/=0 ) CALL read_error(3, 'av_fgrad')

      ALLOCATE ( Hru_glres_melt(Nhru) )
      IF ( declvar(MODNAME, 'hru_glres_melt', 'nhru', Nhru, 'real',     &
     &     'Amount of glacier surface melt (snow, ice, rain) from an HRU that goes into reservoirs', &
     &     'inches', Hru_glres_melt)/=0 ) CALL read_error(3, 'hru_glres_melt')

      ALLOCATE ( Glrette_melt(Nhru) )
      IF ( declvar(MODNAME, 'glrette_melt', 'nhru', Nhru, 'real',     &
     &     'Amount of glacierette surface melt (snow, ice, rain) from an HRU', &
     &     'inches', Glrette_melt)/=0 ) CALL read_error(3, 'glrette_melt')

      ALLOCATE ( Gl_top_melt(Nhru) )
      IF ( declvar(MODNAME, 'gl_top_melt', 'nhru', Nhru, 'real',        &
     &     'Amount of glacier surface melt (snow, ice, rain) coming out of terminus of glacier, indexed by Glacr_tag', &
     &     'inches', Gl_top_melt)/=0 ) CALL read_error(3, 'gl_top_melt')

      ALLOCATE ( Gl_ice_melt(Nhru) )
      IF ( declvar(MODNAME, 'gl_ice_melt', 'nhru', Nhru, 'real',        &
     &     'Amount of glacier ice (firn) melt coming out of terminus of glacier, indexed by Glacr_tag', &
     &     'inches', Gl_ice_melt)/=0 ) CALL read_error(3, 'gl_ice_melt')

      ALLOCATE ( Basal_elev(Nhru) )
      IF ( declvar(MODNAME, 'basal_elev', 'nhru', Nhru, 'real',         &
     &     'Glacier basal elevation mean over HRU',                     &
     &     'elev_units', Basal_elev)/=0 ) CALL read_error(3, 'basal_elev')

      ALLOCATE ( Keep_gl(Nhru,Seven) )
      IF ( declvar(MODNAME, 'keep_gl', 'nhru,seven', Nhru*Seven, 'real', &
     &     'Glacier real variables keeping from first year', &
     &     'none', Keep_gl)/=0 ) CALL read_error(3, 'keep_gl')

      ALLOCATE ( Ikeep_gl(Nhru,Four) )
      IF ( declvar(MODNAME, 'ikeep_gl', 'nhru,four', Nhru*Four, 'integer', &
     &     'Glacier integer variables keeping from first year', &
     &     'none', Ikeep_gl)/=0 ) CALL read_error(3, 'ikeep_gl')

      ALLOCATE ( Basal_slope(Nhru) )
      IF ( declvar(MODNAME, 'basal_slope', 'nhru', Nhru, 'real',        &
     &     'Glacier basal slope down flowline mean over each HRU', &
     &     'decimal fraction', Basal_slope)/=0 ) CALL read_error(3, 'basal_slope')

      ALLOCATE ( Av_basal_slope(Nhru) )
      IF ( declvar(MODNAME, 'av_basal_slope', 'nhru', Nhru, 'real',      &
     &     'Glacier average basal slope at flowline location, indexed by Glacr_tag', &
     &     'decimal fraction', Av_basal_slope)/=0 ) CALL read_error(3, 'av_basal_slope')

      IF ( declvar(MODNAME, 'basin_gl_storage', 'one', 1, 'double', &
     &     'Basin area-weighted average storage change in glacier reservoirs', &
     &     'inches', Basin_gl_storage)/=0 ) CALL read_error(3, 'basin_gl_storage')

      IF ( declvar(MODNAME, 'basin_gl_storstart', 'one', 1, 'double', &
     &     'Basin area-weighted average storage estimated start in glacier reservoirs', &
     &     'inches', Basin_gl_storstart)/=0 ) CALL read_error(3, 'basin_gl_storstart')

      IF ( declvar(MODNAME, 'basin_gl_storvol', 'one', 1, 'double', &
     &     'Basin storage volume in glacier storage reservoirs', &
     &     'acre-inches', Basin_gl_storvol)/=0 ) CALL read_error(3, 'basin_gl_storvol')

      IF ( Init_vars_from_file==0 ) THEN
        ALLOCATE ( Glacr_elev_init(Nhru) )
        IF ( declvar(MODNAME, 'glacr_elev_init', 'nhru', Nhru, 'real',    &
     &     'Glacier surface elevation mean over HRU at initiation extrapolating to 100% glacierized HRU', &
     &     'elev_units', Glacr_elev_init)/=0 ) CALL read_error(3, 'glacr_elev_init')

        ALLOCATE ( Glacr_slope_init(Nhru) )
        IF ( declvar(MODNAME, 'glacr_slope_init', 'nhru', Nhru, 'real',    &
     &     'Glacier surface slope mean over HRU at initiation extrapolating to 100% glacierized HRU', &
     &     'elev_units', Glacr_slope_init)/=0 ) CALL read_error(3, 'glacr_slope_init')
      ENDIF

      ! local arrays
      ALLOCATE ( Hru_area_inch2(Nhru), Gl_mbc_yrend(Nhru) )

! declare parameters
      ALLOCATE ( Tohru(Nhru) )
      IF ( declparam(MODNAME, 'tohru', 'nhru', 'integer', &
     &     '0', 'bounded', 'nhru', &
     &     'The index of the down-flowline HRU for a glacier', &
     &     'Index of down-flowline HRU to which the HRU'// &
     &     ' glacier melt flows, for non-glacier HRUs that do not flow to another HRU enter 0', &
     &     'none')/=0 ) CALL read_error(1, 'tohru')

      ALLOCATE ( Hru_slope(Nhru) )
      IF ( declparam(MODNAME, 'hru_slope', 'nhru', 'real', &
     &     '0.0', '0.0', '10.0', &
     &     'HRU slope', &
     &     'Slope of each HRU, specified as change in vertical length divided by change in horizontal length', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'hru_slope')

      IF ( declparam(MODNAME, 'max_gldepth', 'one', 'real', &
           '1.5', '0.1', '3.0', &
           'Upper bound on glacier thickness, thickest glacier measured is Taku at 1.5 km, ice sheet 3 km', &
           'Upper bound on glacier thickness, thickest glacier measured is Taku at 1.5 km, ice sheet 3 km', &
           'km')/=0 ) CALL read_error(1, 'max_gldepth')

      ALLOCATE ( Glacrva_coef(Nhru) )
      IF ( declparam(MODNAME, 'glacrva_coef', 'nhru', 'real', &
           '0.28', '0.01', '2.0', &
           'Volume area scaling coefficient for glaciers with Luthi 2009', &
           'Volume area scaling coefficient for glaciers, average value by region', &
           'm**(3-2*glacrva_exp)')/=0 ) CALL read_error(1, 'glacrva_coef')

      ALLOCATE ( Glacrva_exp(Nhru) )
      IF ( declparam(MODNAME, 'glacrva_exp', 'nhru', 'real', &
           '1.375', '1.0', '2.0', &
           'Volume area exponential coefficient for glaciers', &
           'Volume area exponential coefficient for glaciers, average value by region', &
           'none')/=0 ) CALL read_error(1, 'glacrva_exp')

      ALLOCATE ( Stor_ice(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'stor_ice', 'nhru,nmonths', 'real', &
           '10.0', '5.0', '29.0', &
           'Monthly Storage coefficient for ice melt on glaciers', &
           'Monthly (January to December) Storage coefficient for ice melt on glaciers', &
           'hours')/=0 ) CALL read_error(1, 'stor_ice')

      ALLOCATE ( Stor_snow(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'stor_snow', 'nhru,nmonths', 'real', &
           '80.0', '30.0', '149.0', &
           'Monthly Storage coefficient for snow melt on glaciers', &
           'Monthly (January to December) Storage coefficient for snow melt on glaciers', &
           'hours')/=0 ) CALL read_error(1, 'stor_snow')

      ALLOCATE ( Stor_firn(Nhru,MONTHS_PER_YEAR) )
      IF ( declparam(MODNAME, 'stor_firn', 'nhru,nmonths', 'real', &
           '400.0', '150.0', '1000.0', &
           'Monthly Storage coefficient for firn melt on glaciers', &
           'Monthly (January to December) Storage coefficient for firn melt on glaciers', &
           'hours')/=0 ) CALL read_error(1, 'stor_firn')

      ALLOCATE ( Hru_length(Nhru) )
      IF ( declparam(MODNAME, 'hru_length', 'nhru', 'real', &
           '0.0', '0.0', '10000.0', &
           'Length of segment covering all of glacier-possible HRU', &
           'Length of segment covering all of glacier-possible HRU', &
           'km')/=0 ) CALL read_error(1, 'hru_length')

      ALLOCATE ( Hru_width(Nhru) )
      IF ( declparam(MODNAME, 'hru_width', 'nhru', 'real', &
           '0.0', '0.0', '10000.0', &
           'Width of glacier-possible HRU', &
           'Width of glacier-possible HRU', &
           'km')/=0 ) CALL read_error(1, 'hru_width')

      ALLOCATE ( Abl_elev_range(Nhru) )
      IF ( declparam(MODNAME, 'abl_elev_range', 'nhru', 'real', &
           '1000.0', '0.0', '17000.0', &
           'Average HRU snowfield ablation zones elevation range',  &
           'Average HRU snowfield ablation zones elevation range or ~ median-min elev', &
           'elev_units')/=0 ) CALL read_error(1, 'abl_elev_range')

      END FUNCTION glacrdecl

!***********************************************************************
!     glacrinit - Initialize glacr module - get parameter values
!***********************************************************************
      INTEGER FUNCTION glacrinit()
      USE PRMS_GLACR
      USE PRMS_BASIN, ONLY: Hru_area, Hru_elev_ts, Active_hrus, Hru_route_order, &
     &    Hru_type, Basin_area_inv, Hru_elev_meters
      USE PRMS_FLOWVARS, ONLY: Glacier_frac, Alt_above_ela, Glrette_frac
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam, get_ftnunit, compute_ela_aar
      INTRINSIC :: ABS, SQRT, SNGL, REAL
      EXTERNAL :: read_error, tag_count, sort5, glacr_restart
! Local Variables
      INTEGER :: i, j, ii, jj, o, p, hru_flowline(Nhru), toflowline(Nhru), doela, termh, len_str
      INTEGER :: iwksp(Nhru), is(Nhru), ie(Nhru), n_inline(Nhru), cell_id(Nhru), str_id(Nhru), prev
      INTEGER :: count, prev0
      REAL :: hru_dcum(Nhru), uraw0(Nhru), xraw0(Nhru), xrawterm(Nhru), urawterm(Nhru), slu
      REAL :: fraw0(Nhru), urawt(Nhru), xrawt(Nhru), hrawt(Nhru), divu, wksp(Nhru), urawm(Nhru+2)
      REAL :: rc(Nhru), rd(Nhru), re(Nhru), ra(Nhru), rb(Nhru), ll(Nhru), u2dd(Nhru)
      REAL :: xrawm(Nhru+2), flowsurf_slope(Nhru), urawtop(Nhru), cell_idm(Nhru), str_idm(Nhru)
      REAL :: glacier_frac_use(Nhru)
      DOUBLE PRECISION :: curr_area(Nhru), add_area(Nhru)
!
!***********************************************************************
      glacrinit = 0

      IF ( Init_vars_from_file>0 ) CALL glacr_restart(1)

      IF ( getparam(MODNAME, 'max_gldepth', 1, 'real', Max_gldepth)/=0 ) CALL read_error(2, 'max_gldepth')
      IF ( getparam(MODNAME, 'glacrva_coef', Nhru, 'real', Glacrva_coef)/=0 ) CALL read_error(2, 'glacrva_coef')
      IF ( getparam(MODNAME, 'glacrva_exp', Nhru, 'real', Glacrva_exp)/=0 ) CALL read_error(2, 'glacrva_exp')
      IF ( getparam(MODNAME, 'stor_ice', Nhru*MONTHS_PER_YEAR, 'real', Stor_ice)/=0 ) CALL read_error(2, 'stor_ice')
      IF ( getparam(MODNAME, 'stor_snow', Nhru*MONTHS_PER_YEAR, 'real', Stor_snow)/=0 ) CALL read_error(2, 'stor_snow')
      IF ( getparam(MODNAME, 'stor_firn', Nhru*MONTHS_PER_YEAR, 'real', Stor_firn)/=0 ) CALL read_error(2, 'stor_firn')
      IF ( getparam(MODNAME, 'hru_length', Nhru, 'real', Hru_length)/=0 ) CALL read_error(2, 'hru_length')
      IF ( getparam(MODNAME, 'hru_width', Nhru, 'real', Hru_width)/=0 ) CALL read_error(2, 'hru_width')
      IF ( getparam(MODNAME, 'abl_elev_range', Nhru, 'real', Abl_elev_range)/=0 ) CALL read_error(2, 'abl_elev_range')
      IF ( getparam(MODNAME, 'tohru', Nhru, 'integer', Tohru)/=0 ) CALL read_error(2, 'tohru')
      IF ( getparam(MODNAME, 'hru_slope', Nhru, 'real', Hru_slope)/=0 ) CALL read_error(2, 'hru_slope')
      IF ( Init_vars_from_file==0 ) THEN
        Alt_above_ela = 0.0
        Prev_out = 0.0
        Prev_outi = 0.0
        Prev_area = 0.0D0
        Hru_glres_melt = 0.0
        Gl_top_melt = 0.0
        Glacr_flow = 0.0
        Gl_ice_melt = 0.0
        Hru_mb_yrend = 0.0
        Top = 0
        Term = 0
        Ela = 0
        Order_flowline = 0
        Ode_glacrva_coef = 0.0
        Hru_mb_yrcumul = 0.0D0
        Delta_volyr = 0.0D0
        Gl_mb_yrcumul = 0.0D0
        Gl_mb_cumul = 0.0D0
        Gl_mbc_yrend = 0.0D0
        Hru_slope_ts = Hru_slope
        Basal_elev = Hru_elev_ts ! Hru_elev_ts always set in basin, need in case of restart
        Basal_slope = Hru_slope_ts
        Av_basal_slope = 0.0
        Glacr_elev_init = Hru_elev_ts
        Glacr_slope_init = Hru_slope_ts
        Av_fgrad = 0.0
        Basin_gl_top_melt = 0.0D0
        Basin_gl_top_gain = 0.0D0
        Basin_gl_ice_melt = 0.0D0
        Glnet_ar_delta = 0.0D0
        Ikeep_gl = 0
        Keep_gl = 0.0
        Basin_gl_storage = 0.0D0
        Basin_gl_storstart = 0.0D0
        Basin_gl_storvol = 0.0D0
      ENDIF

      Glac_HRUnum_down = 1 ! 1 is the way Weasel delineation was designed
      ! 1 is terminus is smallest ID and top is largest. IDs are stacked.
      ! 0 is terminus is smallest ID and top is largest. IDs are stacked.

      hru_flowline = 0
      toflowline = 0
      str_idm = 1.0E15
      cell_idm = 1.0E15
      uraw0 = 1.0E15
      xraw0 = 1.0E15
      fraw0 = 1.0E15
      hrawt = 1.0E15
      ra = 1.0E15
      rb = 1.0E15
      rc = 1.0E15
      rd = 1.0E15
      re = 1.0E15
      divu = 1.0E3 !in meters
      glacier_frac_use = 0.0
      count = 0
      DO jj = 1, Active_hrus
        j = Hru_route_order(jj)
        ! fill all glacier capable hrus to get variables calculated for all possible glaciers with correct branching
        IF ( Hru_type(j)==GLACIER ) THEN
          count = 1 !has at least one glacier
          glacier_frac_use(j) = 1.0
          !should be end of extensions or branches-- will fail if don't set up with indices stacked
          IF ( Glac_HRUnum_down==1) THEN
            IF (Tohru(j)/=j-1 ) glacier_frac_use(j) = 0.999
          ELSEIF ( Glac_HRUnum_down==0) THEN
            IF (Tohru(j)/=j+1 ) glacier_frac_use(j) = 0.999
          ENDIF
        ENDIF
      ENDDO

      IF ( count>0 ) THEN
  ! Number the glaciers and tags parts that belong together
        CALL tag_count(1, hru_flowline, toflowline, glacier_frac_use)
        hru_dcum = 0.0
        DO i = 1, Ntp !do for all glacier capable hrus
        ! will add self and everything above so cumulative dist from top of flowline
          hru_dcum(Top(i)) = Hru_length(Top(i))
          prev = Top(i)
          DO WHILE ( Tohru(prev)>0 )
             IF ( hru_flowline(Tohru(prev))==hru_flowline(Top(i)) )  &
       &          hru_dcum(Tohru(prev)) = Hru_length(Tohru(prev)) + hru_dcum(prev) !in km
             prev = Tohru(prev)
          ENDDO
        ENDDO
  ! Read input from GIS and compute depths
        Nhrugl = 0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)
          IF ( Hru_type(j)==GLACIER ) THEN
            Nhrugl = Nhrugl + 1
            cell_idm(Nhrugl) = REAL(j)
            str_idm(Nhrugl) = REAL(hru_flowline(j))
            uraw0(Nhrugl) = Hru_elev_meters(j) !inital Hru_elev in meters
            xraw0(Nhrugl) = hru_dcum(j) - Hru_length(j)*0.5 !in km, put it at middle
            hrawt(Nhrugl) = Hru_width(j) !in km
          ENDIF
        ENDDO
  ! sort in ascending order by str_id
        CALL sort5(Nhrugl, str_idm(1:Nhrugl), cell_idm(1:Nhrugl), uraw0(1:Nhrugl), xraw0(1:Nhrugl), &
       &           hrawt(1:Nhrugl), wksp(1:Nhrugl), iwksp(1:Nhrugl))
        is(1) = 1
        DO j = 1, Ntp - 1
          n_inline(j) = 0
          DO jj = is(j), Nhrugl
            IF ( str_idm(jj)==j ) THEN
              n_inline(j) = n_inline(j) + 1
            ELSEIF ( str_idm(jj)==j+1 ) THEN !got to end of string
              ie(j) = jj - 1  !end of string
              DO i = is(j), ie(j) !divide up
                ra(i-is(j)+1) = xraw0(i)
                rb(i-is(j)+1) = REAL(str_idm(i))
                rc(i-is(j)+1) = REAL(cell_idm(i))
                rd(i-is(j)+1) = uraw0(i)
                re(i-is(j)+1) = hrawt(i)
              ENDDO
  ! sort in ascending order by xraw0
              CALL sort5(n_inline(j), ra(1:n_inline(j)), rb(1:n_inline(j)), rc(1:n_inline(j)), &
       &                 rd(1:n_inline(j)), re(1:n_inline(j)), wksp(1:n_inline(j)), iwksp(1:n_inline(j)))
              DO i = is(j), ie(j) !put back in
                xraw0(i) = ra(i-is(j)+1)
                str_id(i) = INT(rb(i-is(j)+1))
                cell_id(i) = INT(rc(i-is(j)+1))
                uraw0(i) = rd(i-is(j)+1)
                hrawt(i) = re(i-is(j)+1)
              ENDDO
              is(j+1) = jj !start of next string
              EXIT         !go to next
            ENDIF
          ENDDO
        ENDDO
  !
        IF (Ntp>1) n_inline(Ntp) = Nhrugl - ie(Ntp-1)
        IF (Ntp==1) n_inline(Ntp) = Nhrugl
        ie(Ntp) = Nhrugl
        DO i = is(Ntp), ie(Ntp) !divide up
          ra(i-is(Ntp)+1) = xraw0(i)
          rb(i-is(Ntp)+1) = REAL(str_idm(i))
          rc(i-is(Ntp)+1) = REAL(cell_idm(i))
          rd(i-is(Ntp)+1) = uraw0(i)
          re(i-is(Ntp)+1) = hrawt(i)
        ENDDO
  ! sort in ascending order by xraw0
        CALL sort5(n_inline(Ntp), ra(1:n_inline(Ntp)), rb(1:n_inline(Ntp)), rc(1:n_inline(Ntp)), &
       &           rd(1:n_inline(Ntp)), re(1:n_inline(Ntp)), wksp(1:n_inline(Ntp)), iwksp(1:n_inline(Ntp)))
        DO i = is(Ntp), ie(Ntp) !put back in
          xraw0(i) = ra(i-is(Ntp)+1)
          str_id(i) = INT(rb(i-is(Ntp)+1))
          cell_id(i) = INT(rc(i-is(Ntp)+1))
          uraw0(i) = rd(i-is(Ntp)+1)
          hrawt(i) = re(i-is(Ntp)+1)
        ENDDO
  ! make terminus and top xraw and uraw
  ! if Glacier_frac(termh) <1, then essentially flattening terminus of glacier
  !  Hru_elev_meters determined off glacierized and non-glacierized area and no way
  !   of saying how high above bare ground glacier is, so have to flatten
        DO j = 1, Ntp
          IF ( toflowline(j)==0 ) THEN !terminus is end of last Hru
            xrawterm(j) = xraw0(ie(j))+Hru_length(cell_id(ie(j)))*0.5 !in km
            slu = (uraw0(ie(j))-uraw0(ie(j)-1))/(xraw0(ie(j))-xraw0(ie(j)-1)) !always at least two Hrus in glacier
            urawterm(j) = (uraw0(ie(j)) + slu*Hru_length(cell_id(ie(j)))*0.5)/divu !in km
          ELSE  !terminus is middle of next Hru, 2 segments away and this is a side branch
            termh = Tohru(cell_id(ie(j)))
            urawterm(j) = Hru_elev_meters(termh)/divu !in km
            xrawterm(j) = xraw0(ie(j)) + Hru_length(cell_id(ie(j)))*0.5 + Hru_width(termh)*0.5 !in km
          ENDIF
        ENDDO
  ! normalize xflow, keep H in km (width with xflow)
        DO i = 1, Nhrugl
          DO j = 1, Ntp
            IF ( i>=is(j) .AND. i<=ie(j) ) THEN
              ll(j) = xrawterm(j)*divu  ! in m
              xrawt(i) = xraw0(i)*divu/ll(j)
              urawt(i) = uraw0(i)/divu-urawterm(j) !in km
            ENDIF
          ENDDO
        ENDDO
  ! normalize urawTop, u scaled at term is always 0, x scaled at top is always 0.
        DO j = 1, Ntp
          slu = (urawt(is(j)+1)-urawt(is(j)))/(xrawt(is(j)+1)-xrawt(is(j))) !always at least two Hrus in glacier
          urawtop(j) = urawt(is(j)) - slu*xrawt(is(j)) !in km
        ENDDO
  ! recalculate hru_slope same way as basal_slope, won't agree with GIS slope
        DO j = 1, Ntp
          len_str = ie(j) - is(j) + 1
          DO i = 1, len_str
            urawm(i) = uraw0(is(j)+i-1)! unscaled in m
            xrawm(i) = xraw0(is(j)+i-1)*divu! unscaled in m
          ENDDO
          xrawm(len_str+2) = 0.0
          xrawm(len_str+1) = ll(j)
          urawm(len_str+2) = (urawtop(j)+urawterm(j))*divu !in m, unscaled
          urawm(len_str+1) = urawterm(j)*divu !in m
          u2dd(1)= (urawm(2)-urawm(len_str+2))/(xrawm(2)-xrawm(len_str+2))
          flowSurf_slope(is(j)) = u2dd(1)
          DO i = 2, len_str
            u2dd(i)= (urawm(i+1)-urawm(i-1))/(xrawm(i+1)-xrawm(i-1))
            flowSurf_slope(is(j)+i-1) = u2dd(i)
          ENDDO
        ENDDO
        DO i = 1, Nhrugl
          Keep_gl(i,4) = flowsurf_slope(i) !initiate, if don't do bottom calculations will stay here
          Hru_slope_ts(cell_id(i)) = ABS(flowsurf_slope(i)) !always positive by definition
        ENDDO
        Glacr_slope_init = Hru_slope_ts
        Basal_slope = Hru_slope_ts
  ! Keep stuff
        DO i = 1, Nhrugl
          Keep_gl(i,1) = urawt(i)
          Keep_gl(i,2) = xrawt(i)
          Keep_gl(i,3) = hrawt(i)
          Ikeep_gl(i,1) = cell_id(i)
          Ikeep_gl(i,2) = str_id(i)
        ENDDO
        DO i = 1, Ntp
          Keep_gl(i,5) = urawtop(i)
          Keep_gl(i,6) = urawterm(i)
          Keep_gl(i,7) = xrawterm(i)
          Ikeep_gl(i,3) = is(i)
          Ikeep_gl(i,4) = ie(i)
        ENDDO
        glacier_frac_use = 0.0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)
          Hru_area_inch2(j) = Hru_area(j)*Acre_inch2
          IF ( Hru_type(j)==GLACIER ) THEN
            glacier_frac_use(j)= Glacier_frac(j)
            !should be end of extensions or branches-- will fail if don't set up with indices stacked
            ! making it so has no connected branches because branching bottom calculations don't work
            IF ( Glac_HRUnum_down==1) THEN
              IF (Tohru(j)/=j-1 .AND. glacier_frac_use(j)==1.0 ) glacier_frac_use(j) = 0.999
            ELSEIF ( Glac_HRUnum_down==0) THEN
              IF (Tohru(j)/=j+1 .AND. glacier_frac_use(j)==1.0 ) glacier_frac_use(j) = 0.999
            ENDIF
          ENDIF
        ENDDO
        CALL tag_count(0, hru_flowline, toflowline, glacier_frac_use)
        !CALL tag_count(1, hru_flowline, toflowline, Glacier_frac)
  ! compute area at start
        curr_area = 0.0D0
        add_area = 0.0D0
        DO i = 1, Ntp !do for all glacier capable hrus
        ! will add self and everything above so cumulative area from top of flowline
          curr_area(Top(i)) = DBLE(Glacier_frac(Top(i)))*Hru_area_inch2(Top(i))
          prev = Top(i)
          DO WHILE ( Tohru(prev)>0 )
            IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
            !should be end of extensions or branches-- will fail if Weasel doesn't set up this way,
            ! and then should go off area of branch
            ! making it so has no connected branches because branching bottom calculations don't work
              IF ( Tohru(prev)==prev-1 ) THEN
                curr_area(Tohru(prev)) = DBLE(Glacier_frac(Tohru(prev)))*Hru_area_inch2(Tohru(prev)) &
       &                                 + curr_area(prev)
                prev = Tohru(prev)
              ELSE !a branch join
                prev0 = prev
                DO WHILE ( Tohru(prev)>0 )
                  IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
                    add_area(Tohru(prev)) = curr_area(prev0)
                    prev = Tohru(prev)
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
            ELSE
              EXIT
            ENDIF
          ENDDO
        ENDDO
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)
          IF ( Hru_type(j)==GLACIER ) THEN
            curr_area(j) = curr_area(j) + add_area(j)
            Prev_area(j) = curr_area(j) !need this for ela calcs, actually is current area
          ENDIF
        ENDDO
        !DO j = 1,Nhru
        !    write(9,*) j,',',Glacr_tag(j)
        !ENDDO
        Gl_area = 0.D0
        Basin_gl_area = 0.D0

        DO o = 1, Ngl
          p = Glacr_tag(Term(o)) !index by Glacr_tag
          Basin_gl_area = Basin_gl_area + curr_area(Term(o))
          Gl_area(p) = curr_area(Term(o))/Acre_inch2
          !print*, 'Glacr_tag', p, ', area acres branches=', Gl_area(p), ', terminus HRU=', Term(o)
        ENDDO
        DO i = 1, Active_hrus
          j = Hru_route_order(i)
          IF ( Hru_type(j)==LAND ) Basin_gl_area = Basin_gl_area + DBLE(Glrette_frac(j))*Hru_area_inch2(j)
        ENDDO
  !
        doela = compute_ela_aar() !no previous years MB, get ELA from AAR ratio, need Prev_area
        DO ii = 1, Ntp
          DO i = 1, Active_hrus
            j = Hru_route_order(i)
            IF ( Hru_type(j)==GLACIER ) THEN
              IF ( Top_tag(j)==Top_tag(Top(ii)) ) Alt_above_ela(j) = Hru_elev_ts(j)- Hru_elev_ts(Ela(ii))
            ENDIF
          ENDDO
        ENDDO
  !******Compute basin weighted averages
  ! Basin_area_inv is in 1/acres, Basin_gl_area in inches squared
        Basin_gl_area = (Basin_gl_area/Acre_inch2)*Basin_area_inv
        !print*, 'Basin area acres=', 1.0/Basin_area_inv
      ENDIF ! skip all if no glaciers
!

      END FUNCTION glacrinit

!***********************************************************************
!     glacrrun - Computes surface runoff using contributing area
!                  computations
!***********************************************************************
      INTEGER FUNCTION glacrrun()
      USE PRMS_GLACR
      USE PRMS_BASIN, ONLY: Hru_elev_ts, Active_hrus, Hru_route_order, Hru_type, &
     &                      Elev_units, Hru_elev_feet, Hru_elev_meters
      USE PRMS_FLOWVARS, ONLY: Alt_above_ela, Glrette_frac
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: comp_glsurf, recompute_soltab
! Local Variables
      INTEGER :: dosol, i, j, count
!***********************************************************************
      glacrrun = 0
      count = 0
!
      DO j = 1, Active_hrus
        i = Hru_route_order(j)
        IF ( Hru_type(i)==LAND ) THEN
          IF (Glrette_frac(j)>NEARZERO) THEN
            count=1 !has at least one snowfield
            EXIT
          ENDIF
        ENDIF
      ENDDO

      IF( Ngl==0 ) THEN !no more glaciers, will first happen at 10/1 when size changes
        Gl_area = 0.D0
        Glnet_ar_delta = 0.0D0
        Hru_glres_melt = 0.0
        Gl_top_melt = 0.0
        Gl_ice_melt = 0.0
        Hru_mb_yrend = 0.0
        Hru_mb_yrcumul = 0.0D0
        Delta_volyr = 0.0D0
        Gl_mb_yrcumul = 0.0D0
        Gl_mb_cumul = 0.0D0
        Gl_mbc_yrend = 0.0D0
        Av_basal_slope = 0.0
        Av_fgrad = 0.0
        Hru_slope_ts = Basal_slope
        Alt_above_ela = 0.0 ! doesn't matter if no glaciers
        dosol = recompute_soltab() ! change soltab tables for Hru_slope_ts
        IF (count==0) THEN
          Glacr_flow = 0.0
          Basin_gl_area = 0.D0 !no snowfields either
          Basin_gl_top_melt = 0.0D0
          Basin_gl_top_gain = 0.0D0
          Basin_gl_ice_melt = 0.0D0
          Basin_gl_storage = 0.0D0
          Basin_gl_storstart = 0.0D0
          Basin_gl_storvol = 0.0D0
        ENDIF
        glacrrun = comp_glsurf(0,count) !call with no glaciers
      ELSE ! have glaciers
        glacrrun = comp_glsurf(1,count)
      ENDIF

      ! reset hru_elev variables for glacier HRUs
      DO j = 1, Active_hrus
        i = Hru_route_order(j)
        IF ( Hru_type(i)==GLACIER ) THEN
          IF ( Ngl==0 ) Hru_elev_ts(i) = Basal_elev(i)
          IF ( Elev_units==FEET ) THEN
            Hru_elev_feet(i) = Hru_elev_ts(i)
            Hru_elev_meters(i) = Hru_elev_ts(i)*FEET2METERS
          ELSE
            Hru_elev_meters(i) = Hru_elev_ts(i)
            Hru_elev_feet(i) = Hru_elev_ts(i)*METERS2FEET
          ENDIF
        ENDIF
      ENDDO
!
      END FUNCTION glacrrun

!***********************************************************************
!     function comp_glsurf - Computes surface runoff using contributing area
!                   computations
!***********************************************************************
      INTEGER FUNCTION comp_glsurf(glacr_exist, glrette_exist)
      USE PRMS_GLACR
      USE PRMS_MODULE, ONLY: Nhru, Start_year
      USE PRMS_BASIN, ONLY: Hru_type, Hru_elev_ts, Basin_area_inv, Active_hrus, &
     &    Hru_route_order, Elev_units, Hru_elev
      USE PRMS_SET_TIME, ONLY: Nowyear, Nowmonth, Julwater
      USE PRMS_INTCP, ONLY: Net_rain, Net_snow
      USE PRMS_SNOW, ONLY: Snowcov_area, Snowmelt, Glacrmelt, Glacr_air_deltemp, Glacr_delsnow, &
     &    Glrette_frac_init, Snowcov_area, Basin_snowicecov, Snow_evap, Glacr_evap, Basin_glacrb_melt
      USE PRMS_FLOWVARS, ONLY: Glacier_frac, Alt_above_ela, Glrette_frac
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: get_ftnunit, compute_ela_mb, compute_ela_aar, recompute_soltab
      INTRINSIC :: ABS, EXP, SUM, SQRT, ISNAN, SNGL, DBLE
      EXTERNAL :: tag_count, bottom
! Local Variables
      INTEGER :: i, j, ii, jj, o, p, stact_hrus, endact_hrus, next, curr, keep, doela, dosol
      INTEGER :: thecase(Nhru), toflowline(Nhru), count_delta(Nhru), lowpt(Nhru), lowest(Nhru)
      INTEGER :: cell_id(Nhru), gl_top(Nhru), hru_flowline(Nhru), gln, topn, count_delta2(Nhru)
      INTEGER :: done, gt(Nhru), oldlow, prev, prev0, dobot, botwrite
      REAL :: stor, remain, gl_snow, glacrold, ca, ode_area(Nhru), av_elev(Nhru), gl_evap
      REAL :: fraw0(Nhru), ela_elevt(Nhru), volm(Nhru), aream(Nhru), frawt(Nhru), divu
      REAL :: slope(Nhru), ode_vol(Nhru), flow_slope(Nhru), glacier_frac_use(Nhru), glacier_fracp(Nhru)
      DOUBLE PRECISION :: gl_total(Nhru), tot_reserv(3), delta_vol(Nhru), extra_vol, gl_gain(Nhru)
      DOUBLE PRECISION :: delta_areayr(Nhru), curr_area(Nhru), volresv, volresv_ice, add_areap(Nhru)
      DOUBLE PRECISION :: in_top_melt_tot(3), in_top_melt(3, Nhru), tot_delta_mb(Nhru), add_area(Nhru)
      DOUBLE PRECISION :: tot_reservi(3),in_top_melt_itot(3), in_top_melt_ice(3, Nhru), curr_areap(Nhru)
! Arguments
      INTEGER, INTENT(IN) :: glacr_exist, glrette_exist
!***********************************************************************
      comp_glsurf = 1
      dobot = 1 ! 1 calls bottom calcs, 0 doesn't: Set to 0 for calibrating, then run one extra step with it on
      ! Should change so that saves the basal elevations (or reads in as parameter) and then recalibrating does not change
      botwrite = 0 ! 1 writes bottom calcs, 0 doesn't: Set to 0 for calibrating
! initialize
      ela_elevt = 0.0
      gt = 0
      gl_top = 0
      lowest = 0
      count_delta = 0
      count_delta2 = 0
      delta_areayr = 0.0
      delta_vol = 0.0D0
      lowpt = 0
      stact_hrus = 0
      endact_hrus = 0
      av_elev = 0.0
      flow_slope = 0.0
      ode_area = 0.0
      ode_vol = 0.0
      slope = 0.0
      volm = 0.0
      aream = 0.0
      divu = 1.0E3 !in meters
      Basin_gl_top_melt = 0.0D0
      Basin_gl_top_gain = 0.0D0
      Basin_gl_ice_melt = 0.0D0
      Hru_glres_melt = 0.0
      Glrette_melt = 0.0
      Gl_top_melt = 0.0
      Gl_ice_melt = 0.0
      Glacr_flow = 0.0
      Basin_snowicecov = Basin_snowicecov*Acre_inch2/Basin_area_inv
      gl_gain = 0.D0
!
! Start of year calculations after have a year of data
      IF ( Julwater==1 .AND. Nowyear>=Start_year+1) THEN
        IF (glacr_exist==1 ) THEN !have glaciers
! Save year ending values from previous year
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            IF ( Hru_type(j)==GLACIER ) THEN
              IF ( Glacier_frac(j)>NEARZERO ) THEN
                Hru_mb_yrend(j) = SNGL(Hru_mb_yrcumul(j))
              ENDIF
            ENDIF
          ENDDO
! ELA calculations
          IF ( MBinit_flag==2 ) THEN
            doela = compute_ela_aar() !want steady state ELA estimation for fraw calc
            DO j = 1, Ntp
              ela_elevt(j)=Hru_elev(Ela(j)) !will scale inside subroutine, want initial one without _ts
              IF ( Elev_units==FEET ) ela_elevt(j) = ela_elevt(j)*FEET2METERS !put in meters
            ENDDO
          ENDIF
          doela = compute_ela_mb() !Regular ELA estimation for daily reservoir calcs
          DO i = 1, Nhrugl
            cell_id(i) = Ikeep_gl(i,1)
            fraw0(i) = Hru_mb_yrend(cell_id(i))/39.37 !inches to m
            frawt(i) = fraw0(i)/divu !in km
          ENDDO
! Have a year of mass balance, compute bottom so can do Hru_elev_ts
          IF ( Nowyear==Start_year+1) THEN ! do bottom calcs and get ca
            ALL_unit = get_ftnunit(735)
            Output_unit = get_ftnunit(All_unit)
            IF (botwrite==1) OPEN ( Output_unit, FILE='output.dat' )
! get flowline info
            glacier_frac_use = 0.0
            DO jj = 1, Active_hrus
              j = Hru_route_order(jj)
              IF ( Hru_type(j)==GLACIER ) THEN
                glacier_frac_use(j)= Glacier_frac(j)
                !should be end of extensions or branches-- will fail if don't set up with indices stacked
                ! making it so has no connected branches because branching bottom calculations don't work
                IF ( Glac_HRUnum_down==1) THEN
                  IF (Tohru(j)/=j-1 .AND. glacier_frac_use(j)==1.0 ) glacier_frac_use(j) = 0.999
                ELSEIF ( Glac_HRUnum_down==0) THEN
                  IF (Tohru(j)/=j+1 .AND. glacier_frac_use(j)==1.0 ) glacier_frac_use(j) = 0.999
                ENDIF
              ENDIF
            ENDDO
            CALL tag_count(0, hru_flowline, toflowline, glacier_frac_use)
! order flowline parts of each glacier
            gln = 0 !initalize, end value is Ngl for this timestep
            DO j = 1, Ntp
              IF ( toflowline(j)==0 ) THEN
                gln = gln + 1
                gt(j) = gln
              ELSE
                gt(j) = 0
              ENDIF
            ENDDO
            DO j = 1, Ntp
              IF ( gt(j)/=0 ) THEN !should be gln of these
                done = 0
                DO i = 1, Ntp
                  IF ( toflowline(i)==j .AND. gt(i)==0 ) THEN
                    gt(i) = gt(j)
                  ELSEIF ( gt(i)/=0 ) THEN !how many to go
                    done = done + 1
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
            topn = 0 !initalize, end value is Ntp for this timestep
            DO j = 1, Ntp
              DO i = 1, Ntp
                IF ( j==toflowline(i) ) EXIT !j is not a top stream
              ENDDO
              IF ( i-1==Ntp ) THEN !got to end, j is a top stream
                topn = topn + 1
                gl_top(topn) = j
              ENDIF
            ENDDO
!Find depths and sort these out too, to get av_elev and flow_slope, write stuff for plotting
            IF (dobot ==1) THEN
              CALL bottom(fraw0(1:Nhrugl), gln, gt(1:Ntp), topn, gl_top(1:Ntp), av_elev(1:Nhrugl), &
     &              ela_elevt(1:Ntp), flow_slope(1:Nhrugl), toflowline(1:Ntp), slope(1:Nhrugl), &
     &              ode_area(1:Nhrugl), ode_vol(1:Nhrugl), botwrite)

              IF (botwrite==1) WRITE ( Output_unit, '(A5,6A13)' ) 'HRU', 'Basal_elev', 'Glacr_elev_i','Hru_elev', &
     &         'Basal_slope', 'Hru_slpe_ts', 'Glacr_slpe_i'
              Basin_gl_storstart = 0.D0
              DO i = 1, Nhrugl
                Keep_gl(i,4) = slope(i) !likely negative, assuming glacier HRU downhill
                IF ( Elev_units==FEET ) THEN
                  Basal_elev(cell_id(i)) = av_elev(i)/FEET2METERS !'av_elev in meters, want in elev_units
                ELSE
                  Basal_elev(cell_id(i)) = av_elev(i)
                ENDIF
                Basal_slope(cell_id(i)) = flow_slope(i) !always positive by definition
                volm(cell_id(i)) = ode_vol(i)   !indexed by terminus HRU
                aream(cell_id(i)) = ode_area(i) !indexed by terminus HRU
                Basin_gl_storstart = Basin_gl_storstart+ ode_vol(i) !in km^3
              ENDDO
              Basin_gl_storstart =  (Basin_gl_storstart*(39370.1**3.0)/Acre_inch2)*Basin_area_inv
            ELSE
              DO i = 1, Nhrugl
                volm(cell_id(i)) = 0.0   !indexed by terminus HRU
                aream(cell_id(i)) = 0.0  !indexed by terminus HRU
              ENDDO
            ENDIF
            DO i = 1, Nhrugl
! Need to do this so that Hru_elev_ts is actually the same as Hru_elev before melt in terminus
              IF ( Glacier_frac(cell_id(i))>NEARZERO) THEN !only effects terminus
                Glacr_elev_init(cell_id(i)) = (Hru_elev(cell_id(i)) - (1.0-Glacier_frac(cell_id(i))) &
     &                                 *Basal_elev((cell_id(i))))/Glacier_frac(cell_id(i))
                Glacr_slope_init(cell_id(i)) = (Hru_slope_ts(cell_id(i)) - (1.0-Glacier_frac(cell_id(i))) &
     &                                 *Basal_slope((cell_id(i))))/Glacier_frac(cell_id(i))
              ENDIF
              IF (botwrite==1) WRITE ( Output_unit, '(I5,6F13.5)' ) cell_id(i), Basal_elev(cell_id(i)), &
     &          Glacr_elev_init(cell_id(i)), Hru_elev(cell_id(i)), Basal_slope(cell_id(i)), Hru_slope_ts(cell_id(i)), &
     &          Glacr_slope_init(cell_id(i))
            ENDDO
          ENDIF
!
! Do for all years
          CALL yearly_ca_coef(frawt(1:Nhrugl), ela_elevt(1:Ntp))
!
!After first year estimate Glacrva_coef based on ODE volume and area, ca = Vol/(Area**Glacrva_exp)
! ca = (Glacrva_coef*39.37**(3.0-2.0*Glacrva_exp))*(Av_fgrad(p)**(0.2))*(Av_basal_slope(p)**(-0.4))
          IF ( Nowyear==Start_year+1) THEN ! do bottom calcs and get ca
            IF (botwrite==1) WRITE ( Output_unit, '(A6,A11,A15)' ) 'Glacr', 'Term_HRU_i', 'Gl_area_init'
            DO o = 1, Ngl
              p = Glacr_tag(Term(o)) !index by Glacr_tag
              Ode_glacrva_coef(p) = (volm(Term(o))/(aream(Term(o))**Glacrva_exp(Term(o))))* &
     &             (1000.0**(3.0-2.0*Glacrva_exp(Term(o))))/(Av_fgrad(p)**(0.2))/(Av_basal_slope(p)**(-0.4))
              !Note, this is based off the flowline depth, not the average depth, thus the volume overestimated
              ! so it is an upper estimate of the coefficient
              IF (botwrite==1) WRITE ( Output_unit, '(I6,I10,F15.6)' ) p, Term(o), Gl_area(p)
            ENDDO
            CLOSE ( Output_unit )
            ! Set back to real flowlines without dividing branches
            CALL tag_count(0, hru_flowline, toflowline, Glacier_frac)
! compute area at start
            curr_area = 0.0D0
            add_area = 0.0D0
            DO i = 1, Ntp !do for all glacier capable hrus
            ! will add self and everything above so cumulative area from top of flowline
              curr_area(Top(i)) = DBLE(Glacier_frac(Top(i)))*Hru_area_inch2(Top(i))
              prev = Top(i)
              DO WHILE ( Tohru(prev)>0 )
                IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
                !should be end of extensions or branches-- will fail if Weasel doesn't set up this way,
                ! and then should go off area of branch
                ! making it so has no connected branches because branching bottom calculations don't work
                  IF ( Tohru(prev)==prev-1 ) THEN
                    curr_area(Tohru(prev)) = DBLE(Glacier_frac(Tohru(prev)))*Hru_area_inch2(Tohru(prev)) &
       &                                 + curr_area(prev)
                    prev = Tohru(prev)
                  ELSE !a branch join
                    prev0 = prev
                    DO WHILE ( Tohru(prev)>0 )
                      IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
                        add_area(Tohru(prev)) = curr_area(prev0)
                        prev = Tohru(prev)
                      ELSE
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDDO
            DO jj = 1, Active_hrus
              j = Hru_route_order(jj)
              IF ( Hru_type(j)==GLACIER ) THEN
                curr_area(j) = curr_area(j) + add_area(j)
              ENDIF
            ENDDO
            Gl_area = 0.D0
            DO o = 1, Ngl
              p = Glacr_tag(Term(o)) !index by Glacr_tag
              Gl_area(p) = curr_area(Term(o))/Acre_inch2
              !IF (Term(o)== 1) print*,'Calibrate to Glacr_tag', p, 'of' ,Ngl
              !print*, 'Glacr_tag', p, 'of', Ngl,', area acres initial=', Gl_area(p), ', terminus HRU=', Term(o)
            ENDDO
            Glnet_ar_delta=0.D0 !start at beginning
          ENDIF
          DO o = 1, Ngl
            p = Glacr_tag(Term(o)) !index by Glacr_tag
            !for next year, positive mass outside
            Gl_mbc_yrend(p) = Gl_mbc_yrend(p)+Gl_mb_yrcumul(p)
          ENDDO
!
! Do retreat/advance on whole glacier at end of year
! last year's area/volume is previous area
          Prev_area = 0.D0
          add_area = 0.0D0
          DO i = 1, Ntp !do for all glacier capable hrus
          ! will add self and everything above so cumulative area from top of flowline
            Prev_area(Top(i)) = DBLE(Glacier_frac(Top(i)))*Hru_area_inch2(Top(i))
            prev = Top(i)
            DO WHILE ( Tohru(prev)>0 )
              IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
              !should be end of extensions or branches-- will fail if Weasel doesn't set up this way,
              ! and then should go off area of branch
              ! making it so has no connected branches because branching bottom calculations don't work
                IF ( Tohru(prev)==prev-1 ) THEN
                  Prev_area(Tohru(prev)) = DBLE(Glacier_frac(Tohru(prev)))*Hru_area_inch2(Tohru(prev)) &
       &                                 + Prev_area(prev)
                  prev = Tohru(prev)
                ELSE !a branch join
                  prev0 = prev
                  DO WHILE ( Tohru(prev)>0 )
                    IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
                      add_area(Tohru(prev)) = add_area(Tohru(prev))+Prev_area(prev0)
                      prev = Tohru(prev)
                    ELSE
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ELSE
                EXIT
              ENDIF
            ENDDO
          ENDDO
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            IF ( Hru_type(j)==GLACIER ) THEN
              Prev_area(j) = Prev_area(j) + add_area(j)
            ENDIF
          ENDDO
          glacier_fracp = Glacier_frac !save before recalculate
          DO o = 1, Ngl
            p = Glacr_tag(Term(o)) !index by Glacr_tag
            !from Luthi 09, varying every year may cause trouble as shrinks
            !global average=0.28*39.37**(3-2*1.375)=.70 in^.25 from Arendt 2006
            ca = (Glacrva_coef(Term(o))*39.37**(3.0-2.0*Glacrva_exp(Term(o))))*(Av_fgrad(p)**(0.2))*(Av_basal_slope(p)**(-0.4))
            Prev_vol(p) = DBLE(ca)*(Prev_area(Term(o)))**DBLE(Glacrva_exp(Term(o)))
            IF ( Prev_vol(p)+Delta_volyr(p)<DNEARZERO ) THEN !lost whole glacier
              delta_areayr(o) = -Prev_area(Term(o))
            ELSEIF ( Delta_volyr(p)/=0.0D0 ) THEN
              delta_areayr(o) = ((Prev_vol(p)+Delta_volyr(p))/DBLE(ca))**(1.D0/DBLE(Glacrva_exp(Term(o)))) - Prev_area(Term(o))
            ELSEIF ( ABS(Delta_volyr(p))<DNEARZERO ) THEN
              delta_areayr(o) = 0.0D0
            ENDIF
!
            IF ( Delta_volyr(p)<DNEARZERO ) THEN !glacier retreating, do later
              remain = 0.0
            ELSE !glacier advancing
              remain = SNGL(delta_areayr(o))
            ENDIF
            next = 0
            curr = Term(o)
            DO WHILE ( remain>0.0 )
! in advancing glacier, get rid of total snow till furthest possible
!  terminus of glacier (not letting the glaciers combine) THIS WILL BE DICTATED BY THE HRU MAP
              glacrold = Glacier_frac(curr)
              Glacier_frac(curr) = (Glacier_frac(curr)*SNGL(Hru_area_inch2(curr))+remain) &
     &                         /SNGL(Hru_area_inch2(curr)) !all in inches
              IF ( Glacier_frac(curr)>1.0 ) THEN !glacier can't be more than full, look for next to advance in to
                Glacier_frac(curr) = 1.0
                next = Tohru(curr) !find next to expand in to, could be another glacier then will go to its terminus
                IF ( next==0 ) THEN
                  ! can't advance anymore
                  remain = remain - (1.0-glacrold)*SNGL(Hru_area_inch2(curr))
                  IF ( remain/Hru_area_inch2(curr)<=NEARZERO ) remain = 0.0
                  extra_vol = Prev_vol(Term(o)) + Delta_volyr(p) - &
     &              DBLE(ca)*(Prev_area(Term(o))+delta_areayr(o)-DBLE(remain))**DBLE(Glacrva_exp(Term(o)))
                  !amount not added to area that should, use to thicken later
                  EXIT
                ENDIF
              ENDIF
              remain = remain - (Glacier_frac(curr)-glacrold)*SNGL(Hru_area_inch2(curr))
              IF ( remain/Hru_area_inch2(curr)<=NEARZERO ) THEN !limit of accuracy for reals
                remain = 0.0
              ENDIF
              curr = next
            ENDDO
!
            IF ( Delta_volyr(p)<DNEARZERO ) THEN  !glacier retreating
              remain = SNGL(delta_areayr(o))
            ELSE !glacier advancing, did already
              remain = 0.0
            ENDIF
            lowpt(o) = Term(o)
            DO WHILE ( remain<0.0 )
! in retreating glacier, get rid of total melt till top of glacier.
              IF ( Glacier_frac(lowpt(o))==0.0 ) THEN
                oldlow = lowpt(o)
                DO i = 1, Active_hrus
                  j = Hru_route_order(i)
                  IF ( Hru_type(j)==GLACIER ) THEN
                    ! find a j in glacier
                    IF ( Tohru(j)==oldlow ) THEN
                      lowpt(o) = j
                      EXIT
! if pass a confluence on this iteration, will melt branches cyclic fashion, off of HRU numerical order
!  not in fractional to each or by branch---SHOULD FIX THIS?
                    ENDIF
                  ENDIF
                ENDDO
                IF ( lowpt(o)==oldlow ) EXIT !got to end of glacier since terminus is also top of glacier, so all gone
              ENDIF
              glacrold = Glacier_frac(lowpt(o))
              Glacier_frac(lowpt(o)) = (glacrold*SNGL(Hru_area_inch2(lowpt(o)))+remain) &
     &                              /SNGL(Hru_area_inch2(lowpt(o)))
              IF ( Glacier_frac(lowpt(o))<NEARZERO ) Glacier_frac(lowpt(o)) = 0.0
              !glacier is gone in this hru, look for next glacier-full hru
              remain = remain + (glacrold-Glacier_frac(lowpt(o)))*SNGL(Hru_area_inch2(lowpt(o)))
              IF ( remain/Hru_area_inch2(lowpt(o))>-NEARZERO ) THEN !limit of accuracy for reals
                remain = 0.0
              ENDIF
            ENDDO
          ENDDO !get out of the loop over all the glaciers
! clean up and compute area before start of year
          DO i = 1, Active_hrus
            j = Hru_route_order(i)
            IF ( Hru_type(j)==GLACIER ) THEN
              IF ( Glacier_frac(j)<NEARZERO ) Glacier_frac(j)=0.0
              IF ( Glacier_frac(j)>1.0-NEARZERO ) Glacier_frac(j)=1.0
            ELSE
              Glacier_frac(j)=0.0
            ENDIF
          ENDDO
! Have to renumber the glaciers and tags incase lost some
          CALL tag_count(0, hru_flowline, toflowline, Glacier_frac)
          ! compute new area
          curr_area = 0.0D0
          add_area = 0.0D0
          curr_areap = 0.0D0
          add_areap = 0.0D0
          DO i = 1, Ntp !do for all glacier capable hrus
          ! will add self and everything above so cumulative area from top of flowline
            curr_area(Top(i)) = DBLE(Glacier_frac(Top(i)))*Hru_area_inch2(Top(i))
            curr_areap(Top(i)) = DBLE(glacier_fracp(Top(i)))*Hru_area_inch2(Top(i))
            prev = Top(i)
            DO WHILE ( Tohru(prev)>0 )
              IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
              !should be end of extensions or branches-- will fail if Weasel doesn't set up this way,
              ! and then should go off area of branch
              ! making it so has no connected branches because branching bottom calculations don't work
                IF ( Tohru(prev)==prev-1 ) THEN
                  curr_area(Tohru(prev)) = DBLE(Glacier_frac(Tohru(prev)))*Hru_area_inch2(Tohru(prev)) &
       &                                 + curr_area(prev)
                  curr_areap(Tohru(prev)) = DBLE(glacier_fracp(Tohru(prev)))*Hru_area_inch2(Tohru(prev)) &
       &                                 + curr_areap(prev)
                  prev = Tohru(prev)
                ELSE !a branch join
                  prev0 = prev
                  DO WHILE ( Tohru(prev)>0 )
                    IF ( Glacr_tag(Tohru(prev))==Glacr_tag(Top(i)) ) THEN
                      add_area(Tohru(prev)) = add_area(Tohru(prev))+curr_area(prev0)
                      add_areap(Tohru(prev)) = add_areap(Tohru(prev))+curr_areap(prev0)
                      prev = Tohru(prev)
                    ELSE
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ELSE
                EXIT
              ENDIF
            ENDDO
          ENDDO
          DO jj = 1, Active_hrus
            j = Hru_route_order(jj)
            IF ( Hru_type(j)==GLACIER ) THEN
              curr_area(j) = curr_area(j) + add_area(j)
              curr_areap(j) = curr_areap(j) + add_areap(j)
            ENDIF
          ENDDO
          Gl_area = 0.D0
          Basin_gl_area = 0.D0
          DO o = 1, Ngl
            p = Glacr_tag(Term(o)) !index by Glacr_tag
            Glnet_ar_delta(p) = Glnet_ar_delta(p) + ( curr_area(Term(o)) - curr_areap(Term(o)) )/Acre_inch2
            Gl_area(p) = curr_area(Term(o))/Acre_inch2
            Basin_gl_area = Basin_gl_area + curr_area(Term(o))
            !print*, 'Glacr_tag', p, ', area acres', Nowyear,' =', Gl_area(p), ', terminus HRU=', Term(o)
          ENDDO
! have bottom, compute Hru_elev_ts related stuff
          DO ii = 1, Ntp
            DO i = 1, Active_hrus
              j = Hru_route_order(i)
              IF ( Hru_type(j)==GLACIER ) THEN
                IF ( Top_tag(j)==Top_tag(Top(ii)) ) Alt_above_ela(j) = Hru_elev_ts(j)- Hru_elev_ts(Ela(ii))
                Hru_elev_ts(j) = Glacier_frac(j)*Glacr_elev_init(j) + (1.0-Glacier_frac(j))*Basal_elev(j)
                Hru_slope_ts(j) = Glacier_frac(j)*Glacr_slope_init(j) + (1.0-Glacier_frac(j))*Basal_slope(j)
              ENDIF
            ENDDO
          ENDDO
          dosol = recompute_soltab() ! change soltab tables for Hru_slope_ts
! Beginning of year zero out
          Hru_mb_yrcumul = 0.0D0
          Delta_volyr = 0.0D0
! Clean stuff if gone
          DO i = 1, Active_hrus
            j = Hru_route_order(i)
            IF ( Hru_type(j)==GLACIER ) THEN
              IF ( Glacier_frac(j)==0.0 ) THEN
                Hru_mb_yrend(j) = 0.0
                Hru_mb_yrcumul(j) = 0.0D0
                Alt_above_ela(j) = 0.0
              ENDIF
            ENDIF
          ENDDO
          DO p = 1, Nhru !because could be left over more than Ngl of these
            IF ( Gl_area(p)==0.D0 ) THEN
              Glnet_ar_delta(p) = 0.0D0
              Delta_volyr(p) = 0.0D0
              Gl_mb_yrcumul(p) = 0.0D0
              Gl_mb_cumul(p) = 0.0D0
              Gl_mbc_yrend(p) = 0.0D0
              Av_basal_slope(p) = 0.0
              Av_fgrad(p) = 0.0
            ENDIF
          ENDDO
        ENDIF
!Snowfield area change uses Baumann and Winkler 2010 to change area every 10 years;
! technically each snowfield should have own ablation elevation range.
        IF (glrette_exist==1) THEN !have snowfields,
          IF ( MOD(Nowyear-Start_year,10)==0 ) THEN !change them
            DO i = 1, Active_hrus
              j = Hru_route_order(i)
              IF ( Hru_type(j)==LAND .AND. Glrette_frac(j)>NEARZERO) THEN
                IF ( Elev_units==FEET ) Glrette_frac(j) = ( METERS2FEET*(45.7*Glacr_air_deltemp(j) &
     &               -12.0*Glacr_delsnow(j))/Abl_elev_range(j) +1.0 )*Glrette_frac_init(j)
                IF ( Elev_units==METERS ) Glrette_frac(j) = ( (45.7*Glacr_air_deltemp(j) &
     &               -12.0*Glacr_delsnow(j))/Abl_elev_range(j) +1.0 )*Glrette_frac_init(j)
                IF ( Glrette_frac(j)<0.0 ) Glrette_frac(j)=0.0
                IF ( Glrette_frac(j)>1.0 ) Glrette_frac(j)=1.0
              ENDIF
            ENDDO
          ENDIF
          DO i = 1, Active_hrus !every year
            j = Hru_route_order(i)
            Basin_gl_area = Basin_gl_area + DBLE(Glrette_frac(j))*Hru_area_inch2(j) !keep in inches
          ENDDO
        ENDIF

!
!******Compute basin weighted averages (to units of fraction area)
! Basin_area_inv is in 1/acres, Basin_area_inv is in 1/acres, Basin_gl_area in inches squared
        Basin_gl_area = (Basin_gl_area/Acre_inch2)*Basin_area_inv
      ENDIF !get out of start-of-year computations

!
! Melt runoff calculations, every day
      DO i = 1, Active_hrus
        j = Hru_route_order(i)
        Hru_glres_melt(j) = 0.0
        gl_total(j) = 0.0D0
        IF ( Hru_type(j)==GLACIER ) THEN
          !melting ice + melting snow (energy model), *area = volume
          IF ( Glacier_frac(j)>NEARZERO ) THEN
            Hru_glres_melt(j) = Glacier_frac(j)*(Snowmelt(j) + Glacrmelt(j)/Glacier_frac(j))
            Snowmelt(j) = (1.0 - Glacier_frac(j))*Snowmelt(j) ! this is the snowmelt that is not routed through glacier
            ! all excess rain is included in melt so glacier melt is ( glacier_frac*(snowmelt+net_rain)+ glacrmelt )*hru_area
            gl_snow = Glacier_frac(j)*(Net_rain(j)+Net_Snow(j)) !Pk_precip is zero if no snow, so don't use
            gl_evap = Glacier_frac(j)*(Snow_evap(j) + Glacr_evap(j)/Glacier_frac(j))
            gl_gain(j) = DBLE(gl_snow - gl_evap)
            gl_total(j) = -Hru_glres_melt(j) + gl_gain(j)
            !this is daily mass balance on glacier part of HRU in inches, divide by glacier_frac so averaged over glaciated part of HRU only
            Hru_mb_yrcumul(j) = Hru_mb_yrcumul(j) + gl_total(j)/Glacier_frac(j)
            Basin_gl_top_gain = Basin_gl_top_gain + gl_gain(j)*Hru_area_inch2(j)
            !postive indicates snow, negative indicates melt
          ENDIF
        ENDIF
        IF ( Hru_type(j)==LAND ) THEN
          Glrette_melt(j) = 0.0
          !melting ice + melting snow (energy model), *area = volume
          IF ( Glrette_frac(j)>NEARZERO ) THEN
            Glrette_melt(j) = Glrette_frac(j)*(Snowmelt(j) + Glacrmelt(j)/Glrette_frac(j))
            Snowmelt(j) = (1.0 - Glrette_frac(j))*Snowmelt(j) ! this is the snowmelt that is not included in glacierette melt
            gl_snow = Glrette_frac(j)*(Net_rain(j)+Net_Snow(j)) !Pk_precip is zero if no snow, so don't use
            gl_evap = Glrette_frac(j)*(Snow_evap(j) + Glacr_evap(j)/Glrette_frac(j))
            gl_gain(j) = DBLE(gl_snow - gl_evap)
          ENDIF
        ENDIF
      ENDDO

      in_top_melt = 0.0D0 !initialize
      in_top_melt_ice = 0.0D0 !initialize
      IF (glacr_exist==1) THEN
        DO o = 1, Ngl
          p = Glacr_tag(Term(o)) !index by Glacr_tag
! This code will force at least top of each branch to be above or at ELA
          DO ii = 1, Ntp
            thecase(ii) = 0
            ! find a branch in glacier
            IF ( Glacr_tag(Top(ii))/=p ) CYCLE
            keep = Top(ii)
            DO i = 1, Active_hrus
              j = Hru_route_order(i)
              IF ( Hru_type(j)==GLACIER ) THEN
                !find lowest in branch
                IF ( Glacr_tag(j)==p .AND. &
     &                  (Top_tag(j)==Top_tag(Top(ii)) .OR. (Top_tag(j)==-1.AND.count_delta2(j)==0)) ) THEN
                  count_delta2(j) = 1 !for ones that are end of two branches, will be made part of first branch
                  IF ( Hru_elev_ts(j)<Hru_elev_ts(keep) ) keep = j
                ENDIF
              ENDIF
            ENDDO
            lowest(ii) = keep

            IF ( lowest(ii)==Top(ii) ) THEN !glacier branch only 1 HRU
              thecase(ii) = 1
            ELSEIF ( lowest(ii)==Ela(ii) .OR. Snowcov_area(lowest(ii))>NEARZERO ) THEN !no branch abl.zone
              thecase(ii) = 2
            ENDIF

! Now find reservoirs in each branch: three-- firn, snow, and ice
            DO jj = 1, 3 !firn, snow, ice
              IF ( jj==1 ) THEN !accumulation zone
                IF ( Top(ii)==Ela(ii) .AND. thecase(ii)/=1 ) CYCLE !has no acc. zone
                stact_hrus = Top(ii)
                keep = Top(ii)
                DO i = 1, Active_hrus
                  j = Hru_route_order(i)
                  IF ( Hru_type(j)==GLACIER ) THEN
                    !find lowest in branch above ELA-- there will be one
                    IF ( Top_tag(j)==Top_tag(Top(ii)) .AND.            &
     &                    Hru_elev_ts(j)<Hru_elev_ts(keep) .AND.       &
     &                    Hru_elev_ts(j)>Hru_elev_ts(Ela(ii)) ) keep = j
                  ENDIF
                ENDDO
                endact_hrus = keep
                IF ( thecase(ii)==1 ) endact_hrus = lowest(ii)
              ELSEIF ( jj==2 ) THEN !middle zone
                stact_hrus = Ela(ii)
                keep = Ela(ii)
                DO i = 1, Active_hrus
                  j = Hru_route_order(i)
                  IF ( Hru_type(j)==GLACIER ) THEN
                    !find lowest in branch with snow cover below ELA
                    IF ( Top_tag(j)==Top_tag(Top(ii)) .AND.         &
     &                    Snowcov_area(j)>NEARZERO .AND.        &
     &                    Hru_elev_ts(j)<Hru_elev_ts(keep) ) keep = j
                  ENDIF
                ENDDO
                endact_hrus = keep
                IF ( thecase(ii)==2 ) endact_hrus = lowest(ii)
              ELSEIF ( jj==3 ) THEN !ablation zone
                IF ( thecase(ii)==2 ) CYCLE
                keep = lowest(ii)
                DO i = 1, Active_hrus
                  j = Hru_route_order(i)
                  IF ( Hru_type(j)==GLACIER ) THEN
                    !find highest in branch with no snow cover below ELA
                    IF ( Top_tag(j)==Top_tag(Top(ii)) .AND. Snowcov_area(j)<NEARZERO .AND. &
     &                    Hru_elev_ts(j)>Hru_elev_ts(keep) .AND.                               &
     &                    Hru_elev_ts(j)<Hru_elev_ts(Ela(ii)) ) keep = j
                  ENDIF
                ENDDO
                stact_hrus = keep
                endact_hrus = lowest(ii)
              ENDIF
!
! Now find runoffs for each reservoir of each glacier
              DO i = 1, Active_hrus
                j = Hru_route_order(i)
                volresv = 0.0D0
                volresv_ice = 0.0D0
                IF ( Hru_type(j)==GLACIER ) THEN
                  IF ( Glacr_tag(j)==p .AND. Hru_elev_ts(j)<=Hru_elev_ts(stact_hrus) .AND.  &
     &                   Hru_elev_ts(j)>=Hru_elev_ts(endact_hrus) .AND. &
     &                  (Top_tag(j)==Top_tag(stact_hrus) .OR.           &
     &                  (Top_tag(j)==-1.AND.count_delta(j)==0)) ) THEN
                    count_delta(j) = 1
                    volresv = DBLE(Hru_glres_melt(j))*Hru_area_inch2(j)
                    IF ( volresv>DNEARZERO ) in_top_melt(jj, ii) = in_top_melt(jj, ii)+ volresv
                    ! all excess rain is included in melt, rain on ice goes into reservoirs
                    ! should be true unless Glacrmelt==0
                    IF ( Glacrmelt(j)-Net_rain(j)*Glacier_frac(j)>NEARZERO ) &
     &                     volresv_ice =  DBLE(Glacrmelt(j)-Net_rain(j)*Glacier_frac(j))*Hru_area_inch2(j)
                    IF ( volresv_ice>DNEARZERO ) in_top_melt_ice(jj, ii) = in_top_melt_ice(jj, ii)+ volresv_ice
                    delta_vol(o) = delta_vol(o) + gl_total(j)*Hru_area_inch2(j)/0.917
                    ! divide by density ratio to get in volume, if were all converted to ice (by end of year)
                  ENDIF
                ENDIF
              ENDDO
! skip some reservoirs if glacier branch is too small
              IF ( jj==1 .AND. thecase(ii)==1 ) EXIT
            ENDDO !collected all 3 reservoirs (jj) for this top ii
          ENDDO !collected all branches in glacier
          DO jj = 1, 3
            in_top_melt_tot(jj) = 0.0D0
            in_top_melt_itot(jj) = 0.0D0
            tot_reserv(jj) = 0.0D0
            tot_reservi(jj) = 0.0D0
            stor = 0.0
            IF ( jj==1 ) THEN
              stor = Stor_firn(Term(o),Nowmonth)/24.0 !days
            ELSEIF ( jj==2 ) THEN
              stor = Stor_snow(Term(o),Nowmonth)/24.0 !days
            ELSEIF ( jj==3 ) THEN
              stor = Stor_ice(Term(o),Nowmonth)/24.0  !days
            ENDIF
            DO ii = 1, Ntp
            ! find a branch in glacier
              IF ( Glacr_tag(Top(ii))/=p ) CYCLE
              in_top_melt_tot(jj) = in_top_melt_tot(jj) + in_top_melt(jj, ii)
              in_top_melt_itot(jj) = in_top_melt_itot(jj) + in_top_melt_ice(jj, ii)
              IF ( jj==1 .AND. thecase(ii)==1 ) stor = Stor_snow(Term(o),Nowmonth)/24.0  !days
              IF ( jj==2 .AND. thecase(ii)==2 ) stor = Stor_ice(Term(o),Nowmonth)/24.0  !days
            ENDDO
            tot_reserv(jj) = Prev_out(p, jj)*EXP(-1.0/stor) + in_top_melt_tot(jj)*(1.0-EXP(-1.0/stor))
            tot_reservi(jj) = Prev_outi(p, jj)*EXP(-1.0/stor) + in_top_melt_itot(jj)*(1.0-EXP(-1.0/stor))
            Gl_ice_melt(p) = Gl_ice_melt(p) + SNGL(tot_reservi(jj))
            Gl_top_melt(p) = Gl_top_melt(p) + SNGL(tot_reserv(jj))
            Prev_out(p, jj) = SNGL(in_top_melt_tot(jj))
            Prev_outi(p, jj) = SNGL(in_top_melt_itot(jj))
! On last jj, will be at the terminus and have all out of glacier
! THIS IS THE TOTAL RUNOFF OUT THE ALL THE GLACIERS in inches cubed
            Basin_gl_ice_melt = Basin_gl_ice_melt + tot_reservi(jj)
            Basin_gl_top_melt = Basin_gl_top_melt + tot_reserv(jj)
          ENDDO
          Delta_volyr(p) = Delta_volyr(p) + delta_vol(o)
! only terminus HRUs will have glacier flow, like a flow from a gw reservoir
          Glacr_flow(Term(o)) = Gl_top_melt(p)
        ENDDO
!
! Calculated cumulative mass balance in inches of glaciers to compare to data
! Do off area changed after last day of previous hydrological year
        tot_delta_mb = 0.D0
        DO jj = 1, Active_hrus
          j = Hru_route_order(jj)
          IF ( Hru_type(j)==GLACIER ) THEN
            DO ii = 1, Active_hrus
              i = Hru_route_order(ii)
              IF ( Hru_type(i)==GLACIER ) THEN !find a i in glacier
                IF ( Glacr_tag(i)==Glacr_tag(j) .AND. Hru_elev_ts(i)>=Hru_elev_ts(j) ) THEN
                !will add self (i=j) and everything above
                  tot_delta_mb(j) = tot_delta_mb(j) + Hru_mb_yrcumul(i)*Glacier_frac(i)*Hru_area_inch2(i)
                ENDIF
              ENDIF
            ENDDO
            Basin_snowicecov = Basin_snowicecov + DBLE(( 1.-Snowcov_area(j) )*Glacier_frac(j))*Hru_area_inch2(j)
          ENDIF
        ENDDO
        DO o = 1, Ngl
          p = Glacr_tag(Term(o)) !index by Glacr_tag
          Gl_mb_yrcumul(p) = tot_delta_mb(Term(o))/(Gl_area(p)*Acre_inch2) !if glaciers combine or split will jump to new area and terminus
          Gl_mb_cumul(p) = Gl_mbc_yrend(p)+Gl_mb_yrcumul(p)
        ENDDO
      ENDIF
!
      IF (glrette_exist==1) THEN
        DO i = 1, Active_hrus
          j = Hru_route_order(i)
          IF ( Hru_type(j)==LAND .AND. Glrette_frac(j)>NEARZERO) THEN
            ! all excess rain is included in melt, should be true unless Glacrmelt==0
            IF ( Glacrmelt(j)-Net_rain(j)*Glrette_frac(j)>NEARZERO ) &
     &        Basin_gl_ice_melt = Basin_gl_ice_melt + DBLE(Glacrmelt(j)-Net_rain(j)*Glrette_frac(j))*Hru_area_inch2(j)
            Basin_gl_top_melt = Basin_gl_top_melt + DBLE(Glrette_melt(j))*Hru_area_inch2(j)
            Basin_gl_top_gain = Basin_gl_top_gain + DBLE(gl_gain(j))*Hru_area_inch2(j)
            Basin_snowicecov = Basin_snowicecov + DBLE(( 1.-Snowcov_area(j) )*Glrette_frac(j))*Hru_area_inch2(j)
            Glacr_flow(j) = REAL(Glrette_melt(j)*Hru_area_inch2(j))
          ENDIF
        ENDDO
      ENDIF
!
!******Compute basin weighted averages (to units of inches/dt)
! Basin_area_inv is in 1/acres, Basin_gl_top_* is in inches cubed over all area
!  want in inches over all area
      Basin_gl_ice_melt = (Basin_gl_ice_melt/Acre_inch2)*Basin_area_inv !this will be too small by the amount of rain soaked up by ice
      Basin_gl_top_melt = (Basin_gl_top_melt/Acre_inch2)*Basin_area_inv
      Basin_gl_top_gain = (Basin_gl_top_gain/Acre_inch2)*Basin_area_inv
      Basin_snowicecov = (Basin_snowicecov/Acre_inch2)*Basin_area_inv
! these will be zero when glaciers are gone and negative while receding glaciers
! at zero point could back track it so gl_storage went from postive to zero
      Basin_gl_storage = Basin_gl_storage + Basin_gl_top_gain - Basin_glacrb_melt - Basin_gl_top_melt
      Basin_gl_storvol = Basin_gl_storage/Basin_area_inv
!
      comp_glsurf = 0
!
      END FUNCTION comp_glsurf

!***********************************************************************
!     function compute_ela_mb - Identifies ELA Hru closest to 0 from
!                              last years MB, won't work first year
!***********************************************************************
      INTEGER FUNCTION compute_ela_mb()
      USE PRMS_GLACR, ONLY: Ntp, Ngl, Glacr_tag, Term, Top, Top_tag, Hru_mb_yrend, &
     &    Ela, Nhru, GLACIER
      USE PRMS_BASIN, ONLY: Hru_type, Active_hrus, Hru_route_order
      IMPLICIT NONE
! Functions
      INTRINSIC :: ABS
! Local Variables
      INTEGER :: i, j, ii, o, p, ela2(Nhru)
      REAL :: elamin0(Nhru), elamin20(Nhru), elamin(Nhru), elamin2(Nhru)
!***********************************************************************
      compute_ela_mb = 1
! initialize
      ela2 = 0
      elamin0 = 1.0E15
      elamin20 = 1.0E15
      elamin = 1.0E15
      elamin2 = 1.0E15

      DO o = 1, Ngl
        !find the ELA hru
        p = Glacr_tag(Term(o)) !index by Glacr_tag
        DO i = 1, Active_hrus
          j = Hru_route_order(i)
          IF ( Hru_type(j)/=GLACIER ) CYCLE
          ! find a j in glacier
          IF ( Glacr_tag(j)/=p ) CYCLE
          DO ii = 1, Ntp
            IF ( Top_tag(j)==Top_tag(Top(ii)) ) THEN
            !same glacier branch
              elamin(j) = ABS(Hru_mb_yrend(j))
              IF ( elamin(j)<elamin0(ii) ) THEN
                Ela(ii) = j !j index
                elamin0(ii) = elamin(j)
              ENDIF
            ELSEIF ( Top_tag(j)==-1 ) THEN
            !same glacier confluence branch
                elamin2(j) = ABS(Hru_mb_yrend(j))
                IF ( elamin2(j)<elamin20(ii) ) THEN
                 ela2(ii) = j !j index
                 elamin20(ii) = elamin2(j)
                ENDIF
            ELSE
              CYCLE
            ENDIF
            IF ( elamin20(ii)<elamin0(ii) ) Ela(ii) = ela2(ii)
           !ELA in main branch, divided hrus by elevation
          ENDDO
        ENDDO
      ENDDO
!
      compute_ela_mb = 0
!
      END FUNCTION compute_ela_mb

!***********************************************************************
!     function compute_ela_aar - Identifies ELA Hru by closest to area
!                               from theoretical steady state AAR
!***********************************************************************
      INTEGER FUNCTION compute_ela_aar()
      USE PRMS_GLACR, ONLY: Ntp, Ngl, Glacr_tag, Term, Top, Top_tag, Prev_area, Ela, Nhru, GLACIER
      USE PRMS_BASIN, ONLY: Hru_type, Active_hrus, Hru_route_order
      IMPLICIT NONE
! Functions
      INTRINSIC :: ABS, SNGL
! Local Variables
      INTEGER :: i, j, ii, o, p, ela2(Nhru)
      REAL :: elamin0(Nhru), elamin20(Nhru), elamin(Nhru), elamin2(Nhru)
      REAL :: aar, elaarea
      REAL, PARAMETER :: Convert_units = (39.37*1000)**2
!***********************************************************************
      compute_ela_aar = 1
! initialize
      ela2 = 0
      elamin0 = 1.0E15
      elamin20 = 1.0E15
      elamin = 1.0E15
      elamin2 = 1.0E15
      aar = 0.54 !default value

      DO o = 1, Ngl
        p = Glacr_tag(Term(o)) !index by Glacr_tag
        IF ( Prev_area(Term(o))<1.0*Convert_units ) aar = 0.44
        !for glaciers area <1km^2
        IF ( Prev_area(Term(o))>=4.0*Convert_units ) aar = 0.64
        !for glaciers area >4km^2
        elaarea = SNGL(Prev_area(Term(o)))*aar
!aar is percentage of area from top down, from Kern and Laszlo 2010
!for most glaciers will be 0.64 above ELA, 1 above terminus, 0 above top
!this assumes glacier at steady state
        !find the ELA hru
        DO i = 1, Active_hrus
          j = Hru_route_order(i)
          IF ( Hru_type(j)/=GLACIER ) CYCLE
          ! find a j in glacier
          IF ( Glacr_tag(j)/=p ) CYCLE
          DO ii = 1, Ntp
            IF ( Top_tag(j)==Top_tag(Top(ii)) ) THEN
            !same glacier branch
              elamin(j) = ABS(elaarea-SNGL(Prev_area(j)))
              IF ( elamin(j)<elamin0(ii) ) THEN
                Ela(ii) = j !j index
                elamin0(ii) = elamin(j)
              ENDIF
            ELSEIF ( Top_tag(j)==-1 ) THEN
            !same glacier confluence branch
                elamin2(j) = ABS(elaarea-SNGL(Prev_area(j)))
                IF ( elamin2(j)<elamin20(ii) ) THEN
                 ela2(ii) = j !j index
                 elamin20(ii) = elamin2(j)
                ENDIF
            ELSE
              CYCLE
            ENDIF
            IF ( elamin20(ii)<elamin0(ii) ) Ela(ii) = ela2(ii)
           !ELA in main branch, divided hrus by elevation
          ENDDO
        ENDDO
      ENDDO
!
      compute_ela_aar = 0

      END FUNCTION compute_ela_aar

!***********************************************************************
!     function recompute_soltab - Call the soltab routines again because
! of changing slope
!***********************************************************************
      INTEGER FUNCTION recompute_soltab()
      USE PRMS_GLACR, ONLY: Hru_slope_ts, MAX_DAYS_PER_YEAR, GLACIER
      USE PRMS_SOLTAB, ONLY: Hru_aspect, Hru_cossl, PI, RADIANS, &
     &    Soltab_potsw, Soltab_sunhrs, Solar_declination, &
     &    ECCENTRICY, DEGDAY, DEGDAYRAD
      USE PRMS_BASIN, ONLY: Hru_type, Active_hrus, Hru_route_order, Hru_lat
      IMPLICIT NONE
! Functions
      INTRINSIC :: SIN, COS, FLOAT, SNGL
      EXTERNAL :: compute_soltab
! Local Variables
      INTEGER :: jd, n, nn
      DOUBLE PRECISION :: obliquity(MAX_DAYS_PER_YEAR)
      DOUBLE PRECISION :: y, y2, y3, jddbl
!***********************************************************************
      recompute_soltab = 0
! initialize
      DO jd = 1, MAX_DAYS_PER_YEAR
        jddbl = DBLE(jd)
        obliquity(jd) = 1.0D0 - (ECCENTRICY*COS((jddbl-3.0D0)*DEGDAYRAD))
        y = DEGDAYRAD*(jddbl-1.0D0) ! assume noon
        y2 = 2.0D0*y
        y3 = 3.0D0*y
        Solar_declination(jd) = 0.006918D0 - 0.399912D0*COS(y) + 0.070257D0*SIN(y) &
     &                          - 0.006758D0*COS(y2) + 0.000907D0*SIN(y2) &
     &                          - 0.002697D0*COS(y3) + 0.00148D0*SIN(y3)
      ENDDO
!   Module Variables
      DO nn = 1, Active_hrus
        n = Hru_route_order(nn)
        IF ( Hru_type(n)==GLACIER ) THEN !only call if glacier HRU and could have changed
          Soltab_sunhrs(1, n) = 0.0
          Soltab_potsw(1, n) = 0.0
          CALL compute_soltab(obliquity, Solar_declination, Hru_slope_ts(n), Hru_aspect(n), &
     &                      Hru_lat(n), Hru_cossl(n), Soltab_potsw(1, n), &
     &                      Soltab_sunhrs(1, n), Hru_type(n), n)
          ENDIF
      ENDDO
!
      END FUNCTION recompute_soltab

!***********************************************************************
!     subroutine tag_count - counts glaciers and tops of glaciers, and
! tags them all with markers for which belong together. Also finds
! termini. Assumes there is at least 1 non-glacierized HRU between segments
! that end at top of glacier. Will count tops initially; and count and tag
! glaciers and tag tops each round. Ntp does not change, even if glacier
! disapears.
!***********************************************************************
      SUBROUTINE tag_count(do_init, hru_flowline, toflowline, glacier_frac_use)
      USE PRMS_GLACR, ONLY: Ntp, Ngl, Glacr_tag, Term, Top, Top_tag, Tohru, Nhru, GLACIER
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_type
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: do_init
      INTEGER, INTENT(OUT) :: hru_flowline(Nhru), toflowline(Nhru)
! Local Variables
      INTEGER :: i, j, ii, jj, count, label, curr, next, gltag0
      INTEGER :: numup(Nhru)
      REAL :: glacier_frac_use(Nhru)
!***********************************************************************
      Ngl = 0
      Term = 0
      IF (do_init==1) Ntp = 0
      toflowline = 0
      DO ii = 1, Active_hrus
        j = Hru_route_order(ii)
        IF ( Hru_type(j)==GLACIER ) THEN
          numup(j) = -1
          count = 0
          IF ( glacier_frac_use(j)>0.0 ) THEN ! in glacier HRU that is glacierized
            IF (glacier_frac_use(j)==1.0)  THEN
              IF ( Tohru(j)==0 ) THEN !nowhere else to extend to
                Ngl = Ngl + 1
                Term(Ngl) = j !only one terminus per glacier
              ELSEIF ( glacier_frac_use(Tohru(j))==0.0 ) THEN ! next HRU empty of glacier but could extend
                Ngl = Ngl + 1
                Term(Ngl) = j !only one terminus per glacier
              ENDIF
            ELSE ! glacier_frac_use(j)<1.0 but not 0, so at end, but next HRU could be another glacier
              Ngl = Ngl + 1
              Term(Ngl) = j !only one terminus per glacier
            ENDIF
            DO jj = 1, Active_hrus
              i = Hru_route_order(jj)
              IF ( Hru_type(i)==GLACIER ) THEN
                IF ( glacier_frac_use(i)==1.0 .AND. Tohru(i)==j ) THEN
                  !segment with a glacier does feed to j
                  count = count + 1
                ENDIF
              ENDIF
            ENDDO
            IF ( count==0 ) THEN !j is a top
              IF (do_init==1) THEN
                Ntp = Ntp + 1
                Top(Ntp) = j
              ENDIF
            ELSE !number that feed into j
              numup(j) = count !j not a top
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
      ! zero it all as a switch
      hru_flowline = 0
      Glacr_tag = 0
      Top_tag = 0
      label = 1
      ii = 1

100   DO WHILE ( ii<=Ntp )
        hru_flowline(Top(ii)) = label
        Glacr_tag(Top(ii)) = label
        Top_tag(Top(ii)) = label
        curr = Top(ii)
        DO j = 1, Ngl
          IF ( curr==Term(j) ) THEN
            label = label + 1  !so each term has a label
            ii = ii + 1
            GOTO 100           !end of glacier, go to next top
          ENDIF
        ENDDO
        next = Tohru(curr)
        IF (next==0) THEN !got to end of glacier
          label = label + 1  !so each term has a label
          ii = ii + 1
          GOTO 100           !end of glacier, go to next top
        ENDIF
        DO
          IF ( Glacr_tag(next)==0 ) hru_flowline(next) = label
          IF ( Glacr_tag(next)/=0 .AND. numup(next)>0 ) THEN
            !already been assigned and at confluence
            gltag0 = Glacr_tag(next)
            DO jj = 1, Active_hrus
              i = Hru_route_order(jj)
              !give all glaciers the same tag
              IF ( Hru_type(i)==GLACIER .AND. Glacr_tag(i)==gltag0 ) Glacr_tag(i) = label
             ENDDO
            Top_tag(next) = -1
            toflowline(label) = hru_flowline(next)
            !if single top glacier, will never have a negative part
          ELSE
            Top_tag(next) = label
          ENDIF
          Glacr_tag(next) = label
          curr = next
          DO j = 1, Ngl
            IF ( curr==Term(j) ) THEN
              label = label + 1 !so each term has a label
              ii = ii + 1
              GOTO 100 !end of glacier, go to next top
            ENDIF
          ENDDO
          next = Tohru(curr)
          IF (next==0) THEN !got to end of glacier that no longer exists
            label = label + 1  !so each term has a label
            ii = ii + 1
            GOTO 100           !end of glacier, go to next top
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE tag_count


!***********************************************************************
!     subroutine bottom - calculates bottom topo using Salamatin and Mazo
! equations (1985) without optimization for steady state, instead
! needs a proxy for steady state mass balance. Can do this from mass
! balance calculation with climate data first year (MbInit_flag=1)
! or use max and min balance above and below Ela, respectively and assume
! constant mass balance gradient above and below Ela; e.g. Farinotti (MbInit_flag=2)
! All mass balances are adjust to put glacier in steady state.
!
! The method of Salamatin and Mazo is from conservation of mass,
! solving k*(u-z)*F1(sig)+(u-z)**2*F2(sig)+A(x)=0
! sig=(u-z)*du/dx, u(1)=z(1) (the terminus has height 0)
! F1 =@(sig) sig.^((nn+1)/2); %Kamb law for sliding ice
! F2 =@(sig) (sig.^nn)/(nn+2); %polynomial rheological dependence for
! ice which is sliding + creeping + mass balance distribution =0
!-- out+in=bal
!#of cells=Nhrugl,#of streams=Ntp,#of cells/stream<=Ntp,
!#of glaciers<=Nhru
!***********************************************************************
      SUBROUTINE bottom(Frawt, Gln, Gt, Topn, Gl_top, Av_elev, Ela_elevt, Flow_slope, &
     &                   Toflowline, Slope, Ode_area, Ode_vol,Botwrite)
      USE PRMS_GLACR, ONLY: Ntp, Nhrugl, All_unit, Output_unit, Density, Gravity, Aflow, &
     &   Order_flowline, Max_gldepth, Keep_gl, Ikeep_gl, Hru_length, Hru_slope_ts
      USE PRMS_FLOWVARS, ONLY: Glacier_frac
      IMPLICIT NONE
      INTEGER, PARAMETER :: N = 101
! Functions
      EXTERNAL :: cumtrapzv, spline, splint, solve_poly
      INTRINSIC :: MIN, MAX, ISNAN, SQRT
! Arguments
      INTEGER, INTENT(IN) :: Gln, Topn, Gl_top(Ntp), Gt(Ntp), Toflowline(Ntp), Botwrite
      REAL, INTENT(IN) :: Ela_elevt(Ntp), Frawt(Nhrugl)
      REAL, INTENT(OUT) :: Av_elev(Nhrugl), Flow_slope(Nhrugl), Slope(Nhrugl)
      REAL, INTENT(OUT) :: Ode_area(Nhrugl), Ode_vol(Nhrugl)
! Local Variables
      INTEGER :: i, j, ii, jj, flag, next_do(Ntp), order(Ntp, 2), dothis0, dothis, done
      INTEGER :: thestr, inda, indb, nextstr, ibot, itop, dont, len_str, spg(Nhrugl)
      INTEGER :: str_done(Ntp), len_str_true, cell_id(Nhrugl), str_id(Nhrugl)
      INTEGER :: is(Nhrugl), ie(Nhrugl), iem(Nhrugl)
      REAL :: k, h, hvec(Nhrugl+1), balt(Nhrugl), uraw(Nhrugl), slp, setmax, ela_elev
      REAL :: xraw(Nhrugl), hraw(Nhrugl), fraw(Nhrugl), araw0(Nhrugl), araw(Nhrugl)
      REAL :: sa(Nhrugl+1), arawe0(Nhrugl+1), hvecn(N-1), zv(N), hv(N), xv(N), uv(N)
      REAL :: upv(N), balraw(Nhrugl), add, balv(N), s(N),bal, amt, amt2, bot, top, zpv(N)
      REAL :: hvecn2(N-1), intf3(N), yd, dv(N), dv0(N),dpv(N), dxh(N), areag, volg
      REAL :: alld(N, Ntp), area(Gln), vol(Gln), y2d(N), zraw(Nhrugl+2), plinetop
      REAL :: reparea, repvol, intf5(N), sina, stress, balx, pline(Nhrugl), draw(Nhrugl)
      REAL :: rline(Nhrugl), y2a(Nhrugl), y2b(Nhrugl), y2c(Nhrugl), y2dd0(Nhrugl)
      REAL :: allx(N, Ntp), zraw_av(Nhrugl+2), y2dd(Nhrugl), xrawe(2), xrawm(Nhrugl+2)
      REAL :: balrawe(2), urawe(2), hrawe(2), arawe(2), frawe(2), junk(2), nn, divu, kappa
      REAL :: sh(Nhrugl+2), shf(Nhrugl+2), hf(Nhrugl+2), hraw2(Nhrugl+2), fd(Nhrugl)
      REAL :: kk(Ntp), ll(Ntp), xrawterm(Ntp), urawterm(Ntp), urawtop(Ntp), minterm
      REAL :: allu(N, Ntp), last_frac, check_amt, flowc, dv1k
!***********************************************************************
      divu = 1.E3
      kappa = 0.04 !from Mazo 1995, value after all scaling
      junk(1) = 1.0E36
      junk(2) = 1.0E36
      sa = 1.0E15
      sh = 1.0E15
      shf = 1.0E15
      s = 1.0E15
      area = 0.0  !area of glacier
      vol = 0.0   !vol of glacier
      spg = 0     !streams per glacier
      balt = 0.0  !initialize
      hvec = 0.0
      DO j = 1, Ntp
        next_do(j) = Ntp + 10000
      ENDDO
      DO i = 1, Nhrugl
        cell_id(i) = Ikeep_gl(i,1)
        str_id(i) = Ikeep_gl(i,2)
      ENDDO
      nn = 3.0 ! coefficient of creep, Farinotti used 3, Mazo used 2.2, smaller value makes glacier thicker
      minterm = 1.0E15
      DO j = 1, Ntp
        is(j) = Ikeep_gl(j,3)
        ie(j) = Ikeep_gl(j,4)
        urawtop(j) = Keep_gl(j,5) ! in km, -terminus
        urawterm(j) = Keep_gl(j,6) ! in km,
        IF (urawterm(j)< minterm) minterm = urawterm(j) ! lowest terminus in system
        xrawterm(j) = Keep_gl(j,7) ! in km
        ll(j) = xrawterm(j)*divu  !m
        flowc = 1/(2*Aflow*3.154E7) ! make Pa^nn yr
        !PRINT *,flowc, "1/2A Pa^3/yr", Aflow*3.154E7, "A Pa^-3/yr" !might want to see these
        kk(j) = (( flowc*(ll(j)**(nn+1)) ) /( ((Density*Gravity)**nn) ))**(1.0/(2.0*nn+2)) ! scaling
        len_str = ie(j) - is(j) + 1
        len_str_true = len_str
        DO i = 1, len_str
          IF (Glacier_frac(cell_id(is(j)+i-1)) > 0.0 ) len_str_true = i
          !will write over till end where glacier ends, id = cell_id(is(thestr)+len_str_true-1,1)
        ENDDO
        iem(j) = is(j)+len_str_true -1 !need for all strings before start
      ENDDO!
      IF (Botwrite==1) OPEN ( All_unit, FILE='all.dat' ) !save solutions together
      done = 0
      nextstr = 1
      inda = 1
!start computations with top streams
      dothis = 0
50    DO ii = 1, Topn
        thestr = Gl_top(ii)
        IF ( dothis/=0 ) thestr = dothis !done top streams, not in loop
        k = kappa
        ela_elev=(Ela_elevt(thestr)-urawterm(thestr)*divu)/kk(thestr) !want this to be from from steady state method of ELA (from compute_ela_aar)
        len_str = ie(thestr) - is(thestr) + 1
        len_str_true = iem(thestr) - is(thestr) + 1
        DO i = 1, len_str_true
          uraw(i) = (Keep_gl(is(thestr)+i-1,1)*divu)/kk(thestr) !in - terminus height and fraction of kk(thestr), so O(1) and for scaling
          xraw(i) = Keep_gl(is(thestr)+i-1,2) !in fraction of length ll(thestr), so O(1)
          hraw(i) = Keep_gl(is(thestr)+i-1,3) !in km, so O(1)
          balraw(i) = balt(is(thestr)+i-1)/hraw(i) !need to spread it out
          hraw2(i+1) = hraw(i)
          fraw(i) = Frawt(is(thestr)+i-1) !in m, so O(1)
          IF ( i>1 ) hvec(i) = xraw(i) - xraw(i-1)
        ENDDO
        hvec(1) = xraw(1)
        xrawe(1) = 0.0
        urawe(1) = (urawtop(thestr)*divu)/kk(thestr)
        last_frac = Glacier_frac(cell_id(iem(thestr)))
        xrawe(2) = xraw(len_str_true) + Hru_length(cell_id(iem(thestr)))*last_frac*0.5/xrawterm(thestr) !scale
        hvec(len_str_true+1) = xrawe(2)-xraw(len_str_true)
        urawe(2) = (uraw(len_str_true)*kk(thestr)/divu - Hru_slope_ts(cell_id(iem(thestr)))* &
     &                 Hru_length(cell_id(iem(thestr)))*last_frac*0.5)*divu/kk(thestr) !because slope is set positive but is negative
        ! extrapolate these
        IF (len_str_true>1) THEN
          hrawe(1) = hraw(1) - (hraw(2)-hraw(1))/hvec(2)*hvec(1)
          IF (hrawe(1)<=0.0) hrawe(1) = hraw(1)/10.0 !arbitrary
          hrawe(2) = hraw(len_str_true) + (hraw(len_str_true)-hraw(len_str_true-1)) &
     &         /hvec(len_str_true)*hvec(len_str_true+1)
          IF (hrawe(2)<=0.0) hrawe(2) = hraw(2)/10.0 !arbitrary
        ELSE !make square glacier
          hrawe(1) = hraw(1)
          hrawe(2) = hraw(1)
        ENDIF
        hraw2(1) = hrawe(1)
        hraw2(len_str_true+2) = hrawe(2)
! need frawe and fraw depending on Mbinit_flag
        CALL cumtrapzv(hraw2(1:len_str_true+2), len_str_true+2, hvec(1:len_str_true+1), sh) !sh is cumulative area
        CALL getf_fgrad(ela_elev,len_str_true,sh,urawe,uraw,xraw,hvec,frawe,fraw,fd)
!
        !need to adjust fraw to be steady state, so shift up so int((f(z)+a)*area(z))dz= 0
        ! then add(x)=-int(f(z(x))*area(z(x)))dz/int(area(z(x)))dz
        ! if add is constant in z, then constant in x because no matter what z is, same add
        !If only one Hru on glacier, fraw will all be 0 to be in steady state
        ! add would be zero if did the full integral with Mbinit_flag==2, but wrong at the moment
        hf(1) = hrawe(1)*frawe(1)
        DO i = 2, len_str_true+1
          hf(i) = hraw(i-1)*fraw(i-1)
        ENDDO
        hf(len_str_true+2) = hrawe(2)*frawe(2)
        CALL cumtrapzv(hf(1:len_str_true+2), len_str_true+2, hvec(1:len_str_true+1), shf)
        add = -shf(len_str_true+2)/sh(len_str_true+2) !should be close to 0 if picked Mbinit_flag=2 because adjusted
        DO i = 1, len_str_true
          araw0(i) = hraw(i)*(fraw(i)+add)
          arawe0(i+1) = araw0(i)
        ENDDO
        arawe(2) = 0.0 !boundary condition
        arawe0(1) = hrawe(1)*(frawe(1) + add)
        CALL cumtrapzv(arawe0(1:len_str_true+1), len_str_true+1, hvec(1:len_str_true), sa)
        arawe(1) = (1.0/hrawe(1))*sa(1)
        IF (Botwrite==1) WRITE ( Output_unit,* ) xrawe(1),urawe(1),arawe(1),hrawe(1),frawe(1)+add,frawe(1), len_str_true
        DO i = 1, len_str_true
          araw(i) = (1.0/hraw(i))*sa(i+1)
          IF ( araw(i)<0.0 ) THEN
            ! If A is negative, then won't be able to converge.
            ! Shouldn't happen if MB is steady and glacier gets more negative MB as decrease elevation
            !PRINT *, 'Zeroed integral A (=',araw(i),') at x =',xraw(i)*ll(thestr)/1.0E3,'scaled x =',xraw(i)
            IF (i>1) araw(i) = araw(i-1)!1.0E-5
            IF (i==1) araw(i) = 1.0E-5
          ENDIF
          IF (Botwrite==1) WRITE ( Output_unit,* ) xraw(i),uraw(i),araw(i),hraw(i),fraw(i)+add,fraw(i)
        ENDDO
        IF (Botwrite==1) WRITE ( Output_unit, * ) xrawe(2),urawe(2),arawe(2),hrawe(2),frawe(2)+add,frawe(2),len_str_true, &
     &                           xrawterm(thestr), 0, kk(thestr)/divu, urawterm(thestr)
        CALL spline(xraw(1:len_str_true), uraw(1:len_str_true), xrawe, urawe, len_str_true, 1.0E31, 1.0E31, y2a)
        CALL spline(xraw(1:len_str_true), hraw(1:len_str_true), xrawe, hrawe, len_str_true, 1.0E31, 1.0E31, y2b)
        IF ( dothis/=0 ) THEN
          balrawe(1) = 0.0
          balrawe(2) = 0.0
          CALL spline(xraw(1:len_str_true), balraw(1:len_str_true), xrawe, balrawe, len_str_true, 1.0E31, 1.0E31, y2c)
        ENDIF
        h = xrawe(2)/(N-1)
        DO i = 1, N
          xv(i) = (i-1)*h
          CALL splint(xraw(1:len_str_true), uraw(1:len_str_true), xrawe, urawe, y2a(1:len_str_true), &
     &                len_str_true, xv(i), uv(i), yd)
          CALL splint(xraw(1:len_str_true), hraw(1:len_str_true), xrawe, hrawe, y2b(1:len_str_true), &
     &                len_str_true, xv(i), hv(i), yd)
          balv(i) = 0.0 !once isn't a top stream will be a value
          IF ( i<N ) hvecn(i) = h
          IF ( dothis/=0 ) CALL splint(xraw(1:len_str_true), balraw(1:len_str_true), xrawe, balrawe, &
     &                                 y2c(1:len_str_true), len_str_true, xv(i), balv(i), yd)
        ENDDO
        upv(1)=(uv(2)-uv(1))/(xv(2)-xv(1))
        DO i = 2, N-1
          upv(i)= (uv(i+1)-uv(i-1))/(xv(i+1)-xv(i-1))
        ENDDO
        upv(N)= (uv(N)-uv(N-1))/(xv(N)-xv(N-1))
        setmax = Max_gldepth ! without tight bound things go wrong
        CALL solve_poly(dv0, dpv, flag, xraw(1:len_str_true), araw(1:len_str_true), arawe, xrawe, xv, upv, balv, &
     &                  k, nn, len_str_true, setmax)
        dv(N)=dv0(N)
        dont = 0
        DO i = 1, N ! lots of stability issues here, make better somehow?
          IF ( i<N ) dv(i) = ( dv0(i)+dv0(i+1) )/2.0 !because this polySolve isn't smooth, unstable
          IF ( i>1 .AND. dv(i)<0.0 ) dv(i)=dv(i-1)
          IF ( dv(1)<0.0 ) dv(1)=0.0
          zv(i) = uv(i) - dv(i)
          IF (zv(i) < -2*uv(1)) dont = 1 ! below terminus elevation by twice height elevation,something wrong
          IF (len_str_true==1) dont = 1 !SHOULD MAKE THIS IN WEASEL SO CAN'T HAPPEN only one HRU glacierized
        ENDDO
        IF (dont==1) THEN ! use Li 2012 method with an estimated stress, since having issues with ODE solve
          stress = 78800 !average of Li 2012 Pascal values
          DO i = 1, N
            slp = upv(i)*kk(thestr)/ll(thestr)
            sina = -1.0*slp/SQRT( slp**2.0 + 1.0 )
            dv(i) = (stress/(Density*Gravity*sina))/(1 - (stress/(Density*Gravity*sina))/(0.45*hv(i)*divu))
            dv(i) = dv(i)/kk(thestr)
            zv(i) = uv(i) - dv(i)
          ENDDO
        ENDIF
        DO j = 1, Ntp
          IF ( j/=thestr ) CYCLE
          nextstr = Toflowline(j)  !can only be one
          reparea = 0.0
          repvol = 0.0
          IF ( nextstr==0 ) THEN
            indb = N + 1 !make bigger than any string
            GOTO 200     !skip bal stuff
          ENDIF
          amt = 1.0E10
          DO jj = is(nextstr), iem(nextstr)  !Keep_gl(,1) is uraw, looking for closest HRU in elevation
            check_amt = ABS(Keep_gl(jj,1)*divu+urawterm(nextstr)*divu - &
     &                       Keep_gl(iem(thestr),1)*divu-urawterm(thestr)*divu)
            IF ( amt>=check_amt ) THEN
              amt = check_amt
              inda = jj
            ENDIF
          ENDDO
        ENDDO
        balx = Keep_gl(iem(thestr),2)*ll(thestr) - Keep_gl(inda,2)*0.5*ll(nextstr) !Keep_gl(,2) is xraw in km
        amt = 1.0E10
        indb = N        !default
        DO j = 1, N - 1 !can't be terminus, looking for closest HRU in x distance
          check_amt = ABS(xv(j)*ll(thestr)-balx)
            IF ( amt>=check_amt ) THEN
              amt = check_amt
            indb = j
          ENDIF
        ENDDO
        DO j = 1, N - indb + 1
          intf3(j) = (uv(j+indb-1)-zv(j+indb-1))*hv(j+indb-1)!cross section area
          intf5(j) = hv(j+indb-1)  !width
          IF ( j<N-indb+1 ) hvecn2(j) = h
        ENDDO
        CALL cumtrapzv(intf3(1:N-indb+1), N-indb+1, hvecn2(1:N-indb), s)  !volume
        repvol = s(N-indb+1)*ll(thestr)*kk(thestr)/divu       !in km3 kinda
        bal = repvol*divu/ll(nextstr)/kk(thestr)*0.5
!unscale x and rescale to new x, also divide by 2 since bottom
! approximately triangle
        CALL cumtrapzv(intf5(1:N-indb+1), N-indb+1, hvecn2(1:N-indb), s)  !area
        reparea = s(N-indb+1)*ll(thestr)/divu      !in km2
!spread it hv(indb) onto surrounding inda and save
        amt = 1.0E10
        amt2 = 1.0E10
        bot = Keep_gl(inda,2)*ll(nextstr)/divu - hv(indb)*0.5 !Keep_gl(,2) is xraw, put in km
        top = Keep_gl(inda,2)*ll(nextstr)/divu + hv(indb)*0.5
        ibot = 1 !rsr, make sure these have values
        itop = 0 !rsr, make sure these have values
        DO j = is(nextstr), iem(nextstr) !Keep_gl(,2) is xraw, look for closest x to bot
          check_amt = ABS(Keep_gl(j,2)*ll(nextstr)/divu-bot)
          IF ( amt>=check_amt ) THEN
            amt = check_amt
            ibot = j
          ENDIF !Keep_gl(,2) is xraw, look for closest x to top
          check_amt = ABS(Keep_gl(j,2)*ll(nextstr)-top)
          IF ( amt2>=check_amt ) THEN
            amt2 = check_amt
            itop = j
          ENDIF
        ENDDO
        DO j = ibot, itop
          balt(j) = balt(j) + bal/(itop-ibot+1)
        ENDDO
!
! save solutions together
 200    done = done + 1
        DO i = 1, N
          zv(i) = uv(i) - dv(i)
          IF (Botwrite==1) WRITE ( All_unit, '(4F13.5)' ) xv(i), uv(i), zv(i), hv(i)
          dxh(i) = dv(i)*hv(i)    !in km2
          alld(i, thestr) = dv(i)
          allx(i, thestr) = xv(i)
          allu(i, thestr) = uv(i)
        ENDDO
        CALL cumtrapzv(hv, N, hvecn, s)          !get area, to put in km2
        areag = s(N)*ll(thestr)/divu - reparea  !cut part that will do again
        area(Gt(thestr)) = area(Gt(thestr)) + areag
        CALL cumtrapzv(dxh, N, hvecn, s)         !get vol, to put in km3
        volg = s(N)*ll(thestr)*kk(thestr)/divu - repvol    !cut part that will do again
        vol(Gt(thestr)) = vol(Gt(thestr)) + volg
! N points, first line string, glacier, area
        IF (Botwrite==1) WRITE (All_unit, '(2I5, 2F13.5)') thestr, Gt(thestr), areag, volg
        spg(Gt(thestr)) = spg(Gt(thestr)) + 1
        order(done, 1) = thestr
        order(done, 2) = Gt(thestr)
        str_done(done) = thestr
        next_do(done) = nextstr
!
        IF ( dothis/=0 ) EXIT                  !done top streams, not in loop
!
      ENDDO !end top streams if in loop
!
!if have done all that feed into nextstr then can do nextstr
! this following part is spagetti code. all it is doing is finding the
! streams the top streams feed into, and running them through the
! inverse code next after all feeders are done
!all glaciers have a top stream, so since have all those done, will
!  account for all glaciers
      Q1:DO i = 1, done
        IF ( next_do(i)/=0 ) THEN             !have more to do after top strings
          DO ii = 1, done
            IF ( str_done(ii)==next_do(i) ) THEN  !did it
              IF ( i<=done-1 ) GOTO 300
              IF ( i==done ) GOTO 400  !checked all
            ENDIF
          ENDDO   !got here, found a not done string, or all done
          dothis0 = next_do(i)
          DO
            Q2:DO j = 1, Ntp           !check that did feeders
              IF ( Toflowline(j)==dothis0 ) THEN
                Q3:DO jj = 1, done
                  IF ( j==str_done(jj) ) THEN  !did it
                    IF ( j<=Ntp-1 ) GOTO 210
                    IF ( j==Ntp ) EXIT Q3   !checked all
                  ELSE    !didn't do it
                    dothis0 = j
                    GOTO 220    !so check feeders, back at j=1
                  ENDIF
                ENDDO Q3
              ENDIF     !got here, found next string to do
              dothis = dothis0
              GOTO 50   !the start of the inverse problem for depths
 210        ENDDO Q2
            EXIT
 220      ENDDO
        ENDIF
 300  ENDDO Q1
 400  CLOSE ( All_unit )
!Try to calulate steady surface cross-glacier-profile at every x point
! So want depth function cross-profile (changes with every x), call this d.
! d along the flow line is u-z, =dflow
!Suppose H(x,d)=P(x)*d^R(x) gives a profile for unknown P and R
! d along the flow line is u(x)-z(x), =dflow a function of x
! then H(x)=P(x)*dflow^R(x)
! We know H(x)-- solved for it, rearrange: P(x)=H(x)/dflow^R(x)
! Cross section area S(x) by integrating H(x,d) over depth from 0 to dmax=dflow
! S(x)=(1/(R(x)+1))*P(x)*dflow^(R(x)+1)
!
      DO i = 1, Ntp
        DO j = 1, Ntp
          Order_flowline(j) = order(j, 1)
          IF ( str_id(is(i))==Order_flowline(j) ) THEN  !do this in order
            thestr = Order_flowline(j)
            len_str = ie(thestr) - is(thestr) + 1
            len_str_true = iem(thestr) - is(thestr) + 1
            areag = area(order(j, 2))   !per glacier all in km2
            volg = vol(order(j, 2))     !per glacier all in km3
            Ode_area(iem(thestr)) = areag !store in terminus HRU
            Ode_vol(iem(thestr)) = volg !store in terminus HRU
! this is generous volume because not off average depth, off flowline depth
            DO ii = 1, N
              dv(ii) = alld(ii, thestr) !in per kk
              xv(ii) = allx(ii, thestr) !in km/xterm
              uv(ii) = allu(ii, thestr) !in per kk - uterm*divu/kk
            ENDDO
! uv and zpv is currently not used for anything, but might want to calculate more accurate basal slope
            zpv(1)=(zv(2)-zv(1))/(xv(2)-xv(1))
            DO ii = 2, N-1
              zpv(ii)= (zv(ii+1)-zv(ii-1))/(xv(ii+1)-xv(ii-1)) ! units in per kk
            ENDDO
            zpv(N)= (zv(N)-zv(N-1))/(xv(N)-xv(N-1))
            DO ii = 1, len_str
              uraw(ii) = Keep_gl(is(thestr)+ii-1,1)! in km get whole thing
              xraw(ii) = Keep_gl(is(thestr)+ii-1,2)! scaled get whole thing
              xrawm(ii) = ll(thestr)*xraw(ii) !unscaled in m
              hraw(ii) = Keep_gl(is(thestr)+ii-1,3)  !in km
            ENDDO
            xrawm(len_str+2) = 0.0
            xrawm(len_str+1) = ll(thestr)
            CALL spline(xv, dv, junk, junk, N, 1.0E31, 1.0E31, y2d)
            draw = 0.0
            DO ii = 1, len_str_true
! make R constant =r0, have to because only have one equation for unknown
!  rline denominator from 1.5 at top of each glacier to 2.1 at lowest extent glacier
              rline(ii) = 1.0/(2.1+.6*(uraw(ii)+urawterm(thestr)-minterm)/ &
     &                       (minterm-urawtop(thestr)-urawterm(thestr)))
              CALL splint(xv, dv, junk, junk, y2d, N, xraw(ii), draw(ii), yd)
              IF ( draw(ii)<0.0 ) draw(ii)=0.0 !splint might make some slightly negative
              draw(ii) = draw(ii)*kk(thestr)/divu !make in km
              IF ( draw(ii)==0.0 ) THEN
                pline(ii) = 0.0
              ELSE
                pline(ii) = hraw(ii)/(draw(ii)**rline(ii))
              ENDIF
!
! units of integer,integer,km,km,km, and pline,rline equation will give H in km
! a bottom width for each depth-- center it on flowline
!
! Basal hru_elev
!Calulate average unscaled elevation of bottom and store (zraw_av)
!uraw - average depth= uraw - cross section area divided by H max
              zraw_av(ii) =( uraw(ii) - ( pline(ii)/(rline(ii)+1.0)*draw(ii)**(rline(ii)+1.0) )&
     &          /hraw(ii) + urawterm(thestr))*divu !in m
              zraw(ii) = ( uraw(ii) - draw(ii) + urawterm(thestr))*divu !in m
              Av_elev(is(thestr)+ii-1) = zraw_av(ii)
            ENDDO
            IF ( len_str/=len_str_true) THEN
              DO ii = len_str_true+1, len_str !set all to hru_elev if not glacierized at start
                zraw_av(ii) = (uraw(ii) + urawterm(thestr))*divu !in m
                zraw(ii) = (uraw(ii) + urawterm(thestr))*divu !in m
                Av_elev(is(thestr)+ii-1) = zraw_av(ii)
              ENDDO
            ENDIF
! Basal hru_slope
! Calulate a version of average slope with zraw_av, going down the flowline, scale it and store
! make endpoints at u elevation-- make compatable with way calculate on unglacierized HRUs
! Also calculate centerline slope with zraw. Could more accurately do this by integrating zv over xv
!  but then would change when got off glacier
! Top and bottom of glacier affected by steep negative slope and steep positive slope respectively, may want to exclude
            dv1k = dv(1)*kk(thestr)/divu !in km
            plinetop = hraw(1)/(dv1k**(1.0/2.1))
            IF ( dv1k==0.0 ) plinetop = 0.0
            zraw_av(len_str+2) = ((urawtop(thestr)- ( plinetop/(rline(1)+1.0)*dv1k**(rline(1)+1.0) )/hraw(1)) &
     &                 + urawterm(thestr))*divu
            zraw_av(len_str+1) = urawterm(thestr)*divu
            zraw(len_str+2) = ((urawtop(thestr)-dv1k)+urawterm(thestr))*divu
            zraw(len_str+1) = urawterm(thestr)*divu
            y2dd(1)= (zraw_av(2)-zraw_av(len_str+2))/(xrawm(2)-xrawm(len_str+2))
            y2dd0(1)= (zraw(2)-zraw(len_str+2))/(xrawm(2)-xrawm(len_str+2))
            Flow_slope(is(thestr)) = ABS(y2dd(1))
            Slope(is(thestr)) = y2dd0(1)
            DO ii = 2, len_str
               y2dd(ii)= (zraw_av(ii+1)-zraw_av(ii-1))/(xrawm(ii+1)-xrawm(ii-1))
               y2dd0(ii)= (zraw(ii+1)-zraw(ii-1))/(xrawm(ii+1)-xrawm(ii-1))
               Flow_slope(is(thestr)+ii-1) = ABS(y2dd(ii)) !Basal_slope
               Slope(is(thestr)+ii-1) = y2dd0(ii)
            ENDDO
! Basal hru_aspect: assume aspect the same on surface and base of glacier
          ENDIF
        ENDDO
      ENDDO
!
      END SUBROUTINE bottom

!***********************************************************************
!     subroutine getf_fgrad - extrapolates f and computes fgrad
!***********************************************************************
      SUBROUTINE getf_fgrad(ela_elev,len_str,sh,urawe,uraw,xraw,hvec,frawe,fraw,fd)
      USE PRMS_GLACR, ONLY: Mbinit_flag, Nhrugl
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: len_str
      REAL, INTENT(IN) :: urawe(2), uraw(Nhrugl), xraw(Nhrugl), hvec(Nhrugl), ela_elev
      REAL, INTENT(IN) :: sh(Nhrugl+2)
      REAL, INTENT(INOUT) :: fraw(Nhrugl)
      REAL, INTENT(OUT) :: fd(Nhrugl), frawe(2)
! Local Variables
      INTEGER :: i, ela_i
      REAL :: ela_x, high, low, add0, cons, aarA, sl_acc, sl_abl, frawm(Nhrugl+2)
      REAL :: urawm(Nhrugl+2)
!***********************************************************************
      ela_i = 0
      ela_x = 0.0
      fd = 0.0
      IF (Mbinit_flag==1) THEN !use climate data for mass balance
        IF (len_str>1) THEN
          frawe(1) = fraw(1)- (fraw(2)-fraw(1))/hvec(2)*hvec(1)
          frawe(2) = fraw(len_str) + (fraw(len_str)-fraw(len_str-1))/hvec(len_str)*(1-xraw(len_str))
        ELSE
          frawe(1) = fraw(1)
          frawe(2) = fraw(1)
        ENDIF
      ELSEIF (Mbinit_flag==2) THEN !Farinotti method, using ela_x
        DO i = 1, len_str
          IF (ela_elev-uraw(i)<0.0001) THEN !has rounding errors
            ela_x=xraw(i)
            ela_i=i
          ENDIF
        ENDDO
        high =-1.0E7
        low = 1.0E7
        DO i = 1, len_str
          IF (i<=ela_i) high = MAX(high,fraw(i))
          IF (i>=ela_i) low =  MIN(low,fraw(i))
        ENDDO
        ! say representative of steady state slopes
        ! solve for accumulation and ablation slopes (sl_acc and sl_abl) with integrated f above and below ela_x same
        ! could use 1.8*sl_acc= sl_abl from Furbish 1984
        ! is not quite right since width should be function of altitude
        aarA = sh(ela_i+1)/sh(len_str+2) !aar with this ela_x
        cons = (ela_elev)*(1.0 - aarA)/(urawe(1) - ela_elev)/aarA
        add0 = (-cons*low - high)/(1.0 + cons) ! f value at ela, make it the 0 point for steady state
        frawe(1) = high + add0
        frawe(2) = low  + add0
        sl_acc =  frawe(1)/(urawe(1) - ela_elev)
        sl_abl = -frawe(2)/ela_elev
        DO i = 1, ela_i
          fraw(i) = sl_acc*(uraw(i)- ela_elev)
        ENDDO
      ! fraw(ela_i) is done twice, both times should be = 0
        DO i = ela_i, len_str
          fraw(i) = frawe(2) + sl_abl*uraw(i)
        ENDDO
      ENDIF
      DO i = 1, len_str
        frawm(i) = fraw(i)
        urawm(i) = uraw(i)
      ENDDO
      frawm(len_str+2) = frawe(1)
      frawm(len_str+1) = frawe(2)
      urawm(len_str+2) = urawe(1)
      urawm(len_str+1) = urawe(2)
      ! both scaled by divu (except in bottom routine, but don't use fd there) so slope unscaled
      fd(1)= ( (frawm(2)-frawm(1))/(urawm(2)-urawm(1)) + &
     &          (frawm(1)-frawm(len_str+2))/(urawm(1)-urawm(len_str+2)) )/2.0
      DO i = 2, len_str
        fd(i) = ( (frawm(i+1)-frawm(i))/(urawm(i+1)-urawm(i)) + &
     &         (frawm(i)-frawm(i-1))/(urawm(i)-urawm(i-1)) )/2.0
      ENDDO
      END SUBROUTINE getf_fgrad

!***********************************************************************
!     subroutine yearly_ca_coef - calculates average basal slope at
! center (deepest part of glacier) and average mass balance gradient at
! center for the year (and yearly extent)
! NOTE: The calculations only use the HRUs that are initially glacierized,
! so if glacier grows, will be neglecting extension. MIGHT WANT TO CHANGE
! This wouldn't make a huge difference since extension can't be that long.
! Don't have basal slope calculated there so would need to use something else.
!***********************************************************************
      SUBROUTINE yearly_ca_coef(Frawt, Ela_elevt)
      USE PRMS_GLACR, ONLY: Ntp, Nhrugl, Ngl, Order_flowline, Keep_gl, Ikeep_gl, &
     &    Hru_length, Av_basal_slope, Av_fgrad, Glacr_tag, Term, Glacr_slope_init, &
     &    Hru_length, Nhru, GLACIER
      USE PRMS_BASIN, ONLY: Active_hrus, Hru_route_order, Hru_type
      USE PRMS_FLOWVARS, ONLY: Glacier_frac
      IMPLICIT NONE
! Functions
      INTRINSIC :: ISNAN
      EXTERNAL :: cumtrapzv, getf_fgrad
! Arguments
      REAL, INTENT(IN) :: Frawt(Nhrugl), Ela_elevt(Ntp)
! Local Variables
      INTEGER :: i, j, ii, o, p, thestr, len_str, str_id(Nhrugl), cell_id(Nhrugl), len_str_true
      INTEGER :: is(Ntp), ie(Ntp)
      REAL ::  ela_elev, hvec(Nhrugl+1), uraw(Nhrugl), xraw(Nhrugl), hraw(Nhrugl)
      REAL :: xrawe(2), urawe(2), hrawe(2), frawe(2), fd(Nhrugl), tot_fgrad, tot_length
      REAL :: sh(Nhrugl+2), hraw2(Nhrugl+2), fgrad(Nhrugl), fgradm(Nhru), tot_slope
      REAL :: fraw(Nhrugl), slopem(Nhru), divu
!***********************************************************************
      sh = 1.0E15 !initialize
      hvec = 0.0
      divu = 1.0E3 !in meters
      slopem = 0.0
      fgrad = 0.0
      fgradm = 0.0
      DO i = 1, Nhrugl
        cell_id(i) = Ikeep_gl(i,1)
        str_id(i) = Ikeep_gl(i,2)
        slopem(cell_id(i)) = Keep_gl(i,4) !basal centerline slope
      ENDDO
      DO ii = 1, Ntp
        is(ii) = Ikeep_gl(ii,3)
        ie(ii) = Ikeep_gl(ii,4)
        DO j = 1, Ntp
          IF ( str_id(is(ii))==Order_flowline(j)) THEN  !do this in order
            thestr = Order_flowline(j)
            len_str = ie(thestr) - is(thestr) + 1
 !want this to be from from steady state method of ELA (from compute_ela_aar)
            ela_elev=(Ela_elevt(thestr)-Keep_gl(thestr,6)*divu)/divu !subtract urawterm
            len_str_true = len_str
            DO i = 1, len_str
              IF (Glacier_frac(cell_id(is(thestr)+i-1)) > 0.0 ) len_str_true = i
              !will write over till end where glacier ends, id = cell_id(is(thestr)+len_str_true-1,1)
            ENDDO
            DO i = 1, len_str_true
              uraw(i) = Keep_gl(is(thestr)+i-1,1)
              xraw(i) = Keep_gl(is(thestr)+i-1,2)
              hraw(i) = Keep_gl(is(thestr)+i-1,3)
              hraw2(i+1) = hraw(i)
              fraw(i) = Frawt(is(thestr)+i-1)
              IF ( i>1 ) hvec(i) = xraw(i) - xraw(i-1)
            ENDDO
            hvec(1) = xraw(1)
            xrawe(1) = 0.0
            urawe(1) = Keep_gl(thestr,5) !urawtop(thestr) in km
            IF ( len_str==len_str_true) THEN
              hvec(len_str+1) = 1.0-xraw(len_str)
              xrawe(2) = 1.0
              urawe(2) = 0.0
            ELSE
              xrawe(2) = xraw(len_str_true) + Hru_length(cell_id(is(thestr)+len_str_true-1))*0.5/Keep_gl(thestr,7) !scale
              hvec(len_str_true+1) = xrawe(2)-xraw(len_str_true)
              urawe(2) = uraw(len_str_true) - Glacr_slope_init(cell_id(is(thestr)+len_str_true-1))* &
     &                     Hru_length(cell_id(is(thestr)+len_str_true-1))*0.5 !because slope is set positive but is negative
            ENDIF
           ! extrapolate these
            IF (len_str_true>1) THEN
              hrawe(1) = hraw(1) - (hraw(2)-hraw(1))/hvec(2)*hvec(1)
              IF (hrawe(1)<=0.0) hrawe(1) = hraw(1)/10.0 !arbitrary
              hrawe(2) = hraw(len_str_true) + (hraw(len_str_true)-hraw(len_str_true-1)) &
     &             /hvec(len_str_true)*hvec(len_str_true+1)
              IF (hrawe(2)<=0.0) hrawe(2) = hraw(2)/10.0 !arbitrary
            ELSE !make square glacier
              hrawe(1) = hraw(1)
              hrawe(2) = hraw(1)
            ENDIF
            hraw2(1) = hrawe(1)
            hraw2(len_str_true+2) = hrawe(2)
            CALL cumtrapzv(hraw2(1:len_str_true+2), len_str_true+2, hvec(1:len_str_true+1), sh) !sh is cumulative area
            CALL getf_fgrad(ela_elev,len_str_true,sh,urawe,uraw,xraw,hvec,frawe,fraw,fd)
            DO i = 1, len_str_true
              fgrad(is(thestr)+i-1) = fd(i)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO i = 1, Nhrugl
        fgradm(cell_id(i)) = fgrad(i)
      ENDDO
      DO o = 1, Ngl
        p = Glacr_tag(Term(o)) !index by Glacr_tag
        tot_slope = 0.0
        tot_fgrad = 0.0
        tot_length = 0.0
        DO i = 1, Active_hrus
          j = Hru_route_order(i)
          IF ( Hru_type(i)==GLACIER .AND. Glacr_tag(j)==p ) THEN
            tot_slope = tot_slope + slopem(j)*Glacier_frac(j)*Hru_length(j)
            tot_fgrad = tot_fgrad + fgradm(j)*Glacier_frac(j)*Hru_length(j)
            tot_length = tot_length + Glacier_frac(j)*Hru_length(j)
          ENDIF
        ENDDO
        Av_basal_slope(p) = -tot_slope/tot_length !needs to be the positive version
        Av_fgrad(p) = tot_fgrad/tot_length
        IF (Av_fgrad(p) <= 0.0) Av_fgrad(p) = 0.006 !from Luthi 09
      ENDDO
      END SUBROUTINE yearly_ca_coef

!***********************************************************************
!     subroutine solve_poly - finds roots of polynomial using driver of
! subroutine rtn
!***********************************************************************
      SUBROUTINE solve_poly(Dv, Dpv, Flag, Xraw, Araw, Arawe, Xrawe, Xv, Upv, Balv, K, Nn, Len_str, Setmax)
      IMPLICIT NONE
      INTEGER, PARAMETER :: N = 101
! Arguments
      INTEGER, INTENT(IN) :: Len_str
      INTEGER, INTENT(INOUT) :: Flag
      REAL, INTENT(IN) :: Xraw(Len_str), Araw(Len_str), Xv(N), Upv(N), Balv(N)
      REAL, INTENT(IN) :: Arawe(2), Xrawe(2), K, Nn, Setmax
      REAL, INTENT(OUT) :: Dv(N), Dpv(N)
! Functions
      EXTERNAL :: spline, splint, rtn
! Local Variables
      INTEGER i
      REAL h, dip1, di, upfunciph, xip1, afunciph, bfunciph, yd, xi, xiph
      REAL tol, d(N-1), dp(N-1), rtnewt, y2a(N), y2b(Len_str+2), y2c(N), junk(2)
!***********************************************************************
      h = Xrawe(2)/(N-1)
      DO i = 1, N - 1
        d(i) = 0.0
        dp(i) = 0.0
      ENDDO
      junk(1) = 1.0E36
      junk(2) = 1.0E36
      CALL spline(Xv, Upv, junk, junk, N, 1.0E31, 1.0E31, y2a)
      CALL spline(Xraw(1:Len_str), Araw(1:Len_str), Xrawe, Arawe, Len_str, 1.0E31, 1.0E31, y2b)
      CALL spline(Xv, Balv, junk, junk, N, 1.0E31, 1.0E31, y2c)
      DO i = 0, N - 2
        xip1 = Xrawe(2) - (i+1)*h
        xi = Xrawe(2) - i*h
        xiph = Xrawe(2) - (i+0.5)*h
        CALL splint(Xv, Upv, junk, junk, y2a, N, xiph, upfunciph, yd)
        CALL splint(Xraw(1:Len_str), Araw(1:Len_str), Xrawe, Arawe, y2b(1:Len_str), Len_str, xiph, afunciph, yd)
        CALL splint(Xv, Balv, junk, junk, y2c, N, xiph, bfunciph, yd)
        IF ( i>=1 ) THEN
          di = d(i) !from last iteration
        ELSE  !y at x=1 at initial point
          di = 0 !u0, =z at x=1
          Dv(N) = 0
        ENDIF
        tol = 1.0E-2
!find dip1 root of funcd between [0,Setmax],
        CALL rtn(di, upfunciph, afunciph, bfunciph, 0.0, Setmax, tol, rtnewt, Flag, K, Nn)
        dip1 = rtnewt
        d(i+1) = dip1
        dp(i+1) = (-dip1+di)/h  !dD/dx
      ENDDO
      DO i = 1, N - 1           !could have saved these in earlier do loop
        Dv(i) = d(N-1-i+1)
        Dpv(i) = dp(N-1-i+1)
      ENDDO
      Dpv(N) = Dpv(N-1)         !don't really know it here, end of interval
!
      END SUBROUTINE solve_poly

!***********************************************************************
!     subroutine funcd - defines function
!***********************************************************************
      SUBROUTINE funcd(Dip1, Di, Upfunciph, Afunciph, Bfunciph, F1, Df, K, Nn)
      IMPLICIT NONE
      INTRINSIC :: ABS
! Arguments
      REAL, INTENT(IN) :: Dip1, Di, Upfunciph, Afunciph, Bfunciph, K, Nn
      REAL, INTENT(OUT) :: F1, Df
! Local Variables
      REAL sigiph, ds
!***********************************************************************
      sigiph = 0.5*(Dip1+Di)*ABS(Upfunciph)
      ds = 0.5*ABS(Upfunciph)
      F1 = -( K*0.5*(Dip1+Di)*(sigiph**((Nn+1.0)/2.0))         &
     &       + ((0.5*(Dip1+Di))**2.0)*(sigiph**Nn)/(Nn+2.0) )  &
     &      + Afunciph + Bfunciph

      Df = -( K*0.5*(sigiph**((Nn+1)/2.0))                             &
     &       + K*0.5*(Dip1+Di)*(Nn+1)*0.5*(sigiph**((Nn-1.0)/2.0))*ds  &
     &       + .5*(Dip1+Di)*(sigiph**Nn)/(Nn+2.0)                         &
     &       + ((0.5*(Dip1+Di))**2.0)*Nn*(sigiph**(Nn-1.0))*ds/(Nn+2.0) )

      END SUBROUTINE FUNCD

!***********************************************************************
!     subroutine rtn - Using the Newton-Raphson method, and the root of
!a function known to lie in the interval [x1; x2]. The root rtnewt will
!be refined until its accuracy is known within xacc. funcd is a
!user-supplied subroutine that returns both the function value and the
!first derivative of the function at the point x.
!***********************************************************************
      SUBROUTINE rtn(Di, Upfunciph, Afunciph, Bfunciph, X1, X2, Xacc, Rtnewt, Flag, K, Nn)
      IMPLICIT NONE
      INTEGER, PARAMETER :: JMAX = 1000 !Set to maximum number of iterations.
! Arguments
      INTEGER, INTENT(OUT) :: Flag
      REAL, INTENT(IN) :: Di, Upfunciph, Afunciph, K, Nn, Bfunciph, X1, X2, Xacc
      REAL, INTENT(OUT) :: Rtnewt
! Functions
      EXTERNAL :: funcd
! Local Variables
      INTEGER j
      REAL df, dx, f
!***********************************************************************
      Flag = 0
      Rtnewt = 0.5*(X1+X2) !Initial guess.
!      Rtnewt = 1.0 !???
      DO j = 1, JMAX
        CALL funcd(Rtnewt, Di, Upfunciph, Afunciph, Bfunciph, f, df, K, Nn)
        dx = f/df
        Rtnewt = Rtnewt - dx
!            PRINT*, 'Rtnewt jumped out of brackets, Rtnewt =', Rtnewt
        IF ( Rtnewt<X1 ) Rtnewt = X1
        IF ( Rtnewt>X2 ) Rtnewt = X2
        IF ( ABS(dx)<Xacc ) RETURN  !Convergence.
      ENDDO
!      PRINT*, 'Rtnewt exceeded maximum iterations, Rtnewt = ', Rtnewt
      Flag = -2
      END SUBROUTINE rtn

!***********************************************************************
!     subroutine cumtrapzv - computes the trapezoidal rule on a vector
!of length n and spacing h, cumulative
!***********************************************************************
      SUBROUTINE cumtrapzv(Vec, N, H, S)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: H(N-1), Vec(N)
      REAL, INTENT(OUT) :: S(N)
! Local Variable
      INTEGER i
!***********************************************************************
      S(1) = 0.0
      DO i = 2, N
        S(i) = S(i-1) + 0.5*H(i-1)*(Vec(i)+Vec(i-1))
      ENDDO
      END SUBROUTINE cumtrapzv

!***********************************************************************
!     subroutine spline - given arrays x0(1:n) and y0(1:n)
!containing a tabulated function with
!ex1 the right endpoint of x (contained in array if >1.0E35),ex2 the
!left endpoint of x0 (contained in array if >1.0E35) i.e., yi = f(xi),
!with ex1 < x01 < x02 < ::: < x0N < ex2, and given values yp1 and ypn for
!the first derivative of the interpolating function at points 1 and n,
!respectively, this routine returns an array y2(1:n) of length n which
!contains the second derivatives of the interpolating function at the
!tabulated points x0i. If yp1 and/or ypn are equal to 10**30 or larger,
!the routine is signaled to set the corresponding boundary condition
!for a natural spline, with zero second derivative on that boundary.
!Parameter: NMAX is the largest anticipated value of n.
!***********************************************************************
      SUBROUTINE spline(X0, Y0, Ex, Ey, N, Yp1, Ypn, Y2)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: Yp1, Ypn, X0(N), Y0(N), Ex(2), Ey(2)
      REAL, INTENT(OUT) :: Y2(N)
! Local Variables
      INTEGER, PARAMETER :: NMAX = 1000
      INTEGER i, k
      REAL  x(N), y(N), p, qn, sig, un, u(NMAX)
!***********************************************************************
      IF ( Ex(1)>0.99E35 ) THEN
        DO i = 1, N-1
          x(i) = X0(i)
        ENDDO
      ELSE
        x(1) = Ex(1)
        DO i = 1, N-2
          x(i+1) = X0(i)
        ENDDO
      ENDIF
      IF ( Ex(2)>0.99E35 ) THEN
        x(N) = X0(N)
      ELSE
        x(N) = Ex(2)
      ENDIF
      IF ( Ey(1)>0.99E35 ) THEN
        DO i = 1, N-1
          y(i) = Y0(i)
        ENDDO
      ELSE
        y(1) = Ey(1)
        DO i = 1, N-2
          y(i+1) = Y0(i)
        ENDDO
      ENDIF
      IF ( Ey(2)>0.99E35 ) THEN
        y(N) = Y0(N)
      ELSE
        y(N) = Ey(2)
      ENDIF
!
      IF ( Yp1>0.99E30 ) THEN
!The lower boundary condition is set either to be natural
        Y2(1) = 0.0
        u(1) = 0.0
      ELSE  !or else to have a specified first derivative.
        Y2(1) = -0.5
        u(1) = (3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-Yp1)
      ENDIF
      DO i = 2, N-1
!This is the decomposition loop of the tridiagonal algorithm.
!y2 and u are used for temporary storage of the decomposed factors.
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*Y2(i-1) + 2.0
        Y2(i) = (sig-1.0)/p
        u(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
     &          /(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO
      IF ( Ypn>0.99E30 ) THEN
!The upper boundary condition is set either to be natural
        qn = 0.0
        un = 0.0
      ELSE !or else to have a specified first derivative.
        qn = 0.5
        un = (3.0/(x(N)-x(N-1)))*(Ypn-(y(N)-y(N-1))/(x(N)-x(N-1)))
      ENDIF
      Y2(N) = (un-qn*u(N-1))/(qn*Y2(N-1)+1.0)
      DO k = N-1, 1, -1
!This is the backsubstitution loop of the tridiagonal algorithm.
        Y2(k) = Y2(k)*Y2(k+1) + u(k)
      ENDDO
!
      END SUBROUTINE spline

!***********************************************************************
!     subroutine splint - given the arrays xa0(1:n) and ya0(1:n) of
!length n, which tabulate acfunction (with the xa0i 's in order),
!with ex1 the right endpoint of x (contained in array if >1.0E35)
!ex2 the left endpoint of x0 (contained in array if >1.0E35),
!and given the array y2a(1:n), which is the output from spline,
!and given a value of x, this routine returns a
!cubic-spline interpolated value y and its derivative yd.
!***********************************************************************
      SUBROUTINE splint(Xa0, Ya0, Ex, Ey, Y2a, N, X, Y, Yd)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: Xa0(N), Ya0(N), Y2a(N), X, Ex(2), Ey(2)
      REAL, INTENT(INOUT) :: Yd, Y
! Local Variables
      INTEGER i, k, khi, klo
      REAL xa(N), ya(N), a, b, h
!***********************************************************************
      IF ( Ex(1)>0.99E35 ) THEN
        DO i = 1, N - 1
          xa(i) = Xa0(i)
        ENDDO
      ELSE
        xa(1) = Ex(1)
        DO i = 1, N - 2
          xa(i+1) = Xa0(i)
        ENDDO
      ENDIF
      IF ( Ex(2)>0.99E35 ) THEN
        xa(N) = Xa0(N)
      ELSE
        xa(N) = Ex(2)
      ENDIF
      IF ( Ey(1)>0.99E35 ) THEN
        DO i = 1, N - 1
          ya(i) = Ya0(i)
        ENDDO
      ELSE
        ya(1) = Ey(1)
        DO i = 1, N - 2
          ya(i+1) = Ya0(i)
        ENDDO
      ENDIF
      IF ( Ey(2)>0.99E35 ) THEN
        ya(N) = Ya0(N)
      ELSE
        ya(N) = Ey(2)
      ENDIF

!We will find the right place in the table by means of bisection.
!This is optimal if sequential calls to this routine are at random
!values of x. If sequential calls are in order, and closely
!spaced, one would do better to store previous values of
!klo and khi and test if they remain appropriate on the next call.
!
      klo = 1
      khi = N
      DO WHILE ( khi-klo>1 )
        k = (khi+klo)/2
        IF ( xa(k)>X ) THEN
          khi = k
        ELSE
          klo = k
        ENDIF
      ENDDO
      !klo and khi now bracket the input value of x.
      h = xa(khi) - xa(klo)
      IF ( h==0.0 ) PRINT *, 'Bad xa input in splint.'
      !The xa's must be distinct.
      a = (xa(khi)-X)/h  !Cubic spline polynomial is now evaluated.
      b = (X-xa(klo))/h
      Y = a*ya(klo) + b*ya(khi) + ((a**3-a)*Y2a(klo)+(b**3-b)*Y2a(khi))*(h**2)/6.0
      Yd = -ya(klo) + ya(khi) + (-(a**2-1.0/3.0)*Y2a(klo)+(b**2-1.0/3.0)*Y2a(khi))*(h**2)/2.0
!
      END SUBROUTINE splint

!***********************************************************************
!     subroutine indexx - indexes an array Arr(1:n), i.e.,
!outputs the array Indx(1:n) such that
!Arr(Indx(j))is in ascending order for j = 1, 2, . . . ,N. The input
!quantities N and Arr are not changed.
!***********************************************************************
      SUBROUTINE indexx(N, Arr, Indx)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(OUT) :: Indx(N)
      REAL, INTENT(IN) :: Arr(N)
! Local Variables
      INTEGER, PARAMETER :: M = 7, NSTACK = 100
      INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK), ind
      REAL a
!***********************************************************************
      DO j = 1, N
        Indx(j) = j
      ENDDO
      jstack = 0
      l = 1
      ir = N
      DO
        IF ( ir-l<M ) THEN
          DO j = l+1, ir
            indxt = Indx(j)
            a = Arr(indxt)
            ind = l
            DO i = j-1, l, -1
              IF ( Arr(Indx(i))<=a ) THEN
                ind = i + 1
                EXIT
              ENDIF
              Indx(i+1) = Indx(i)
            ENDDO
            Indx(ind) = indxt
          ENDDO
          IF ( jstack==0 ) RETURN
          ir = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack - 2
        ELSE
          k = (l+ir)/2
          itemp = Indx(k)
          Indx(k) = Indx(l+1)
          Indx(l+1) = itemp
          IF ( Arr(Indx(l))>Arr(Indx(ir)) ) THEN
            itemp = Indx(l)
            Indx(l) = Indx(ir)
            Indx(ir) = itemp
          ENDIF
          IF ( Arr(Indx(l+1))>Arr(Indx(ir)) ) THEN
            itemp = Indx(l+1)
            Indx(l+1) = Indx(ir)
            Indx(ir) = itemp
          ENDIF
          IF ( Arr(Indx(l))>Arr(Indx(l+1)) ) THEN
            itemp = Indx(l)
            Indx(l) = Indx(l+1)
            Indx(l+1) = itemp
          ENDIF
          i = l + 1
          j = ir
          indxt = Indx(l+1)
          a = Arr(indxt)
          DO
            i = i + 1
            IF ( Arr(Indx(i))<a ) CYCLE
            DO
              j = j - 1
              IF ( Arr(Indx(j))>a ) CYCLE
              IF ( j<i ) THEN
                Indx(l+1) = Indx(j)
                Indx(j) = indxt
                jstack = jstack + 2
                IF ( jstack>NSTACK ) PRINT *, 'NSTACK too small in indexx.'
                IF ( ir-i+1>=j-l ) THEN
                  istack(jstack) = ir
                  istack(jstack-1) = i
                  ir = j - 1
                ELSE
                  istack(jstack) = j - 1
                  istack(jstack-1) = l
                  l = i
                ENDIF
                GOTO 100
              ELSE
                itemp = Indx(i)
                Indx(i) = Indx(j)
                Indx(j) = itemp
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
 100  ENDDO
      END SUBROUTINE indexx

!***********************************************************************
!     subroutine sort5 - sorts an array ra(1:n) into
!ascending numerical order while making the
!corresponding rearrangements of the other arrays rb to re. An
!index table is constructed via the routine indexx.
!***********************************************************************
      SUBROUTINE sort5(N, Ra, Rb, Rc, Rd, Re, Wksp, Iwksp)
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(OUT) :: Iwksp(N)
      REAL, INTENT(INOUT) :: Ra(N), Rb(N), Rc(N), Rd(N), Re(N)
      REAL, INTENT(OUT) :: Wksp(N)
! Functions
      EXTERNAL :: indexx
! Local Variables
      INTEGER j
!***********************************************************************
      CALL indexx(N, Ra, Iwksp) !Make the index table.
      DO j = 1, N               !Save the array ra.
        Wksp(j) = Ra(j)
      ENDDO
      DO j = 1, N               !Copy it back in the rearranged order.
        Ra(j) = Wksp(Iwksp(j))
      ENDDO
      DO j = 1, N               !Ditto rb.
        Wksp(j) = Rb(j)
      ENDDO
      DO j = 1, N
        Rb(j) = Wksp(Iwksp(j))
      ENDDO
      DO j = 1, N               !Ditto rc.
        Wksp(j) = Rc(j)
      ENDDO
      DO j = 1, N
        Rc(j) = Wksp(Iwksp(j))
      ENDDO
      DO j = 1, N               !Ditto rd.
        Wksp(j) = Rd(j)
      ENDDO
      DO j = 1, N
        Rd(j) = Wksp(Iwksp(j))
      ENDDO
      DO j = 1, N               !Ditto re.
        Wksp(j) = Re(j)
      ENDDO
      DO j = 1, N
        Re(j) = Wksp(Iwksp(j))
      ENDDO
      END SUBROUTINE sort5

!***********************************************************************
!     glacr_restart- write or read glacrrestart file
!***********************************************************************
      SUBROUTINE glacr_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit
      USE PRMS_GLACR
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL :: check_restart
      ! Local Variable
      CHARACTER(LEN=10) :: module_name
!***********************************************************************
      IF ( In_out==0 ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Nhrugl, Basin_gl_top_melt, Gl_mb_yrcumul
        WRITE ( Restart_outunit ) Gl_mb_cumul, Glnet_ar_delta, Gl_mbc_yrend
        WRITE ( Restart_outunit ) Basin_gl_top_gain
        WRITE ( Restart_outunit ) Basin_gl_area, Gl_area, Basin_gl_ice_melt
        WRITE ( Restart_outunit ) Hru_glres_melt, Basin_gl_storstart
        WRITE ( Restart_outunit ) Gl_top_melt, Basin_gl_storage, Basin_gl_storvol
        WRITE ( Restart_outunit ) Basal_elev
        WRITE ( Restart_outunit ) Keep_gl
        WRITE ( Restart_outunit ) Ikeep_gl
        WRITE ( Restart_outunit ) Basal_slope
        WRITE ( Restart_outunit ) Av_basal_slope
        WRITE ( Restart_outunit ) Av_fgrad
        WRITE ( Restart_outunit ) Prev_out, Prev_outi
        WRITE ( Restart_outunit ) Hru_mb_yrend
        WRITE ( Restart_outunit ) Glacr_flow
        WRITE ( Restart_outunit ) Gl_ice_melt
        WRITE ( Restart_outunit ) Top
        WRITE ( Restart_outunit ) Term
        WRITE ( Restart_outunit ) Ela
        WRITE ( Restart_outunit ) Order_flowline
        WRITE ( Restart_outunit ) Ode_glacrva_coef
        WRITE ( Restart_outunit ) Top_tag
        WRITE ( Restart_outunit ) Glacr_tag
        WRITE ( Restart_outunit ) Delta_volyr
        WRITE ( Restart_outunit ) Hru_mb_yrcumul
        WRITE ( Restart_outunit ) Hru_slope_ts
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Nhrugl, Basin_gl_top_melt, Gl_mb_yrcumul
        READ ( Restart_inunit ) Gl_mb_cumul, Glnet_ar_delta, Gl_mbc_yrend
        READ ( Restart_inunit ) Basin_gl_top_gain
        READ ( Restart_inunit ) Basin_gl_area, Gl_area, Basin_gl_ice_melt
        READ ( Restart_inunit ) Hru_glres_melt, Basin_gl_storstart
        READ ( Restart_inunit ) Gl_top_melt, Basin_gl_storage, Basin_gl_storvol
        READ ( Restart_inunit ) Basal_elev
        READ ( Restart_inunit ) Keep_gl
        READ ( Restart_inunit ) Ikeep_gl
        READ ( Restart_inunit ) Basal_slope
        READ ( Restart_inunit ) Av_basal_slope
        READ ( Restart_inunit ) Av_fgrad
        READ ( Restart_inunit ) Prev_out, Prev_outi
        READ ( Restart_inunit ) Hru_mb_yrend
        READ ( Restart_inunit ) Glacr_flow
        READ ( Restart_inunit ) Gl_ice_melt
        READ ( Restart_inunit ) Top
        READ ( Restart_inunit ) Term
        READ ( Restart_inunit ) Ela
        READ ( Restart_inunit ) Order_flowline
        READ ( Restart_inunit ) Ode_glacrva_coef
        READ ( Restart_inunit ) Top_tag
        READ ( Restart_inunit ) Glacr_tag
        READ ( Restart_inunit ) Delta_volyr
        READ ( Restart_inunit ) Hru_mb_yrcumul
        READ ( Restart_inunit ) Hru_slope_ts
      ENDIF
      END SUBROUTINE glacr_restart
