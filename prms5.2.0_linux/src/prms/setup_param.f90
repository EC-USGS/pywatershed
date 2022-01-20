!***********************************************************************
! Declares necessary parameters that are not used by PRMS
!  These are Lumen, GIS, and ncdf parameters that PRMS does not use.
!  They need to come out in the .par_name and default parameter file.
!***********************************************************************

!***********************************************************************
!     Main setup routine
!***********************************************************************
      INTEGER FUNCTION setup()
      USE PRMS_MODULE, ONLY: Process_flag, DECL
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: setupdecl
!***********************************************************************
      setup = 0

      IF ( Process_flag==DECL ) setup = setupdecl()
      END FUNCTION setup

!***********************************************************************
!     setupdecl - set up parameters
!   Declared Parameters
!     poi_gage_id, poi_gage_segment, poi_type, parent_poigages
!     parent_gw, parent_ssr, parent_segment, parent_hru, hru_lon
!***********************************************************************
      INTEGER FUNCTION setupdecl()
      IMPLICIT NONE
      character(len=*), parameter :: MODDESC = 'Parameter Setup'
      character(len=*), parameter :: MODNAME = 'setup_param'
      character(len=*), parameter :: Version_setup = '2020-08-03'
! Functions
      INTEGER, EXTERNAL :: declparam
      EXTERNAL :: read_error, print_module
!***********************************************************************
      setupdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_setup)

      IF (declparam(MODNAME, 'tosegment_nhm', 'nsegment', 'integer', &
     &     '0', '0', '9999999', &
     &     'National Hydrologic Model downstream segment ID', 'National Model Hydrologic downstream segment ID', &
     &     'none') /= 0 ) CALL read_error(1, 'tosegment_nhm')

      IF (declparam(MODNAME, 'nhm_seg', 'nsegment', 'integer', &
     &     '0', '0', '9999999', &
     &     'National Hydrologic Model segment ID', 'National Hydrologic Model segment ID', &
     &     'none') /= 0 ) CALL read_error(1, 'nhm_seg')

      IF (declparam(MODNAME, 'nhm_id', 'nhru', 'integer', &
     &     '1', '1', '9999999', &
     &     'National Hydrologic Model HRU ID', 'National Hydrologic Model HRU ID', &
     &     'none') /= 0 ) CALL read_error(1, 'nhm_id')

      IF ( declparam(MODNAME, 'poi_gage_id', 'npoigages', 'string', &
     &     '0', '0', '9999999', &
     &     'POI Gage ID', 'USGS stream gage for each POI gage', &
     &     'none')/=0 ) CALL read_error(1, 'poi_gage_id')

      IF ( declparam(MODNAME, 'poi_gage_segment', 'npoigages', 'integer', &
     &     '0', 'bounded', 'nsegment', &
     &     'Segment index for each POI gage', &
     &     'Segment index for each POI gage', &
     &     'none')/=0 ) CALL read_error(1, 'poi_gage_segment')

      ! Type code for each point-of-interest (POI) gage (1=USGS streamgage; 2=NWS River Forecast Center Node 
      ! (http://water.weather.gov/ahps/download.php); 3=USGS SPARROW node (https://water.usgs.gov/nawqa/sparrow/);
      ! 4=outlet from or inlet to NHD-defined waterbody; 5=confluence of NHDPlus flowlines exceeding a Strahler stream order (Strahler, 1952, 1957) of four;
      ! 6=connectivity segment not in NHDPlus; 7=split segment from NHDPlus flowline because of travel time > 24 hours;
      ! 8=split segment from NHMPlus flowline because of change in elevation > 500 meters)
      IF ( declparam(MODNAME, 'poi_type', 'npoigages', 'integer', &
     &     '1', '1', '1', &
     &     'Type code for each POI gage', 'Type code for each POI gage', &
     &     'none')/=0 ) CALL read_error(1, 'poi_type')

      IF ( declparam(MODNAME, 'parent_poigages', 'npoigages', 'integer', &
     &     '0', '0', '9999999', &
     &     'Lumen index in parent model for each POI gage','Lumen index in parent model for each POI gage', &
     &     'none')/=0 ) CALL read_error(1, 'parent_poigages')

      IF ( declparam(MODNAME, 'parent_gw', 'ngw', 'integer', &
     &     '1', '1', '9999999', &
     &     'Lumen index in parent model for each GWR','Lumen index in parent model for each GWR', &
     &     'none')/=0 ) CALL read_error(1, 'parent_gw')

      IF ( declparam(MODNAME, 'parent_ssr', 'nssr', 'integer', &
     &     '1', '1', '99999999', &
     &     'Lumen index in parent model for each SSR','Lumen index in parent model  for each SSR', &
     &     'none')/=0 ) CALL read_error(1, 'parent_ssr')

      IF ( declparam(MODNAME, 'parent_segment', 'nsegment', 'integer', &
     &     '0', '0', '9999999', &
     &     'Lumen index in parent model for each segment','Lumen index in parent model for each segment', &
     &     'none')/=0 ) CALL read_error(1, 'parent_segment')

      IF ( declparam(MODNAME, 'parent_hru', 'nhru', 'integer', &
     &     '1', '1', '9999999', &
     &     'Lumen index in parent model for each HRU','Lumen index in parent model for each HRU', &
     &     'none')/=0 ) CALL read_error(1, 'parent_hru')

      IF ( declparam(MODNAME, 'hru_lon', 'nhru', 'real', &
     &     '-105.0', '-360.0', '360.0', &
     &     'HRU longitude', 'Longitude of each HRU', &
     &     'degrees East')/=0 ) CALL read_error(1, 'hru_lon')

      IF ( declparam(MODNAME, 'hru_x', 'nhru', 'real', &
     &     '0.0', '-1.0E7', '1.0E7', &
     &     'X for each HRU (albers)', &
     &     'Longitude (X) for each HRU in albers projection', &
     &     'meters')/=0 ) CALL read_error(1, 'hru_x')

      IF ( declparam(MODNAME, 'hru_y', 'nhru', 'real', &
     &     '0.0', '-1.0E7', '1.0E7', &
     &     'Y for each HRU (albers)', &
     &     'Latitude (Y) for each HRU in albers projection', &
     &     'meters')/=0 ) CALL read_error(1, 'hru_y')

      END FUNCTION setupdecl
