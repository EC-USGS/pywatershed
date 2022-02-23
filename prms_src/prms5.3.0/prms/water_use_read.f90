!***********************************************************************
! Read and makes available water-use data (diversions and gains)
! from files pre-processed Data Files available for other PRMS modules
!***********************************************************************
      MODULE PRMS_WATER_USE
      USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH
      IMPLICIT NONE
      ! Local Variables
      character(len=*), parameter :: MODDESC = 'Time Series Data'
      character(len=*), parameter :: MODNAME = 'water_use_read'
      character(len=*), parameter :: Version_water_use_read = '2021-11-19'

      ! transfer type
      integer, parameter :: STREAM = 1
      integer, parameter :: GROUNDWATER = 2
      integer, parameter :: DPRST = 3
      integer, parameter :: EXTRNAL = 4
      integer, parameter :: LAKE = 5
      integer, parameter :: CAPILLARY = 6
      integer, parameter :: CONSUMPTIVE = 7
      integer, parameter :: CANOPY = 8

      ! Declared Variables
      DOUBLE PRECISION, SAVE :: Total_external_transfer, Total_external_gain
      REAL, ALLOCATABLE, SAVE :: External_transfer(:), External_gain(:), External_transfer_tot(:), External_gain_tot(:)
      DOUBLE PRECISION, SAVE ::  Total_consumed_gain
      REAL, ALLOCATABLE, SAVE :: Consumed_gain(:), Consumed_gain_tot(:)
      DOUBLE PRECISION, SAVE ::  Total_soilzone_gain
      REAL, ALLOCATABLE, SAVE :: Soilzone_gain(:), Soilzone_gain_tot(:), Soilzone_gain_hru(:)
      DOUBLE PRECISION, SAVE :: Total_dprst_transfer, Total_dprst_gain
      REAL, ALLOCATABLE, SAVE :: Dprst_transfer(:), Dprst_gain(:), Dprst_transfer_tot(:), Dprst_gain_tot(:)
      DOUBLE PRECISION, SAVE :: Total_gwr_transfer, Total_gwr_gain
      REAL, ALLOCATABLE, SAVE :: Gwr_transfer(:), Gwr_gain(:), Gwr_transfer_tot(:), Gwr_gain_tot(:)
      DOUBLE PRECISION, SAVE :: Total_segment_transfer, Total_segment_gain
      REAL, ALLOCATABLE, SAVE :: Segment_transfer(:), Segment_gain(:), Segment_transfer_tot(:), Segment_gain_tot(:)
      DOUBLE PRECISION, SAVE :: Total_lake_transfer, Total_lake_gain
      REAL, ALLOCATABLE, SAVE :: Lake_transfer(:), Lake_gain(:), Lake_transfer_tot(:), Lake_gain_tot(:)
      DOUBLE PRECISION, SAVE ::  Total_canopy_gain
      REAL, SAVE, ALLOCATABLE :: Canopy_gain(:), Canopy_gain_tot(:)
      DOUBLE PRECISION, SAVE ::  Total_transfers
      REAL, ALLOCATABLE, SAVE :: Transfer_rate(:)
      ! Local Variables
      INTEGER, SAVE :: Outunit, Ndiversions
      INTEGER, SAVE :: Consumed_transfers_on, Lake_transfers_on, Segment_transfers_on
      INTEGER, SAVE :: External_transfers_on, Dprst_transfers_on, Gwr_transfers_on
      INTEGER, ALLOCATABLE, SAVE :: Source_id(:), Destination_id(:), Source_type(:), Destination_type(:)
! Control Parameters
      CHARACTER(LEN=MAXFILE_LENGTH) :: Segment_transfer_file, Gwr_transfer_file, Dprst_transfer_file
      CHARACTER(LEN=MAXFILE_LENGTH) :: External_transfer_file, Lake_transfer_file
      END MODULE PRMS_WATER_USE

      INTEGER FUNCTION water_use_read()
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ACTIVE, OFF, &
     &    RUN, DECL, INIT, CLEAN, ERROR_water_use, strmflow_noroute_module, strmflow_muskingum_lake_module
     use PRMS_CONTROL_FILE, only: control_string
     use PRMS_MMFAPI, only: declvar_dble, declvar_real
     use PRMS_READ_PARAM_FILE, only: decldim, getdim
      USE PRMS_MODULE, ONLY: Process_flag, Nhru, Nsegment, Nwateruse, Nexternal, Nconsumed, &
     &    Segment_transferON_OFF, Gwr_transferON_OFF, Lake_transferON_OFF, &
     &    External_transferON_OFF, Dprst_transferON_OFF, Dprst_flag, Strmflow_flag, &
     &    Model, Inputerror_flag, Start_year, Start_month, Start_day, Soilzone_add_water_use, &
     &    End_year, End_month, End_day, Dprst_transfer_water_use, Dprst_add_water_use, &
     &    Gwr_transfer_water_use, Gwr_add_water_use, Lake_transfer_water_use, Lake_add_water_use
      USE PRMS_WATER_USE
      use prms_utils, only: error_stop, find_current_file_time, find_header_end, print_module, PRMS_open_module_file, read_error
      IMPLICIT NONE
! Functions
      INTRINSIC :: SNGL, DBLE
      EXTERNAL :: read_event
! Local Variables
      INTEGER, SAVE :: external_unit, external_next_year, external_next_month, external_next_day
      INTEGER, SAVE :: segment_unit, dprst_unit, gwr_unit, lake_unit
      INTEGER :: year, month, day, ierr, istop, i, id_src, id_dest
      INTEGER, SAVE :: dprst_next_year, dprst_next_month, dprst_next_day
      INTEGER, SAVE :: gwr_next_year, gwr_next_month, gwr_next_day
      INTEGER, SAVE :: segment_next_year, segment_next_month, segment_next_day
      INTEGER, SAVE :: lake_next_year, lake_next_month, lake_next_day
      DOUBLE PRECISION :: transfer_rate_dble !, factor
!***********************************************************************
      ! Types
      ! (1) stream segments; (2) groundwater reservoirs; (3) surface-depression storage;
      ! (4) external locations; (5) lakes; (6) capillary reservoir of the soil zone;
      ! (7) internal consumptive-use locations; and (8) plant canopy.
!***********************************************************************
      water_use_read = 0

      IF ( Process_flag==RUN ) THEN
        IF ( External_transferON_OFF==ACTIVE ) THEN
          CALL read_event(external_unit, EXTRNAL, external_next_year, external_next_month, external_next_day)
          Total_external_transfer = 0.0D0
          External_transfer = 0.0
        ENDIF
        IF ( External_transfers_on==ACTIVE ) THEN
          Total_external_gain = 0.0D0
          External_gain = 0.0
        ENDIF

        IF ( Gwr_transferON_OFF==ACTIVE ) THEN
          CALL read_event(gwr_unit, GROUNDWATER, gwr_next_year, gwr_next_month, gwr_next_day)
          Total_gwr_transfer = 0.0D0
          Gwr_transfer = 0.0
        ENDIF
        IF ( Gwr_transfers_on==ACTIVE ) THEN
          Total_gwr_gain = 0.0D0
          Gwr_gain = 0.0
        ENDIF

        IF ( Dprst_transferON_OFF==ACTIVE ) THEN
          CALL read_event(dprst_unit, DPRST, dprst_next_year, dprst_next_month, dprst_next_day)
          Total_dprst_transfer = 0.0D0
          Dprst_transfer = 0.0
        ENDIF
        IF ( Dprst_transfers_on==ACTIVE ) THEN
          Total_dprst_gain = 0.0D0
          Dprst_gain = 0.0
        ENDIF

        IF ( Segment_transferON_OFF==ACTIVE ) THEN
          CALL read_event(segment_unit, STREAM, segment_next_year, segment_next_month, segment_next_day)
          Total_segment_transfer = 0.0D0
          Segment_transfer = 0.0
        ENDIF
        IF ( Segment_transfers_on==ACTIVE ) THEN
          Total_segment_gain = 0.0D0
          Segment_gain = 0.0
        ENDIF

        IF ( Lake_transferON_OFF==ACTIVE ) THEN
          CALL read_event(lake_unit, LAKE, lake_next_year, lake_next_month, lake_next_day)
          Total_lake_transfer = 0.0D0
          Lake_transfer = 0.0
        ENDIF
        IF ( Lake_transfers_on==ACTIVE ) THEN
          Total_lake_gain = 0.0D0
          Lake_gain = 0.0
        ENDIF

        IF ( Consumed_transfers_on==ACTIVE ) THEN
          Total_consumed_gain = 0.0D0
          Consumed_gain = 0.0
        ENDIF

        Total_soilzone_gain = 0.0D0
        Soilzone_gain = 0.0

        Total_canopy_gain = 0.0D0
        Canopy_gain = 0.0

        Total_transfers = 0.0D0
        DO i = 1, Ndiversions
          id_src = Source_id(i)
          id_dest = Destination_id(i)
          transfer_rate_dble = DBLE( Transfer_rate(i) )
          Total_transfers = Total_transfers + transfer_rate_dble

          IF ( Gwr_transfers_on==ACTIVE ) THEN
            IF ( Source_type(i)==GROUNDWATER ) THEN
              Gwr_transfer(id_src) = Gwr_transfer(id_src) + Transfer_rate(i)
              Gwr_transfer_tot(id_src) = Gwr_transfer_tot(id_src) + Transfer_rate(i)
              Total_gwr_transfer = Total_gwr_transfer + transfer_rate_dble
              Gwr_transfer_water_use = ACTIVE
            ENDIF
            IF ( Destination_type(i)==GROUNDWATER ) THEN
              Gwr_gain(id_dest) = Gwr_gain(id_dest) + Transfer_rate(i)
              Gwr_gain_tot(id_dest) = Gwr_gain_tot(id_dest) + Transfer_rate(i)
              Total_gwr_gain = Total_gwr_gain + transfer_rate_dble
              Gwr_add_water_use = ACTIVE
            ENDIF
          ENDIF

          IF ( Dprst_transfers_on==ACTIVE ) THEN
            IF ( Source_type(i)==DPRST ) THEN
              Dprst_transfer(id_src) = Dprst_transfer(id_src) + Transfer_rate(i)
              Dprst_transfer_tot(id_src) = Dprst_transfer_tot(id_src) + Transfer_rate(i)
              Total_dprst_transfer = Total_dprst_transfer + transfer_rate_dble
              Dprst_transfer_water_use = ACTIVE
            ENDIF
            IF ( Destination_type(i)==DPRST ) THEN
              Dprst_gain(id_dest) = Dprst_gain(id_dest) + Transfer_rate(i)
              Dprst_gain_tot(id_dest) = Dprst_gain_tot(id_dest) + Transfer_rate(i)
              Total_dprst_gain = Total_dprst_gain + transfer_rate_dble
              Dprst_add_water_use = ACTIVE
            ENDIF
          ENDIF

          IF ( Segment_transfers_on==ACTIVE ) THEN
            IF ( Source_type(i)==STREAM ) THEN
              Segment_transfer(id_src) = Segment_transfer(id_src) + Transfer_rate(i)
              Segment_transfer_tot(id_src) = Segment_transfer_tot(id_src) + Transfer_rate(i)
              Total_segment_transfer = Total_segment_transfer + transfer_rate_dble
            ENDIF
            IF ( Destination_type(i)==STREAM ) THEN
              Segment_gain(id_dest) = Segment_gain(id_dest) + Transfer_rate(i)
              Segment_gain_tot(id_dest) = Segment_gain_tot(id_dest) + Transfer_rate(i)
              Total_segment_gain = Total_segment_gain + transfer_rate_dble
            ENDIF
          ENDIF

          IF ( Lake_transfers_on==ACTIVE ) THEN
            IF ( Source_type(i)==LAKE ) THEN
              Lake_transfer(id_src) = Lake_transfer(id_src) + Transfer_rate(i)
              Lake_transfer_tot(id_src) = Lake_transfer_tot(id_src) + Transfer_rate(i)
              Total_lake_transfer = Total_lake_transfer + transfer_rate_dble
              Lake_transfer_water_use = ACTIVE
            ENDIF
            IF ( Destination_type(i)==STREAM ) THEN
              Lake_gain(id_dest) = Lake_gain(id_dest) + Transfer_rate(i)
              Lake_gain_tot(id_dest) = Lake_gain_tot(id_dest) + Transfer_rate(i)
              Total_lake_gain = Total_lake_gain + transfer_rate_dble
              Lake_add_water_use = ACTIVE
            ENDIF
          ENDIF

          IF ( External_transfers_on==ACTIVE ) THEN
            IF ( Source_type(i)==EXTRNAL ) THEN
              External_transfer(id_src) = External_transfer(id_src) + Transfer_rate(i)
              External_transfer_tot(id_src) = External_transfer_tot(id_src) + Transfer_rate(i)
              Total_external_transfer = Total_external_transfer + transfer_rate_dble
            ENDIF
            IF ( Destination_type(i)==EXTRNAL ) THEN
              External_gain(id_dest) = External_gain(id_dest) + Transfer_rate(i)
              External_gain_tot(id_dest) = External_gain_tot(id_dest) + Transfer_rate(i)
              Total_external_gain = Total_external_gain + transfer_rate_dble
            ENDIF
          ENDIF

          IF ( Consumed_transfers_on==ACTIVE ) THEN
            IF ( Destination_type(i)==CONSUMPTIVE ) THEN
              Consumed_gain(id_dest) = Consumed_gain(id_dest) + Transfer_rate(i)
              Consumed_gain_tot(id_dest) = Consumed_gain_tot(id_dest) + Transfer_rate(i)
              Total_consumed_gain = Total_consumed_gain + transfer_rate_dble
            ENDIF
          ENDIF

          IF ( Destination_type(i)==CAPILLARY ) THEN
            Soilzone_gain(id_dest) = Soilzone_gain(id_dest) + Transfer_rate(i)
            Soilzone_gain_tot(id_dest) = Soilzone_gain_tot(id_dest) + Transfer_rate(i)
            Total_soilzone_gain = Total_soilzone_gain + transfer_rate_dble
            Soilzone_add_water_use = ACTIVE
          ENDIF

          IF ( Destination_type(i)==CANOPY ) THEN
            Canopy_gain(id_dest) = Canopy_gain(id_dest) + Transfer_rate(i)
            Canopy_gain_tot(id_dest) = Canopy_gain_tot(id_dest) + Transfer_rate(i)
            Total_canopy_gain = Total_canopy_gain + transfer_rate_dble
          ENDIF
        ENDDO

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_water_use_read)

        Dprst_transfers_on = OFF
        IF ( Dprst_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          IF ( Dprst_transferON_OFF==ACTIVE .OR. Model==DOCUMENTATION ) THEN
            Dprst_transfers_on = ACTIVE
            ALLOCATE ( Dprst_transfer(Nhru) )
            CALL declvar_real(MODNAME, 'dprst_transfer', 'nhru', Nhru, &
     &           'Transfer flow rate from surface-depression storage for each HRU for each time step', &
     &           'cfs', Dprst_transfer)
            ALLOCATE ( Dprst_transfer_tot(Nhru) )
            CALL declvar_real(MODNAME, 'dprst_transfer_tot', 'nhru', Nhru, &
     &           'Transfer flow rate from surface-depression storage for each HRU for the simulation', &
     &           'cfs', Dprst_transfer_tot)
            CALL declvar_dble(MODNAME, 'total_dprst_transfer', 'one', 1, &
     &           'Transfer flow rates from all surface-depression storage for each time step', &
     &           'cfs', Total_dprst_transfer)
          ENDIF
          ALLOCATE ( Dprst_gain(Nhru) )
          CALL declvar_real(MODNAME, 'dprst_gain', 'nhru', Nhru, &
     &         'Transfer gains to surface-depression storage for each HRU for each time step', &
     &         'cfs', Dprst_gain)
          ALLOCATE ( Dprst_gain_tot(Nhru) )
          CALL declvar_real(MODNAME, 'dprst_gain_tot', 'nhru', Nhru, &
     &         'Transfer gains to surface-depression storage for each HRU for the simulation', &
     &         'cfs', Dprst_gain_tot)
          CALL declvar_dble(MODNAME, 'total_dprst_gain', 'one', 1, &
     &         'Transfer gains to all surface-depression storage for each time step', &
     &         'cfs', Total_dprst_gain)
        ELSEIF ( Dprst_transferON_OFF==ACTIVE .AND. Model/=DOCUMENTATION ) THEN
          PRINT *, 'ERROR, specified to transfer water from surface-depression storage when dprst_flag = 0'
          Inputerror_flag = 1
        ENDIF

        Segment_transfers_on = OFF
        IF ( Strmflow_flag>strmflow_noroute_module .OR. Model==DOCUMENTATION ) THEN
          IF ( Segment_transferON_OFF==ACTIVE .OR. Model==DOCUMENTATION ) THEN
            Segment_transfers_on = ACTIVE
            ALLOCATE ( Segment_transfer(Nsegment) )
            CALL declvar_real(MODNAME, 'segment_transfer', 'nsegment', Nsegment, &
     &           'Transfer flow rate from each stream segment for each time step', &
     &           'cfs', Segment_transfer)
            ALLOCATE ( Segment_transfer_tot(Nsegment) )
            CALL declvar_real(MODNAME, 'segment_transfer_tot', 'nsegment', Nsegment, &
     &           'Transfer flow rate from each stream segment for the simulation', &
     &           'cfs', Segment_transfer_tot)
            CALL declvar_dble(MODNAME, 'total_segment_transfer', 'one', 1, &
     &           'Transfer flow rates from all stream segments for each time step', &
     &           'cfs', Total_segment_transfer)
          ENDIF
          ALLOCATE ( Segment_gain(Nsegment) )
          CALL declvar_real(MODNAME, 'segment_gain', 'nsegment', Nsegment, &
     &         'Transfer gains for each stream segment for each time step', &
     &         'cfs', Segment_gain)
          ALLOCATE ( Segment_gain_tot(Nsegment) )
          CALL declvar_real(MODNAME, 'segment_gain_tot', 'nsegment', Nsegment, &
     &         'Transfer gains for each stream segment for the simulation', &
     &         'cfs', Segment_gain_tot)
          CALL declvar_dble(MODNAME, 'total_segment_gain', 'one', 1, &
     &         'Transfer gains to all stream segments for each time step', &
     &         'cfs', Total_segment_gain)
        ELSEIF ( Segment_transferON_OFF==ACTIVE .AND. Model/=DOCUMENTATION ) THEN
          PRINT *, 'ERROR, specified to transfer water from stream segments when they are not present'
          Inputerror_flag = 1
        ENDIF

        Gwr_transfers_on = OFF
        IF ( Gwr_transferON_OFF==ACTIVE .OR. Model==DOCUMENTATION ) THEN
          Gwr_transfers_on = ACTIVE
          ALLOCATE ( Gwr_transfer(Nhru) )
          CALL declvar_real(MODNAME, 'gwr_transfer', 'nhru', Nhru, &
     &         'Transfer flow rate from the groundwater reservoir of each HRU for each time step', &
     &         'cfs', Gwr_transfer)
          ALLOCATE ( Gwr_transfer_tot(Nhru) )
          CALL declvar_real(MODNAME, 'gwr_transfer_tot', 'nhru', Nhru, &
     &         'Transfer flow rate from the groundwater reservoir of each HRU for the simulation', &
     &         'cfs', Gwr_transfer_tot)
          CALL declvar_dble(MODNAME, 'total_gwr_transfer', 'one', 1, &
     &         'Transfer flow rates from all groundwater reservoirs for each time step', &
     &         'cfs', Total_gwr_transfer)
        ENDIF
        ALLOCATE ( Gwr_gain(Nhru) )
        CALL declvar_real(MODNAME, 'gwr_gain', 'nhru', Nhru, &
     &       'Transfer gains to the groundwater reservoir of each HRU for each time step', &
     &       'cfs', Gwr_gain)
        ALLOCATE ( Gwr_gain_tot(Nhru) )
        CALL declvar_real(MODNAME, 'gwr_gain_tot', 'nhru', Nhru, &
     &       'Transfer gains to the groundwater reservoir of each HRU for the simulation', &
     &       'cfs', Gwr_gain_tot)
        CALL declvar_dble(MODNAME, 'total_gwr_gain', 'one', 1, &
     &       'Flow to all groundwater reservoirs for each time step', &
     &       'cfs', Total_gwr_gain)

        Lake_transfers_on = OFF
        IF ( Strmflow_flag==strmflow_muskingum_lake_module .OR. Model==DOCUMENTATION ) THEN
          IF ( Lake_transferON_OFF==ACTIVE .OR. Model==DOCUMENTATION ) THEN
            Lake_transfers_on = ACTIVE
            ALLOCATE ( Lake_transfer(Nhru) )
            CALL declvar_real(MODNAME, 'lake_transfer', 'nhru', Nhru, &
     &           'Transfer flow rate from each lake HRU for each time step', &
     &           'cfs', Lake_transfer)
            ALLOCATE ( Lake_transfer_tot(Nhru) )
            CALL declvar_real(MODNAME, 'lake_transfer_tot', 'nhru', Nhru, &
     &           'Transfer flow rate from each lake HRU for the simulation', &
     &           'cfs', Lake_transfer_tot)
            CALL declvar_dble(MODNAME, 'total_lake_transfer', 'one', 1, &
     &           'Transfer flow rates from all lake HRUs for each time step', &
     &           'cfs', Total_lake_transfer)
          ENDIF
          ALLOCATE ( Lake_gain(Nhru) )
          CALL declvar_real(MODNAME, 'lake_gain', 'nhru', Nhru, &
     &         'Transfer gains to each lake HRU for each time step', &
     &         'cfs', Lake_gain)
          ALLOCATE ( Lake_gain_tot(Nhru) )
          CALL declvar_real(MODNAME, 'lake_gain_tot', 'nhru', Nhru, &
     &         'Transfer gains to each lake HRU for the simulation', &
     &         'cfs', Lake_gain_tot)
          CALL declvar_dble(MODNAME, 'total_lake_gain', 'one', 1, &
     &         'Transfer gains to all lake HRUs for each time step', &
     &         'cfs', Total_lake_gain)
        ELSEIF ( Lake_transferON_OFF==ACTIVE .AND. Model/=DOCUMENTATION ) THEN
          PRINT *, 'ERROR, specified to transfer water from lakes when lake module is not active'
          Inputerror_flag = 1
        ENDIF

        External_transfers_on = OFF
        IF ( (External_transferON_OFF==ACTIVE.AND.Nexternal>0) .OR. Model==DOCUMENTATION ) THEN
          External_transfers_on = ACTIVE
          ALLOCATE ( External_transfer(Nexternal) )
          CALL declvar_real(MODNAME, 'external_transfer', 'nexternal', Nexternal, &
     &         'Transfer flow rate from each external source for each time step', &
     &         'cfs', External_transfer)
          ALLOCATE ( External_transfer_tot(Nexternal) )
          CALL declvar_real(MODNAME, 'external_transfer_tot', 'nexternal', Nexternal, &
     &         'Transfer flow rate from each external source for the simulation', &
     &         'cfs', External_transfer_tot)
          CALL declvar_dble(MODNAME, 'total_external_transfer', 'one', 1, &
     &         'Transfer flow rates from all external sources for each time step', &
     &         'cfs', Total_external_transfer)
        ENDIF
        IF ( Nexternal>0 ) THEN
          ALLOCATE ( External_gain(Nexternal) )
          CALL declvar_real(MODNAME, 'external_gain', 'nexternal', Nexternal, &
     &         'Transfer gains to each external location for each time step', &
     &         'cfs', External_gain)
          ALLOCATE ( External_gain_tot(Nexternal) )
          CALL declvar_real(MODNAME, 'external_gain_tot', 'nexternal', Nexternal, &
     &         'Transfer gains to each external location for each time step', &
     &         'cfs', External_gain_tot)
        ENDIF
        CALL declvar_dble(MODNAME, 'total_external_gain', 'one', 1, &
     &       'Transfer gains to all external locations for each time step', &
     &       'cfs', Total_external_gain)

        Consumed_transfers_on = OFF
        IF ( Nconsumed>0 ) THEN
          Consumed_transfers_on = ACTIVE
          ALLOCATE ( Consumed_gain(Nconsumed) )
          CALL declvar_real(MODNAME, 'consumed_gain', 'nconsumed', Nconsumed, &
     &         'Transfer flow rate to each water-use comsumption destination for each time step', &
     &         'cfs', Consumed_gain)
          ALLOCATE ( Consumed_gain_tot(Nconsumed) )
          CALL declvar_real(MODNAME, 'consumed_gain_tot', 'nconsumed', Nconsumed, &
     &         'Transfer flow rate to each water-use comsumption destination for the simulation', &
     &         'cfs', Consumed_gain_tot)
        ENDIF
        CALL declvar_dble(MODNAME, 'total_consumed_gain', 'one', 1, &
     &       'Transfer flow rates to all water-use comsumption destinations for each time step', &
     &       'cfs', Total_consumed_gain)

        ALLOCATE ( Soilzone_gain(Nhru) )
        CALL declvar_real(MODNAME, 'soilzone_gain', 'nhru', Nhru, &
     &       'Transfer gains to the capillary reservoir within the soilzone for each HRU for each time step', &
     &       'cfs', Soilzone_gain)
        ALLOCATE ( Soilzone_gain_tot(Nhru) )
        CALL declvar_real(MODNAME, 'soilzone_gain_tot', 'nhru', Nhru, &
     &       'Transfer gains to the capillary reservoir within the soilzone for each HRU for the simulation', &
     &       'cfs', Soilzone_gain_tot)
        ALLOCATE ( Soilzone_gain_hru(Nhru) )
        CALL declvar_real(MODNAME, 'soilzone_gain_hru', 'nhru', Nhru, &
     &       'Irrigation added to soilzone as depth over each HRU', &
     &       'inches', Soilzone_gain_hru)
        CALL declvar_dble(MODNAME, 'total_soilzone_gain', 'one', 1, &
     &       'Transfer gains to all capillary reservoirs for each time step', &
     &       'cfs', Total_soilzone_gain)

        ALLOCATE ( Canopy_gain(Nhru) )
        CALL declvar_real(MODNAME, 'canopy_gain', 'nhru', Nhru, &
     &       'Transfer gains to the canopy reservoir for each HRU for each time step', &
     &       'cfs', Canopy_gain)
        ALLOCATE ( Canopy_gain_tot(Nhru) )
        CALL declvar_real(MODNAME, 'canopy_gain_tot', 'nhru', Nhru, &
     &       'Transfer gains to the canopy reservoir for each HRU for the simulation', &
     &       'cfs', Canopy_gain_tot)
        CALL declvar_dble(MODNAME, 'total_canopy_gain', 'one', 1, &
     &       'Transfer gains to all canopy reservoirs for each time step', &
     &       'cfs', Total_canopy_gain)

        ALLOCATE ( Transfer_rate(Nwateruse), Source_id(Nwateruse), Source_type(Nwateruse) )
        ALLOCATE ( Destination_id(Nwateruse), Destination_type(Nwateruse) )
        CALL declvar_dble(MODNAME, 'total_transfers', 'one', 1, &
     &           'Transfer of all water-use transfers for each time step', &
     &           'cfs', Total_transfers)
        CALL declvar_real(MODNAME, 'transfer_rate', 'nwateruse', nwateruse, &
     &           'Transfer of each water-use transfer for each time step', &
     &           'cfs', Transfer_rate)

      ELSEIF ( Process_flag==INIT ) THEN
        Ndiversions = 0
        year = Start_year
        month = Start_month
        day = Start_day
        ierr = 0

        CALL PRMS_open_module_file(Outunit, 'water_use.out')
        WRITE ( Outunit, 10 ) 'Simulation Start Date:', year, month, day, '   End Date:', &
     &                         End_year, End_month, End_day
10      FORMAT ( 'Water Use Summary File', /, 2(A, I5, 2('/',I2.2)), / )

        istop = 0
        IF ( Segment_transferON_OFF==ACTIVE ) THEN ! type STREAM
          IF ( control_string(Segment_transfer_file, 'segment_transfer_file')/=0 ) &
     &         CALL read_error(5, 'segment_transfer_file')
          CALL find_header_end(segment_unit, Segment_transfer_file, 'segment_transfer_file', ierr, 0, 0)
          IF ( ierr==0 ) THEN
            CALL find_current_file_time(segment_unit, year, month, day, segment_next_year, segment_next_month, segment_next_day)
            Total_segment_transfer = 0.0D0
            Segment_transfer = 0.0
            Segment_transfer_tot = 0.0
          ELSE
            istop = 1
          ENDIF
        ENDIF
        IF ( Strmflow_flag>strmflow_noroute_module ) THEN
          Total_segment_gain = 0.0D0
          Segment_gain = 0.0
          Segment_gain_tot = 0.0
        ENDIF

        IF ( Gwr_transferON_OFF==ACTIVE ) THEN ! type GROUNDWATER
          IF ( control_string(Gwr_transfer_file, 'gwr_transfer_file')/=0 ) &
     &         CALL read_error(5, 'gwr_transfer_file')
          CALL find_header_end(gwr_unit, Gwr_transfer_file, 'gwr_transfer_file', ierr, 0, 0)
          IF ( ierr==0 ) THEN
            CALL find_current_file_time(gwr_unit, year, month, day, &
     &                                  gwr_next_year, gwr_next_month, gwr_next_day)
            Total_gwr_transfer = 0.0D0
            Gwr_transfer = 0.0
            Gwr_transfer_tot = 0.0
          ELSE
            istop = 1
          ENDIF
        ENDIF
        Gwr_gain = 0.0
        Gwr_gain_tot = 0.0
        Total_gwr_gain = 0.0D0

        IF ( Dprst_transferON_OFF==ACTIVE ) THEN ! type DPRST
          IF ( control_string(Dprst_transfer_file, 'dprst_transfer_file')/=0 ) CALL read_error(5, 'dprst_transfer_file')
          CALL find_header_end(dprst_unit, Dprst_transfer_file, 'dprst_transfer_file', ierr, 0, 0)
          IF ( ierr==0 ) THEN
            CALL find_current_file_time(dprst_unit, year, month, day, dprst_next_year, dprst_next_month, dprst_next_day)
            Total_dprst_transfer = 0.0D0
            Dprst_transfer = 0.0
            Dprst_transfer_tot = 0.0
          ELSE
            istop = 1
          ENDIF
        ENDIF
        IF ( Dprst_flag==ACTIVE ) THEN
          Dprst_gain = 0.0
          Dprst_gain_tot = 0.0
          Total_dprst_gain = 0.0D0
        ENDIF

        IF ( External_transferON_OFF==ACTIVE ) THEN ! type EXTRNAL
          IF ( control_string(External_transfer_file, 'external_transfer_file')/=0 ) &
     &         CALL read_error(5, 'external_transfer_file')
          CALL find_header_end(external_unit, External_transfer_file, 'external_transfer_file', ierr, 0, 0)
          IF ( ierr==0 ) THEN
            CALL find_current_file_time(external_unit, year, month, day, &
     &                                  external_next_year, external_next_month, external_next_day)
            Total_external_transfer = 0.0D0
            External_transfer = 0.0
            External_transfer_tot = 0.0
          ELSE
            istop = 1
          ENDIF
        ENDIF
        IF ( Nexternal>0 ) THEN
          External_gain = 0.0
          External_gain_tot = 0.0
          Total_external_gain = 0.0D0
        ENDIF

        IF ( Lake_transferON_OFF==ACTIVE ) THEN ! Type LAKE
          IF ( control_string(Lake_transfer_file, 'lake_transfer_file')/=0 ) CALL read_error(5, 'lake_transfer_file')
          CALL find_header_end(lake_unit, Lake_transfer_file, 'lake_transfer_file', ierr, 0, 0)
          IF ( ierr==0 ) THEN
            CALL find_current_file_time(lake_unit, year, month, day, lake_next_year, lake_next_month, lake_next_day)
            Total_lake_transfer = 0.0D0
            Lake_transfer = 0.0
            Lake_transfer_tot = 0.0
          ELSE
            istop = 1
          ENDIF
        ENDIF
        IF ( Strmflow_flag==strmflow_muskingum_lake_module ) THEN
          Total_lake_gain = 0.0D0
          Lake_gain = 0.0
          Lake_gain_tot = 0.0
        ENDIF

        IF ( istop==1 ) CALL error_stop('in water_use_read module', ERROR_water_use)

        ! type CAPILLARY
        Soilzone_gain = 0.0
        Soilzone_gain_tot = 0.0
        Soilzone_gain_hru = 0.0
        Total_soilzone_gain = 0.0D0

        IF ( Consumed_transfers_on==ACTIVE ) THEN ! type CONSUMPTIVE
          Consumed_gain = 0.0
          Consumed_gain_tot = 0.0
          Total_consumed_gain = 0.0D0
        ENDIF

        ! type CANOPY
        Canopy_gain = 0.0
        Canopy_gain_tot = 0.0
        Total_canopy_gain = 0.0D0

        Total_transfers = 0.0D0
        Source_id = 0
        Source_type = 0
        Destination_id = 0
        Destination_type = 0
        Transfer_rate = 0.0
      ENDIF

      END FUNCTION water_use_read

!*****************************
! Read event for a source type
!*****************************
      SUBROUTINE read_event(Iunit, Src_type, Next_yr, Next_mo, Next_day)
      USE PRMS_CONSTANTS, ONLY: ERROR_water_use, ACTIVE, OFF
      USE PRMS_WATER_USE, ONLY: CANOPY
      USE PRMS_MODULE, ONLY: Nowyear, Nowmonth, Nowday
      use prms_utils, only: is_eof
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Iunit, Src_type
      INTEGER, INTENT (INOUT) :: Next_yr, Next_mo, Next_day
! Funcions
      EXTERNAL :: check_event, set_transfers
! Local Variables
      INTEGER src_id, dest_type, dest_id, keep_reading, ignore
      REAL transfer
!*******************************************************************************
      IF ( Next_mo==0 ) RETURN ! already found end of file
      keep_reading = ACTIVE
      DO WHILE ( keep_reading==ACTIVE )
        IF ( Next_yr==Nowyear .AND. Next_mo==Nowmonth .AND. Next_day==Nowday ) THEN
          READ ( Iunit, * ) Next_yr, Next_mo, Next_day, src_id, dest_type, dest_id, transfer
          IF ( dest_type>CANOPY ) THEN
            PRINT *, 'ERROR, destination flag>8 for date:', Next_yr, Next_mo, Next_day, ' destination:', dest_type
            ERROR STOP ERROR_water_use
          ENDIF
          CALL check_event(Src_type, dest_type, src_id, dest_id, ignore)
          IF ( ignore==0 ) CALL set_transfers(Src_type, src_id, dest_type, dest_id, transfer)
          CALL is_eof(Iunit, Next_yr, Next_mo, Next_day)
          IF ( Next_mo==0 ) keep_reading = 0
        ELSE
          keep_reading = OFF
        ENDIF
      ENDDO
      END SUBROUTINE read_event

! ****************************
! check event for validity
! ****************************
      SUBROUTINE check_event(Src_type, Dest_type, Src_id, Dest_id, Ignore)
      USE PRMS_CONSTANTS, ONLY: ERROR_water_use, ACTIVE, OFF, &
          &    strmflow_muskingum_lake_module, strmflow_noroute_module
      USE PRMS_WATER_USE, ONLY: Outunit, Segment_transfers_on, Dprst_transfers_on, Lake_transfers_on, &
     &    Consumed_transfers_on, External_transfers_on, Gwr_transfers_on, &
     &    STREAM, GROUNDWATER, DPRST, EXTRNAL, LAKE, CAPILLARY, CONSUMPTIVE, CANOPY
      USE PRMS_MODULE, ONLY: Segment_transferON_OFF, Gwr_transferON_OFF, Lake_transferON_OFF, &
     &    Dprst_transferON_OFF, External_transferON_OFF, Strmflow_flag, Nexternal, Dprst_flag, Nowyear, Nowmonth, Nowday
      use prms_utils, only: error_stop
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Src_type, Dest_type, Src_id, Dest_id
      INTEGER, INTENT(OUT) :: Ignore
!*******************************************************************************
      WRITE ( Outunit, '(/,A,I5,2("/",I2.2))' ) 'Event date:', Nowyear, Nowmonth, Nowday

      Ignore = 0
       IF ( Src_type>LAKE ) THEN
        PRINT *, 'Invalid src_type:', Src_type
        CALL error_stop('Valid src_type values are 1 to 5', ERROR_water_use)
      ENDIF
      IF ( Src_type==STREAM ) THEN
        IF ( Segment_transferON_OFF==OFF ) THEN
          PRINT *, 'Warning, specified a stream segment transfer, but segment_transferON_OFF=0, transfer ignored'
          Ignore = 1
        ELSEIF ( Segment_transfers_on==OFF ) THEN
          CALL error_stop('specified a stream segment transfer, but stream segments not present in model', ERROR_water_use)
        ENDIF
      ENDIF
      IF ( Dest_type==STREAM ) THEN
        IF ( Strmflow_flag==strmflow_noroute_module ) THEN
          CALL error_stop('specified a transfer to stream segment, but stream segments not present in model', ERROR_water_use)
        ELSE
         Segment_transfers_on = ACTIVE
        ENDIF
      ENDIF
      IF ( Src_type==EXTRNAL ) THEN
        IF ( External_transferON_OFF==OFF ) THEN
          PRINT *, 'Warning, specified a external transfer, but external_transferON_OFF=0, transfer ignored'
          Ignore = 1
        ELSEIF ( External_transfers_on==OFF ) THEN
          CALL error_stop('specified a external transfer, but external locations not present in model', ERROR_water_use)
        ENDIF
      ENDIF
      IF ( Dest_type==EXTRNAL ) THEN
        IF ( Nexternal==0 ) THEN
          CALL error_stop('specified a transfer to external location, but nexternal = 0', ERROR_water_use)
        ELSE
          External_transfers_on = ACTIVE
        ENDIF
      ENDIF
      IF ( Dest_type==CONSUMPTIVE ) THEN
        IF ( Consumed_transfers_on==OFF ) &
     &       CALL error_stop('specified a consumption-use transfer, but consumption locations not present in model', &
     &       ERROR_water_use)
      ENDIF
      IF ( Src_type==GROUNDWATER ) THEN
        IF ( Gwr_transferON_OFF==OFF ) THEN
          PRINT *, 'Warning, specified a groundwater transfer, but gwr_transferON_OFF=0, transfer ignored'
          Ignore = 1
        ENDIF
      ENDIF
      IF ( Dest_type==GROUNDWATER ) Gwr_transfers_on = ACTIVE
      IF ( Src_type==DPRST ) THEN
        IF ( Dprst_transferON_OFF==OFF ) THEN
          PRINT *, 'Warning, specified a external transfer, but dprst_transferON_OFF=0, transfer ignored'
          Ignore = 1
        ELSEIF ( Dprst_transfers_on==OFF ) THEN
          CALL error_stop('specified a surface-depression transfer, but dprst_flag=0', ERROR_water_use)
        ENDIF
      ENDIF
      IF ( Dest_type==DPRST ) THEN
        IF ( Dprst_flag==OFF ) THEN
          CALL error_stop('specified a transfer to depression storage, but dprst_flag = 0', ERROR_water_use)
        ELSE
          Dprst_transfers_on = ACTIVE
        ENDIF
      ENDIF
      IF ( Src_type==LAKE ) THEN
        IF ( Lake_transferON_OFF==OFF ) THEN
          PRINT *, 'Warning, specified a lake transfer, but lake_transferON_OFF=0, transfer ignored'
          Ignore = 1
        ELSEIF ( Lake_transfers_on==OFF ) THEN
          CALL error_stop('specified a lake transfer, but lake module is not active', ERROR_water_use)
        ENDIF
      ENDIF
      IF ( Dest_type==LAKE ) THEN
        IF ( Strmflow_flag/=strmflow_muskingum_lake_module ) THEN
          CALL error_stop('specified a transfer to lake, but lake simulation is inactive', ERROR_water_use)
        ELSE
          Lake_transfers_on = ACTIVE
        ENDIF
      ENDIF
      IF ( Src_type==Dest_type .AND. Dest_id==Src_id ) THEN
        PRINT *, 'Warning, specified source and destination as equal, transfer ignored'
        Ignore = 1
      ENDIF
      IF ( Ignore==1 ) PRINT *, 'src_type=', Src_type, '; dest_type=', Dest_type
      END SUBROUTINE check_event

! ****************************
! reset transfers when new event is found
! ****************************
      SUBROUTINE set_transfers(Src_type, Src_id, Dest_type, Dest_id, Diversion)
      USE PRMS_CONSTANTS, ONLY: ERROR_water_use
      USE PRMS_WATER_USE
      IMPLICIT NONE
! Arguments
      INTEGER, INTENT(IN) :: Src_type, Src_id, Dest_type, Dest_id
      REAL, INTENT(IN) :: Diversion
! Functions
      EXTERNAL :: check_transfer
!*******************************************************************************
      IF ( Src_type==EXTRNAL ) THEN
        CALL check_transfer('external location', Src_type, Src_id, Dest_type, Dest_id, Diversion)
      ELSEIF ( Src_type==GROUNDWATER ) THEN
        CALL check_transfer('groundwater reservior', Src_type, Src_id, Dest_type, Dest_id, Diversion )
      ELSEIF ( Src_type==DPRST ) THEN
        CALL check_transfer('open surface-depression storage', Src_type, Src_id, Dest_type, Dest_id, Diversion)
      ELSEIF ( Src_type==STREAM ) THEN
        CALL check_transfer('stream segment', Src_type, Src_id, Dest_type, Dest_id, Diversion)
      ELSEIF ( Src_type==LAKE ) THEN
        CALL check_transfer('lake storage', Src_type, Src_id, Dest_type, Dest_id, Diversion)
      ELSE
        PRINT *, 'ERROR, invalid src_type:', Src_type
        ERROR STOP ERROR_water_use
      ENDIF
      END SUBROUTINE set_transfers

! ****************************
! check transfer when new event is found
! ****************************
      SUBROUTINE check_transfer(Ctype, Src_type, Src_id, Dest_type, Dest_id, Diversion)
      USE PRMS_WATER_USE
      USE PRMS_MODULE, ONLY: Nwateruse
      IMPLICIT NONE
      ! Functions
      EXTERNAL :: nwateruse_error
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: Ctype
      INTEGER, INTENT(IN) :: Src_type, Src_id, Dest_type, Dest_id
      REAL, INTENT(IN) :: Diversion
      ! Local Variables
      INTEGER :: i, done
      ! ***********************************
      i = 0 ! index of all transfers
      done = 0
      DO WHILE ( done==0 )
        i = i + 1
        IF ( i>Nwateruse ) CALL nwateruse_error(Ctype)
        ! is this a new diversion
        IF ( Source_type(i)==0 ) THEN ! new
          Source_id(i) = Src_id
          Source_type(i) = Src_type
          Destination_id(i) = Dest_id
          Destination_type(i) = Dest_type
          Transfer_rate(i) = Diversion
          Ndiversions = Ndiversions + 1
          done = 1
        ELSEIF ( Source_type(i)==Src_type .AND. Source_id(i)==Src_id .AND. &
     &           Destination_type(i)==Dest_type .AND. Destination_id(i)==Dest_id ) THEN
          ! replace old value with new value
          Transfer_rate(i) = Diversion
          done = 1
        ENDIF
      ENDDO
      WRITE ( Outunit, '(3A,1X,I0)' ) 'Source: ', Ctype, ':', Src_id
      IF ( Dest_type==STREAM ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Destination: stream segment:', Dest_id
      ELSEIF ( Dest_type==GROUNDWATER ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Destination: groundwater reservior, HRU:', Dest_id
      ELSEIF ( Dest_type==DPRST ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Destination: open surface-depression storage, HRU:', Dest_id
      ELSEIF ( Dest_type==EXTRNAL ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Destination: external location:', Dest_id
      ELSEIF ( Dest_type==LAKE ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Lake storage, HRU:', Dest_id
      ELSEIF ( Dest_type==CAPILLARY ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Capillary reservoir storage, HRU:', Dest_id
      ELSEIF ( Dest_type==CONSUMPTIVE ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Consumptive-use location:', Dest_id
      ELSEIF ( Dest_type==CANOPY ) THEN
        WRITE ( Outunit, '(A,1X,I0)' ) 'Canopy storage, HRU:', Dest_id
      ENDIF
      WRITE ( Outunit, '(A,1X,F0.5)' ) 'Transfer flow rate:', Diversion
      END SUBROUTINE check_transfer

      SUBROUTINE nwateruse_error(ctype)
      USE PRMS_CONSTANTS, ONLY: ERROR_water_use
      USE PRMS_MODULE, ONLY: Nwateruse, Nowyear, Nowmonth, Nowday
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN) :: ctype
      ! ***********************************
      PRINT ('(3A, I0)'), 'ERROR, too many water-use ', ctype, ': specified nwateruse= ', Nwateruse
      PRINT *, '       increase nwateruse to number of unique transfers'
      PRINT *, Nowyear, Nowmonth, Nowday
      ERROR STOP ERROR_water_use
      END SUBROUTINE nwateruse_error
