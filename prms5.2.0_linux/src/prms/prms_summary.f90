!***********************************************************************
!     WRITES NHM CSV SUMMARY FILE
!***********************************************************************
      MODULE PRMS_PRMS_SUMMARY
        USE PRMS_CONSTANTS, ONLY: MAXFILE_LENGTH, RUN, DECL, INIT, CLEAN, ACTIVE, OFF, &
     &      DOCUMENTATION, ERROR_open_out, ERROR_OPEN_IN, ERROR_READ, MAXDIM
        USE PRMS_MODULE, ONLY: Model, Process_flag, Nobs, Nsegment, Npoigages, &
     &      Csv_output_file, Inputerror_flag, Parameter_check_flag, CsvON_OFF
        IMPLICIT NONE
        ! Local Variables
        character(len=*), parameter :: MODDESC = 'Output Summary'
        character(len=*), parameter :: MODNAME = 'prms_summary'
        character(len=*), parameter :: Version_prms_summary = '2020-12-02'
        INTEGER, PARAMETER :: NVARS = 51
        INTEGER, SAVE :: Iunit
        INTEGER, SAVE, ALLOCATABLE :: Gageid_len(:)
        DOUBLE PRECISION, SAVE, ALLOCATABLE :: Segmentout(:) !, Gageout(:)
        CHARACTER(LEN=40), ALLOCATABLE :: Streamflow_pairs(:)
        CHARACTER(LEN=4), ALLOCATABLE :: Cfs_strings(:)
        CHARACTER(LEN=24) :: Fmt
        CHARACTER(LEN=40), SAVE :: Fmt2
        ! Declared Variables
        DOUBLE PRECISION, SAVE :: Basin_total_storage, Basin_surface_storage
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Poi_gage_segment(:)
!        INTEGER, SAVE, ALLOCATABLE :: Parent_poigages(:)
        CHARACTER(LEN=16), SAVE, ALLOCATABLE :: Poi_gage_id(:)
      END MODULE PRMS_PRMS_SUMMARY

      SUBROUTINE prms_summary()
      USE PRMS_PRMS_SUMMARY
      USE PRMS_CLIMATEVARS, ONLY: Basin_potet, Basin_tmax, Basin_tmin, Basin_swrad, Basin_ppt
      USE PRMS_FLOWVARS, ONLY: Basin_soil_moist, Basin_ssstor, Basin_soil_to_gw, &
     &    Basin_lakeevap, Basin_perv_et, Basin_actet, Basin_lake_stor, &
     &    Basin_gwflow_cfs, Basin_sroff_cfs, Basin_ssflow_cfs, Basin_cfs, Basin_stflow_in, &
     &    Basin_stflow_out, Seg_outflow
      USE PRMS_SET_TIME, ONLY : Nowyear, Nowmonth, Nowday
      USE PRMS_OBS, ONLY: Streamflow_cfs
      USE PRMS_INTCP, ONLY: Basin_intcp_evap, Basin_intcp_stor
      USE PRMS_SNOW, ONLY: Basin_pweqv, Basin_snowevap, Basin_snowmelt, Basin_snowcov, Basin_pk_precip
      USE PRMS_SRUNOFF, ONLY: Basin_imperv_stor, Basin_dprst_evap, Basin_imperv_evap, Basin_dprst_seep, &
     &    Basin_dprst_volop, Basin_dprst_volcl, Basin_hortonian
      USE PRMS_SOILZONE, ONLY: Basin_capwaterin, Basin_pref_flow_infil, Basin_prefflow, Basin_recharge, Basin_slowflow, &
     &    Basin_pref_stor, Basin_slstor, Basin_soil_rechr, Basin_sz2gw, Basin_dunnian
      USE PRMS_GWFLOW, ONLY: Basin_gwstor, Basin_gwin, Basin_gwsink, Basin_gwflow, &
     &    Basin_gwstor_minarea_wb, Basin_dnflow
      IMPLICIT NONE
! Functions
      INTRINSIC :: CHAR, INDEX, MAX
      INTEGER, EXTERNAL :: declparam, declvar, getparam !, control_integer
      EXTERNAL :: read_error, PRMS_open_output_file, print_module, statvar_to_csv, checkdim_bounded_limits
      INTEGER, EXTERNAL :: getparamstring, control_string
! Local Variables
      INTEGER :: i, ios, foo, idim !, statsON_OFF
      DOUBLE PRECISION :: gageflow
      CHARACTER(LEN=10) :: chardate
!***********************************************************************
      IF ( Process_flag==RUN ) THEN
        DO i = 1, Npoigages
          Segmentout(i) = Seg_outflow(Poi_gage_segment(i))
!          Gageout(i) = Streamflow_cfs(Parent_poigages(i))
        ENDDO

        gageflow = 0.0D0
        IF ( Nobs>0 ) gageflow = Streamflow_cfs(1)
        Basin_total_storage = Basin_soil_moist + Basin_intcp_stor + Basin_gwstor + Basin_ssstor + Basin_pweqv + &
     &                        Basin_imperv_stor + Basin_lake_stor + Basin_dprst_volop + Basin_dprst_volcl
        Basin_surface_storage = Basin_intcp_stor + Basin_pweqv + Basin_imperv_stor + Basin_lake_stor + &
     &                          Basin_dprst_volop + Basin_dprst_volcl
        IF ( CsvON_OFF==ACTIVE ) THEN
          WRITE ( chardate, '(I4.4,2("-",I2.2))' ) Nowyear, Nowmonth, Nowday
          WRITE ( Iunit, Fmt2 ) chardate, &
     &            Basin_potet, Basin_actet, Basin_dprst_evap, Basin_imperv_evap, Basin_intcp_evap, Basin_lakeevap, &
     &            Basin_perv_et, Basin_snowevap, Basin_swrad, Basin_ppt, Basin_pk_precip, &
     &            Basin_tmax, Basin_tmin, Basin_snowcov, &
     &            Basin_total_storage, Basin_surface_storage, &
     &            Basin_dprst_volcl, Basin_dprst_volop, Basin_gwstor, Basin_imperv_stor, Basin_intcp_stor, Basin_lake_stor, &
     &            Basin_pweqv, Basin_soil_moist, Basin_ssstor, &
     &            Basin_pref_stor, Basin_slstor, Basin_soil_rechr, &
     &            Basin_capwaterin, Basin_dprst_seep, Basin_gwin, Basin_pref_flow_infil, Basin_recharge, Basin_snowmelt, &
     &            Basin_soil_to_gw, Basin_sz2gw, &
     &            Basin_gwsink, Basin_prefflow, Basin_slowflow, Basin_hortonian, Basin_dunnian, &
     &            Basin_stflow_in, Basin_stflow_out, Basin_gwflow, Basin_dnflow, &
     &            Basin_gwstor_minarea_wb, &
     &            Basin_cfs, Basin_gwflow_cfs, Basin_sroff_cfs, Basin_ssflow_cfs, gageflow, &
     &            (Segmentout(i), i = 1, Npoigages)
!     &            (Segmentout(i), Gageout(i), i = 1, Npoigages)
        ELSE
          WRITE ( chardate, '(I4.4,2(1X,I2.2))' ) Nowyear, Nowmonth, Nowday
          WRITE ( Iunit, Fmt2 ) chardate, (Segmentout(i), i = 1, Npoigages)
        ENDIF

! Declare procedure
      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_prms_summary)

!       Open summary file
        IF ( control_string(Csv_output_file, 'csv_output_file')/=0 ) CALL read_error(5, 'csv_output_file')
        IF ( Model/=DOCUMENTATION ) THEN
          CALL PRMS_open_output_file(Iunit, Csv_output_file, 'csv_output_file', 0, ios)
          IF ( ios/=0 ) ERROR STOP ERROR_open_out
        ENDIF

        IF ( declvar(MODNAME, 'basin_total_storage', 'one', 1, 'double', &
     &       'Basin area-weighted average storage in all water storage reservoirs', &
     &       'inches', Basin_total_storage)/=0 ) CALL read_error(3, 'basin_total_storage')
        IF ( declvar(MODNAME, 'basin_surface_storage', 'one', 1, 'double', &
     &       'Basin area-weighted average storage in all surface water storage reservoirs', &
     &       'inches', Basin_surface_storage)/=0 ) CALL read_error(3, 'basin_surface_storage')

        IF ( Npoigages>0 .OR. Model==DOCUMENTATION ) THEN
!          ALLOCATE ( Parent_poigages(Npoigages) )
!          IF ( declparam(MODNAME, 'parent_poigages', 'npoigages', 'integer', &
!     &         '1', '1', '1000000', &
!     &         'Lumen index in parent model','Lumen index in parent model',&
!     &         'none')/=0 ) CALL read_error(1, 'parent_poigages')
          ALLOCATE ( Poi_gage_segment(Npoigages) )
          IF ( declparam(MODNAME, 'poi_gage_segment', 'npoigages', 'integer', &
     &         '0', 'bounded', 'nsegment', &
     &         'Segment index for each POI gage', &
     &         'Segment index for each POI gage', &
     &         'none')/=0 ) CALL read_error(1, 'poi_gage_segment')
          ALLOCATE ( Poi_gage_id(Npoigages) )
          IF ( declparam(MODNAME, 'poi_gage_id', 'npoigages', 'string', &
     &         '0', '0', '9999999', &
     &         'POI Gage ID', 'USGS stream gage for each POI gage', &
     &         'none')/=0 ) CALL read_error(1, 'poi_gage_id')
        ENDIF

! Initialize Procedure
      ELSEIF ( Process_flag==INIT ) THEN
        idim = MAX(1, Npoigages)
        ALLOCATE ( Streamflow_pairs(idim), Cfs_strings(idim), Segmentout(idim), Gageid_len(idim) )
!        ALLOCATE ( Gageout(idim) )
        Streamflow_pairs = ' '
!        Cfs_strings = ',cfs,cfs'
        IF ( CsvON_OFF==ACTIVE ) THEN
          Cfs_strings = ',cfs'
        ELSE
          Cfs_strings = ' cfs'
        ENDIF

        Basin_total_storage = 0.0D0
        Basin_surface_storage = 0.0D0

        IF ( Npoigages>0 ) THEN
!          IF ( getparam(MODNAME, 'parent_poigages', Npoigages, 'integer', Parent_poigages)/=0 ) &
!     &         CALL read_error(2, 'parent_poigages')
          IF ( getparam(MODNAME, 'poi_gage_segment', Npoigages, 'integer', Poi_gage_segment)/=0 ) &
     &         CALL read_error(2, 'poi_gage_segment')
          IF ( Parameter_check_flag>0 ) &
     &      CALL checkdim_bounded_limits('poi_gage_segment', 'nsegment', Poi_gage_segment, Npoigages, 1, Nsegment, Inputerror_flag)
          DO i = 1, Npoigages
            Poi_gage_id(i) = '                '
          ENDDO

          DO i = 1, Npoigages
            foo = getparamstring(MODNAME, 'poi_gage_id', Npoigages, 'string', &
     &            i-1, Poi_gage_id(i))
          ENDDO

          DO i = 1, Npoigages
            IF ( Poi_gage_segment(i)<1 .OR. Poi_gage_segment(i)>Nsegment ) CYCLE
            Gageid_len(i) = INDEX( Poi_gage_id(i), ' ' ) - 1
            IF ( Gageid_len(i)<0 ) Gageid_len(i) = INDEX( Poi_gage_id(i), CHAR(0) ) - 1
!            PRINT *, 'gageid_len ', Gageid_len(i), ' :', Poi_gage_id(i), ':'
            IF ( Gageid_len(i)<1 ) Gageid_len(i) = 0
            IF ( Gageid_len(i)>0 ) THEN
              IF ( Gageid_len(i)>15 ) Gageid_len(i) = 15
              IF ( CsvON_OFF==ACTIVE ) THEN
                WRITE (Streamflow_pairs(i), '(A,I0,2A)' ) ',seg_outflow_', Poi_gage_segment(i), '_gage_', &
     &                                                    Poi_gage_id(i)(:Gageid_len(i))
              ELSE
                WRITE (Streamflow_pairs(i), '(A,I0,2A)' ) ' seg_outflow_', Poi_gage_segment(i), '_gage_', &
     &                                                    Poi_gage_id(i)(:Gageid_len(i))
              ENDIF
              IF ( Poi_gage_segment(i)>9 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>99 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>999 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>9999 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>99999 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>999999 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>9999999 ) Gageid_len(i) = Gageid_len(i) + 1
              IF ( Poi_gage_segment(i)>99999999 ) Gageid_len(i) = Gageid_len(i) + 1
            ELSE
              Gageid_len(i) = -6
              WRITE (Streamflow_pairs(i), '(A,I0)' ) ',seg_outflow_', Poi_gage_segment(i)
            ENDIF
          ENDDO
        ENDIF

        IF ( CsvON_OFF==ACTIVE ) THEN
          WRITE ( Fmt, '(A,I0,A)' ) '( ', Npoigages+14, 'A )'
          WRITE ( Iunit, Fmt ) 'Date,', &
     &            'basin_potet,basin_actet,basin_dprst_evap,basin_imperv_evap,basin_intcp_evap,basin_lakeevap,', &
     &            'basin_perv_et,basin_snowevap,basin_swrad,basin_ppt,basin_pk_precip,', &
     &            'basin_tmax,basin_tmin,basin_snowcov,', &
     &            'basin_total_storage,basin_surface_storage,', &
     &            'basin_dprst_volcl,basin_dprst_volop,basin_gwstor,basin_imperv_stor,basin_intcp_stor,basin_lake_stor,', &
     &            'basin_pweqv,basin_soil_moist,basin_ssstor,', &
     &            'basin_pref_stor,basin_slstor,basin_soil_rechr,', &
     &            'basin_capwaterin,basin_dprst_seep,basin_gwin,basin_pref_flow_in,basin_recharge,basin_snowmelt,', &
     &            'basin_soil_to_gw,basin_sz2gw,', &
     &            'basin_gwsink,basin_prefflow,basin_slowflow,basin_hortonian,basin_dunnian,', &
     &            'basin_stflow_in,basin_stflow_out,basin_gwflow,basin_dnflow,', &
     &            'basin_gwstor_minarea_wb,', &
     &            'basin_cfs,basin_gwflow_cfs,basin_sroff_cfs,basin_ssflow_cfs,runoff_cfs', &
     &            (Streamflow_pairs(i)(:Gageid_len(i)+23), i = 1, Npoigages)

          WRITE ( Iunit, Fmt ) 'year-month-day,', &
     &            'inches/day,inches/day,inches/day,inches/day,inches/day,inches/day,', &
     &            'inches/day,inches/day,Langleys,inches/day,inches/day,', &
     &            'degrees,degrees,fraction,', &
     &            'inches,inches,', &
     &            'inches,inches,inches,inches,inches,inches,', &
     &            'inches,inches,inches,', &
     &            'inches,inches,inches,', &
     &            'inches/day,inches/day,inches/day,inches/day,inches/day,inches/day,', &
     &            'inches/day,inches/day,', &
     &            'inches/day,inches/day,inches/day,inches/day,inches/day,', &
     &            'inches/day,inches/day,inches/day,inches/day,', &
     &            'inches,', &
     &            'cfs,cfs,cfs,cfs,cfs', &
     &            (Cfs_strings(i), i = 1, Npoigages)

          WRITE ( Fmt2, '(A,I0,A)' )  '( A,', Npoigages+NVARS, '(",",F0.4) )'
!        WRITE ( Fmt2, '(A,I0,A)' )  '( A,', 2*Npoigages+NVARS, '(",",SPES10.3) )'
        ELSE
          WRITE ( Fmt, '(A,I0,A)' ) '( ', Npoigages+1, 'A )'
          WRITE ( Iunit, Fmt ) 'Date', &
     &            (Streamflow_pairs(i)(:Gageid_len(i)+20), i = 1, Npoigages)
          WRITE ( Iunit, Fmt ) 'year month day', (Cfs_strings(i), i = 1, Npoigages)
          WRITE ( Fmt2, '(A,I0,A)' )  '( A,', Npoigages, '(1X,F0.4) )'
        ENDIF

      ELSEIF ( Process_flag==CLEAN ) THEN
        !IF ( control_integer(statsON_OFF, 'statsON_OFF')/=0 ) statsON_OFF = 1
        !IF ( statsON_OFF==ACTIVE ) CALL statvar_to_csv()
        CLOSE ( Iunit )
      ENDIF

      END SUBROUTINE prms_summary

!***********************************************************************
!     statvar_to_csv - write a CSV file based on the statvar file
!***********************************************************************
      SUBROUTINE statvar_to_csv()
      USE PRMS_PRMS_SUMMARY
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: control_string, numchars
      EXTERNAL PRMS_open_input_file, PRMS_open_output_file, error_stop
      ! Local Variable
      INTEGER :: inunit, numvariables, ios, i, outunit, ts, yr, mo, day, hr, mn, sec, num
      INTEGER, ALLOCATABLE :: varindex(:), nc(:)
      REAL, ALLOCATABLE :: values(:)
      CHARACTER(LEN=32), ALLOCATABLE :: varname(:)
      CHARACTER(LEN=MAXFILE_LENGTH), SAVE :: statvar_file, statvar_file_csv
      CHARACTER(LEN=10) :: chardate
      CHARACTER(LEN=13) :: fmt5
      CHARACTER(LEN=17) :: fmt3
      CHARACTER(LEN=27) :: fmt6
!***********************************************************************
      IF ( control_string(statvar_file, 'stat_var_file')/=0 ) CALL read_error(5, 'stat_var_file')
      CALL PRMS_open_input_file(inunit, statvar_file, 'stat_var_file', 0, ios)
      IF ( ios/=0 ) CALL error_stop('opening statvar file', ERROR_open_in)
      statvar_file_csv = statvar_file(:numchars(statvar_file))//'.csv'
      CALL PRMS_open_output_file(outunit, statvar_file_csv, 'statvar_csv', 0, ios)
      IF ( ios/=0 ) CALL error_stop('opening statvar CSV file', ERROR_open_out)
      READ ( inunit, * ) numvariables
      ALLOCATE ( varname(numvariables), varindex(numvariables), values(numvariables), nc(numvariables) )
      DO i = 1, numvariables
        READ ( inunit, '(A)', IOSTAT=ios ) varname(i)
        IF ( ios/=0 ) CALL error_stop('reading statvar file', ERROR_read)
        num = numchars(varname(i))
        READ ( varname(i)(num+1:32), '(I5)' ) varindex(i)
        WRITE ( varname(i), '(A,I0)' ) varname(i)(:num)//'_', varindex(i)
        nc(i) = num + 6
      ENDDO
      WRITE ( fmt5, '(A,I0,A)' ) '( A, ', 2*numvariables, 'A )'
      WRITE ( outunit, fmt5 ) 'Date,', ( varname(i)(:nc(i)), ',', i = 1, numvariables )
      WRITE ( fmt3, '(A,I0,A)' ) '(A, ', 2*numvariables, '(I0,A))'
      WRITE ( outunit, fmt3 ) 'date,', ( varindex(i), ',', i = 1, numvariables )
      WRITE ( fmt6, '(A,I0,A)' ) '( A, ', numvariables, '(",",E14.6) )'
      DO WHILE ( ios/=-1 )
        READ ( inunit, *, IOSTAT=ios ) ts, yr, mo, day, hr, mn, sec, (values(i), i = 1, numvariables )
        IF ( ios==-1 ) EXIT
        IF ( ios/=0 ) THEN
          PRINT *, 'ERROR, reading statvar file values, IOSTAT:', ios
          PRINT *, ts, yr, mo, day, hr, 'number of variables:', numvariables
          PRINT *, (values(i), i = 1, numvariables )
          ERROR STOP ERROR_read
        ENDIF
        WRITE ( chardate, '(I0,2("-",I2.2))' )  yr, mo, day
        WRITE ( outunit, fmt6 ) chardate, (values(i), i = 1, numvariables )
      ENDDO
      CLOSE ( outunit )
      END SUBROUTINE statvar_to_csv
