!***********************************************************************
! Routes water between segments in the system as inflow equals outflow
!***********************************************************************
      INTEGER FUNCTION strmflow_in_out()
      USE PRMS_CONSTANTS, ONLY: RUN, DECL, DEBUG_less, OUTFLOW_SEGMENT, CFS2CMS_CONV
      USE PRMS_MODULE, ONLY: Nsegment, Process_flag, Print_debug
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_BASIN, ONLY: Active_area
      USE PRMS_GWFLOW, ONLY: Basin_gwflow
      USE PRMS_FLOWVARS, ONLY: Basin_ssflow, Basin_cfs, Basin_cms, Basin_stflow_in, &
     &    Basin_sroff_cfs, Basin_ssflow_cfs, Basin_gwflow_cfs, Basin_stflow_out, &
     &    Seg_inflow, Seg_outflow, Seg_upstream_inflow, Seg_lateral_inflow, Flow_out, Basin_sroff
      USE PRMS_ROUTING, ONLY: Obsin_segment, Segment_order, Tosegment, Obsout_segment, Segment_type, &
     &    Flow_to_lakes, Flow_to_ocean, Flow_to_great_lakes, Flow_out_region, Flow_replacement, &
     &    Flow_out_NHM, Flow_terminus, Flow_in_region, Flow_in_nation, Flow_headwater, Flow_in_great_lakes
      USE PRMS_OBS, ONLY: Streamflow_cfs
      IMPLICIT NONE
! Functions
      EXTERNAL :: print_module
! Local Variables
      character(len=*), parameter :: MODDESC = 'Streamflow Routing'
      character(len=*), parameter :: MODNAME = 'strmflow_in_out'
      character(len=*), parameter :: Version_strmflow = '2020-08-03'
      INTEGER :: i, iorder, toseg, segtype
      DOUBLE PRECISION :: area_fac, segout
!***********************************************************************
      strmflow_in_out = 0

      IF ( Process_flag==RUN ) THEN
        Seg_inflow = 0.0D0
        Seg_outflow = 0.0D0
        Seg_upstream_inflow = 0.0D0
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
        DO i = 1, Nsegment
          iorder = Segment_order(i)
          toseg = Tosegment(iorder)
          segtype = Segment_type(iorder)
          IF ( Obsin_segment(iorder)>0 ) Seg_upstream_inflow(iorder) = Streamflow_cfs(Obsin_segment(iorder))
          Seg_inflow(iorder) = Seg_upstream_inflow(iorder) + Seg_lateral_inflow(iorder)
          IF ( Obsout_segment(iorder)>0 ) THEN
            Seg_outflow(iorder) = Streamflow_cfs(Obsout_segment(iorder))
          ELSE
            Seg_outflow(iorder) = Seg_inflow(iorder)
          ENDIF

          IF ( Seg_outflow(iorder) < 0.0 ) THEN
            IF ( Print_debug>DEBUG_less ) THEN
              PRINT *, 'WARNING, negative flow from segment:', iorder, ' flow:', Seg_outflow(iorder)
              PRINT *, '         likely a water-use specification or replacement flow issue'
            ENDIF
          ENDIF

          segout = Seg_outflow(iorder)
! Flow_out is the total flow out of the basin, which allows for multiple outlets
! includes closed basins (tosegment=0)
          IF ( segtype==1 ) THEN
            Flow_headwater = Flow_headwater + segout
          ELSEIF ( segtype==2 ) THEN
            Flow_to_lakes = Flow_to_lakes + segout
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
          IF ( toseg==OUTFLOW_SEGMENT ) THEN
            Flow_out = Flow_out + segout
          ELSE
            Seg_upstream_inflow(toseg) = Seg_upstream_inflow(toseg) + segout
          ENDIF
        ENDDO

        area_fac = Cfs_conv*Active_area
        Basin_stflow_in = Basin_sroff + Basin_gwflow + Basin_ssflow ! not equal to basin_stflow_out if replacement flows
        Basin_cfs = Flow_out
        Basin_stflow_out = Basin_cfs/area_fac
        Basin_cms = Basin_cfs*CFS2CMS_CONV
        Basin_sroff_cfs = Basin_sroff*area_fac
        Basin_ssflow_cfs = Basin_ssflow*area_fac
        Basin_gwflow_cfs = Basin_gwflow*area_fac

      ELSEIF ( Process_flag==DECL ) THEN
        CALL print_module(MODDESC, MODNAME, Version_strmflow)
      ENDIF

      END FUNCTION strmflow_in_out
