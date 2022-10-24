! Muskingum routing function that calculates the upstream inflow and
! outflow for each segment

! Args:
!     segment_order: segment routing order
!     to_segment: downstream segment for each segment
!     seg_lateral_inflow: segment lateral inflow
!     seg_inflow0: previous segment inflow variable (internal calculations)
!     outflow_ts: outflow timeseries variable (internal calculations)
!     tsi: integer flood wave travel time
!     ts: float version of integer flood wave travel time
!     c0: Muskingum c0 variable
!     c1: Muskingum c1 variable
!     c2: Muskingum c2 variable

! Returns:
!     seg_upstream_inflow: inflow for each segment for the current day
!     seg_inflow0: segment inflow variable
!     seg_inflow: segment inflow variable
!     seg_outflow0: outflow for each segment for the current day
!     seg_outflow: outflow for each segment for the current day
!     inflow_ts: inflow timeseries variable
!     outflow_ts: outflow timeseries variable (internal calculations)
!     seg_current_sum: summation variable

pure subroutine calc_muskingum_mann( &
    nseg, &  ! in
    segment_order, &  ! in
    to_segment, &  ! in
    seg_lateral_inflow, &  ! in
    seg_inflow0_in, &  ! in
    outflow_ts_in, &  ! in
    tsi, &  ! in
    ts, &  ! in
    c0, &  ! in
    c1, &  ! in
    c2, &  ! in
    seg_upstream_inflow, &  ! out
    seg_inflow0, &  ! out
    seg_inflow, &  ! out
    seg_outflow, &  ! out
    inflow_ts, &  ! out
    outflow_ts, &  ! out
    seg_current_sum &  ! out
)

    implicit none
    
    ! Inputs
    integer(kind=4), intent(in) :: nseg
    integer(kind=8), intent(in), dimension(nseg) :: &
        segment_order, to_segment, tsi
    real(kind=8), intent(in), dimension(nseg) :: &
        seg_lateral_inflow, seg_inflow0_in, outflow_ts_in, ts, &
        c0, c1, c2
    
    ! Outputs
    real(kind=8), intent(out), dimension(nseg) :: &
        seg_upstream_inflow,  seg_inflow0,  seg_inflow, seg_outflow, inflow_ts, &
        outflow_ts,  seg_current_sum

    ! Locals
    real(kind=8), dimension(nseg) :: seg_outflow0
    integer(kind=4) :: ihr, jj
    integer(kind=8) :: jseg, to_seg, remainder
    real(kind=8) :: seg_current_inflow

    real(kind=8), parameter :: zero = 0.0000000000000000000000000000000

    ! Initialize variables
   
    ! In to out copies.
    seg_inflow0 = seg_inflow0_in
    outflow_ts = outflow_ts_in
    
    ! initialize variables for the day
    seg_inflow = zero
    seg_outflow = zero
    inflow_ts = zero
    seg_current_sum = zero
    
    do ihr = 0, 23
        seg_upstream_inflow = zero

        do jj = 1, nseg
            jseg = segment_order(jj) + 1  ! convert from zero based in python
            
            ! current inflow to the segment is the time-weighted average
            ! of the outflow of the upstream segments and the lateral HRU
            ! inflow plus any gains
            seg_current_inflow = seg_lateral_inflow(jseg) + seg_upstream_inflow(jseg)

            ! todo: evaluate if obsin_segment needs to be implemented -
            !  would be needed needed if headwater basins are not included
            !  in a simulation
            ! seg_current_inflow += seg_upstream_inflow(jseg)

            seg_inflow(jseg) = seg_inflow(jseg) + seg_current_inflow
            inflow_ts(jseg) = inflow_ts(jseg) + seg_current_inflow
            seg_current_sum(jseg) = seg_current_sum(jseg) + seg_upstream_inflow(jseg)

            remainder = modulo(ihr + 1, tsi(jseg))
            if (remainder == 0) then
                ! segment routed on current hour
                inflow_ts(jseg) = inflow_ts(jseg) / ts(jseg)

                if (tsi(jseg) > 0) then

                    ! todo: evaluated if denormal results should be dealt with
                    
                    ! Muskingum routing equation
                    outflow_ts(jseg) = ( &
                        inflow_ts(jseg) * c0(jseg) &
                        + seg_inflow0(jseg) * c1(jseg) & 
                        + outflow_ts(jseg) * c2(jseg) )
                    
                else
                    ! travel time is 1 hour or less so outflow is set
                    ! equal to the inflow - outflow_ts is the value for
                    ! the previous hour
                    outflow_ts(jseg) = inflow_ts(jseg)
                    
                end if
                        
                ! previous inflow is equal to inflow_ts from the previous
                ! routed time step
                seg_inflow0(jseg) = inflow_ts(jseg)

                ! upstream inflow is used, reset it to zero so a new
                ! average can be calculated next routing time step
                inflow_ts(jseg) = zero

            end if
            
            ! todo: evaluate if obsout_segment needs to be implemented -
            !  would be needed needed fixing ourflow to observed data is
            !  required in a simulation

            ! todo: water use

            ! segment outflow (the mean daily flow rate for each segment)
            ! will be the average of hourly values
            seg_outflow(jseg) = seg_outflow(jseg) + outflow_ts(jseg)

            ! previous segment outflow is equal to the inflow_ts on the
            ! previous routed timestep
            seg_outflow0(jseg) = outflow_ts(jseg)
            
            ! add current time step flow rate to the upstream flow rate
            ! for the segment this segment is connected to
            to_seg = to_segment(jseg) + 1 ! convert from zero based in python
            if (to_seg >= 1) then
                seg_upstream_inflow(to_seg) = &
                    seg_upstream_inflow(to_seg) + outflow_ts(jseg)
            endif

        end do
    end do

    seg_outflow = seg_outflow / 24.0
    seg_inflow = seg_inflow / 24.0
    seg_upstream_inflow = seg_current_sum / 24.0
        
end subroutine calc_muskingum_mann
