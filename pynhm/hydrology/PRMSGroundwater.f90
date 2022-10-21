pure subroutine calc_groundwater( &
    nhru, &  ! in
    gwarea, &  ! in
    soil_to_gw, &  ! in
    ssr_to_gw, &   ! in
    dprst_seep_hru, &  ! in
    gwres_stor_in, &  ! in
    gwflow_coef, &  ! in
    gwsink_coef, &  ! in
    gwres_stor_old, &  ! in
    hru_in_to_cf, &  ! in
    gwres_stor_out, &  ! out
    gwres_flow, &  ! out
    gwres_sink, &  ! out
    gwres_stor_change, &  ! out
    gwres_flow_vol &  ! out
)

    implicit none
    integer(kind=4), intent(in) :: nhru
    real(kind=8), intent(in), dimension(nhru) :: &
        gwarea, soil_to_gw, ssr_to_gw, dprst_seep_hru, &
        gwflow_coef, gwsink_coef, gwres_stor_old, hru_in_to_cf, &
        gwres_stor_in
    real(kind=8), intent(out), dimension(nhru) :: &
        gwres_flow, gwres_sink, &
        gwres_stor_change, gwres_flow_vol, &
        gwres_stor_out

    ! local
    real(kind=8), dimension(nhru) :: &
        l_gwres_stor, l_gwres_flow, l_gwres_sink, &
        soil_to_gw_vol, ssr_to_gw_vol, dprst_seep_hru_vol

    soil_to_gw_vol = soil_to_gw * gwarea
    ssr_to_gw_vol = ssr_to_gw * gwarea
    dprst_seep_hru_vol = dprst_seep_hru * gwarea

    ! todo: what about route order
    l_gwres_stor = gwres_stor_in * gwarea
    l_gwres_stor = l_gwres_stor + &
         soil_to_gw_vol + ssr_to_gw_vol + dprst_seep_hru_vol

    l_gwres_flow = l_gwres_stor * gwflow_coef
    l_gwres_stor = l_gwres_stor - l_gwres_flow

    l_gwres_sink = l_gwres_stor * gwsink_coef
    where(l_gwres_sink > l_gwres_stor)
      l_gwres_sink = l_gwres_stor
    endwhere
    l_gwres_stor = l_gwres_stor - l_gwres_sink

    ! convert most units back to self variables
    ! output variables
    gwres_stor_out = l_gwres_stor / gwarea
    ! for some stupid reason this is left in acre-inches
    gwres_flow = l_gwres_flow / gwarea
    gwres_sink = l_gwres_sink / gwarea

    gwres_stor_change = gwres_stor_out - gwres_stor_old
    gwres_flow_vol = gwres_flow * hru_in_to_cf

end subroutine calc_groundwater
