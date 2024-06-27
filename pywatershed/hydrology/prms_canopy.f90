module canopy

    implicit none
    real(kind=8), parameter :: zero = 0.0000000000000000000000000000000
    real(kind=8), parameter :: one = 1.0000000000000000000000000000000
    
    contains

    pure subroutine calc_canopy( &
        nhru, &  ! in
        cov_type, &  ! in
        covden_sum, &  ! in
        covden_win, &  ! in
        hru_intcpstor, &  ! in
        hru_intcpevap, &  ! in
        hru_ppt, &  ! in
        hru_rain, &  ! in
        hru_snow, &  ! in
        intcp_changeover, &  ! in
        intcp_evap, &  ! in
        intcp_stor, &  ! in
        intcp_transp_on, &  ! in
        net_ppt, &  ! in
        net_rain, &  ! in
        net_snow, &  ! in
        pptmix, & !in
        pk_ice_prev, &  ! in
        freeh2o_prev, &  ! in
        potet, &  ! in
        potet_sublim, &  ! in
        snow_intcp, &  ! in
        srain_intcp, &  ! in
        transp_on, &  ! in
        wrain_intcp, &  ! in
        time_length, &  ! in
        hru_type, &  ! in
        NEARZERO, &  ! in
        DNEARZERO, &  ! in
        BARESOIL, &  ! in
        GRASSES, &  ! in
        LAND, &  ! in
        LAKE, &  ! in
        RAIN, &  ! in
        SNOW, &  ! in
        OFF, &  ! in
        ACTIVE, &  ! in
        intcp_evap_out, &  ! out
        intcp_form_out, &  ! out
        intcp_stor_out, &  ! out
        net_rain_out, &  ! out
        net_snow_out, &  ! out
        pptmix_out, & ! out
        net_ppt_out, &  ! out
        hru_intcpstor_out, &  ! out
        hru_intcpevap_out, &  ! out
        intcp_changeover_out, &  ! out
        intcp_transp_on_out &  ! out
    )

        implicit none
        ! Input Scalars
        integer(kind=4), intent(in) :: &
            nhru, BARESOIL, GRASSES, LAND, LAKE, &
            RAIN, SNOW, OFF, ACTIVE
        real(kind=8), intent(in) :: &
            time_length, NEARZERO, DNEARZERO

        ! Input Vectors
        integer(kind=4), intent(in), dimension(nhru) :: &
            intcp_transp_on, pptmix
        integer(kind=8), intent(in), dimension(nhru) :: &
            cov_type, hru_type

        real(kind=8), intent(in), dimension(nhru) :: &
            covden_sum, covden_win, hru_intcpstor, hru_intcpevap, hru_ppt, &
            hru_rain, hru_snow, intcp_changeover, intcp_evap, intcp_stor,&
            net_ppt, net_rain, net_snow, pk_ice_prev, freeh2o_prev, potet, &
            potet_sublim, snow_intcp, srain_intcp, transp_on, wrain_intcp

        ! Output vectors
        integer(kind=4), intent(out), dimension(nhru) :: &
            intcp_transp_on_out, intcp_form_out, pptmix_out
        real(kind=8), intent(out), dimension(nhru) :: &
            intcp_evap_out, intcp_stor_out, net_rain_out, net_snow_out, net_ppt_out, &
            hru_intcpstor_out, hru_intcpevap_out, intcp_changeover_out

        ! locals
        integer(kind=4) :: ii

        real(kind=8) :: &
            netrain, netsnow, cov, stor_max_rain,  &
            intcpstor, intcpevap, changeover, extra_water, diff, &
            epan_coef, evrn, evsn, zz, dd, &
            last, intcpstor_in, netrain_in, netsnow_in
        
        ! initialize output values ?
        ! intcp_evap_out = intcpevap
        ! intcp_stor_out = intcpstor
        ! net_rain_out = netrain
        ! net_snow_out = netsnow
        ! net_ppt_out = netrain + netsnow
        ! hru_intcpstor_out = intcpstor * cov
        ! hru_intcpevap_out = intcpevap * cov
        intcp_transp_on_out = intcp_transp_on
        pptmix_out = pptmix
        
        ! Begin the HRU LOOP !
        do ii = 1, nhru
            netrain = hru_rain(ii)
            netsnow = hru_snow(ii)

            if (transp_on(ii) == ACTIVE) then
                cov = covden_sum(ii)
                stor_max_rain = srain_intcp(ii)
            else
                cov = covden_win(ii)
                stor_max_rain = wrain_intcp(ii)
            end if

            intcp_form_out(ii) = RAIN
            if (hru_snow(ii) > zero) then
                intcp_form_out(ii) = SNOW
            end if

            intcpstor = intcp_stor(ii)
            intcpevap = zero
            changeover = zero
            extra_water = zero

            ! lake or bare ground hrus
            if ((hru_type(ii) == LAKE) .or. (cov_type(ii) == BARESOIL)) then
                if ((cov_type(ii) == BARESOIL) .and. (intcpstor > zero)) then
                    extra_water = intcp_stor(ii)
                end if
                intcpstor = zero
            end if

            if ((transp_on(ii) == OFF) .and. (intcp_transp_on_out(ii) == ACTIVE)) then
                ! ***** go from summer to winter cover density
                intcp_transp_on_out(ii) = OFF
                if (intcpstor > zero) then
                    diff = covden_sum(ii) - cov
                    changeover = intcpstor * diff
                    if (cov > zero) then
                        if (changeover < zero) then
                            intcpstor = intcpstor * covden_sum(ii) / cov
                            changeover = zero
                        end if
                    else
                        intcpstor = zero
                    end if
                end if

            else if ((transp_on(ii) == ACTIVE) .and. (intcp_transp_on_out(ii) == OFF)) then
                ! ****** go from winter to summer cover density, excess = throughfall
                intcp_transp_on_out(ii) = ACTIVE
                if (intcpstor > zero) then
                    diff = covden_win(ii) - cov
                    changeover = intcpstor * diff
                    if (cov > zero) then
                        if (changeover < zero) then
                            intcpstor = intcpstor * covden_win(ii) / cov
                            changeover = zero
                        end if
                    else
                        intcpstor = zero
                    end if
                end if
            end if

            ! *****Determine the amount of interception from rain
            if ((hru_type(ii) /= LAKE) .and. (cov_type(ii) /= BARESOIL)) then
                if (hru_rain(ii) > zero) then
                    if (cov > zero) then
                        if (cov_type(ii) > GRASSES) then
                            intcpstor_in = intcpstor
                            netrain_in = netrain
                            call intercept( &
                                hru_rain(ii), stor_max_rain, cov, &  ! in
                                intcpstor_in, netrain_in, &  ! in
                                intcpstor, netrain)  ! out
                        else if (cov_type(ii) == GRASSES) then
                            ! if there is no snowpack and no snowfall, then apparently,
                            ! grasses can intercept rain.
                            if ((pk_ice_prev(ii) + freeh2o_prev(ii) < DNEARZERO) &
                                .and. (netsnow < NEARZERO)) then
                                intcpstor_in = intcpstor
                                netrain_in = netrain
                                call intercept( &
                                    hru_rain(ii), stor_max_rain, cov, &  ! in
                                    intcpstor_in, netrain_in, &  ! in
                                    intcpstor, netrain)
                            end if
                        end if
                    end if
                end if
            end if

            ! ******Determine amount of interception from snow
            if (hru_snow(ii) > zero) then
                if (cov > zero) then
                    if (cov_type(ii) > GRASSES) then
                        intcpstor_in = intcpstor
                        netsnow_in = netsnow
                        call intercept( &
                            hru_snow(ii), snow_intcp(ii), cov, &  ! in
                            intcpstor_in, netsnow_in, &  ! in
                            intcpstor, netsnow)  ! out
                        if (netsnow < NEARZERO) then
                            netrain = netrain + netsnow
                            netsnow = zero
                            ! todo: deal with newsnow and pptmix?
                            ! Newsnow(i) = OFF
                            pptmix_out(ii) = 0
                        end if
                    end if
                end if
            end if

            ! todo: canopy application of irrigation water based on irr_type

            ! ******compute evaporation or sublimation of interception

            ! if precipitation assume no evaporation or sublimation
            if (intcpstor > zero) then
                if (hru_ppt(ii) < NEARZERO) then
                    epan_coef = one
                    evrn = potet(ii) / epan_coef
                    evsn = potet(ii) * potet_sublim(ii)

                    ! todo: pan et
                    ! IF (( Use_pandata==ACTIVE ) THEN
                    ! evrn = Pan_evap(Hru_pansta(i))
                    ! IF (( evrn<zero ) evrn = zero
                    ! ENDIF

                    if (intcp_form_out(ii) == SNOW) then
                        zz = intcpstor - evsn
                        if (zz > 0) then
                            intcpstor = zz
                            intcpevap = evsn
                        else
                            intcpevap = intcpstor
                            intcpstor = zero
                        end if
                    else
                        dd = intcpstor - evrn
                        if (dd > zero) then
                            intcpstor = dd
                            intcpevap = evrn
                        else
                            intcpevap = intcpstor
                            intcpstor = zero
                        end if
                    end if
                end if
            end if

            if (intcpevap * cov > potet(ii)) then
                last = intcpevap
                if (cov > zero) then
                    intcpevap = potet(ii) / cov
                else
                    intcpevap = zero
                end if
                intcpstor = intcpstor + last - intcpevap
            end if

            ! Store calculated values in output variables
            intcp_evap_out(ii) = intcpevap
            intcp_stor_out(ii) = intcpstor
            net_rain_out(ii) = netrain
            net_snow_out(ii) = netsnow
            net_ppt_out(ii) = netrain + netsnow
            hru_intcpstor_out(ii) = intcpstor * cov
            hru_intcpevap_out(ii) = intcpevap * cov
            intcp_changeover_out(ii) = changeover + extra_water


        end do

    end subroutine calc_canopy


    pure subroutine intercept( &
        precip, stor_max, cov, intcp_stor, net_precip, &  ! in
        intcp_stor_out, net_precip_out & ! out
    )
        implicit none

        real(kind=8), intent(in) :: &
            precip, stor_max, cov, intcp_stor, net_precip
        real(kind=8), intent(out) :: &
            intcp_stor_out, net_precip_out
        
        net_precip_out = precip * (one - cov)
        intcp_stor_out = intcp_stor + precip
        if (intcp_stor_out > stor_max) then
            net_precip_out = net_precip_out + (intcp_stor_out - stor_max) * cov
            intcp_stor_out = stor_max
        end if

    end subroutine intercept

end module canopy
