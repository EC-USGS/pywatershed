import numpy as np

from pynhm.base.storageUnit import StorageUnit

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import ETType, HruType, SoilType, epsilon, nan, one, zero

ONETHIRD = 1 / 3
TWOTHIRDS = 2 / 3


class PRMSSoilzone(StorageUnit):
    """PRMS soil zone

    Args:
    """

    def __init__(
        self,
        control: Control,
        dprst_evap_hru: adaptable,
        dprst_seep_hru: adaptable,
        hru_impervevap: adaptable,
        hru_intcpevap: adaptable,
        infil_hru: adaptable,  # file by /pynhm/analysis/budget_soilzone.py
        sroff: adaptable,
        potet: adaptable,
        transp_on: adaptable,
        snow_evap: adaptable,
        snowcov_area: adaptable,
        budget_type: str = None,
        verbose: bool = False,
    ) -> "PRMSSoilzone":

        super().__init__(
            control=control,
            verbose=verbose,
        )
        self.name = "PRMSSoilzone"

        self.set_inputs(locals())
        self.set_budget(budget_type)

        return

    @staticmethod
    def get_parameters() -> tuple:
        """Get soil zone parameters

        Returns:
            parameters: input parameters

        """
        return (
            "nhru",
            "nssr",
            "dprst_frac",
            "cov_type",
            "fastcoef_lin",
            "fastcoef_sq",
            "hru_area",
            "hru_percent_imperv",
            "hru_type",
            "pref_flow_den",
            "sat_threshold",
            "slowcoef_lin",
            "slowcoef_sq",
            "soil_moist_max",
            "soil_moist_init_frac",
            "soil_rechr_init_frac",
            "soil_rechr_max_frac",
            "soil_type",
            "soil2gw_max",
            "ssr2gw_exp",
            "ssr2gw_rate",
            "ssstor_init_frac",
        )

    @staticmethod
    def get_inputs() -> tuple:
        """Get soil zone input variables

        Returns:
            variables: input variables

        """
        return (
            "dprst_evap_hru",  # JLM ?? needs this stuff to calculate evap?
            "dprst_seep_hru",
            "hru_impervevap",  # JLM ??
            "hru_intcpevap",  # JLM ???
            "infil_hru",
            # soil_moist_chg => model_runoff%soil_moist_chg, &
            # soil_rechr_chg => model_runoff%soil_rechr_chg, &
            "sroff",
            "potet",
            # hru_ppt => model_precip%hru_ppt, & # JLM ??
            "transp_on",
            "snow_evap",
            "snowcov_area",
        )

    @staticmethod
    def get_restart_variables() -> tuple:
        """Get soil zone restart variables

        Returns:
            variables: restart variables
        """

        return ()

    @staticmethod
    def get_init_values() -> dict:
        """Get soil zone inital values

        Returns:
            dict: inital values for named variables
        """

        return {
            "cap_infil_tot": zero,
            "cap_waterin": zero,
            "dunnian_flow": zero,
            "hru_actet": zero,
            "infil": nan,
            "perv_actet": zero,
            "perv_actet_hru": zero,
            "potet_lower": zero,
            "potet_rechr": zero,
            "pref_flow": zero,
            "pref_flow_in": zero,
            "pref_flow_infil": zero,
            "pref_flow_max": zero,
            "pref_flow_stor": zero,
            "pref_flow_stor_change": nan,
            "pref_flow_stor_prev": nan,
            "pref_flow_thrsh": zero,
            "recharge": zero,
            "slow_flow": zero,
            "slow_stor": zero,
            "slow_stor_change": zero,
            "slow_stor_prev": nan,
            "snow_free": nan,  # sm_soilzone
            "soil_lower": nan,  # completely set later
            "soil_lower_change": nan,
            "soil_lower_change_hru": nan,
            "soil_lower_prev": nan,
            "soil_lower_ratio": zero,
            "soil_lower_max": nan,  # completely set later
            "soil_moist": nan,  # sm_climateflow
            "soil_moist_prev": nan,  # sm_climateflow
            "soil_moist_tot": nan,  # completely set later
            "soil_rechr": nan,  # sm_climateflow
            "soil_rechr_change": nan,  # sm_climateflow
            "soil_rechr_change_hru": nan,  # sm_climateflow
            "soil_rechr_prev": nan,  # sm_climateflow
            "soil_to_gw": zero,
            "soil_to_ssr": zero,
            "soil_zone_max": nan,  # this is completely later
            "ssr_to_gw": zero,
            "ssres_flow": zero,  # todo: privatize keep vol public
            "ssres_flow_vol": nan,
            "ssres_in": zero,
            "ssres_stor": nan,  # sm_soilzone
            "swale_actet": zero,
            "unused_potet": zero,
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "infil_hru",
            ],
            "outputs": [
                "perv_actet_hru",
                "soil_to_gw",
                "ssr_to_gw",
                "slow_flow"
                # "pref_flow",
            ],
            "storage_changes": [
                "soil_rechr_change_hru",
                "soil_lower_change_hru",
                "slow_stor_change",
                # "pref_flow_stor_change",
            ],
        }

    def set_initial_conditions(self):
        """Initialize PRMSSoilzone variables."""

        # Derived parameters
        # JLM: is this awkward here?
        # JLM: it's definitely awkward to edit a parameter. maybe
        #      changes/checks on params should throw errors so
        #      parameter values remain the responsibility of the users
        #      and their deficiencies are transparent?
        self.hru_area_imperv = self.hru_percent_imperv * self.hru_area
        self.hru_area_perv = self.hru_area - self.hru_area_imperv

        self.soil_rechr_max = self.soil_rechr_max_frac * self.soil_moist_max
        # asdf

        # apparent issues with hru_frac_perv
        # this%hru_frac_perv = 1.0 - this%hru_percent_imperv
        self.hru_frac_perv = one - self.hru_percent_imperv

        dprst_area_max = self.dprst_frac * self.hru_area

        wh_active = np.where(self.hru_type != HruType.INACTIVE.value)
        self.hru_area_perv[wh_active] = (
            self.hru_area_perv[wh_active] - dprst_area_max[wh_active]
        )
        # Recompute hru_frac_perv to reflect the depression storage area
        self.hru_frac_perv[wh_active] = (
            self.hru_area_perv[wh_active] / self.hru_area[wh_active]
        )

        self.snow_free = one - self.snowcov_area
        # edit a param
        wh_inactive_or_lake = np.where(
            (self.hru_type == HruType.INACTIVE.value)
            | (self.hru_type == HruType.LAKE.value)
        )
        self.sat_threshold[wh_inactive_or_lake] = zero
        # edit a param
        wh_not_land = np.where(self.hru_type != HruType.LAND.value)
        self.pref_flow_den[wh_not_land] = zero

        # variables
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            # these are set in sm_climateflow
            self.soil_moist[:] = (
                self.soil_moist_init_frac * self.soil_moist_max
            )
            self.soil_rechr[:] = (
                self.soil_rechr_init_frac * self.soil_rechr_max
            )
        else:
            # call ctl_data%read_restart_variable(
            #    'soil_moist', this%soil_moist)
            # call ctl_data%read_restart_variable(
            #    'soil_rechr', this%soil_rechr)
            raise RuntimeError("Soilzone restart capability not implemented")

        # Note that the following is often editing a copy of the parameters,
        # which is a bit odd in that it might be contrary to the users
        # expectations. Move this parameter business to __init__

        # Do these need set on self or added to variables?
        self._grav_dunnian_flow = np.full(self.nhru, zero, "float64")
        self._pfr_dunnian_flow = np.full(self.nhru, zero, "float64")

        # ssres_stor
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            self.ssres_stor = self.ssstor_init_frac * self.sat_threshold
            wh_inactive_or_lake = np.where(
                (self.hru_type == HruType.INACTIVE.value)
                | (self.hru_type == HruType.LAKE.value)
            )
            self.ssres_stor[wh_inactive_or_lake] = zero
            del self.ssstor_init_frac

        else:
            # call ctl_data%read_restart_variable(
            #     'ssres_stor', this%ssres_stor)
            raise RuntimeError("Soilzone restart capability not implemented")

        # Parameter edits in climateflow
        # JLM: some of these should just be runtime/value errors
        # JLM: These are for "ACTIVE and non-lake" hrus....
        # JLM check that.
        self.soil_moist_max = np.where(
            self.soil_moist_max < 0.00001, 0.00001, self.soil_moist_max
        )
        self.soil_rechr_max = np.where(
            self.soil_rechr_max < 0.00001, 0.00001, self.soil_rechr_max
        )
        self.soil_rechr_max = np.where(
            self.soil_rechr_max > self.soil_moist_max,
            self.soil_moist_max,
            self.soil_rechr_max,
        )
        self.soil_rechr = np.where(
            self.soil_rechr > self.soil_rechr_max,
            self.soil_rechr_max,
            self.soil_rechr,
        )
        self.soil_moist = np.where(
            self.soil_moist > self.soil_moist_max,
            self.soil_moist_max,
            self.soil_moist,
        )
        self.soil_rechr = np.where(
            self.soil_rechr > self.soil_moist, self.soil_moist, self.soil_rechr
        )
        self.ssres_stor = np.where(
            self.ssres_stor > self.sat_threshold,
            self.sat_threshold,
            self.ssres_stor,
        )

        # <
        # need to set on swale_limit self? move to variables?
        self._swale_limit = np.full(self.nhru, zero, "float64")
        wh_swale = np.where(self.hru_type == HruType.SWALE.value)
        self._swale_limit[wh_swale] = 3.0 * self.sat_threshold[wh_swale]

        self.pref_flow_thrsh[wh_swale] = self.sat_threshold[wh_swale]
        wh_land = np.where(self.hru_type == HruType.LAND.value)
        self.pref_flow_thrsh[wh_land] = self.sat_threshold[wh_land] * (
            one - self.pref_flow_den[wh_land]
        )
        self.pref_flow_max[wh_land] = (
            self.sat_threshold[wh_land] - self.pref_flow_thrsh[wh_land]
        )

        # Need to set pref_flow_flag on self? or add to variables?
        wh_land_and_prf_den = np.where(
            (self.hru_type == HruType.LAND.value) & (self.pref_flow_den > zero)
        )
        self._pref_flow_flag = np.full(self.nhru, False, dtype=int)
        self._pref_flow_flag[wh_land_and_prf_den] = True

        # can this one be combined with the restart read logic above?
        if self.control.config["init_vars_from_file"] in [0, 2, 5]:
            wh_land_or_swale = np.where(
                (self.hru_type == HruType.LAND.value)
                | (self.hru_type == HruType.SWALE.value)
            )
            # JLM: to verify this
            self.slow_stor[wh_land_or_swale] = np.minimum(
                self.ssres_stor[wh_land_or_swale],
                self.pref_flow_thrsh[wh_land_or_swale],
            )
            self.pref_flow_stor[wh_land_or_swale] = (
                self.ssres_stor[wh_land_or_swale]
                - self.slow_stor[wh_land_or_swale]
            )

        else:
            # call ctl_data%read_restart_variable(
            #    'pref_flow_stor', self.pref_flow_stor)
            # call ctl_data%read_restart_variable('slow_stor', self.slow_stor)
            raise RuntimeError("Soilzone restart capability not implemented")

        # <
        # need to set soil2gw_flag on self? move to variables?
        self._soil2gw_flag = np.full(self.nhru, False, dtype=bool)
        wh_soil2gwmax = np.where(self.soil2gw_max > zero)
        self._soil2gw_flag[wh_soil2gwmax] = True

        self.soil_zone_max = (
            self.sat_threshold + self.soil_moist_max * self.hru_frac_perv
        )
        self.soil_moist_tot = (
            self.ssres_stor + self.soil_moist * self.hru_frac_perv
        )

        self.soil_lower = self.soil_moist - self.soil_rechr
        self.soil_lower_max = self.soil_moist_max - self.soil_rechr_max

        wh_soil_lower_stor = np.where(self.soil_lower_max > zero)
        self.soil_lower_ratio[wh_soil_lower_stor] = (
            self.soil_lower[wh_soil_lower_stor]
            / self.soil_lower_max[wh_soil_lower_stor]
        )

        self._gvr2pfr = np.full(self.nhru, zero, dtype=float)

        return

    def _advance_variables(self) -> None:
        self.pref_flow_stor_prev[:] = self.pref_flow_stor
        self.soil_rechr_prev[:] = self.soil_rechr
        self.soil_lower_prev[:] = self.soil_lower
        self.slow_stor_prev[:] = self.slow_stor
        return

    def _calculate(self, simulation_time):
        """Calculate soil zone for a time step"""

        # JLM: not clear we need this / for GSFlow
        # if self.srunoff_updated_soil:
        #     self.soil_moist = self.soil_moist_chg
        #     self.soil_rechr = self.soil_rechr_chg
        # # <

        # <
        gwin = zero
        # update_potet = 0

        # diagnostic state resets
        self.soil_to_gw[:] = zero
        self.soil_to_ssr[:] = zero
        self.ssr_to_gw[:] = zero
        self.slow_flow[:] = zero
        self.ssres_flow[:] = zero
        self.potet_rechr[:] = zero
        self.potet_lower[:] = zero

        self.snow_free = one - self.snowcov_area

        # Do this here and not in advance as this is not an individual storage
        self.soil_moist_prev[:] = self.soil_moist

        # JLM: ET calculations to be removed from soilzone.
        self.hru_actet = (
            self.hru_impervevap + self.hru_intcpevap + self.snow_evap
        )

        # This is obnoxious. i guess this should be an
        # optional input? should default to zero?
        if self.control.config["dprst_flag"] == 1:
            self.hru_actet = self.hru_actet + self.dprst_evap_hru

        # <
        for hh in range(self.nhru):

            dunnianflw = zero
            dunnianflw_pfr = zero
            dunnianflw_gvr = zero
            interflow = zero

            # JLM: ET calculation to be removed from soilzone.
            avail_potet = np.maximum(zero, self.potet[hh] - self.hru_actet[hh])

            # Add infiltration to soil and compute excess
            # Whole HRU area:
            #   * cap_infil_tot is the depth in whole HRU
            #   * preferential flow reservoir for whole HRU
            #   * gravity reservoir for whole HRU
            #   * upslope flow for whole HRU
            # Pervious HRU area:
            #   * infil: for pervious portion of HRU
            #   * capwatermaxin: capillary reservoir for pervious area

            # note: hru_area_perv[hh] has to be > zero for infil > 0
            #       hru_frac_perv[hh] has to be > zero
            # IE capwater_maxin is just on the pervious area

            #  this is a bit ridiculous until cleaned up... should only
            # be passing HRU quantities not a mix... eventually volumes
            self.infil = self.infil_hru / self.hru_frac_perv[hh]
            capwater_maxin = self.infil[hh]

            # Compute preferential flow and storage, and any dunnian flow
            prefflow = zero
            # if self._pref_flow_flag[hh]:
            #     self.pref_flow_infil[hh] = zero
            #     pref_flow_maxin = zero

            #     if capwater_maxin > zero:
            #         # PRMSIV Step 1 - Partition infil between capillary and
            #         #                 preferential-flow (eqn 1-121)
            #         # pref_flow for whole HRU but capwater is pervious area
            #         # calculations on pervious area
            #         pref_flow_maxin = capwater_maxin * self.pref_flow_den[hh]

            #         # PRMSIV Step 3: no cascades and already normalized to
            #         #                pervious area. (eqn 1-124)
            #         capwater_maxin = capwater_maxin - pref_flow_maxin

            #         # renormalize pref_flow to whole HRU from pervious area
            #         pref_flow_maxin = pref_flow_maxin * self.hru_frac_perv[hh]

            #         # PRMSIV Step 2 - Compute PFR storage, excess to Dunnian
            #         #                 (eqns 1-122 and 1-123)
            #         # Compute contribution to preferential-flow reservoir
            #         # storage
            #         self.pref_flow_stor[hh] = (
            #             self.pref_flow_stor[hh] + pref_flow_maxin
            #         )
            #         dunnianflw_pfr = max(
            #             zero, self.pref_flow_stor[hh] - self.pref_flow_max[hh]
            #         )

            #         if dunnianflw_pfr > zero:
            #             self.pref_flow_stor[hh] = self.pref_flow_max[hh]

            #         # <
            #         self.pref_flow_infil[hh] = pref_flow_maxin - dunnianflw_pfr

            #     # <
            #     self._pfr_dunnian_flow[hh] = dunnianflw_pfr

            # <
            # whole HRU
            self.cap_infil_tot[hh] = capwater_maxin * self.hru_frac_perv[hh]

            # ****** Add infiltration to soil and compute excess
            # Step 3 above
            self.cap_waterin[hh] = capwater_maxin

            # PRMSIV Steps 4, 5, 6 (see compute_soilmoist)
            if (capwater_maxin + self.soil_moist[hh]) > zero:

                # JLM: not sure why the function returns cap_waterin
                (
                    self.cap_waterin[hh],
                    self.soil_moist[hh],
                    self.soil_rechr[hh],
                    self.soil_to_gw[hh],
                    self.soil_to_ssr[hh],
                ) = self.compute_soilmoist(
                    self._soil2gw_flag[hh],
                    self.hru_frac_perv[hh],
                    self.soil_moist_max[hh],
                    self.soil_rechr_max[hh],
                    self.soil2gw_max[hh],
                    self.cap_waterin[hh],
                    self.soil_moist[hh],
                    self.soil_rechr[hh],
                    self.soil_to_gw[hh],
                    self.soil_to_ssr[hh],
                )

                self.cap_waterin[hh] = (
                    self.cap_waterin[hh] * self.hru_frac_perv[hh]
                )

            # <

            topfr = zero

            # soil_to_ssr also known as grv_maxin
            # grv_availh2o = grv_stor + soil_to_ssr
            availh2o = self.slow_stor[hh] + self.soil_to_ssr[hh]

            if self.hru_type[hh] == HruType.LAND.value:
                # PRMSIV Step 7
                # PRMSIV eqn 1-132? not necessary?
                # PRMSIV eqn 1-133:
                #     gvr_excess = max(
                #         0, slow_stor_old + gvr_maxin - pref_flow_thresh)
                topfr = max(zero, availh2o - self.pref_flow_thrsh[hh])

                # topfr = max(
                #    zero,
                #    self.slow_stor[hh]
                #    + self.soil_to_ssr[hh]
                #    - self.pref_flow_thrsh[hh],
                # )

                # PRMSIV eqn 1-134: ssres_in gvr_maxin - gvr_excess
                ssresin = self.soil_to_ssr[hh] - topfr
                # ssresin = self.soil_to_ssr[hh] - max(
                #     zero,
                #     self.slow_stor[hh]
                #     + self.soil_to_ssr[hh]
                #     - self.pref_flow_thrsh[hh],
                # )

                # JLM: what equation is this?
                self.slow_stor[hh] = max(zero, availh2o - topfr)
                # JLM: This is expansion is absurd
                # self.slow_stor[hh] = max(
                #     zero,
                #     self.slow_stor[hh]
                #     + self.soil_to_ssr[hh]
                #     - max(
                #         zero,
                #         self.slow_stor[hh]
                #         + self.soil_to_ssr[hh]
                #         - self.pref_flow_thrsh[hh],
                #     ),
                # )
                # if pref_flow_thrsh < availh2o, then
                # slow_stor = pref_flow_thrsh ???

                # JLM: the order of the PRMSIV theory is not maintained here
                # This skips over 8 to step 9

                # PRMSIV Step 9
                # Compute slow contribution to interflow, if any
                if self.slow_stor[hh] > epsilon:
                    (
                        self.slow_stor[hh],
                        self.slow_flow[hh],
                    ) = self.compute_interflow(
                        self.slowcoef_lin[hh],
                        self.slowcoef_sq[hh],
                        ssresin,
                        self.slow_stor[hh],
                        self.slow_flow[hh],
                    )

                # <
            elif self.hru_type[hh] == HruType.SWALE:
                self.slow_stor[hh] = availh2o

            # <
            if (self.slow_stor[hh] > epsilon) and (
                self.ssr2gw_rate[hh] > zero
            ):
                (
                    self.ssr_to_gw[hh],
                    self.slow_stor[hh],
                ) = self.compute_gwflow(
                    self.ssr2gw_rate[hh],
                    self.ssr2gw_exp[hh],
                    self.slow_stor[hh],
                )

            # <
            # JLM: right about here the code gets difficult to follow
            #      compared to the theory/description.
            #      change variables to match eqns, reuse intermediate
            #      variables less.

            # Compute contribution to Dunnian flow from PFR, if any
            if self._pref_flow_flag[hh]:
                # PRMSIV eqn 1-135 (? - value of topfr is not obvious)
                availh2o = self.pref_flow_stor[hh] + topfr
                dunnianflw_gvr = max(zero, availh2o - self.pref_flow_max[hh])

                if dunnianflw_gvr > zero:
                    # topfr = topfr - dunnianflw_gvr
                    # JLM: maybe using this variable too much?
                    # PRMSIV eqn 1-136
                    topfr = max(zero, topfr - dunnianflw_gvr)

                # <
                #
                self.pref_flow_in[hh] = self.pref_flow_infil[hh] + topfr
                self.pref_flow_stor[hh] = self.pref_flow_stor[hh] + topfr

                if self.pref_flow_stor[hh] > zero:
                    self.compute_interflow(
                        self.fastcoef_lin[hh],
                        self.fastcoef_sq[hh],
                        self.pref_flow_in[hh],
                        self.pref_flow_stor[hh],
                        prefflow,
                    )
                # <
            elif self.hru_type[hh] == HruType.LAND:
                dunnianflw_gvr = topfr  # ?? is this right

            # <
            self._gvr2pfr[hh] = topfr
            pervactet = zero

            # Compute actual evapotranspiration
            if self.soil_moist[hh] > zero:
                (
                    self.soil_moist[hh],
                    self.soil_rechr[hh],
                    avail_potet,
                    self.potet_rechr[hh],
                    self.potet_lower[hh],
                    pervactet,
                ) = self.compute_szactet(
                    self.transp_on[hh],
                    self.cov_type[hh],
                    self.soil_type[hh],
                    self.soil_moist_max[hh],
                    self.soil_rechr_max[hh],
                    self.snow_free[hh],
                    self.soil_moist[hh],
                    self.soil_rechr[hh],
                    avail_potet,
                    self.potet_rechr[hh],
                    self.potet_lower[hh],
                )

            # <
            self.hru_actet[hh] = (
                self.hru_actet[hh] + pervactet * self.hru_frac_perv[hh]
            )
            avail_potet = self.potet[hh] - self.hru_actet[hh]
            self.perv_actet[hh] = pervactet
            self.soil_lower[hh] = self.soil_moist[hh] - self.soil_rechr[hh]

            if self.hru_type[hh] == HruType.LAND.value:
                interflow = self.slow_flow[hh] + prefflow

                dunnianflw = dunnianflw_gvr + dunnianflw_pfr
                self.dunnian_flow[hh] = dunnianflw

                # Treat pref_flow as interflow
                self.ssres_flow[hh] = self.slow_flow[hh]

                if self._pref_flow_flag[hh]:
                    self.pref_flow[hh] = prefflow
                    self.ssres_flow[hh] = self.ssres_flow[hh] + prefflow

                # <
                # Treat dunnianflw as surface runoff to streams
                # WARNING: PAN This is modifying sroff from the srunoff module
                self.sroff[hh] = self.sroff[hh] + self.dunnian_flow[hh]
                self.ssres_stor[hh] = (
                    self.slow_stor[hh] + self.pref_flow_stor[hh]
                )

            else:
                # For swales
                availh2o = self.slow_stor[hh] - self.sat_threshold[hh]
                self.swale_actet[hh] = zero

                if availh2o > zero:
                    # If ponding, as storage > sat_threshold
                    unsatisfied_et = self.potet[hh] - self.hru_actet[hh]

                    if unsatisfied_et > zero:
                        availh2o = min(availh2o, unsatisfied_et)
                        self.swale_actet[hh] = availh2o
                        self.hru_actet[hh] = (
                            self.hru_actet[hh] + self.swale_actet[hh]
                        )
                        self.slow_stor[hh] = (
                            self.slow_stor[hh] - self.swale_actet[hh]
                        )
                    # <

                # <
                self.ssres_stor[hh] = self.slow_stor[hh]

            # <
            self.ssres_in[hh] = (
                self.soil_to_ssr[hh] + self.pref_flow_infil[hh] + gwin
            )
            self._grav_dunnian_flow[hh] = dunnianflw_gvr
            self.unused_potet[hh] = self.potet[hh] - self.hru_actet[hh]

        # <

        # refactor with np.where
        wh_lower_stor_max_gt_zero = np.where(self.soil_lower_max > zero)
        self.soil_lower_ratio[wh_lower_stor_max_gt_zero] = (
            self.soil_lower[wh_lower_stor_max_gt_zero]
            / self.soil_lower_max[wh_lower_stor_max_gt_zero]
        )
        # if self.control.current_time == np.datetime64("1979-03-18T00:00:00"):
        # asdf

        self.soil_moist_tot = (
            self.ssres_stor + self.soil_moist * self.hru_frac_perv
        )
        self.recharge = self.soil_to_gw + self.ssr_to_gw

        if self.control.config["dprst_flag"] == 1:
            self.recharge = self.recharge + self.dprst_seep_hru

        self.pref_flow_stor_change[:] = (
            self.pref_flow_stor - self.pref_flow_stor_prev
        )
        self.soil_lower_change[:] = self.soil_lower - self.soil_lower_prev
        self.soil_rechr_change[:] = self.soil_rechr - self.soil_rechr_prev
        self.slow_stor_change[:] = self.slow_stor - self.slow_stor_prev
        # Apparently the following are sums of the above and not actual
        # inddividual storage changes
        # self.soil_moist_change[:] = self.soil_moist - self.soil_moist_prev
        # self.ssres_stor_change[:] = self.ssres_stor - self.ssres_stor_prev

        self.soil_lower_change_hru[:] = (
            self.soil_lower_change * self.hru_frac_perv
        )
        self.soil_rechr_change_hru[:] = (
            self.soil_rechr_change * self.hru_frac_perv
        )
        self.perv_actet_hru[:] = self.perv_actet * self.hru_frac_perv

        self.ssres_flow_vol[:] = (
            self.ssres_flow * self.control.params.hru_in_to_cf
        )

        return

    @staticmethod
    def compute_soilmoist(
        soil2gw_flag,
        perv_frac,
        soil_moist_max,
        soil_rechr_max,
        soil2gw_max,
        infil,
        soil_moist,
        soil_rechr,
        soil_to_gw,
        soil_to_ssr,
    ) -> tuple:

        # JLM: i dont see any real advantage of using this function.

        # PRMSIV Step 4 (eqn 1-125)
        # prognostic
        soil_rechr = np.minimum(soil_rechr + infil, soil_rechr_max)

        # PRMSIV Step 5
        # soil_moist_max from previous time step or soil_moist_max has
        # changed for a restart simulation
        # prognostic
        excess = soil_moist + infil
        soil_moist = np.minimum(excess, soil_moist_max)  # PRMSIV eqn 1-126
        # JLM: not yet sure why eqn 1-127 is skipped here

        # PRMSIV Step 6
        # The following are are eqns 1-128 and 1-129
        # not prognostic
        excess = (excess - soil_moist_max) * perv_frac

        if excess > zero:
            if soil2gw_flag:
                # PRMSIV eqn 1-130
                soil_to_gw = np.minimum(soil2gw_max, excess)
                # PRMSIV eqn 1-131 (start)
                # not prognostic
                excess = excess - soil_to_gw
                # this "excess" is also known as gvr_maxin which is used in
                # calculations in later sections but it is local here.

            # <
            # JLM: why tampering with infil/cap_waterin and returning that?
            #      this is being set on self and renormalized.
            # are theyse just SANITY CHECK?
            if excess > (infil * perv_frac):
                # Probably dynamic
                infil = zero

            else:
                # ???? what if infil<0 ??? might happen with dynamic and
                #      small values,
                # ???? maybe ABS < NEARZERO = zero
                infil = infil - (excess / perv_frac)

            # <
            # PRMSIV eqn 1-131 (finish)
            soil_to_ssr = np.maximum(zero, excess)
            # THis is also known as gvr_maxin

        # <
        return (
            infil,
            soil_moist,
            soil_rechr,
            soil_to_gw,
            soil_to_ssr,
        )

    @staticmethod
    def compute_interflow(
        coef_lin, coef_sq, ssres_in, storage, inter_flow
    ) -> tuple:
        # inter_flow is in inches for the timestep
        # JLM: this is being way too clever. I am not sure there's a need for
        # this function. The theory shows 3 oneline equations, this is
        # a wild departure which makes me assume that there are things going
        # on with the coefficients that I dont know about.

        if (coef_lin <= zero) and (ssres_in <= zero):
            # print("interflow: 1")
            c1 = coef_sq * storage
            inter_flow = storage * (c1 / (one + c1))

        elif (coef_lin > zero) and (coef_sq <= zero):
            # print("interflow: 2")
            c2 = one - np.exp(-coef_lin)
            inter_flow = ssres_in * (one - c2 / coef_lin) + storage * c2

        elif coef_sq > zero:
            # print("interflow: 3")
            c3 = np.sqrt(coef_lin**2.0 + 4.0 * coef_sq * ssres_in)
            sos = storage - ((c3 - coef_lin) / (2.0 * coef_sq))
            if c3 == zero:
                msg = (
                    "ERROR, in compute_interflow sos=0, "
                    "please contact code developers"
                )
                raise RuntimeError(msg)

            # <
            c1 = coef_sq * sos / c3
            c2 = one - np.exp(-c3)

            if one + c1 * c2 > zero:
                inter_flow = ssres_in + (
                    (sos * (one + c1) * c2) / (one + c1 * c2)
                )
            else:
                inter_flow = ssres_in

            # <
        else:
            inter_flow = zero

        # <

        # TODO: Check if inter_flow < zero; if so, throw warning, reset to zero
        inter_flow = min(inter_flow, storage)
        # ifinter_flow > storage:
        #   inter_flow = storage
        # # <
        storage = storage - inter_flow
        return storage, inter_flow

    @staticmethod
    def compute_gwflow(ssr2gw_rate, ssr2gw_exp, slow_stor) -> tuple:
        # Compute flow to groundwater
        ssr_to_gw = max(0.0, ssr2gw_rate * slow_stor**ssr2gw_exp)
        ssr_to_gw = min(ssr_to_gw, slow_stor)
        slow_stor = slow_stor - ssr_to_gw
        return ssr_to_gw, slow_stor

    @staticmethod
    def compute_szactet(
        transp_on,
        cov_type,
        soil_type,
        soil_moist_max,
        soil_rechr_max,
        snow_free,
        soil_moist,
        soil_rechr,
        avail_potet,
        potet_rechr,
        potet_lower,
    ) -> tuple:

        # Determine type of evapotranspiration.
        #   et_type=2    - evaporation only
        #   et_type=3    - transpiration plus evaporation
        #   et_type=1    - default... JLM: which is ??

        if avail_potet < epsilon:
            et_type = ETType.ET_DEFAULT  # 1
            avail_potet = zero

        elif not transp_on:

            if snow_free < 0.01:
                et_type = ETType.ET_DEFAULT  # 1

            else:
                et_type = ETType.EVAP_ONLY  # 2

            # <
        elif cov_type > 0:
            et_type = ETType.EVAP_PLUS_TRANSP  # 3

        elif snow_free < 0.01:
            et_type = ETType.ET_DEFAULT  # 1

        else:
            et_type = ETType.EVAP_ONLY  # 2

        # <

        # ifet_type > 1:
        if et_type in [ETType.EVAP_ONLY, ETType.EVAP_PLUS_TRANSP]:
            pcts = soil_moist / soil_moist_max
            pctr = soil_rechr / soil_rechr_max
            potet_lower = avail_potet
            potet_rechr = avail_potet

            if soil_type == SoilType.SAND.value:

                if pcts < 0.25:
                    potet_lower = 0.5 * pcts * avail_potet

                # <
                if pctr < 0.25:
                    potet_rechr = 0.5 * pctr * avail_potet

                # <
            elif soil_type == SoilType.LOAM.value:

                if pcts < 0.5:
                    potet_lower = pcts * avail_potet

                # <
                if pctr < 0.5:
                    potet_rechr = pctr * avail_potet

                # <
            elif soil_type == SoilType.CLAY.value:

                if (pcts < TWOTHIRDS) and (pcts > ONETHIRD):
                    potet_lower = pcts * avail_potet

                elif pcts <= ONETHIRD:
                    potet_lower = 0.5 * pcts * avail_potet

                # <
                if (pctr < TWOTHIRDS) and (pctr > ONETHIRD):
                    potet_rechr = pctr * avail_potet

                elif pctr <= ONETHIRD:
                    potet_rechr = 0.5 * pctr * avail_potet

                # <
            else:
                pass  ## JLM?

            # <
            # ****** Soil moisture accounting
            if et_type == ETType.EVAP_ONLY:
                potet_rechr = potet_rechr * snow_free

            # <
            if potet_rechr > soil_rechr:
                potet_rechr = soil_rechr
                soil_rechr = zero

            else:
                soil_rechr = soil_rechr - potet_rechr

            # <
            if (et_type == ETType.EVAP_ONLY) or (potet_rechr >= potet_lower):

                if potet_rechr > soil_moist:
                    potet_rechr = soil_moist
                    soil_moist = zero

                else:
                    soil_moist = soil_moist - potet_rechr

                # <
                et = potet_rechr

            elif potet_lower > soil_moist:
                et = soil_moist
                soil_moist = zero

            else:
                soil_moist = soil_moist - potet_lower
                et = potet_lower

            # <
            if soil_rechr > soil_moist:
                soil_rechr = soil_moist

            # <
        else:
            et = zero

        # <
        return (
            soil_moist,
            soil_rechr,
            avail_potet,
            potet_rechr,
            potet_lower,
            et,
        )
