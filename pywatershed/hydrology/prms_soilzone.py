from typing import Literal
from warnings import warn

import numpy as np
from numba import prange

from ..base.adapter import adaptable
from ..base.conservative_process import ConservativeProcess
from ..base.control import Control
from ..constants import (
    ETType,
    HruType,
    SoilType,
    nan,
    nearzero,
    numba_num_threads,
    one,
    zero,
)
from ..parameters import Parameters

ONETHIRD = 1 / 3
TWOTHIRDS = 2 / 3


class PRMSSoilzone(ConservativeProcess):
    """PRMS soil zone.

    Implementation based on PRMS 5.2.1 with theoretical documentation given in
    the PRMS-IV documentation:

    `Markstrom, S. L., Regan, R. S., Hay, L. E., Viger, R. J., Webb, R. M.,
    Payn, R. A., & LaFontaine, J. H. (2015). PRMS-IV, the
    precipitation-runoff modeling system, version 4. US Geological Survey
    Techniques and Methods, 6, B7.
    <https://pubs.usgs.gov/tm/6b7/pdf/tm6-b7.pdf>`__

    Args:
        control: a Control object
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters
        dprst_evap_hru: Evaporation from surface-depression storage for each
            HRU
        dprst_seep_hru: Seepage from surface-depression storage to associated
            GWR for each HRU
        hru_impervevap: HRU area-weighted average evaporation from impervious
            area for each HRU
        hru_intcpevap: HRU area-weighted average evaporation from the
            canopy for each HRU
        infil_hru: Infiltration to the capillary and preferential-flow
            reservoirs, depth on HRU area
        sroff: Surface runoff to the stream network for each HRU
        potet: Potential ET for each HRU
        transp_on: Flag indicating whether transpiration is occurring
            (0=no;1=yes)
        snow_evap: Evaporation and sublimation from snowpack on each HRU
        snowcov_area: Snow-covered area on each HRU prior to melt and
            sublimation unless snowpack
        budget_type: one of [None, "warn", "error"]
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        adjust_parameters: one of ["warn", "error", "no"]. Default is "warn",
            the code edits the parameters and issues a warning. If "error" is
            selected the the code issues warnings about all edited parameters
            before raising the error to give you information. If "no" is
            selected then no parameters are adjusted and there will be no
            warnings or errors.
        verbose: Print extra information or not?
    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        dprst_evap_hru: adaptable,
        dprst_seep_hru: adaptable,
        hru_impervevap: adaptable,
        hru_intcpevap: adaptable,
        infil_hru: adaptable,  # in /pywatershed/analysis/budget_soilzone.py
        sroff: adaptable,
        potet: adaptable,
        transp_on: adaptable,
        snow_evap: adaptable,
        snowcov_area: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        calc_method: Literal["numba", "numpy"] = None,
        adjust_parameters: Literal["warn", "error", "no"] = "warn",
        verbose: bool = None,
    ) -> "PRMSSoilzone":
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "PRMSSoilzone"

        self._set_inputs(locals())
        self._set_options(locals())

        # This uses options
        self._initialize_soilzone_data()

        self._set_budget()
        self._init_calc_method()

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "dprst_frac",
            "cov_type",
            "fastcoef_lin",
            "fastcoef_sq",
            "hru_area",
            "hru_in_to_cf",
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
        return (
            "dprst_evap_hru",  # JLM ?? needs this stuff to calculate evap?
            "dprst_seep_hru",
            "hru_impervevap",  # JLM ??
            "hru_intcpevap",  # JLM ???
            "infil_hru",
            "sroff",
            "potet",
            # hru_ppt => model_precip%hru_ppt, & # JLM ??
            "transp_on",
            "snow_evap",
            "snowcov_area",
        )

    @staticmethod
    def get_restart_variables() -> tuple:
        return ()

    @staticmethod
    def get_init_values() -> dict:
        return {
            "cap_infil_tot": zero,
            "cap_waterin": zero,
            "dunnian_flow": zero,
            "hru_actet": zero,
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
            "soil_lower": nan,  # completely set later
            "soil_lower_change": nan,
            "soil_lower_change_hru": nan,
            "soil_lower_prev": nan,
            "soil_lower_ratio": zero,
            "soil_lower_max": nan,  # completely set later
            "soil_moist": nan,  # sm_climateflow
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
                "slow_flow",
                # "pref_flow",
            ],
            "storage_changes": [
                "soil_rechr_change_hru",
                "soil_lower_change_hru",
                "slow_stor_change",
                # "pref_flow_stor_change",
            ],
        }

    def _set_initial_conditions(self):
        # this is called in the super before options are set on self
        pass

    def _initialize_soilzone_data(self):
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

        self._snow_free = one - self.snowcov_area
        # edit a param
        wh_inactive_or_lake = np.where(
            (self.hru_type == HruType.INACTIVE.value)
            | (self.hru_type == HruType.LAKE.value)
        )
        self._sat_threshold = self.sat_threshold.copy()
        self._sat_threshold[wh_inactive_or_lake] = zero
        # edit a param
        wh_not_land = np.where(self.hru_type != HruType.LAND.value)
        self._pref_flow_den = self.pref_flow_den.copy()
        self._pref_flow_den[wh_not_land] = zero

        # variables
        if True:
            # For now there is no restart capability. we'll use the following
            # line when there is
            # if self.control.options["restart"] in [0, 2, 5]:

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

        # ssres_stor
        if True:
            # For now there is no restart capability. we'll use the following
            # line when there is
            # if self.control.options["restart"] in [0, 2, 5]:

            self.ssres_stor = self.ssstor_init_frac * self._sat_threshold
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

        throw_error = False
        mask = self.soil_moist_max < 1.0e-5
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_moist_max < 1.0e-5, set to 1.0e-5 at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_moist_max = np.where(mask, 1.0e-5, self.soil_moist_max)

        mask = self.soil_rechr_max < 1.0e-5
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_rechr_max < 1.0e-5, set to 1.0e-5 at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_rechr_max = np.where(mask, 1.0e-5, self.soil_rechr_max)

        mask = self.soil_rechr_max > self.soil_moist_max
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_rechr_max > soil_moist_max, "
                "soil_rechr_max set to soil_moist_max at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_rechr_max = np.where(
                mask,
                self.soil_moist_max,
                self.soil_rechr_max,
            )

        mask = self.soil_rechr > self.soil_rechr_max
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_rechr_init > soil_rechr_max, "
                "setting soil_rechr_init to soil_rechr_max at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_rechr = np.where(
                mask,
                self.soil_rechr_max,
                self.soil_rechr,
            )

        mask = self.soil_moist > self.soil_moist_max
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_moist_init > soil_moist_max, "
                "setting soil_moist to soil_moist max at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_moist = np.where(
                mask,
                self.soil_moist_max,
                self.soil_moist,
            )

        mask = self.soil_rechr > self.soil_moist
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "soil_rechr > soil_moist, "
                "setting soil_rechr to soil_moist at indices: "
                f"{np.where(mask)[0]}"
            )
            warn(msg, UserWarning)
            self.soil_rechr = np.where(mask, self.soil_moist, self.soil_rechr)

        mask = self.ssres_stor > self._sat_threshold
        if mask.any() and self._adjust_parameters != "no":
            if self._adjust_parameters == "error":
                throw_error = True
            msg = (
                "ssres_stor > _sat_threshold, "
                "setting ssres_stor to _sat_threshold at indices: "
                f"{np.where(mask)[0]}"
            )
            self.ssres_stor = np.where(
                mask,
                self._sat_threshold,
                self.ssres_stor,
            )

        if throw_error:
            raise ValueError(
                "Some parameter values were edited and an error was requested."
                " See warnings for additional details."
            )

        # <
        # need to set on swale_limit self? move to variables?
        self._swale_limit = np.full(self.nhru, zero, "float64")
        wh_swale = np.where(self.hru_type == HruType.SWALE.value)
        self._swale_limit[wh_swale] = 3.0 * self._sat_threshold[wh_swale]

        self.pref_flow_thrsh[wh_swale] = self._sat_threshold[wh_swale]
        wh_land = np.where(self.hru_type == HruType.LAND.value)
        self.pref_flow_thrsh[wh_land] = self._sat_threshold[wh_land] * (
            one - self._pref_flow_den[wh_land]
        )
        self.pref_flow_max[wh_land] = (
            self._sat_threshold[wh_land] - self.pref_flow_thrsh[wh_land]
        )

        # Need to set pref_flow_flag on self? or add to variables?
        wh_land_and_prf_den = np.where(
            (self.hru_type == HruType.LAND.value)
            & (self._pref_flow_den > zero)
        )
        self._pref_flow_flag = np.full(self.nhru, False, dtype=int)
        self._pref_flow_flag[wh_land_and_prf_den] = True

        # can this one be combined with the restart read logic above?
        if True:
            # For now there is no restart capability. we'll use the following
            # line when there is
            # if self.control.options["restart"] in [0, 2, 5]:

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
            self._sat_threshold + self.soil_moist_max * self.hru_frac_perv
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

        return

    def _init_calc_method(self):
        if self._calc_method is None:
            self._calc_method = "numba"

        if self._calc_method.lower() not in ["numpy", "numba"]:
            msg = (
                f"Invalid calc_method={self._calc_method} for {self.name}. "
                f"Setting calc_method to 'numba' for {self.name}"
            )
            warn(msg, UserWarning)
            self._calc_method = "numba"

        if self._calc_method.lower() == "numba":
            import numba as nb

            numba_msg = f"{self.name} jit compiling with numba "
            nb_parallel = (numba_num_threads is not None) and (
                numba_num_threads > 1
            )
            if nb_parallel:
                numba_msg += f"and using {numba_num_threads} threads"
            print(numba_msg, flush=True)

            self._calculate_soilzone = nb.njit(
                fastmath=True, parallel=nb_parallel
            )(self._calculate_numpy)

            self._compute_gwflow = nb.njit(fastmath=True)(self._compute_gwflow)
            self._compute_interflow = nb.njit(fastmath=True)(
                self._compute_interflow
            )
            self._compute_soilmoist = nb.njit(fastmath=True)(
                self._compute_soilmoist
            )
            self._compute_szactet = nb.njit(fastmath=True)(
                self._compute_szactet
            )

        else:
            self._calculate_soilzone = self._calculate_numpy

        return

    def _advance_variables(self) -> None:
        self.pref_flow_stor_prev[:] = self.pref_flow_stor
        self.soil_rechr_prev[:] = self.soil_rechr
        self.soil_lower_prev[:] = self.soil_lower
        self.slow_stor_prev[:] = self.slow_stor
        return

    def _calculate(self, simulation_time):
        (
            self.soil_to_gw[:],
            self.soil_to_ssr[:],
            self.ssr_to_gw[:],
            self.slow_flow[:],
            self.ssres_flow[:],
            self.potet_rechr[:],
            self.potet_lower[:],
            self.cap_waterin[:],
            self.soil_moist[:],
            self.soil_rechr[:],
            self.hru_actet[:],
            self.cap_infil_tot[:],
            self.slow_stor[:],
            self.pref_flow_in[:],
            self.pref_flow_stor[:],
            self.perv_actet[:],
            self.soil_lower[:],
            self.dunnian_flow[:],
            self.perv_actet_hru[:],
            self.pref_flow[:],
            self.pref_flow_stor_change[:],
            self.recharge[:],
            self.slow_stor_change[:],
            self.soil_lower_change[:],
            self.soil_lower_change_hru[:],
            self.soil_lower_ratio[:],
            self.soil_moist_tot[:],
            self.soil_rechr_change[:],
            self.soil_rechr_change_hru[:],
            self.sroff[:],
            self.ssres_flow_vol[:],
            self.ssres_in[:],
            self.ssres_stor[:],
            self.swale_actet[:],
            self.unused_potet[:],
        ) = self._calculate_soilzone(
            _pref_flow_flag=self._pref_flow_flag,
            _snow_free=self._snow_free,
            _soil2gw_flag=self._soil2gw_flag,
            cap_infil_tot=self.cap_infil_tot,
            cap_waterin=self.cap_waterin,
            compute_gwflow=self._compute_gwflow,
            compute_interflow=self._compute_interflow,
            compute_soilmoist=self._compute_soilmoist,
            compute_szactet=self._compute_szactet,
            cov_type=self.cov_type,
            current_time=self.control.current_time,
            dprst_evap_hru=self.dprst_evap_hru,
            dprst_flag=True,  # self.control.options["dprst_flag"],
            dprst_seep_hru=self.dprst_seep_hru,
            dunnian_flow=self.dunnian_flow,
            fastcoef_lin=self.fastcoef_lin,
            fastcoef_sq=self.fastcoef_sq,
            hru_actet=self.hru_actet,
            hru_frac_perv=self.hru_frac_perv,
            hru_impervevap=self.hru_impervevap,
            hru_in_to_cf=self.hru_in_to_cf,
            hru_intcpevap=self.hru_intcpevap,
            hru_type=self.hru_type,
            infil_hru=self.infil_hru,
            nhru=self.nhru,
            perv_actet=self.perv_actet,
            perv_actet_hru=self.perv_actet_hru,
            potet=self.potet,
            potet_lower=self.potet_lower,
            potet_rechr=self.potet_rechr,
            pref_flow=self.pref_flow,
            pref_flow_in=self.pref_flow_in,
            pref_flow_infil=self.pref_flow_infil,
            pref_flow_max=self.pref_flow_max,
            pref_flow_stor=self.pref_flow_stor,
            pref_flow_stor_change=self.pref_flow_stor_change,
            pref_flow_stor_prev=self.pref_flow_stor_prev,
            pref_flow_thrsh=self.pref_flow_thrsh,
            recharge=self.recharge,
            sat_threshold=self._sat_threshold,
            slow_flow=self.slow_flow,
            slow_stor=self.slow_stor,
            slow_stor_change=self.slow_stor_change,
            slow_stor_prev=self.slow_stor_prev,
            slowcoef_lin=self.slowcoef_lin,
            slowcoef_sq=self.slowcoef_sq,
            snow_evap=self.snow_evap,
            snowcov_area=self.snowcov_area,
            soil2gw_max=self.soil2gw_max,
            soil_lower=self.soil_lower,
            soil_lower_change=self.soil_lower_change,
            soil_lower_change_hru=self.soil_lower_change_hru,
            soil_lower_max=self.soil_lower_max,
            soil_lower_prev=self.soil_lower_prev,
            soil_lower_ratio=self.soil_lower_ratio,
            soil_moist=self.soil_moist,
            soil_moist_max=self.soil_moist_max,
            soil_moist_tot=self.soil_moist_tot,
            soil_rechr=self.soil_rechr,
            soil_rechr_change=self.soil_rechr_change,
            soil_rechr_change_hru=self.soil_rechr_change_hru,
            soil_rechr_max=self.soil_rechr_max,
            soil_rechr_prev=self.soil_rechr_prev,
            soil_to_gw=self.soil_to_gw,
            soil_to_ssr=self.soil_to_ssr,
            soil_type=self.soil_type,
            sroff=self.sroff,
            ssr2gw_exp=self.ssr2gw_exp,
            ssr2gw_rate=self.ssr2gw_rate,
            ssr_to_gw=self.ssr_to_gw,
            ssres_flow=self.ssres_flow,
            ssres_flow_vol=self.ssres_flow_vol,
            ssres_in=self.ssres_in,
            ssres_stor=self.ssres_stor,
            swale_actet=self.swale_actet,
            transp_on=self.transp_on,
            unused_potet=self.unused_potet,
        )

        return

    @staticmethod
    def _calculate_numpy(
        _pref_flow_flag,
        _snow_free,
        _soil2gw_flag,
        cap_infil_tot,
        cap_waterin,
        compute_gwflow,
        compute_interflow,
        compute_soilmoist,
        compute_szactet,
        cov_type,
        current_time,
        dprst_evap_hru,
        dprst_flag,
        dprst_seep_hru,
        dunnian_flow,
        fastcoef_lin,
        fastcoef_sq,
        hru_actet,
        hru_frac_perv,
        hru_impervevap,
        hru_in_to_cf,
        hru_intcpevap,
        hru_type,
        infil_hru,
        nhru,
        potet,
        potet_lower,
        potet_rechr,
        slow_flow,
        slow_stor,
        soil_moist,
        soil_moist_max,
        soil_rechr,
        soil_to_gw,
        soil_to_ssr,
        soil_type,
        ssr_to_gw,
        ssres_flow,
        perv_actet,
        perv_actet_hru,
        pref_flow,
        pref_flow_in,
        pref_flow_infil,
        pref_flow_max,
        pref_flow_stor,
        pref_flow_stor_change,
        pref_flow_stor_prev,
        pref_flow_thrsh,
        recharge,
        sat_threshold,
        slow_stor_change,
        slow_stor_prev,
        slowcoef_lin,
        slowcoef_sq,
        snow_evap,
        snowcov_area,
        soil2gw_max,
        soil_lower,
        soil_lower_change,
        soil_lower_change_hru,
        soil_lower_max,
        soil_lower_prev,
        soil_lower_ratio,
        soil_moist_tot,
        soil_rechr_change,
        soil_rechr_change_hru,
        soil_rechr_max,
        soil_rechr_prev,
        sroff,
        ssr2gw_exp,
        ssr2gw_rate,
        ssres_flow_vol,
        ssres_in,
        ssres_stor,
        swale_actet,
        transp_on,
        unused_potet,
    ):
        """Calculate soil zone for a time step"""

        # JLM: not clear we need this / for GSFlow
        # if srunoff_updated_soil:
        #     soil_moist = soil_moist_change
        #     soil_rechr = soil_rechr_change
        # # <

        # <
        gwin = zero
        # update_potet = 0

        # diagnostic state resets
        soil_to_gw[:] = zero
        soil_to_ssr[:] = zero
        ssr_to_gw[:] = zero
        slow_flow[:] = zero
        ssres_flow[:] = zero
        potet_rechr[:] = zero
        potet_lower[:] = zero

        _snow_free[:] = one - snowcov_area

        # we dont track soil_moist_prev as it's not prognostic
        # soil_moist_prev = soil_rechr and soil_lower
        # soil_moist_prev[:] = soil_moist

        # JLM: ET calculations to be removed from soilzone.
        hru_actet = hru_impervevap + hru_intcpevap + snow_evap

        # This is obnoxious. i guess this should be an
        # optional input? should default to zero?
        if dprst_flag:
            hru_actet = hru_actet + dprst_evap_hru

        # <
        for hh in prange(nhru):
            dunnianflw = zero
            dunnianflw_pfr = zero
            dunnianflw_gvr = zero
            # interflow = zero

            # JLM: ET calculation to be removed from soilzone.
            avail_potet = np.maximum(zero, potet[hh] - hru_actet[hh])

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
            capwater_maxin = infil_hru[hh] / hru_frac_perv[hh]

            # Compute preferential flow and storage, and any dunnian flow
            prefflow = zero
            # if _pref_flow_flag[hh]:
            #     pref_flow_infil[hh] = zero
            #     pref_flow_maxin = zero

            #     if capwater_maxin > zero:
            #         # PRMSIV Step 1 - Partition infil between capillary and
            #         #                 preferential-flow (eqn 1-121)
            #         # pref_flow for whole HRU but capwater is pervious area
            #         # calculations on pervious area
            #         pref_flow_maxin = capwater_maxin * pref_flow_den[hh]

            #         # PRMSIV Step 3: no cascades and already normalized to
            #         #                pervious area. (eqn 1-124)
            #         capwater_maxin = capwater_maxin - pref_flow_maxin

            #         # renormalize pref_flow to whole HRU from pervious area
            #         pref_flow_maxin = pref_flow_maxin * hru_frac_perv[hh]

            #         # PRMSIV Step 2 - Compute PFR storage, excess to Dunnian
            #         #                 (eqns 1-122 and 1-123)
            #         # Compute contribution to preferential-flow reservoir
            #         # storage
            #         pref_flow_stor[hh] = (
            #             pref_flow_stor[hh] + pref_flow_maxin
            #         )
            #         dunnianflw_pfr = max(
            #             zero, pref_flow_stor[hh] - pref_flow_max[hh]
            #         )

            #         if dunnianflw_pfr > zero:
            #             pref_flow_stor[hh] = pref_flow_max[hh]

            #         # <
            #         pref_flow_infil[hh] = pref_flow_maxin - dunnianflw_pfr

            #     # <
            #     _pfr_dunnian_flow[hh] = dunnianflw_pfr

            # <
            # whole HRU
            cap_infil_tot[hh] = capwater_maxin * hru_frac_perv[hh]

            # ****** Add infiltration to soil and compute excess
            # Step 3 above
            cap_waterin[hh] = capwater_maxin

            # PRMSIV Steps 4, 5, 6 (see compute_soilmoist)
            if (capwater_maxin + soil_moist[hh]) > zero:
                # JLM: not sure why the function returns cap_waterin
                (
                    cap_waterin[hh],
                    soil_moist[hh],
                    soil_rechr[hh],
                    soil_to_gw[hh],
                    soil_to_ssr[hh],
                ) = compute_soilmoist(
                    _soil2gw_flag[hh],
                    hru_frac_perv[hh],
                    soil_moist_max[hh],
                    soil_rechr_max[hh],
                    soil2gw_max[hh],
                    cap_waterin[hh],
                    soil_moist[hh],
                    soil_rechr[hh],
                    soil_to_gw[hh],
                    soil_to_ssr[hh],
                )

                cap_waterin[hh] = cap_waterin[hh] * hru_frac_perv[hh]

            # <

            topfr = zero

            # soil_to_ssr also known as grv_maxin
            # grv_availh2o = grv_stor + soil_to_ssr
            availh2o = slow_stor[hh] + soil_to_ssr[hh]

            if hru_type[hh] == HruType.LAND.value:
                # PRMSIV Step 7
                # PRMSIV eqn 1-132? not necessary?
                # PRMSIV eqn 1-133:
                #     gvr_excess = max(
                #         0, slow_stor_old + gvr_maxin - pref_flow_thresh)
                topfr = max(zero, availh2o - pref_flow_thrsh[hh])

                # topfr = max(
                #    zero,
                #    slow_stor[hh]
                #    + soil_to_ssr[hh]
                #    - pref_flow_thrsh[hh],
                # )

                # PRMSIV eqn 1-134: ssres_in gvr_maxin - gvr_excess
                ssresin = soil_to_ssr[hh] - topfr
                # ssresin = soil_to_ssr[hh] - max(
                #     zero,
                #     slow_stor[hh]
                #     + soil_to_ssr[hh]
                #     - pref_flow_thrsh[hh],
                # )

                # JLM: what equation is this?
                slow_stor[hh] = max(zero, availh2o - topfr)
                # JLM: This is expansion is absurd
                # slow_stor[hh] = max(
                #     zero,
                #     slow_stor[hh]
                #     + soil_to_ssr[hh]
                #     - max(
                #         zero,
                #         slow_stor[hh]
                #         + soil_to_ssr[hh]
                #         - pref_flow_thrsh[hh],
                #     ),
                # )
                # if pref_flow_thrsh < availh2o, then
                # slow_stor = pref_flow_thrsh ???

                # JLM: the order of the PRMSIV theory is not maintained here
                # This skips over 8 to step 9

                # PRMSIV Step 9
                # Compute slow contribution to interflow, if any
                if slow_stor[hh] > zero:  # single precision?
                    (
                        slow_stor[hh],
                        slow_flow[hh],
                    ) = compute_interflow(
                        slowcoef_lin[hh],
                        slowcoef_sq[hh],
                        ssresin,
                        slow_stor[hh],
                        slow_flow[hh],
                    )

                # <
            elif hru_type[hh] == HruType.SWALE.value:
                slow_stor[hh] = availh2o

            # <
            if (slow_stor[hh] > zero) and (ssr2gw_rate[hh] > zero):
                (
                    ssr_to_gw[hh],
                    slow_stor[hh],
                ) = compute_gwflow(
                    ssr2gw_rate[hh],
                    ssr2gw_exp[hh],
                    slow_stor[hh],
                )

            # <
            # JLM: right about here the code gets difficult to follow
            #      compared to the theory/description.
            #      change variables to match eqns, reuse intermediate
            #      variables less.

            # Compute contribution to Dunnian flow from PFR, if any
            if _pref_flow_flag[hh]:
                # PRMSIV eqn 1-135 (? - value of topfr is not obvious)
                availh2o = pref_flow_stor[hh] + topfr
                dunnianflw_gvr = max(zero, availh2o - pref_flow_max[hh])

                if dunnianflw_gvr > zero:
                    # topfr = topfr - dunnianflw_gvr
                    # JLM: maybe using this variable too much?
                    # PRMSIV eqn 1-136
                    topfr = max(zero, topfr - dunnianflw_gvr)

                # <
                #
                pref_flow_in[hh] = pref_flow_infil[hh] + topfr
                pref_flow_stor[hh] = pref_flow_stor[hh] + topfr

                if pref_flow_stor[hh] > zero:
                    compute_interflow(
                        fastcoef_lin[hh],
                        fastcoef_sq[hh],
                        pref_flow_in[hh],
                        pref_flow_stor[hh],
                        prefflow,
                    )
                # <
            elif hru_type[hh] == HruType.LAND.value:
                dunnianflw_gvr = topfr  # ?? is this right

            # <
            perv_actet[hh] = zero

            # Compute actual evapotranspiration
            if soil_moist[hh] > zero:
                (
                    soil_moist[hh],
                    soil_rechr[hh],
                    avail_potet,
                    potet_rechr[hh],
                    potet_lower[hh],
                    perv_actet[hh],
                ) = compute_szactet(
                    transp_on[hh],
                    cov_type[hh],
                    soil_type[hh],
                    soil_moist_max[hh],
                    soil_rechr_max[hh],
                    _snow_free[hh],
                    soil_moist[hh],
                    soil_rechr[hh],
                    avail_potet,
                    potet_rechr[hh],
                    potet_lower[hh],
                )

            # <
            hru_actet[hh] = hru_actet[hh] + perv_actet[hh] * hru_frac_perv[hh]
            avail_potet = potet[hh] - hru_actet[hh]
            soil_lower[hh] = soil_moist[hh] - soil_rechr[hh]

            if hru_type[hh] == HruType.LAND.value:
                # interflow = slow_flow[hh] + prefflow

                dunnianflw = dunnianflw_gvr + dunnianflw_pfr
                dunnian_flow[hh] = dunnianflw

                # Treat pref_flow as interflow
                ssres_flow[hh] = slow_flow[hh]

                if _pref_flow_flag[hh]:
                    pref_flow[hh] = prefflow
                    ssres_flow[hh] = ssres_flow[hh] + prefflow

                # <
                # Treat dunnianflw as surface runoff to streams
                # WARNING: PAN This is modifying sroff from the srunoff module
                sroff[hh] = sroff[hh] + dunnian_flow[hh]
                ssres_stor[hh] = slow_stor[hh] + pref_flow_stor[hh]

            else:
                # For swales
                availh2o = slow_stor[hh] - sat_threshold[hh]
                swale_actet[hh] = zero

                if availh2o > zero:
                    # If ponding, as storage > sat_threshold
                    unsatisfied_et = potet[hh] - hru_actet[hh]

                    if unsatisfied_et > zero:
                        availh2o = min(availh2o, unsatisfied_et)
                        swale_actet[hh] = availh2o
                        hru_actet[hh] = hru_actet[hh] + swale_actet[hh]
                        slow_stor[hh] = slow_stor[hh] - swale_actet[hh]
                    # <

                # <
                ssres_stor[hh] = slow_stor[hh]

            # <
            ssres_in[hh] = soil_to_ssr[hh] + pref_flow_infil[hh] + gwin
            unused_potet[hh] = potet[hh] - hru_actet[hh]

        # <
        # Could mo move the remaining code to _calculate
        # refactor with np.where
        wh_lower_stor_max_gt_zero = np.where(soil_lower_max > zero)
        soil_lower_ratio[wh_lower_stor_max_gt_zero] = (
            soil_lower[wh_lower_stor_max_gt_zero]
            / soil_lower_max[wh_lower_stor_max_gt_zero]
        )
        # if current_time == np.datetime64("1979-03-18T00:00:00"):
        # asdf

        soil_moist_tot = ssres_stor + soil_moist * hru_frac_perv
        recharge = soil_to_gw + ssr_to_gw

        if dprst_flag:
            recharge = recharge + dprst_seep_hru

        pref_flow_stor_change[:] = pref_flow_stor - pref_flow_stor_prev
        soil_lower_change[:] = soil_lower - soil_lower_prev
        soil_rechr_change[:] = soil_rechr - soil_rechr_prev
        slow_stor_change[:] = slow_stor - slow_stor_prev
        # Apparently the following are sums of the above and not actual
        # inddividual storage changes
        # soil_moist_change[:] = soil_moist - soil_moist_prev
        # ssres_stor_change[:] = ssres_stor - ssres_stor_prev

        soil_lower_change_hru[:] = soil_lower_change * hru_frac_perv
        soil_rechr_change_hru[:] = soil_rechr_change * hru_frac_perv
        perv_actet_hru[:] = perv_actet * hru_frac_perv

        ssres_flow_vol[:] = ssres_flow * hru_in_to_cf

        return (
            soil_to_gw,
            soil_to_ssr,
            ssr_to_gw,
            slow_flow,
            ssres_flow,
            potet_rechr,
            potet_lower,
            cap_waterin,
            soil_moist,
            soil_rechr,
            hru_actet,
            cap_infil_tot,
            slow_stor,
            pref_flow_in,
            pref_flow_stor,
            perv_actet,
            soil_lower,
            dunnian_flow,
            perv_actet_hru,
            pref_flow,
            pref_flow_stor_change,
            recharge,
            slow_stor_change,
            soil_lower_change,
            soil_lower_change_hru,
            soil_lower_ratio,
            soil_moist_tot,
            soil_rechr_change,
            soil_rechr_change_hru,
            sroff,
            ssres_flow_vol,
            ssres_in,
            ssres_stor,
            swale_actet,
            unused_potet,
        )

    @staticmethod
    def _compute_soilmoist(
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
    def _compute_interflow(
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
    def _compute_gwflow(ssr2gw_rate, ssr2gw_exp, slow_stor) -> tuple:
        # Compute flow to groundwater
        ssr_to_gw = max(0.0, ssr2gw_rate * slow_stor**ssr2gw_exp)
        ssr_to_gw = min(ssr_to_gw, slow_stor)
        slow_stor = slow_stor - ssr_to_gw
        return ssr_to_gw, slow_stor

    @staticmethod
    def _compute_szactet(
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

        if avail_potet < nearzero:
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
                pass  # JLM

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
            et,  # -> perv_actet
        )
