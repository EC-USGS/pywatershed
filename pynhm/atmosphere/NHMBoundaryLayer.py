import pathlib as pl
from pprint import pprint
import warnings
from itertools import chain
from typing import Union

import netCDF4 as nc4

import numpy as np

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import inch2cm, nan, one, zero
from .solar_constants import solf

from ..utils.parameters import PrmsParameters
from ..utils.time_utils import datetime_month

from pynhm.base.storageUnit import StorageUnit
from .PRMSSolarGeometry import PRMSSolarGeometry

from pynhm.utils.netcdf_utils import NetCdfRead


fileish = Union[str, pl.Path]


# may not use this if they cant be called with jit
# if it can, put it in a common utility.
def tile_space_to_time(arr: np.ndarray, n_time) -> np.ndarray:
    return np.tile(arr, (n_time, 1))


def tile_time_to_space(arr: np.ndarray, n_space) -> np.ndarray:
    return np.transpose(np.tile(arr, (n_space, 1)))


class PRMSBoundaryLayer(StorageUnit):
    def __init__(
        self,
        control: Control,
        params: PrmsParameters,
        prcp: fileish,
        tmax: fileish,
        tmin: fileish,
        budget_type: str = None,
        verbose: bool = False,
        netcdf_output_dir: fileish = None,
    ):
        """PRMS atmospheric boundary layer models

        This representation uses precipitation and temperature inputs.
        Relative humidity could be added as well."""

        self.name = "PRMSBoundaryLayer"
        super().__init__(
            control=control,
            params=params,
            verbose=verbose,
            subclass_name=self.name,
        )

        # Override self.set_inputs(locals())
        # Get all data at all times: do all forcings up front
        # There will be no inputs to advance
        self._input_variables_dict = {}
        self._datetime = None
        for input in self.get_inputs():
            nc_data = NetCdfRead(locals()[input])
            self[input] = nc_data.dataset[input][:].data
            # Get the datetimes or check against the first
            if self._datetime is None:
                self._datetime = nc_data._datetime
            else:
                assert (nc_data._datetime == self._datetime).all()

        start_time_ind = np.where(self._datetime == self.control._start_time)
        start_time_ind = start_time_ind[0]
        if not len(start_time_ind):
            msg = "Control start_time is not in the input data datetime"
            raise ValueError(msg)
        self._init_time_ind = start_time_ind[0]

        # atm/boundary layer has a solar geometry
        self.solar_geom = PRMSSolarGeometry(control, params)

        # Solve all variables for all time
        self._month_ind_12 = datetime_month(self._datetime) - 1  # (time)
        self._month_ind_1 = np.zeros(self._datetime.shape, dtype=int)  # (time)

        self.adjust_temperature()
        self.adjust_precip()
        self.calculate_sw_rad_degree_day()
        self.calculate_potential_et_jh()

        # Budget is not for all time
        # self.set_budget(budget_type)
        # This is a nonstandard budget.
        self.budget = None
        # self.budget = {
        #     "inputs": {"potet": self.potet},
        #     "outputs": {"hru_actet": self.hru_actet},
        #     "storage_changes": {"available_potet": self.available_potet},
        # }

        if netcdf_output_dir:
            # write out all variables
            # check that netcdf_output_dir is a directory
            pass

        return

    @staticmethod
    def get_inputs() -> tuple:
        return ("prcp", "tmax", "tmin")

    @staticmethod
    def get_variables() -> tuple:
        return (
            "tmaxf",
            "tminf",
            "prmx",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "swrad",
            "potet",
            # "hru_actet",
            # "available_potet",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "tmaxf": nan,
            "tminf": nan,
            "prmx": nan,
            "hru_ppt": nan,
            "hru_rain": nan,
            "hru_snow": nan,
            "swrad": nan,
            "potet": nan,
            # "hru_actet": zero,
            # "available_potet": nan,
        }

    @staticmethod
    def get_parameters():
        return (
            "radadj_intcp",
            "radadj_slope",
            "tmax_index",
            "dday_slope",
            "dday_intcp",
            "radmax",
            "ppt_rad_adj",
            "tmax_allsnow",
            "tmax_allrain_offset",
            "hru_slope",
            "radj_sppt",
            "radj_wppt",
            "hru_lat",
            "hru_area",
            "nhru",
            "hru_aspect",
            "jh_coef",
            "jh_coef_hru",
            "tmax_cbh_adj",
            "tmin_cbh_adj",
            "tmax_allsnow",
            "tmax_allrain_offset",
            "snow_cbh_adj",
            "rain_cbh_adj",
            "adjmix_rain",
        )

    def set_initial_conditions(self):
        return

    def _advance_variables(self):
        """Advance the PRMSBoundaryLayer

        Returns:
            None

        """
        return

    def _calculate(self, time_length):
        """Calculate PRMSBoundaryLayer"

        Sets the current value to the precalculated values.

        Args:
            time_length: time step length

        Returns:
            None

        """
        current_ind = self._init_time_ind + self._itime_step
        for var in self.variables:
            self[var][:] = self[f"_{var}"][current_ind, :]
        return

    def adjust_temperature(self):
        """Input temperature adjustments using calibrated parameters."""

        # throw an error if these have different shapes
        if self["tmax_cbh_adj"].shape != self["tmin_cbh_adj"].shape:
            msg = (
                "Not implemented: tmin/tmax cbh adj parameters "
                "with different shapes"
            )
            raise NotImplementedError(msg)

        if self["tmax_cbh_adj"].shape[0] == 12:
            month_ind = self._month_ind_12
        elif self["tmax_cbh_adj"].shape[0] == 1:
            month_ind = self._month_ind_1
        else:
            msg = (
                "Unexpected month dimension for cbh "
                "temperature adjustment params"
            )
            raise ValueError(msg)

        # (time, space) dimensions on these variables
        self._tmaxf = np.zeros(self.tmax.shape, dtype=self.tmax.dtype)
        self._tminf = np.zeros(self.tmin.shape, dtype=self.tmin.dtype)
        self._tmaxf = self["tmax"] + self["tmax_cbh_adj"][month_ind]
        self._tminf = self["tmin"] + self["tmin_cbh_adj"][month_ind]

        return

    def adjust_precip(self):
        """Input precipitation adjustments using calibrated parameters.

        Snow/rain partitioning of total precip depends on adjusted temperature
        in addition to depending on additonal parameters.
        """

        # throw an error shapes are inconsistent
        shape_list = np.array(
            [
                self.tmax_allsnow.shape[0],
                self.tmax_allrain_offset.shape[0],
                self.snow_cbh_adj.shape[0],
                self.rain_cbh_adj.shape[0],
                self.adjmix_rain.shape[0],
            ]
        )
        if not (shape_list == 12).all():
            msg = (
                "Not implemented: cbh precip adjustment parameters with "
                " different shapes"
            )
            raise NotImplementedError(msg)

        tmax_allrain = self.tmax_allsnow + self.tmax_allrain_offset

        if self.tmax_allsnow.shape[0] == 12:
            month_ind = self._month_ind_12
        elif self.tmax_allsnow.shape[0] == 1:
            month_ind = self._month_ind_1
        else:
            msg = "Unexpected month dimension for cbh precip adjustment params"
            raise ValueError(msg)

        # (time, space) dimensions
        self._hru_ppt = np.zeros(self.prcp.shape, dtype=self.prcp.dtype)
        self._hru_rain = np.zeros(self.prcp.shape, dtype=self.prcp.dtype)
        self._hru_snow = np.zeros(self.prcp.shape, dtype=self.prcp.dtype)

        self._prmx = np.zeros(self.prcp.shape, dtype=self.prcp.dtype)
        # Order MATTERS in calculating the prmx mask
        # The logic in PRMS is if(all_snow),elif(all_rain),else(mixed)
        # so we set the mask in the reverse order
        # Calculate the mix everywhere, then set the precip/rain/snow amounts from the conditions.
        tdiff = self._tmaxf - self._tminf
        self._prmx = (
            (self._tmaxf - self.tmax_allsnow[month_ind]) / tdiff
        ) * self.adjmix_rain[month_ind]
        del tdiff

        wh_all_snow = np.where(self._tmaxf <= self.tmax_allsnow[month_ind])
        wh_all_rain = np.where(
            np.logical_or(
                self._tminf > self.tmax_allsnow[month_ind],
                self._tmaxf >= tmax_allrain[month_ind],
            )
        )
        self._prmx[wh_all_rain] = one
        self._prmx[wh_all_snow] = zero

        # Recalculate/redefine these now based on prmx instead of the temperature logic
        wh_all_snow = np.where(self._prmx <= zero)
        wh_all_rain = np.where(self._prmx >= one)

        # Mixed case (everywhere, to be overwritten by the all-snow/rain-fall cases)
        self._hru_ppt = self.prcp * self.snow_cbh_adj[month_ind]
        self._hru_rain = self._prmx * self._hru_ppt
        self._hru_snow = self._hru_ppt - self._hru_rain

        # All precip is snow case
        # The condition to be used later:
        self._hru_ppt[wh_all_snow] = (
            self.prcp * self.snow_cbh_adj[month_ind]
        )[wh_all_snow]
        self._hru_snow[wh_all_snow] = self._hru_ppt[wh_all_snow]
        self._hru_rain[wh_all_snow] = zero

        # All precip is rain case
        # The condition to be used later:
        self._hru_ppt[wh_all_rain] = (
            self.prcp * self.rain_cbh_adj[month_ind]
        )[wh_all_rain]
        self._hru_rain[wh_all_rain] = self._hru_ppt[wh_all_rain]
        self._hru_snow[wh_all_rain] = zero
        return

    def calculate_sw_rad_degree_day(self) -> None:
        """Calculate shortwave radiation using the degree day method."""

        solar_params = {}
        solar_param_names = (
            "radadj_intcp",
            "radadj_slope",
            "tmax_index",
            "dday_slope",
            "dday_intcp",
            "radmax",
            "ppt_rad_adj",
            "tmax_allsnow",
            "tmax_allrain_offset",
            "hru_slope",
            "radj_sppt",
            "radj_wppt",
            "hru_lat",
            "hru_area",
        )
        for name in solar_param_names:
            solar_params[name] = self.solar_geom[name]

        self._swrad = self._ddsolrad_run(
            dates=self._datetime,
            tmax_hru=self._tmaxf,
            hru_ppt=self._hru_ppt,
            soltab_potsw=self.solar_geom._soltab_potsw,
            **solar_params,
        )

        return

    # @jit
    @staticmethod
    def _ddsolrad_run(
        dates: np.ndarray,  # [n_time]
        tmax_hru: np.ndarray,  # [n_time, n_hru]
        hru_ppt: np.ndarray,  # [n_time, n_hru]
        soltab_potsw: np.ndarray,  # [n_time, n_hru] ??
        radadj_intcp,  # param [12, n_hru]
        radadj_slope,  # "    "
        tmax_index,
        dday_slope,
        dday_intcp,
        radmax,
        ppt_rad_adj,
        tmax_allsnow,
        tmax_allrain_offset,
        hru_slope,  # pram [n_hru]
        radj_sppt,
        radj_wppt,
        hru_lat,
        hru_area,
    ) -> np.ndarray:  # [n_time, n_hru]

        # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation_degday.f90
        n_time, n_hru = tmax_hru.shape

        # Transforms of time
        doy = (dates - dates.astype("datetime64[Y]")).astype(
            "timedelta64[h]"
        ).astype(int) / 24 + 1
        doy = tile_time_to_space(doy, n_hru).astype(int)
        month = (dates.astype("datetime64[M]").astype(int) % 12) + 1
        month = tile_time_to_space(month, n_hru)

        # is_summer
        # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/misc/sm_prms_time.f90#L114
        is_summer = (doy >= 79) & (doy <= 265)
        northern_hemisphere = hru_lat > zero
        if not any(northern_hemisphere):
            northern_hemisphere = tile_space_to_time(
                northern_hemisphere, n_time
            )
            msg = "Implementation not checked"
            raise NotImplementedError(msg)
            is_summer = is_summer & northern_hemisphere

        hru_cossl = tile_space_to_time(np.cos(np.arctan(hru_slope)), n_time)

        def monthly_to_daily(arr):
            return arr[month[:, 0] - 1, :]

        def doy_to_daily(arr):
            return arr[doy[:, 0] - 1, :]

        # This is the start of the loop over hru & time
        # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_solar_radiation_degday.f90#L103

        # dday
        dday_slope_day = monthly_to_daily(dday_slope)
        dday_intcp_day = monthly_to_daily(dday_intcp)
        dday = (dday_slope_day * tmax_hru) + dday_intcp_day + one
        dday[np.where(dday < one)] = one
        del dday_slope_day, dday_intcp_day

        # radadj
        radmax_day = monthly_to_daily(radmax)
        radadj = monthly_to_daily(radmax)  # the else condition
        wh_dday_lt_26 = np.where(dday < 26.0)
        if len(wh_dday_lt_26[0]):
            kp = dday.astype(int)[wh_dday_lt_26]  # dddayi = float(kp)
            radadj[wh_dday_lt_26] = solf[kp - 1] + (
                (solf[kp] - solf[kp - 1]) * (dday[wh_dday_lt_26] - kp)
            )
            wh_radadj_gt_max = np.where(radadj > radmax_day)
            if len(wh_radadj_gt_max[0]):
                radadj[wh_radadj_gt_max] = radmax_day[wh_radadj_gt_max]
            del kp, wh_radadj_gt_max
        del radmax_day, wh_dday_lt_26

        # pptadj
        pptadj = np.ones((n_time, n_hru))
        radadj_intcp_day = monthly_to_daily(radadj_intcp)
        radadj_slope_day = monthly_to_daily(radadj_slope)
        tmax_index_day = monthly_to_daily(tmax_index)

        ppt_rad_adj_day = monthly_to_daily(ppt_rad_adj)

        # * if
        # This is the outer if. All of the following conditions work on this
        cond_ppt_gt_rad_adj = hru_ppt > ppt_rad_adj_day
        cond_if = cond_ppt_gt_rad_adj
        wh_if = np.where(cond_if)
        if len(wh_if[0]):

            # * if.else
            pptadj[wh_if] = (
                radadj_intcp_day
                + radadj_slope_day * (tmax_hru - tmax_index_day)
            )[wh_if]

            # * if.else.if
            cond_ppt_adj_gt_one = pptadj > one
            cond_if_else_if = cond_if & cond_ppt_adj_gt_one
            wh_if_else_if = np.where(cond_if_else_if)
            if len(wh_if_else_if[0]):
                pptadj[wh_if_else_if] = one
            del radadj_intcp_day, radadj_slope_day, tmax_index_day

            # * if.if
            tmax_index_day = monthly_to_daily(tmax_index)
            cond_tmax_lt_index = tmax_hru < tmax_index_day
            cond_if_if = cond_if & cond_tmax_lt_index
            wh_if_if = np.where(cond_if_if)
            if len(wh_if_if[0]):

                # The logic equiv to but changed from the original
                radj_sppt_day = tile_space_to_time(radj_sppt, n_time)
                pptadj[wh_if_if] = radj_sppt_day[wh_if_if]

                # if.if.else: negate the if.if.if
                radj_wppt_day = tile_space_to_time(radj_wppt, n_time)
                tmax_allrain_day = monthly_to_daily(
                    tmax_allrain_offset + tmax_allsnow
                )
                cond_tmax_lt_allrain = tmax_hru < tmax_allrain_day
                cond_if_if_else = cond_if_if & cond_tmax_lt_allrain
                wh_if_if_else = np.where(cond_if_if_else)
                if len(wh_if_if_else[0]):
                    pptadj[wh_if_if_else] = radj_wppt_day[wh_if_if_else]

                # if.if.if(.if): the cake is taken!
                cond_tmax_gt_allrain_and_not_summer = (
                    tmax_hru >= tmax_allrain_day
                ) & (~is_summer)
                cond_if_if_if = (
                    cond_if_if & cond_tmax_gt_allrain_and_not_summer
                )
                wh_if_if_if = np.where(cond_if_if_if)
                if len(wh_if_if_if[0]):
                    pptadj[wh_if_if_if] = radj_wppt_day[wh_if_if_if]

        radadj = radadj * pptadj
        wh_radadj_lt_2_tenths = np.where(radadj < 0.2)
        if len(wh_radadj_lt_2_tenths[0]):
            radadj[wh_radadj_lt_2_tenths] = 0.2

        swrad = doy_to_daily(soltab_potsw) * radadj / hru_cossl
        return swrad

    # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_potet_jh.f90

    def calculate_potential_et_jh(self) -> None:
        """Calculate potential evapotranspiration following Jensen and Haise
        (1963)."""

        self._potet = self._potet_jh_run(
            dates=self._datetime,
            tmax_hru=self._tmaxf,
            tmin_hru=self._tminf,
            swrad=self._swrad,
            jh_coef=self.jh_coef,
            jh_coef_hru=self.jh_coef_hru,
        )

        return

    @staticmethod
    def _potet_jh_run(dates, tmax_hru, tmin_hru, swrad, jh_coef, jh_coef_hru):
        n_time, n_hru = tmax_hru.shape

        # The vector
        month = (dates.astype("datetime64[M]").astype(int) % 12) + 1
        month = tile_time_to_space(month, n_hru)
        tavg_f = (tmax_hru + tmin_hru) / 2.0
        tavg_c = (tavg_f - 32.0) * 5 / 9  # f_to_c
        elh = (597.3 - (0.5653 * tavg_c)) * inch2cm

        # JLM: move this out?
        def monthly_to_daily(arr):
            return arr[month[:, 0] - 1, :]

        jh_coef_day = monthly_to_daily(jh_coef)
        jh_coef_hru_day = tile_space_to_time(jh_coef_hru, n_time)

        potet = jh_coef_day * (tavg_f - jh_coef_hru_day) * swrad / elh
        cond_potet_lt_zero = potet < zero
        wh_cond = np.where(cond_potet_lt_zero)
        if len(wh_cond[0]):
            potet[wh_cond] = zero

        return potet

    # def advance(self, itime_step, current_date):
    #     self.precip_current = self.precip[itime_step]
    #     self.pot_et_current = self.pot_et[itime_step]
    #     self.pot_et_consumed = 0.0
    #     self.current_date = current_date

    # # Track the amount of potential ET used at a given timestep
    # # JLM: is this strange to track here? I suppose not.
    # def consume_pot_et(self, requested_et):
    #     et = requested_et
    #     available_et = self.pot_et_current - self.pot_et_consumed
    #     if et > available_et:
    #         et = available_et
    #     self.pot_et_consumed += et
    #     return et
