import pathlib as pl
import warnings
from copy import deepcopy
from itertools import chain
from typing import Union

import netCDF4 as nc4

# from numba import jit
import numpy as np

from ..preprocess.cbh_utils import cbh_adjust
from ..utils.parameters import PrmsParameters
from .AtmBoundaryLayer import AtmBoundaryLayer
from .NHMSolarGeometry import NHMSolarGeometry

fileish = Union[str, pl.Path]

zero = np.zeros(1)[0]
one = np.ones(1)[0]

inch_to_cm = 2.54

# https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/c_solar_radiation_degday.f90#L19
# define this somehow
solf = np.array(
    [
        0.20,
        0.35,
        0.45,
        0.51,
        0.56,
        0.59,
        0.62,
        0.64,
        0.655,
        0.67,
        0.682,
        0.69,
        0.70,
        0.71,
        0.715,
        0.72,
        0.722,
        0.724,
        0.726,
        0.728,
        0.73,
        0.734,
        0.738,
        0.742,
        0.746,
        0.75,
    ]
)


parameters_swrad = [
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
]

# JLM: the swrad parameters should be available from that module.

params_required = {
    "swrad": parameters_swrad,
    "solar_geom": ["nhru", "hru_aspect"] + parameters_swrad,
    "pot_et_jh": ["jh_coef", "jh_coef_hru"],
    "adj": [
        "nhru",
        "tmax_cbh_adj",
        "tmin_cbh_adj",
        "tmax_allsnow",
        "tmax_allrain_offset",
        "snow_cbh_adj",
        "rain_cbh_adj",
        "adjmix_rain",
    ],
}


# JLM: where is metadata


# may not use this if they cant be called with jit
# if it can, put it in a common utility.
def tile_space_to_time(arr: np.ndarray, n_time) -> np.ndarray:
    return np.tile(arr, (n_time, 1))


def tile_time_to_space(arr: np.ndarray, n_space) -> np.ndarray:
    return np.transpose(np.tile(arr, (n_space, 1)))


class NHMBoundaryLayer(AtmBoundaryLayer):
    def __init__(
        self,
        state_dict: dict,
        parameters: PrmsParameters,
        *args,
        **kwargs,
    ):
        """The atmospheric boundary layer of the NHM model."""
        super().__init__(*args, **kwargs)

        self.name = "NHMBoundaryLayer"

        # Dimensions and dimension variables
        self.spatial_id = None
        self.spatial_id_is_nhm = None
        # property hrus or space & nspace, time &ntime

        # JLM: these names are historical, I dont love them
        self._potential_variables = [
            "prcp",
            "rhavg",
            "rainfall",
            "snowfall",
            "tmax",
            "tmin",
            "swrad",
            "potet",
        ]

        self._parameter_list = list(
            set(
                chain.from_iterable(
                    [val for key, val in params_required.items()]
                )
            )
        )
        self.parameters = parameters.get_parameters(self._parameter_list)

        self._n_states_adj = 0
        self._allow_param_adjust = None
        self._state_adjusted_here = False

        # To be calculated (diagnostic?)
        self.pot_et = None
        self.pot_et_consumed = None

        # These are set in super(), delete them now
        # to be able to set
        del self.datetime
        del self.spatial_id

        for key, val in state_dict.items():
            self[key] = val

        return

    @classmethod
    def load_netcdf(
        cls,
        nc_file: fileish,
        parameters: PrmsParameters,
        *args,
        nc_read_vars: list = None,
        **kwargs,
    ) -> "NHMBoundaryLayer":
        obj = cls({}, *args, parameters=parameters, **kwargs)

        # netcdf handling. consolidate these?
        obj.dataset = None
        obj.ds_var_list = None
        obj.ds_var_chunking = None

        obj._self_file_vars = {}
        obj._nc_read_vars = nc_read_vars
        obj._optional_read_vars = ["rhavg", "swrad", "potet"]
        obj._nc_file = nc_file
        obj._open_nc_file()
        obj._read_nc_file()
        # obj._close_nc_file()

        return obj

    # move this down? create helper functions? it is too long to be near top.
    def _open_nc_file(self):
        self.dataset = nc4.Dataset(self._nc_file, "r")
        self.ds_var_list = list(self.dataset.variables.keys())
        self.ds_var_chunking = {
            vv: self.dataset.variables[vv].chunking()
            for vv in self.ds_var_list
        }
        # Set dimension variables which are not chunked
        self["datetime"] = (
            nc4.num2date(
                self.dataset.variables["datetime"][:],
                units=self.dataset.variables["datetime"].units,
                calendar=self.dataset.variables["datetime"].calendar,
                only_use_cftime_datetimes=False,
            )
            .filled()
            .astype("datetime64[s]")
            # JLM: the global time type as in cbh_utils, define somewhere
        )
        spatial_id_name = (
            "nhm_id" if "nhm_id" in self.ds_var_list else "hru_ind"
        )
        self.spatial_id = self.dataset.variables[spatial_id_name][:]

        if self._nc_read_vars is None:

            # If specific variables are not specified for reading,
            # then set the list of variables in the file to use.
            # Assumption is: use adjusted variables over unadjusted ones
            # The logic below just sets this perference for adjusted
            # variables over unadjusted vars which are warned about when
            # used as backup
            n_states_adj = 0
            for self_var in self._potential_variables:
                adj_var = f"{self_var}_adj"
                if adj_var in self.ds_var_list:
                    self._self_file_vars[self_var] = adj_var
                    n_states_adj += 1
                elif self_var in self.ds_var_list:
                    self._self_file_vars[self_var] = self_var
                    msg = (
                        f"Adjusted variable '{adj_var}' not found in dataset, "
                        f"using apparently unadjusted variable instead: "
                        f"'{self_var}'"
                    )
                    warnings.warn(msg)
                else:
                    msg = (
                        f"Neither variable '{self_var}' nor {adj_var}"
                        f"found in the file: '{self._nc_file}'"
                    )
                    if self_var in self._optional_read_vars:
                        warnings.warn(msg)
                    else:
                        raise ValueError(msg)

        else:

            # Find the specified/requested variables to read from file.
            n_states_adj = 0
            for self_var in self._potential_variables:
                adj_var = f"{self_var}_adj"
                if adj_var in self._nc_read_vars:
                    self._self_file_vars[self_var] = adj_var
                    n_states_adj += 1
                elif self_var in self._nc_read_vars:
                    self._self_file_vars[self_var] = self_var
                else:
                    continue

            # Post-mortem on what requested varaibles were not found:
            for read_var in self._nc_read_vars:
                if read_var not in list(self._self_file_vars.values()):
                    msg = (
                        f"Requested variable '{read_var}' not found in file "
                        f"'{self._nc_file}'"
                    )
                    raise ValueError(msg)

            # finally reduce the state vars to what we have
            for state_var in deepcopy(self.variables):
                if state_var not in self._self_file_vars.keys():
                    del self.__dict__[state_var]
                    _ = self.variables.remove(state_var)

        # If it looks like we read adjusted parameters, try to block double
        # adjustments on parameter_adjust
        if n_states_adj > 0:
            self._allow_param_adjust = False
        else:
            self._allow_param_adjust = True

        return

    def _read_nc_file(self) -> None:
        # JLM: "front load" option vs "load as you go"
        for self_var, file_var in self._self_file_vars.items():
            # JLM: Later use chunking here to get data
            # JLM: have to update datetime in that case?
            data = self.dataset.variables[file_var][:]
            if data.mask:
                msg = (
                    "NetCDF file contains missing values "
                    "not currently handled"
                )
                raise ValueError(msg)
            self[self_var] = data.data
        return

    def _close_nc_file(self) -> None:
        self.dataset.close()
        return

    def param_adjust(self) -> None:
        """Adjust atmospheric variables using PRMS parameters.

        Args:
            None

        Returns:
            None
        """
        if not self._allow_param_adjust:
            msg = (
                "Parameter adjustments not permitted when any state "
                " variables are already adjusted. \n {self._self_file_vars}"
            )
            raise ValueError(msg)

        # construct a dict of references
        var_dict_adj = {key: self[key] for key in self.variables}
        var_dict_adj["datetime"] = self["datetime"]
        # This modified var_dict_adj, which is sus
        params = self.parameters.get_parameters(params_required["adj"])
        _ = cbh_adjust(var_dict_adj, params)
        self._state_adjusted_here = True
        self._allow_param_adjust = False

        map_adj_dict = {
            "prcp": "prcp_adj",
            "tmax": "tmax_adj",
            "tmin": "tmin_adj",
            "rainfall": "rainfall_adj",
            "snowfall": "snowfall_adj",
        }

        for self_var, adj_var in map_adj_dict.items():
            self[self_var] = var_dict_adj[adj_var]

        return

    # JLM: make this available for preprocess with cbh
    def calculate_sw_rad_degree_day(self) -> None:
        """Calculate shortwave radiation using the degree day method."""

        solar_geom = NHMSolarGeometry(self.parameters)
        params = self.parameters.get_parameters(
            params_required["swrad"]
        ).parameters

        self["swrad"] = self._ddsolrad_run(
            dates=self["datetime"],
            tmax_hru=self["tmax"],
            hru_ppt=self["prcp"],
            soltab_potsw=solar_geom["potential_sw_rad"],
            **params,
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

        # JLM: could get this from solar_geom object or delete it there
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
            kp = dday.astype(int)  # dddayi = float(kp)
            radadj[wh_dday_lt_26] = (
                solf[kp - 1] + ((solf[kp] - solf[kp - 1]) * (dday - kp))
            )[wh_dday_lt_26]
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

    def calculate_potential_et_jh(
        self,
    ) -> None:
        """Calculate potential evapotranspiration following Jensen and Haise (1963)."""
        params = self.parameters.get_parameters(
            params_required["pot_et_jh"]
        ).parameters

        self["potet"] = self._potet_jh_run(
            dates=self["datetime"],
            tmax_hru=self["tmax"],
            tmin_hru=self["tmin"],
            swrad=self["swrad"],
            **params,
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
        elh = (597.3 - (0.5653 * tavg_c)) * inch_to_cm

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
