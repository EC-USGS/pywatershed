import warnings
from copy import deepcopy

import netCDF4 as nc4

# from numba import jit
import numpy as np

from ..preprocess.cbh_utils import cbh_adjust
from ..utils.parameters import PrmsParameters
from .AtmBoundaryLayer import AtmBoundaryLayer
from .NHMSolarGeometry import NHMSolarGeometry

# JLM: where is metadata


class NHMBoundaryLayer(AtmBoundaryLayer):
    def __init__(self, nc_file, *args, nc_read_vars=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.name = "NHMBoundaryLayer"
        self._nc_file = nc_file
        self._nc_read_vars = nc_read_vars

        # Dimensions and dimension variables
        self.datetime = None
        self.spatial_id = None
        self.spatial_id_is_nhm = None
        # property hrus or space & nspace, time &ntime

        # The states
        # JLM: these names are historical, I dont love them
        self._potential_variables = [
            "prcp",
            "rhavg",
            "rainfall",
            "snowfall",
            "tmax",
            "tmin",
            "swrad",
        ]
        self.variables = []
        self._self_file_vars = {}
        self._optional_read_vars = ["swrad", "rhavg"]
        self._n_states_adj = 0
        self.prcp = None
        self.rhavg = None
        self.rainfall = None
        self.snowfall = None
        self.tmax = None
        self.tmin = None

        self._allow_param_adjust = None
        self._state_adjusted_here = False

        # To be calculated (diagnostic?)
        self.pot_et = None
        self.pot_et_consumed = None

        # netcdf handling. consolidate these?
        self.dataset = None
        self.ds_var_list = None
        self.ds_var_chunking = None
        self._open_nc_file()
        self._read_nc_file()
        # self._close_nc_file()
        return

    # move this down? create helper functions? it is too long to be near top.
    def _open_nc_file(self):
        self.dataset = nc4.Dataset(self._nc_file, "r")
        self.ds_var_list = list(self.dataset.variables.keys())
        self.ds_var_chunking = {
            vv: self.dataset.variables[vv].chunking()
            for vv in self.ds_var_list
        }
        # Set dimension variables which are not chunked
        self.datetime = (
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
            _ = self[self_var] = self.dataset.variables[file_var][:]
        return

    def _close_nc_file(self) -> None:
        self.dataset.close()
        return

    def param_adjust(self, parameters: PrmsParameters) -> None:
        if not self._allow_param_adjust:
            msg = (
                "Parameter adjustments not permitted when any state "
                " variables are already adjusted. \n {self._self_file_vars}"
            )
            raise ValueError(msg)

        # construct a dict of references
        var_dict_adj = {key: self[key] for key in self.variables}
        var_dict_adj["datetime"] = self.datetime
        self._parameters_for_adj = parameters
        # This modified var_dict_adj, which is sus
        _ = cbh_adjust(var_dict_adj, self._parameters_for_adj)
        self._state_adjusted_here = True
        self._allow_param_adjust = False

        map_adj_dict = {
            "datetime": "datetime",
            "prcp": "prcp_adj",
            "tmax": "tmax_adj",
            "tmin": "tmin_adj",
            "rainfall": "rainfall_adj",
            "snowfall": "snowfall_adj",
        }

        for self_var, adj_var in map_adj_dict.items():
            self[self_var] = var_dict_adj[adj_var]
            self.variables = list(set(self.variables))

        return

    # JLM: make this available for preprocess with cbh
    def calculate_sw_rad_degree_day(
        self,
        parameters: PrmsParameters,
    ) -> None:

        solar_geom = NHMSolarGeometry(parameters)
        asdf
        self["swrad"] = self.ddsolrad_run(
            self["datetime"],
            self["tmax"],  # adj?
            self["prcp"],
            solar_geom.soltab_potsw,
        )
        return

    # @jit
    @staticmethod
    def ddsolrad_run(
        dates: np.ndarray,  # [n_time]
        tmax_hru: np.ndarray,  # [n_time, n_hru]
        hru_ppt: np.ndarray,  # [n_time, n_hru]
        soltab_potsw: np.ndarray,  # [n_time, n_hru] ??
        radadj_intcp,  # param
        radadj_slope,  # param
        tmax_index,  # param
        dday_slope,  # param
        dday_intcp,  # param
        radmax,  # param
        hru_slope,  # param
        ppt_rad_adj,  # param
        radj_sppt,  # param
        radj_wppt,  # param
        tmax_allrain,  # params
        hemisphere,  # make internal based on hru lat?
    ) -> np.ndarray:  # [n_time, n_hru]

        # PRMS6 disagrees with Markstrom's code
        # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation_degday.f90

        swrad = np.zeros(len(dates))

        ihru = 0
        hru_cossl = math.cos(math.atan(hru_slope[ihru]))

        ii = 0
        for date in dates:
            jday = day_of_year(date)
            imon = date.month - 1

            # ! set degree day and radiation adjustment limited by radmax
            dday = (
                dday_slope[imon, ihru] * tmax_hru[ii]
                + dday_intcp[imon, ihru]
                + 1.0
            )
            if dday < 1.0:
                dday = 1.0

            if dday < 26.0:
                kp = int(dday)
                ddayi = float(kp)
                kp1 = kp + 1
                radadj = solf[kp - 1] + (
                    (solf[kp1 - 1] - solf[kp - 1]) * (dday - ddayi)
                )
                if radadj > radmax[imon, ihru]:
                    radadj = radmax[imon, ihru]
            else:
                radadj = radmax[imon, ihru]

            #           ! Set precipitation adjument factor based on temperature
            #           ! and amount of precipitation

            pptadj = 1.0
            if hru_ppt[ii] > ppt_rad_adj[imon, ihru]:
                if tmax_hru[ii] < tmax_index[imon, ihru]:
                    pptadj = radj_sppt[ihru]
                    if tmax_hru[ii] >= tmax_allrain[imon]:
                        if not is_summer(jday, hemisphere):
                            pptadj = radj_wppt[ihru]
                    else:
                        pptadj = radj_wppt[ihru]
                else:
                    pptadj = radadj_intcp[imon, ihru] + radadj_slope[
                        imon, ihru
                    ] * (tmax_hru[ii] - tmax_index[imon, ihru])
                    if pptadj > 1.0:
                        pptadj = 1.0

            radadj = radadj * pptadj
            if radadj < 0.2:
                radadj = 0.2

            swrad[ii] = soltab_potsw[jday - 1] * radadj / hru_cossl
            ii += 1

        return swrad

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

    # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_potet_jh.f90
    # is this part of the full class?
    # # potet_jh run code:
    # INCH2CM = 2.54

    # def potet_jh_run(dates, tmax_hru, tmin_hru, swrad, jh_coef, jh_coef_hru):
    #     potet = np.zeros(len(dates))

    #     ihru = 0
    #     ii = 0
    #     for date in dates:
    #         imon = date.month - 1
    #         tavgf = (tmax_hru[ii] + tmin_hru[ii]) / 2.0
    #         tavgc = (tavgf - 32.0) * 5 / 9
    #         elh = (597.3 - (0.5653 * tavgc)) * INCH2CM
    #         potet[ii] = jh_coef[imon, ihru] * (tavgf - jh_coef_hru[ihru]) * swrad[ii] / elh
    #         if potet[ii] < 0.0:
    #             potet[ii] = 0.0
    #         ii += 1
    #     return potet
