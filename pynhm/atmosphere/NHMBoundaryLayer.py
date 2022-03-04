import netCDF4 as nc4
import numpy as np
import warnings

# JLM: "front load" option vs "load as you go"
# JLM: metadata
from .AtmBoundaryLayer import AtmBoundaryLayer


class NHMBoundaryLayer(AtmBoundaryLayer):
    def __init__(self, nc_file, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._nc_file = nc_file

        # Dimensions and dimension variables
        self.datetime = None
        self.hru_id = None
        self.hru_id_is_nhm = None
        # property hrus or space & nspace, time &ntime

        # The states
        self.state_vars_in_file = [
            "prcp",
            "rainfall",
            "snowfall",
            "tmax",
            "tmin",
        ]
        self.prcp = None
        self.rainfall = None
        self.snowfall = None
        self.tmax = None
        self.tmin = None

        # To be calculated (diagnostic?)
        self.pot_et = None
        self.pot_et_consumed = None

        # netcdf handling. consolidate these?
        self.dataset = None
        self.ds_var_list = None
        self.ds_var_chunking = None
        self._open_nc_file()
        self._read_nc_file()
        return

    # Set state from the netcdf file
    # get the adjusted states? depends on how these are written out

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
        hru_id_name = "nhm_id" if "nhm_id" in self.ds_var_list else "hru_ind"
        self.hru_id = self.dataset.variables[hru_id_name][:]
        return

    def _read_nc_file(self):
        for self_var in self.state_vars_in_file:
            file_var = f"{self_var}_adj"
            msg = None
            if file_var not in self.ds_var_list:
                msg = (
                    f"Adjusted variable '{file_var}' not found in dataset, "
                    f"using apparently unadjusted variable instead: '{self_var}'"
                )
                file_var = self_var
            # JLM: Later use chunking here
            _ = setattr(self, self_var, self.dataset.variables[file_var][:])
            if msg is not None:
                warnings.warn(msg)
        return

    # def param_adjust(self):
    #     msg = (
    #         "Base AtmosphericForcings class does"
    #         "not provide forcing adjustment"
    #     )
    #     raise NotImplementedError(msg)

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

    # def adjust(self):
    #     self._adjust_temp()
    #     self._adjust_precip()


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
