import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np

fileish = Union[str, pl.Path]


class NetCdfRead:
    def __init__(
        self,
        nc_file: fileish,
        nc_read_vars: list = None,
    ) -> "NetCdfRead":

        self._nc_file = nc_file
        self._nc_read_vars = nc_read_vars
        self._itime_step = 0
        self._open_nc_file()

    def __del__(self):
        if self.dataset.isopen():
            self.dataset.close()

    def _open_nc_file(self):
        self.dataset = nc4.Dataset(self._nc_file, "r")
        self.ds_var_list = list(self.dataset.variables.keys())
        if self._nc_read_vars is None:
            self._nc_read_vars = self.ds_var_list
        self.ds_var_chunking = {
            vv: self.dataset.variables[vv].chunking()
            for vv in self.ds_var_list
        }
        # Set dimension variables which are not chunked
        self._datetime = (
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
        self._ntimes = self._datetime.shape[0]
        spatial_id_name = (
            "nhm_id" if "nhm_id" in self.ds_var_list else "hru_ind"
        )
        self.spatial_id = self.dataset.variables[spatial_id_name][:]

    @property
    def ntimes(self):
        return self._ntimes

    def get_data(
        self,
        variable: str,
        itime_step: int = None,
    ) -> np.ndarray:
        if variable not in self._nc_read_vars:
            raise ValueError(
                f"'{variable}' not in list of available variables"
            )
        if itime_step is None:
            return self.dataset[variable][:, :]
        else:
            if itime_step >= self._ntimes:
                raise ValueError(
                    f"requested time step {itime_step} but only "
                    + f"{self._ntimes} time steps are available."
                )
            return self.dataset[variable][itime_step, :]

    def advance(
        self,
        variable: str,
    ) -> np.ndarray:
        arr = self.get_data(
            variable,
            itime_step=self._itime_step,
        )
        self._itime_step += 1
        return arr
