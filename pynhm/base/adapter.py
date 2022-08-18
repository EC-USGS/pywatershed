import pathlib as pl
from typing import Union

import numpy as np

from pynhm.base.control import Control
from pynhm.base.timeseries import TimeseriesArray
from pynhm.constants import fileish
from pynhm.utils.netcdf_utils import NetCdfRead


class Adapter:
    def __init__(
        self,
        variable: str,
    ) -> None:
        self.name = "Adapter"
        self._variable = variable
        return None

    def advance(self):
        raise NotImplementedError("Must be overridden")

    def sanity_check(self):
        raise RuntimeError("contact code developers (sos=0); not really")

    @property
    def current(self):
        return self._current_value


class AdapterNetcdf(Adapter):
    def __init__(
        self,
        fname: fileish,
        variable: str,
        control: Control,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterNetcdf"

        self._fname = fname
        self._nc_read = NetCdfRead(fname)

        self.data = self._nc_read.dataset[self._variable][:].data
        self.time = self._nc_read.times

        self.control = control
        self._start_time = self.control.start_time
        self._current_value = control.get_var_nans(
            self._variable, drop_time_dim=True
        )
        return

    def advance(self):
        if self._nc_read._itime_step[self._variable] > self.control.itime_step:
            return
        self._current_value[:] = self._nc_read.advance(
            self._variable, self.control.current_time
        )
        return None


class AdapterOnedarray(Adapter):
    def __init__(
        self,
        data: np.ndarray,
        variable: str,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterOnedarray"
        self._current_value = data
        return

    def advance(self, *args):
        return None


adaptable = Union[str, np.ndarray, Adapter]


def adapter_factory(
    var: adaptable,
    variable_name: str = None,
    control: Control = None,
):
    if isinstance(var, Adapter):
        """Adapt an adapter"""
        return var

    elif isinstance(var, (str, pl.Path)):
        """Paths and strings are considered paths to netcdf files"""
        if pl.Path(var).suffix == ".nc":
            return AdapterNetcdf(
                var,
                variable=variable_name,
                control=control,
            )

    elif isinstance(var, np.ndarray) and len(var.shape) == 1:
        """Adapt 1-D np.ndarrays"""
        return AdapterOnedarray(var, variable=variable_name)

    elif isinstance(var, TimeseriesArray):
        """Adapt TimeseriesArrays as is."""
        return var

    elif var is None:
        """var is specified as None so return None"""
        return None

    else:
        raise TypeError("oops you screwed up")
