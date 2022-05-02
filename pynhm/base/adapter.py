import pathlib as pl
from typing import Union

import numpy as np

from pynhm.utils.netcdf_utils import NetCdfRead

from ..utils.time_utils import datetime_doy

fileish = Union[str, pl.Path]

# these names will change and reduce with amt/geom refactor
# could get these from SolarGeom but there are other non-soltab vars
soltab_vars = [
    "soltab_horad_potsw",
    "potential_sw_rad_flat",
    "potential_sw_rad",
]


def adapter_factory(
    var,
    variable_name: str = None,
    start_time: np.datetime64 = None,
):
    if isinstance(var, (str, pl.Path)):
        if pl.Path(var).suffix == ".nc":
            return AdapterNetcdf(
                var,
                variable=variable_name,
                start_time=start_time,
            )
    elif isinstance(var, np.ndarray) and len(var.shape) == 1:
        """One-D np.ndarrays"""
        return AdapterOnedarray(var, variable=variable_name)
    elif (
        isinstance(var, np.ndarray)
        and len(var.shape) == 2
        and variable_name in soltab_vars
    ):
        """Two-D np.ndarrays that are Soltabs"""
        return AdapterTwodarraySoltab(var, variable=variable_name)
    else:
        raise TypeError("oops you screwed up")


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


class AdapterNetcdf(Adapter):
    def __init__(
        self,
        fname: fileish,
        variable: str,
        start_time: np.datetime64 = None,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterNetcdf"
        self._fname = fname
        self._dataset = NetCdfRead(fname)
        self._datetime = self._dataset.date_times
        if start_time is None:
            self._start_time = self._datetime[0]
        else:
            self._start_time = start_time

        self._current_value = None

    def advance(self):
        self._current_value = self._dataset.advance(self._variable)
        return None

    @property
    def current(self):
        return self._current_value


class AdapterOnedarray(Adapter):
    def __init__(
        self,
        data: np.ndarray,
        variable,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterOnedarray"
        self._current_value = data
        return

    def advance(self):
        return None

    @property
    def current(self):
        return self._current_value


class AdapterTwodarraySoltab(Adapter):
    def __init__(
        self,
        data: np.ndarray,
        variable: str,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterTwodarraySoltab"
        self._data = data
        self._current_value = None
        return

    def advance(self, datetime: np.datetime64):
        index = datetime_doy(datetime) - 1
        self._current_value = self._data[index]
        return None

    @property
    def current(self):
        return self._current_value
