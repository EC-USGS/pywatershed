import pathlib as pl
from typing import Union

import numpy as np

from pynhm.base.control import Control
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
        control: Control,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterNetcdf"

        self._fname = fname
        self._dataset = NetCdfRead(fname)
        self._datetime = self._dataset.date_times

        if control is None:
            raise ValueError(f"{self.name} requires a valid Control object.")
        self.control = control
        self._start_time = self.control.start_time

        self._nhru = self._dataset.dataset.dimensions["nhm_id"].size
        self._dtype = self._dataset.dataset.variables[self._variable].dtype
        self._current_value = np.full(self._nhru, np.nan, self._dtype)

    def advance(self):
        # JLM: Seems like the time of the ncdf dataset or variable
        # should be public
        if self._dataset._itime_step[self._variable] > self.control.itime_step:
            return
        self._current_value[:] = self._dataset.advance(self._variable)
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
        """One-D np.ndarrays"""
        return AdapterOnedarray(var, variable=variable_name)
    elif (
        isinstance(var, np.ndarray)
        and len(var.shape) == 2
        and variable_name in soltab_vars
    ):
        """Two-D np.ndarrays that are Soltabs"""
        return AdapterTwodarraySoltab(var, variable=variable_name)
    elif var is None:
        """var is specified as None so return None"""
        # using the variable_name, get the dimensions from metadata
        # and their size from control
        # this requires "params" to be available
        # I propose adding them to control (makes sense to me,
        # control already has config)
        var_dims = control.meta.get_dimensions(variable_name)[variable_name]
        var_dim_sizes = [control.params.parameters[vv] for vv in var_dims]
        var_type = control.meta.get_numpy_types(variable_name)[variable_name]
        # return None
        return np.full(var_dim_sizes, np.nan, var_type)

    else:
        raise TypeError("oops you screwed up")
