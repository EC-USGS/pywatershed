import pathlib as pl
from typing import Union

import numpy as np

from pynhm.base.control import Control
from pynhm.utils.netcdf_utils import NetCdfRead

fileish = Union[str, pl.Path]


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
        self._dataset = NetCdfRead(fname)
        self._datetime = self._dataset.date_times

        if control is None:
            raise ValueError(f"{self.name} requires a valid Control object.")
        self.control = control
        self._start_time = self.control.start_time

        self._nhru = self._dataset.dataset.dimensions["nhm_id"].size
        self._dtype = self._dataset.dataset.variables[self._variable].dtype
        self._current_value = np.full(self._nhru, np.nan, self._dtype)
        # Adopting the next line requires minor changes to most tests
        # self._current_value = control.get_var_nans(self._variable)

    def advance(self):
        # JLM: Seems like the time of the ncdf dataset or variable
        # should be public
        if self._dataset._itime_step[self._variable] > self.control.itime_step:
            return
        self._current_value[:] = self._dataset.advance(self._variable)
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
        """One-D np.ndarrays"""
        print(f"1darray adapter: {variable_name}")
        return AdapterOnedarray(var, variable=variable_name)

    elif var is None:
        """var is specified as None so return None"""
        return None

    else:
        raise TypeError("oops you screwed up")
