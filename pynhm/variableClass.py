import pathlib as pl
from typing import Union

import numpy as np

from .utils.netcdf_utils import NetCdfRead

fileish = Union[str, pl.Path]


def variable_factory(
    var,
    variable_name: str = None,
    start_time: np.datetime64 = None,
):
    if isinstance(var, (str, pl.Path)):
        if pl.Path(var).suffix == ".nc":
            return VariableFromNetcdf(
                var,
                variable=variable_name,
                start_time=start_time,
            )
    else:
        raise TypeError("oops you screwed up")


class Variable:
    def __init__(
        self,
        variable: str,
    ) -> None:
        self._variable = variable
        return None

    def advance(self):
        raise NotImplementedError("Must be overridden")

    def sanity_check(self):
        raise RuntimeError("contact code developers (sos=0); not really")


class VariableFromNetcdf(Variable):
    def __init__(
        self,
        fname: fileish,
        variable: str,
        start_time: np.datetime64 = None,
    ) -> None:
        super().__init__(variable)
        self._fname = fname
        self._dataset = NetCdfRead(fname)
        self._datetime = self._dataset.date_times
        if start_time is None:
            self._start_time = self._datetime[0]
            self._itime_step = 0
        else:
            self._start_time = start_time
            self._itime_step = np.where(self._datetime == start_time)[0]

        self._previous_value = None
        self._current_value = None

    def advance(self):
        self._previous_value = self._current_value
        self._current_value = self._dataset.advance(self._variable)
        return None

    @property
    def current(self):
        return self._current_value
