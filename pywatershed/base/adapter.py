import pathlib as pl
from typing import Union

import numpy as np

from ..base.control import Control
from ..base.timeseries import TimeseriesArray
from ..constants import fileish
from ..utils.netcdf_utils import NetCdfRead


class Adapter:
    """Adapter base class for getting data from a variety of sources.

    Args:
        variable: string name ov variable
    """

    def __init__(
        self,
        variable: str,
    ) -> None:
        self.name = "Adapter"
        self._variable = variable
        return None

    def advance(self):
        "Advance the adapter in time"
        raise NotImplementedError("Must be overridden")

    @property
    def current(self):
        return self._current_value


class AdapterNetcdf(Adapter):
    """Adapter class for a netcdf file

    Args:
        fname: filename of netcdf as string or Path
        variable: variable name string
        dim_sizes: a tuple of dimension sizes
        type: a variable dtype
        control: a Control object
        load_n_time_batches: number of times to read from file.

    """

    def __init__(
        self,
        fname: fileish,
        variable: str,
        dim_sizes: tuple,
        type: str,
        control: Control,
        load_n_time_batches: int = 1,
    ) -> None:
        super().__init__(variable)
        self.name = "AdapterNetcdf"

        self._dim_sizes = dim_sizes
        self._type = type
        self._fname = fname
        self.control = control
        self._start_time = self.control.start_time
        self._end_time = self.control.end_time

        self._nc_read = NetCdfRead(
            fname,
            start_time=self._start_time,
            end_time=self._end_time,
            load_n_time_batches=load_n_time_batches,
        )

        # would like to make this a check if dim_sizes and type are available
        nc_type = self._nc_read.dataset[variable].dtype
        nc_shape = list(self._nc_read.dataset[variable].shape)
        nc_dims = list(self._nc_read.dataset[variable].dimensions)
        for time_dim in ["time", "doy"]:
            if time_dim in nc_dims:
                _ = nc_shape.pop(nc_dims.index(time_dim))

        self.time = self._nc_read.times
        # self._current_value = np.full(self._dim_sizes, np.nan, self._type)
        self._current_value = np.full(nc_shape, np.nan, nc_type)

        return

    def advance(self):
        if self._nc_read._itime_step[self._variable] > self.control.itime_step:
            return
        self._current_value[:] = self._nc_read.advance(
            self._variable, self.control.current_time
        )
        return None

    @property
    def data(self):
        # TODO JLM: seems like we'd want to cache this data if we invoke once
        return self._nc_read.all_time(self._variable).data


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


adaptable = Union[str, pl.Path, np.ndarray, Adapter]


def adapter_factory(
    var: adaptable,
    variable_name: str = None,
    control: Control = None,
    variable_dim_sizes: tuple = None,
    variable_type: str = None,
    load_n_time_batches: int = 1,
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
                dim_sizes=variable_dim_sizes,
                type=variable_type,
                load_n_time_batches=load_n_time_batches,
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
