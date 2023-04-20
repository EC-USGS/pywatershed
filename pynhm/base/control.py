"""The control class."""
import datetime

import numpy as np

from ..base import meta
from ..constants import fileish
from ..utils import ControlVariables
from .parameters import Parameters
from ..utils.time_utils import (
    datetime_dowy,
    datetime_doy,
    datetime_epiweek,
    datetime_month,
    datetime_year,
)
from .accessor import Accessor


class Control(Accessor):
    """The control class."""

    def __init__(
        self,
        start_time: np.datetime64,
        end_time: np.datetime64,
        time_step: np.timedelta64,
        init_time: np.datetime64 = None,
        config: dict = None,
        params: Parameters = None,
        verbosity: int = 0,
        **kwargs,
    ):
        """Initialize time with data and parameters.

        Args:
            start_time: this is the first time of integration NOT the restart
                time
            end_time: the last integration time
            time_step: the length fo the time step
            config: a PRMS config file to read and use for contorl
            verbosity: the level of verbosity in [0,10]
        """
        super().__init__(**kwargs)
        self.name = "Control"

        self.verbosity = verbosity

        if end_time <= start_time:
            raise ValueError("end_time <= start_time")

        n_times_m1 = (end_time - start_time) / time_step
        if n_times_m1 != int(n_times_m1):
            raise ValueError("time_step does not divide end_time - start_time")

        self._start_time = start_time
        self._end_time = end_time
        self._time_step = time_step
        self._n_times = int(n_times_m1) + 1

        if init_time:
            self._init_time = init_time
        else:
            self._init_time = self._start_time - time_step

        self._current_time = None
        self._previous_time = None
        self._itime_step = -1

        self.config = config
        self.params = params

        self.meta = meta
        # This will have the time dimension name
        # This will have the time coordimate name

    @classmethod
    def load(
        cls,
        control_file: fileish,
        params: Parameters = None,
        verbosity: int = 0,
    ) -> "Control":
        """Initialize a control object from a PRMS control file

        Args:
            control_file: PRMS control file
            verbosity: output verbosity level

        Returns:
            Time: Time object initialized from a PRMS control file

        """
        control = ControlVariables.load(control_file)

        return cls(
            control.control["start_time"],
            control.control["end_time"],
            control.control["initial_deltat"],
            config=control.control,
            params=params,
            verbosity=verbosity,
        )

    @property
    def current_time(self):
        """Get the current time."""
        return self._current_time

    @property
    def current_datetime(self):
        """Get the current time as a datetime.datetime object"""
        return self._current_time.astype(datetime.datetime)

    @property
    def current_year(self):
        """Get the current year."""
        return datetime_year(self._current_time)

    @property
    def current_month(self):
        """Get the current month."""
        return datetime_month(self._current_time)

    @property
    def current_doy(self):
        """Get the current day of year."""
        return datetime_doy(self._current_time)

    @property
    def current_dowy(self):
        """Get the current day of water year."""
        return datetime_dowy(self._current_time)

    @property
    def current_epiweek(self):
        """Get the current epiweek [1, 53]."""
        return datetime_epiweek(self._current_time)

    @property
    def previous_time(self):
        """The previous time."""
        return self._previous_time

    @property
    def itime_step(self):
        """The counth of the current time [0, self.n_times-1]"""
        return self._itime_step

    @property
    def time_step(self):
        """The time step"""
        return self._time_step

    @property
    def init_time(self):
        """Get the simulation initialization time"""
        return self._init_time

    @property
    def start_time(self):
        """The simulation start time"""
        return self._start_time

    @property
    def start_doy(self):
        """The simulation start day of year"""
        return datetime_doy(self._start_time)

    @property
    def start_month(self):
        """The simulation start month"""
        return datetime_month(self._start_time)

    @property
    def end_time(self):
        """The simulation end time"""
        return self._end_time

    @property
    def n_times(self):
        """The number of time steps"""
        return self._n_times

    @property
    def time_step_seconds(self):
        return self.time_step / np.timedelta64(1, "s")

    def advance(self):
        """Advance time"""
        if self._current_time == self._end_time:
            raise ValueError("End of time reached")

        self._itime_step += 1

        if self._current_time is None:
            self._current_time = self._start_time
            return

        self._previous_time = self._current_time
        self._current_time += self.time_step

        return None

    def get_var_nans(self, var_name: str, drop_time_dim: bool = None):
        """Get an array filled with nans for a given variable"""
        var_dims = self.meta.get_dimensions(var_name)[var_name]
        if drop_time_dim:
            # This accomodates Timeseries like objects that need to init both
            # full rank and reduced rank versions of their data
            # this is pretty adhoc
            check_list = ["time", "doy"]
            if len([mm for mm in check_list if mm in var_dims[0]]):
                del var_dims[0]

        var_dim_sizes = self.params.get_dim_values(var_dims)
        var_dim_shape = [var_dim_sizes[vv] for vv in var_dims]
        var_type = self.meta.get_numpy_types(var_name)[var_name]
        return np.full(var_dim_shape, np.nan, var_type)
