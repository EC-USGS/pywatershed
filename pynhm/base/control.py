"""The control class."""

import pathlib as pl
from typing import Union

import numpy as np

from ..base.meta import Meta
from ..utils import ControlVariables
from ..utils.time_utils import (
    datetime_dowy,
    datetime_doy,
    datetime_month,
    datetime_year,
)
from .accessor import Accessor

fileish = Union[str, pl.PosixPath]


class Control(Accessor):
    """The control class."""

    def __init__(
        self,
        start_time: np.datetime64,
        end_time: np.datetime64,
        time_step: np.timedelta64,
        config: dict = None,
        verbosity: int = 0,
        **kwargs,
    ):
        """Initialize time with data and parameters.

        Args:
            start_time: this is the first time of integration NOT the restart time
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

        self._init_time = None  # get from control file if restart?
        self._current_time = None
        self._previous_time = None
        self._i_time = 0

        self.config = config

        self.meta = Meta()
        # This will have the time dimension name
        # This will have the time coordimate name

    @classmethod
    def load(
        cls,
        control_file: fileish,
        verbosity: int = 0,
    ) -> "Time":
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
            verbosity=verbosity,
        )

    @property
    def current_time(self):
        """Get the current time."""
        return self._current_time

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
    def previous_time(self):
        """Get the previous time."""
        return self._previous_time

    @property
    def i_time(self):
        """The counth of the current time [0, self.n_times-1]"""
        return self._i_time

    @property
    def time_step(self):
        """Get the time step."""
        return self._time_step

    @property
    def start_time(self):
        """Get the simulation start time"""
        return self._start_time

    @property
    def end_time(self):
        """Get the simulation end time"""
        return self._end_time

    @property
    def n_times(self):
        """Get the number of times"""
        return self._n_times

    def advance(self):
        """Advance time."""

        if self._current_time is None:
            self._current_time = self._start_time
            return

        if self._current_time == self._end_time:
            raise ValueError("End of time reached")

        self._previous_time = self._current_time
        self._current_time += self.time_step
        self._i_time += 1

        return None
