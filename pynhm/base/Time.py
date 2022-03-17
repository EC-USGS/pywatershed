"""The base time class.

Currently, all objects with a time dimension inherit from this. (Consider
composition instead).
"""

import pathlib as pl
from typing import Tuple, Union

import numpy as np

from ..utils import ControlVariables
from .StateAccess import StateAccess

fileish = Union[str, pl.PosixPath]


class Time(StateAccess):
    """The time base class.

    To manage and access time the same way across the model.

    You may initialize with or with out datetime coordinate.

    * datetime supplied: The Time object behaves like "timeseries" where it
      keeps track of time relative to the datetime data.

    * datetime NOT supplied: The Time object behaves like a Markov Model in
      the sense that only the current and previous times are kept.

    Dimension:
        'time'
    Coordinate data variable (optional):
        'datetime': a one-dimensional np.ndarray of type np.datetime64
        The length of the datetime data is the number of times
        to maintain in state. (This could be changed in the future
        with an optional argument).

    Parameters:
        start_time: np.datetime64 scalar for the simulation start time.
        end_time: np.datetime64 scalar for the simulation end time.
        time_step: nptimedelta64 for the distance between times.
        datetime: optional np.ndarray of type np.datetime64 to mark
            all discrete times available for this object.
    """

    def __init__(
        self,
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        time_step: np.timedelta64 = None,
        datetime: np.ndarray = None,
        verbosity: int = 0,
    ):
        "Initialize time with data and parameters."
        super().__init__()
        self.name = "Time"
        self._coords = ["datetime"]
        self._potential_variables = []
        self.verbosity = verbosity

        if (datetime) is not None:
            # JLM check that it's an np.ndarray of type np.datetime64
            if start_time is None:
                self._start_time = datetime[0]
            else:
                self._start_time = start_time
            if end_time is None:
                self._end_time = datetime[-1]
            else:
                self._end_time = end_time

            if time_step is None:
                self._time_step = datetime[1] - datetime[0]
            else:
                self._time_step = time_step

            # Check that all the time deltas in datetime match.
            assert (np.diff(datetime) == self.time_step).all()

            # set this after setting self._start_time and self._time_step
            self["datetime"] = datetime

        else:

            self.datetime = None
            self._start_time = start_time
            self._end_time = end_time
            self._current_time_index = 1
            self._time_step = time_step

        self._current_time = self._start_time
        self._previous_time = None
        self._previous_time_index = None

        # JLM: should probably make the coordinate data private once it is set?
        # JLM: is that using a different setitem option or _set_coord

        # self._current_time_index = start_time

        # self.attrs
        # self.global_attrs

    def __setitem__(self, name: str, value: np.ndarray) -> None:
        super().__setitem__(name, value)
        if name == "datetime":
            wh_start = np.where(self.datetime == self._start_time)[0]
            if len(wh_start) == 0:
                msg = (
                    f"start_time '{self._start_time}' is not in "
                    "supplied datetime array."
                )
                raise ValueError(msg)
            self._current_time_index = wh_start.tolist()[0]
            _ = self.n_time
        return

    @property
    def n_time(self):
        """Get the length of the time dimension."""
        if self.datetime is None:
            return 2
        else:
            len_datetime = len(self.datetime)
            if len_datetime < 2:
                msg = "Time's datetime data must be at least 2 times"
                raise ValueError(msg)
            return len_datetime

    @property
    def current_time(self):
        """Get the current time."""
        return self._current_time

    @property
    def current_time_index(self):
        """Get the current time."""
        return self._current_time_index

    @property
    def previous_time(self):
        """Get the previous time index."""
        return self._previous_time

    @property
    def previous_time_index(self):
        """Get the previous time index."""
        return self._previous_time_index

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

    @staticmethod
    def load(
        control_file: fileish,
        verbosity: int = 0,
    ) -> "Time":
        """Initialize Time object from a PRMS control file

        Args:
            control_file: PRMS control file
            verbosity: output verbosity level

        Returns:
            Time: Time object initialized from a PRMS control file

        """
        control = ControlVariables.load(control_file)
        return Time(
            control.control.start_time,
            control.control.end_time,
            control.control.initial_deltat,
            verbosity=verbosity,
        )

    def advance(self):
        """Advance time."""
        if self.datetime is not None:
            if self._current_time_index + 1 == self.n_time:
                msg = f"End of timeseries reached, can not advance {self.name}"
                raise ValueError(msg)
            self._previous_time_index = self._current_time_index
            self._current_time_index += 1

        if self._previous_time_index is None:
            self._previous_time_index = 0

        self._previous_time = self._current_time
        self._current_time += self.time_step
        return
