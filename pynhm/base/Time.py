import numpy as np

from .DataAccess import DataAccess


class Time(DataAccess):
    """The time object base class.

    Objects can access and manage time the same way.

    Dimension:
        'time'
    Coordinate data variable:
        'datetime', a one-dimensional
        np.ndarray of type np.datetime64
        The length of the datetime data is the number of times
        to maintain in state. (This could be changed in the future
        with an optional argument).

    Attributes:
        verbosity: The verbosity level [0-10].
    """

    def __init__(
        self,
        start_time: np.datetime64 = None,
        time_step: np.timedelta64 = None,
        datetime: np.ndarray = None,
        verbosity: int = 0,
    ):
        "Initialize time with data and parameters."
        super().__init__()
        self.name = "Time"
        self._coords = ["datetime"]
        self._variables = []
        self._potential_variables = []
        self.verbosity = verbosity

        if (datetime) is not None:
            if len(datetime) < 2:
                msg = "Time's datetime data must be at least 2 times"
                raise ValueError(msg)

            # JLM check that it's an np.ndarray of type np.datetime64
            self["datetime"] = datetime

            if start_time is None:
                self._start_time = self["datetime"][0]
            else:
                self._start_time = start_time

            if time_step is None:
                self._time_step = self["datetime"][1] - self["datetime"][0]
            else:
                self._time_step = time_step

            # Check that all the time deltas in datetime match.
            assert (np.diff(self.datetime) == self.time_step).all()

            wh_start = np.where(self.datetime == self._start_time)[0]
            if len(wh_start) == 0:
                msg = (
                    f"start_time '{start_time}' is not in "
                    "supplied datetime array."
                )
                raise ValueError(msg)
            self._current_time_index = wh_start.tolist()[0]

        else:

            self.datetime = None
            self._start_time = start_time
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

    def advance(self):
        """Advance time."""
        self._previous_time = self._current_time
        self._current_time += self.time_step
        if self.datetime is not None:
            self._previous_time_index = self._current_time_index
            self._current_time_index += 1
        if self._previous_time_index is None:
            self._previous_time_index = 0
        return
