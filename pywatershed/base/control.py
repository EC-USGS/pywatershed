"""The control class."""
import datetime
import pathlib as pl

import numpy as np

from ..base import meta
from ..constants import fileish
from ..utils import ControlVariables
from ..utils.path import assert_exists, path_rel_to_yml
from ..utils.time_utils import (
    datetime_dowy,
    datetime_doy,
    datetime_epiweek,
    datetime_month,
    datetime_year,
)
from .accessor import Accessor


class Control(Accessor):
    def __init__(
        self,
        start_time: np.datetime64,
        end_time: np.datetime64,
        time_step: np.timedelta64,
        init_time: np.datetime64 = None,
        config: dict = None,
        verbosity: int = 0,
        **kwargs,
    ):
        """Initialize the control class

        The Control class manages input passed at run time, metadata, and
        keeps track of time.

        Args:
            start_time: this is the first time of integration NOT the restart
                time
            end_time: the last integration time
            time_step: the length fo the time step
            config: a PRMS config file to read and use for control parameters
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
        self.meta = meta
        # This will have the time dimension name
        # This will have the time coordimate name

    @classmethod
    def load(
        cls,
        control_file: fileish,
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

    def edit_end_time(self, new_end_time: np.datetime64):
        "Supply a new end time for the simulation."

        self._end_time = new_end_time
        assert self._end_time - self._start_time > 0
        self._n_times = (
            (self._end_time - self._start_time) / self._time_step
        ) + 1
        return

    def edit_n_time_steps(self, new_n_time_steps: int):
        "Supply a new number of timesteps to change the simulation end time."
        self._n_times = new_n_time_steps
        self._end_time = (
            self._start_time + (self._n_times - 1) * self._time_step
        )
        return

    @staticmethod
    def from_yml(yml_file):
        """Instantate a Control object from a yml file

        Required key:value pairs:
            start_time: ISO8601 string for numpy datetime64,
                e.g. 1979-01-01T00:00:00
            end_time: ISO8601 string for numpy datetime64,
                e.g. 1980-12-31T00:00:00
            time_step: The first argument to get a numpy.timedelta64, e.g. 24
            time_step_units: The second argument to get a numpy.timedelta64,
                e.g. 'h'
            verbosity: integer 0-10
            input_dir: path relative to this file.
            budget_type: None | "warn" | "error"
            calc_method: None | "numpy" | "numba" | "fortran" (
                depending on availability)
            load_n_time_batches: integer < total number of timesteps (
                optionalize?)
            init_vars_from_file: False (optionalize, rename restart)
            dprst_flag: True (optionalize) only for PRMSSoilzone (should
                be supplied in its process dictionary)

        Optional key, value pairs:
            netcdf_output: boolean
            netcdf_output_var_names: list of variable names to output, e.g.
                - albedo
                - cap_infil_tot
                - contrib_fraction

        Returns:
            Control object
        """
        import yaml

        with pl.Path(yml_file).open("r") as file_stream:
            control_dict = yaml.load(file_stream, Loader=yaml.Loader)

        start_time = np.datetime64(control_dict["start_time"])
        end_time = np.datetime64(control_dict["end_time"])
        time_step = np.timedelta64(
            control_dict["time_step"], control_dict["time_step_units"]
        )
        verbosity = control_dict["verbosity"]

        paths_to_convert = ["input_dir"]
        for path_name in paths_to_convert:
            control_dict[path_name] = path_rel_to_yml(
                control_dict[path_name], yml_file
            )
            assert_exists(control_dict[path_name])

        control = Control(
            start_time,
            end_time,
            time_step,
            config=control_dict,
            verbosity=verbosity,
        )
        return control
