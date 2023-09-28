import datetime
import pathlib as pl
from warnings import warn

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

# This is the list of control variables currently used by pywatershed
# It is important to maintain this list to issue warnings about what
# variables are unrecognized/ignored in legacy and non-legacy control
# files
# TODO: where should these be documented?
# TODO: identify which are PRMS-legacy?
pws_control_options_avail = [
    "budget_type",
    "calc_method",
    "dprst_flag",  # to remove?
    "restart",
    "input_dir",
    "load_n_time_batches",
    "netcdf_output_dir",
    "netcdf_output_var_names",
    # "netcdf_output_separate_files",
    # "netcdf_budget_args",
    "start_time",
    "time_step_units",
    "verbosity",
]

prms_legacy_options_avail = [
    "dprst_flag",
    "end_time",
    "init_vars_from_file",
    "initial_deltat",
    "nhruOutBaseFileName",
    "nhruOutVar_names",
    "nsegmentOutBaseFileName",
    "nsegmentOutVar_names",
    "start_time",
    "print_debug",
]

prms_to_pws_option_map = {
    "init_vars_from_file": "restart",
    "initial_deltat": "time_step",
    "nhruOutBaseFileName": "netcdf_output_dir",
    "nhruOutVar_names": "netcdf_output_var_names",
    "nsegmentOutBaseFileName": "netcdf_output_dir",
    "nsegmentOutVar_names": "netcdf_output_var_names",
    "print_debug": "verbosity",
}

assert (
    len(set(prms_to_pws_option_map.keys()) - set(prms_legacy_options_avail))
    == 0
)


class Control(Accessor):
    """Control manages global time and options, and provides metadata.

    Args:
        start_time: this is the first time of integration NOT the restart
            time
        end_time: the last integration time
        time_step: the length fo the time step
        options: a dictionary of global Process options.

    """

    def __init__(
        self,
        start_time: np.datetime64,
        end_time: np.datetime64,
        time_step: np.timedelta64,
        init_time: np.datetime64 = None,
        options: dict = None,
    ):
        super().__init__()
        self.name = "Control"

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

        if options is None:
            options = {}
        self.options = {}
        self._set_options(options)
        self.meta = meta
        # This will have the time dimension name
        # This will have the time coordimate name

    @classmethod
    def load(
        cls,
        control_file: fileish,
    ) -> "Control":
        msg = "Control.load will be deprecated for Control.load_prms"
        warn(msg, PendingDeprecationWarning)
        return Control.load_prms(control_file)

    @classmethod
    def load_prms(
        cls,
        control_file: fileish,
        warn_unused_options: bool = True,
    ) -> "Control":
        """Initialize a control object from a PRMS control file

        Args:
            control_file: PRMS control file
            warn_unused_options: bool if warnings are to be issued for unused
                options from the PRMS control file. Recommended and True by
                default. See below for a list of used/available legacy options.

        Returns:
            Time: Time object initialized from a PRMS control file



        Available PRMS legacy options :
            nhruOutVar_names: mapped to netcdf_output_var_names
            nsegmentOutVar_names: mapped to netcdf_output_var_names


        """
        control = ControlVariables.load(control_file)

        if warn_unused_options:
            for vv in control.control.keys():
                if vv not in prms_legacy_options_avail:
                    msg = (
                        f"Option '{vv}' in supplied control file is not used "
                        "by pywatershed"
                    )
                    warn(msg, RuntimeWarning)

        opts = control.control
        opt_names = list(opts.keys())

        for oo in opt_names:
            if oo not in prms_legacy_options_avail:
                del opts[oo]
            if oo in prms_to_pws_option_map.keys():
                pws_option_key = prms_to_pws_option_map[oo]
                val = opts[oo]
                del opts[oo]
                if pws_option_key in opts.keys():
                    # combine to a list with only unique entries
                    # use value instead of list if only one value in list
                    opts[pws_option_key] = list(
                        set(opts[pws_option_key].tolist() + val.tolist())
                    )
                    if len(opts[pws_option_key]) == 1:
                        opts[pws_option_key] = opts[pws_option_key][0]
                else:
                    opts[pws_option_key] = val

        start_time = control.control["start_time"]
        end_time = control.control["end_time"]
        time_step = control.control["time_step"]
        del control.control["start_time"]
        del control.control["end_time"]
        del control.control["time_step"]

        return cls(
            start_time=start_time,
            end_time=end_time,
            time_step=time_step,
            options=control.control,
        )

    def _set_options(self, options):
        for okey, oval in options.items():
            if okey not in pws_control_options_avail:
                msg = f"'{okey}' is not an available control option"
                raise ValueError(msg)
            self.options[okey] = oval

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
    def current_month(self, zero_based: bool = False):
        """Get the current month in 1-12 (unless zero based)."""
        return datetime_month(self._current_time, zero_based=zero_based)

    @property
    def current_doy(self, zero_based: bool = False):
        """Get the current day of year in 1-366 (unless zero based)."""
        return datetime_doy(self._current_time, zero_based=zero_based)

    @property
    def current_dowy(self, zero_based: bool = False):
        """Get the current day of water year in 1-366 (unless zero-based)."""
        return datetime_dowy(self._current_time, zero_based=zero_based)

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
            int((self._end_time - self._start_time) / self._time_step) + 1
        )
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
        del control_dict["start_time"]
        del control_dict["end_time"]
        del control_dict["time_step"]
        del control_dict["time_step_units"]

        paths_to_convert = ["input_dir"]
        for path_name in paths_to_convert:
            if path_name in control_dict.keys():
                control_dict[path_name] = path_rel_to_yml(
                    control_dict[path_name], yml_file
                )
                assert_exists(control_dict[path_name])

        control = Control(
            start_time,
            end_time,
            time_step,
            options=control_dict,
        )
        return control
