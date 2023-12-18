import datetime
import pathlib as pl
from collections import UserDict
from copy import deepcopy
from typing import Union
from warnings import warn

import numpy as np
import yaml

from ..base import meta
from ..constants import fileish
from ..utils import ControlVariables
from ..utils.path import assert_exists, dict_pl_to_str, path_rel_to_yaml
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
# The following are duplicated in the Control docstring below and that
# docstring needs updated whenever any of these change.
pws_control_options_avail = [
    "budget_type",
    "calc_method",
    # "restart",
    "input_dir",
    # "load_n_time_batches",
    "netcdf_output_dir",
    "netcdf_output_var_names",
    "netcdf_output_separate_files",
    "netcdf_budget_args",
    "start_time",
    "time_step_units",
    "verbosity",
]

prms_legacy_options_avail = [
    "end_time",
    # "init_vars_from_file",
    "initial_deltat",
    "nhruOutBaseFileName",
    "nhruOutVar_names",
    "nsegmentOutBaseFileName",
    "nsegmentOutVar_names",
    "start_time",
    "print_debug",
]

prms_to_pws_option_map = {
    # "init_vars_from_file": "restart",
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

    Available pywatershed options:
      * budget_type: one of [None, "warn", "error"]
      * calc_method: one of ["numpy", "numba", "fortran"]
      * input_dir: str or pathlib.path directory to search for input data
      * netcdf_output_dir: str or pathlib.Path directory for output
      * netcdf_output_var_names: a list of variable names to output
      * netcdf_output_separate_files: bool if output is grouped by Process or
        if each variable is written to an individual file
      * netcdf_budget_args:
      * start_time: np.datetime64
      * end_time: np.datetime64
      * time_step_units: str containing single character code for
        np.timedelta64
      * verbosity: 0-10

    Available PRMS legacy options:
      Either used as-is or mapped to pywatershed options as indicated below.

      * start_time
      * end_time
      * initial_deltat: translates to "time_step"
      * init_vars_from_file: translates to "restart"
      * nhruOutBaseFileName: translates to "netcdf_output_dir"
      * nhruOutVar_names: translates to a subset of "netcdf_output_var_names"
      * nsegmentOutBaseFileName: translates to "netcdf_output_dir"
      * nsegmentOutVar_names: translates to a subset of
        "netcdf_output_var_names"
      * print_debug: translates to "verbosity"


    Examples:
    ---------

    >>> import pathlib as pl
    >>>
    >>> import numpy as np
    >>> import pywatershed as pws
    >>>
    >>> control = pws.Control(
    ...     start_time=np.datetime64("2023-01-01T00:00:00"),
    ...     end_time=np.datetime64("2023-01-02T00:00:00"),
    ...     time_step=np.timedelta64(24, "h"),
    ...     options={"input_dir": pl.Path("./input")},
    ... )
    >>>
    >>> # A more interesting example reads from a PRMS control file
    >>> pws_root = pws.constants.__pywatershed_root__
    >>> drb_control_file = pws_root / "data/drb_2yr/control.test"
    >>> control_drb = pws.Control.load_prms(
    ...     drb_control_file,
    ...     warn_unused_options=False,
    ... )
    >>> control_drb.current_time
    numpy.datetime64('1978-12-31T00:00:00')
    >>> control_drb.previous_time
    >>> control_drb.init_time
    numpy.datetime64('1978-12-31T00:00:00')
    >>> control_drb.time_step
    numpy.timedelta64(24,'h')
    >>> control_drb.time_step_seconds
    86400.0
    >>> control_drb.start_time
    numpy.datetime64('1979-01-01T00:00:00')
    >>> control_drb.advance()
    >>> control_drb.current_time
    numpy.datetime64('1979-01-01T00:00:00')
    >>> control_drb.current_doy
    1
    >>> control_drb.current_dowy
    93
    >>> control_drb.previous_time
    numpy.datetime64('1978-12-31T00:00:00')



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

        self._current_time = self._init_time
        self._previous_time = None
        self._itime_step = -1

        if options is None:
            options = OptsDict()
        self.options = options
        self.meta = meta
        # This will have the time dimension name
        # This will have the time coordimate name

    @classmethod
    def load(
        cls,
        control_file: fileish,
        warn_unused_options: bool = True,
    ) -> "Control":
        msg = "Control.load will be deprecated for Control.load_prms"
        warn(msg, PendingDeprecationWarning)
        return Control.load_prms(
            control_file, warn_unused_options=warn_unused_options
        )

    @classmethod
    def load_prms(
        cls,
        control_file: fileish,
        warn_unused_options: bool = True,
    ) -> "Control":
        """Initialize a control object from a PRMS control file.

        Args:
            control_file: PRMS control file
            warn_unused_options: bool if warnings are to be issued for unused
                options from the PRMS control file. Recommended and True by
                default. See below for a list of used/available legacy options.

        Returns:
            An instance of a Control object.
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

    def _set_options(self, options: dict):
        if not isinstance(options, (OptsDict, dict)):
            raise ValueError("control.options must be a dictionary")
        valid_options = OptsDict()
        for key, val in options.items():
            valid_options[key] = val

        return valid_options

    def __setitem__(self, key, value) -> None:
        if key == "options":
            value = self._set_options(value)

        super().__setitem__(key, value)
        return None

    def __setattr__(self, name, value) -> None:
        if name == "options":
            value = self._set_options(value)

        super().__setattr__(name, value)
        return None

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        del self.meta
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))

        self.meta = meta
        result.meta = meta
        return result

    @property
    def current_time(self) -> np.datetime64:
        """Get the current time."""
        return self._current_time

    @property
    def current_datetime(self) -> datetime.datetime:
        """Get the current time as a datetime.datetime object"""
        return self._current_time.astype(datetime.datetime)

    @property
    def current_year(self) -> int:
        """Get the current year."""
        return datetime_year(self._current_time)

    @property
    def current_month(self) -> int:
        """Get the current month in 1-12 (unless zero based)."""
        return datetime_month(self._current_time)

    @property
    def current_doy(self) -> int:
        """Get the current day of year in 1-366 (unless zero based)."""
        return datetime_doy(self._current_time)

    @property
    def current_dowy(self) -> int:
        """Get the current day of water year in 1-366 (unless zero-based)."""
        return datetime_dowy(self._current_time)

    @property
    def current_epiweek(self) -> int:
        """Get the current epiweek [1, 53]."""
        return datetime_epiweek(self._current_time)

    @property
    def previous_time(self) -> np.datetime64:
        """The previous time."""
        return self._previous_time

    @property
    def itime_step(self) -> int:
        """The counth of the current time [0, self.n_times-1]"""
        return self._itime_step

    @property
    def time_step(self) -> int:
        """The time step"""
        return self._time_step

    @property
    def init_time(self) -> np.datetime64:
        """Get the simulation initialization time"""
        return self._init_time

    @property
    def start_time(self) -> np.datetime64:
        """The simulation start time"""
        return self._start_time

    @property
    def start_doy(self) -> int:
        """The simulation start day of year."""
        return datetime_doy(self._start_time)

    @property
    def start_month(self) -> int:
        """The simulation start month."""
        return datetime_month(self._start_time)

    @property
    def end_time(self) -> np.datetime64:
        """The simulation end time."""
        return self._end_time

    @property
    def n_times(self) -> int:
        """The number of time steps."""
        return self._n_times

    @property
    def time_step_seconds(self) -> np.float64:
        """The timestep length in units of seconds."""
        return self.time_step / np.timedelta64(1, "s")

    def advance(self) -> None:
        """Advance time."""
        if self._current_time == self._end_time:
            raise ValueError("End of time reached")

        self._itime_step += 1

        if self._current_time is None:
            self._current_time = self._start_time
            return

        self._previous_time = self._current_time
        self._current_time += self.time_step

        return None

    def edit_end_time(self, new_end_time: np.datetime64) -> None:
        """Supply a new end time for the simulation.

        Args:
            new_end_time: the new time at which to end the simulation.
        """

        self._end_time = new_end_time
        assert self._end_time - self._start_time > 0
        self._n_times = (
            int((self._end_time - self._start_time) / self._time_step) + 1
        )
        return None

    def edit_n_time_steps(self, new_n_time_steps: int) -> None:
        """Supply a new number of timesteps to change the simulation end time.

        Args:
            new_n_time_steps: The new number of timesteps.
        """
        self._n_times = new_n_time_steps
        self._end_time = (
            self._start_time + (self._n_times - 1) * self._time_step
        )
        return None

    def __str__(self):
        from pprint import pformat

        return pformat(self.to_dict())

    def __repr__(self):
        # TODO: this is not really an object representation
        return self.__str__()

    def to_dict(self, deep_copy=True) -> dict:
        """Export a control object to a dictionary

        Args:
            deep_copy: If the dictionary should be a deep copy or not.
        """

        control_dict = {}

        # I suppose this list could grow with time but these are
        # the only non .option items in __dict__ required to reconstitute a
        # Control instance
        control_dict["start_time"] = str(self.start_time)
        control_dict["end_time"] = str(self.end_time)
        control_dict["time_step"] = str(self.time_step)[0:2]
        control_dict["time_step_units"] = str(self.time_step)[3:4]

        if deep_copy:
            control = deepcopy(self)
        else:
            control = self

        control_dict["options"] = {}
        for kk, vv in control.options.items():
            control_dict["options"][kk] = control.options[kk]

        return control_dict

    def to_yaml(self, yaml_file: Union[pl.Path, str]) -> None:
        """Export to a yaml file

        Args:
            yaml_file: The file to write to.

        Note: This flattens .options to the top level of the yaml/dict
            so that option keys are all at the same level as "start_time",
            "end_time", "time_step", and "time_step_units". Using .from_yaml
            will restore options to a nested dictionary.

        Args:
            yaml_file: pl.Path or str to designate the output path/file.
        """
        control_dict = dict_pl_to_str(self.to_dict())
        opts = control_dict["options"]
        for kk, vv in opts.items():
            if kk in control_dict.keys():
                msg = "Control option keys collide with non-option keys"
                raise ValueError(msg)
            control_dict[kk] = vv

        del control_dict["options"]

        yaml_file = pl.Path(yaml_file)
        with open(yaml_file, "w") as file:
            _ = yaml.dump(control_dict, file)

        assert yaml_file.exists()
        return None

    @staticmethod
    def from_yaml(yaml_file: Union[str, pl.Path]) -> "Control":
        """Instantate a Control object from a yaml file

        Args:
            yaml_file: a yaml file to parse.


        Required key:value pairs:
            * budget_type: None | "warn" | "error"
            * calc_method: None | "numpy" | "numba" | "fortran" (depending on
              availability)
            * end_time: ISO8601 string for numpy datetime64, e.g.
              1980-12-31T00:00:00.
            * input_dir: path relative to this file.
            * start_time: ISO8601 string for numpy datetime64, e.g.
              1979-01-01T00:00:00.
            * time_step: The first argument to get a numpy.timedelta64, e.g. 24
            * time_step_units: The second argument to get a numpy.timedelta64,
              e.g. 'h'.
            * verbosity: integer 0-10

        Optional key, value pairs:
            * netcdf_output: boolean
            * netcdf_output_var_names: list of variable names to output, e.g.
              - albedo
              - cap_infil_tot
              - contrib_fraction

        Returns:
            Control object
        """
        import yaml

        with pl.Path(yaml_file).open("r") as file_stream:
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
                control_dict[path_name] = path_rel_to_yaml(
                    control_dict[path_name], yaml_file
                )
                assert_exists(control_dict[path_name])

        control = Control(
            start_time,
            end_time,
            time_step,
            options=control_dict,
        )
        return control


class OptsDict(UserDict):
    def __setitem__(self, key, value):
        if key not in pws_control_options_avail:
            msg = f"'{key}' is not an available control option"
            raise NameError(msg)
        super().__setitem__(key, value)
        return None
