from copy import deepcopy

import numpy as np

# This is a base AtmBoundaryLayer class
# It has no state, only time.

# JLM: where do we keep the metadata attributes?

# JLM: The state is concieved here as a timeseries, not just an instant or
# JLM: markov (previous and current).

# JLM: Are the forcings valid over the interval [current_time, current_time + time_step]?
# JLM: This should go in the attributes

# JLM: state or FLUX? differentiate? just in metadata? change nomenclature?


class AtmBoundaryLayer:
    def __init__(
        self,
        start_time: np.datetime64,
        time_step: np.timedelta64,  # could we infer timestep?
        height_m: int = None,
        verbose: int = 0,
    ):
        self.name = "AtmBoundaryLayer"
        self.coords = ["datetime", "spatial_id"]
        self.variables = []
        self._potential_variables = deepcopy(self.variables)
        self.start_time = start_time
        self._current_time = start_time
        self._time_step = time_step
        self.height_m = height_m
        self.verbose = verbose

        # method for initing this... need time to get it
        self._current_time_index = 0

        self.datetime = None
        # spatial dimension?
        self.spatial_id = None

        return

    # The set/get methods seem like they could be inherited from
    # a base state class.
    def _has_potential_state(self, state_name, fail=True) -> bool:
        if state_name in self._potential_variables:
            return True
        return False

    def _has_state(self, state_name, fail=True) -> bool:
        if state_name in self.variables + self.coords:
            return True
        return False

    def _has_coord(self, coord_name, fail=True) -> bool:
        if coord_name in self.coords:
            return True
        return False

    def __setitem__(self, name: str, value: np.ndarray) -> None:
        if not isinstance(value, np.ndarray):
            raise ValueError(
                f"Value to set must be an np.ndarray not {type(value)}"
            )
        if self._has_potential_state(name):
            setattr(self, name, value)
            self.variables = list(set(self.variables + [name]))
        elif self._has_coord(name):
            setattr(self, name, value)
        else:
            msg = (
                f"{self.name} does not have potential state or coordinate "
                f"'{name}'"
            )
            raise KeyError(msg)
        return None

    def __getitem__(self, name: str) -> np.ndarray:
        if self._has_state(name) or self._has_coord(name):
            return getattr(self, name)
        else:
            msg = f"'{name}' not in state or coordinate variables"
            raise KeyError(msg)

    def get_current_state(self, state_name: str) -> np.ndarray:
        return self[state_name].take(indices=self.current_time_index, axis=0)

    def __delitem__(self, name: str):
        if self._has_state(name):
            return getattr(self, name)
        elif self._has_coord(name):
            msg = f"Can not delete coordinate '{name}'"
            raise KeyError(msg)
        else:
            msg = f"'{name}' not in state or coordinate variables"
            raise KeyError(msg)
        return None

    # Time tracking could also be handled by state base class.
    # JLM all this stuff might go in a time or datetime subdir and be common to
    # all objects that need time?
    @property
    def current_time_index(self) -> np.datetime64:
        if self.datetime is None:
            return None
        else:
            # JLM this is going to be slow, there are ways to ensure
            # continuity with more speed.
            ind = np.where(self.datetime == self._current_time)[0]
            if len(ind) == 0:
                msg = (
                    f"Current/new datetime is not in the time data: "
                    f"{self._current_time}"
                )
                raise ValueError(msg)
            return ind

    @property
    def current_time(self) -> np.datetime64:
        if self.datetime is None:
            return None
        else:
            return self.datetime[self.current_time_index]

    @property
    def time_step(self) -> np.timedelta64:
        return self._time_step

    # JLM check that time and the time dimensions on the states match on set or inits
    # JLM ensure the spatial dimension matches.

    def advance(self) -> None:
        self._current_time += self._time_step
        # this ensures sanity
        _ = self.current_time
        return

    def calculate(self) -> None:
        pass
