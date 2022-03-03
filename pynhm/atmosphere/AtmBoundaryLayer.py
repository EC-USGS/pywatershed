import numpy as np
import pandas as pd

# JLM: where do we keep the metadata attributes?

# This is a base AtmBoundaryLayer class
# It has no state, only time.
# How do we define the start time?
# Are the forcings valid over the interval [current_time, current_time + time_step]?


class AtmBoundaryLayer:
    def __init__(
        self,
        start_time: np.datetime64,
        time_step: np.timedelta64,  # could we infer timestep?
        height_m: int = None,
        verbose: int = 0,
    ):
        self.name = "AtmBoundaryLayer"
        self.start_time = start_time
        self._current_time = start_time
        self._time_step = time_step
        self.height_m = height_m
        self.verbose = verbose

        # method for initing this... need time to get it
        self._current_time_index = 0

        self.datetime = None

        return

    def _has_state(self, state_name, fail=True) -> bool:
        if hasattr(self, state_name):
            return True
        else:
            if fail:
                msg = f"{self.name} does not have state '{state_name}'"
                raise ValueError(msg)
            else:
                return False

    # JLM: do want set_state public? or part of init?
    def set_state(self, state_name: str, value: np.array) -> None:
        if self._has_state(state_name):
            setattr(self, state_name, value)
        return None

    def get_state(self, state_name: str) -> np.array:
        if self._has_state(state_name):
            return getattr(self, state_name)

    # JLM: Getattr is nice but i assume it is copying.
    def get_current_state(self, state_name: str) -> np.array:
        if self._has_state(state_name):
            return getattr(self, state_name).take(
                indices=self.current_time_index, axis=0
            )

    # JLM all this stuff might go in a time or datetime subdir and be common to
    # all objects that need time?
    @property
    def current_time_index(self) -> np.datetime64:
        if self.datetime is None:
            return None
        else:
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
