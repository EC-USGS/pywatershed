import numpy as np

# JLM this seems treating the symptom
try:
    import pandas as pd
    from pandas import DataFrame
except ModuleNotFoundError:
    pd = None
    DataFrame = None


# JLM: where do we keep the metadata attributes?

# This is a base AtmBoundaryLayer class
# It has no state.
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

        return

    def _has_state(self, state_name, fail=True):
        if hasattr(self, state_name):
            return True
        else:
            if fail:
                msg = f"{self.name} does not have state '{state_name}'"
                raise ValueError(msg)
            else:
                return False

    # do we really want set_state public??
    def set_state(self, state_name: str, value: np.array) -> None:
        if self._has_state(state_name):
            setattr(self, state_name, value)
        return None

    # give a time argument for getting the state at some/current time?
    def get_state(self, state_name: str) -> np.array:
        if self._has_state(state_name):
            return getattr(self, state_name)

    def get_current_state(self, state_name: str) -> np.array:
        if self._has_state(state_name):
            return getattr(self, state_name)[self.current_time_index, :]

    @property
    def current_time(self) -> np.datetime64:
        return self.time[self.current_time_index]

    @property
    def time_step(self) -> np.timedelta64:
        return self._time_step

    # JLM check that time and the time dimensions on the states match on set or inits
    # JLM ensure the spatial dimension matches.

    def advance(
        self,
        itime_step: int,
        current_date: np.datetime64,
    ):
        self.itime_step = itime_step
        self.current_date = current_date
        dt = pd.to_datetime(current_date)
        for key, df in self.forcings.items():
            self.current[key] = df.iloc[df.index.get_loc(dt, method="nearest")]
