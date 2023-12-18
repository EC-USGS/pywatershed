import numpy as np

from ..base.control import Control

# should this be subclassed to DoyArray and TimeseriesArray?


class TimeseriesArray:
    def __init__(
        self,
        control: Control,
        var_name: str,
        array: np.ndarray,
        time: np.ndarray = None,
    ):
        # check that array is > 1d
        # check that first dimension of array is the same length as time array
        self.name = "TimeseriresArray"
        self.control = control
        self.variable_name = var_name
        self.time = time
        self.data = array

        self._current = self.data[0, :].copy()  # copy is necessary

        # check if control.start_time is in time
        # integer is doy type and np.datetime64 is time type
        self._time_type = self.time.dtype

        return

    def advance(self):
        self._set_current()
        return

    def _set_current(self):
        if not self.control._current_time:
            # then data is all nans, just use 0
            time_ind = 0

        else:
            # This is for DOY time which is integer
            # Windows represents this as int32, others as int64
            if self._time_type in [np.int32, np.int64]:
                time_ind = self.control.current_doy - 1

            else:
                # This is for datetime, which is float
                if self.control.itime_step == 0:
                    start_time_ind = np.where(
                        self.time == self.control._start_time
                    )[0]
                    if not len(start_time_ind):
                        msg = "Control start_time is not in the input time"
                        raise ValueError(msg)
                    self._init_time_ind = start_time_ind[0]

                time_ind = self._init_time_ind + self.control.itime_step

        self._current[:] = self.data[time_ind, :]
        return

    @property
    def current(self):
        # return self.get_current()
        return self._current
