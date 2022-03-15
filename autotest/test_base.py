from datetime import datetime, timedelta

import numpy as np
import pytest

from pynhm.base.DataAccess import DataAccess
from pynhm.base.Time import Time


class TestDataAccess:
    def test_init(self):
        da = DataAccess()
        assert da.name == "DataAccess"
        assert len(da.coords) == 0
        assert len(da.variables) == 0
        assert len(da._potential_variables) == 0
        return

    @pytest.mark.parametrize(
        "data", [np.arange(4), list(range(4))], ids=["valid", "invalid"]
    )
    def test_coord(self, data):
        da = DataAccess()
        da._coords = ["foo"]  # strictly verboten!
        try:
            da["foo"] = data
            assert isinstance(data, np.ndarray)
            assert np.isclose(da["foo"], data).all()
            # can not delete a coord
        except TypeError:
            assert not isinstance(data, np.ndarray)
            return

        # can not set or delete coords directly
        try:
            da.coords = ["foo"]
            assert False
        except AttributeError:
            assert True

        try:
            da.coords["foo"] = data
            assert False
        except TypeError:
            assert True

        # can not delete a coord
        try:
            del da["foo"]
            assert False
        except KeyError:
            assert True

        return

    @pytest.mark.parametrize(
        "data", [np.arange(4), list(range(4))], ids=["valid", "invalid"]
    )
    def test_variable(self, data):
        da = DataAccess()
        da._potential_variables = ["foo"]  # strictly verboten!
        try:
            da["foo"] = data
            assert "foo" in da.variables
            assert isinstance(data, np.ndarray)
            assert np.isclose(da["foo"], data).all()
            del da["foo"]
            assert len(da.variables) == 0
        except TypeError:
            assert not isinstance(data, np.ndarray)
            return

        # can not set or delete coords directly
        try:
            da.variables = ["foo"]
            assert False
        except AttributeError:
            assert True

        return


time_data = np.arange(
    datetime(1979, 1, 1), datetime(1979, 1, 5), timedelta(days=1)
).astype(np.datetime64)

start_times = [time_data[0], time_data[3], np.datetime64(datetime(1980, 1, 1))]
time_step = np.timedelta64(24, "h")


class TestTime:
    @pytest.mark.parametrize(
        "start_time", start_times, ids=["valid0", "valid1", "invalid"]
    )
    def test_init_markov(self, start_time):
        time = Time(start_time=start_time, time_step=time_step)
        assert time.current_time == start_time
        assert time.current_time_index == 1
        assert time.previous_time is None
        assert time.previous_time_index is None
        assert time.time_step == time_step

        time.advance()
        assert time.current_time_index == 1
        assert time.current_time == start_time + time_step
        assert time.previous_time == start_time
        assert time.previous_time_index == 0
        assert time.time_step == time_step

        # fail re-setting datetime
        try:
            time["datetime"] = np.array(np.datetime64(start_time))
            assert False
        except KeyError:
            assert True

        return

    @pytest.mark.parametrize(
        "start_time", start_times, ids=["valid0", "valid1", "invalid"]
    )
    def test_init_timeseries(self, start_time):
        try:
            time = Time(
                start_time=start_time, time_step=time_step, datetime=time_data
            )
            # make sure we have the right case, as in the except
            wh_start = np.where(time_data == start_time)[0]
            assert len(wh_start) > 0
            assert time.current_time_index == wh_start
            assert time.current_time == start_time
            assert time.previous_time is None
            assert time.previous_time_index is None
            assert time.time_step == time_step
        except ValueError:
            assert len(np.where(time_data == start_time)[0]) == 0
            return

        time.advance()
        assert time.current_time_index == wh_start + 1
        assert time.current_time == start_time + time_step
        assert time.previous_time == start_time
        assert time.previous_time_index == wh_start
        assert time.time_step == time_step

        # fail re-setting datetime
        try:
            time["datetime"] += np.timedelta64(24, "h")
            assert False
        except KeyError:
            assert True

        return
