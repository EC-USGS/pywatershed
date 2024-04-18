from datetime import datetime
import math
import random

import epiweeks as ew
import numpy as np
from pywatershed.utils.time_utils import (
    dt64_to_dt,
    datetime_year,
    datetime_month,
    datetime_day_of_month,
    datetime_doy,
    datetime_dowy,
)

# should probably test the offset too


def random_datetime_datetime64():
    start = datetime.fromisoformat("1980-01-01")
    end = datetime.fromisoformat("2100-01-01")
    dt = start + (end - start) * random.uniform(0, 1)
    dt64 = np.datetime64(str(dt))
    return (dt, dt64)


def test_dt64_to_dt():
    dt, dt64 = random_datetime_datetime64()
    assert dt == dt64_to_dt(dt64)


def test_datetime_year():
    dt, dt64 = random_datetime_datetime64()
    assert dt.year == datetime_year(dt64)


def test_datetime_month():
    dt, dt64 = random_datetime_datetime64()
    assert dt.month == datetime_month(dt64)


def test_datetime_day_of_month():
    dt, dt64 = random_datetime_datetime64()
    assert dt.day == datetime_day_of_month(dt64)


def test_datetime_doy():
    dt, dt64 = random_datetime_datetime64()
    assert dt.timetuple().tm_yday == datetime_doy(dt64)


def test_datetime_dowy():
    dt, dt64 = random_datetime_datetime64()
    year = dt.year
    month = dt.month
    if month < 10:
        year -= 1
    diff = dt - datetime.fromisoformat(f"{year}-10-01")
    diff = math.ceil(diff.total_seconds() / 60 / 60 / 24)
    assert diff == datetime_dowy(dt64)


def test_epiweek():
    dt, dt64 = random_datetime_datetime64()
    assert ew.Week.fromdate(dt) == ew.Week.fromdate(dt64_to_dt(dt64))
