import datetime
import numpy as np

# JLM these np time manipulations need to be tested b/c i fear np might change
# These should be built into numpy but it is stalled: https://github.com/numpy/numpy/pull/14276


def dt64_to_dt(dt64: np.datetime64) -> datetime.datetime:
    """np.datetime64 to datetime.datetime, only works for scalars"""
    # This is because I always forget this exists. And to take a vector
    # solution at some point.
    return dt64.astype(datetime.datetime)


def datetime_year(dt64: np.datetime64) -> int:
    """Get the year from np.datetime64"""
    return dt64.astype("datetime64[Y]").astype(int) + 1970


def datetime_month(dt64: np.datetime64) -> int:
    """Get the month from np.datetime64"""
    return dt64.astype("datetime64[M]").astype(int) % 12 + 1


def datetime_doy(dt64: np.datetime64) -> int:
    """Get day of year from np.datetime64"""
    return (dt64 - dt64.astype("datetime64[Y]")).astype(
        "timedelta64[D]"
    ).astype(int) + 1


def datetime_dowy(dt64: np.datetime64) -> int:
    """Get day of water year from np.datetime64"""
    year_start = datetime_year(dt64)
    if datetime_month(dt64) < 10:
        year_start -= 1
    diff = dt64 - np.datetime64(f"{year_start}-10-01")
    return diff.astype("timedelta64[D]").astype(int)
