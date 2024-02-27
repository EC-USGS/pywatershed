import datetime

import epiweeks as ew
import numpy as np

# JLM these np time manipulations need to be tested b/c i fear np might change
# These should be built into numpy but it is stalled:
# https://github.com/numpy/numpy/pull/14276


def dt64_to_dt(dt64: np.datetime64) -> datetime.datetime:
    """np.datetime64 to datetime.datetime"""
    # This is because I always forget this exists.
    return dt64.astype(datetime.datetime)


def datetime_year(dt64: np.datetime64) -> int:
    """Get the year from np.datetime64"""
    return dt64.astype("datetime64[Y]").astype(int) + 1970


def datetime_month(dt64: np.datetime64, zero_based=False) -> int:
    """Get the month from np.datetime64
    Args:
         zero_based: count as a zero-based index.
    """
    return dt64.astype("datetime64[M]").astype(int) % 12 + _offset(zero_based)


def datetime_doy(dt64: np.datetime64, zero_based=False) -> int:
    """Get day of year from np.datetime64
    Args:
        zero_based: count as a zero-based index.
    """
    return (dt64 - dt64.astype("datetime64[Y]")).astype(
        "timedelta64[D]"
    ).astype(int) + _offset(zero_based)


def datetime_dowy(dt64: np.datetime64, zero_based=False) -> int:
    """Get day of water year from np.datetime64
    Args:
        zero_based: count as a zero-based index.
    """
    year_start = datetime_year(dt64)
    if datetime_month(dt64) < 10:
        year_start -= 1
    diff = dt64 - np.datetime64(f"{year_start}-10-01")
    return diff.astype("timedelta64[D]").astype(int) + _offset(zero_based)


def datetime_epiweek(dt64: np.datetime64) -> int:
    """Get CDC eipweek [1, 53] from np.datetime64"""
    return ew.Week.fromdate(dt64_to_dt(dt64)).week


def _offset(zero_based: bool):
    if zero_based:
        return 0
    else:
        return 1
