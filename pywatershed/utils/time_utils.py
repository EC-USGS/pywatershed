import datetime

import epiweeks as ew
import numpy as np

# These should be built into numpy but it is stalled:
# https://github.com/numpy/numpy/pull/14276


def dt64_to_dt(dt64: np.datetime64, dt_precision="s") -> datetime.datetime:
    """np.datetime64 to datetime.datetime"""
    # This is because I always forget this exists.
    val = dt64.astype(datetime.datetime)
    # If the precision of dt64 is such that it is an integer too large to be
    # coerced to datetime.datetime, then it just returns an integer and
    # not a datetime.datetime object, which breaks other stuff downstream.
    # default to 's' precision in case where this happens.
    if isinstance(val, int):
        dt = dt64.astype(f"datetime64[{dt_precision}]")
        val = dt.astype(datetime.datetime)
    # <
    return val


def datetime_year(dt64: np.datetime64) -> int:
    """Get the year from np.datetime64"""
    return dt64.astype("datetime64[Y]").astype(int) + 1970


def datetime_month(dt64: np.datetime64, zero_based=False) -> int:
    """Get the month from np.datetime64
    Args:
         zero_based: count as a zero-based index.
    """
    return dt64.astype("datetime64[M]").astype(int) % 12 + _offset(zero_based)


def datetime_day_of_month(dt64: np.datetime64, zero_based=False) -> int:
    """Get the month from np.datetime64
    Args:
         zero_based: count as a zero-based index.
    """
    dfmt = "datetime64[D]"
    mfmt = "datetime64[M]"
    diff = (dt64.astype(dfmt) - dt64.astype(mfmt).astype(dfmt)).astype(int)
    return diff + _offset(zero_based)


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
