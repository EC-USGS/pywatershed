import numpy as np

# JLM these np time manipulations need to be tested b/c i fear np might change
# These should be built into numpy but it is stalled: https://github.com/numpy/numpy/pull/14276


def datetime_year(datetime: np.datetime64):
    """Get the year from np.datetime64"""
    return datetime.astype("datetime64[Y]").astype(int) + 1970


def datetime_month(datetime: np.datetime64):
    """Get the month from np.datetime64"""
    return datetime.astype("datetime64[M]").astype(int) % 12 + 1


def datetime_doy(datetime: np.datetime64):
    """Get day of year from np.datetime64"""
    return (datetime - datetime.astype("datetime64[Y]")).astype(
        "timedelta64[D]"
    ).astype(int) + 1


def datetime_dowy(datetime: np.datetime64):
    """Get day of water year from np.datetime64"""
    year_start = datetime_year(datetime)
    if datetime_month(datetime) < 10:
        year_start -= 1
    diff = datetime - np.datetime64(f"{year_start}-10-01")
    return diff.astype("timedelta64[D]").astype(int)
