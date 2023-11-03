import os
import pathlib as pl
from enum import Enum
from typing import Union

import numpy as np

# Environment variables
numba_num_threads = os.getenv("NUMBA_NUM_THREADS")
if numba_num_threads is None:
    numba_num_threads = 0
else:
    numba_num_threads = int(numba_num_threads)

# Typing constants
fileish = Union[str, pl.Path]
listish = Union[str, list, tuple]  # Todo deprecate to Typing.Iterable

# PRMS6 Constants module:
# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/misc/m_constants.f90

__pywatershed_root__ = pl.Path(__file__).parent

zero = np.zeros([1])[0]
one = np.ones([1])[0]
nan = np.nan

epsilon = np.finfo(zero).eps
# https://en.wikipedia.org/wiki/Machine_epsilon
# use values slightly larger than the informal definition
epsilon64 = 2.23e-16  # epsilon
epsilon32 = 1.20e-07  # np.finfo(zero.astype("float32")).eps

# These are PRMS conventions, should not be used elsewhere
nearzero = 1.0e-6
dnearzero = epsilon64
closezero = epsilon32

fill_value_f4 = 9.96921e36

# work in progress...
fill_values_dict = {
    np.dtype("float64"): np.nan,
    np.dtype("float32"): np.nan,
    np.dtype("float16"): np.nan,
    np.dtype("int64"): None,
    np.dtype("int32"): None,
    np.dtype("int16"): None,
    np.dtype("int8"): None,
    np.dtype("bool"): None,
}

np_type_to_netcdf_type_dict = {
    np.dtype("float64"): "f8",
    np.dtype("float32"): "f4",
    np.dtype("int64"): "i8",
    np.dtype("int32"): "i4",
    np.dtype("int16"): "i2",
    np.dtype("int8"): "i1",
    np.dtype("uint64"): "u8",
    np.dtype("uint32"): "u4",
    np.dtype("uint16"): "u2",
    np.dtype("uint8"): "u1",
    np.dtype("bool"): None,
}

inch2cm = 2.54
ft2_per_acre = 43560.0
inches_per_foot = 12.0

ndoy = 366
nmonth = 12


class HruType(Enum):
    INACTIVE = 0
    LAND = 1
    LAKE = 2
    SWALE = 3


class CovType(Enum):
    BARESOIL = 0
    GRASSES = 1
    SHRUBS = 2
    TREES = 3
    CONIFEROUS = 4


class SoilType(Enum):
    SAND = 1
    LOAM = 2
    CLAY = 3


class ETType(Enum):
    ET_DEFAULT = 1
    EVAP_ONLY = 2
    EVAP_PLUS_TRANSP = 3


class SegmentType(Enum):
    SEGMENT = 0
    HEADWATER = 1
    LAKE = 2
    REPLACEINFLOW = 3
    INBOUNDNHM = 4
    OUTBOUNDNHM = 5
    INBOUNDREGION = 6
    OUTBOUNDREGION = 7
    OUTBOUNDOCEAN = 8
    SINK = 9
    INBOUNDGREATLAKES = 10
    OUTBOUNDGREATLAKES = 11
