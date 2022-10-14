import pathlib as pl
from enum import Enum
from typing import Union

import numpy as np

fileish = Union[str, pl.Path]

# PRMS6 Constants module:
# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/misc/m_constants.f90

__pynhm_root__ = pl.Path(__file__).parent

zero = np.zeros([1])[0]
one = np.ones([1])[0]
nan = np.nan

epsilon = np.finfo(zero).eps
epsilon64 = epsilon
epsilon32 = np.finfo(zero.astype("float32")).eps

fill_value_f4 = 9.96921e36

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
