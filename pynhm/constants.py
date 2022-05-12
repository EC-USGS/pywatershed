import pathlib as pl
from enum import Enum

import numpy as np

__pynhm_root__ = pl.Path(__file__).parent

zero = np.zeros([1])[0]
one = np.ones([1])[0]
nan = np.nan

fill_value_f4 = 9.96921e36


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
