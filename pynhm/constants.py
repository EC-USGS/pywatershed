import pathlib as pl
from enum import Enum

import numpy as np

# PRMS6 Constatns module:
# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/misc/m_constants.f90

__pynhm_root__ = pl.Path(__file__).parent

zero = np.zeros([1])[0]
one = np.ones([1])[0]
nan = np.nan

epsilon = np.finfo(zero).eps

fill_value_f4 = 9.96921e36


class HruType(Enum):
    INACTIVE = 0
    LAND = 1
    LAKE = 2
    SWALE = 3


inch2cm = 2.54
