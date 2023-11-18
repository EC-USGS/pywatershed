import math

import numpy as np

from ..constants import ndoy, one

# SolarGeom Constants: Lots of them
# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/c_solar_radiation.f90

pi = math.pi
n_days_per_year = ndoy
n_days_per_year_flt = 365.242
eccentricy = 0.01671
two_pi = 2 * pi
pi_12 = 12 / pi
# JLM: only place prms6 uses 365.242 and commented value is wrong
# rad day is ~0.0172028 not 0.00143356672
rad_day = two_pi / n_days_per_year_flt

julian_days = np.arange(n_days_per_year) + 1

obliquity = 1 - (eccentricy * np.cos((julian_days - 3) * rad_day))

# integer julian day values will yield local noon in the sunrise equation
yy = (julian_days - 1) * rad_day
yy2 = yy * 2
yy3 = yy * 3
solar_declination = (
    0.006918
    - 0.399912 * np.cos(yy)
    + 0.070257 * np.sin(yy)
    - 0.006758 * np.cos(yy2)
    + 0.000907 * np.sin(yy2)
    - 0.002697 * np.cos(yy3)
    + 0.00148 * np.sin(yy3)
)


# Solar constant cal/cm2/min (r0 could also be 1.95 (Drummond, et al 1968))
r0 = 2 * one
# Solar constant for 60 minutes
r1 = (60.0 * r0) / (obliquity**2)

# https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/c_solar_radiation_degday.f90#L19
solf = np.array(
    [
        0.20,
        0.35,
        0.45,
        0.51,
        0.56,
        0.59,
        0.62,
        0.64,
        0.655,
        0.67,
        0.682,
        0.69,
        0.70,
        0.71,
        0.715,
        0.72,
        0.722,
        0.724,
        0.726,
        0.728,
        0.73,
        0.734,
        0.738,
        0.742,
        0.746,
        0.75,
    ]
)
