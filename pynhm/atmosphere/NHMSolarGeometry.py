import math
import warnings
from typing import Tuple

# from numba import jit
import numpy as np

from ..utils.parameters import PrmsParameters
from ..base.StateAccess import StateAccess

# The solar geometry model for NHM/PRMS
# Primary reference
# * Appendix E of Dingman, S. L., 1994,
#   Physical Hydrology. Englewood Cliffs, NJ: Prentice Hall, 575 p.

# Lots of constants for this model
# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/c_solar_radiation.f90

# These could be in a constants module
zero = np.zeros(1)[0]
one = np.ones(1)[0]
pi = math.pi


# This could be a util
def epsilon(array: np.ndarray):
    return np.finfo(array.dtype).eps


n_days_per_year = 366
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


# Dimensions
# The time dimension is n_days_per_year (known apriori)
# The spatial dimension n_hru (only known on init from parameters)

# may not use this if they cant be called with jit
def tile_space_to_time(arr: np.ndarray) -> np.ndarray:
    return np.tile(arr, (n_days_per_year, 1))


def tile_time_to_space(arr: np.ndarray, n_hru) -> np.ndarray:
    return np.transpose(np.tile(r1, (n_hru, 1)))


# JLM metadata ?


class NHMSolarGeometry(StateAccess):
    def __init__(
        self,
        parameters: PrmsParameters,
    ):
        # JLM: Document. these names are bad, fix them when it's tested.
        # JLM: It would be nice to inherit state accessors here
        super().__init__()
        self._potential_variables = [
            "hru_cossl",
            "potential_sw_rad_flat",
            "potential_sw_rad",
            "sun_hrs",
        ]
        self.parameters = parameters
        asdf
        # JLM: This should be in the base class for handling space
        if "nhm_id" in parameters._parameter_data.keys():
            space
        self._coords = ["julian_day"]
        self["julian_day"] = julian_days
        self._compute_solar_geometry()
        # dimensions

        return None

    def _compute_solar_geometry(self, parameters: PrmsParameters):
        params = parameters.parameters
        n_hru = parameters.parameters.nhru

        self["hru_cossl"] = np.cos(np.arctan(params["hru_slope"]))

        # The potential radiation on horizontal surfce
        self["potential_sw_rad_flat"], _ = self.compute_soltab(
            np.zeros(n_hru),
            np.zeros(n_hru),
            params["hru_lat"],
            self.compute_t,
            self.func3,
        )

        # The potential radiaton given slope and aspect
        self["potential_sw_rad"], self["sun_hrs"] = self.compute_soltab(
            params["hru_slope"],
            params["hru_aspect"],
            params["hru_lat"],
            self.compute_t,
            self.func3,
        )

        return

    # @jit
    @staticmethod
    def compute_soltab(
        slopes: np.ndarray,
        aspects: np.ndarray,
        lats: np.ndarray,
        compute_t: callable,
        func3: callable,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Swift's daily potential solar radiation and number of hours
        of duration on a sloping surface in [cal/cm2/day].
        Swift, 1976, equation 6.

        Arguments:
          cossl: cos(atan(hru_slope)) [n_hru] ?
          slope: slope [n_hru]
          aspect: aspect [n_hru]
          latitude: latitude [n_hru]

        Return Values: (solt, sunh)
          solt: Swift's potential solar radiation on a sloping surface in
            [cal/cm2/day]. Swift, 1976, equation 6.
            Dimensions: [n_days_per_year, n_hru]
          sunh: The number of hours of direct sunlight
            Dimensions: [n_days_per_year, n_hru]
        """
        n_hru = len(slopes)

        # Slope derived quantities
        sl = np.arctan(slopes)
        sl_sin = np.sin(sl)
        sl_cos = np.cos(sl)

        # Aspect derived quantities
        aspects_rad = np.radians(aspects)
        aspects_cos = np.cos(aspects_rad)

        # Latitude derived quantities
        x0 = np.radians(lats)
        x0_cos = np.cos(x0)

        # x1 latitude of equivalent slope
        # This is equation 13 from Lee, 1963
        x1 = np.arcsin(sl_cos * np.sin(x0) + sl_sin * x0_cos * aspects_cos)

        # d1 is the denominator of equation 12, Lee, 1963
        d1 = sl_cos * x0_cos - sl_sin * np.sin(x0) * aspects_cos
        eps_d1 = epsilon(d1)
        wh_d1_lt_eps = np.where(d1 < eps_d1)
        if len(wh_d1_lt_eps[0]) > 0:
            d1[wh_d1_lt_eps] = eps_d1

        # x2 is the difference in longitude between the location of
        # the HRU and the equivalent horizontal surface expressed in angle hour
        # This is equation 12 from Lee, 1963
        x2 = np.arctan(sl_sin * np.sin(aspects_rad) / d1)
        wh_d1_lt_zero = np.where(d1 < zero)
        if len(wh_d1_lt_zero[0]) > 0:
            x2[wh_d1_lt_zero] = x2[wh_d1_lt_zero] + pi

        # The hour angle from the local meridian (local solar noon) to the
        # sunrise (negative) or sunset (positive)
        # t6: is the hour angle of sunrise on the equivalent slope
        # t7: is the hour angle of sunset on the equivalent slope
        tt = compute_t(x1, solar_declination)
        t6 = (-1 * tt) - x2
        t7 = tt - x2

        # Hours of sunrise and sunset on a horizontal surface at lat
        # t0: is the hour angle of sunrise on a hroizontal surface at the HRU
        # t1: is the hour angle of sunset on a hroizontal surface at the HRU
        tt = compute_t(x0, solar_declination)
        t0 = -1 * tt
        t1 = tt

        # For HRUs that have an east or west direction component to their
        # aspect, the longitude adjustment (moving the effective slope east
        # or west) will cause either:
        # (1) sunrise to be earlier than at the horizontal plane at the HRU
        # (2) sunset to be later than at the horizontal plane at the HRU
        # This is not possible. The if statements below check for this and
        # adjust the sunrise/sunset angle hours on the equivalent slopes as
        # necessary.

        # t2: is the hour angle of sunset on the slope at the HRU
        # t3: is the hour angle of sunrise on the slope at the HRU
        t3 = t7
        wh_t7_gt_t1 = np.where(t7 > t1)
        if len(wh_t7_gt_t1[0]) > 0:
            t3[wh_t7_gt_t1] = t1[wh_t7_gt_t1]

        t2 = t6
        wh_t6_lt_t0 = np.where(t6 < t0)
        if len(wh_t6_lt_t0[0]) > 0:
            t2[wh_t6_lt_t0] = t0[wh_t6_lt_t0]

        # JLM: vectorizing requires are reverse looped order

        t6 = t6 + two_pi
        t7 = t7 - two_pi
        wh_t3_lt_t2 = np.where(t3 < t2)
        if len(wh_t3_lt_t2[0]):
            t2[wh_t3_lt_t2] = zero
            t3[wh_t3_lt_t2] = zero

        # This is if no other conditions are met
        solt = func3(x2, x1, t3, t2)
        sunh = (t3 - t2) * pi_12

        # t7 > t0
        wh_t7_gt_t0 = np.where(t7 > t0)
        if len(wh_t7_gt_t0[0]):
            solt[wh_t7_gt_t0] = (
                func3(x2, x1, t3, t2)[wh_t7_gt_t0]
                + func3(x2, x1, t7, t0)[wh_t7_gt_t0]
            )
            sunh[wh_t7_gt_t0] = (t3 - t2 + t7 - t0)[wh_t7_gt_t0] * pi_12

        # t6 < t1
        wh_t6_lt_t1 = np.where(t6 < t1)
        if len(wh_t6_lt_t1[0]):
            solt[wh_t6_lt_t1] = (
                func3(x2, x1, t3, t2)[wh_t6_lt_t1]
                + func3(x2, x1, t1, t6)[wh_t6_lt_t1]
            )
            sunh[wh_t6_lt_t1] = (t3 - t2 + t1 - t6)[wh_t6_lt_t1] * pi_12

        # The first condition checked
        wh_sl_zero = np.where(tile_space_to_time(np.abs(sl)) < epsilon(sl))
        if len(wh_sl_zero[0]):
            solt[wh_sl_zero] = func3(np.zeros(n_hru), x0, t1, t0)[wh_sl_zero]
            sunh[wh_sl_zero] = (t1 - t0)[wh_sl_zero] * pi_12

        wh_sunh_lt_zero = np.where(sunh < epsilon(sunh))
        if len(wh_sunh_lt_zero[0]):
            sunh[wh_sunh_lt_zero] = zero

        wh_solt_lt_zero = np.where(solt < epsilon(solt))
        if len(wh_solt_lt_zero[0]):
            solt[wh_solt_lt_zero] = zero
            warnings.warn(
                f"{len(wh_solt_lt_zero[0])}/{np.product(solt.shape)} "
                f"loacations-times with negative "
                f"potential solar radiation."
            )

        return solt, sunh

    # @jit
    @staticmethod
    def compute_t(
        lats: np.ndarray, solar_declination: np.ndarray
    ) -> np.ndarray:
        # This function is
        #   * named as in prms for historical purposes
        #   * is for numba compilation
        # JLM: why is division by earth's angular velocity not done here?
        """
        The "sunrise" equation
        lats - latitudes of the hrus
        solar_declination: the declination of the sun on all julian days
        result: the angle hour from the local meridian (local solar noon) to
            the sunrise (negative) or sunset (positive).  The Earth rotates at
            the angular speed of 15 degrees/hour (2 pi / 24 hour in radians)
            and, therefore, result*(24/(2pi) radians gives the time of sunrise
            as the number of hours before the local noon, or the time of sunset
            as the number of hours after the local noon. Here the term local
            noon indicates the local time when the sun is exactly to the
            south or north or exactly overhead.

        https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation.f90
        """
        n_hru = len(lats)
        lats_mat = np.tile(-1 * np.tan(lats), (n_days_per_year, 1))
        sol_dec_mat = np.transpose(
            np.tile(np.tan(solar_declination), (n_hru, 1))
        )
        tx = lats_mat * sol_dec_mat
        result = np.copy(tx)
        result[np.where((tx >= (-1 * one)) & (tx <= one))] = np.arccos(
            tx[np.where((tx >= (-1 * one)) & (tx <= one))]
        )
        result[np.where(tx < (-1 * one))] = pi
        result[np.where(tx > one)] = zero
        return result

    # @jit
    @staticmethod
    def func3(
        v: np.ndarray,
        w: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
    ) -> np.ndarray:
        """
        This is the radian angle version of FUNC3 (eqn 6) from Swift, 1976
        or Lee, 1963 equation 5.
        result: (R4) is potential solar radiation on the surface cal/cm2/day [n_days_per_year, n_hru]
        v: (L2) latitude angle hour offset between actual and equivalent slope [n_hru]
        w: (L1) latitude of the equivalent slope [n_hru]
        x: (T3) hour angle of sunset on equivalent slope [n_days_per_year, n_hru]
        y: (T2) hour angle of sunrise on equivalent slope [n_days_per_year, n_hru]

        # Constants
        r1: solar constant for 60 minutes [n_days_per_year]
        solar_declination: declination of sun [n_days_per_year]

        https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation.f90
        """
        # Must alter the dimensions of the inputs to get the outputs.
        n_hru = len(v)
        vv = np.tile(v, (n_days_per_year, 1))
        ww = np.tile(w, (n_days_per_year, 1))
        # These are known at init time, not sure they are worth saving in self
        # and passing
        rr = np.transpose(np.tile(r1, (n_hru, 1)))
        dd = np.transpose(np.tile(solar_declination, (n_hru, 1)))

        f3 = (
            rr
            * pi_12
            * (
                np.sin(dd) * np.sin(ww) * (x - y)
                + np.cos(dd) * np.cos(ww) * (np.sin(x + vv) - np.sin(y + vv))
            )
        )
        return f3
