import pathlib as pl
import warnings
from typing import Tuple, Union

# from numba import jit
import numpy as np

# would like to not subclass storageUnit but it is much simpler to do so
from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.netcdf_utils import NetCdfWrite

from ..base.control import Control
from ..constants import epsilon32, nan, one, zero
from ..utils.prms5util import load_soltab_debug
from .solar_constants import (
    n_days_per_year,
    pi,
    pi_12,
    r1,
    solar_declination,
    two_pi,
)

fileish = Union[str, pl.Path]
epsilon = epsilon32

doy = np.arange(366) + 1

# The solar geometry model for NHM/PRMS
# Primary reference
# * Appendix E of Dingman, S. L., 1994,
#   Physical Hydrology. Englewood Cliffs, NJ: Prentice Hall, 575 p.

# Dimensions
# The time dimension is n_days_per_year (known apriori)
# The spatial dimension n_hru (only known on init from parameters)


def tile_space_to_time(arr: np.ndarray) -> np.ndarray:
    return np.tile(arr, (n_days_per_year, 1))


# def tile_time_to_space(arr: np.ndarray, n_hru) -> np.ndarray:
#    return np.transpose(np.tile(arr, (n_hru, 1)))


class PRMSSolarGeometry(StorageUnit):
    def __init__(
        self,
        control: Control,
        verbose: bool = False,
        netcdf_output_dir: fileish = None,
        from_prms_file: fileish = None,
        from_nc_files_dir: fileish = None,
    ):
        """PRMS Solar Geometry."""

        self.set_inputs(locals())
        budget_type = None
        self.set_budget(budget_type)
        self.netcdf_output_dir = netcdf_output_dir
        self._time = doy

        super().__init__(control=control, verbose=verbose)
        self.name = "PRMSSolarGeometry"

        if from_prms_file:
            (
                self.soltab_potsw.data[:],
                self.soltab_horad_potsw.data[:],
                self.soltab_sunhrs.data[:],
            ) = load_soltab_debug(from_prms_file)

        elif from_nc_files_dir:
            raise NotImplementedError()

        self._calculated = False

        if self.netcdf_output_dir:
            self._output_netcdf = True
            self._calculate_all_time()
            self.netcdf_output_dir = pl.Path(netcdf_output_dir)
            assert self.netcdf_output_dir.exists()
            self._write_netcdf_timeseries()

        else:
            self._output_netcdf = False

        return

    @staticmethod
    def get_inputs() -> tuple:
        return ()

    @staticmethod
    def get_variables():
        """Get solar geometry variables

        Returns:
            variables: output variables

        """
        return (
            "soltab_horad_potsw",
            "soltab_potsw",
            "soltab_sunhrs",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "soltab_horad_potsw": nan,
            "soltab_potsw": nan,
            "soltab_sunhrs": nan,
        }

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "doy",
            "nhru",
            "radadj_intcp",
            "radadj_slope",
            "tmax_index",
            "dday_slope",
            "dday_intcp",
            "radmax",
            "ppt_rad_adj",
            "tmax_allsnow",
            "tmax_allrain_offset",
            "hru_slope",
            "radj_sppt",
            "radj_wppt",
            "hru_lat",
            "hru_area",
            "hru_aspect",
        )

    def set_initial_conditions(self):
        return

    def _calculate_all_time(self):
        self._hru_cossl = np.cos(np.arctan(self["hru_slope"]))
        # The potential radiation on horizontal surfce
        self.soltab_horad_potsw.data[:], _ = self.compute_soltab(
            np.zeros(self["nhru"]),
            np.zeros(self["nhru"]),
            self["hru_lat"],
            self.compute_t,
            self.func3,
        )

        # The potential radiaton given slope and aspect
        (
            self.soltab_potsw.data[:],
            self.soltab_sunhrs.data[:],
        ) = self.compute_soltab(
            self["hru_slope"],
            self["hru_aspect"],
            self["hru_lat"],
            self.compute_t,
            self.func3,
        )

        self._calculated = True
        return

    def _advance_variables(self):
        """Advance solar geometry

        Returns:
            None

        """
        if not self._calculated:
            self._calculate_all_time()
            self._write_netcdf_timeseries()
        return

    def _calculate(self, time_length):
        """Calculate solar_geom"

        Sets the current value to the precalculated values.

        Args:
            time_length: time step length

        Returns:
            None

        """
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
          cossl: cos(atan(hru_slope)) [nhru] ?
          slope: slope [nhru]
          aspect: aspect [nhru]
          latitude: latitude [nhru]

        Return Values: (solt, sunh)
          solt: Swift's potential solar radiation on a sloping surface in
            [cal/cm2/day]. Swift, 1976, equation 6.
            Dimensions: [n_days_per_year, nhru]
          sunh: The number of hours of direct sunlight
            Dimensions: [n_days_per_year, nhru]
        """
        nhru = len(slopes)

        # Slope derived quantities
        sl = np.arctan(slopes)
        sl_sin = np.sin(sl)
        sl_cos = np.cos(sl)

        # Aspect derived quantities
        aspects_rad = np.radians(aspects)
        aspects_cos = np.cos(aspects_rad)

        # Latitude derived quantities
        x0 = np.radians(lats)
        # errors in PRMS ingesting latitudes can be seen here (or earlier)
        x0_cos = np.cos(x0)

        # x1 latitude of equivalent slope
        # This is equation 13 from Lee, 1963
        x1 = np.arcsin(sl_cos * np.sin(x0) + sl_sin * x0_cos * aspects_cos)

        # d1 is the denominator of equation 12, Lee, 1963
        d1 = sl_cos * x0_cos - sl_sin * np.sin(x0) * aspects_cos
        eps_d1 = epsilon
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
        wh_sl_zero = np.where(tile_space_to_time(np.abs(sl)) < epsilon)
        if len(wh_sl_zero[0]):
            solt[wh_sl_zero] = func3(np.zeros(nhru), x0, t1, t0)[wh_sl_zero]
            sunh[wh_sl_zero] = (t1 - t0)[wh_sl_zero] * pi_12

        wh_sunh_lt_zero = np.where(sunh < epsilon)
        if len(wh_sunh_lt_zero[0]):
            sunh[wh_sunh_lt_zero] = zero

        wh_solt_lt_zero = np.where(solt < epsilon)
        if len(wh_solt_lt_zero[0]):
            solt[wh_solt_lt_zero] = zero
            warnings.warn(
                f"{len(wh_solt_lt_zero[0])}/{np.product(solt.shape)} "
                f"locations-times with negative "
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
        nhru = len(lats)
        lats_mat = np.tile(-1 * np.tan(lats), (n_days_per_year, 1))
        sol_dec_mat = np.transpose(
            np.tile(np.tan(solar_declination), (nhru, 1))
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
        result: (R4) is potential solar radiation on the surface cal/cm2/day [n_days_per_year, nhru]
        v: (L2) latitude angle hour offset between actual and equivalent slope [nhru]
        w: (L1) latitude of the equivalent slope [nhru]
        x: (T3) hour angle of sunset on equivalent slope [n_days_per_year, nhru]
        y: (T2) hour angle of sunrise on equivalent slope [n_days_per_year, nhru]

        # Constants
        r1: solar constant for 60 minutes [n_days_per_year]
        solar_declination: declination of sun [n_days_per_year]

        https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation.f90
        """
        # Must alter the dimensions of the inputs to get the outputs.
        nhru = len(v)
        vv = np.tile(v, (n_days_per_year, 1))
        ww = np.tile(w, (n_days_per_year, 1))
        # These are known at init time, not sure they are worth saving in self
        # and passing
        rr = np.transpose(np.tile(r1, (nhru, 1)))
        dd = np.transpose(np.tile(solar_declination, (nhru, 1)))

        f3 = (
            rr
            * pi_12
            * (
                np.sin(dd) * np.sin(ww) * (x - y)
                + np.cos(dd) * np.cos(ww) * (np.sin(x + vv) - np.sin(y + vv))
            )
        )

        return f3

    def _write_netcdf_timeseries(self) -> None:
        if not self._output_netcdf:
            return
        for var in self.variables:
            nc_path = self.netcdf_output_dir / f"{var}.nc"
            nc = NetCdfWrite(
                nc_path,
                self.params.nhm_coordinates,
                [var],
                {var: self.var_meta[var]},
            )
            nc.add_all_data(
                var,
                self[var].data,
                self._time,
                time_coord="doy",
            )
            nc.close()
            assert nc_path.exists()
            print(f"Wrote file: {nc_path}")

        self._output_netcdf = False
        return

    def initialize_netcdf(self, dir):
        self.netcdf_output_dir = dir
        self._output_netcdf = True
        return

    def output(self):
        if self._output_netcdf:
            if self.verbose:
                print(f"writing FULL timeseries output for: {self.name}")
            self._write_netcdf_timeseries()
            self._output_netcdf = False
        return
