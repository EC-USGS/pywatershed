from typing import Literal
from warnings import warn

import numpy as np
from numba import prange

from ..base.adapter import adaptable
from ..base.conservative_process import ConservativeProcess
from ..base.control import Control
from ..constants import numba_num_threads, one, zero
from ..parameters import Parameters

in2mm = 2.54 * 10

# Stefan-Boltzman constant (mm/K/hr)
stefan = 6.12 * (10 ** (-10))


class Snow17(ConservativeProcess):
    """PRMS snow pack.

    The Snow-17 snow representation based on

    Anderson, E. A. "SNOW-17 a model based on temperature and precipitation."
    NOAA, publication SW-102 (1973) and the 2006 update
    https://www.weather.gov/media/owp/oh/hrl/docs/22snow17.pdf

    Implementation based on
    * First: https://github.com/UW-Hydro/tonic/blob/master/tonic/models/
        snow17/snow17.py
    * With ambition to expand to: https://github.com/NOAA-OWP/snow17

    The first being a point representation already in python. The second
    appearing to also include a depletion curve and written in Fotran.

    Args:
        control: a Control object
        discretization: a discretization of class Parameters
        parameters: a parameter object of class Parameters

        tavgc: Average air temperature distributed to each HRU
        net_ppt: Precipitation (rain and/or snow) that falls through the
            canopy for each HRU
        budget_type: one of [None, "warn", "error"]
        calc_method: one of ["fortran", "numba", "numpy"]. None defaults to
            "numba".
        verbose: Print extra information or not?

    """

    def __init__(
        self,
        control: Control,
        discretization: Parameters,
        parameters: Parameters,
        tavgc: adaptable,
        net_ppt: adaptable,
        budget_type: Literal[None, "warn", "error"] = None,
        calc_method: Literal["numba", "numpy"] = "numpy",
        verbose: bool = None,
    ) -> "Snow17":
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "Snow17"

        self._set_inputs(locals())
        self._set_options(locals())

        self._set_budget()
        self._init_calc_method()

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nhru",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "lat",
            "elevation",
            "scf",
            "rvs",
            "uadj",
            "mbase",
            "mfmax",
            "mfmin",
            "tipm",
            "nmf",
            "plwhc",
            "pxtemp",
            "pxtemp1",
            "pxtemp2",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return (
            "tavgc",
            "net_ppt",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "freeh2o": zero,
            "freeh2o_prev": zero,
            "freeh2o_change": zero,
            "pk_def": zero,
            "pk_ice": zero,
            "pk_ice_prev": zero,
            "pk_ice_change": zero,
            "pkwater_equiv": zero,  # model_swe
            "snowmelt": zero,  # outflow
        }

    @staticmethod
    def get_mass_budget_terms():
        return {
            "inputs": [
                "net_ppt",
            ],
            "outputs": [
                "snowmelt",
            ],
            "storage_changes": ["freeh2o_change", "pk_ice_change"],
        }

    @staticmethod
    def get_restart_variables() -> tuple:
        return ("pkwater_equiv",)

    def _set_initial_conditions(self):
        # define self variables without needing metadata.
        # these wont be output
        zero_array = self.pkwater_equiv * zero
        # Liquid water capacity
        self.freeh2o_capacity = zero_array.copy()  # w_qx
        # Antecedent Temperature Index, deg C
        self.ait = zero_array.copy()

        self._dt_hrs = self.control.time_step_seconds / (60 * 60)

    def _init_calc_method(self):
        if self._calc_method is None:
            self._calc_method = "numba"

        if self._calc_method.lower() not in ["numpy", "numba"]:
            msg = (
                f"Invalid calc_method={self._calc_method} for {self.name}. "
                f"Setting calc_method to 'numba' for {self.name}"
            )
            warn(msg, UserWarning)
            self._calc_method = "numba"

        if self._calc_method.lower() == "numba":
            import numba as nb

            numba_msg = f"{self.name} jit compiling with numba "
            nb_parallel = (numba_num_threads is not None) and (
                numba_num_threads > 1
            )
            if nb_parallel:
                numba_msg += f"and using {numba_num_threads} threads"
            print(numba_msg, flush=True)

            self._calculate_snow = nb.njit(
                fastmath=True, parallel=nb_parallel
            )(self._calculate_numpy)

            fns = [
                "_melt_function",
            ]
            for fn in fns:
                setattr(
                    self,
                    fn,
                    nb.njit(fastmath=True)(getattr(self, fn)),
                )

        else:
            self._calculate_snow = self._calculate_numpy

        return

    def _advance_variables(self) -> None:
        self.freeh2o_prev[:] = self.freeh2o
        self.pk_ice_prev[:] = self.pk_ice
        return

    def _calculate(self, simulation_time):
        (
            self.pkwater_equiv[:],
            self.snowmelt[:],
            self.freeh2o[:],
            self.freeh2o_capacity[:],
            self.freeh2o_change[:],
            self.pk_ice[:],
            self.pk_ice_change[:],
            self.pk_def[:],
            self.ait[:],
        ) = self._calculate_snow(
            current_doy=self.control.current_doy,
            dt_hrs=self._dt_hrs,
            nhru=self.nhru,
            tavgc=self.tavgc,
            ait_vec=self.ait,
            net_ppt=self.net_ppt,
            freeh2o=self.freeh2o,
            freeh2o_capacity=self.freeh2o_capacity,
            freeh2o_prev=self.freeh2o_prev,
            pk_ice=self.pk_ice,
            pk_ice_prev=self.pk_ice_prev,
            pk_def=self.pk_def,
            pkwater_equiv=self.pkwater_equiv,
            snowmelt=self.snowmelt,
            elevation=self.elevation,
            lat=self.lat,
            mbase=self.mbase,
            mfmax=self.mfmax,
            mfmin=self.mfmin,
            nmf=self.nmf,
            plwhc=self.plwhc,
            pxtemp=self.pxtemp,
            pxtemp1=self.pxtemp1,
            pxtemp2=self.pxtemp2,
            rvs=self.rvs,
            scf=self.scf,
            tipm=self.tipm,
            uadj=self.uadj,
            melt_function=self._melt_function,
            verbose=self._verbose,
        )

        return

    @staticmethod
    def _calculate_numpy(
        current_doy,
        dt_hrs,
        nhru,
        tavgc,
        ait_vec,
        net_ppt,
        freeh2o,
        freeh2o_capacity,
        freeh2o_prev,
        pk_ice,
        pk_ice_prev,
        pk_def,
        pkwater_equiv,
        snowmelt,
        elevation,
        lat,
        mbase,
        mfmax,
        mfmin,
        nmf,
        plwhc,
        pxtemp,
        pxtemp1,
        pxtemp2,
        rvs,
        scf,
        tipm,
        uadj,
        melt_function,
        verbose,
    ):
        """Calculate snow pack terms for a time step

        Args:
            simulation_time: current simulation time

        Returns:
            None
        """

        for hh in prange(nhru):
            # convert PRMS names to orig names going from array to indiv hru
            # air temperature at this time step (deg C)
            t_air_mean = tavgc[hh]
            ait = ait_vec[hh]
            # convert inches to mm
            precip = net_ppt[hh] * in2mm
            w_qx = freeh2o_capacity[hh] * in2mm
            w_q = freeh2o[hh] * in2mm
            w_i = pk_ice[hh] * in2mm
            # what are the units??
            deficit = pk_def[hh]

            # atmospheric pressure (mb) where elevation is in HUNDREDS of
            # meters (this is incorrectly stated in the manual)
            p_atm = 33.86 * (
                29.9
                - (0.335 * elevation[hh] / 100)
                + (0.00022 * ((elevation[hh] / 100) ** 2.4))
            )

            transitionx = [pxtemp1[hh], pxtemp2[hh]]
            transitiony = [one, zero]
            tipm_dt = one - ((one - tipm[hh]) ** (dt_hrs / 6))

            mf = melt_function(
                current_doy, dt_hrs, lat[hh], mfmax[hh], mfmin[hh]
            )

            # Divide rain and snow
            if rvs[hh] == 0:
                if t_air_mean <= pxtemp[hh]:
                    # then the air temperature is cold enough for snow to occur
                    fracsnow = one
                else:
                    # then the air temperature is warm enough for rain
                    fracsnow = zero
            elif rvs[hh] == 1:
                if t_air_mean <= pxtemp1[hh]:
                    fracsnow = one
                elif t_air_mean >= pxtemp2[hh]:
                    fracsnow = zero
                else:
                    fracsnow = np.interp(t_air_mean, transitionx, transitiony)
            elif rvs[hh] == 2:
                fracsnow = one
            else:
                raise ValueError("Invalid rain vs snow option")

            fracrain = one - fracsnow

            # Snow Accumulation
            # water equivalent of new snowfall (mm)
            pn = precip * fracsnow * scf[hh]
            # w_i: accumulated water equivalent of the ice portion of the snow
            # cover (mm)
            w_i += pn
            # ee: snowmelt in mm
            ee = zero
            # amount of precip (mm) that is rain during this time step
            rain = fracrain * precip

            # Temperature and Heat deficit from new Snow
            if t_air_mean < zero:
                t_snow_new = t_air_mean
                # delta_hd_snow: mm change in the heat deficit due to snowfall
                delta_hd_snow = -(t_snow_new * pn) / (80 / 0.5)
                t_rain = pxtemp[hh]
            else:
                t_snow_new = zero
                delta_hd_snow = zero
                t_rain = t_air_mean

            # Antecedent temperature Index
            if pn > (1.5 * dt_hrs):
                ait = t_snow_new
            else:
                # Antecedent temperature index
                ait = ait + tipm_dt * (t_air_mean - ait)
            if ait > 0:
                ait = 0

            # Heat Exchange when no Surface melt
            # delta_hd_t = change in heat deficit due to a temperature
            # gradient(mm)
            delta_hd_t = (
                nmf[hh]
                * (dt_hrs / 6.0)
                * ((mf) / mfmax[hh])
                * (ait - t_snow_new)
            )

            # Rain-on-snow melt
            # saturated vapor pressure at t_air_mean (mb)
            e_sat = (
                2.7489
                * (10**8)
                * np.exp((-4278.63 / (t_air_mean + 242.792)))
            )
            # 1.5 mm/ 6 hrs
            if rain > (0.25 * dt_hrs):
                # melt (mm) during rain-on-snow periods is:
                m_ros1 = np.maximum(
                    stefan * dt_hrs * (((t_air_mean + 273) ** 4) - (273**4)),
                    zero,
                )
                m_ros2 = np.maximum((0.0125 * rain * t_rain), zero)
                m_ros3 = np.maximum(
                    (
                        8.5
                        * uadj[hh]
                        * (dt_hrs / 6.0)
                        * (
                            ((0.9 * e_sat) - 6.11)
                            + (0.00057 * p_atm * t_air_mean)
                        )
                    ),
                    zero,
                )
                m_ros = m_ros1 + m_ros2 + m_ros3
            else:
                m_ros = zero

            # Non-Rain melt
            if rain <= (0.25 * dt_hrs) and (t_air_mean > mbase[hh]):
                # melt during non-rain periods is:
                m_nr = (mf * (t_air_mean - mbase[hh])) + (
                    0.0125 * rain * t_rain
                )
            else:
                m_nr = zero

            # Ripeness of the snow cover
            melt = m_ros + m_nr
            if melt <= zero:
                melt = zero

            if melt < w_i:
                w_i = w_i - melt
            else:
                melt = w_i + w_q
                w_i = zero

            # qw: liquid water available melted/rained at the snow surface (mm)
            qw = melt + rain
            # w_qx = liquid water capacity (mm)
            w_qx = plwhc[hh] * w_i
            # deficit = heat deficit (mm)
            deficit += delta_hd_snow + delta_hd_t

            # limits of heat deficit
            if deficit < 0:
                deficit = zero
            elif deficit > 0.33 * w_i:
                deficit = 0.33 * w_i

            # Snow cover is ripe when both (deficit=0) & (w_q = w_qx)
            if w_i > zero:
                if (qw + w_q) > ((deficit * (1 + plwhc[hh])) + w_qx):
                    # THEN the snow is RIPE
                    # Excess liquid water (mm)
                    ee = qw + w_q - w_qx - (deficit * (1 + plwhc[hh]))
                    # fills liquid water capacity
                    w_q = w_qx
                    # w_i increases because water refreezes as heat deficit is
                    # decreased
                    w_i = w_i + deficit
                    deficit = zero
                elif (qw >= deficit) and (
                    (qw + w_q) <= ((deficit * (1 + plwhc[hh])) + w_qx)
                ):
                    # THEN the snow is NOT yet ripe, but ice is being melted
                    ee = zero
                    w_q = w_q + qw - deficit
                    # w_i increases because water refreezes as heat deficit is
                    # decreased
                    w_i = w_i + deficit
                    deficit = zero
                else:
                    # (qw < deficit) %elseif ((qw + w_q) < deficit):
                    # THEN the snow is NOT yet ripe
                    ee = zero
                    # w_i increases because water refreezes as heat deficit is
                    # decreased
                    w_i = w_i + qw
                    deficit = deficit - qw
                swe = w_i + w_q
            else:
                ee = qw
                swe = 0

            if deficit == 0:
                ait = 0

            # End of model execution
            pkwater_equiv[hh] = swe / in2mm  # total swe (mm) at this time step
            snowmelt[hh] = ee / in2mm

            # precip only on RHS
            freeh2o_capacity[hh] = w_qx / in2mm
            freeh2o[hh] = w_q / in2mm
            freeh2o_change = freeh2o - freeh2o_prev

            pk_ice[hh] = w_i / in2mm
            pk_ice_change = pk_ice - pk_ice_prev

            pk_def[hh] = deficit  # TODO: convert to langelys?

            ait_vec[hh] = ait  # rename?

        return (
            pkwater_equiv,
            snowmelt,
            freeh2o,
            freeh2o_capacity,
            freeh2o_change,
            pk_ice,
            pk_ice_change,
            pk_def,
            ait_vec,
        )

    @staticmethod
    def _melt_function(doy, dt, lat, mfmax, mfmin):
        """
        Seasonal variation calcs - indexed for Non-Rain melt

        Parameters
        ----------
        dt : float
            Timestep in hours.
        lat : float
            Latitude of simulation point or grid cell.
        mfmax : float
            Maximum melt factor during non-rain periods (mm/deg C 6 hr) - in
            western facing slope assumed to occur on June 21. Default value of
            1.05 is from the American River Basin fromShamir & Georgakakos
            2007.
        mfmin : float
            Minimum melt factor during non-rain periods (mm/deg C 6 hr) - in
            western facing slope assumed to occur on December 21. Default
            value of 0.60 is from the American River Basin from Shamir &
            Georgakakos 2007.

        Returns
        ----------
        meltf : float
            Melt function for current timestep.
        """
        # tt = t.timetuple()
        jday = doy  # tt[-2]
        n_mar21 = jday - 80
        days = 365

        # seasonal variation
        sv = (0.5 * np.sin((n_mar21 * 2 * np.pi) / days)) + 0.5
        if lat < 54:
            # latitude parameter, av=1.0 when lat < 54 deg N
            av = one
        else:
            if jday <= 77 or jday >= 267:
                # av = 0.0 from September 24 to March 18,
                av = zero
            elif jday >= 117 and jday <= 227:
                # av = 1.0 from April 27 to August 15
                av = one
            elif jday >= 78 and jday <= 116:
                # av varies linearly between 0.0 and 1.0 from 3/19-4/26 and
                # between 1.0 and 0.0 from 8/16-9/23.
                av = np.interp(jday, [78, 116], [0, 1])
            elif jday >= 228 and jday <= 266:
                av = np.interp(jday, [228, 266], [1, 0])

        meltf = (dt / 6) * ((sv * av * (mfmax - mfmin)) + mfmin)
        return meltf
