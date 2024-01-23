#
# use function from tonic as reference
# tonic: two independent basins
# pws: run two hrus simultaneously.

# SHOULD USE DRB FORCINGS...

import pathlib as pl

import numpy as np
import xarray as xr
from pandas import Timestamp

import pywatershed as pws

dom_dir = pl.Path("../test_data/snow17")

forcing_files = sorted(
    [
        dom_dir / "snow17_forcing_data_HHWM8IL.nc",
        dom_dir / "snow17_forcing_data_HHWM8IU.nc",
    ]
)
forcing_dict = {
    ff.name.split("_")[-1].split(".")[0]: ff for ff in forcing_files
}

in2mm = 2.54 * 10.0
mm2in = 1.0 / in2mm


def test_compare_tonic():
    forcings = {kk: xr.open_dataset(vv) for kk, vv in forcing_dict.items()}

    tonic_results = {}
    for kk, vv in forcings.items():
        tonic_results[kk] = tonic_snow17(
            time=vv.time.values,
            prec=vv.net_ppt.values * in2mm,
            tair=vv.tavgc.values,
        )

    tonic_results = list(tonic_results.values())

    input_dir = pl.Path(dom_dir / "two_forcing_data")
    control = pws.Control(
        start_time=vv.time.values[0],
        end_time=vv.time.values[-1],
        time_step=np.timedelta64(24, "h"),
        options={"input_dir": input_dir},
    )
    parameters = pws.Parameters.from_netcdf(
        dom_dir / "two_snow17_params_default.nc"
    )
    discretization = pws.Parameters.from_netcdf(
        dom_dir / "two_snow17_params_default.nc"
    )
    input_variables = {}
    for key in pws.Snow17.get_inputs():
        nc_path = input_dir / f"{key}.nc"
        input_variables[key] = nc_path

    pws_snow17 = pws.Snow17(
        control, discretization, parameters, **input_variables
    )

    tol = 1e-13
    for istep in range(control.n_times):
        control.advance()
        pws_snow17.advance()
        pws_snow17.calculate(1.0)
        pws_snow17.output()

        for hh in range(len(forcing_files)):
            assert (
                abs(
                    pws_snow17.pkwater_equiv[hh]
                    - tonic_results[hh][0][istep] * mm2in
                )
                < tol
            )
            assert abs(
                abs(
                    pws_snow17.snowmelt[hh]
                    - tonic_results[hh][1][istep] * mm2in
                )
                < tol
            )

    pws_snow17.finalize()

    return


def tonic_snow17(
    time,
    prec,
    tair,
    lat=50,
    elevation=0,
    dt=24,
    scf=1.0,
    rvs=1,
    uadj=0.04,
    mbase=1.0,
    mfmax=1.05,
    mfmin=0.6,
    tipm=0.1,
    nmf=0.15,
    plwhc=0.04,
    pxtemp=1.0,
    pxtemp1=-1.0,
    pxtemp2=3.0,
):
    """
    Snow-17 accumulation and ablation model. This version of Snow-17 is
    intended for use at a point location.
    The time steps for precipitation and temperature must be equal for this
    code.

    Parameters
    ----------
    time : 1d numpy.ndarray or scalar
        Array of datetime objects.
    prec : 1d numpy.ndarray or scalar
        Array of precipitation forcings, size of `time`, units mm.
    tair : 1d numpy.ndarray or scalar
        Array of air temperature forcings, size of `time`, units deg C.
    lat : float, optional
        Latitude of simulation point or grid cell.
    elevation : float, optional
        Elevation of simulation point or grid cell. Default is 0.
    dt : float, optional
        Timestep in hours, default is 24 hours but should always match the
        timestep in `time`.
    scf : float, optional
        Gauge under-catch snow correction factor. Default is 1.0.
    rvs : {0, 1, 2}, optional
        Rain vs. Snow option. Default value of 1 is a linear transition between
        2 temperatures (pxtemp1 and pxtemp2).
    uadj : float, optional
        Average wind function during rain on snow (mm/mb). Default value of
        0.04 is from the American River Basin from Shamir & Georgakakos 2007.
    mbase : float, optional
        Base temperature above which melt typically occurs (deg C). Default
        value of 1 is from the American River Basin from Shamir &
        Georgakakos 2007.
        Note:  must be greater than 0 deg C.
    mfmax : float, optional
        Maximum melt factor during non-rain periods (mm/deg C 6 hr) - in
        western facing slope assumed to occur on June 21. Default value of 1.05
        is from the American River Basin fromShamir & Georgakakos 2007.
    mfmin : float, optional
        Minimum melt factor during non-rain periods (mm/deg C 6 hr) - in
        western facing slope assumed to occur on December 21. Default value of
        0.60 is from the American River Basin from Shamir & Georgakakos 2007.
    tipm : float, optional
        Model parameter (>0.0 and <1.0) - Anderson Manual recommends 0.1 to 0.2
        for deep snowpack areas.  Default is 0.1.
    nmf : float, optional
        Percent liquid water holding capacity of the snow pack - max is 0.4.
        Default value of 0.04 is from the American River Basin from Shamir &
        Georgakakos 2007.
    plwhc : float, optional
        percent liquid water holding capacity of the snow pack - max is 0.4.
        Default value of 0.04 for the American River Basin from Shamir &
        Georgakakos 2007.
    pxtemp : float, optional
        Temperature dividing rain from snow, deg C - if temp is less than or
        equal to pxtemp, all precip is snow.  Otherwise it is rain.  Default is
        1.0.
    pxtemp1 : float, optional
        Lower Limit Temperature dividing tranistion from snow, deg C - if temp
        is less than or equal to pxtemp1, all precip is snow.  Otherwise it is
        mixed linearly. Default is -1.0.
    pxtemp2 : float, optional
        Upper Limit Temperature dividing rain from transition, deg C - if temp
        is greater than or equal to pxtemp2, all precip is rain.  Otherwise it
        is mixed linearly. Default is 3.0.

    Returns
    ----------
    model_swe : numpy.ndarray
        Simulated snow water equivalent.
    outflow : numpy.ndarray
        Simulated runoff outflow.
    """

    # Convert to numpy array if scalars
    time = np.asarray(time)
    prec = np.asarray(prec)
    tair = np.asarray(tair)
    assert time.shape == prec.shape == tair.shape

    # Initialization
    # Antecedent Temperature Index, deg C
    ait = 0.0
    # Liquid water capacity
    w_qx = 0.0
    # Liquid water held by the snow (mm)
    w_q = 0.0
    # accumulated water equivalent of the iceportion of the snow cover (mm)
    w_i = 0.0
    # Heat deficit, also known as NEGHS, Negative Heat Storage
    deficit = 0.0

    # number of time steps
    nsteps = len(time)
    model_swe = np.zeros(nsteps)
    outflow = np.zeros(nsteps)

    # Stefan-Boltzman constant (mm/K/hr)
    stefan = 6.12 * (10 ** (-10))
    # atmospheric pressure (mb) where elevation is in HUNDREDS of meters
    # (this is incorrectly stated in the manual)
    p_atm = 33.86 * (
        29.9
        - (0.335 * elevation / 100)
        + (0.00022 * ((elevation / 100) ** 2.4))
    )

    transitionx = [pxtemp1, pxtemp2]
    transitiony = [1.0, 0.0]

    tipm_dt = 1.0 - ((1.0 - tipm) ** (dt / 6))

    # Model Execution
    for i, dt64 in enumerate(time):
        t = Timestamp(dt64).to_pydatetime()
        mf = melt_function(t, dt, lat, mfmax, mfmin)

        # air temperature at this time step (deg C)
        t_air_mean = tair[i]
        # precipitation at this time step (mm)
        precip = prec[i]

        # Divide rain and snow
        if rvs == 0:
            if t_air_mean <= pxtemp:
                # then the air temperature is cold enough for snow to occur
                fracsnow = 1.0
            else:
                # then the air temperature is warm enough for rain
                fracsnow = 0.0
        elif rvs == 1:
            if t_air_mean <= pxtemp1:
                fracsnow = 1.0
            elif t_air_mean >= pxtemp2:
                fracsnow = 0.0
            else:
                fracsnow = np.interp(t_air_mean, transitionx, transitiony)
        elif rvs == 2:
            fracsnow = 1.0
        else:
            raise ValueError("Invalid rain vs snow option")

        fracrain = 1.0 - fracsnow

        # Snow Accumulation
        # water equivalent of new snowfall (mm)
        pn = precip * fracsnow * scf
        # w_i = accumulated water equivalent of the ice portion of the snow
        # cover (mm)
        w_i += pn
        e = 0.0
        # amount of precip (mm) that is rain during this time step
        rain = fracrain * precip

        # Temperature and Heat deficit from new Snow
        if t_air_mean < 0.0:
            t_snow_new = t_air_mean
            # delta_hd_snow = change in the heat deficit due to snowfall (mm)
            delta_hd_snow = -(t_snow_new * pn) / (80 / 0.5)
            t_rain = pxtemp
        else:
            t_snow_new = 0.0
            delta_hd_snow = 0.0
            t_rain = t_air_mean

        # Antecedent temperature Index
        if pn > (1.5 * dt):
            ait = t_snow_new
        else:
            # Antecedent temperature index
            ait = ait + tipm_dt * (t_air_mean - ait)
        if ait > 0:
            ait = 0

        # Heat Exchange when no Surface melt
        # delta_hd_t = change in heat deficit due to a temperature gradient(mm)
        delta_hd_t = nmf * (dt / 6.0) * ((mf) / mfmax) * (ait - t_snow_new)

        # Rain-on-snow melt
        # saturated vapor pressure at t_air_mean (mb)
        e_sat = (
            2.7489 * (10**8) * np.exp((-4278.63 / (t_air_mean + 242.792)))
        )
        # 1.5 mm/ 6 hrs
        if rain > (0.25 * dt):
            # melt (mm) during rain-on-snow periods is:
            m_ros1 = np.maximum(
                stefan * dt * (((t_air_mean + 273) ** 4) - (273**4)), 0.0
            )
            m_ros2 = np.maximum((0.0125 * rain * t_rain), 0.0)
            m_ros3 = np.maximum(
                (
                    8.5
                    * uadj
                    * (dt / 6.0)
                    * (((0.9 * e_sat) - 6.11) + (0.00057 * p_atm * t_air_mean))
                ),
                0.0,
            )
            m_ros = m_ros1 + m_ros2 + m_ros3
        else:
            m_ros = 0.0

        # Non-Rain melt
        if rain <= (0.25 * dt) and (t_air_mean > mbase):
            # melt during non-rain periods is:
            m_nr = (mf * (t_air_mean - mbase)) + (0.0125 * rain * t_rain)
        else:
            m_nr = 0.0

        # Ripeness of the snow cover
        melt = m_ros + m_nr
        if melt <= 0:
            melt = 0.0

        if melt < w_i:
            w_i = w_i - melt
        else:
            melt = w_i + w_q
            w_i = 0.0

        # qw = liquid water available melted/rained at the snow surface (mm)
        qw = melt + rain
        # w_qx = liquid water capacity (mm)
        w_qx = plwhc * w_i
        # deficit = heat deficit (mm)
        deficit += delta_hd_snow + delta_hd_t

        # limits of heat deficit
        if deficit < 0:
            deficit = 0.0
        elif deficit > 0.33 * w_i:
            deficit = 0.33 * w_i

        # Snow cover is ripe when both (deficit=0) & (w_q = w_qx)
        if w_i > 0.0:
            if (qw + w_q) > ((deficit * (1 + plwhc)) + w_qx):
                # THEN the snow is RIPE
                # Excess liquid water (mm)
                e = qw + w_q - w_qx - (deficit * (1 + plwhc))
                # fills liquid water capacity
                w_q = w_qx
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i = w_i + deficit
                deficit = 0.0
            elif (qw >= deficit) and (
                (qw + w_q) <= ((deficit * (1 + plwhc)) + w_qx)
            ):
                # THEN the snow is NOT yet ripe, but ice is being melted
                e = 0.0
                w_q = w_q + qw - deficit
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i = w_i + deficit
                deficit = 0.0
            else:
                # (qw < deficit) %elseif ((qw + w_q) < deficit):
                # THEN the snow is NOT yet ripe
                e = 0.0
                # w_i increases because water refreezes as heat deficit is
                # decreased
                w_i = w_i + qw
                deficit = deficit - qw
            swe = w_i + w_q
        else:
            e = qw
            swe = 0

        if deficit == 0:
            ait = 0

        model_swe[i] = swe  # total swe (mm) at this time step
        outflow[i] = e

    return model_swe, outflow


def melt_function(t, dt, lat, mfmax, mfmin):
    """
    Seasonal variation calcs - indexed for Non-Rain melt

    Parameters
    ----------
    t : datetime object
        Datetime ojbect for current timestep.
    dt : float
        Timestep in hours.
    lat : float
        Latitude of simulation point or grid cell.
    mfmax : float
        Maximum melt factor during non-rain periods (mm/deg C 6 hr) - in
        western facing slope assumed to occur on June 21. Default value of 1.05
        is from the American River Basin fromShamir & Georgakakos 2007.
    mfmin : float
        Minimum melt factor during non-rain periods (mm/deg C 6 hr) - in
        western facing slope assumed to occur on December 21. Default value of
        0.60 is from the American River Basin from Shamir & Georgakakos 2007.

    Returns
    ----------
    meltf : float
        Melt function for current timestep.
    """
    tt = t.timetuple()
    jday = tt[-2]
    n_mar21 = jday - 80
    days = 365

    # seasonal variation
    sv = (0.5 * np.sin((n_mar21 * 2 * np.pi) / days)) + 0.5
    if lat < 54:
        # latitude parameter, av=1.0 when lat < 54 deg N
        av = 1.0
    else:
        if jday <= 77 or jday >= 267:
            # av = 0.0 from September 24 to March 18,
            av = 0.0
        elif jday >= 117 and jday <= 227:
            # av = 1.0 from April 27 to August 15
            av = 1.0
        elif jday >= 78 and jday <= 116:
            # av varies linearly between 0.0 and 1.0 from 3/19-4/26 and
            # between 1.0 and 0.0 from 8/16-9/23.
            av = np.interp(jday, [78, 116], [0, 1])
        elif jday >= 228 and jday <= 266:
            av = np.interp(jday, [228, 266], [1, 0])

    meltf = (dt / 6) * ((sv * av * (mfmax - mfmin)) + mfmin)
    return meltf
