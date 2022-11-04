import pathlib as pl

import numpy as np

from pynhm.base.storageUnit import StorageUnit
from pynhm.utils.netcdf_utils import NetCdfWrite

from ..base.adapter import adaptable
from ..base.control import Control
from ..constants import epsilon, fileish, inch2cm, nan, one, zero
from ..utils.time_utils import datetime_doy, datetime_month
from .solar_constants import solf


# may not use this if they cant be called with jit
# if it can, put it in a common utility.
def tile_space_to_time(arr: np.ndarray, n_time) -> np.ndarray:
    return np.tile(arr, (n_time, 1))


def tile_time_to_space(arr: np.ndarray, n_space) -> np.ndarray:
    return np.transpose(np.tile(arr, (n_space, 1)))


class PRMSAtmosphere(StorageUnit):
    """PRMS atmospheric boundary layer model.

    Implementation based on PRMS 5.2.1 with theoretical documentation based on
    PRMS-IV:

    Markstrom, S. L., Regan, R. S., Hay, L. E., Viger, R. J., Webb, R. M.,
    Payn, R. A., & LaFontaine, J. H. (2015). PRMS-IV, the
    precipitation-runoff modeling system, version 4. US Geological Survey
    Techniques and Methods, 6, B7.
    https://pubs.usgs.gov/tm/6b7/pdf/tm6-b7.pdf

    This representation uses precipitation and temperature inputs. Relative
    humidity could be added as well.

    The boundary layer calculates and manages the following variables (given
    by PRMSAtmosphere.get_variables()):

        tmaxf,
        tminf,
        prmx,
        hru_ppt,
        hru_rain,
        hru_snow,
        swrad,
        potet,
        transp_on,

    PRMS adjustments to temperature and precipitation are applied here to
    the inputs. Shortwave radiation (using degree day method) and potential
    evapotranspiration (Jensen and Haise ,1963) and a temperature based
    transpiration flag (transp_on) are also calculated.

    Note that all variables are calculate for all time upon initialization. The
    full time version of a variable is given by the "private" version of the
    variable which is named with a single-leading underscore (eg tmaxf for all
    time is _tmaxf).

    This full-time initialization ma not be tractable for large domains and/or
    long periods of time and require changes to batch the processing of the
    variables. The benefits of full-time initialization are 1) the code is
    vectorized and fast for such a large calculation, 2) the initialization of
    this class effectively preprocess all the inputs to the rest of the model
    and can then be skipped in subsequent model calls (unless the parameters
    are changing).

    Args:
        control: control object
        prcp: precipitation cbh netcdf file
        tmax: maximum daily temperature cbh netcdf file
        tmin: minimum daily temperature cbh netcdf file
        budget_type: [None | "warn" |  "error"].
        verbose: bool indicating amount of output to terminal.
        netcdf_output_dir: an existing directory to which to write all
            variables for all time.

    """

    def __init__(
        self,
        control: Control,
        prcp: fileish,
        tmax: fileish,
        tmin: fileish,
        soltab_potsw: adaptable,
        soltab_horad_potsw: adaptable,
        budget_type: str = None,
        verbose: bool = False,
        netcdf_output_dir: fileish = None,
        from_file_dir: fileish = None,
        n_time_chunk: int = -1,
        load_n_time_batches: int = 1,
    ):

        # This could be used to subclass storageUnit or Process classes to have
        # timeseries. Solar geom bas doy dimension not actual simulation times

        # Defering handling batch handling of time chunks but self.n_time_chunk
        # is a dimension used in the metadata/variables dimensions.
        if n_time_chunk <= 0:
            self.n_time_chunk = control.n_times
        else:
            self.n_time_chunk = n_time_chunk

        # Initialize full time with nans
        self._time = np.full(self.n_time_chunk, nan, dtype="datetime64[s]")

        self.netcdf_output_dir = netcdf_output_dir

        super().__init__(
            control=control,
            verbose=verbose,
            load_n_time_batches=load_n_time_batches,
        )
        self._set_inputs(locals())

        self.name = "PRMSAtmosphere"
        self.budget = None

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

    def _calculate_all_time(self):
        # Eventually refactor to work at specific chunks of time
        for input in ["prcp", "tmax", "tmin"]:
            # # this is a bit of a mess: ._dataset.dataset
            input_time = self._input_variables_dict[input]._nc_read.times
            # dont do this...
            # self[input][:] = ds.dataset[input][:].data
            if np.isnan(self._time[0]):
                self._time[:] = input_time
            else:
                assert (input_time == self._time).all()

        start_time_ind = np.where(self._time == self.control._start_time)[0]
        if not len(start_time_ind):
            msg = "Control start_time is not in the input data time"
            raise ValueError(msg)
        self._init_time_ind = start_time_ind[0]

        # Solve all variables for all time
        self._month_ind_12 = datetime_month(self._time) - 1  # (time)
        self._month_ind_1 = np.zeros(self._time.shape, dtype=int)  # (time)
        self._month = datetime_month(self._time)  # (time)
        self._doy = datetime_doy(self._time)  # (time)

        self.adjust_temperature()
        self.adjust_precip()
        self.calculate_sw_rad_degree_day()
        self.calculate_potential_et_jh()
        self.calculate_transp_tindex()

        # Budget is not for all time
        # self.set_budget(budget_type)
        # This is a nonstandard budget.
        self.budget = None
        # self.budget = {
        #     "inputs": {"potet": self.potet},
        #     "outputs": {"hru_actet": self.hru_actet},
        #     "storage_changes": {"available_potet": self.available_potet},
        # }

        # JLM todo: delete large variables on self for memory management
        self._calculated = True
        return

    @staticmethod
    def get_inputs() -> tuple:
        return (
            "prcp",
            "tmax",
            "tmin",
            "soltab_potsw",
            "soltab_horad_potsw",
        )

    @staticmethod
    def get_variables() -> tuple:
        return (
            "tmaxf",
            "tminf",
            "prmx",
            "hru_ppt",
            "hru_rain",
            "hru_snow",
            "swrad",
            "potet",
            # "hru_actet",
            # "available_potet",
            "transp_on",
            "tmaxc",
            "tavgc",
            "tminc",
            "pptmix",
            "orad_hru",
        )

    @staticmethod
    def get_init_values() -> dict:
        return {
            "tmaxf": nan,
            "tminf": nan,
            "prmx": nan,
            "hru_ppt": nan,
            "hru_rain": nan,
            "hru_snow": nan,
            "swrad": nan,
            "potet": nan,
            # "hru_actet": zero,
            # "available_potet": nan,
            "transp_on": 0,
            "tmaxc": nan,
            "tavgc": nan,
            "tminc": nan,
            "pptmix": nan,
            "orad_hru": nan,
        }

    @staticmethod
    def get_parameters():
        return (
            "nmonth",
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
            "nhru",
            "hru_aspect",
            "jh_coef",
            "jh_coef_hru",
            "tmax_cbh_adj",
            "tmin_cbh_adj",
            "tmax_allsnow",
            "tmax_allrain_offset",
            "snow_cbh_adj",
            "rain_cbh_adj",
            "adjmix_rain",
            "transp_beg",
            "transp_end",
            "transp_tmax",
            "radadj_intcp",  # below are solar params used by Atmosphere
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
            "temp_units",
        )

    def _set_initial_conditions(self):
        return

    def _advance_variables(self):
        """Advance the PRMSAtmosphere

        Returns:
            None

        """
        if not self._calculated:
            self._calculate_all_time()

        return

    def _calculate(self, time_length):
        """Calculate PRMSAtmosphere"

        Sets the current value to the precalculated values.

        Args:
            time_length: time step length

        Returns:
            None

        """
        return

    def adjust_temperature(self):
        """Input temperature adjustments using calibrated parameters."""

        # throw an error if these have different shapes
        if self["tmax_cbh_adj"].shape != self["tmin_cbh_adj"].shape:
            msg = (
                "Not implemented: tmin/tmax cbh adj parameters "
                "with different shapes"
            )
            raise NotImplementedError(msg)

        if self["tmax_cbh_adj"].shape[0] == self.nmonth:
            month_ind = self._month_ind_12
        elif self["tmax_cbh_adj"].shape[0] == 1:
            month_ind = self._month_ind_1
        else:
            msg = (
                "Unexpected month dimension for cbh "
                "temperature adjustment params"
            )
            raise ValueError(msg)

        # (time, space) dimensions on these variables
        ivd = self._input_variables_dict
        self.tmaxf.data[:] = ivd["tmax"].data + self.tmax_cbh_adj[month_ind]
        self.tminf.data[:] = ivd["tmin"].data + self.tmin_cbh_adj[month_ind]
        self.tminc.data[:] = (self["tminf"].data - 32.0) * (5 / 9)
        self.tmaxc.data[:] = (self["tmaxf"].data - 32.0) * (5 / 9)
        self.tavgc.data[:] = (self["tmaxc"].data + self["tminc"].data) / 2.0

        return

    def adjust_precip(self):
        """Input precipitation adjustments using calibrated parameters.

        Snow/rain partitioning of total precip depends on adjusted temperature
        in addition to depending on additonal parameters.
        """

        ivd = self._input_variables_dict

        # throw an error shapes are inconsistent
        shape_list = np.array(
            [
                self.tmax_allsnow.shape[0],
                self.tmax_allrain_offset.shape[0],
                self.snow_cbh_adj.shape[0],
                self.rain_cbh_adj.shape[0],
                self.adjmix_rain.shape[0],
            ]
        )
        if not (shape_list == self.nmonth).all():
            msg = (
                "Not implemented: cbh precip adjustment parameters with "
                " different shapes"
            )
            raise NotImplementedError(msg)

        tmax_allrain = self.tmax_allsnow + self.tmax_allrain_offset

        if self.tmax_allsnow.shape[0] == self.nmonth:
            month_ind = self._month_ind_12
        elif self.tmax_allsnow.shape[0] == 1:
            month_ind = self._month_ind_1
        else:
            msg = "Unexpected month dimension for cbh precip adjustment params"
            raise ValueError(msg)

        # Order MATTERS in calculating the prmx mask
        # The logic in PRMS is if(all_snow),elif(all_rain),else(mixed)
        # so we set the mask in the reverse order
        # Calculate the mix everywhere, then set the precip/rain/snow amounts
        # from the conditions.
        tdiff = self.tmaxf.data - self.tminf.data
        tdiff = np.where(tdiff < epsilon, 0.000100, tdiff)
        self.prmx.data[:] = (
            (self.tmaxf.data - self.tmax_allsnow[month_ind]) / tdiff
        ) * self.adjmix_rain[month_ind]

        self.prmx.data[:] = np.where(
            self.prmx.data < zero, zero, self.prmx.data
        )
        self.prmx.data[:] = np.where(self.prmx.data > one, one, self.prmx.data)
        del tdiff

        wh_all_rain = np.where(
            np.logical_or(
                self.tminf.data > self.tmax_allsnow[month_ind],
                self.tmaxf.data >= tmax_allrain[month_ind],
            )
        )
        self.prmx.data[wh_all_rain] = one

        wh_all_snow = np.where(self.tmaxf.data <= self.tmax_allsnow[month_ind])
        self.prmx.data[wh_all_snow] = zero

        # This is in climate_hru as a condition of calling climateflow
        # (eye roll)
        self.prmx.data[:] = np.where(
            ivd["prcp"].data <= zero, zero, self.prmx.data
        )

        self.pptmix.data[:] = np.where(
            (self.prmx.data > zero) & (self.prmx.data < one), 1, 0
        )

        # Recalculate/redefine these now based on prmx instead of the
        # temperature logic
        wh_all_snow = np.where(self.prmx.data <= zero)
        wh_all_rain = np.where(self.prmx.data >= one)

        # Mixed case (everywhere, to be overwritten by the all-snow/rain-fall
        # cases)
        self.hru_ppt.data[:] = ivd["prcp"].data * self.snow_cbh_adj[month_ind]
        self.hru_rain.data[:] = self.prmx.data * self.hru_ppt.data
        self.hru_snow.data[:] = self.hru_ppt.data - self.hru_rain.data

        # All precip is snow case
        # The condition to be used later:
        self.hru_ppt.data[wh_all_snow] = (
            ivd["prcp"].data * self.snow_cbh_adj[month_ind]
        )[wh_all_snow]
        self.hru_snow.data[wh_all_snow] = self.hru_ppt.data[wh_all_snow]
        self.hru_rain.data[wh_all_snow] = zero

        # All precip is rain case
        # The condition to be used later:
        self.hru_ppt.data[wh_all_rain] = (
            ivd["prcp"].data * self.rain_cbh_adj[month_ind]
        )[wh_all_rain]
        self.hru_rain.data[wh_all_rain] = self.hru_ppt.data[wh_all_rain]
        self.hru_snow.data[wh_all_rain] = zero
        return

    def calculate_sw_rad_degree_day(self) -> None:
        """Calculate shortwave radiation using the degree day method."""

        solar_params = {}
        solar_param_names = (
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
        )
        for name in solar_param_names:
            solar_params[name] = self[name]

        ivd = self._input_variables_dict
        self.swrad.data[:], self.orad_hru.data[:] = self._ddsolrad_run(
            dates=self._time,
            tmax_hru=self.tmaxf.data,
            hru_ppt=self.hru_ppt.data,
            soltab_potsw=ivd["soltab_potsw"].data,
            soltab_horad_potsw=ivd["soltab_horad_potsw"].data,
            **solar_params,
            nmonth=self.nmonth,
        )

        return

    # @jit
    @staticmethod
    def _ddsolrad_run(
        dates: np.ndarray,  # [n_time]
        tmax_hru: np.ndarray,  # [n_time, n_hru]
        hru_ppt: np.ndarray,  # [n_time, n_hru]
        soltab_potsw: np.ndarray,  # [n_time, n_hru] ??
        soltab_horad_potsw: np.ndarray,  # [n_time, n_hru] ??
        radadj_intcp,  # param [12, n_hru]
        radadj_slope,  # "    "
        tmax_index,
        dday_slope,
        dday_intcp,
        radmax,
        ppt_rad_adj,
        tmax_allsnow,
        tmax_allrain_offset,
        hru_slope,  # pram [n_hru]
        radj_sppt,
        radj_wppt,
        hru_lat,
        hru_area,
        nmonth,
    ) -> (np.ndarray, np.ndarray):  # [n_time, n_hru]

        # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_solar_radiation_degday.f90
        n_time, n_hru = tmax_hru.shape

        # Transforms of time
        doy = (dates - dates.astype("datetime64[Y]")).astype(
            "timedelta64[h]"
        ).astype(int) / 24 + 1
        doy = tile_time_to_space(doy, n_hru).astype(int)
        month = (dates.astype("datetime64[M]").astype(int) % nmonth) + 1
        month = tile_time_to_space(month, n_hru)

        # is_summer
        # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/misc/sm_prms_time.f90#L114
        is_summer = (doy >= 79) & (doy <= 265)
        northern_hemisphere = hru_lat > zero
        if not any(northern_hemisphere):
            northern_hemisphere = tile_space_to_time(
                northern_hemisphere, n_time
            )
            msg = "Implementation not checked"
            raise NotImplementedError(msg)
            is_summer = is_summer & northern_hemisphere

        hru_cossl = tile_space_to_time(np.cos(np.arctan(hru_slope)), n_time)

        def monthly_to_daily(arr):
            return arr[month[:, 0] - 1, :]

        def doy_to_daily(arr):
            return arr[doy[:, 0] - 1, :]

        # This is the start of the loop over hru & time
        # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_solar_radiation_degday.f90#L103

        # dday
        dday_slope_day = monthly_to_daily(dday_slope)
        dday_intcp_day = monthly_to_daily(dday_intcp)
        dday = (dday_slope_day * tmax_hru) + dday_intcp_day + one
        dday[np.where(dday < one)] = one
        del dday_slope_day, dday_intcp_day

        # radadj
        radmax_day = monthly_to_daily(radmax)
        radadj = monthly_to_daily(radmax)  # the else condition
        wh_dday_lt_26 = np.where(dday < 26.0)
        if len(wh_dday_lt_26[0]):
            kp = dday.astype(int)[wh_dday_lt_26]  # dddayi = float(kp)
            radadj[wh_dday_lt_26] = solf[kp - 1] + (
                (solf[kp] - solf[kp - 1]) * (dday[wh_dday_lt_26] - kp)
            )
            wh_radadj_gt_max = np.where(radadj > radmax_day)
            if len(wh_radadj_gt_max[0]):
                radadj[wh_radadj_gt_max] = radmax_day[wh_radadj_gt_max]
            del kp, wh_radadj_gt_max
        del radmax_day, wh_dday_lt_26

        # pptadj
        pptadj = np.ones((n_time, n_hru))
        radadj_intcp_day = monthly_to_daily(radadj_intcp)
        radadj_slope_day = monthly_to_daily(radadj_slope)
        tmax_index_day = monthly_to_daily(tmax_index)

        ppt_rad_adj_day = monthly_to_daily(ppt_rad_adj)

        # * if
        # This is the outer if. All of the following conditions work on this
        cond_ppt_gt_rad_adj = hru_ppt > ppt_rad_adj_day
        cond_if = cond_ppt_gt_rad_adj
        wh_if = np.where(cond_if)
        if len(wh_if[0]):

            # * if.else
            pptadj[wh_if] = (
                radadj_intcp_day
                + radadj_slope_day * (tmax_hru - tmax_index_day)
            )[wh_if]

            # * if.else.if
            cond_ppt_adj_gt_one = pptadj > one
            cond_if_else_if = cond_if & cond_ppt_adj_gt_one
            wh_if_else_if = np.where(cond_if_else_if)
            if len(wh_if_else_if[0]):
                pptadj[wh_if_else_if] = one
            del radadj_intcp_day, radadj_slope_day, tmax_index_day

            # * if.if
            tmax_index_day = monthly_to_daily(tmax_index)
            cond_tmax_lt_index = tmax_hru < tmax_index_day
            cond_if_if = cond_if & cond_tmax_lt_index
            wh_if_if = np.where(cond_if_if)
            if len(wh_if_if[0]):

                # The logic equiv to but changed from the original
                radj_sppt_day = tile_space_to_time(radj_sppt, n_time)
                pptadj[wh_if_if] = radj_sppt_day[wh_if_if]

                # if.if.else: negate the if.if.if
                radj_wppt_day = tile_space_to_time(radj_wppt, n_time)
                tmax_allrain_day = monthly_to_daily(
                    tmax_allrain_offset + tmax_allsnow
                )
                cond_tmax_lt_allrain = tmax_hru < tmax_allrain_day
                cond_if_if_else = cond_if_if & cond_tmax_lt_allrain
                wh_if_if_else = np.where(cond_if_if_else)
                if len(wh_if_if_else[0]):
                    pptadj[wh_if_if_else] = radj_wppt_day[wh_if_if_else]

                # if.if.if(.if): the cake is taken!
                cond_tmax_gt_allrain_and_not_summer = (
                    tmax_hru >= tmax_allrain_day
                ) & (~is_summer)
                cond_if_if_if = (
                    cond_if_if & cond_tmax_gt_allrain_and_not_summer
                )
                wh_if_if_if = np.where(cond_if_if_if)
                if len(wh_if_if_if[0]):
                    pptadj[wh_if_if_if] = radj_wppt_day[wh_if_if_if]

        radadj = radadj * pptadj
        wh_radadj_lt_2_tenths = np.where(radadj < 0.2)
        if len(wh_radadj_lt_2_tenths[0]):
            radadj[wh_radadj_lt_2_tenths] = 0.2

        swrad = doy_to_daily(soltab_potsw) * radadj / hru_cossl
        orad_hru = radadj * doy_to_daily(soltab_horad_potsw)
        return swrad, orad_hru

    # https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_potet_jh.f90

    def calculate_potential_et_jh(self) -> None:
        """Calculate potential evapotranspiration following Jensen and Haise
        (1963)."""

        self.potet.data[:] = self._potet_jh_run(
            dates=self._time,
            tavgc=self.tavgc.data,
            swrad=self.swrad.data,
            jh_coef=self.jh_coef,
            jh_coef_hru=self.jh_coef_hru,
            nmonth=self.nmonth,
        )
        return

    @staticmethod
    def _potet_jh_run(dates, tavgc, swrad, jh_coef, jh_coef_hru, nmonth):
        """This code mixes celcius and fahrenheit units in its calculation"""

        n_time, n_hru = tavgc.shape

        # The vector
        month = (dates.astype("datetime64[M]").astype(int) % nmonth) + 1
        month = tile_time_to_space(month, n_hru)
        tavgf = (tavgc * 9 / 5) + 32

        elh = (597.3 - (0.5653 * tavgc)) * inch2cm

        # JLM: move this out?
        def monthly_to_daily(arr):
            return arr[month[:, 0] - 1, :]

        jh_coef_day = monthly_to_daily(jh_coef)
        jh_coef_hru_day = tile_space_to_time(jh_coef_hru, n_time)

        potet = jh_coef_day * (tavgf - jh_coef_hru_day) * swrad / elh
        cond_potet_lt_zero = potet < zero
        wh_cond = np.where(cond_potet_lt_zero)
        if len(wh_cond[0]):
            potet[wh_cond] = zero

        return potet

    # # Track the amount of potential ET used at a given timestep
    # # JLM: is this strange to track here? I suppose not.
    # def consume_pot_et(self, requested_et):
    #     et = requested_et
    #     available_et = self.pot_et_current - self.pot_et_consumed
    #     if et > available_et:
    #         et = available_et
    #     self.pot_et_consumed += et
    #     return et

    def calculate_transp_tindex(self):

        # INIT: Process_flag==INIT
        # transp_on inited to 0 everywhere above

        # candidate for worst code lines
        if self.control.params.get_parameters("temp_units")["temp_units"] == 0:
            transp_tmax_f = self.transp_tmax
        else:
            transp_tmax_f = (self.transp_tmax * (9.0 / 5.0)) + 32.0

        transp_check = self.transp_on.current.copy()  # dim nhrus only
        tmax_sum = self.transp_on.current.copy().astype(
            "float64"
        )  # dim nhrus only
        start_day = self.control.start_doy
        start_month = self.control.start_month

        motmp = start_month + self.nmonth

        for hh in range(self.nhru):
            if start_month == self.transp_beg[hh]:
                # rsr, why 10? if transp_tmax < 300, should be < 10
                if start_day > 10:
                    self.transp_on.data[0, hh] = 1
                else:
                    transp_check[hh] = 1

            elif self.transp_end[hh] > self.transp_beg[hh]:
                if (start_month > self.transp_beg[hh]) and (
                    start_month < self.transp_end[hh]
                ):
                    self.transp_on.data[0, hh] = 1
            else:
                if (start_month > self.transp_beg[hh]) or (
                    motmp < self.transp_end[hh] + self.nmonth
                ):
                    self.transp_on.data[0, hh] = 1

        # vectorize
        ntime = self.transp_on.data.shape[0]
        # _transp_check = self.transp_on.data.copy()
        # self._transp_beg = tile_space_to_time(self.transp_beg, ntime)
        # self._transp_end = tile_space_to_time(self.transp_end, ntime)

        # RUN: Process_flag == RUN
        # Set switch for active transpiration period
        for tt in range(ntime):
            for hh in range(self.nhru):
                if tt > 0:
                    self.transp_on.data[tt, hh] = self.transp_on.data[
                        tt - 1, hh
                    ]

                # if tt == 2 and hh == 980:
                #      asdf

                # check for month to turn check switch on or
                # transpiration switch off

                if self._doy[tt] == 1:
                    # check for end of period
                    if self._month[tt] == self.transp_end[hh]:
                        self.transp_on.data[tt, hh] = 0
                        transp_check[hh] = 0
                        tmax_sum[hh] = zero

                    # <
                    # check for month to turn transpiration switch on or off
                    if self._month[tt] == self.transp_beg[hh]:
                        transp_check[hh] = 1
                        tmax_sum[hh] = zero

                # <<
                # If in checking period, then for each day
                # sum maximum temperature until greater than temperature index
                # parameter, at which time, turn transpiration switch on, check
                # switch off freezing temperature assumed to be 32 degrees
                # Fahrenheit
                if transp_check[hh] == 1:
                    if self.tmaxf.data[tt, hh] > 32.0:
                        tmax_sum[hh] = tmax_sum[hh] + self.tmaxf.data[tt, hh]

                    # <
                    if tmax_sum[hh] > transp_tmax_f[hh]:
                        self.transp_on.data[tt, hh] = 1
                        transp_check[hh] = 0
                        tmax_sum[hh] = 0.0

        # <<<
        return

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
            nc.add_all_data(var, self[var].data, self._time)
            nc.close()
            assert nc_path.exists()
            print(f"Wrote preprocessed forcing file: {nc_path}")

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
