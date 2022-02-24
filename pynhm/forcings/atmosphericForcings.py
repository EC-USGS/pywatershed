try:
    import numpy as np
except ModuleNotFoundError:
    np = None

# JLM: "front load" option vs "load as you go"


class AtmosphericForcings:
    def __init__(self, forcings_file, verbose=0):

        # netcdf open
        # Set state from the netcdf file
        # get the adjusted states? depends on how these are written out
        if forcings_file is None:
            forcing_df = _read_cbh_individual_new(
                input_data_dir=input_data_dir,
                precip_file=precip_file,
                temp_min_file=temp_min_file,
                temp_max_file=temp_max_file,
                verbose=verbose,
            )
        else:
            raise NotImplementedError

        # self.datetime =
        # self.prcp =
        # self.rainfall =
        # self.snowfall =
        # self.tmax =
        # self.tmin =

        #
        # self.pot_et = pot_et
        self.pot_et_consumed = None
        self.current_date = None

    def adjust(self):
        msg = (
            "Base AtmosphericForcings class does"
            "not provide forcing adjustment"
        )
        raise NotImplementedError(msg)

    def advance(self, itime_step, current_date):
        self.precip_current = self.precip[itime_step]
        self.pot_et_current = self.pot_et[itime_step]
        self.pot_et_consumed = 0.0
        self.current_date = current_date

    # Track the amount of potential ET used at a given timestep
    # JLM: is this strange to track here? I suppose not.
    def consume_pot_et(self, requested_et):
        et = requested_et
        available_et = self.pot_et_current - self.pot_et_consumed
        if et > available_et:
            et = available_et
        self.pot_et_consumed += et
        return et


# JLM: There could be alternative adjustment approaches
# JLM: subclassing is one way to handle. Could also use decorators or just functions
# JLM: In the case of NHM, seems like we just want a particular subclass
class AtmForcingsNHM(AtmosphericForcings):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def adjust(self):
        self._adjust_temp()
        self._adjust_precip()

    # Compute the adjusted HRU temperature
    # Need parameters here
    def _adjust_temp(dates, tmax, tmin, tmaxadj, tminadj):
        tmax_hru = np.zeros(len(dates))
        tmin_hru = np.zeros(len(dates))
        ii = 0
        # JLM can date->month  be "vectorized"?
        for date in dates:
            # jday = int(date.strftime('%j'))
            imon = date.month - 1
            tmax_hru[ii] = tmax[ii] + tmaxadj[imon]
            tmin_hru[ii] = tmin[ii] + tminadj[imon]
            ii += 1

        return tmin_hru, tmax_hru


# https://github.com/nhm-usgs/prms/blob/6.0.0_dev/src/prmslib/physics/sm_potet_jh.f90
# is this part of the full class?
# # potet_jh run code:
# INCH2CM = 2.54

# def potet_jh_run(dates, tmax_hru, tmin_hru, swrad, jh_coef, jh_coef_hru):
#     potet = np.zeros(len(dates))

#     ihru = 0
#     ii = 0
#     for date in dates:
#         imon = date.month - 1
#         tavgf = (tmax_hru[ii] + tmin_hru[ii]) / 2.0
#         tavgc = (tavgf - 32.0) * 5 / 9
#         elh = (597.3 - (0.5653 * tavgc)) * INCH2CM
#         potet[ii] = jh_coef[imon, ihru] * (tavgf - jh_coef_hru[ihru]) * swrad[ii] / elh
#         if potet[ii] < 0.0:
#             potet[ii] = 0.0
#         ii += 1
#     return potet
