from io import StringIO
from functools import reduce
import numpy as np
import os
import pandas as pd
from .utils import timer
from .prms5util import unit_conversion

# JLM: "front load" option vs "load as you go"

class AtmosphericForcings:
    def __init__(
            self,
            #precip, tmax, tmin, pot_et,  # This will always be init'd from a file
            input_data_dir,
            forcings_file=None,
            precip_file=None,
            temp_min_file=None,
            temp_max_file=None,
            convert=True,
            verbose=False):

        if forcings_file is None:
            forcing_df = _read_cbh_individual_new(
                input_data_dir=input_data_dir,
                precip_file=precip_file,
                temp_min_file=temp_min_file,
                temp_max_file=temp_max_file,
                convert=convert, verbose=verbose)
        else:
            raise NotImplementedError

        self.datetime = forcing_df.precipitation.reset_index().date.values
        self.precipitation = forcing_df.precipitation.values
        self.temp_max = forcing_df.temp_max.values
        self.temp_min = forcing_df.temp_min.values
        del forcing_df

        #
        # self.pot_et = pot_et
        self.pot_et_consumed = None
        self.current_date = None

    def adjust(self):
        msg = ("Base AtmosphericForcings class does"
               "not provide forcing adjustment")
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


@timer
def _read_cbh_individual(
        input_data_dir,
        precip_file,
        temp_min_file,
        temp_max_file,
        convert=True, verbose=False):

    col_file_dict = {
        'precipitation': precip_file,
        'temp_min': temp_min_file,
        'temp_max': temp_max_file}

    templist = []
    for colname, filename in col_file_dict.items():
        fpath = os.path.join(input_data_dir, filename)
        print(f"Loading {fpath}")
        filelist = f"date,{colname}\n"
        with open(fpath) as f:
            for i, line in enumerate(f):
                if i > 2:
                    ll = line.strip().split()
                    yr = int(ll[0])
                    mo = int(ll[1])
                    da = int(ll[2])
                    data = float(ll[-1])
                    # Note: A fast-path exists for iso8601-formatted dates.
                    # filelist += f"{da:02d}-{mo:02d}-{yr:04d}, {data}\n"
                    filelist += f"{yr:04d}-{mo:02d}-{da:02d}, {data}\n"
        tdf = pd.read_csv(
            StringIO(filelist),
            parse_dates=["date"],
            index_col=["date"],
            dtype=float,
            float_precision="high",
        )
        templist.append(tdf)

    df = pd.concat([v for v in templist], axis=1)

    if convert:
        unit_conversion(df, verbose=verbose)

    return df


@timer
def _read_cbh_individual_new(
        input_data_dir,
        precip_file,
        temp_min_file,
        temp_max_file,
        convert=True,
        verbose=False):

    col_file_dict = {
        'precipitation': precip_file,
        'temp_min': temp_min_file,
        'temp_max': temp_max_file}

    df_list = []
    for colname, filename in col_file_dict.items():
        fpath = os.path.join(input_data_dir, filename)
        print(f"Loading {fpath}")

        # JLM: These files appear inconsistent across prcp vs others and also
        # markstroms code didnt handle this extra column in the temp files
        if colname == 'precipitation':
            cols = ["Y", "m", "d", "H", "M", "S", colname]
        else:
            cols = ["Y", "m", "d", "H", "M", "S", '-', colname]

        dtype_dict = dict()
        for cc in cols:
            type = 'str'
            if cc == colname:
                type = 'float64'
            dtype_dict[cc] = type

        data = pd.read_csv(
            fpath, names=cols, skiprows=3, delim_whitespace=True, dtype=dtype_dict)
        data['date'] = pd.to_datetime(
            data.Y + '-'
            + data.m.str.zfill(2) + '-'
            + data.d.str.zfill(2))
        data = data.drop(columns=set(cols).difference(set([colname])))
        data = data.set_index('date')
        df_list += [data]

    df = pd.concat(df_list, axis=1)

    if convert:
        unit_conversion(df, verbose=verbose)

    return df


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
