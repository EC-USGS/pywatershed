import math
import pandas as pd
import pathlib as pl
from pynhm import PrmsParameters
from typing import Union

file_type = Union[str, pl.PosixPath]

created_line = 'Created '
written_line = 'Written '
hash_line = '##############'
hash_line_official = '########################################'

cbh_std_name_from_dict = {
    'prcp': ['precip', 'prcp'],
    'tmax': ['tmax'],
    'tmin': ['tmin'],
    'orad': ['orad'],
    'ptet': ['ptet'],
    'runoff': ['runoff']
}

required_vars = ['precp', 'tmax', 'tmin']


def get_std_name(name: str) -> str:
    for key, val in cbh_std_name_from_dict.items():
        if name in val:
            return key
    msg = f"The variable name '{name}' can not be mapped to a standard name."
    raise ValueError(msg)


def cbh_file_to_df(the_file:file_type) -> pd.DataFrame:
    # This is attempting to handle as many input formats for these kinds of files
    # as we can find. It may not be comprehensive. See tests for what is currently
    # handled
    # Handle netcdf?
    meta_lines = []
    with open(the_file, 'r') as file_open:
        wh_hash_line = -1
        the_line = ''
        while hash_line not in the_line:
            wh_hash_line += 1
            the_line = file_open.readline()
            if (
                    ('//' not in the_line)
                    and (created_line not in the_line)
                    and (written_line not in the_line)
                    and (hash_line not in the_line)):
                meta_lines += [the_line.strip()]

    col_names = []
    var_count_dict = {}
    for posn in range(len(meta_lines)):
        key, count = meta_lines[posn].split(' ')
        count = int(count)
        zs = math.ceil(math.log(count, 10))
        std_name = get_std_name(key)
        col_names += [f"{std_name}{str(ii).zfill(zs)}" for ii in range(count)]
        # col_names += [f"{key}{str(ii).zfill(zs)}" for ii in range(count)]        
        var_count_dict[key] = count

    # Can we get the hru info? is that standarized at all?

    dtypes = (['str'] * 6) + (['float64'] * len(col_names))
    date_cols = ["Y", "m", "d", "H", "M", "S"]
    col_names = date_cols + col_names
    assert len(dtypes) == len(col_names)
    dtype_dict = dict(zip(col_names, dtypes))
    assert len(dtype_dict) == len(col_names)

    # Some files erroneously specify the number of columns
    # catch that here
    data = pd.read_csv(
        the_file,
        nrows=1,
        # names=col_names,  # dont specify names
        header=None,  # dont use default header from first line
        index_col=False,
        skiprows=wh_hash_line+1,
        delim_whitespace=True, dtype=dtype_dict)
    msg = f"Number of actual data columns does not match metadata info: {meta_lines}"
    assert len(data.columns) == len(col_names), msg
    # JLM: is the above sufficient?

    data = pd.read_csv(
        the_file,
        names=col_names,
        index_col=False,
        skiprows=wh_hash_line+1,
        delim_whitespace=True, dtype=dtype_dict)

    data['date'] = pd.to_datetime(
        data.Y + '-'
        + data.m.str.zfill(2) + '-'
        + data.d.str.zfill(2))
    data = data.drop(columns=set(date_cols))
    data = data.set_index('date')

    # Convert to a multi-index? For getting to a dict of np arrays
    # This might go elsewhere
    def col_name_split(string):
        char_list = list(string)
        wh_digit = [ii for ii in range(len(char_list)) if char_list[ii].isdigit()]
        return string[:wh_digit[0]], string[wh_digit[0]:]
    col_name_tuples = [col_name_split(kk) for kk,vv in data.iteritems()]
    # data.columns = pd.MultiIndex.from_tuples(col_name_tuples)

    return data


def cbh_files_to_df(file_dict):
    dfs = [cbh_file_to_df(val) for val in file_dict.values()]
    return pd.concat(dfs, axis=1)


def cbh_adjust(df, params: PrmsParameters):
    # adjust temp
    param_data = params._parameter_data
    nhru = params._dimensions['nhru']

    tmax_allsnow = param_data["tmax_allsnow"]
    snow_cbh_adj = param_data["snow_cbh_adj"]
    rain_cbh_adj = param_data["rain_cbh_adj"]
    tmax_cbh_adj = param_data["tmax_cbh_adj"]
    tmin_cbh_adj = param_data["tmin_cbh_adj"]
    tmax_allrain_offset = param_data["tmax_allrain_offset"]
    tmax_allrain = tmax_allsnow + tmax_allrain_offset
    adjmix_rain = param_data["adjmix_rain"]

    # Have to handle this data having lens 1 or 12

    # markstrom code for adjusted temperatures: >>>
    # tmax_hru = np.zeros(len(dates))
    # tmin_hru = np.zeros(len(dates))

    ii = 0
    for date in dates:
        jday = day_of_year(date)
        imon = date.month - 1
        tmax_hru[ii] = tmax[ii] + tmax_cbh_adj[imon]
        tmin_hru[ii] = tmin[ii] + tmin_cbh_adj[imon]
        ii += 1
    # << markstrom code

    asdf

    # adjust precip

    pass


def cbh_convert_units():
    pass


def cbh_check():

    # Markstroms code
    # tiff = tmax - tmin
    # print("there are", np.sum(np.array(tiff)) < 0, "negative values")
    # print("there are", np.sum(np.array(tiff)) == 0, "equal values")

    # Parkers comments
    #    tmax is never less than tmin
    #    prcp is never negative
    #    any missing data/missing date is filled with (?? avg of bracketing dates??)


    pass
