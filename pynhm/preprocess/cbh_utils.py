import math
import numpy as np
import pandas as pd
import pathlib as pl
from pynhm import PrmsParameters
from typing import Union

zero = np.zeros((1))[0]
one = np.ones((1))[0]

# Compound types
file_type = Union[str, pl.PosixPath]
fileish = Union[str, pl.PosixPath, dict]

# Absorb some slop in the naming conventions
cbh_std_name_from_dict = {
    'prcp': ['precip', 'prcp'],
    'rhavg': ['rhavg'],
    'tmax': ['tmax'],
    'tmin': ['tmin'],
    #'orad': ['orad'],  # these are in some cbh for now only take NHM inputs above
    #'ptet': ['ptet'],
    #'runoff': ['runoff']
}

required_vars = ['precp', 'tmax', 'tmin']

cbh_metadata = {
    'time': {
        'long_name':'time',
        'standard_name':'time',
        'calendar':'standard',                       # Depends, revisit
        'units':'days since 1979-01-01 00:00:00',},  # Depends, may not be correct
    'hruid': {
        'long_name':'Hydrologic Response Unit ID (HRU)',
        'cf_role':'timeseries_id',},
    'tmax': {
        '_FillValue':9.96921e36,
        'long_name':'Maximum daily air temperature',
        'units':'degree_fahrenheit',
        'standard_name':'maximum_daily_air_temperature',},
    'tmin': {
        '_FillValue':9.96921e36,
        'long_name':'Minimum daily air temperature',
        'units':'degree_fahrenheit',
        'standard_name':'minimum_daily_air_temperature',},
    'prcp': {
        '_FillValue':9.96921e36,
        'long_name':'daily total precipitation',
        'units':'in',
        'standard_name':'daily_total_precipitation',},
    'rhavg': {
        '_FillValue':9.96921e36,
        'long_name':'Daily mean relative humidity',
        'units':'percent',
        'standard_name':'rhavg',}}

cbh_units = {key:val['units'] for key, val in cbh_metadata.items() if 'units' in list(val.keys())}

created_line = 'Created '
written_line = 'Written '
hash_line = '##############'
hash_line_official = '########################################'


def _get_std_name(name: str) -> str:
    for key, val in cbh_std_name_from_dict.items():
        if name in val:
            return key
    msg = f"The variable name '{name}' can not be mapped to a standard name."
    raise ValueError(msg)


def _cbh_file_to_df(the_file: file_type) -> pd.DataFrame:
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
        std_name = _get_std_name(key)
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
    # JLM TODO: Set datetime resolution to hours? or mins. Could do days but might look forward a bit.
    data = data.drop(columns=set(date_cols))
    data = data.set_index('date')

    return data


def _cbh_files_to_df(file_dict: dict) -> pd.DataFrame:
    dfs = [_cbh_file_to_df(val) for val in file_dict.values()]
    return pd.concat(dfs, axis=1)


def cbh_files_to_df(files: fileish) -> pd.DataFrame:
    if isinstance(files, (str, pl.PosixPath)):
        df = _cbh_file_to_df(files)
    elif isinstance(files, (dict)):
        df = _cbh_files_to_df(files)
    else:
        raise ValueError(f'"files" argument of type {type(files)} not accepted.')
    return df


def _col_name_split(string: str) -> tuple:
    char_list = list(string)
    wh_digit = [ii for ii in range(len(char_list)) if char_list[ii].isdigit()]
    return string[:wh_digit[0]], string[wh_digit[0]:]


def cbh_df_to_np_dict(df: pd.DataFrame) -> dict:
    # Convert to a multi-index? For getting to a dict of np arrays
    # This might go elsewhere
    col_name_tuples = [_col_name_split(kk) for kk,vv in df.iteritems()]
    df.columns = pd.MultiIndex.from_tuples(col_name_tuples)
    var_names = df.columns.unique(level=0)
    np_dict = {}
    np_dict['datetime'] = df.index.to_numpy(copy=True).astype('datetime64[D]')
    for vv in var_names:
        np_dict[vv] = df[vv].to_numpy(copy=True)
    return np_dict


def cbh_files_to_np_dict(files: fileish) -> dict:
    np_dict = cbh_df_to_np_dict(cbh_files_to_df(files))
    return np_dict


def cbh_adjust(cbh_dict: dict, params: PrmsParameters) -> dict:
    # Param object has no defined interface at this time.
    param_data = params._parameter_data
    nhru = params._dimensions['nhru']

    # I dislike using pd for something that seems like it should exist in np
    month_ind_12 = pd.to_datetime(cbh_dict['datetime']).month - 1  # (time)
    month_ind_1 = np.zeros(cbh_dict['datetime'].shape, dtype=int) # (time)

    # adjust temp ---------------------------------
    # Temperature bias corrections are calibrated by PRMS/HNM
    # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_temperature_hru.f90#L91
    tmax_param = param_data["tmax_cbh_adj"]  # (12 or 1, space)
    tmin_param = param_data["tmin_cbh_adj"]

    # i suppose thse could have different shapes.
    # throw an error if that happens
    if tmax_param.shape != tmin_param.shape:
        msg = "Not implemented: tmin/tmax cbh adj parameters with different shapes"
        raise NotImplementedError(msg)

    if tmax_param.shape[0] == 12:
        month_ind = month_ind_12
    elif tmax_param.shape[0] == 1:
        month_ind = month_ind_1
    else:
        msg = "Unexpected month dimension for cbh temperature adjustment params"
        raise ValueError(msg)

    cbh_dict['tmax_adj'] = np.zeros(cbh_dict['tmax'].shape, dtype=cbh_dict['tmax'].dtype)  # (time, space)
    cbh_dict['tmin_adj'] = np.zeros(cbh_dict['tmin'].shape, dtype=cbh_dict['tmin'].dtype)
    cbh_dict['tmax_adj'] = cbh_dict['tmax'] + tmax_param[month_ind]
    cbh_dict['tmin_adj'] = cbh_dict['tmin'] + tmin_param[month_ind]

    # To check the resulting shape used above
    # tmax_time_params = tmax_param[month_ind]
    # for tt in range(len(month_ind)):
    #     assert (tmax_time_params[tt,:] == tmax_param[month_ind[tt]]).all()

    # adjust precip/snowfall ---------------------------------
    # Snow/rain partitioning of total precip is (adjusted) temperature dependent in addition
    # to being dependent on the following parameters (not sure if all are calibrated).
    # https://github.com/nhm-usgs/prms/blob/92f3c470bbf10e37ee23f015f60d42f6a028cf48/src/prmslib/physics/sm_precipitation.f90#L227
    tmax_allsnow_param = param_data["tmax_allsnow"]
    tmax_allrain_offset_param = param_data["tmax_allrain_offset"]
    snow_cbh_adj_param = param_data["snow_cbh_adj"]
    rain_cbh_adj_param = param_data["rain_cbh_adj"]
    adjmix_rain_param = param_data["adjmix_rain"]

    # i suppose thse could have different shapes.
    # throw an error if that happens
    shape_list = np.array(
        [tmax_allsnow_param.shape[0], tmax_allrain_offset_param.shape[0],
         snow_cbh_adj_param.shape[0], rain_cbh_adj_param.shape[0],
         adjmix_rain_param.shape[0]])
    if not (shape_list == 12).all():
        msg = "Not implemented: tmin/tmax cbh adj parameters with different shapes"
        raise NotImplementedError(msg)

    tmax_allrain_param = tmax_allsnow_param + tmax_allrain_offset_param

    if tmax_allsnow_param.shape[0] == 12:
        month_ind = month_ind_12
    elif tmax_allsnow_param.shape[0] == 1:
        month_ind = month_ind_1
    else:
        msg = "Unexpected month dimension for cbh precip adjustment params"
        raise ValueError(msg)

    cbh_dict['prcp_adj'] = np.zeros(cbh_dict['prcp'].shape, dtype=cbh_dict['prcp'].dtype)  # (time, space)
    cbh_dict['rainfall_adj'] = np.zeros(cbh_dict['prcp'].shape, dtype=cbh_dict['prcp'].dtype)
    cbh_dict['snowfall_adj'] = np.zeros(cbh_dict['prcp'].shape, dtype=cbh_dict['prcp'].dtype)

    # JLM: in my test vectorization was 45x faster for drb_2yr: 3.728s:.083s for loop:vectorized
    prmx = np.zeros(cbh_dict['prcp'].shape, dtype=cbh_dict['prcp'].dtype)

    # Calculate the mix everywhere, then set the precip/rain/snow amounts from the conditions.
    tdiff = cbh_dict['tmax_adj'] - cbh_dict['tmin_adj']
    prmx = ((cbh_dict['tmax_adj'] - tmax_allsnow_param[month_ind]) / tdiff) * adjmix_rain_param[month_ind]
    del tdiff

    wh_all_snow = np.where(cbh_dict['tmax_adj'] <= tmax_allsnow_param[month_ind])
    wh_all_rain = np.where(
        np.logical_or(
            cbh_dict['tmin_adj'] > tmax_allsnow_param[month_ind],
            cbh_dict['tmax_adj'] >= tmax_allrain_param[month_ind]))
    # This order MATTERS per the logic in PRMS: if(all_snow),elif(all_rain),else(mixed)
    prmx[wh_all_rain] = one
    prmx[wh_all_snow] = zero

    # Recalculate these
    wh_all_snow = np.where(prmx <= zero)
    wh_all_rain = np.where(prmx >= one)
    wh_mixed = np.where(np.logical_and(prmx < one, prmx > zero))

    # Mixed case (to be over written in the all snow/rain fall cases)
    wh_prmx_mixed = np.where(np.logical_and(prmx < one, prmx > zero))
    cbh_dict['prcp_adj'] = (cbh_dict['prcp'] * snow_cbh_adj_param[month_ind])
    cbh_dict['rainfall_adj'] = (prmx * cbh_dict['prcp_adj'])
    cbh_dict['snowfall_adj'] = (cbh_dict['prcp_adj'] - cbh_dict['rainfall_adj'])
    del prmx

    # All precip is snow case
    # The condition to be used later:
    all_snow_prcp = cbh_dict['prcp'] * snow_cbh_adj_param[month_ind]
    cbh_dict['prcp_adj'][wh_all_snow] = all_snow_prcp[wh_all_snow]
    cbh_dict['rainfall_adj'][wh_all_snow] = zero
    cbh_dict['snowfall_adj'][wh_all_snow] = all_snow_prcp[wh_all_snow]
    del all_snow_prcp

    # All precip is rain case
    # The condition to be used later:
    all_rain_prcp = cbh_dict['prcp'] * rain_cbh_adj_param[month_ind]
    cbh_dict['prcp_adj'][wh_all_rain] = all_rain_prcp[wh_all_rain]
    cbh_dict['rainfall_adj'][wh_all_rain] = all_rain_prcp[wh_all_rain]
    cbh_dict['snowfall_adj'][wh_all_rain] = zero
    del all_rain_prcp

    return None


def cbh_convert_units():
    pass


def cbh_check(cbh_dict: dict, verbosity: int = 0) -> None:

    assert (~np.isnat(cbh_dict['datetime'])).all()

    for adj in ['', '_adj']:

        # assume one variable represents if there are any adjustments
        if f'tmax{adj}' not in cbh_dict.keys():
            continue

        tmaxvar = f'tmax{adj}'
        tminvar = f'tmin{adj}'
        if not (cbh_dict[tmaxvar] > cbh_dict[tminvar]).all():
            msg = f"{tmaxvar} !> {tminvar}: strictly greater maybe too stringent"
            raise ValueError(msg)

        # assert (cbh_dict['tmax'] >= zero).all()  # use units for checking max/minimums?

        prcpvar = f'prcp{adj}'
        if not (cbh_dict[prcpvar] >= zero).all():
            msg = f'{prcpvar} contains negative values'
            raise ValueError(msg)

    if verbosity >= 1:
        print('cbh state check successful')

    return
