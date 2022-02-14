import math
import pandas as pd
import pathlib as pl
from typing import Union
import xarray as xr

fileish = Union[str, pl.PosixPath, dict]


class CBH:
    def __init__(
            self,
            files: fileish,
            convert_units: bool = False,
            adjust: bool = False,
            output_file: fileish = None,
            verbosity: bool = 0) -> None:

        # Data -----
        self.files = files
        self.convert_units = convert_units
        self.adjust = adjust
        self.output_file = output_file
        self.verbosity = verbosity
        self.df = None

        # Methods -----
        # self.to_csv = _df_to_csv  # pre-deprecated?
        self.convert_units_cbh = _convert_units_cbh
        self.adjust_cbh = _adjust_cbh
        self.check_cbh = _check_cbh

        # How to get a data frame depends on the files argument
        if isinstance(self.files, (str, pl.PosixPath)):
            self.to_df = _cbh_file_to_df
        elif isinstance(self.files, (dict)):
            self.to_df = _cbh_files_to_df
        else:
            raise ValueError(f'"files" argument of type {type(files)} not accepted.')

        # Action -----
        self.df = self.to_df(self.files)  # do I want to do this at init?

        # Special action -----
        # If the output_file is specified: get it all done
        # (if you dont want the object and you just want the written file side-effect)
        if output_file is not None:
            if convert_units:
                self.df = self.convert_units(self.df)
            if adjust:
                self.df = self._
                # print the path that was written

    # Method
    def to_netcdf(self, nc_file):
        return _df_to_netcdf(self.df, nc_file)


        # self.checks

        # accept output file to just write the converted files?

        return None


def _df_to_netcdf(df, file):
    asdf
    pass


def _convert_units_cbh():
    pass


def _adjust_cbh():
    pass


def _check_cbh():
    pass


created_line = 'Created '
written_line = 'Written '
hash_line = '##############'
hash_line_official = '########################################'


def _cbh_file_to_df(the_file):
    # This is trying to handle as many input formats for these kinds of files
    # as we can find. It may not be comprehensive. See tests for what is currently
    # handled
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
    for posn in range(len(meta_lines)):
        key, count = meta_lines[posn].split(' ')
        count = int(count)
        zs = math.ceil(math.log(count, 10))
        col_names += [f"{key}{str(ii).zfill(zs)}" for ii in range(count)]

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

    return data


def _cbh_files_to_df(file_dict):
    dfs = [_cbh_file_to_df(val) for val in file_dict.values()]
    return pd.concat(dfs, axis=1)
