import pathlib as pl
from typing import Union
# import xarray as xr

from .utils import cbh_file_to_df, cbh_files_to_df

fileish = Union[str, pl.PosixPath, dict]


class CBH:
    def __init__(
            self,
            files: fileish,
            convert_units: bool = False,
            adjust: bool = False,
            output_file: fileish = None,
            verbosity: bool = 0) -> None:

        self.files = files
        self.convert_units = convert_units
        self.adjust = adjust
        self.output_file = output_file
        self.verbosity = verbosity
        self.df = None

        self.to_df()

        if convert_units:
            self.convert_units()

        if adjust:
            self.adjust()

        if output_file is not None:
            self.write(output_file)

        return None

    def to_df(self):
        if isinstance(self.files, (str, pl.PosixPath)):
            self.df = cbh_file_to_df(self.files)
        elif isinstance(self.files, (dict)):
            self.df = cbh_files_to_df(self.files)
        else:
            raise ValueError(f'"files" argument of type {type(self.files)} not accepted.')

    def convert_units_cbh(self):
        pass

    def adjust_cbh(self):
        pass

    def check_cbh(self):
        pass

    def to_netcdf(self, nc_file):
        pass
