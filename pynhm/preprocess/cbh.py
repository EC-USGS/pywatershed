import pathlib as pl
from typing import Union
# import xarray as xr

from .utils import cbh_file_to_df, cbh_files_to_df, cbh_adjust
from pynhm import PrmsParameters

fileish = Union[str, pl.PosixPath, dict]

# JLM: Should also take a file describing the HRUs


class CBH:
    def __init__(
            self,
            files: fileish,
            params: PrmsParameters = False,
            convert_units: bool = False,
            adjust: bool = False,
            output_file: fileish = None,
            verbosity: bool = 0) -> None:

        self.files = files
        self.convert_units_flag = convert_units
        self.adjust_flag = adjust
        self.output_file = output_file
        self.verbosity = verbosity
        self.df = None

        self.to_df()
        self.check()

        if self.convert_units_flag:
            self.convert_units()

        if self.adjust_flag:
            self.adjust(params)
            self.check()

        if self.output_file is not None:
            self.write(output_file)

        return None

    def to_df(self):
        if isinstance(self.files, (str, pl.PosixPath)):
            self.df = cbh_file_to_df(self.files)
        elif isinstance(self.files, (dict)):
            self.df = cbh_files_to_df(self.files)
        else:
            raise ValueError(f'"files" argument of type {type(self.files)} not accepted.')

    def convert_units(self):
        pass

    def adjust(self, parameters):
        self.df = cbh_adjust(self.df, parameters)

    def check(self):
        pass

    def to_netcdf(self, nc_file):
        pass
