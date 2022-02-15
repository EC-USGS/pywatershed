import pathlib as pl
from typing import Union
# import xarray as xr

from .cbh_utils import cbh_files_to_np_dict, cbh_adjust
from pynhm import PrmsParameters

fileish = Union[str, pl.PosixPath, dict]

# JLM: Should also take a file describing the HRUs
# CBH dosent have a precisesly defined sate, so use a dictonary for state.

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
        self.state = None

        self.set_state()
        self.check()

        if self.convert_units_flag:
            self.convert_units()

        if self.adjust_flag:
            self.adjust(params)
            self.check()

        if self.output_file is not None:
            self.write(output_file)

        return None

    def set_state(self):
        self.state = cbh_files_to_np_dict(self.files)
        return

    def convert_units(self):
        pass

    def adjust(self, parameters):
        self.df = cbh_adjust(self.df, parameters)

    def check(self):
        pass

    def to_netcdf(self, nc_file):
        pass
