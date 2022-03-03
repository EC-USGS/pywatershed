import pathlib as pl
from copy import deepcopy
from typing import Union

from pynhm.utils import PrmsParameters

from .cbh_utils import (
    cbh_adjust,
    cbh_check,
    cbh_files_to_np_dict,
    cbh_n_hru,
    cbh_n_time,
    cbh_to_netcdf,
)

fileish = Union[str, pl.PosixPath, dict]

# Notes:
# * Could/Should also take a file describing the HRUs by ID (to go in that dim in the netcdf)
# * CBH dosent have a precisesly defined state, so use a dictonary for state.
# * I've relegated the code to cbh_utils becase there is potential reuse by the
#   AtmForcings object in model.


class CBH:
    def __init__(
        self,
        files: fileish,
        parameters: PrmsParameters = None,
        adjust: bool = False,
        units: dict = None,
        new_units: dict = None,
        output_vars: list = None,
        output_file: fileish = None,
        verbosity: bool = 0,
    ) -> None:

        self.files = files
        self.params = parameters
        self.adjust = adjust
        self.units = units
        self.new_units = new_units
        self.state = None
        self.output_vars = output_vars
        self.output_file = output_file
        self.verbosity = verbosity

        self.set_state()
        # option to convert state to pd or xr for plotting capabilities (prior to nc write)?

        if self.new_units is not None:
            self.convert_units()

        if (self.adjust) and (self.params is not None):
            self._adjust()

        if self.output_file is not None:
            self.to_netcdf()

        return None

    # Is this private?
    # other setters
    def set_state(self):
        self.state = cbh_files_to_np_dict(self.files, self.params)
        self.check()
        return

    # @property
    # def get_df
    # @property
    # def get_variable_names
    # @property
    # def get_variable

    @property
    def n_rhu(self):
        return cbh_n_hru(self.state)

    @property
    def n_time(self):
        return cbh_n_time(self.state)

    def convert_units(self):
        pass

    def _adjust(self):
        _ = cbh_adjust(self.state, self.params)
        # self.check()  # JLM: put back, but this is also wrong in the prms/nhm outputs
        return

    def check(self):
        _ = cbh_check(self.state, verbosity=self.verbosity)
        return

    def to_netcdf(self, *args, **kwargs):
        _ = cbh_to_netcdf(self.state, *args, **kwargs)
        pass
