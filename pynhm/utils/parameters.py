import pathlib as pl
from typing import Union

import netCDF4 as nc4
import numpy as np
import pynhm

from ..constants import ft2_per_acre, inches_per_foot, ndoy, nmonth
from .prms5_file_util import PrmsFile

fileish = Union[str, pl.PosixPath, dict]
listish = Union[str, list, tuple]


class PrmsParameters:
    """
    PRMS parameter class

    Parameters
    ----------
    parameter_dict : dict
        parameters dictionary
    parameter_dimensions_dict : dict
        parameters dimensions dictionary

    """

    def __init__(
        self,
        parameter_dict: dict,
        parameter_dimensions_dict: dict = None,
    ) -> "PrmsParameters":

        self.parameters = parameter_dict

        # build dimensions from data
        # todo: this could be done from metadata. probably fewer lines of code
        #       and would ensure consistency
        if parameter_dimensions_dict is None:
            dimensions = self.dimensions
            parameter_dimensions_dict = {}
            for key, value in parameter_dict.items():
                if isinstance(value, int):
                    parameter_dimensions_dict[key] = None
                elif isinstance(value, np.ndarray):
                    shape = value.shape
                    temp_dims = []
                    for isize in shape:
                        found_dim = False
                        for dim_key, dim_value in dimensions.items():
                            if dim_value == isize:
                                found_dim = True
                                temp_dims.append(dim_key)
                                break
                        if isize == 1:
                            found_dim = True
                            temp_dims.append("scalar")
                        if not found_dim:
                            temp_dims.append("unknown")
                    parameter_dimensions_dict[key] = tuple(temp_dims)

        self.parameter_dimensions = parameter_dimensions_dict

    def get_parameters(self, keys: listish) -> "PrmsParameters":
        """Get a subset of keys in the parameter dictionary

        Args:
            keys: keys to retrieve from the full PRMS parameter object

        Returns:
            PrmsParameters : subset of full parameter dictionary
                Passed keys that do not exist in the full parameter
                dictionary are skipped.

        """
        if isinstance(keys, str):
            keys = [keys]

        return PrmsParameters(
            {
                key: self.parameters.get(key)
                for key in keys
                if key in self.parameters.keys()
            },
            {
                key: self.parameter_dimensions.get(key)
                for key in keys
                if key in self.parameter_dimensions.keys()
            },
        )

    @property
    def dimensions(self) -> dict:
        """Get the dimensions from the parameters

        Returns:
            dimensions in the PRMS parameter dictionary

        """
        dimensions = {}
        for key, value in self.parameters.items():
            if isinstance(value, int):
                dimensions[key] = value
        return dimensions

    @property
    def nhm_hru_coordinate(self) -> np.ndarray:
        """Get the nhm hru coordinate

        Returns:
            id: nhm coordinate for each hru

        """
        if "nhm_id" in self.parameters.keys():
            nhm_id = self.parameters["nhm_id"]
        else:
            nhm_id = np.arange(1, self.parameters["nhru"] + 1)
        return nhm_id

    @property
    def nhm_segment_coordinate(self) -> np.ndarray:
        """Get the nhm segment coordinate

        Returns:
            id: nhm coordinate for each segment

        """
        if "nhm_seg" in self.parameters.keys():
            nhm_seg = self.parameters["nhm_seg"]
        else:
            nhm_seg = np.arange(1, self.parameters["nsegment"] + 1)
        return nhm_seg

    @property
    def nhm_coordinates(self) -> dict:
        """Get the nhm coordinates

        Returns:
            id: nhm coordinates

        """
        return {
            "nhm_id": self.nhm_hru_coordinate,
            "nhm_seg": self.nhm_segment_coordinate,
        }

    @staticmethod
    def load(parameter_file: fileish) -> "PrmsParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        data = PrmsFile(parameter_file, "parameter").get_data()

        # could insert dimenion data here. going to add this one for now.
        # it's a little unclear if constants like this are parameters or not
        data["parameter"]["parameters"]["ndoy"] = ndoy
        data["parameter"]["parameter_dimensions"]["ndoy"] = None
        data["parameter"]["parameters"]["nmonth"] = 12
        data["parameter"]["parameter_dimensions"]["nmonth"] = None

        return PrmsParameters(
            data["parameter"]["parameters"],
            data["parameter"]["parameter_dimensions"],
        )

    @staticmethod
    def from_nc_file(parameter_nc_file: fileish) -> "PRMSParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        ds = nc4.Dataset(parameter_nc_file)
        param_dict = {kk: vv[:].data for kk, vv in ds.variables.items()}
        param_dict_dimensions = {
            kk: vv.dimensions for kk, vv in ds.variables.items()
        }

        scalar_params = {
            kk: len(vv) for kk, vv in ds.dimensions.items() if kk != "scalar"
        }
        scalar_param_dims = {kk: None for kk in scalar_params.keys()}

        param_dict = {**param_dict, **scalar_params}
        param_dict_dimensions = {**param_dict_dimensions, **scalar_param_dims}

        ds_attrs = ds.__dict__
        if ("nhm_process" in ds_attrs) and (
            ds_attrs["nhm_process"] == "PRMSSolarGeometry"
        ):
            param_dict["ndoy"] = ndoy
            param_dict_dimensions["ndoy"] = None

        return PrmsParameters(param_dict, param_dict_dimensions)

    def hru_in_to_cfs(self, time_step: np.timedelta64) -> np.ndarray:
        "Derived parameter converting inches to cfs on hrus."
        time_step_seconds = time_step / np.timedelta64(1, "s")
        return self.hru_in_to_cf / time_step_seconds

    @property
    def hru_in_to_cf(self) -> np.ndarray:
        "Derived parameter converting inches to cubic feet on hrus."
        return self.parameters["hru_area"] * ft2_per_acre / (inches_per_foot)
