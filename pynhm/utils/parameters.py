import pathlib as pl
from typing import Union

import numpy as np

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
                        if not found_dim:
                            temp_dims.append("unknown")
                    parameter_dimensions_dict[key] = temp_dims

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
    def nhm_coordinate(self) -> np.ndarray:
        """Get the nhm coordinate

        Returns:
            id: nhm coordinate for each hru

        """
        if "nhm_id" in self.parameters.keys():
            id = self.parameters["nhm_id"]
        else:
            id = np.arange(1, self.parameters["nhru"] + 1)
        return id

    @staticmethod
    def load(parameter_file: fileish) -> "PrmsParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        data = PrmsFile(parameter_file, "parameter").get_data()
        return PrmsParameters(
            data["parameter"]["parameters"],
            data["parameter"]["parameter_dimensions"],
        )
