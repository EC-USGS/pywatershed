import pathlib as pl
from typing import Union

import numpy as np

from .dictionary_as_properties import DictionaryAsProperties
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

        self.parameters = DictionaryAsProperties(parameter_dict)

        # build dimensions from data
        if parameter_dimensions_dict is None:
            dimensions = self.get_dimensions
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

        self.parameter_dimensions = DictionaryAsProperties(
            parameter_dimensions_dict
        )

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
    def get_dimensions(self) -> dict:
        """Get the dimensions from the parameters

        Returns:
            dimensions in the PRMS parameter dictionary

        """
        dimensions = {}
        for key, value in self.parameters.items():
            if isinstance(value, int):
                dimensions[key] = value
        return DictionaryAsProperties(dimensions)

    @staticmethod
    def load(parameter_file: fileish) -> "PrmsParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        # (
        #     dimensions,
        #     parameter_data,
        #     parameter_dimensions,
        #     parameter_types,
        # ) = _load_prms_parameters(parameter_file)
        # parameters = parameter_data.copy()
        # for key, value in dimensions.items():
        #     parameters[key] = value
        data = PrmsFile(parameter_file, "parameter").get_data()
        return PrmsParameters(
            data["parameter"]["parameters"],
            data["parameter"]["parameter_dimensions"],
        )


def _load_prms_parameters(parameter_file: fileish):
    """Read a PRMS parameter file

    :param parameter_file:
    :return:
    """
    line_num = 0
    vals = {}
    dims = {}
    param_dims = {}
    param_type = {}

    with open(parameter_file) as f:
        reading_dims = False
        for line in f:
            try:
                line = line.rstrip()  # remove '\n' at end of line
                line_num += 1
                if line == "** Dimensions **":
                    reading_dims = True
                    line = f.readline().rstrip()
                    line_num += 1

                if line == "** Parameters **":
                    reading_dims = False
                    break

                if reading_dims:
                    line = f.readline().rstrip()
                    line_num += 1
                    dim_name = line

                    line = f.readline().rstrip()
                    line_num += 1
                    size = line

                    if dim_name in dims.keys():
                        pass
                    else:
                        dims[dim_name] = int(size)
            except:
                msg = (
                    f"read parameters exception line = {line}\n"
                    + f"read parameters exception line_num = {str(line_num)}\n"
                )
                raise ValueError(msg)

        #        read params
        for line in f:
            try:
                line = line.rstrip()  # remove '\n' at end of line
                line_num += 1

                if line == "####":
                    line = f.readline().rstrip()
                    line = line.split(" ", 1)[
                        0
                    ]  # old format parameter files have a blank (' ') and then a width format value. Strip this off.
                    param_name = line
                    line_num += 1

                    line = f.readline().rstrip()
                    line_num += 1
                    num_dims = int(line)
                    pd = [None] * num_dims
                    for ii in range(num_dims):
                        line = f.readline().rstrip()
                        pd[ii] = line
                        line_num += 1

                    param_dims[param_name] = pd

                    line = f.readline().rstrip()
                    line_num += 1
                    num_vals = int(line)
                    line = f.readline().rstrip()
                    line_num += 1
                    tp = int(line)
                    param_type[param_name] = tp

                    if tp == 2:
                        vs = np.zeros(num_vals, dtype=float)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = float(line)

                    elif tp == 1:
                        vs = np.zeros(num_vals, dtype=int)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = int(line)

                    else:
                        vs = np.zeros(num_vals, dtype=np.chararray)
                        for jj in range(num_vals):
                            line = f.readline().rstrip()
                            line_num += 1
                            vs[jj] = line

                    if num_dims == 2:
                        vs.shape = (dims[pd[1]], dims[pd[0]])

                    if param_name in vals.keys():
                        print(
                            "parameter ",
                            param_name,
                            " is already in ",
                            parameter_file,
                        )
                    else:
                        vals[param_name] = vs

            except:
                raise ValueError(
                    f"read parameters exception line_num = {line_num}"
                )

    return dims, vals, param_dims, param_type
