import pathlib as pl
from enum import Enum
from typing import Union

import numpy as np

from ..base import meta

fileish = Union[str, pl.PosixPath]


def rename_dims(dim_name):
    if dim_name == "nmonths":
        return "nmonth"
    elif dim_name == "one":
        return "scalar"
    else:
        return dim_name


def expand_scalar_to_dims(param_dict, param_dim_dict):
    from pywatershed.utils.separate_nhm_params import (
        params_expand_scalar_to_dims,
    )

    # sometimes a scalar is allowed to represent a uniform values for
    # the full dimensions of a parameter. Going to handle those on a
    # case-by-case basis for now: just expand them to full size
    for param_name in param_dict.keys():
        if param_name in params_expand_scalar_to_dims:
            dims = meta.find_variables(param_name)[param_name]["dims"]
            exp_dim = params_expand_scalar_to_dims[param_name]
            param_shape = param_dict[param_name].shape
            param_val = param_dict[param_name]
            param_dims = param_dim_dict[param_name]
            for ii, dim_name in enumerate(dims):
                if dim_name == exp_dim:
                    param_shape = list(param_shape)
                    param_dims = list(param_dims)
                    param_shape[ii] = param_dict[exp_dim]
                    param_dims[ii] = exp_dim
                    param_shape = tuple(param_shape)
                    param_vals = np.zeros(param_shape) + param_val
                    param_dict[param_name] = param_vals
                    param_dim_dict[param_name] = tuple(param_dims)

    return param_dict, param_dim_dict


class PrmsDataType(Enum):
    INTEGER = 1
    FLOAT = 2
    CHARACTER = 4


class PrmsFileType(Enum):
    CONTROL = "control"
    PARAMETER = "parameter"


class PrmsFileSection(Enum):
    CONTROL = "control"
    DIMENSIONS = "dimensions"
    PARAMETER = "parameter"
    UNDEFINED = "undefined"


class PrmsFile:
    def __init__(
        self,
        file_path: fileish,
        file_type: str = None,
    ) -> "PrmsFile":
        self.file_path = file_path
        self.file_type = file_type
        self.file_object = None
        self.line_number = None
        self.section = PrmsFileSection.UNDEFINED
        self.eof = None
        self.set_file_type(file_type)
        self.dimensions = None

    def set_file_type(
        self,
        file_type: str = None,
    ) -> None:
        """

        Args:
            file_type:

        Returns:

        """
        self.line_number = None
        if file_type.lower() == PrmsFileType.CONTROL.value:
            self.file_type = PrmsFileType.CONTROL
            self.section = PrmsFileSection.CONTROL
        elif file_type.lower() == PrmsFileType.PARAMETER.value:
            self.file_type = PrmsFileType.PARAMETER
            self.section = PrmsFileSection.DIMENSIONS
        else:
            raise ValueError(
                f"file_type ({file_type}) "
                + f"must be '{PrmsFileType.CONTROL.value}' or "
                + f"'{PrmsFileType.PARAMETER.value}'"
            )

    def get_data(self) -> dict:
        self._get_file_object()
        if self.file_type == PrmsFileType.CONTROL:
            return {self.file_type.value: self._get_control_variables()}
        elif self.file_type == PrmsFileType.PARAMETER:
            (
                parameters,
                parameter_dimensions,
            ) = self._get_dimensions_parameters()
            return {
                self.file_type.value: {
                    "parameters": parameters,
                    "parameter_dimensions": parameter_dimensions,
                }
            }

    def _get_file_object(
        self,
    ) -> None:
        """Get an open file object"""
        if isinstance(self.file_path, (str, pl.Path)):
            self.file_object = open(self.file_path, "r")
            self.line_number = 0
            self.eof = False
        else:
            raise TypeError("file_path must be a file path")
        return

    def _get_control_variables(
        self,
    ) -> dict:
        """Get control file variables from PRMS control file

        Args:
            fpth: control file object

        Returns:
            variable_dict: dictionary with control file variables

        """
        variable_dict = {}
        while True:
            var_temp = self._get_next_variable()
            if var_temp is None:
                break
            else:
                for key, value in var_temp.items():
                    if key in (
                        "start_time",
                        "end_time",
                    ):
                        value = np.datetime64(
                            f"{value[0]:04d}-{value[1]:02d}-{value[2]:02d} "
                            + f"{value[3]:02d}:{value[4]:02d}:{value[5]:02d}"
                        )
                    elif key in ("initial_deltat",):
                        value = np.timedelta64(int(value[0]), "h")
                    variable_dict[key] = value
        self.file_object.close()
        return variable_dict

    def _get_dimensions_parameters(self):
        dimensions_dict = {}
        parameters_dict = {}
        parameter_dimensions_dict = {}
        parameters_full_dict = {}
        parameter_dimensions_full_dict = {}
        while True:
            current_position = self.file_object.tell()
            line = self._get_line()
            if line == "** Dimensions **":
                dimensions_start = current_position
            elif line == "** Parameters **":
                parameters_start = current_position
                break

        # read dimensions data
        self.file_object.seek(dimensions_start)
        while self.file_object.tell() < parameters_start:
            dim_temp = self._get_next_variable()
            if dim_temp is None:
                break
            else:
                for key, value in dim_temp.items():
                    dimensions_dict[rename_dims(key)] = value

        # read parameter data
        self.file_object.seek(parameters_start)
        self.section = PrmsFileSection.PARAMETER
        self.dimensions = dimensions_dict
        while True:
            par_temp = self._get_next_variable()
            if par_temp is None:
                break
            else:
                for key, value in par_temp.items():
                    parameters_dict[key] = value[0]
                    parameter_dimensions_dict[key] = value[1]

        # fill dictionaries that will be returned
        for key, value in dimensions_dict.items():
            parameters_full_dict[key] = value
            parameter_dimensions_full_dict[key] = None

        for key, value in parameters_dict.items():
            parameters_full_dict[key] = value

        for key, value in parameter_dimensions_dict.items():
            parameter_dimensions_full_dict[key] = value

        (
            parameters_full_dict,
            parameter_dimensions_full_dict,
        ) = expand_scalar_to_dims(
            parameters_full_dict, parameter_dimensions_full_dict
        )

        parameter_dimensions_full_dict = {
            kk: {"dims": vv}
            for kk, vv in parameter_dimensions_full_dict.items()
        }

        return parameters_full_dict, parameter_dimensions_full_dict

    def _get_parameters(self):
        parameters_dict = {}
        return parameters_dict

    def _get_line(self) -> str:
        line = self.file_object.readline()
        if line == "":
            self.eof = True
        else:
            self.line_number += 1
        return line.rstrip()

    def _get_next_variable(
        self,
    ) -> dict:
        """Get the next variable from a file

        Returns:
            var_dict: variable dict for one variable

        """
        while True:
            line = self._get_line()
            if self.eof:
                return None
            if line == "####":
                name = self._get_line().split()[0]
                if self.section == PrmsFileSection.CONTROL:
                    data = self._parse_variable()
                elif self.section == PrmsFileSection.DIMENSIONS:
                    data = self._parse_dimension()
                elif self.section == PrmsFileSection.PARAMETER:
                    data = self._parse_parameter()
                else:
                    raise NotImplementedError(
                        "reader not implemented for "
                        + f"'{self.section.value}' section type."
                    )

                return {name: data}

    def _get_variable(
        self,
        name: str,
    ) -> dict:
        """Get a variable from a file

        Args:
            name: variable name

        Returns:
            var_dict: variable dict for one variable

        """
        while True:
            line = self._get_line()
            if self.eof:
                return None
            if line == "####":
                var_name = self._get_line().split()[0]
                if var_name == name:
                    arr = self._parse_variable()
                    return {name: arr}

    def _parse_variable(
        self,
    ) -> dict:
        """Parse the data in a file to create variable data

        Returns:
            arr: numpy array that contains the data for a variable

        """
        try:
            num_values = int(self._get_line().split()[0])
            data_type = int(self._get_line().rstrip().split()[0])
            if data_type == PrmsDataType.INTEGER.value:
                arr = np.zeros(num_values, dtype=int)
                for idx in range(num_values):
                    arr[idx] = int(self._get_line().split()[0])
            elif data_type == PrmsDataType.FLOAT.value:
                arr = np.zeros(num_values, dtype=float)
                for idx in range(num_values):
                    arr[idx] = float(self._get_line().split()[0])
            elif data_type == PrmsDataType.CHARACTER.value:
                arr = np.zeros(num_values, dtype=np.chararray)
                for idx in range(num_values):
                    arr[idx] = self._get_line().split()[0]
            else:
                raise TypeError(
                    f"data type ({data_type}) can only be "
                    + f"int ({PrmsDataType.INTEGER.value}), "
                    + f"float ({PrmsDataType.FLOAT.value}), "
                    + f"or character ({PrmsDataType.CHARACTER.value}). "
                    + f"Error on line {self.line_number} in PRMS "
                    + f"input file '{self.file_object.name}'."
                )
        except TypeError:
            raise ValueError(
                f"Error on line {self.line_number} in PRMS "
                + f"input file '{self.file_object.name}'."
            )
        return arr

    def _parse_dimension(
        self,
    ) -> dict:
        """Parse the data in a file to create dimension data

        Returns:
            arr: int that contains the dimension of a variable

        """
        try:
            dimension = int(self._get_line().split()[0])
        except:
            raise ValueError(
                f"Error on line {self.line_number} in PRMS "
                + f"input file '{self.file_object.name}'."
            )
        return dimension

    def _parse_parameter(
        self,
    ) -> dict:
        """Parse the data in a file to create parameter data

        Returns:
            arr: numpy array that contains the data for a parameter

        """
        try:
            num_dims = int(self._get_line().split()[0])
            dim_names = []
            dims = []
            for idx in range(num_dims):
                dim_name = rename_dims(self._get_line().split()[0])
                dim_names.append(dim_name)
                dims.append(self.dimensions[dim_name])
            # if len(dims) == 1:
            #     dim_names = [1] + dim_names
            #     dims = [1] + dims
            shape = tuple(dims[::-1])
            dim_names = tuple(dim_names[::-1])
            len_array = int(self._get_line().split()[0])
            data_type = int(self._get_line().rstrip().split()[0])
            if data_type == PrmsDataType.INTEGER.value:
                arr = np.zeros(len_array, dtype=int)
                for idx in range(len_array):
                    arr[idx] = int(self._get_line().split()[0])
            elif data_type == PrmsDataType.FLOAT.value:
                arr = np.zeros(len_array, dtype=float)
                for idx in range(len_array):
                    arr[idx] = float(self._get_line().split()[0])
            elif data_type == PrmsDataType.CHARACTER.value:
                arr = np.zeros(len_array, dtype=np.chararray)
                for idx in range(len_array):
                    arr[idx] = self._get_line().split()[0]
            else:
                raise TypeError(
                    f"data type ({data_type}) can only be "
                    + f"int ({PrmsDataType.INTEGER.value}), "
                    + f"float ({PrmsDataType.FLOAT.value}), "
                    + f"or character ({PrmsDataType.CHARACTER.value}). "
                    + f"Error on line {self.line_number} in PRMS "
                    + f"input file '{self.file_object.name}'."
                )
        except TypeError:
            raise ValueError(
                f"Error on line {self.line_number} in PRMS "
                + f"input file '{self.file_object.name}'."
            )
        if len(shape) == 2:
            arr = arr.reshape(shape)
        return arr, dim_names
