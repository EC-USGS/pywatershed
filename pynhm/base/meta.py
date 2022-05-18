import pathlib as pl
from typing import Iterable, Union

import numpy as np
import yaml

from ..constants import __pynhm_root__

fileish = Union[str, pl.Path]
varoptions = Union[str, list, tuple]


def load_yaml_file(the_file):
    with pl.Path(the_file).open("r") as file_stream:
        data = yaml.load(file_stream, Loader=yaml.Loader)
    return data


def meta_netcdf_type(meta_item: dict) -> str:
    type_str = meta_type(meta_item)
    if type_str == "I":
        netcdf_type_str = "i4"
    elif type_str == "B":
        netcdf_type_str = "i4"
    elif type_str == "F":
        netcdf_type_str = "f4"
    elif type_str == "D":
        netcdf_type_str = "f8"
    else:
        raise ValueError(f"invalid metadata type '{type_str}'")
    return netcdf_type_str


def meta_numpy_type(meta_item: dict) -> str:
    type_str = meta_type(meta_item)
    if type_str == "I":
        numpy_type_str = int
    elif type_str == "B":
        numpy_type_str = int
    elif type_str == "F":
        numpy_type_str = np.float32
    elif type_str == "D":
        numpy_type_str = float
    else:
        raise ValueError(f"invalid metadata type '{type_str}'")
    return numpy_type_str


def meta_type(meta_item: dict) -> str:
    return meta_item["type"]


def meta_dimensions(meta_item: dict) -> list:
    """Get the dimensions from a meta dictionary for a variable name

    Args:
        meta_item: meta dictionary for a variable name

    Returns:
        dimensions: tuple with dimension strings for a variable name

    """
    return [value for key, value in meta_item["dimensions"].items()]


class Meta:
    def __init__(
        self,
        dimensions: fileish = (
            __pynhm_root__ / "static/metadata/dimensions.yaml"
        ),
        control: fileish = (__pynhm_root__ / "static/metadata/control.yaml"),
        parameters: fileish = (
            __pynhm_root__ / "static/metadata/parameters.yaml"
        ),
        variables: fileish = (
            __pynhm_root__ / "static/metadata/variables.yaml"
        ),
    ):

        self.name = "Meta"

        print(dimensions)

        self.dimensions = load_yaml_file(dimensions)
        self.control = load_yaml_file(control)
        self.parameters = load_yaml_file(parameters)
        self.variables = load_yaml_file(variables)

        return

    def is_available(self, variable_name: str) -> bool:
        if (
            variable_name in self.variables.keys()
            or variable_name in self.dimensions.keys()
            or variable_name in self.control.keys()
            or variable_name in self.parameters.keys()
        ):
            avail = True
        else:
            avail = False
        return avail

    def _get_meta_in_list(self, meta: str, the_list: Iterable) -> dict:
        return {
            key: value
            for key, value in getattr(self, meta).items()
            if key in the_list
        }

    def get_dims(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("dimensions", var_list)

    def get_control(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("control", var_list)

    def get_params(self, inputs_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("parameters", inputs_list)

    def get_vars(self, var_list: Iterable) -> dict:
        """Get the variable metadata for the supplied subclass."""
        return self._get_meta_in_list("variables", var_list)

    def get_dimensions(
        self,
        variables: varoptions,
    ) -> dict:
        """
        Get the dimensions strings for one or more variables for a specific
        metadata data type.

        Args:
            variables: variable name(s)

        Returns:
            dimensions: dictionary with dimension strings for each variable

        """
        if isinstance(variables, str):
            variables = [variables]
        variable_dict = self.find_variables(variables)
        return {
            key: [dimension for tag, dimension in value["dimensions"].items()]
            for key, value in variable_dict.items()
        }

    def find_variables(self, variables: varoptions) -> dict:
        """

        Args:
            variables:

        Returns:

        """
        if isinstance(variables, str):
            variables = [variables]
        variable_dict = {}
        for variable_name in variables:
            if variable_name in self.variables.keys():
                variable_dict[variable_name] = self.variables[variable_name]
            elif variable_name in self.dimensions.keys():
                variable_dict[variable_name] = self.dimensions[variable_name]
            elif variable_name in self.control.keys():
                variable_dict[variable_name] = self.control[variable_name]
            elif variable_name in self.parameters.keys():
                variable_dict[variable_name] = self.parameters[variable_name]
        return variable_dict
