import pathlib as pl
from typing import Iterable, Union

import numpy as np
import yaml

from ..constants import __pynhm_root__

fileish = Union[str, pl.Path]
varoptions = Union[str, list, tuple]

"""Metadata are static, there is no point in a class"""


dims_file = __pynhm_root__ / "static/metadata/dimensions.yaml"
control_file = __pynhm_root__ / "static/metadata/control.yaml"
params_file = __pynhm_root__ / "static/metadata/parameters.yaml"
vars_file = __pynhm_root__ / "static/metadata/variables.yaml"


def load_yaml_file(the_file: pl.Path) -> dict:
    """Load a yaml file

    Args:
        the_file: metadata yaml file

    Returns:
        data: dictionary with metadata in a metadata yaml file

    """
    with pl.Path(the_file).open("r") as file_stream:
        data = yaml.load(file_stream, Loader=yaml.Loader)
    return data


dimensions = load_yaml_file(dims_file)
control = load_yaml_file(control_file)
parameters = load_yaml_file(params_file)
variables = load_yaml_file(vars_file)


def meta_netcdf_type(meta_item: dict) -> str:
    """Get the NetCDF data type for a variable

    Args:
        meta_item: variable meta data

    Returns:
        netcdf_type_str: NetCDF data type

    """
    type_str = meta_type(meta_item)
    if type_str == "int32":
        netcdf_type_str = "i4"
    elif type_str == "bool":
        netcdf_type_str = "i4"
    elif type_str == "float32":
        netcdf_type_str = "f4"
    elif type_str == "float64":
        netcdf_type_str = "f8"
    else:
        raise ValueError(f"invalid metadata type '{type_str}'")
    return netcdf_type_str


def meta_numpy_type(meta_item: dict) -> str:
    """Get the numpy dtype for a variable

    Args:
        meta_item: variable meta data

    Returns:
        numpy_type_str: numpy dtype

    """
    type_str = meta_type(meta_item)
    if type_str == "int32":
        numpy_type_str = int
    elif type_str == "bool":
        numpy_type_str = bool
    elif type_str == "float32":
        numpy_type_str = np.float32
    elif type_str == "float64":
        numpy_type_str = np.float64
    else:
        raise ValueError(f"invalid metadata type '{type_str}'")
    return numpy_type_str


def meta_type(meta_item: dict) -> str:
    """Get variable data type from meta data

    Args:
        meta_item: variable meta data

    Returns:
        type_str: variable type from meta data

    """
    return meta_item["type"]


def meta_dimensions(meta_item: dict) -> list:
    """Get the dimensions from a meta dictionary for a variable name

    Args:
        meta_item: meta dictionary for a variable name

    Returns:
        dimensions: tuple with dimension strings for a variable name

    """
    return list(meta_item["dimensions"].values())


# TODO: this class is pointless since the data is always the same. it
# could just be the module. Evenutally move that way.
def is_available(variable_name: str) -> bool:
    """Determine if a variable is available in the meta data

    Args:
        variable_name: variable name

    Returns:
        avail: boolean indicating of the variable name is available

    """
    if (
        variable_name in variables.keys()
        or variable_name in dimensions.keys()
        or variable_name in control.keys()
        or variable_name in parameters.keys()
    ):
        avail = True
    else:
        avail = False
    return avail


def _get_meta_in_list(meta_dict: dict, the_list: Iterable) -> dict:
    return {key: value for key, value in meta_dict.items() if key in the_list}


def get_dims(var_list: Iterable) -> dict:
    """Get the variable metadata for the supplied subclass."""
    return _get_meta_in_list(dimensions, var_list)


def get_control(var_list: Iterable) -> dict:
    """Get the variable metadata for the supplied subclass."""
    return _get_meta_in_list(control, var_list)


def get_params(inputs_list: Iterable) -> dict:
    """Get the variable metadata for the supplied subclass."""
    return _get_meta_in_list(parameters, inputs_list)


def get_vars(var_list: Iterable) -> dict:
    """Get the variable metadata for the supplied subclass."""
    return _get_meta_in_list(variables, var_list)


def get_dimensions(
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
    variable_dict = find_variables(variables)
    return {
        key: [dimension for tag, dimension in value["dimensions"].items()]
        for key, value in variable_dict.items()
    }


def find_variables(vars: varoptions) -> dict:
    """
    Find metadata for variables from variable, dimensions, control, and
    parameter metadata types.

    Args:
        variables: select variable names

    Returns:
        variable_dict: metadata for select variable names. An empty
        dictionary will be returned if the select variable names are not
        found.

    """
    if isinstance(vars, str):
        vars = [vars]
    variable_dict = {}
    for variable_name in vars:
        if variable_name in variables.keys():
            variable_dict[variable_name] = variables[variable_name]
        elif variable_name in dimensions.keys():
            variable_dict[variable_name] = dimensions[variable_name]
        elif variable_name in control.keys():
            variable_dict[variable_name] = control[variable_name]
        elif variable_name in parameters.keys():
            variable_dict[variable_name] = parameters[variable_name]
    return variable_dict


def get_types(variables: Iterable) -> dict:
    """Get the types for the supplied variables."""
    vars = find_variables(variables)
    return {kk: meta_type(vv) for kk, vv in vars.items()}


def get_numpy_types(variables: Iterable) -> dict:
    """Get the types for the supplied variables."""
    vars = find_variables(variables)
    return {kk: meta_numpy_type(vv) for kk, vv in vars.items()}
