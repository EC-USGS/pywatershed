"""The metadata module.

Metadata for variables, parameters and dimensions are accssed by the functions
in this module.

The metadata are static (pywatershed/static/metadata) so this is a module and
not a class.

"""

import pathlib as pl
from typing import Iterable, Union

import numpy as np
import yaml

from ..constants import __pywatershed_root__

varoptions = Union[str, list, tuple]

dims_file = __pywatershed_root__ / "static/metadata/dimensions.yaml"
control_file = __pywatershed_root__ / "static/metadata/control.yaml"
params_file = __pywatershed_root__ / "static/metadata/parameters.yaml"
vars_file = __pywatershed_root__ / "static/metadata/variables.yaml"


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


# change the representaiton of dimensions in the metadata to tuples
def _dims_to_tuples(dd):
    result = {}
    for kk, vv in dd.items():
        if kk == "dims":
            result[kk] = tuple(vv.values())
        elif isinstance(vv, dict):
            result[kk] = _dims_to_tuples(vv)
        else:
            result[kk] = vv

    return result


dimensions = _dims_to_tuples(load_yaml_file(dims_file))
control = _dims_to_tuples(load_yaml_file(control_file))
parameters = _dims_to_tuples(load_yaml_file(params_file))
variables = _dims_to_tuples(load_yaml_file(vars_file))


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
    return list(meta_item["dims"])


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


# repeat this for dims, control, params as needed. or refactor to DRY.
def filter_vars(var_list: Iterable, key, vals=None) -> dict:
    return _filter_dict(get_vars(var_list), key, vals)


def _filter_dict(meta_dict, key, vals=None) -> dict:
    if not isinstance(vals, list):
        vals = [vals]
    result = {}
    for name, meta in meta_dict.items():
        if key in meta:
            if not vals:
                result[name] = meta
            elif meta[key] in vals:
                result[name] = meta
    return result


def get_dimensions(
    vars: varoptions,
) -> dict:
    """
    Get the dimensions strings for one or more variables for a specific
    metadata data type.

    Args:
        vars: variable name(s)

    Returns:
        dimensions: dictionary with dimension strings for each variable

    """
    if isinstance(vars, str):
        vars = [vars]
    variable_dict = find_variables(vars)
    return {key: value["dims"] for key, value in variable_dict.items()}


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


def get_units(vars: varoptions, to_pint: bool = False) -> dict:
    """For names, return units

    Args:
        vars: select variable names in variables, dimensions, control and/or
            parameters

    Returns:
        units_dict: units for select variable names. An empty
        dictionary will be returned if the select variable names are not
        found.
    """
    units_dict = {
        key: val["units"] for key, val in find_variables(vars).items()
    }
    if to_pint:
        units_dict = prms_to_pint(units_dict)
    return units_dict


def get_types(variables: Iterable) -> dict:
    """Get the types for the supplied variables."""
    vars = find_variables(variables)
    return {kk: meta_type(vv) for kk, vv in vars.items()}


def get_numpy_types(variables: Iterable) -> dict:
    """Get the types for the supplied variables."""
    vars = find_variables(variables)
    return {kk: meta_numpy_type(vv) for kk, vv in vars.items()}


def prms_to_pint(var_units_dict: dict) -> dict:
    """Convert PRMS units to pint units

    A work in progress to keep track of PRMS units/dimensions that dont
    work with pint and their.

    Args: A dictionary of {var_name: prms_units, ...}

    Returns: A dictonary of {var_name: pint_units, ...}
    """
    prms_to_pint = {
        "decimal fraction": "dimensionless",
        "cfs": "feet **3 / seconds",
        "cubicfeet": "feet ** 3",
        # more to come!
    }

    return {
        key: (val if val not in prms_to_pint.keys() else prms_to_pint[val])
        for key, val in var_units_dict.items()
    }
