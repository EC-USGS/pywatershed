import io
import pathlib
from enum import Enum

import numpy as np


class prms_type(Enum):
    INTEGER = 1
    FLOAT = 2
    CHARACTER = 4


def get_file_object(fobj, rewind=False):
    if isinstance(fobj, io.IOBase):
        if fobj.closed:
            open(fobj.name).read()
        if rewind:
            fobj.seek(0)
    elif isinstance(fobj, (str, pathlib.Path)):
        fobj = open(fobj, "r")
    else:
        raise TypeError("fobj must be a file object or a file path")
    return fobj


def get_all_variables(fpth):
    variable_dict = {}
    fobj = get_file_object(fpth)
    while True:
        fobj, var_temp = get_next_variable(fobj)
        if var_temp is None:
            break
        else:
            for key, value in var_temp.items():
                variable_dict[key] = value
    fobj.close()
    return variable_dict


def get_next_variable(fobj, rewind=False):
    fobj = get_file_object(fobj, rewind=rewind)
    name = None
    arr = None
    while True:
        line = fobj.readline()
        if line == "":
            return fobj, None
        if line.rstrip() == "####":
            name = fobj.readline().rstrip().split()[0]
            fobj, arr = _parse_variable(fobj)
            return fobj, {name: arr}


def get_variable(name, fobj, rewind=False):
    fobj = get_file_object(fobj, rewind=rewind)
    arr = None
    while True:
        line = fobj.readline()
        if line == "":
            return fobj, None
        if line.rstrip() == "####":
            var_name = fobj.readline().rstrip().split()[0]
            if var_name == name:
                fobj, arr = _parse_variable(fobj)
                return fobj, {name: arr}


def _parse_variable(fobj):
    arr = None
    try:
        num_values = int(fobj.readline().rstrip().split()[0])
        data_type = int(fobj.readline().rstrip().split()[0])
        if data_type == prms_type.INTEGER.value:
            arr = np.zeros(num_values, dtype=int)
            for idx in range(num_values):
                arr[idx] = int(fobj.readline().rstrip().split()[0])
        elif data_type == prms_type.FLOAT.value:
            arr = np.zeros(num_values, dtype=float)
            for idx in range(num_values):
                arr[idx] = float(fobj.readline().rstrip().split()[0])
        elif data_type == prms_type.CHARACTER.value:
            arr = np.zeros(num_values, dtype=np.chararray)
            for idx in range(num_values):
                arr[idx] = fobj.readline().rstrip().split()[0]
        else:
            raise TypeError(
                f"data type ({data_type}) can only be "
                + f"int ({prms_type.INTEGER.value}), "
                + f"float ({prms_type.FLOAT.value}), "
                + f"or character ({prms_type.CHARACTER.value})"
            )
    except:
        raise ValueError("oops...add better information")
    return fobj, arr
