import io
import pathlib as pl
from enum import Enum
from typing import Tuple, Union

import numpy as np

fileish = Union[str, pl.PosixPath, io.IOBase]


class PrmsDataType(Enum):
    INTEGER = 1
    FLOAT = 2
    CHARACTER = 4


def _get_file_object(
    fobj: fileish,
    rewind: bool = False,
) -> io.IOBase:
    """Get an open file object

    Args:
        fobj: file path or file object
        rewind: boolean indicating if the file object should
            be rewound.

    Returns:
        fobj: open file object

    """
    if isinstance(fobj, io.IOBase):
        if fobj.closed:
            open(fobj.name).read()
        if rewind:
            fobj.seek(0)
    elif isinstance(fobj, (str, pl.Path)):
        fobj = open(fobj, "r")
    else:
        raise TypeError("fobj must be a file object or a file path")
    return fobj


def _get_control_variables(
    fpth: fileish,
) -> dict:
    """Get control file variables from PRMS control file

    Args:
        fpth: control file path or file object

    Returns:
        variable_dict: dictionary with control file variables

    """
    variable_dict = {}
    fobj = _get_file_object(fpth)
    while True:
        fobj, var_temp = _get_next_variable(fobj)
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
    fobj.close()
    return variable_dict


def _get_next_variable(
    fobj: fileish,
    rewind: bool = False,
) -> Tuple[io.IOBase, dict]:
    """Get the next variable from a file

    Args:
        fobj: file path or file object
        rewind: boolean indicating if the file object should
            be rewound.

    Returns:
        fobj: control file object
        var_dict: variable dict for one variable

    """
    fobj = _get_file_object(fobj, rewind=rewind)
    while True:
        line = fobj.readline()
        if line == "":
            return fobj, None
        if line.rstrip() == "####":
            name = fobj.readline().rstrip().split()[0]
            fobj, arr = _parse_variable(fobj)
            return fobj, {name: arr}


def _get_variable(
    name: str,
    fobj: fileish,
    rewind: bool = False,
) -> Tuple[io.IOBase, dict]:
    """Get a variable from a file

    Args:
        name: variable name
        fobj: file path or file object
        rewind: boolean indicating if the file object should
            be rewound.

    Returns:
        fobj: control file object
        var_dict: variable dict for one variable

    """
    fobj = _get_file_object(fobj, rewind=rewind)
    while True:
        line = fobj.readline()
        if line == "":
            return fobj, None
        if line.rstrip() == "####":
            var_name = fobj.readline().rstrip().split()[0]
            if var_name == name:
                fobj, arr = _parse_variable(fobj)
                return fobj, {name: arr}


def _parse_variable(
    fobj: io.IOBase,
) -> Tuple[io.IOBase, dict]:
    """Parse the data in a file to create variable data

    Args:
        fobj: file object

    Returns:
        fobj: control file object
        arr: numpy array that contains the data for a variable

    """
    try:
        num_values = int(fobj.readline().rstrip().split()[0])
        data_type = int(fobj.readline().rstrip().split()[0])
        if data_type == PrmsDataType.INTEGER.value:
            arr = np.zeros(num_values, dtype=int)
            for idx in range(num_values):
                arr[idx] = int(fobj.readline().rstrip().split()[0])
        elif data_type == PrmsDataType.FLOAT.value:
            arr = np.zeros(num_values, dtype=float)
            for idx in range(num_values):
                arr[idx] = float(fobj.readline().rstrip().split()[0])
        elif data_type == PrmsDataType.CHARACTER.value:
            arr = np.zeros(num_values, dtype=np.chararray)
            for idx in range(num_values):
                arr[idx] = fobj.readline().rstrip().split()[0]
        else:
            raise TypeError(
                f"data type ({data_type}) can only be "
                + f"int ({PrmsDataType.INTEGER.value}), "
                + f"float ({PrmsDataType.FLOAT.value}), "
                + f"or character ({PrmsDataType.CHARACTER.value})"
            )
    except:
        raise ValueError("oops...add better information")
    return fobj, arr
