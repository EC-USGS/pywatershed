from collections.abc import Mapping

import numpy as np


class DotDict(dict):
    """
    Quick and dirty implementation of a dot-able dict, which allows access and
    assignment via object properties rather than dict indexing.
    source: https://stackoverflow.com/questions/13520421/recursive-dotdict
    modified OrderdDict to standard dict
    """

    def __init__(self, *args, **kwargs):
        od = dict(*args, **kwargs)
        for key, val in od.items():
            if isinstance(val, Mapping):
                value = DotDict(val)
            else:
                value = val
            self[key] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError as ex:
            raise AttributeError(f"No attribute called: {name}") from ex

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError as ex:
            raise AttributeError(f"No attribute called: {k}") from ex

    __setattr__ = dict.__setitem__


class PrmsParameters(dict):
    """
    PRMS parameter class

    Parameters
    ----------
    dict : dict
        parameters dictionary

    """

    def __init__(self, parameter_dict):
        self.parameters = DotDict(parameter_dict)

    def get_parameters(self, keys):
        """
        Get a subset of keys in the parameter dictionary

        Parameters
        ----------
        keys : str list or tuple
            parameters to extract from the full parameter dictionary

        Returns
        -------
        subset : dict
            Subset of full parameter dictionary with the passed parameter
            keys. Passed keys that do not exist in the full parameter
            dictionary are set to None

        """
        if isinstance(keys, str):
            keys = [keys]

        return PrmsParameters(
            {
                key: self.parameters.get(key)
                for key in keys
                if key in self.parameters.keys()
            }
        )

    @staticmethod
    def load(parameter_file):
        """
        Load parameters from a PRMS parameter file

        Parameters
        ----------
        parameter_file : str or pathlib.Path


        Returns
        -------
        Parameters : PrmsParameters
            PRMS parameter object


        """
        (
            dimensions,
            parameter_data,
            parameter_dimensions,
            parameter_types,
        ) = _load_prms_parameters(parameter_file)
        parameters = parameter_data.copy()
        for key, value in dimensions.items():
            parameters[key] = value
        return PrmsParameters(parameters)


def _load_prms_parameters(parameter_file):
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
