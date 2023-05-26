from types import MappingProxyType
from typing import Union

import numpy as np

from .data_model import DatasetDict


class Parameters(DatasetDict):
    def __init__(
        self,
        dims: dict = {},
        coords: dict = {},
        data_vars: dict = {},
        metadata: dict = {},
        encoding: dict = {},
        validate: bool = True,
    ) -> "Parameters":
        """Parameter class

        This is a subclass of data_model.DatasetDict, but that all the data are
        read-only by design.
        Parameters has all the same methods as DatasetDict plus several new
        ones that map to DatasetDict as follows:
            parameters: dd.variables
            get_param_values: get values from dd.variables
            get_dim_values: get values from dd.dims

        """
        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
            metadata=metadata,
            encoding=encoding,
            validate=validate,
        )

        for kk in self._keys():
            self[f"_{kk}"] = _set_dict_read_only(self[kk])

        return

    @property
    def parameters(self) -> dict:
        return self.variables

    def get_param_values(
        self,
        keys: Union[list, str] = None,
    ) -> Union[dict, np.ndarray]:
        """Get the values of the parameters (coords or data_vars) by keys

        Also see:
            subset() method is a Parameter object is desired.
        """
        if not isinstance(keys, list):
            return self.variables[keys]
        else:
            return {kk: self.variables[kk] for kk in keys}

    def get_dim_values(
        self,
        keys: Union[list, tuple, str] = None,
    ) -> Union[dict, np.ndarray]:
        """Get the values of the dimensions by keys."""
        if not isinstance(keys, (list, tuple)):
            return self.dims[keys]
        else:
            return {kk: self.dims[kk] for kk in keys}


def _set_dict_read_only(dd: dict):
    for kk, vv in dd.items():
        if isinstance(vv, dict):
            _set_dict_read_only(vv)
            dd[kk] = MappingProxyType(vv)
        elif isinstance(vv, np.ndarray):
            vv.flags.writeable = False

    return MappingProxyType(dd)
