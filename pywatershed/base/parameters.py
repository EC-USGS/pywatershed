from types import MappingProxyType
from typing import Union

import numpy as np
import xarray as xr

from .data_model import DatasetDict, dd_to_nc4_ds, dd_to_xr_ds

# MappingProxyType used as per
# https://adamj.eu/tech/2022/01/05/how-to-make-immutable-dict-in-python/


class Parameters(DatasetDict):
    """Parameter base class

    This is a subclass of data_model.DatasetDict, but all the data are
    read-only by design.
    Parameters has all the same methods as DatasetDict plus several new
    ones that map to DatasetDict as follows:

    * parameters: dd.variables
    * get_param_values: get values from dd.variables
    * get_dim_values: get values from dd.dims

    Args:
        dims: A dictionary of pairs of `dim_names: dim_len` where `dim_len` is
            an integer value.
        coords: A dictionary of pairs of `coord_names: coord_data` where
            `coord_data` is an np.ndarray.
        data_vars: A dictionary of pairs of `var_names: var_data` where
            `coord_data` is an np.ndarray.
        metadata: For all names in `coords` and `data_vars`, metadata entries
            with the required fields:

            - dims: tuple of names in dim,
            - attrs: dictionary whose values may be strings, ints, floats

            The metadata argument may also contain a special `global` key
            paired with a dictionary of global metadata of arbitrary name and
            values of string, integer, or float types.
        encoding: The encoding attributes to/from file when reading/writing.
        validate: A bool that defaults to True, enforcing the consistency
            of the supplied dictionaries
    """

    def __init__(
        self,
        dims: dict = None,
        coords: dict = None,
        data_vars: dict = None,
        metadata: dict = None,
        encoding: dict = None,
        validate: bool = True,
    ) -> None:
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

    def to_xr_ds(self) -> xr.Dataset:
        # must pass as dictionary not a mapping proxy: dict = mp | {}
        return dd_to_xr_ds(_set_dict_read_write(self.data))

    def to_nc4_ds(self, filename) -> None:
        # must pass as dictionary not a mapping proxy: dict = mp | {}
        dd_to_nc4_ds(_set_dict_read_write(self.data), filename)
        return

    def to_dd(self) -> DatasetDict:
        # must pass as dictionary not a mapping proxy: dict = mp | {}
        return DatasetDict(_set_dict_read_write(self.data))

    @classmethod
    def merge(cls, *param_list, copy=True, del_global_src=True):
        dd_list = [
            DatasetDict.from_dict(_set_dict_read_write(pp.data))
            for pp in param_list
        ]
        merged = super().merge(*dd_list, copy=copy, del_global_src=True)
        return merged


# TODO: test that these dont modify in place
def _set_dict_read_write(mp: MappingProxyType):
    dd = mp | {}
    for kk, vv in dd.items():
        if isinstance(vv, MappingProxyType):
            dd[kk] = _set_dict_read_write(dd[kk] | {})
        elif isinstance(vv, np.ndarray):
            vv.flags.writeable = True

    return dd


def _set_dict_read_only(dd: dict):
    for kk, vv in dd.items():
        if isinstance(vv, dict):
            _set_dict_read_only(vv)
            dd[kk] = MappingProxyType(vv)
        elif isinstance(vv, np.ndarray):
            vv.flags.writeable = False

    return MappingProxyType(dd)
