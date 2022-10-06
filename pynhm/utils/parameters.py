import pathlib as pl
from typing import Union

import numpy as np

from ..constants import ft2_per_acre, inches_per_foot
from .prms5_file_util import PrmsFile

fileish = Union[str, pl.PosixPath, dict]
listish = Union[str, list, tuple]


class PRMSParameters:
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
    ) -> "PRMSParameters":

        self.parameters = parameter_dict
        if parameter_dimensions_dict is None:
            self.parameter_dimensions = self._parameter_dimensions()
        else:
            self.parameter_dimensions = parameter_dimensions_dict

        # could insert dimension data here. going to add this one for now.
        # it's a little unclear if constants like this are parameters or not
        self.parameters["doy"] = 366
        self.parameter_dimensions["doy"] = None

        return

    def _parameter_dimensions(self) -> dict:
        # build dimensions from data
        dimensions = self.dimensions
        parameter_dimensions_dict = {}
        for key, value in self.parameters.items():
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

        return parameter_dimensions_dict

    def get_parameters(self, keys: listish) -> "PRMSParameters":
        """Get a subset of keys in the parameter dictionary

        Args:
            keys: keys to retrieve from the full PRMS parameter object

        Returns:
            PRMSParameters : subset of full parameter dictionary
                Passed keys that do not exist in the full parameter
                dictionary are skipped.

        """
        if isinstance(keys, str):
            keys = [keys]

        return PRMSParameters(
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
    def dimensions(self) -> dict:
        """Get the dimensions from the parameters

        Returns:
            dimensions in the PRMS parameter dictionary

        """
        dimensions = {}
        for key, value in self.parameters.items():
            if isinstance(value, int):
                dimensions[key] = value
        return dimensions

    @property
    def nhm_hru_coordinate(self) -> np.ndarray:
        """Get the nhm hru coordinate

        Returns:
            id: nhm coordinate for each hru

        """
        if "nhm_id" in self.parameters.keys():
            nhm_id = self.parameters["nhm_id"]
        else:
            nhm_id = np.arange(1, self.parameters["nhru"] + 1)
        return nhm_id

    @property
    def nhm_segment_coordinate(self) -> np.ndarray:
        """Get the nhm segment coordinate

        Returns:
            id: nhm coordinate for each segment

        """
        if "nhm_seg" in self.parameters.keys():
            nhm_seg = self.parameters["nhm_seg"]
        else:
            nhm_seg = np.arange(1, self.parameters["nsegment"] + 1)
        return nhm_seg

    @property
    def nhm_coordinates(self) -> dict:
        """Get the nhm coordinates

        Returns:
            id: nhm coordinates

        """
        return {
            "nhm_id": self.nhm_hru_coordinate,
            "nhm_seg": self.nhm_segment_coordinate,
        }

    @staticmethod
    def load(parameter_file: fileish) -> "PRMSParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PRMSParameters: full PRMS parameter dictionary

        """
        data = PrmsFile(parameter_file, "parameter").get_data()
        return PRMSParameters(
            data["parameter"]["parameters"],
            data["parameter"]["parameter_dimensions"],
        )

    def hru_in_to_cfs(self, time_step: np.timedelta64) -> np.ndarray:
        "Derived parameter converting inches to cfs on hrus."
        time_step_seconds = time_step / np.timedelta64(1, "s")
        return self.hru_in_to_cf / time_step_seconds

    @property
    def hru_in_to_cf(self) -> np.ndarray:
        "Derived parameter converting inches to cubic feet on hrus."
        return self.parameters["hru_area"] * ft2_per_acre / (inches_per_foot)

    def subset_to_ids(
        self, nhm_id: np.ndarray, nhm_seg: np.ndarray
    ) -> "PRMSParameters":
        wh_nhm_id = np.where(np.isin(self.parameters["nhm_id"], nhm_id))
        wh_nhm_seg = np.where(np.isin(self.parameters["nhm_seg"], nhm_seg))

        n_nhm_id = len(wh_nhm_id[0])
        n_nhm_seg = len(wh_nhm_seg[0])
        assert n_nhm_id == len(nhm_id)
        assert n_nhm_seg == len(nhm_seg)

        n_nhm_id_orig = self.parameters["nhru"]
        n_nhm_seg_orig = self.parameters["nsegment"]

        parameters_sub = {}
        for pkey, pval in self.parameters.items():
            # print(f"\n{pkey}")
            if isinstance(pval, np.ndarray):
                shp = list(pval.shape)
                if n_nhm_id_orig in shp:
                    wh_shp = shp.index(n_nhm_id_orig)
                    indices = {wh_shp: wh_nhm_id}
                elif n_nhm_seg_orig in shp:
                    wh_shp = shp.index(n_nhm_seg_orig)
                    indices = {wh_shp: wh_nhm_seg}
                else:
                    print(pkey, shp)
                    parameters_sub[pkey] = pval
                    continue

                subset_inds = [
                    indices.get(dim, slice(None)) for dim in range(pval.ndim)
                ]
                parameters_sub[pkey] = np.squeeze(pval[subset_inds][0])

            else:
                print(pkey)
                parameters_sub[pkey] = pval
                continue

        # return wh_nhm_id, wh_nhm_seg?
        # override some parameters?

        return PRMSParameters(parameters_sub)
