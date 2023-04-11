import json
from copy import deepcopy

import netCDF4 as nc4
import numpy as np
import xarray as xr

from ..base import meta
from ..constants import fileish, ft2_per_acre, inches_per_foot, listish, ndoy
from .prms5_file_util import PrmsFile
from ..utils import data_model as dm


class JSONParameterEncoder(json.JSONEncoder):
    """
    Simple encoder to cast numpy objects to json-friendly formats
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(JSONParameterEncoder, self).default(obj)


def _json_load(json_filename):
    pars = json.load(open(json_filename))
    # need to convert lists to numpy arrays
    for k, v in pars.items():
        if isinstance(v, list):
            pars[k] = np.array(v)
    return pars


class Parameters:
    def __init__(
        self,
        dims: dict,
        coords: dict,
        data_vars: dict,
        attrs: dict,
        metadata: dict,
        encoding: dict = None,
    ) -> "Parameters":
        self._data_vars = data_vars
        self._dims = dims
        self._coords = coords
        self._attrs = attrs

    @property
    def dims(self) -> dict:
        """Return the dimensions over all parameters"""
        return self._dims

    @property
    def coords(self) -> dict:
        """Return the coordinates over all parameters"""
        return self._coords

    @property
    def parameters(self) -> dict:
        """Return all parameters"""
        return self._data_vars

    @property
    def attrs(self) -> dict:
        """Return the parameter (global) attributes"""
        return self._data_vars

    @property
    def spatial_coord_names(self) -> dict:
        """Return the spatial coordinate names."""
        return {kk: vv for kk, vv in self._attrs.items() if "spatial" in kk}

    @staticmethod
    def from_netcdf(nc_file_list) -> "Parameters":
        """Load parameters object from a netcdf file(s?)"""
        # Provide in base class
        # handle more than one file? # see prms ?
        raise NotImplementedError

    def write_netcdf(self, filename) -> None:
        """Write parameters to a json file"""
        xr.Dataset.from_dict(
            {
                "dims": self._dims,
                "coords": self._coords,
                "data_vars": self._data_vars,
                "attrs": self._attrs,
            }
        ).to_netcdf(filename)
        return

    @staticmethod
    def from_json(json_filename: fileish) -> "Parameters":
        """Load parameters object from a json file"""
        raise NotImplementedError

    def write_json(self, json_filename) -> None:
        """Write parameters to a json file"""
        raise NotImplementedError

    def subset(
        self,
        keys: listish = None,
        process: str = None,
        keep_attrs=True,
    ) -> "Parameters":
        """Returns a Parameters object with a subset of the data."""
        # Provide in base class
        raise NotImplementedError

    def get_parameters(
        self, keys: listish, process: str = None, dims: bool = False
    ) -> dict:
        """Returns Parameter data for requested keys."""
        # Provide in base class
        raise NotImplementedError


class PrmsParameters(Parameters):
    """
    PRMS parameter class

    Parameters
    ----------
    parameter_dict : dict
        parameters dictionary: either structure
          * param: value
          * process: {param: value ... }
        where the later is a parameter dictionary grouped by process.
        The keys for process should be either the class itself, class.name, or
        type(class.__name__).
    parameter_dimensions_dict : dict
        parameters dimensions dictionary with a structure mirring the parameter
        dict as described above but with shape tuples in place of parameter
        value data.

    Returns:
        PrmsParameter object

    """

    def __init__(
        self,
        parameter_dict: dict,
        parameter_dimensions_dict: dict = None,
    ) -> "PrmsParameters":
        super().__init__(
            dims={},
            coords=parameter_dimensions_dict,
            data_vars=parameter_dict,
            attrs={},
            metadata={},
        )

        self._params_sep_procs = all(
            [isinstance(pp, dict) for pp in self.parameters.values()]
        )

        # build dimensions from data
        if parameter_dimensions_dict is None:
            # todo: handle self._params_sep_procs
            parameter_dimensions_dict = {}
            for key, value in parameter_dict.items():
                param_dim_names = meta.get_params(key)
                if len(param_dim_names):
                    param_dim_names = list(
                        param_dim_names[key]["dimensions"].values()
                    )
                else:
                    param_dim_names = [key]

                common_params = set(param_dim_names) & set(
                    self.dimensions.keys()
                )
                if not len(common_params):
                    parameter_dimensions_dict[key] = tuple(["unknown"])
                    continue

                param_dims = {kk: self.dimensions[kk] for kk in common_params}

                if isinstance(value, int):
                    parameter_dimensions_dict[key] = None
                elif isinstance(value, np.ndarray):
                    shape = value.shape
                    temp_dims = []
                    for isize in shape:
                        found_dim = False
                        for dim_key, dim_value in param_dims.items():
                            if dim_value == isize:
                                found_dim = True
                                temp_dims.append(dim_key)
                                break

                        if isize == 1 and not found_dim:
                            found_dim = True
                            temp_dims.append("scalar")

                    parameter_dimensions_dict[key] = tuple(temp_dims)

        self.parameter_dimensions = parameter_dimensions_dict

        return

    def subset(
        self, keys: listish = None, process: str = None
    ) -> "PrmsParameters":
        """Get a parameter object subset to passed keys and processes"""

        if (process is not None) and (keys is None):
            if not self._params_sep_procs:
                return PrmsParameters(
                    parameter_dict=deepcopy(self.parameters),
                    parameter_dimensions_dict=deepcopy(
                        self.parameter_dimensions
                    ),
                )
            return PrmsParameters(
                parameter_dict=deepcopy(self.parameters[process]),
                parameter_dimensions_dict=deepcopy(
                    self.parameter_dimensions[process]
                ),
            )

        elif (process is None) and (keys is not None):
            if not self._params_sep_procs:
                return PrmsParameters(
                    parameter_dict={
                        kk: self.parameters[kk]
                        for kk in keys
                        if kk in self.parameters.keys()
                    },
                    parameter_dimensions_dict={
                        kk: self.parameter_dimensions[kk]
                        for kk in keys
                        if kk in self.parameter_dimensions.keys()
                    },
                )
            else:
                raise NotImplementedError
            pass

        else:
            # This is still a work in progress
            raise NotImplementedError

    def get_parameters(
        self, keys: listish, process: str = None, dims: bool = False
    ) -> dict:
        """Get a subset of keys in the parameter dictionary

        Args:
            keys: keys to retrieve from the full PRMS parameter object
            process: return a single process (all keys)
            dims: return dims instead of values (default is False)

        Returns:
            dict: containing requested key:val pairs

        """
        if isinstance(keys, str):
            keys = [keys]

        if self._params_sep_procs:
            # If a parameter dict, get the values from any dict entry
            # there should never be a param with different values
            return_params = {}
            for proc_key in self.parameters.keys():
                if (process is not None) and (proc_key != type(process)):
                    continue
                # only allow a single process and do not return a process key
                var_keys = self.parameters[proc_key].keys()
                for key in keys:
                    if key in var_keys:
                        if dims:
                            return_params[key] = self.parameters_dimensions[
                                proc_key
                            ][key]
                        else:
                            return_params[key] = self.parameters[proc_key][key]

            return return_params

        else:
            if dims:
                return {
                    key: self.parameter_dimensions.get(key)
                    for key in keys
                    if key in self.parameter_dimensions.keys()
                }
            else:
                return {
                    key: self.parameters.get(key)
                    for key in keys
                    if key in self.parameters.keys()
                }

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

    def parameters_to_json(self, json_filename):
        """write the parameters dictionary out to a json file"""
        json.dump(
            self.parameters,
            open(json_filename, "w"),
            indent=4,
            cls=JSONParameterEncoder,
        )

        return

    @staticmethod
    def load_from_json(json_filename: fileish) -> "PrmsParameters":
        """Load parameters from a json file

        Args:
            : json file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        pars = _json_load(json_filename)
        params = PrmsParameters(pars)

        # add implied dimensions
        params.parameters["ndoy"] = ndoy
        params.parameters["nmonth"] = 12
        params.parameter_dimensions["ndoy"] = None
        params.parameter_dimensions["nmonth"] = None

        return params

    @staticmethod
    def load(parameter_file: fileish) -> "PrmsParameters":
        """Load parameters from a PRMS parameter file

        Args:
            parameter_file: parameter file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        data = PrmsFile(parameter_file, "parameter").get_data()
        params = PrmsParameters(
            data["parameter"]["parameters"],
            data["parameter"]["parameter_dimensions"],
        )

        # add implied dimensions
        params.parameters["ndoy"] = ndoy
        params.parameters["nmonth"] = 12
        params.parameter_dimensions["ndoy"] = None
        params.parameter_dimensions["nmonth"] = None

        return params

    @staticmethod
    def _from_nc_file(parameter_nc_file) -> dict:
        import importlib

        ds = nc4.Dataset(parameter_nc_file)
        attrs = ds.__dict__
        process_name = attrs["nhm_process"]
        process = importlib.import_module("pynhm").__dict__[process_name]
        proc_params = process.get_parameters()

        param_dict = {kk: vv[:].data for kk, vv in ds.variables.items()}
        param_dict_dimensions = {
            kk: vv.dimensions for kk, vv in ds.variables.items()
        }

        dimension_params = {
            kk: len(vv) for kk, vv in ds.dimensions.items() if kk != "scalar"
        }
        dimension_params_dims = {kk: None for kk in dimension_params.keys()}

        param_dict = {**param_dict, **dimension_params}
        param_dict_dimensions = {
            **param_dict_dimensions,
            **dimension_params_dims,
        }

        # Additional scalar parameters are in the golbal attributes
        neglected_params = list(set(proc_params).difference(param_dict.keys()))
        for nn in neglected_params:
            param_dict[nn] = ds.__dict__[nn]
            param_dict_dimensions[nn] = None

        return {"parameters": param_dict, "dimensions": param_dict_dimensions}

    @staticmethod
    def from_nc_files(proc_param_nc_file_dict: dict) -> "PrmsParameters":
        """Load parameters from a PRMS parameter file

        Args:
            proc_param_nc_file_dict: dict of proc: param_nc_file pairs
            for all processes in the model

        Returns:
            param_dict: A dictionary containing PRMSParameters objects for
            each process.

        """
        param_dim = {
            proc: PrmsParameters._from_nc_file(file)
            for proc, file in proc_param_nc_file_dict.items()
        }
        parameters = {key: val["parameters"] for key, val in param_dim.items()}
        dims = {key: val["dimensions"] for key, val in param_dim.items()}
        return PrmsParameters(parameters, dims)

    def hru_in_to_cfs(self, time_step: np.timedelta64) -> np.ndarray:
        "Derived parameter converting inches to cfs on hrus."
        time_step_seconds = time_step / np.timedelta64(1, "s")
        return self.hru_in_to_cf / time_step_seconds

    @property
    def hru_in_to_cf(self) -> np.ndarray:
        "Derived parameter converting inches to cubic feet on hrus."
        return (
            self.get_parameters("hru_area")["hru_area"]
            * ft2_per_acre
            / (inches_per_foot)
        )


class StarfitParameters(Parameters):
    """
    Starfit parameter class

    Parameters
    ----------
    parameter_dict : dict
        parameters dictionary: either structure
          * param: value
          * process: {param: value ... }
        where the later is a parameter dictionary grouped by process.
        The keys for process should be either the class itself, class.name, or
        type(class.__name__).
    parameter_dimensions_dict : dict
        parameters dimensions dictionary with a structure mirring the parameter
        dict as described above but with shape tuples in place of parameter
        value data.

    Returns:
        StarfitParameters object

    """

    def __init__(
        self,
        dims: dict,
        coords: dict,
        data_vars: dict,
        attrs: dict,
        metadata: dict,
    ) -> "StarfitParameters":
        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
            attrs=attrs,
            metadata=metadata,
        )
        # remove this throughout, no prms specific parameter methods should
        # be used in netcdf utils
        self.nhm_coordinates = {"grand_id": self.parameters["grand_id"]}

    @staticmethod
    def from_netcdf(
        resops_domain: fileish,
        istarf_conus: fileish,
        grand_dams: fileish,
        grand_ids: fileish = None,
        param_names: list = None,
    ) -> dict:
        # istarf_conus_ds = dm.nc4_ds_to_dd(istarf_conus)

        istarf_conus_ds = nc4.Dataset(istarf_conus)
        resops_ds = nc4.Dataset(resops_domain)
        grand_dams_ds = nc4.Dataset(grand_dams)

        if grand_ids is None:
            domain_grand_ids = resops_ds["grand_id"][:].data
        else:
            domain_grand_ids = nc4.Dataset(grand_ids)["grand_id"]

        # As a convenience, subset the full reservoir datasets to
        # the domain instead of requiring a preprocess.
        param_dict = {}
        param_dim_dict = {}

        # implied parameters
        param_dict["nreservoirs"] = len(domain_grand_ids)
        param_dict["nhru"] = len(domain_grand_ids)  # to remove

        def get_params(data, grand_id_name):
            wh_domain = np.where(
                np.isin(data[grand_id_name], domain_grand_ids)
            )
            assert len(wh_domain[0]) == len(domain_grand_ids)

            for vv in data.variables:
                if param_names is not None and vv not in param_names:
                    if vv.lower() != "grand_id":
                        continue

                if vv in param_dict.keys():
                    raise KeyError(
                        f"the key {vv} is already in the param_dict"
                    )

                # param_dim_dict[vv] = list(data[vv].dimensions)
                param_dict[vv] = data[vv][:][wh_domain]

                if hasattr(data[vv], "units") and "since" in data[vv].units:
                    param_dict[vv] = (
                        nc4.num2date(
                            param_dict[vv],
                            units=data[vv].units,
                            calendar=data[vv].calendar,
                            only_use_cftime_datetimes=False,
                        )
                        .filled()
                        .astype("datetime64[s]")
                    )

                if isinstance(param_dict[vv], np.ma.core.MaskedArray):
                    param_dict[vv] = param_dict[vv].data

            return

        # subset the GRanD data provided
        # https://sedac.ciesin.columbia.edu/data/set/grand-v1-dams-rev01
        get_params(grand_dams_ds, "GRAND_ID")

        # subset the istarf data provided
        # https://zenodo.org/record/4602277#.ZCtYj-zMJqs
        get_params(istarf_conus_ds, "GRanD_ID")

        # subset the resops data provided
        # https://zenodo.org/record/5893641#.ZCtakuzMJqs
        get_params(resops_ds, "grand_id")

        # make sure the three subsets are colocated
        # remove the duplicated grand_id (ridiculous)
        assert (param_dict["GRAND_ID"] == param_dict["GRanD_ID"]).all()
        assert (param_dict["GRAND_ID"] == param_dict["grand_id"]).all()
        del param_dict["GRAND_ID"], param_dict["GRanD_ID"]

        return StarfitParameters(
            attrs={}, dims={}, coords={}, data_vars=param_dict, metadata={}
        )

    def subset(self, **kwargs):
        return self

    def get_parameters(self, keys, **kwargs):
        if not isinstance(keys, list):
            keys = [keys]
        return {
            key: self.parameters[key]
            for key in keys
            if key in self.parameters.keys()
        }

    def write(self, file):
        # it's desirable to have a netcdf write method for these parameters
        pass
