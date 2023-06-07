from copy import deepcopy
import json

import numpy as np

from ..base import meta
from ..base.parameters import Parameters
from ..constants import fileish, ft2_per_acre, inches_per_foot, listish, ndoy
from ..utils.prms5_file_util import PrmsFile

prms_dim_names = (
    "nhru",
    "nsegment",
    "nssr",
    "ngw",
    "npoigages",
    "nobs",
    "ndeplval",
    "ndepl",
    "nmonth",
    "ndoy",
    "scalar",
)


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
        dims: dict = None,
        coords: dict = None,
        data_vars: dict = None,
        metadata: dict = None,
        encoding: dict = None,
        validate: bool = True,
    ) -> "PrmsParameters":
        if dims is None:
            dims = {}
        if coords is None:
            coords = {}
        if data_vars is None:
            data_vars = {}
        if metadata is None:
            metadata = {}
        if encoding is None:
            encoding = {}

        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
            metadata=metadata,
            encoding=encoding,
            validate=validate,
        )

        return

    def get_parameters(self, keys: listish, dims: bool = False) -> dict:
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
        for key, value in self.dims.items():
            if isinstance(value, int):
                dimensions[key] = value
        return dimensions

    @property
    def nhm_hru_coordinate(self) -> np.ndarray:
        """Get the nhm hru coordinate

        Returns:
            id: nhm coordinate for each hru

        """
        if "nhru" not in self.dims.keys():
            return None

        if "nhm_id" in self.parameters.keys():
            nhm_id = self.parameters["nhm_id"]
        else:
            nhm_id = np.arange(1, self.dims["nhru"] + 1)
        return nhm_id

    @property
    def nhm_segment_coordinate(self) -> np.ndarray:
        """Get the nhm segment coordinate

        Returns:
            id: nhm coordinate for each segment

        """
        if "nsegment" not in self.dims.keys():
            return None

        if "nhm_seg" in self.parameters.keys():
            nhm_seg = self.parameters["nhm_seg"]
        else:
            nhm_seg = np.arange(1, self.dims["nsegment"] + 1)
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
            {**self.dims, **self.parameters},
            open(json_filename, "w"),
            indent=4,
            cls=JSONParameterEncoder,
        )

        return

    @staticmethod
    def _from_dict(param_dict: dict) -> "PrmsParameters":
        """Load parameters from a dictionary of just parameters"""
        return PrmsParameters._after_load(param_dict)

    @staticmethod
    def load_from_json(json_filename: fileish) -> "PrmsParameters":
        """Load parameters from a json file

        Args:
            : json file path

        Returns:
            PrmsParameters: full PRMS parameter dictionary

        """
        pars = _json_load(json_filename)
        params = PrmsParameters._process_file_input(pars)

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
        params = PrmsParameters._process_file_input(
            data["parameter"]["parameters"],
            # data["parameter"]["parameter_dimensions"],
        )

        return params

    def _process_file_input(
        parameter_dict: dict,
        parameter_dimensions_dict: dict = None,
    ) -> "PrmsParameters":
        # move dims from params to dims
        if parameter_dimensions_dict is None:
            parameter_dimensions_dict = {}

        dims = {}
        for dd in prms_dim_names:
            if dd in parameter_dict.keys():
                dims[dd] = parameter_dict.pop(dd)
                if dd in parameter_dimensions_dict.keys():
                    _ = parameter_dimensions_dict.pop(dd)

        # add implied dimensions
        dims["ndoy"] = ndoy
        dims["nmonth"] = 12

        # build dimension metadata from data
        if len(parameter_dimensions_dict) == 0:
            for key, value in parameter_dict.items():
                param_dim_names = meta.get_params(key)[key]["dims"]
                parameter_dimensions_dict[key] = {"dims": param_dim_names}

                common_params = set(param_dim_names) & set(dims)
                if not len(common_params):
                    parameter_dimensions_dict[key] = {
                        "dims": tuple(["unknown"])
                    }
                    continue

                param_dims = {kk: dims[kk] for kk in common_params}

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

                parameter_dimensions_dict[key] = {"dims": tuple(temp_dims)}

        # build coords, only some dims have coords that are not indexes
        coords = {}
        # edge case of dims that have only index or implied coords
        # this is justified, as this data should probably be in the parameter
        # files. netcdf neccisitates this to even carry the dimension data.
        # alternative would be to have a scalar, but in fact this is a
        # dimension of the data.
        coords["doy"] = np.arange(dims["ndoy"], dtype="int32") + 1

        for cc in ["nhm_id", "nhm_seg", "poi_gage_id", "doy"]:
            coord_to_dim = {
                "nhm_id": "nhru",
                "nhm_seg": "nsegment",
                "poi_gage_id": "npoigages",
                "doy": "ndoy",
            }
            dim_name = coord_to_dim[cc]
            if dim_name not in dims:
                continue

            dim_val = dims[dim_name]

            if cc in parameter_dict.keys():
                coords[cc] = parameter_dict.pop(cc)
            else:
                coord_to_dim = {"nhm_id": "nhru", "nhm_seg": "nsegment"}
                coords[cc] = np.arange(1, dim_val + 1)

            parameter_dimensions_dict[cc] = {"dims": (dim_name,)}

        for key, value in parameter_dimensions_dict.items():
            key_meta = deepcopy(meta.get_params(key)[key])
            _ = key_meta.pop("dims")
            parameter_dimensions_dict[key]["attrs"] = key_meta

        parameter_dimensions_dict["global"] = {
            "Description": "Parameter data for PRMS"
        }

        prms_params = PrmsParameters(
            dims=dims,
            coords=coords,
            data_vars=parameter_dict,
            metadata=parameter_dimensions_dict,
            validate=True,
        )

        return prms_params

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
