from copy import deepcopy

import netCDF4 as nc4  # remove
import numpy as np

from ..base.parameters import Parameters
from ..constants import fileish
from ..base import data_model as dm


class StarfitParameters(Parameters):
    """
    Starfit parameter class

    The GRanD data
    https://sedac.ciesin.columbia.edu/data/set/grand-v1-dams-rev01

    The istarf data
    https://zenodo.org/record/4602277#.ZCtYj-zMJqs

    The resops data
    https://zenodo.org/record/5893641#.ZCtakuzMJqs

    # add citiatons. add this information to the starfit model too
    # add a working example

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
        metadata: dict,
    ) -> "StarfitParameters":
        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
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
        istarf_conus_dd = dm.nc4_ds_to_dd(istarf_conus)  # from_nc
        resops_dd = dm.nc4_ds_to_dd(resops_domain)  # from_netcdf
        params_dd = dm.DatasetDict.merge(istarf_conus_dd, resops_dd)

        asdf

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
            dims={}, coords={}, data_vars=param_dict, metadata={}
        )

    def subset(self, **kwargs) -> Parameters:
        return deepcopy(self)

    def get_parameters(self, names, **kwargs) -> Parameters:
        if not isinstance(names, list):
            names = [names]
        return {name: self._data_vars[name] for name in list(names)}
