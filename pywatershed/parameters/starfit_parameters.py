import numpy as np

from ..base.data_model import DatasetDict
from ..base.parameters import Parameters
from ..constants import fileish


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
        encoding: dict = {},
    ) -> "StarfitParameters":
        super().__init__(
            dims=dims,
            coords=coords,
            data_vars=data_vars,
            metadata=metadata,
            encoding=encoding,
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
        resops_dd = DatasetDict.from_netcdf(resops_domain)
        istarf_conus_dd = DatasetDict.from_netcdf(istarf_conus)
        istarf_rename = {"GRanD_ID": "grand_id"}
        istarf_conus_dd.rename_dim(istarf_rename)
        istarf_conus_dd.rename_var(istarf_rename)
        wh_subset = np.where(
            np.isin(
                istarf_conus_dd.coords["grand_id"],
                resops_dd.coords["grand_id"],
            )
        )
        istarf_conus_dd.subset_on_coord("grand_id", wh_subset)

        grand_dams_dd = DatasetDict.from_netcdf(grand_dams)
        grand_rename = {"GRAND_ID": "grand_id"}
        grand_dams_dd.rename_dim(grand_rename)
        grand_dams_dd.rename_var(grand_rename)
        wh_subset = np.where(
            np.isin(
                grand_dams_dd.coords["grand_id"],
                resops_dd.coords["grand_id"],
            )
        )
        grand_dams_dd.subset_on_coord("grand_id", wh_subset)

        istarf_conus_dd.drop_var("subset_inds")
        grand_dams_dd.drop_var("subset_inds")

        params_dd = DatasetDict.merge(
            resops_dd, istarf_conus_dd, grand_dams_dd
        )
        # probably should have named the spatial dim nreservoirs when
        # I created the netcdf file
        _ = params_dd.rename_dim({"grand_id": "nreservoirs"})

        return StarfitParameters.from_dict(params_dd.data)
