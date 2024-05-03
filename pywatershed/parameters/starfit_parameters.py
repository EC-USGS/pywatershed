import pathlib as pl
import tempfile
import urllib
from typing import Union
from warnings import warn

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from ..base.data_model import DatasetDict
from ..base.parameters import Parameters
from ..constants import fileish

# from ..hydrology.starfit import Starfit
# from pywatershed import Starfit is circular so copy the needed info
starfit_param_names = (
    "grand_id",
    "GRanD_NAME",
    "initial_storage",
    "start_time",
    "end_time",
    "inflow_mean",
    "NORhi_min",
    "NORhi_max",
    "NORhi_alpha",
    "NORhi_beta",
    "NORhi_mu",
    "NORlo_min",
    "NORlo_max",
    "NORlo_alpha",
    "NORlo_beta",
    "NORlo_mu",
    "Release_min",
    "Release_max",
    "Release_alpha1",
    "Release_alpha2",
    "Release_beta1",
    "Release_beta2",
    "Release_p1",
    "Release_p2",
    "Release_c",
    "GRanD_CAP_MCM",
    "Obs_MEANFLOW_CUMECS",
    "GRanD_MEANFLOW_CUMECS",
)


class StarfitParameters(Parameters):
    """
    Starfit parameter class

    The GRanD data
    https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/default/447819520013  # noqa

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

        if param_names:
            params_dd = params_dd.subset(param_names)

        return StarfitParameters.from_dict(params_dd.data)

    @staticmethod
    def from_istarf_conus_grand(
        grand_file: Union[pl.Path, str],
        istarf_file: Union[pl.Path, str] = None,
        files_directory: Union[pl.Path, str] = None,
        grand_ids: list = None,
    ):
        """Build parameter object from istarf-conus and the GRanD v1.3 sources.

        Args:
        istarf_file: a path to an existing file. If file does not exist or is
            None then the file will be dowladed to files_directory
        grand_file: a path to an existing dbf or shp file. If the file does not
            exist, an error will be thrown and you must download it manually
            from
             https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/default/447819520013  # noqa

        """
        grand_ds = _get_grand(grand_file)
        istarf_ds = _get_istarf_conus(istarf_file, files_directory)

        common_grand_ids = list(
            set(grand_ds.grand_id.values).intersection(
                set(istarf_ds.grand_id.values)
            )
        )
        grand_ds = grand_ds.where(
            grand_ds.grand_id.isin(common_grand_ids), drop=True
        )
        istarf_ds = istarf_ds.where(
            istarf_ds.grand_id.isin(common_grand_ids), drop=True
        )

        ds = xr.combine_by_coords([grand_ds, istarf_ds])

        params_dd = DatasetDict.from_ds(ds)
        return StarfitParameters.from_dict(params_dd.data)


def _get_grand(grand_file):
    if not pl.Path(grand_file).exists():
        msg = f"the GRanD file {grand_file} does not exist."
        raise ValueError(msg)
    # check that it's a dbf or a shp file?

    cols_keep = ["GRAND_ID", "LONG_DD", "LAT_DD"]
    grand_ds = gpd.read_file(grand_file)[cols_keep].to_xarray().drop("index")
    grand_ds = grand_ds.rename(
        {"GRAND_ID": "grand_id", "index": "nreservoirs"}
    ).set_coords("grand_id")
    return grand_ds


def _get_istarf_conus(istarf_file, files_directory):
    files_directory_in = files_directory
    istarf_file_in = istarf_file

    files_directory_in_exists = (
        files_directory_in is not None and pl.Path(files_directory_in).exists()
    )
    if files_directory_in is None or not files_directory_in_exists:
        files_directory = pl.Path(tempfile.mkdtemp())

    istarf_file_in_exists = (
        istarf_file_in is not None and pl.Path(istarf_file_in).exists()
    )
    if istarf_file_in is None or not istarf_file_in_exists:
        # dowload source to files_directory
        istarf_url = (
            "https://zenodo.org/records/4602277/files/"
            "ISTARF-CONUS.csv?download=1"
        )
        istarf_file = files_directory / "ISTARF-CONUS.csv"
        if not istarf_file_in_exists:
            warn(
                "The specified istarf_file does not exist: "
                f"{istarf_file_in}"
            )
        if files_directory_in is not None and not files_directory_in_exists:
            warn(
                "The specified files_directory does not exist: "
                f" {files_directory_in}"
            )
        print(f"Downloading and saving ISTARF-CONUS.csv to {istarf_file}")
        urllib.request.urlretrieve(istarf_url, istarf_file)

    # <
    istarf_ds = pd.read_csv(istarf_file).to_xarray()
    rename_map = {
        "GRanD_ID": "grand_id",
        "index": "nreservoirs",
        "GRanD_MEANFLOW_CUMECS": "inflow_mean",
    }
    istarf_ds = istarf_ds.rename(rename_map)
    istarf_ds = istarf_ds.set_coords(["nreservoirs", "grand_id"])
    # drop variables not in the starfit parameters
    data_vars = list(istarf_ds.variables)
    for vv in data_vars:
        if vv not in starfit_param_names:
            # print(vv)  # fit, match, nreservoirs, GRanD_MEANFLOW_CUMECS
            del istarf_ds[vv]

    return istarf_ds
