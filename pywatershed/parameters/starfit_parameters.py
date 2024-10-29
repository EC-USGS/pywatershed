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
from ..constants import nan, nat

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
    """Starfit parameter class.

    This parameter class provides STARFIT parameters to for modeling. This
    class does NOT calculate the parameters from inputs (e.g. as ISTARF-CONUS
    did using ResOpsUS), it simply provides the format for the model to get the
    the parameter data.

    The data supplied can come from whatever means. The method
    `from_istarf_conus_grand` uses existing ISTARF-CONUS and GRanD data to
    create a parameter object for the user.

    References:

    **ISTARF-CONUS (Inferred Storage Targets and Release Functions - Continental
    US)**: Sean W.D. Turner, Jennie Clarice Steyaert, Laura Condon,
    Nathalie Voisin, Water storage and release policies for all large
    reservoirs of conterminous United States, Journal of Hydrology,
    Volume 603, Part A, 2021, 126843, ISSN 0022-1694,
    https://doi.org/10.1016/j.jhydrol.2021.126843.
    https://zenodo.org/records/4602277

    **GRanD (Global Reservoir and Dam) database**: Lehner, Bernhard, Catherine
    Reidy Liermann, Carmen Revenga, Charles
    Vörösmarty, Balazs Fekete, Philippe Crouzet, Petra Döll et al. "High‐
    resolution mapping of the world's reservoirs and dams for sustainable
    river‐flow management." Frontiers in Ecology and the Environment 9, no. 9
    (2011): 494-502.
    https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/default/447819520013

    **ResOpsUS**: Steyaert, Jennie C., Laura E. Condon, Sean WD Turner, and
    Nathalie Voisin. "ResOpsUS, a dataset of historical reservoir operations
    in the contiguous United States." Scientific Data 9, no. 1 (2022): 34.
    https://zenodo.org/records/6612040

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

    """  # noqa: E501

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
        resops_domain: Union[str, pl.Path],
        istarf_conus: Union[str, pl.Path],
        grand_dams: Union[str, pl.Path],
        grand_ids: Union[str, pl.Path] = None,
        param_names: list = None,
    ) -> dict:
        """
        TODO: what are the netcdf parameter files? describe their format
        """
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
        files_directory: Union[pl.Path, str] = pl.Path("."),
        grand_ids: list = None,
    ):
        """Build parameter object from istarf-conus and the GRanD v1.3 sources.

        This returns the parameters for the STARFIT method. The parameters
        are in the original units of the method.

        Note that this method returns nan for the fields of start_time,
        end_time, and initial_storage. The user can edit the parameter set if
        she would like to change these with the following basic steps (outlined
        in an example notebook) export the parameters to and xarray data set
        via params.to_xr_ds(), then edit params using xarray, finally
        instantiate a parameter object from the xarray dataset using
        params = StarfitParameters.from_ds(param_ds). The units of
        initial_storage supplied should match the units of flow input to
        Starfit.

        Args:
        grand_file: a path to an existing dbf or shp file. If the file does not
            exist, an error will be thrown and you must download it manually
            at https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/de
        istarf_file: a path to an existing file. If file does not exist or is
            None then the file will be dowladed to files_directory. You can
            download the file yourself here
            https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/default/447819520013
        files_directory: A local directory where to download the file.
        grand_ids: a subset of grand_ids to keep.

        Examples:
        ---------

        Read the full ISTART-CONUS dataset, identify the "big sandy" reservoir
        by name to get its grand_id, then subset the parameters to this
        grand_id. This requires downloading the GRanD and ISTARF-CONUS datasets
        in advance and specifying the paths to those files.

        >>> import pywatershed as pws
        >>> grand_file = (
        ...     your_data_dir / "GRanD_Version_1_3/GRanD_reservoirs_v1_3.dbf"
        ... )
        >>> istarf_file = your_data_dir / "ISTARF-CONUS.csv"
        >>> sf_params = (
        ...     pws.parameters.StarfitParameters.from_istarf_conus_grand(
        ...         grand_file=grand_file, istarf_file=istarf_file
        ...     )
        ... )
        >>> grand_names = sf_params.parameters["GRanD_NAME"].tolist()
        >>> # where is Big Sandy?
        >>> big_sandy_index = [
        ...     ii
        ...     for ii, nn in enumerate(grand_names)
        ...     if "big sandy" in str(nn).lower()
        ... ][0]
        >>> big_sandy_grand_id = sf_params.parameters["grand_id"][
        ...     big_sandy_index
        ... ]
        >>> # get a parameter set with just the big sandy dike
        >>> sf_params = (
        ...     pws.parameters.StarfitParameters.from_istarf_conus_grand(
        ...         grand_file=grand_file,
        ...         istarf_file=istarf_file,
        ...         grand_ids=[big_sandy_grand_id],
        ...     )
        ... )

        """  # noqa: E501

        grand_ds = _get_grand(grand_file)
        istarf_ds = _get_istarf_conus(istarf_file, files_directory)

        common_grand_ids = list(
            set(grand_ds.grand_id.values).intersection(
                set(istarf_ds.grand_id.values)
            )
        )
        if grand_ids is not None:
            # make sure all requested grand_ids are in common_grand_ids
            avail_grand_ids = set(grand_ids).intersection(
                set(common_grand_ids)
            )
            if avail_grand_ids != set(grand_ids):
                unavail_grand_ids = set(grand_ids) - set(avail_grand_ids)
                msg = (
                    "The following requested grand_ids wer not available "
                    f"in the data sources: {unavail_grand_ids}"
                )
                raise ValueError(msg)
            else:
                common_grand_ids = list(avail_grand_ids)

        # <<
        common_grand_ids = np.array(common_grand_ids)
        grand_ds = grand_ds.where(
            grand_ds.grand_id.isin(common_grand_ids), drop=True
        )
        istarf_ds = istarf_ds.where(
            istarf_ds.grand_id.isin(common_grand_ids), drop=True
        )
        ds = xr.combine_by_coords([grand_ds, istarf_ds])

        nreservoirs = len(ds.nreservoirs)
        for vv in ["start_time", "end_time"]:
            ds[vv] = xr.Variable(
                "nreservoirs", np.array([nat] * nreservoirs, "<M8[ns]")
            )
        # <
        vv = "initial_storage"
        ds[vv] = xr.Variable("nreservoirs", np.array([nan] * nreservoirs))

        params_dd = DatasetDict.from_ds(ds)
        return StarfitParameters.from_dict(params_dd.data)


def _get_grand(grand_file):
    if grand_file is None:
        msg = (
            "You must acquire the GRanD file manually at\n"
            "https://ln.sync.com/dl/bd47eb6b0/anhxaikr-62pmrgtq-k44xf84f-pyz4atkm/view/default/447819520013"  # noqa: E501
        )
        raise IOError(msg)
    if not pl.Path(grand_file).exists():
        msg = f"the GRanD file {grand_file} does not exist."
        raise ValueError(msg)
    # check that it's a dbf or a shp file?

    cols_keep = ["GRAND_ID", "LONG_DD", "LAT_DD"]
    grand_ds = (
        gpd.read_file(grand_file)[cols_keep].to_xarray().drop_vars("index")
    )
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
        if not istarf_file_in_exists and istarf_file_in is not None:
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
