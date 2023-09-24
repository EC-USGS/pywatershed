import pathlib as pl

import pywatershed as pws
from pywatershed import CsvFile, Soltab

"""This module is intended to aid/make consistent conversions of PRMS files"""


def convert_csv_to_nc(
    var_name: str, data_dir: pl.Path, output_dir: pl.Path = None
):
    """Convert PRMS CSV files to netcdf.

    Args:
        var_name: str name of the variable to create
        data_dir: where the csv file is found and the netcdf file will be
            written (could add argument to output to a differnt dir)
    """
    if output_dir is None:
        output_dir = data_dir
    csv_path = data_dir / f"{var_name}.csv"
    nc_path = output_dir / f"{var_name}.nc"
    CsvFile(csv_path).to_netcdf(nc_path)
    assert nc_path.exists()


def convert_soltab_to_nc(
    soltab_file: pl.Path, output_dir: pl.Path, domain_dir: pl.Path
):
    """Convert soltab files to NetCDF, one file for each variable

    The inputs files are "soltab_debug" specifically written by the PRMS5.2.1
    in the pywatershed repository.

    Args:
        soltab_file: the pl.Path of the  soltab_debug file.
        output_dir: pl.Path where the netcdf file output will be written
        domain_dir: defaults to the parent dir of soltab_file, the pl.Path
            where domain files (parameter & control) for the domain can be
            found
    """
    if domain_dir is None:
        domain_dir = soltab_file.parent

    params = pws.parameters.PrmsParameters.load(domain_dir / "myparam.param")
    nhm_ids = params.parameters["nhm_id"]

    soltab = Soltab(soltab_file, nhm_ids=nhm_ids)
    soltab.to_netcdf(output_dir=output_dir)
    print(f"Created NetCDF files from soltab file {soltab_file}:")

    for var in soltab.variables:
        nc_path = output_dir / f"{var}.nc"
        assert nc_path.exists()
