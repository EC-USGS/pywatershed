from filelock import FileLock
from warnings import warn

import pytest
from filelock import FileLock
from prms_convert_to_netcdf import convert_csv_to_nc, convert_soltab_to_nc
from prms_diagnostic_variables import (
    diagnose_final_vars_to_nc,
    diagnose_simple_vars_to_nc,
)

import pywatershed as pws


@pytest.fixture
def netcdf_file(control_csv_file):
    """Convert CSV files from model output to NetCDF"""
    control_file = control_csv_file[0]
    csv_file = control_csv_file[1]

    var_name = csv_file.stem
    data_dir = csv_file.parent
    convert_csv_to_nc(var_name, data_dir)

    success = diagnose_simple_vars_to_nc(var_name, data_dir, control_file)

    if not success:
        # should this fail or not?
        # pytest.skip(f"Unable to diagnose {var_name}")
        assert False, "Unable to diagnose {var_name}"

    return f"{var_name}.nc"


def make_netcdf_files(netcdf_file):
    print(f"Created NetCDF from CSV: {netcdf_file}")


@pytest.fixture()
def soltab_netcdf_file(tmp_path_factory, control_soltab_file):
    """Convert soltab files to NetCDF, one file for each variable"""
    control_file = control_soltab_file[0]
    soltab_file = control_soltab_file[1]
    domain_dir = soltab_file.parent
    control = pws.Control.load_prms(control_file, warn_unused_options=False)
    output_dir = control_file.parent / control.options["netcdf_output_dir"]

    convert_soltab_to_nc(soltab_file, output_dir, control_file, domain_dir)


def make_soltab_netcdf_files(soltab_netcdf_file):
    print(f"Creating NetCDF files for soltab file {soltab_netcdf_file}")


@pytest.fixture(scope="session")
def final_netcdf_file(tmp_path_factory, control_final_file):
    """Create NetCDF files that depend on multiple other NetCDFs"""

    control_file = control_final_file[0]
    final_file = control_final_file[1]
    domain_dir = final_file.parent
    var_name = final_file.name
    control = pws.Control.load_prms(control_file, warn_unused_options=False)
    output_dir = control_file.parent / control.options["netcdf_output_dir"]

    root_tmpdir = tmp_path_factory.getbasetemp().parent
    with FileLock(root_tmpdir / "final_nc.lock"):
        yield final_file  # do this in session cleanup

        success = diagnose_final_vars_to_nc(
            var_name, output_dir, control_file, domain_dir
        )

        if not success:
            warn(
                "make_final_netcdf_files False PASS above: "
                f"unable to diagnose {final_file}"
            )

    return final_file  # dosent really matter


def make_final_netcdf_files(final_netcdf_file):
    print(f"Creating final NetCDF file {final_netcdf_file}")
