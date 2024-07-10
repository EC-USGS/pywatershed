import pathlib as pl

import pytest
import pywatershed as pws
import xarray as xr

file_types = (
    pl.Path("../test_data/drb_2yr/parameters_dis_hru.nc"),
    pl.Path("../test_data/drb_2yr/tmax.nc"),
    pl.Path("../test_data/drb_2yr/output/sroff.nc"),
)


@pytest.fixture(scope="function")
def nhm_ids():
    param_file = pl.Path("../test_data/drb_2yr/parameters_dis_hru.nc")
    param_ds = xr.open_dataset(param_file)
    return param_ds.nhm_id.values[(0, 10, 100),]


@pytest.mark.domainless
@pytest.mark.parametrize("file_type", file_types)
def test_subset_netcdf_file(file_type, nhm_ids, tmp_path):
    new_file = tmp_path / "new.nc"
    pws.utils.netcdf_utils.subset_file(
        file_type,
        new_file,
        coord_dim_name="nhm_id",
        coord_dim_values_keep=nhm_ids,
    )

    old_ds = xr.open_dataset(file_type)
    new_ds = xr.open_dataset(new_file)

    assert set(new_ds.nhm_id.values) == set(nhm_ids)

    assert set(new_ds.variables) == set(old_ds.variables)

    for vv in new_ds:
        assert new_ds[vv].dims == old_ds[vv].dims
