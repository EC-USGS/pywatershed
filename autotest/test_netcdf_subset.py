import pathlib as pl

import pytest
import xarray as xr

import pywatershed as pws

this_dir = pl.Path(__file__).parent

file_types = (
    this_dir / "../test_data/drb_2yr/parameters_dis_hru.nc",
    this_dir / "../test_data/drb_2yr/tmax.nc",
    this_dir / "../test_data/drb_2yr/output/sroff.nc",
)


@pytest.fixture(scope="function")
def nhm_ids():
    param_file = pl.Path("../test_data/drb_2yr/parameters_dis_hru.nc")
    param_ds = xr.open_dataset(param_file)
    return param_ds.nhm_id.values[(0, 10, 100),]


@pytest.mark.domainless
@pytest.mark.parametrize("file_type", file_types)
def test_subset_netcdf_file(file_type, nhm_ids, tmp_path):
    # do time in the test on the dataset not the file
    new_file = tmp_path / "new.nc"
    pws.utils.netcdf_utils.subset_netcdf_file(
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


@pytest.mark.domainless
@pytest.mark.parametrize("file_type", file_types)
def test_subset_xr_ds(file_type, nhm_ids, tmp_path):
    start_time = end_time = None
    old_ds = xr.open_dataset(file_type)
    if "time" in old_ds.variables:
        start_time = old_ds.time[100]
        end_time = old_ds.time[125]
        # create a variables without time and nhm_id dimensions
        data_var_names = list(old_ds.data_vars)
        old_ds["dum_var"] = old_ds[data_var_names[0]].isel(time=0)
        old_ds["some_var"] = old_ds[data_var_names[0]][:, 0].squeeze()

    new_ds = pws.utils.netcdf_utils.subset_xr(
        old_ds,
        start_time=start_time,
        end_time=end_time,
        coord_dim_name="nhm_id",
        coord_dim_values_keep=nhm_ids,
    )

    assert set(new_ds.nhm_id.values) == set(nhm_ids)

    assert set(new_ds.variables) == set(old_ds.variables)

    for vv in new_ds:
        assert new_ds[vv].dims == old_ds[vv].dims

    if start_time is not None:
        assert len(new_ds.time == 26)


@pytest.mark.domainless
@pytest.mark.parametrize("file_type", file_types[1:])
def test_subset_xr_da(file_type, nhm_ids, tmp_path):
    start_time = end_time = None
    old_da = xr.open_dataarray(file_type)
    start_time = old_da.time[100]
    end_time = old_da.time[125]

    new_da = pws.utils.netcdf_utils.subset_xr(
        old_da,
        start_time=start_time,
        end_time=end_time,
        coord_dim_name="nhm_id",
        coord_dim_values_keep=nhm_ids,
    )

    assert set(new_da.nhm_id.values) == set(nhm_ids)

    assert new_da.dims == old_da.dims

    assert len(new_da.time == 26)
