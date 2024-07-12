import pytest
import xarray as xr

import pywatershed as pws

file_types = (
    "parameters_dis_hru.nc",
    "tmax.nc",
    "output/sroff.nc",
)


@pytest.fixture(scope="function")
def nhm_ids(simulation):
    domain_name = simulation["name"].split(":")[0]
    if domain_name not in ["hru_1", "drb_2yr"]:
        pytest.skip("Only test_netcdf_subset hru_1 and drb_2yr")
    if domain_name == "hru_1":
        subset_inds = (0,)
    else:
        subset_inds = (0, 10, 100)

    param_file = simulation["dir"] / "parameters_dis_hru.nc"
    param_ds = xr.open_dataset(param_file)
    return param_ds.nhm_id.values[subset_inds,]


@pytest.mark.parametrize("file_type", file_types)
def test_subset_netcdf_file(simulation, file_type, nhm_ids, tmp_path):
    # do time in the test on the dataset not the file
    old_file = simulation["dir"] / file_type
    new_file = tmp_path / "new.nc"
    pws.utils.netcdf_utils.subset_netcdf_file(
        old_file,
        new_file,
        coord_dim_name="nhm_id",
        coord_dim_values_keep=nhm_ids,
    )

    old_ds = xr.open_dataset(old_file)
    new_ds = xr.open_dataset(new_file)

    assert set(new_ds.nhm_id.values) == set(nhm_ids)

    assert set(new_ds.variables) == set(old_ds.variables)

    for vv in new_ds:
        assert new_ds[vv].dims == old_ds[vv].dims


@pytest.mark.parametrize("file_type", file_types)
def test_subset_xr_ds(simulation, file_type, nhm_ids, tmp_path):
    start_time = end_time = None
    old_file = simulation["dir"] / file_type
    old_ds = xr.open_dataset(old_file)
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


@pytest.mark.parametrize("file_type", file_types[1:])
def test_subset_xr_da(simulation, file_type, nhm_ids, tmp_path):
    start_time = end_time = None
    old_file = simulation["dir"] / file_type
    old_da = xr.open_dataarray(old_file)
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
