from pywatershed.utils.netcdf_utils import NetCdfRead


def test_netcdf(domain):
    # output_files = list(domain["prms_outputs"].values())
    variable = "gwres_stor"
    output_dir = domain["prms_output_dir"]
    nc_pth = output_dir / f"{variable}.nc"

    nc_data = NetCdfRead(nc_pth)
    ntimes, nhru = nc_data.ntimes, nc_data.nhru
    shape = (nhru,)

    for idx in range(ntimes):
        arr = nc_data.advance(variable)
        assert arr.shape == shape, f"shape is {arr.shape} but should be {shape}"

    shape = (ntimes, nhru)
    arr = nc_data.get_data(variable)
    assert arr.shape == shape, f"shape is {arr.shape} but should be {shape}"
