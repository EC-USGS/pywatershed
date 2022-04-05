from pynhm.utils.netcdf_utils import NetCdfRead


def test_netcdf(domain):
    output_files = list(domain["prms_outputs"].values())
    output_path = output_files[-1]
    nc_pth = output_path.with_suffix(".nc")
    variable = nc_pth.stem

    nc_data = NetCdfRead(nc_pth)
    ntimes, nhru = nc_data.ntimes, nc_data.nhru
    shape = (nhru,)

    for idx in range(ntimes):
        arr = nc_data.advance(variable)
        assert (
            arr.shape == shape
        ), f"shape is {arr.shape} but should be {shape}"

    shape = (ntimes, nhru)
    arr = nc_data.get_data(variable)
    assert arr.shape == shape, f"shape is {arr.shape} but should be {shape}"
