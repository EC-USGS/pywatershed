import numpy as np

from pynhm import CsvFile


def test_csv_to_netcdf(csv_files):
    nc_path = csv_files.with_suffix(".nc")
    CsvFile(csv_files).to_netcdf(nc_path)
    assert nc_path.exists()


def test_csv_to_previous_netcdf(csv_files_prev):
    nc_name = csv_files_prev.with_suffix("").name + "_prev.nc"
    nc_path = csv_files_prev.parent / nc_name
    csv = CsvFile(csv_files_prev)
    csv._get_data(variable_name=nc_path.stem)  # why so private?
    orig_dates = csv._data["date"].copy()
    csv._data = np.roll(csv._data, 1, axis=0)
    csv._data["date"] = orig_dates
    # Here we will eventually want to supply the desired initial conditions
    for hh in range(len(csv._data[0])):
        if hh == 0:
            continue
        csv._data[0][hh] = np.zeros([1])[0]
    csv.to_netcdf(nc_path)
    assert nc_path.exists()
