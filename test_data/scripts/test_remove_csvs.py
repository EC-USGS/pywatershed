from .conftest import previous_vars

csvs_to_keep = ["hru_ppt", "intcp_stor", "potet", "gwres_stor"]


def test_csv_to_netcdf(csv_files):
    print((csvs_to_keep + previous_vars))
    if csv_files.with_suffix("").name not in (csvs_to_keep + previous_vars):
        csv_files.unlink()
