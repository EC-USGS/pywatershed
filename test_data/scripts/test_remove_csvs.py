csvs_to_keep = ["hru_ppt", "intcp_stor", "potet", "gwres_stor"]


def test_csv_to_netcdf(csv_file):
    if csv_file.with_suffix("").name not in (csvs_to_keep):
        csv_file.unlink()
