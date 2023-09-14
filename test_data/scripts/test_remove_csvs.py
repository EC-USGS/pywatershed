csvs_to_keep = ["hru_ppt", "intcp_stor", "potet", "gwres_stor"]


def test_csv_to_netcdf(csv_files):
    if csv_files.with_suffix("").name not in (csvs_to_keep):
        csv_files.unlink()
