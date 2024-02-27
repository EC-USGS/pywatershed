csvs_to_keep = ["hru_ppt", "intcp_stor", "potet", "gwres_stor"]


def test_remove_csv_files(control_csv_file):
    csv_file = control_csv_file[1]
    if csv_file.with_suffix("").name not in (csvs_to_keep):
        csv_file.unlink()
