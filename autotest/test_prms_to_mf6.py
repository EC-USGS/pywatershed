# import hashlib

from pynhm.utils.prms_to_mf6 import MMRToMF6

# not currently in repo or being used but might be good to add
# shape_file = (
#     "/Users/jamesmcc/usgs/data/pynhm/20220209_gm_delaware_river"
#     "/GIS_simple/HRU_subset.shp"
# )


def test_mmr_to_mf6(domain, tmp_path):

    param_file = domain["param_file"]
    control_file = domain["control_file"]

    mmr = MMRToMF6(
        param_file=param_file,
        control_file=control_file,
        output_dir=tmp_path,
        inflow_dir=control_file.parent / "output",
        # , hru_shapefile=shape_file)
    )

    # mmr.write(tmp_dir)

    output_files = sorted(tmp_path.glob("*"))
    assert len(output_files) == 8

    # a lightweight regression test to make sure we are aware of changes.
    # with open(out_file, "rb") as ff:
    #     md5sum_result = hashlib.md5(ff.read()).hexdigest()

    # md5sum_answer = "eedac6070a6073b793c3c27d5a11c9f9"
    # assert md5sum_result == md5sum_answer

    return


# def test_dis_hru(tmp_path):
#     dis = DisHru(param_file=param_file)  # , hru_shapefile=shape_file)
#     out_file = tmp_path / "dis_hru"
#     dis.write(out_file)
#     assert out_file.exists()
#     print(out_file)

#     # a lightweight regression test to make sure we are aware of changes.
#     with open(out_file, "rb") as ff:
#         md5sum_result = hashlib.md5(ff.read()).hexdigest()

#     md5sum_answer = "eedac6070a6073b793c3c27d5a11c9f9"
#     assert md5sum_result == md5sum_answer

#     return
