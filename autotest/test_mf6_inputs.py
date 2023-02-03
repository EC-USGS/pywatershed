import hashlib

from pynhm.utils import BcSeg, DisHru, DisSeg
from pynhm.constants import __pynhm_root__

# not currently in repo or being used but might be good to add
# shape_file = (
#     "/Users/jamesmcc/usgs/data/pynhm/20220209_gm_delaware_river"
#     "/GIS_simple/HRU_subset.shp"
# )
param_file = __pynhm_root__ / "../test_data/drb_2yr/myparam.param"


def test_dis_hru(tmp_path):
    dis = DisHru(param_file=param_file)  # , hru_shapefile=shape_file)
    out_file = tmp_path / "dis_hru"
    dis.write(out_file)
    assert out_file.exists()
    print(out_file)

    # a lightweight regression test to make sure we are aware of changes.
    with open(out_file, "rb") as ff:
        md5sum_result = hashlib.md5(ff.read()).hexdigest()

    md5sum_answer = "eedac6070a6073b793c3c27d5a11c9f9"
    assert md5sum_result == md5sum_answer

    return


def test_dis_seg(tmp_path):
    dis = DisSeg(param_file=param_file)  # , hru_shapefile=shape_file)
    out_file = tmp_path / "dis_seg"
    dis.write(out_file)
    assert out_file.exists()

    # a lightweight regression test to make sure we are aware of changes.
    with open(out_file, "rb") as ff:
        md5sum_result = hashlib.md5(ff.read()).hexdigest()

    # md5sum_answer = "eedac6070a6073b793c3c27d5a11c9f9"
    # assert md5sum_result == md5sum_answer

    return


def test_bc_seg(tmp_path):
    netcdf_bc_file = (
        __pynhm_root__ / "../test_data/drb_2yr/output/gwres_flow_vol.nc"
    )
    bc_seg = BcSeg(netcdf_bc_file, param_file=param_file, print_input=True)
    out_file = tmp_path / "bc_seg_gwres_flow_vol"
    bc_seg.write(out_file)
    assert out_file.exists()

    # a lightweight regression test to make sure we are aware of changes.
    with open(out_file, "rb") as ff:
        md5sum_result = hashlib.md5(ff.read()).hexdigest()

    md5sum_answer = "808382240973b425b5cbe9b1ce72684d"
    assert md5sum_result == md5sum_answer

    return
