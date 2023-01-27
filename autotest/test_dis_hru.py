import hashlib

from pynhm.utils import DisHru
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

    # a lightweight regression test to make sure we are aware of changes.
    with open(out_file, "rb") as ff:
        md5sum_result = hashlib.md5(ff.read()).hexdigest()

    md5sum_answer = "eedac6070a6073b793c3c27d5a11c9f9"
    assert md5sum_result == md5sum_answer

    return
