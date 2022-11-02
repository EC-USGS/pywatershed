from pynhm.utils import DisHru

shape_file = (
    "/Users/jamesmcc/usgs/data/pynhm/20220209_gm_delaware_river"
    "/GIS_simple/HRU_subset.shp"
)
param_file = "/Users/jamesmcc/usgs/pynhm/test_data/drb_2yr/myparam.param"


def test_dis_hru(tmp_path):
    dis = DisHru(param_file=param_file, hru_shapefile=shape_file)

    out_file = tmp_path / "dis_hru"
    dis.write(out_file)
    assert out_file.exists()

    return
