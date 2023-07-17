import urllib
import zipfile
from shutil import rmtree

import pywatershed as pws

pkg_root_dir = pws.constants.__pywatershed_root__
gis_dir = pkg_root_dir / "data/pywatershed_gis"


def download(force=False):
    if force and gis_dir.exists():
        rmtree(gis_dir)

    if not gis_dir.exists():
        gis_url = (
            "https://github.com/EC-USGS/pywatershed/"
            "releases/download/v2022.0.1/pynhm_gis.zip"
        )
        gis_file = pkg_root_dir / "data/pynhm_gis.zip"
        urllib.request.urlretrieve(gis_url, gis_file)

        with zipfile.ZipFile(gis_file, "r") as zz:
            zz.extractall(pkg_root_dir / "data")

        (pkg_root_dir / "data/pynhm_gis").rename(gis_dir)

    return
