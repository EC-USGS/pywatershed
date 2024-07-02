# This utility gets GIS files needed for certain visualizations and
# calculations. These GIS files are not part of the pws repo.

import urllib.request as request
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
            "releases/download/1.1.0/pywatershed_gis.zip"
        )
        gis_file = pkg_root_dir / "data/pywatershed_gis.zip"
        request.urlretrieve(gis_url, gis_file)

        with zipfile.ZipFile(gis_file, "r") as zz:
            zz.extractall(pkg_root_dir / "data")

        assert gis_dir.exists()

    return


if __name__ == "__main__":
    download()
