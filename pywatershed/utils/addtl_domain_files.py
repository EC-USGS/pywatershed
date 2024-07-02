# This utility gets additional doain files these files are not part of the pws
# repo.

import urllib.request as request
import zipfile
from shutil import rmtree

import pywatershed as pws

pkg_root_dir = pws.constants.__pywatershed_root__
addtl_domains_dir = pkg_root_dir / "data/pywatershed_addtl_domains"


def download(force=False):
    if force and addtl_domains_dir.exists():
        rmtree(addtl_domains_dir)

    if not addtl_domains_dir.exists():
        addtl_domains_url = (
            "https://github.com/EC-USGS/pywatershed/"
            "releases/download/1.1.0/pywatershed_addtl_domains.zip"
        )
        addtl_domains_file = (
            pkg_root_dir / "data/pywatershed_addtl_domains.zip"
        )
        request.urlretrieve(addtl_domains_url, addtl_domains_file)

        with zipfile.ZipFile(addtl_domains_file, "r") as zz:
            zz.extractall(pkg_root_dir / "data")

        assert addtl_domains_dir.exists()

    return


if __name__ == "__main__":
    download()
