import os
import platform
import warnings

import setuptools  # noqa
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

pyws_fortran = str(os.getenv("PYWS_FORTRAN"))

if pyws_fortran.lower() == "true":
    pyws_fortran = True
    if platform.system() == "Windows":
        pyws_fortran = False
        warnings.warn(
            "Fortran source compilation not enabled on Windows", Warning
        )

else:
    pyws_fortran = False


config = Configuration("pywatershed")

if pyws_fortran:
    source_dir_names = {
        "hydrology": ["prms_groundwater", "prms_canopy", "prms_channel"]
    }
    for dir_name, names in source_dir_names.items():
        for name in names:
            config.add_extension(
                f"{name}_f",
                sources=[
                    f"pywatershed/{dir_name}/{name}.pyf",
                    f"pywatershed/{dir_name}/{name}.f90",
                ],
            )

setup(**config.todict(), packages=setuptools.find_packages())
