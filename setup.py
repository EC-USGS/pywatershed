import os
import platform
import setuptools  # noqa
import warnings

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

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
        "hydrology": ["PRMSGroundwater", "PRMSCanopy", "PRMSChannel"]
    }
    for dir, names in source_dir_names.items():
        for name in names:
            config.add_extension(
                f"{name}_f",
                sources=[
                    f"pywatershed/{dir}/{name}.pyf",
                    f"pywatershed/{dir}/{name}.f90",
                ],
            )


setup(**config.todict(), packages=setuptools.find_packages())
