import os
import setuptools  # noqa

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

pynhm_fortran = str(os.getenv("PYNHM_FORTRAN"))
if pynhm_fortran.lower() == "false":
    pynhm_fortran = False
else:
    pynhm_fortran = True

config = Configuration("pynhm")

if pynhm_fortran:
    config.add_extension(
        "PRMSGroundwater_f",
        sources=[
            "pynhm/hydrology/PRMSGroundwater.pyf",
            "pynhm/hydrology/PRMSGroundwater.f90",
        ],
    )
    # add more f2py extensions here

setup(**config.todict(), packages=setuptools.find_packages())
