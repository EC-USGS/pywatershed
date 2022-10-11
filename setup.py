import setuptools  # noqa

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

config = Configuration("pynhm")
config.add_extension(
    "PRMSGroundwater_f",
    sources=[
        "pynhm/hydrology/PRMSGroundwater.pyf",
        "pynhm/hydrology/PRMSGroundwater.f90",
    ],
)
# add more f2py extensions here

setup(**config.todict(), packages=setuptools.find_packages())
