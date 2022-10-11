import setuptools


def configuration():
    from numpy.distutils.misc_util import Configuration

    config = Configuration()
    config.add_extension(
        "PRMSGroundwater_f",
        sources=[
            "pynhm/hydrology/PRMSGroundwater.pyf",
            "pynhm/hydrology/PRMSGroundwater.f90",
        ],
    )
    # add more extensions here
    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration().todict())
