import os
import sys

import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = ".."
elif cwd.endswith("prmsNHMpy"):
    sys.path.append(".")
    rel_path = "."

from pynhm import PrmsParameters

base_path = os.path.join(rel_path, "test_models")
parameter_files = (
    os.path.join(base_path, "merced", "input", "mercdIDE.param"),
    os.path.join(base_path, "merced", "input", "mercdXYZ.param"),
    os.path.join(base_path, "acf", "input", "acf.params"),
)
parameter_keys = [
    os.path.basename(fpth.replace(".param", "").replace(".params", ""))
    for fpth in parameter_files
]

test_dimensions = {
    parameter_keys[0]: {
        "nhru": 90,
        "ngw": 90,
        "nssr": 90,
        "ndays": 366,
        "nrain": 26,
        "ntemp": 26,
    },
    parameter_keys[1]: {
        "nhru": 90,
        "ngw": 90,
        "nssr": 90,
        "ndays": 366,
        "nrain": 26,
        "ntemp": 26,
    },
}


@pytest.mark.parametrize(
    "parameter_key, parameter_file", zip(parameter_keys, parameter_files)
)
def test_parameter_read(parameter_key, parameter_file):
    print(f"parsing...'{parameter_file}'")

    p_parameters = PrmsParameters(parameter_file)
    dimensions = p_parameters.get_dimensions
    param_data = p_parameters.get_parameter_data

    if param_data in tuple(test_dimensions.keys()):
        test_dimension = test_dimensions[parameter_key]
        for key in test_dimensions[parameter_key].keys():
            assert dimensions[key] == test_dimension[key], (
                f"{parameter_key}: parsed value '{dimensions[key]}' "
                + f"does not equal '{test_dimension[key]}'"
            )
    print(f"success parsing...'{parameter_file}'")

    return


if __name__ == "__main__":

    for parameter_key, parameter_file in zip(parameter_keys, parameter_files):
        test_parameter_read(parameter_key, parameter_file)
