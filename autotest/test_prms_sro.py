import os
import sys

import pytest

cwd = os.getcwd()
if cwd.endswith("autotest"):
    sys.path.append("..")
    rel_path = ".."
elif cwd.endswith("pynhm"):
    sys.path.append(".")
    rel_path = "."

import pynhm

base_path = os.path.join(rel_path, "test_data")

parameter_file_dict = {
    "box_01": os.path.join(base_path, "box_01", "box_01.params"),
}


@pytest.mark.parametrize("domain_key", parameter_file_dict.keys())
def test_prms_surface_runoff(domain_key):
    parameter_file = parameter_file_dict[domain_key]
    print(f"parsing...'{parameter_file}'")

    p_parameters = pynhm.utils.PrmsParameters(parameter_file)
    dimensions = p_parameters.get_dimensions
    param_data = p_parameters.get_parameter_data

    sro_vol = pynhm.runoff.prmsSurfaceRunoff(1, None, param_data)

    return


if __name__ == "__main__":

    for parameter_key in list(parameter_file_dict.keys()):
        test_prms_surface_runoff(parameter_key)
