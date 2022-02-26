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

from pynhm.utils import PrmsParameters


def test_parameter_read(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")
    ans_dict_domain = domain["test_ans"]["parameter_read"]

    p_parameters = PrmsParameters(parameter_file)
    dimensions = p_parameters.get_dimensions
    # param_data = p_parameters.get_parameter_data

    # check dimensions
    for param_key in dimensions.keys():
        if param_key in ans_dict_domain.keys():
            msg = (
                f"{param_key}: parsed value '{dimensions[param_key]}' "
                + f"does not equal '{ans_dict_domain[param_key]}'"
            )
            # print(f'"{param_key}": {dimensions[param_key]},')
            assert dimensions[param_key] == ans_dict_domain[param_key], msg

    print(f"success parsing...'{parameter_file}'")

    return
