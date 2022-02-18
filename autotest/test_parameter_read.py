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

base_path = os.path.join(rel_path, "test_data")

parameter_file_dict = {
    'drb_2yr': os.path.join(base_path, "drb_2yr", "myparam.param"),}

ans_dict = {
    'drb_2yr': {
        "nhru": 765,
        "ngw": 765,
        "nssr": 765, }, }


@pytest.mark.parametrize("domain_key", parameter_file_dict.keys())
def test_parameter_read(domain_key):
    parameter_file = parameter_file_dict[domain_key]
    print(f"parsing...'{parameter_file}'")
    ans_dict_domain = ans_dict[domain_key]

    p_parameters = PrmsParameters(parameter_file)
    dimensions = p_parameters.get_dimensions
    param_data = p_parameters.get_parameter_data
    
    # check dimensions
    for param_key in dimensions.keys():
        if param_key in ans_dict_domain.keys():
            msg = (
                f"{param_key}: parsed value '{dimensions[param_key]}' "
                + f"does not equal '{ans_dict_domain[param_key]}'")
            # print(f'"{param_key}": {dimensions[param_key]},')
            assert dimensions[param_key] == ans_dict_domain[param_key], msg
            
    print(f"success parsing...'{parameter_file}'")

    return


if __name__ == "__main__":

    for parameter_key, parameter_file in zip(parameter_keys, parameter_files):
        test_parameter_read(parameter_key, parameter_file)
