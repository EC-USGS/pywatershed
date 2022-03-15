import os
import sys

import pytest

from pynhm.utils import PrmsParameters
from utils import assert_or_print


@pytest.fixture
def canopy_parameters():
    return tuple(("unknown", "srain_intcp", "wrain_intcp", "snow_intcp"))


def test_parameter_read(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    p_parameters = PrmsParameters(parameter_file)
    dimensions = p_parameters.get_dimensions
    # param_data = p_parameters.get_parameter_data

    # check dimensions
    answers = domain["test_ans"]["parameter_read"]
    results = {
        key: val for key, val in dimensions.items() if key in answers.keys()
    }
    assert_or_print(results, answers, print_ans=domain["print_ans"])

    print(f"success parsing...'{parameter_file}'")

    return


def test_parameter_canopy_subset(domain, canopy_parameters):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    p_parameters = PrmsParameters(parameter_file)

    canopy_subset = p_parameters.get_parameters(canopy_parameters)

    assert (
        canopy_subset["unknown"] is None
    ), f"'unknown' parameter should be None"

    for key in canopy_parameters:
        if key == "unknown":
            assert (
                canopy_subset[key] is None
            ), f"'{key}' parameter should be None"
        else:
            assert (
                canopy_subset[key] is not None
            ), f"'{key}' parameter should not be None"

    print(f"success parsing...'{parameter_file}'")


def test_parameter_access(domain, canopy_parameters):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    p_parameters = PrmsParameters(parameter_file)

    assert (
        p_parameters.parameters.srain_intcp is not None
    ), "'srain_intcp' should not return None"

    assert p_parameters.parameters.unknown is None, "'unknown' should return None"

    del p_parameters.parameters.srain_intcp
    assert p_parameters.parameters.srain_intcp is None, "'srain_intcp' should return None"
