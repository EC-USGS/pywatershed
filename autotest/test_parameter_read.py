import os
import sys

import pytest

from pynhm.utils import PrmsParameters
from utils import assert_or_print


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
