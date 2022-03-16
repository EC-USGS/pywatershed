import pytest

from pynhm.utils import PrmsParameters
from utils import assert_or_print


@pytest.fixture
def canopy_parameters():
    return tuple(("unknown", "srain_intcp", "wrain_intcp", "snow_intcp"))


def test_parameter_read(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    # check dimensions
    answers = domain["test_ans"]["parameter_read"]
    results = {
        key: val
        for key, val in parameters.parameters.items()
        if key in answers.keys()
    }
    assert_or_print(results, answers, print_ans=domain["print_ans"])

    print(f"success parsing...'{parameter_file}'")

    return


def test_parameter_canopy_subset(domain, canopy_parameters):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    canopy_subset = parameters.get_parameters(canopy_parameters)

    with pytest.raises(AttributeError):
        v = canopy_subset.parameters.unknown

    for key in canopy_parameters:
        if key != "unknown":
            assert (
                canopy_subset.parameters[key] is not None
            ), f"'{key}' parameter should not be None"

    print(f"success parsing...'{parameter_file}'")


def test_parameter_access(domain, canopy_parameters):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    assert (
        parameters.parameters.srain_intcp is not None
    ), "'srain_intcp' should not return None"

    with pytest.raises(AttributeError):
        v = parameters.parameters.unknown

    del parameters.parameters.srain_intcp
    with pytest.raises(AttributeError):
        v = parameters.parameters.srain_intcp

    with pytest.raises(AttributeError):
        del parameters.parameters.srain_intcp
