import pathlib as pl

import numpy as np
import pytest

from pynhm import Parameters
from pynhm.parameters import PrmsParameters
from utils import assert_or_print


@pytest.fixture
def canopy_parameters():
    return tuple(("unknown", "srain_intcp", "wrain_intcp", "snow_intcp"))


def test_parameter_init():
    # TODO: this is now an invalid parameter object, fix?
    parameters = {
        "nhru": 2,
        "abc": np.zeros(2, dtype=float),
        "nmonths": 12,
        "nsegment": 10,
        "xyz": np.arange(12 * 10, dtype=float).reshape(12, 10),
    }
    param_obj = Parameters(data_vars=parameters, validate=False)

    answers = {"nhru": 2, "nmonths": 12, "nsegment": 10}

    results = {
        key: val
        for key, val in param_obj.parameters.items()
        if key in answers.keys()
    }
    assert_or_print(results, answers)

    print("success initializing PrmsParameters object")


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

    canopy_subset = parameters.subset(canopy_parameters)

    with pytest.raises(KeyError):
        v = canopy_subset.parameters["unknown"]

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
        parameters.parameters["srain_intcp"] is not None
    ), "'srain_intcp' should not return None"

    with pytest.raises(KeyError):
        v = parameters.parameters["unknown"]

    del parameters.parameters["srain_intcp"]
    with pytest.raises(KeyError):
        v = parameters.parameters["srain_intcp"]

    with pytest.raises(KeyError):
        del parameters.parameters["srain_intcp"]


def test_parameter_json(domain, tmp_path):
    # read a myparams.param file to a parameter object,
    # write it to json,
    # read the json back in
    # check that the read json gives an identical parameter object.
    parameter_file = domain["param_file"]
    parameters = PrmsParameters.load(parameter_file)

    json_file = pl.Path(tmp_path) / "params.json"
    parameters.parameters_to_json(json_file)
    assert json_file.exists()

    params_from_json = PrmsParameters.load_from_json(json_file)

    param_obj_keys = params_from_json.__dict__.keys()
    assert sorted(param_obj_keys) == sorted(parameters.__dict__.keys())

    for kk in param_obj_keys:
        vv_result = params_from_json.__dict__["parameter_dimensions"]
        vv_ans = parameters.__dict__["parameter_dimensions"]
        assert len(vv_result) == len(vv_ans)
        for kk1 in vv_result.keys():
            assert vv_result[kk1] == vv_ans[kk1]

    return
