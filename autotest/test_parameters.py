import pathlib as pl

import numpy as np
import pytest
from utils import assert_or_print

from pywatershed import Parameters, PRMSCanopy
from pywatershed.parameters import PrmsParameters


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
        key: val for key, val in param_obj.parameters.items() if key in answers.keys()
    }
    assert_or_print(results, answers)

    print("success initializing PrmsParameters object")


def test_parameter_read(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    # check dimensions
    answers = domain["test_ans"]["parameter_read"]
    results = parameters.dims
    assert_or_print(results, answers, print_ans=domain["print_ans"])

    print(f"success parsing...'{parameter_file}'")
    return


def test_parameter_canopy_subset(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")
    parameters = PrmsParameters.load(parameter_file)

    canopy_params = PRMSCanopy.get_parameters()
    canopy_subset = parameters.subset(canopy_params)

    canopy_params_2 = tuple([*canopy_params, "bar"])
    with pytest.raises(KeyError):
        canopy_subset = parameters.subset(canopy_params_2)

    for key in canopy_params:
        if key != "unknown":
            assert (
                canopy_subset.parameters[key] is not None
            ), f"'{key}' parameter should not be None"

    print(f"success parsing...'{parameter_file}'")


def test_parameter_access(domain):
    parameter_file = domain["param_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    assert (
        parameters.parameters["srain_intcp"] is not None
    ), "'srain_intcp' should not return None"

    with pytest.raises(KeyError):
        _ = parameters.parameters["unknown"]

    # can not delete parameters this way so will not raise a keyerror
    del parameters.parameters["srain_intcp"]
    _ = parameters.parameters["srain_intcp"]


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
    assert parameters.data.keys() == params_from_json.data.keys()
    for kk in parameters.data.keys():
        np.testing.assert_equal(
            dict(parameters.data[kk]), dict(params_from_json.data[kk])
        )
    return
