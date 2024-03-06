import pathlib as pl

import numpy as np
import pytest

from pywatershed import Control, Parameters, PRMSCanopy
from pywatershed.parameters import PrmsParameters
from utils import assert_or_print

test_ans = {
    "drb_2yr": {
        "ndepl": 5,
        "ndeplval": 55,
        "ndoy": 366,
        "ngw": 765,
        "nhru": 765,
        "nmonth": 12,
        "nobs": 168,
        "npoigages": 168,
        "nsegment": 456,
        "nssr": 765,
        "scalar": 1,
    },
    "hru_1": {
        "nhru": 1,
        "nsegment": 1,
        "nssr": 1,
        "ngw": 1,
        "nobs": 1,
        "ndeplval": 11,
        "ndepl": 1,
        "nmonth": 12,
        "scalar": 1,
        "ndoy": 366,
    },
    "ucb_2yr": {
        "nhru": 3851,
        "nsegment": 1942,
        "nssr": 3851,
        "ngw": 3851,
        "npoigages": 347,
        "nobs": 347,
        "ndeplval": 55,
        "ndepl": 5,
        "nmonth": 12,
        "scalar": 1,
        "ndoy": 366,
    },
}


def test_parameter_init():
    # TODO: this test should be moved to tests for the Parameters class
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

    print("success initializing Parameters object")


def test_parameter_read(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    parameter_file = simulation["dir"] / ctl.options["parameter_file"]
    print(f"parsing...'{parameter_file}'")

    parameters = PrmsParameters.load(parameter_file)

    # check dimensions
    # this only depends on domain and not simulation AFAICT
    test_domain = simulation["name"].split(":")[0]
    answers = test_ans[test_domain]
    results = parameters.dims
    assert_or_print(results, answers, print_ans=simulation["print_ans"])

    print(f"success parsing...'{parameter_file}'")
    return


def test_parameter_canopy_subset(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    parameter_file = simulation["dir"] / ctl.options["parameter_file"]
    print(f"parsing...'{parameter_file}'")
    parameters = PrmsParameters.load(parameter_file)

    canopy_params = PRMSCanopy.get_parameters()
    canopy_subset = parameters.subset(canopy_params)

    canopy_params_2 = tuple([*canopy_params, "bar"])
    with pytest.raises(KeyError):
        canopy_subset = parameters.subset(canopy_params_2, strict=True)

    for key in canopy_params:
        if key != "unknown":
            assert (
                canopy_subset.parameters[key] is not None
            ), f"'{key}' parameter should not be None"

    print(f"success parsing...'{parameter_file}'")


def test_parameter_access(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    parameter_file = simulation["dir"] / ctl.options["parameter_file"]
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


def test_parameter_json(simulation, tmp_path):
    # read a myparams.param file to a parameter object,
    # write it to json,
    # read the json back in
    # check that the read json gives an identical parameter object.
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    parameter_file = simulation["dir"] / ctl.options["parameter_file"]
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
