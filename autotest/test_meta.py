import pytest
from pywatershed.base import meta
from pywatershed.hydrology.prms_groundwater import PRMSGroundwater


@pytest.mark.domainless
def test_init():
    # some random checks
    assert "nhru" in meta.dimensions.keys()
    assert "transp_on_dynamic" in meta.control.keys()
    assert "covden_win" in meta.parameters.keys()
    assert "transp_on" in meta.variables.keys()
    return


@pytest.mark.domainless
def test_get_in_list():
    gw_vars = PRMSGroundwater.get_variables()
    gw_var_meta = meta.get_vars(gw_vars)
    assert set(gw_var_meta.keys()) == set(gw_vars)

    gw_inputs = PRMSGroundwater.get_inputs()
    gw_input_meta = meta.get_vars(gw_inputs)
    assert set(gw_input_meta.keys()) == set(gw_inputs)

    # There are dimensions in the PRMSGroundwater parameters
    # so these are not in the returned set. Not clear the best
    # way to go
    # gw_params = PRMSGroundwater.get_parameters()
    # gw_param_meta = meta.get_params(gw_params)
    # assert set(gw_param_meta.keys()) == set(gw_params)

    return
