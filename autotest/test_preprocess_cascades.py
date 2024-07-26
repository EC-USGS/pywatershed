import pytest

from pywatershed.base.control import Control
from pywatershed.parameters import Parameters, PrmsParameters
from pywatershed.utils.preprocess_cascades import calc_hru_route_order


@pytest.fixture(scope="function")
def control(simulation):
    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    return control


@pytest.fixture(scope="function")
def parameters(simulation, control):
    param_file = simulation["dir"] / control.options["parameter_file"]
    params = PrmsParameters.load(param_file)
    return params


def test_preprocess_cascades(parameters):
    new_params = calc_hru_route_order(parameters)
    assert "hru_route_order" in new_params.variables
    assert isinstance(new_params, Parameters)
