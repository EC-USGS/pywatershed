from pprint import pprint

from numpy.testing import assert_equal

from pywatershed.base.parameters import _set_dict_read_write
from pywatershed.parameters import Parameters, PrmsParameters
from utils import assert_dicts_equal


def test_param_dd_param(domain):
    # round trip from read-only to read-write to read-only
    # use a PRMS Parameter file for now
    domain_dir = domain["param_file"].parent
    params = Parameters.from_netcdf(
        domain_dir / "parameters_PRMSGroundwater.nc"
    )
    # These both copy by default
    param_dd = params.to_dd()
    params2 = Parameters(**param_dd.data)

    assert_dicts_equal(params.data, params2.data)

    param_dd.data_vars["gwflow_coef"] *= 4
    params3 = Parameters(**param_dd.data)

    try:
        assert_dicts_equal(params.data, params3.data)
        assert False
    except AssertionError:
        pass

    return
