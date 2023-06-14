import pathlib as pl

import numpy as np
import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.PRMSCanopy import PRMSCanopy, has_prmscanopy_f
from pywatershed.parameters import PrmsParameters

calc_methods = ("numpy", "numba", "fortran")


def test_simple():
    time_dict = {
        "start_time": np.datetime64("1979-01-03T00:00:00.00"),
        "end_time": np.datetime64("1979-01-04T00:00:00.00"),
        "time_step": np.timedelta64(1, "D"),
    }

    nhru = 2
    prms_params = {
        "dims": {"nhru": nhru},
        "data_vars": {
            "hru_area": np.array(nhru * [1.0]),
            "covden_sum": np.array(nhru * [0.5]),
            "covden_win": np.array(nhru * [0.5]),
            "srain_intcp": np.array(nhru * [1.0]),
            "wrain_intcp": np.array(nhru * [1.0]),
            "snow_intcp": np.array(nhru * [1.0]),
            "epan_coef": np.array(nhru * [1.0]),
            "potet_sublim": np.array(nhru * [1.0]),
            "cov_type": np.array(nhru * [1]),
        },
    }
    prms_params["metadata"] = {}
    for kk in prms_params["data_vars"].keys():
        prms_params["metadata"][kk] = {"dims": ("nhru",)}

    prms_params = PrmsParameters(**prms_params)

    control = Control(**time_dict, params=prms_params)

    input_variables = {}
    for key in PRMSCanopy.get_inputs():
        input_variables[key] = np.ones([nhru])

    # todo: this is testing instantiation, but not physics
    cnp = PRMSCanopy(
        control=control,
        **input_variables,
        budget_type="error",
    )
    control.advance()
    cnp.advance()
    cnp.calculate(time_length=1.0)

    return


@pytest.fixture(scope="function", params=["params_sep", "params_one"])
def params(domain, request):
    if request.param == "params_one":
        params = PrmsParameters.load(domain["param_file"])
    else:
        params = PrmsParameters.from_netcdf(
            domain["dir"] / "parameters_PRMSCanopy.nc"
        )

    return params


@pytest.fixture(scope="function")
def control(domain, params):
    return Control.load(domain["control_file"], params=params)


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(domain, control, tmp_path, calc_method):
    if not has_prmscanopy_f and calc_method == "fortran":
        pytest.skip(
            "PRMSCanopy fortran code not available, skipping its test."
        )

    tmp_path = pl.Path(tmp_path)
    output_dir = domain["prms_output_dir"]

    # get the answer data
    comparison_var_names = [
        "net_rain",
        "net_snow",
        "net_ppt",
        "intcp_stor",
        "intcp_evap",
        "hru_intcpstor",
        "hru_intcpevap",
        "intcp_changeover",
    ]
    ans = {}
    for key in comparison_var_names:
        nc_pth = output_dir / f"{key}.nc"
        ans[key] = adapter_factory(nc_pth, variable_name=key, control=control)

    # setup the canopy
    input_variables = {}
    for key in PRMSCanopy.get_inputs():
        nc_pth = output_dir / f"{key}.nc"
        input_variables[key] = nc_pth

    cnp = PRMSCanopy(
        control=control,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    all_success = True
    for istep in range(control.n_times):
        control.advance()
        cnp.advance()
        cnp.calculate(1.0)

        # compare along the way
        atol = 1.0e-5
        for key, val in ans.items():
            val.advance()
        for key in ans.keys():
            a1 = ans[key].current
            a2 = cnp[key]
            success = np.isclose(a2, a1, atol=atol).all()
            if not success:
                all_success = False
                diff = a1 - a2
                diffmin = diff.min()
                diffmax = diff.max()
                print(f"time step {istep}")
                print(f"output variable {key}")
                print(f"prms   {a1.min()}    {a1.max()}")
                print(f"pywatershed  {a2.min()}    {a2.max()}")
                print(f"diff   {diffmin}  {diffmax}")

    cnp.finalize()

    if not all_success:
        raise Exception("pywatershed results do not match prms results")

    return
