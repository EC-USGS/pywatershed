import pathlib as pl

import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.hydrology.prms_soilzone import PRMSSoilzone
from pywatershed.parameters import Parameters, PrmsParameters
from utils_compare import compare_in_memory, compare_netcdfs

# compare in memory (faster) or full output files? or both!
do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 5.0e-6

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(simulation):
    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    control.options["netcdf_output_var_names"] = PRMSSoilzone.get_variables()

    return control


@pytest.fixture(scope="function")
def discretization(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    return Parameters.from_netcdf(dis_hru_file, encoding=False)


@pytest.fixture(scope="function", params=params)
def parameters(simulation, control, request):
    if request.param == "params_one":
        param_file = simulation["dir"] / control.options["parameter_file"]
        params = PrmsParameters.load(param_file)
    else:
        param_file = simulation["dir"] / "parameters_PRMSSoilzone.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)

    comparison_var_names = list(
        set(PRMSSoilzone.get_variables())
        # These are not prms variables per se.
        # The _hru ones have non-hru equivalents being checked.
        # soil_zone_max and soil_lower_max would be nice to check but
        # prms5.2.1 wont write them as hru variables.
        - set(
            [
                "perv_actet_hru",
                "soil_lower_change_hru",
                "soil_lower_max",
                "soil_rechr_change_hru",
                "soil_zone_max",  # not a prms variable?
            ]
        )
    )
    # TODO: this is hacky, improve the design
    if not control.options["dprst_flag"]:
        comparison_var_names = {
            vv for vv in comparison_var_names if "dprst" not in vv
        }

    output_dir = simulation["output_dir"]

    input_variables = {}
    for key in PRMSSoilzone.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        # TODO: this is hacky for accommodating dprst_flag, improve the design
        # so people dont have to pass None for dead options.
        if not nc_path.exists():
            nc_path = None
        input_variables[key] = nc_path

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"]
        control.options["netcdf_output_dir"] = nc_parent

    soil = PRMSSoilzone(
        control=control,
        discretization=discretization,
        parameters=parameters,
        **input_variables,
        budget_type="error",
        calc_method=calc_method,
    )

    if do_compare_output_files:
        soil.initialize_netcdf()

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        soil.advance()
        soil.calculate(1.0)
        soil.output()
        if do_compare_in_memory:
            compare_in_memory(
                soil, answers, atol=atol, rtol=rtol, skip_missing_ans=True
            )

    soil.finalize()

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path / simulation["name"],
            output_dir,
            atol=atol,
            rtol=rtol,
            # fail_after_all_vars=False,
            verbose=True,
        )

    return
