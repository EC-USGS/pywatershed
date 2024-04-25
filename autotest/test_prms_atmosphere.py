from warnings import warn

import numpy as np
import pytest
from utils_compare import compare_in_memory, compare_netcdfs

from pywatershed.atmosphere.prms_atmosphere import PRMSAtmosphere
from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.parameters import PrmsParameters

# pptmix is altered by PRMSCanopy
# the answers on those times,locations are set to the true answers below
# only for compare_in_memory. comparing output files will fail if
# pptmix is changed by PRMSCanopy

# compare in memory (faster) or full output files? or both!
# do_compare_output_files=True results in failure on pptmix in certain domains
do_compare_output_files = False  # see above
do_compare_in_memory = True
rtol = 1.0e-5
atol = 1.0e-5  # why is this relatively low accuracy?

params = ["params_sep", "params_one"]


@pytest.fixture(scope="function")
def control(simulation):
    ctl = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    del ctl.options["netcdf_output_dir"]
    return ctl


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
        param_file = simulation["dir"] / "parameters_PRMSAtmosphere.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.domain
def test_compare_prms(
    simulation, control, discretization, parameters, tmp_path
):
    comparison_var_names = PRMSAtmosphere.get_variables()

    output_dir = simulation["output_dir"]
    cbh_dir = simulation["dir"]
    # cbh_dir = simulation["cbh_inputs"]["prcp"].parent.resolve()

    input_variables = {}
    for key in PRMSAtmosphere.get_inputs():
        if "soltab" in key:
            nc_path = simulation["output_dir"] / f"{key}.nc"
        else:
            nc_path = cbh_dir / f"{key}.nc"
        # <
        print(nc_path, nc_path.exists())
        input_variables[key] = nc_path

    atm = PRMSAtmosphere(
        control=control,
        discretization=discretization,
        parameters=parameters,
        **input_variables,
    )

    # check the advance/calculate the state
    tmaxf_id = id(atm.tmaxf)

    atm.initialize_netcdf(output_dir=tmp_path)

    # this is the condition which may modify pptmix, see below
    new_snow = adapter_factory(
        output_dir / "newsnow.nc",
        variable_name="newsnow",
        control=control,
    )

    if do_compare_in_memory:
        answers = {}
        for var in comparison_var_names:
            var_pth = output_dir / f"{var}.nc"
            answers[var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for ii in range(control.n_times):
        control.advance()
        atm.advance()
        if ii == 0:
            atm.output()
        atm.calculate(1.0)

        if do_compare_in_memory:
            # this is because pptmix is altered by canopy using this
            # net_snow condition from canopy
            new_snow.advance()
            # answers are behind, a second advance in compare_in_memory will
            # keep it at the same time, even with control
            answers["pptmix"].advance()
            # this is the condition that zeros pptmix
            cond = new_snow.current == 0
            wh_cond = np.where(cond)
            if (
                cond.any()
                and (
                    answers["pptmix"].current[wh_cond]
                    != atm.pptmix.data[ii, wh_cond]
                ).any()
            ):
                msg = (
                    "Editing answers for pptmix to match results where "
                    "PRMSCanopy edits the value later, this pywatershed "
                    "functionality is properly tested by "
                    "test_prms_above_snow."
                )
                warn(msg)
                atm.pptmix._current = np.where(
                    cond, answers["pptmix"].current, atm.pptmix.data[ii, :]
                )

            compare_in_memory(
                atm,
                answers,
                atol=atol,
                rtol=rtol,
                fail_after_all_vars=False,
                verbose=False,
            )

            assert id(atm.tmaxf) == tmaxf_id

    if do_compare_output_files:
        compare_netcdfs(
            comparison_var_names,
            tmp_path,
            output_dir,
            atol=atol,
            rtol=rtol,
            print_var_max_errs=False,
        )

    return
