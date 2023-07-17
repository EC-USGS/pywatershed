import pathlib as pl

import numpy as np
import pytest

from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.constants import epsilon32, zero
from pywatershed.hydrology.prms_snow import PRMSSnow
from pywatershed.parameters import Parameters, PrmsParameters

calc_methods = ("numpy", "numba")
params = ("params_sep", "params_one")


@pytest.fixture(scope="function")
def control(domain):
    return Control.load(domain["control_file"])


@pytest.fixture(scope="function")
def discretization(domain):
    dis_hru_file = domain["dir"] / "parameters_dis_hru.nc"
    return Parameters.from_netcdf(dis_hru_file, encoding=False)


@pytest.fixture(scope="function", params=params)
def parameters(domain, request):
    if request.param == "params_one":
        params = PrmsParameters.load(domain["param_file"])
    else:
        param_file = domain["dir"] / "parameters_PRMSSnow.nc"
        params = PrmsParameters.from_netcdf(param_file)

    return params


@pytest.mark.xfail
@pytest.mark.parametrize("calc_method", calc_methods)
def test_compare_prms(
    domain, control, discretization, parameters, tmp_path, calc_method
):
    tmp_path = pl.Path(tmp_path)
    output_dir = domain["prms_output_dir"]

    # get the answer data
    comparison_var_names = [
        # # "ai",
        # "albedo",
        # # "frac_swe",
        # # "freeh2o",
        # # "iasw",
        # # "int_alb",
        "iso",
        # # "lso",
        # # "lst",
        # # "mso",
        # # "pk_def",
        # # "pk_den",
        # # "pk_depth",
        # # "pk_ice",
        # # "pk_precip",
        # # "pk_temp",
        # # "pksv",
        # "pkwater_ante", -- dosent respect the iso memory eval
        "pkwater_equiv",
        # "pptmix_nopack",
        # # "pss",
        # "pst",
        # # "salb",
        # # "scrv",
        # # "slst",
        "snow_evap",
        # "snowcov_area",  # issues on 12-23-1979 in drb
        # "snowcov_areasv",
        # "snowmelt",  # issues on 12-22-1979 in drb
        # #"snsv",
        "tcal",
    ]

    ans = {}
    for key in comparison_var_names:
        nc_pth = output_dir / f"{key}.nc"
        ans[key] = adapter_factory(nc_pth, variable_name=key, control=control)

    # setup the snow
    input_variables = {}
    for key in PRMSSnow.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = nc_path

    snow = PRMSSnow(
        control,
        discretization,
        parameters,
        **input_variables,
        budget_type="warn",
        calc_method=calc_method,
    )

    all_success = True
    iso_censor_mask = None
    for istep in range(control.n_times):
        # for istep in range(10):

        print("\n")
        print(f"solving time: {control.current_time}")

        control.advance()
        snow.advance()

        snow.calculate(float(istep))

        # compare along the way
        for key, val in ans.items():
            val.advance()

        # pkwater_equiv ---
        # first check pkwater_equiv and assert they are with
        # their own tolerance (approx float error)
        key = "pkwater_equiv"
        pkwe_ans = ans[key].current
        pkwe = snow[key]
        pkwater_absdiff = abs(pkwe - pkwe_ans)
        atol = 5e-2  # 5 hundreds of an inch

        close = np.isclose(pkwe, pkwe_ans, atol=atol, rtol=zero)
        if iso_censor_mask is None:
            assert close.all()
        else:
            assert (close[~iso_censor_mask]).all()

        # if one of them is "zero" and the other is not,
        # then there can be differences in other variables.
        # so we'll exclude this case from the comparison of
        # other variables as long as the SWEs pass the above test
        pkwater_equiv_one_zero = (pkwe <= atol) | (pkwe_ans <= atol)
        if iso_censor_mask is not None:
            pkwater_equiv_one_zero = np.where(
                iso_censor_mask, True, pkwater_equiv_one_zero
            )
        print(f"pkwater_equiv_one_zero.sum(): " f"{pkwater_equiv_one_zero.sum()}")

        # iso ---
        # iso has a memory of packwater equiv which can cause about
        # a 20% change in albedo based on pkwater differences in the
        # noise... AND those differences in albedo can translate into
        # longer term differences in other variables. The solution is
        # to flag (memory) when iso gets out of snyc over reasonable/small
        # pkwater equiv differences no not evaluate those points and
        # after pk water has gone to zero for both pywatershed and prms.
        # we'll want to report the number of points being censored this
        # way at each time. this censor mask has to be used for pkwater
        # above
        key = "iso"
        iso_ans = ans[key].current
        iso = snow[key]
        iso_absdiff = abs(iso - iso_ans)

        # anytime at least one of them is zero AND there
        # are iso differences, then we add these points to the
        # iso_censor list
        iso_censor_mask_new = pkwater_equiv_one_zero & (iso_absdiff > 0)

        if iso_censor_mask is None:
            iso_censor_mask = iso_censor_mask_new
        else:
            iso_censor_mask = iso_censor_mask | iso_censor_mask_new

        # reset points where near zero and have matching iso
        pkwater_both_zero_and_iso = (
            (pkwe <= epsilon32) & (pkwe_ans <= epsilon32) & (iso == iso_ans)
        )
        iso_censor_mask = np.where(pkwater_both_zero_and_iso, False, iso_censor_mask)
        print(f"iso_censor_mask.sum(): {iso_censor_mask.sum()}")
        assert np.where(~iso_censor_mask, iso - iso_ans == 0, True).all()

        censor_comp = iso_censor_mask | pkwater_equiv_one_zero

        # The rest of the variables to be checked
        atol = 1e-5
        for key in ans.keys():
            if key in ["pkwater_equiv", "iso"]:
                continue

            print(key)
            a1 = ans[key].current
            a2 = snow[key]

            success = np.allclose(a1, a2, atol=atol, rtol=zero)
            atol = 1e-3
            success_prelim = np.isclose(a1, a2, atol=atol)
            success = np.where(~censor_comp, success_prelim, True).all()

            if not success:
                all_success = False
                diff = a1 - a2
                diffmin = diff.min()
                diffmax = diff.max()
                print(f"time step {istep}")
                print(f"output variable: {key}")
                print(f"prms   {a1.min()}    {a1.max()}")
                print(f"pywatershed  {a2.min()}    {a2.max()}")
                print(f"diff   {diffmin}  {diffmax}")

                zz = abs(a2 - a1)
                max_diff = np.where(~censor_comp, zz, zero).max()
                wh_max_diff = np.where(zz == max_diff)
                print(f"wh_max_diff: {wh_max_diff}")
                print(f"max diff: {zz[wh_max_diff]}")
                print(f"prms: {a1[wh_max_diff]}")
                print(f"pywatershed: {a2[wh_max_diff]}")
                print(
                    f"pkwater_absdiff[wh_max_diff]: " f"{pkwater_absdiff[wh_max_diff]}"
                )
                print(f"pkwe[wh_max_diff]: {pkwe[wh_max_diff]}")
                print(f"pkwe_ans[wh_max_diff]: {pkwe_ans[wh_max_diff]}")
                assert success

    snow.finalize()

    if not all_success:
        raise Exception("pywatershed results do not match prms results")

    return
