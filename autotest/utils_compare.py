import pathlib as pl

import numpy as np
import xarray as xr

import pywatershed as pws


def assert_allclose(
    actual: np.ndarray,
    desired: np.ndarray,
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    print_max_errs: bool = False,
    var_name: str = "",
):
    """Reinvent np.testing.assert_allclose to get useful diagnostincs in debug

    Args:
        actual: Array obtained.
        desired: Array desired.
        rtol: Relative tolerance.
        atol: Absolute tolerance.
        equal_nan: If True, NaNs will compare equal.
        strict: If True, raise an ``AssertionError`` when either the shape or
            the data type of the arguments does not match. The special
            handling of scalars mentioned in the Notes section is disabled.
        also_check_w_np: first check using np.testing.assert_allclose using
            the same options.
        print_max_errs: bool=False. Print max abs and rel err for each var.
    """

    if also_check_w_np:
        np.testing.assert_allclose(
            actual,
            desired,
            rtol=rtol,
            atol=atol,
            equal_nan=equal_nan,
            # strict=strict,  # to add to newer versions of numpy
        )

    if strict:
        assert actual.shape == desired.shape
        assert isinstance(actual, type(desired))
        assert isinstance(desired, type(actual))

    if equal_nan:
        actual_nan = np.where(np.isnan(actual), True, False)
        desired_nan = np.where(np.isnan(desired), True, False)
        assert (actual_nan == desired_nan).all()

    abs_diff = abs(actual - desired)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_abs_diff = abs_diff / desired

    abs_close = abs_diff < atol
    rel_close = rel_abs_diff < rtol
    rel_close = np.where(np.isnan(rel_close), False, rel_close)

    close = abs_close | rel_close

    if print_max_errs:
        sp = "" if len(var_name) == 0 else " "
        print(f"{var_name}{sp}max abs err: {abs_diff.max()}")
        print(f"{var_name}{sp}max rel err: {rel_abs_diff.max()}")

    assert close.all()


def compare_in_memory(
    process: pws.base.Process,
    answers: dict[pws.base.adapter.AdapterNetcdf],
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    skip_missing_ans: bool = False,
    fail_after_all_vars: bool = True,
    verbose: bool = False,
):
    # TODO: docstring

    fail_list = []
    for var in process.get_variables():
        if var not in answers.keys():
            if skip_missing_ans:
                continue
            else:
                msg = f"Variable '{var}' not found in the answers provided."
                raise KeyError(msg)

        if verbose:
            print(f"checking {var}")
        answers[var].advance()

        if isinstance(process[var], pws.base.timeseries.TimeseriesArray):
            actual = process[var].current
        else:
            actual = process[var]

        if isinstance(answers[var], pws.base.adapter.Adapter):
            desired = answers[var].current.data
        else:
            desired = answers[var]

        if not fail_after_all_vars:
            assert_allclose(
                actual,
                desired,
                atol=atol,
                rtol=rtol,
                equal_nan=equal_nan,
                strict=strict,
                also_check_w_np=also_check_w_np,
                var_name=var,
            )

        else:
            try:
                assert_allclose(
                    actual,
                    desired,
                    atol=atol,
                    rtol=rtol,
                    equal_nan=equal_nan,
                    strict=strict,
                    also_check_w_np=also_check_w_np,
                    var_name=var,
                )
                if verbose:
                    print(f"compare_netcdfs all close for variable: {var}")

            except AssertionError:
                fail_list += [var]
                if verbose:
                    print(f"compare_netcdfs NOT all close for variable: {var}")

    if len(fail_list) > 0:
        assert False, f"compare_netcdfs failed for variables: {fail_list}"


def compare_netcdfs(
    var_list: list,
    results_dir: pl.Path,
    answers_dir: pl.Path,
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    print_var_max_errs: bool = False,
    fail_after_all_vars: bool = True,
    verbose: bool = False,
):
    # TODO: docstring
    # TODO: improve error message
    # TODO: collect failures in a try and report at end
    fail_list = []
    for var in var_list:
        answer = xr.open_dataarray(
            answers_dir / f"{var}.nc", decode_timedelta=False
        )
        result = xr.open_dataarray(
            results_dir / f"{var}.nc", decode_timedelta=False
        )

        if not fail_after_all_vars:
            assert_allclose(
                actual=result.values,
                desired=answer.values,
                rtol=rtol,
                atol=atol,
                equal_nan=equal_nan,
                strict=strict,
                also_check_w_np=also_check_w_np,
                print_max_errs=print_var_max_errs,
                var_name=var,
            )
            if verbose:
                print(f"compare_netcdfs all close for variable: {var}")

        else:
            try:
                assert_allclose(
                    actual=result.values,
                    desired=answer.values,
                    rtol=rtol,
                    atol=atol,
                    equal_nan=equal_nan,
                    strict=strict,
                    also_check_w_np=also_check_w_np,
                    print_max_errs=print_var_max_errs,
                    var_name=var,
                )
                if verbose:
                    print(f"compare_netcdfs all close for variable: {var}")

            except AssertionError:
                fail_list += [var]
                if verbose:
                    print(f"compare_netcdfs NOT all close for variable: {var}")

    if len(fail_list) > 0:
        assert False, f"compare_netcdfs failed for variables: {fail_list}"
