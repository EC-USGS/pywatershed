import pathlib as pl

import numpy as np
import pywatershed as pws
import xarray as xr


def assert_allclose(
    actual: np.ndarray,
    desired: np.ndarray,
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    error_message: str = "Comparison unsuccessful (default message)",
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
    rel_abs_diff = abs_diff / desired

    abs_close = abs_diff < atol
    rel_close = rel_abs_diff < rtol
    rel_close = np.where(np.isnan(rel_close), False, rel_close)

    close = abs_close | rel_close

    assert close.all()


def compare_in_memory(
    process: pws.base.Process,
    answers: dict[pws.base.adapter.AdapterNetcdf],
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    error_message: str = None,
):
    # TODO: docstring
    for var in process.get_variables():
        answers[var].advance()
        assert_allclose(
            process[var],
            answers[var].current.data,
            atol=atol,
            rtol=rtol,
            equal_nan=equal_nan,
            strict=strict,
            also_check_w_np=also_check_w_np,
            error_message=error_message,
        )


def compare_netcdfs(
    var_list: list,
    results_dir: pl.Path,
    answers_dir: pl.Path,
    rtol: float = 1.0e-15,
    atol: float = 1.0e-15,
    equal_nan: bool = True,
    strict: bool = False,
    also_check_w_np: bool = True,
    error_message: str = None,
):
    # TODO: docstring
    # TODO: improve error message
    # TODO: collect failures in a try and report at end
    for var in var_list:
        answer = xr.open_dataarray(answers_dir / f"{var}.nc")
        result = xr.open_dataarray(results_dir / f"{var}.nc")

        if error_message is None:
            error_message = f"Comparison of variable '{var}' was unsuccessful"

        assert_allclose(
            actual=result.values,
            desired=answer.values,
            rtol=rtol,
            atol=atol,
            equal_nan=equal_nan,
            strict=strict,
            also_check_w_np=also_check_w_np,
            error_message=error_message,
        )
