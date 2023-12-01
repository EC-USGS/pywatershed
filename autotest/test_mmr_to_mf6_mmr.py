import flopy
import numpy as np
import pint
import pytest
import xarray as xr

from pywatershed.utils.mmr_to_mf6_mmr import MMRToMF6MMR

start_time = np.datetime64("1979-01-01T00:00:00")
end_time = np.datetime64("1979-01-07T00:00:00")


def lateral_flow_ans_ds(domain):
    control_file = domain["control_file"]
    inflow_dir = control_file.parent / "output"
    ds = xr.open_dataset(inflow_dir / "seg_lateral_inflow.nc")
    ds = ds.sel(time=slice(start_time, end_time))["seg_lateral_inflow"]
    return ds


# When we can point at modflow6 develop we'll un-xfail this
@pytest.mark.xfail
@pytest.mark.parametrize("bc_binary_files", [True, False])
@pytest.mark.parametrize("bc_flows_combine", [True, False])
def test_mmr_to_mf6(domain, tmp_path, bc_binary_files, bc_flows_combine):
    units = pint.UnitRegistry()
    ans = lateral_flow_ans_ds(domain)
    times = ans.time.values

    param_file = domain["param_file"]
    control_file = domain["control_file"]
    domain_name = domain["domain_name"]

    _ = MMRToMF6MMR(
        param_file=param_file,
        control_file=control_file,
        output_dir=tmp_path,
        inflow_dir=control_file.parent / "output",
        sim_name=domain_name,
        start_time=start_time,
        end_time=end_time,
        bc_binary_files=bc_binary_files,
        bc_flows_combine=bc_flows_combine,
    )

    if bc_flows_combine:
        flow_names = ["combined"]
    else:
        flow_names = ["sroff", "ssres_flow", "gwres_flow"]

    if bc_binary_files:
        for tt in times:
            flow_dict = {}
            for flow_name in flow_names:
                np_file = (
                    tmp_path / "snf_flw_bc"
                ) / f"flw_{flow_name}_{str(tt)[0:10]}T00_00_00.bin"
                assert np_file.exists()
                flow_dict[flow_name] = np.fromfile(
                    np_file, dtype=[("irch", "<i4"), ("q", "<f8")]
                )

            if bc_flows_combine:
                result = flow_dict["combined"]["q"]
            else:
                result = sum([ra["q"] for ra in flow_dict.values()])

            ans_tt = (
                (ans.sel(time=tt).values * units("feet ** 3 / s"))
                .to("meter ** 3 / s")
                .magnitude
            )

    else:
        sim = flopy.mf6.MFSimulation.load(
            "sim", "mfsim.nam", sim_ws=str(tmp_path)
        )
        model = sim.get_model(domain_name)
        flw_spds = {}
        for flow_name in flow_names:
            flw_spds[flow_name] = model.get_package(
                flow_name
            ).stress_period_data.get_data()  # all the way DOWN

        # can we check the times?
        for ii, tt in enumerate(times):
            result = sum([flw[ii]["q"] for flw in flw_spds.values()])
            ans_tt = (
                (ans.sel(time=tt).values * units("feet ** 3 / s"))
                .to("meter ** 3 / s")
                .magnitude
            )

    # <<
    # Compare
    abs_diff = abs(result - ans_tt)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = abs_diff / ans_tt

    abs_tol = 1.0e-5
    rel_tol = 1.0e-5

    abs_close = abs_diff < abs_tol
    rel_close = rel_diff < rel_tol
    rel_close = np.where(np.isnan(rel_close), False, rel_close)

    close = abs_close | rel_close

    assert close.all()

    return
