import pathlib as pl

import flopy
import numpy as np
import pytest
import xarray as xr

import pywatershed as pws
from pywatershed.utils.mmr_to_mf6_dfw import MmrToMf6Dfw

# The point of this is to reproduce the swf_dfw test run in mf6
# the answers of this test on 12/4/2023 are below.

# THis test calls mf6 via flopy and wont be tested in CI
# until swf merges to develop.

# needs mf6 in PATH on mac/linux

# Using modflow6 remote:branch langevin-usgs:feat-snf
# mf6/autotest: python update_flopy.py
# compile mf6 using gfortran and export LDFLAGS="$LDFLAGS -Wl,-ld_classic"

# This DOES check the output of mf6_mmr and needs a binary from a specialized
# mf6 branch. When feat-snf merges to develop we can un-xfail.

# Here we supply "seg_mid_elevation" in the parameter data and
# "stress_period_data" in chd options. This bypasses the calculation of these
# given PRMS parameter information. In this case, it would be impossible
# to supply PRMS-like data that would reproduce this test.

answers_swf_dfw = {
    "ia": np.array([0, 2, 5, 7]),
    "ja": np.array([0, 1, 1, 0, 2, 2, 1], dtype=np.int32),
    "stage": np.array([[[[1.00196123, 1.00003366, 1.0]]]]),
    "flowja": np.array(
        [
            [
                [
                    5.49604806e-12,
                    -1.93841991e00,
                    2.45701848e-10,
                    1.93841991e00,
                    -2.55350330e-01,
                    0.00000000e00,
                    2.55350330e-01,
                ]
            ]
        ]
    ),
    "qstorage": np.array([[[-98.06158009, -1.68306958, 0.0]]]),
    "qflw": np.array([1.0, 1.0, 100.0]),
    "qchd": np.array([3.0, 1.0, -0.25535033]),
    "qresidual": 0.0,
}


@pytest.mark.fail
@pytest.mark.parametrize("binary_flw", [True, False])
def test_mmr_to_mf6_dfw(tmp_path, binary_flw):
    name = "swf-dfw01"
    output_dir = tmp_path / "test_swf_dfw"
    save_flows = True

    # This example shows how to formulate pywatershed inputs to match MF6
    # inputs unfortunately, these have to be translated through some of the
    # assumptions of PRMS to go through MMRToMF6DFW
    # One is that there are flows coming from potentially multiple HRUs to
    # stream segments, the other is about the units of volume/flow being in
    # cubicfeet.

    # create a dummy control with
    control = pws.Control(
        start_time=np.datetime64("2023-01-01T00:00:00"),
        end_time=np.datetime64("2023-01-01T00:00:00"),
        time_step=np.timedelta64(1, "s"),
    )

    # create a dummy domain with
    nreach = 3
    params = pws.parameters.PrmsParameters(
        dims={
            "nsegment": nreach,
            "nhru": nreach,
        },
        coords={
            "seg_id": np.array(range(nreach)),  # todo: rename seg_id or such
        },
        data_vars={
            "tosegment": np.array([2, 3, 0]),  # one-based index, 0 is outflow
            "seg_length": np.ones(nreach) * 1.0e3,
            "mann_n": np.ones(nreach) * 0.035,
            "seg_slope": np.ones(nreach) * 0.001,
            "seg_width": np.ones(nreach) * 50.0,
            "seg_mid_elevation": np.zeros(nreach),
            "hru_segment": np.array(range(nreach)) + 1,
        },
        metadata={
            "seg_id": {"dims": ["nsegment"]},
            "tosegment": {"dims": ["nsegment"]},
            "seg_length": {"dims": ["nsegment"]},
            "mann_n": {"dims": ["nsegment"]},
            "seg_slope": {"dims": ["nsegment"]},
            "seg_width": {"dims": ["nsegment"]},
            "seg_mid_elevation": {"dims": ["nsegment"]},
            "hru_segment": {"dims": ["nsegment"]},
        },
        validate=True,
    )

    xfraction = [0.0, 0.0, 1.0, 1.0]
    height = [100.0, 0.0, 0.0, 100.0]
    mannfraction = [1.0, 1.0, 1.0, 1.0]
    cxsdata = list(zip(xfraction, height, mannfraction))
    cxs_options = {
        "nsections": 1,
        "npoints": 4,
        "packagedata": [(0, 4)],
        "crosssectiondata": cxsdata,
    }

    ims_options = {
        "print_option": "all",
        "linear_acceleration": "BICGSTAB",
        "outer_dvclose": 1.0e-7,
        "inner_dvclose": 1.0e-8,
    }

    dfw_options = {
        "print_flows": True,
        "save_flows": save_flows,
        "idcxs": 0,  # zero based in flopy, None for hydraulically wide
    }

    sto_options = {"save_flows": save_flows}

    oc_options = {
        "saverecord": [
            ("STAGE", "ALL"),
            ("BUDGET", "ALL"),
            # ("QOUTFLOW", "ALL"),
        ],
        "printrecord": [
            ("STAGE", "LAST"),
            ("BUDGET", "ALL"),
            # ("QOUTFLOW", "ALL"),
        ],
    }

    ic_options = {"strt": 1.0}

    # boundary conditions / FLW
    # if we use the volumes from pws rather than the inches from
    # PRMS, we dont have to supply 'hru_area' in the parameters for
    # depth to volume conversion.
    inflow_from_PRMS = False
    # this will aggregate to a single flw variable when mf6 writes to file
    bc_flows_combined = True
    # we have to create files with the boundary conditions in them
    inflow_list = ["sroff_vol", "ssres_flow_vol", "gwres_flow_vol"]
    inflow_dir = output_dir / "inflow_dir"
    inflow_dir.mkdir(parents=True)
    for inflw in inflow_list:
        if inflw == "sroff_vol":
            flw_vol = 3531.4666721489  # 100.0 m3/s in ft^3/s
        else:
            flw_vol = 0.0

        _ = xr.Dataset(
            coords=dict(
                time=np.array([control.start_time]),
                nsegment=params.parameters["seg_id"],
            ),
            data_vars={
                f"{inflw}": (
                    ["time", "nsegment"],
                    np.array([[flw_vol, 0.0, 0.0]]),
                ),
            },
        ).to_netcdf(inflow_dir / f"{inflw}.nc")

    flw_options = {
        "print_input": True,
        "print_flows": True,
    }

    chd_options = {
        "maxbound": 1,
        "print_input": True,
        "print_flows": True,
        "stress_period_data": [(2, 1.0)],
    }

    dfw = MmrToMf6Dfw(
        control=control,
        params=params,
        output_dir=output_dir,
        sim_name=name,
        save_flows=save_flows,
        ims_options=ims_options,
        dfw_options=dfw_options,
        sto_options=sto_options,
        ic_options=ic_options,
        oc_options=oc_options,
        chd_options=chd_options,
        cxs_options=cxs_options,
        inflow_from_PRMS=inflow_from_PRMS,
        bc_flows_combined=bc_flows_combined,
        inflow_dir=inflow_dir,
        flw_options=flw_options,
    )

    dfw.write()

    success, buff = dfw.run(silent=False, report=True)

    # ======================
    # Checks

    # check binary grid file
    grb = flopy.mf6.utils.MfGrdFile(output_dir / f"{name}.disl.grb")
    ia = grb.ia
    ja = grb.ja
    assert (answers_swf_dfw["ia"] == ia).all()
    assert (answers_swf_dfw["ja"] == ja).all()
    assert ia.shape[0] == grb.nodes + 1, "ia in grb file is not correct size"

    # check stage file
    qobj = flopy.utils.HeadFile(
        output_dir / f"{name}.stage", precision="double", text="STAGE"
    )
    stage = qobj.get_alldata()
    assert (answers_swf_dfw["stage"] - stage < 1.0e-7).all()

    # read the budget file
    budobj = flopy.utils.binaryfile.CellBudgetFile(output_dir / f"{name}.bud")
    flowja = budobj.get_data(text="FLOW-JA-FACE")[0]
    qstorage = budobj.get_data(text="STORAGE")[0]
    qflw = np.array(budobj.get_data(text="FLW")[0].tolist()[0])
    qchd = np.array(budobj.get_data(text="CHD")[0].tolist()[0])
    qresidual = np.zeros(grb.nodes)[0]

    assert (answers_swf_dfw["flowja"] - flowja < 1.0e-7).all()
    assert (answers_swf_dfw["qstorage"] - qstorage < 1.0e-7).all()
    assert (answers_swf_dfw["qflw"] - qflw < 1.0e-7).all()
    assert (answers_swf_dfw["qchd"] - qchd < 1.0e-7).all()
    assert (answers_swf_dfw["qresidual"] - qresidual < 1.0e-7).all()
