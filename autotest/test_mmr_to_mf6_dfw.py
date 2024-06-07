import pathlib as pl

import flopy
import numpy as np
import pytest
import xarray as xr

import pywatershed as pws
from pywatershed.utils.mmr_to_mf6_dfw import MmrToMf6Dfw

# The point of this is to reproduce the modflow6/autotest/test_swf_dfw.py
# Needs mf6 in PATH on mac/linux


# See below to check if these answers are still up-to-date with
# mf6/autotest/test_swf_dfw.py
# Note the answers are from the binary FLW files in
# mf6/autotest/test_swf_dfw.py, if you switch to text files, the line should be
# flw_list = [
#     (int(binary), 100),
# ]  # one-based cell numbers here if binary, zero-based if text

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


@pytest.mark.parametrize("binary_flw", [True, False])
def test_mmr_to_mf6_dfw(tmp_path, binary_flw):
    # Here we supply "seg_mid_elevation" in the parameter data and
    # "stress_period_data" in chd options. This bypasses the calculation of
    # these given PRMS parameter information. In this case, it would be
    # impossible to supply PRMS-like data that would reproduce this test.

    # This example shows how to formulate pywatershed inputs to match MF6
    # inputs unfortunately, these have to be translated through some of the
    # assumptions of PRMS to go through MMRToMF6DFW
    # One is that there are flows coming from potentially multiple HRUs to
    # stream segments, the other is about the units of volume/flow being in
    # cubicfeet.

    name = "swf-dfw01"
    output_dir = tmp_path / "test_swf_dfw01"
    save_flows = True
    print_flows = True

    control = pws.Control(
        start_time=np.datetime64("2023-01-01T00:00:00"),
        end_time=np.datetime64("2023-01-01T00:00:00"),
        time_step=np.timedelta64(1, "s"),
    )

    # ims
    ims_options = {
        "print_option": "all",
        "linear_acceleration": "BICGSTAB",
        "outer_dvclose": 1.0e-7,
        "inner_dvclose": 1.0e-8,
    }

    # disv1d
    dx = 1000.0
    nreach = 3

    # mannings coefficients are actually considered part of the dfw package
    params = pws.parameters.PrmsParameters(
        dims={
            "nsegment": nreach,
            "nhru": nreach,
        },
        coords={
            "seg_id": np.array(range(nreach)),  # todo: rename seg_id or such
        },
        data_vars={
            # "tosegment" no longer used
            "seg_length": np.ones(nreach) * dx,
            "mann_n": np.ones(nreach) * 0.035,
            # "seg_slope": np.ones(nreach) * 0.001,
            "seg_width": np.ones(nreach) * 50.0,
            "seg_mid_elevation": np.zeros(nreach),
            "hru_segment": np.array(range(nreach)) + 1,
        },
        metadata={
            "seg_id": {"dims": ["nsegment"]},
            "seg_length": {"dims": ["nsegment"]},
            "mann_n": {"dims": ["nsegment"]},
            "seg_width": {"dims": ["nsegment"]},
            "seg_mid_elevation": {"dims": ["nsegment"]},
            "hru_segment": {"dims": ["nsegment"]},
        },
        validate=True,
    )

    # vertices could also be supplied by shapefiles
    vertices = []
    vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
    cell2d = []
    for j in range(nreach):
        cell2d.append([j, 0.5, 2, j, j + 1])
    nodes = len(cell2d)
    nvert = len(vertices)
    disv1d_options = {
        "nodes": nodes,
        "nvert": nvert,
        "vertices": vertices,
        "cell2d": cell2d,
    }

    # dfw
    dfw_options = {
        "print_flows": print_flows,
        "save_flows": save_flows,
        "idcxs": 0,  # zero based in flopy, None for hydraulically wide
    }

    # sto
    sto_options = {"save_flows": save_flows}

    # ic
    ic_options = {"strt": 1.0}

    # CXS
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

    # oc
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
            flw_vol = 3531.4666721489 * 86400  # 100.0 m3/s in ft^3/day
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
        tdis_perlen=1,  # does this matter
        tdis_nstp=1,
        output_dir=output_dir,
        bc_binary_files=binary_flw,
        sim_name=name,
        save_flows=save_flows,
        ims_options=ims_options,
        dfw_options=dfw_options,
        sto_options=sto_options,
        ic_options=ic_options,
        oc_options=oc_options,
        chd_options=chd_options,
        disv1d_options=disv1d_options,
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
    assert success

    # one can verify the answers match the current mf6 results for test_swf_dfw
    # by placing the path to its output here
    # output_dir = pl.Path(
    #     "../../modflow6/autotest/keepers/test_mf6model[0-swf-dfw01]0"
    # )

    # check binary grid file
    grb = flopy.mf6.utils.MfGrdFile(output_dir / f"{name}.disv1d.grb")
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
    assert ((answers_swf_dfw["stage"] - stage) < 1.0e-7).all()

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
