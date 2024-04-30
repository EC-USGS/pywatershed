import numpy as np
import pytest
import xarray as xr
from utils_compare import compare_in_memory

from pywatershed import PRMSChannel
from pywatershed.base.adapter import AdapterNetcdf, adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowGraph
from pywatershed.base.model import Model
from pywatershed.base.parameters import Parameters
from pywatershed.constants import zero
from pywatershed.hydrology.prms_channel_flow_graph import (
    HruSegmentFlowAdapter,
    HruSegmentFlowExchange,
    PRMSChannelFlowNodeMaker,
)
from pywatershed.parameters import PrmsParameters

# NB: THere is no real comparison of output files because the answer files
#     have different units. Could create a class to manage this but
#     I've checked that the memory values match the file values.
#     do_compare_output_files write files but does not check them
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7

calc_methods = ("numpy", "numba")

rename_vars = {
    "channel_outflow_vol": "outflows",
    "seg_upstream_inflow": "node_upstream_inflows",
    "seg_outflow": "node_outflows",
    "seg_stor_change": "node_storage_changes",
}

# the above values in rename_vars need converted to volume in these cases
convert_to_vol = ["outflows", "node_storage_changes"]


@pytest.fixture(scope="function")
def control(simulation):
    if "obsin" in simulation["name"]:
        pytest.skip(
            "Not testing prms_channel_flow_graph for drb_2yr:nhm_obsin"
        )
    if "hru_1" in simulation["name"]:
        pytest.skip("Not testing prms_channel_flow_graph for hru_1")

    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    if control.options["streamflow_module"] == "strmflow":
        pytest.skip(
            f"PRMSChannel not present in simulation {simulation['name']}"
        )

    del control.options["netcdf_output_dir"]
    del control.options["netcdf_output_var_names"]
    return control


@pytest.fixture(scope="function")
def discretization(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = simulation["dir"] / "parameters_dis_seg.nc"
    return Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )


@pytest.fixture(scope="function")
def parameters(simulation, control):
    param_file = simulation["dir"] / "parameters_PRMSChannel.nc"
    return PrmsParameters.from_netcdf(param_file)


@pytest.mark.parametrize("calc_method", calc_methods)
def test_prms_channel_flow_graph_compare_prms(
    simulation, control, discretization, parameters, tmp_path, calc_method
):
    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization, parameters, calc_method=calc_method
        ),
    }

    # combine lateral inflows to a single non-volumetric inflow
    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = AdapterNetcdf(nc_path, key, control)

    inflow_prms = HruSegmentFlowAdapter(parameters, **input_variables)

    # FlowGraph
    nnodes = parameters.dims["nsegment"]
    params_flow_graph = Parameters(
        dims={
            "nnodes": nnodes,
            # "nnode_types": 1,
        },
        coords={
            "node_coord": np.arange(nnodes),
        },
        data_vars={
            "node_maker_name": ["prms_channel"] * nnodes,
            "node_maker_index": np.arange(nnodes),
            "to_graph_index": discretization.parameters["tosegment"] - 1,
        },
        metadata={
            "node_coord": {"dims": ["nnodes"]},
            "node_maker_name": {"dims": ["nnodes"]},
            "node_maker_index": {"dims": ["nnodes"]},
            "to_graph_index": {"dims": ["nnodes"]},
        },
        validate=True,
    )

    flow_graph = FlowGraph(
        control,
        discretization=None,
        parameters=params_flow_graph,
        inflows=inflow_prms,
        node_maker_dict=node_maker_dict,
        budget_type="error",
    )

    if do_compare_output_files:
        nc_parent = tmp_path / simulation["name"].replace(":", "_")
        flow_graph.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            if var not in rename_vars.keys():
                continue
            new_var = rename_vars[var]
            var_pth = output_dir / f"{var}.nc"
            answers[new_var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

        # answers for the exchange
        var = "seg_lateral_inflow"
        var_pth = output_dir / f"{var}.nc"
        lateral_inflow_answers = adapter_factory(
            var_pth, variable_name=var, control=control
        )

    for istep in range(control.n_times):
        control.advance()
        flow_graph.advance()
        flow_graph.calculate(1.0)
        flow_graph.output()

        # check exchange
        lateral_inflow_answers.advance()
        np.testing.assert_allclose(
            inflow_prms.current,
            lateral_inflow_answers.current,
            rtol=rtol,
            atol=atol,
        )

        for var in answers.values():
            var.advance()

        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers.items():
                if key in convert_to_vol:
                    answers_conv_vol[key] = val.current / (24 * 60 * 60)
                else:
                    answers_conv_vol[key] = val.current

            # <<
            # there are no expected sources or sinks in this test
            answers_conv_vol["sink_source"] = val.current * zero

            compare_in_memory(
                flow_graph,
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=False,
                fail_after_all_vars=False,
            )

    flow_graph.finalize()


def test_hru_segment_flow_exchange(
    simulation,
    control,
    discretization,
    parameters,
    tmp_path,
):
    control.options["netcdf_output_var_names"] = ["outflow"]
    # run_dir = tmp_path / "test_hru_segment_exchange"

    control.options = control.options | {
        "input_dir": simulation["output_dir"],
        "budget_type": "error",
        "calc_method": "numba",
        "netcdf_output_dir": None,
    }

    # FlowGraph arguments
    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization, parameters, calc_method="numpy"
        ),
    }

    nnodes = parameters.dims["nsegment"]
    params_flow_graph = Parameters(
        dims={
            "nnodes": nnodes,
        },
        coords={
            "node_coord": np.arange(nnodes),
        },
        data_vars={
            "node_maker_name": ["prms_channel"] * nnodes,
            "node_maker_index": np.arange(nnodes),
            "to_graph_index": discretization.parameters["tosegment"] - 1,
        },
        metadata={
            "node_coord": {"dims": ["nnodes"]},
            "node_maker_name": {"dims": ["nnodes"]},
            "node_maker_index": {"dims": ["nnodes"]},
            "to_graph_index": {"dims": ["nnodes"]},
        },
        validate=True,
    )

    # Exchange parameters
    # TODO: this is funky, can we make this more elegant?
    params_ds = parameters.to_xr_ds().copy()
    params_ds["node_coord"] = xr.Variable(
        dims="nnodes",
        data=np.arange(parameters.dims["nsegment"]),
    )
    params_ds = params_ds.set_coords("node_coord")
    params_exchange = Parameters.from_ds(params_ds)

    model_dict = {
        "control": control,
        "dis_both": discretization,
        "model_order": ["hru_seg_exchange", "prms_channel_graph"],
        "hru_seg_exchange": {
            "class": HruSegmentFlowExchange,
            "parameters": params_exchange,
            "dis": "dis_both",
        },
        "prms_channel_graph": {
            "class": FlowGraph,
            "node_maker_dict": node_maker_dict,
            "parameters": params_flow_graph,
            "dis": None,
            "budget_type": "error",
        },
    }

    # if do_compare_output_files:
    #     nc_parent = tmp_path / simulation["name"].replace(":", "_")
    #     flow_graph.initialize_netcdf(nc_parent)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            if var not in rename_vars.keys():
                continue
            new_var = rename_vars[var]
            var_pth = simulation["output_dir"] / f"{var}.nc"
            answers[new_var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

        # answers for the exchange
        var = "seg_lateral_inflow"
        var_pth = simulation["output_dir"] / f"{var}.nc"
        lateral_inflow_answers = adapter_factory(
            var_pth, variable_name=var, control=control
        )

    model = Model(model_dict)
    for istep in range(control.n_times):
        model.advance()
        model.calculate()
        model.output()

        # check exchange
        lateral_inflow_answers.advance()
        np.testing.assert_allclose(
            model.processes["hru_seg_exchange"]["inflows"],
            lateral_inflow_answers.current,
            rtol=rtol,
            atol=atol,
        )

        for var in answers.values():
            var.advance()

        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers.items():
                if key in convert_to_vol:
                    answers_conv_vol[key] = val.current / (24 * 60 * 60)
                else:
                    answers_conv_vol[key] = val.current

            # <<
            # there are no expected sources or sinks in this test
            answers_conv_vol["sink_source"] = val.current * zero

            compare_in_memory(
                model.processes["prms_channel_graph"],
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=False,
                fail_after_all_vars=False,
            )

    model.finalize()
