import pathlib as pl

import numpy as np
import pytest
import xarray as xr
from utils_compare import compare_in_memory

from pywatershed import (
    Model,
    PRMSChannel,
    PRMSGroundwater,
    PRMSRunoff,
    PRMSSoilzone,
)
from pywatershed.base.adapter import adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.parameters import Parameters
from pywatershed.constants import nan, zero
from pywatershed.hydrology.pass_through_node import PassThroughNodeMaker
from pywatershed.hydrology.prms_channel_flow_graph import (
    prms_channel_flow_graph_postprocess,
    prms_channel_flow_graph_to_model_dict,
)
from pywatershed.hydrology.starfit import StarfitFlowNodeMaker
from pywatershed.parameters import PrmsParameters

# Purpose: Test that STARTFIT in a PRMSChannel FlowGraph gives correct results.
# Below the inserted reservoir, we dont have answers, but we'll check
# that the PRMSChannel results are not altered elsewhere on the graph.
# Tests below use "postprocess" and "to_model_dict" helper functions

# Todo: add the model_dict run of flowgraph

# NB: THere is no real comparison of output files because the answer files
#     have different units. Could create a class to manage this but
#     I've checked that the memory values match the file values.
#     do_compare_output_files write files but does not check them
do_compare_output_files = True
do_compare_in_memory = True
rtol = atol = 1.0e-7


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
    if simulation["name"][0:7] != "ucb_2yr":
        pytest.skip("Only testing prms_channel_flow_graph for ucb_2yr")

    control = Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )
    control.edit_n_time_steps(180)
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


@pytest.fixture(scope="function")
def big_sandy_parameters():
    param_file = pl.Path("../test_data/starfit/istarf_conus_grand.nc")
    bs_grand_id = 419
    params = xr.open_dataset(param_file)
    params = params.where(params.grand_id == bs_grand_id, drop=True)
    return Parameters.from_ds(params)


# TODO fixture for daily/hourly
@pytest.mark.parametrize(
    "compute_daily", [True, False], ids=("daily", "hourly")
)
def test_starfit_flow_graph_postprocess(
    simulation,
    control,
    discretization,
    parameters,
    big_sandy_parameters,
    compute_daily,
    tmp_path,
):
    input_dir = simulation["output_dir"]

    # are we testing flowgraph netcdf output elsewhere?
    # i dont think so, so do it here.
    # cant really test answers per above note
    control.options["input_dir"] = input_dir
    control.options["budget_type"] = "error"
    control.options["verbosity"] = True
    # control.options["netcdf_output_dir"] = tmp_path  # TODO
    control.options["netcdf_output_var_names"] = [
        "node_outflows",
        "node_upstream_inflows",
        "node_storages",
    ]

    # Currently this is the same as notebook 06
    flow_graph = prms_channel_flow_graph_postprocess(
        control=control,
        input_dir=input_dir,
        prms_channel_dis=discretization,
        prms_channel_params=parameters,
        new_nodes_maker_dict={
            "starfit": StarfitFlowNodeMaker(
                None,
                big_sandy_parameters,
                budget_type="error",
                compute_daily=compute_daily,
            ),
            "pass_through": PassThroughNodeMaker(),
        },
        new_nodes_maker_names=["starfit", "pass_through"],
        new_nodes_maker_indices=[0, 0],
        new_nodes_flow_to_nhm_seg=[
            44426,
            44435,
        ],
    )

    # get the segments un affected by flow, where the PRMS solutions should
    # match
    wh_44426 = np.where(discretization.parameters["nhm_seg"] == 44426)[0][0]
    upind = wh_44426
    wh_ignore = []
    toseg = discretization.parameters["tosegment"] - 1
    while upind >= 0:
        downind = toseg[upind]
        wh_ignore += [upind]
        upind = downind

    if do_compare_output_files:
        flow_graph.initialize_netcdf(tmp_path)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            if var not in rename_vars.keys():
                # this drops the inflow-vars and renames the rest
                continue
            new_var = rename_vars[var]
            var_pth = input_dir / f"{var}.nc"
            answers[new_var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    for istep in range(control.n_times):
        control.advance()
        flow_graph.advance()
        flow_graph.calculate(1.0)
        flow_graph.output()

        for var in answers.values():
            var.advance()

        # use the results/actual as answers on wh_ingore and on the new nodes
        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers.items():
                if key in convert_to_vol:
                    current = val.current / (24 * 60 * 60)
                else:
                    current = val.current

                # <
                # Fill on the downstream reaches of PRMS, where we "ingnore"
                current[(wh_ignore,)] = flow_graph[key][(wh_ignore,)]
                # Fill in the last two nodes
                answers_conv_vol[key] = np.concatenate(
                    [current, flow_graph[key][-2:]]
                )

            # <<
            # there are no expected sources or sinks in this test
            answers_conv_vol["node_sink_source"] = np.concatenate(
                [val.current * zero, flow_graph["node_sink_source"][-2:]]
            )
            answers_conv_vol["node_negative_sink_source"] = np.concatenate(
                [
                    val.current * zero,
                    flow_graph["node_negative_sink_source"][-2:],
                ]
            )

            answers_conv_vol["node_storages"] = np.concatenate(
                [val.current * nan, flow_graph["node_storages"][-2:]]
            )

            compare_in_memory(
                flow_graph,
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=False,
                fail_after_all_vars=False,
            )

    # this checks that the budget was actually active for the starfit node
    assert flow_graph._nodes[-2].budget is not None

    flow_graph.finalize()


@pytest.mark.parametrize(
    "compute_daily", [True, False], ids=("daily", "hourly")
)
def test_starfit_flow_graph_model_dict(
    simulation,
    control,
    discretization,
    parameters,
    big_sandy_parameters,
    compute_daily,
    tmp_path,
):
    domain_dir = simulation["dir"]
    input_dir = simulation["output_dir"]

    # are we testing flowgraph netcdf output elsewhere?
    # i dont think so, so do it here.
    # cant really test answers per above note
    control.options["input_dir"] = input_dir
    control.options["budget_type"] = "error"
    control.options["verbosity"] = True
    # control.options["netcdf_output_dir"] = tmp_path  # TODO
    control.options["netcdf_output_var_names"] = [
        "node_outflows",
        "node_upstream_inflows",
        "node_storages",
    ]

    nhm_processes = [
        # PRMSSnow,  #  snow introduces significant error
        PRMSRunoff,
        PRMSSoilzone,
        PRMSGroundwater,
    ]

    model_dict = {
        "control": control,
        "dis_both": discretization,
        "dis_hru": discretization,
        "model_order": [],
    }

    # As in notebook 01
    for proc in nhm_processes:
        proc_name = proc.__name__
        proc_rename = "prms_" + proc_name[4:].lower()
        model_dict["model_order"] += [proc_rename]
        model_dict[proc_rename] = {}
        proc_dict = model_dict[proc_rename]
        proc_dict["class"] = proc
        proc_param_file = domain_dir / f"parameters_{proc_name}.nc"
        proc_dict["parameters"] = Parameters.from_netcdf(proc_param_file)
        proc_dict["dis"] = "dis_hru"

    # <
    new_nodes_maker_dict = {
        "starfit": StarfitFlowNodeMaker(
            None,
            big_sandy_parameters,
            budget_type="error",
            compute_daily=compute_daily,
        ),
        "pass_through": PassThroughNodeMaker(),
    }
    new_nodes_maker_names = ["starfit", "pass_through"]
    new_nodes_maker_indices = [0, 0]
    new_nodes_flow_to_nhm_seg = [
        44426,
        44435,
    ]

    model_dict = prms_channel_flow_graph_to_model_dict(
        model_dict=model_dict,
        prms_channel_dis=discretization,
        prms_channel_dis_name="dis_both",
        prms_channel_params=parameters,
        new_nodes_maker_dict=new_nodes_maker_dict,
        new_nodes_maker_names=new_nodes_maker_names,
        new_nodes_maker_indices=new_nodes_maker_indices,
        new_nodes_flow_to_nhm_seg=new_nodes_flow_to_nhm_seg,
        graph_budget_type="error",  # move to error
    )
    model = Model(model_dict)

    # on the segments unaffected by flow, the PRMS solutions should match
    # so "ignore" the affected segments, we'll use the results as answers later
    wh_44426 = np.where(discretization.parameters["nhm_seg"] == 44426)[0][0]
    upind = wh_44426
    wh_ignore = []
    toseg = discretization.parameters["tosegment"] - 1
    while upind >= 0:
        downind = toseg[upind]
        wh_ignore += [upind]
        upind = downind

    if do_compare_output_files:
        # not really feasible as noted by NB section at top
        model.initialize_netcdf(tmp_path)

    if do_compare_in_memory:
        answers = {}
        for var in PRMSChannel.get_variables():
            if var not in rename_vars.keys():
                continue
            new_var = rename_vars[var]
            var_pth = input_dir / f"{var}.nc"
            answers[new_var] = adapter_factory(
                var_pth, variable_name=var, control=control
            )

    # short cut
    flow_graph = model.processes["prms_channel_flow_graph"]

    for istep in range(control.n_times):
        model.advance()
        model.calculate()
        model.output()

        for var in answers.values():
            var.advance()

        # use the results as answers on wh_ingore and on the new nodes
        if do_compare_in_memory:
            answers_conv_vol = {}
            for key, val in answers.items():
                if key in convert_to_vol:
                    current = val.current / (24 * 60 * 60)
                else:
                    current = val.current

                # <
                # Fill on the downstream reaches of PRMS
                current[(wh_ignore,)] = flow_graph[key][(wh_ignore,)]
                # Fill in the last two nodes
                answers_conv_vol[key] = np.concatenate(
                    [current, flow_graph[key][-2:]]
                )

            # <<
            # there are no expected sources or sinks in this test
            answers_conv_vol["node_sink_source"] = np.concatenate(
                [val.current * zero, flow_graph["node_sink_source"][-2:]]
            )
            answers_conv_vol["node_negative_sink_source"] = np.concatenate(
                [
                    val.current * zero,
                    flow_graph["node_negative_sink_source"][-2:],
                ]
            )

            answers_conv_vol["node_storages"] = np.concatenate(
                [val.current * nan, flow_graph["node_storages"][-2:]]
            )

            compare_in_memory(
                flow_graph,
                answers_conv_vol,
                atol=atol,
                rtol=rtol,
                skip_missing_ans=False,
                fail_after_all_vars=False,
                verbose=True,
            )

    # this checks that the budget was actually active for the starfit node
    assert flow_graph._nodes[-2].budget is not None

    flow_graph.finalize()
