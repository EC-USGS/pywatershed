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
from pywatershed.hydrology.pass_through_flow_node import (
    PassThroughFlowNodeMaker,
)
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


@pytest.fixture(
    params=[True, False], ids=("daily", "hourly"), scope="function"
)
def compute_daily(request):
    return request.param


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

    # We'll test adding multiple new nodes in-series into a FlowGraph. Above
    # and below the starfit we'll add pass-through nodes and check these
    # match.

    # We add a random passthrough node with no upstream node, to test if that
    # works.

    # in this test we'll cover the case of disconnected nodes and allowing
    # them or not
    def add_disconnected_node_data(ds: xr.Dataset) -> xr.Dataset:
        # This just takes data from the first segment and repeats it
        # at the end of the array. Editing the values is done outside this
        # function.
        seg_data = {}
        for vv in list(ds.variables):
            if "nsegment" in ds[vv].dims:
                data = ds[vv].values
                seg_data[vv] = (["nsegment"], np.append(data, data[0]))
                del ds[vv]

        nhm_seg = seg_data.pop("nhm_seg")
        ds_seg = xr.Dataset(
            data_vars=seg_data,
            coords={"nhm_seg": nhm_seg},
        )
        return xr.merge([ds, ds_seg])

    dis_ds = add_disconnected_node_data(discretization.to_xr_ds())
    params_ds = add_disconnected_node_data(parameters.to_xr_ds())
    dis_ds.tosegment[-1] = 0
    dis_ds.tosegment_nhm[-1] = 99999
    dis_ds.nhm_seg[-1] = 99998
    params_ds.nhm_seg[-1] = 99998

    discretization = Parameters.from_ds(dis_ds)
    parameters = Parameters.from_ds(params_ds)

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

    new_nodes_maker_names = ["starfit"] + ["pass_through"] * 3
    new_nodes_maker_indices = [0, 0, 1, 2]
    new_nodes_maker_ids = [-2, -1, -100, -1000]
    # The starfit node flows to the third passthrough node, in index 3.
    # The first passthrough node flows to some random nhm_seg, not connected to
    # the other new nodes.
    # The second passthrough flows to the starfit node in index 0.
    # The last passthrough node flows to the seg above which the reservoir
    # is placed.
    new_nodes_flow_to_nhm_seg = [-3, 44409, 0, 44426]

    # the first in the list is for the disconnected node
    check_names = ["prms_channel"] + new_nodes_maker_names
    check_indices = [dis_ds.dims["nsegment"] - 1] + new_nodes_maker_indices
    check_ids = [dis_ds.nhm_seg[-1].values.tolist()] + new_nodes_maker_ids

    # This warning should say: TODO
    with pytest.warns(
        UserWarning, match="Disconnected nodes present in FlowGraph."
    ):
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
                "pass_through": PassThroughFlowNodeMaker(),
            },
            new_nodes_maker_names=new_nodes_maker_names,
            new_nodes_maker_indices=new_nodes_maker_indices,
            new_nodes_maker_ids=new_nodes_maker_ids,
            new_nodes_flow_to_nhm_seg=new_nodes_flow_to_nhm_seg,
            addtl_output_vars=["spill", "release"],
            allow_disconnected_nodes=True,
        )

    # <
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

        # use the results/actual as answers on wh_ingore and on the new nodes
        if do_compare_in_memory:
            for var in answers.values():
                var.advance()

            answers_conv_vol = {}
            # The -5 is not -4 because of the synthetic disconnected node also
            # present on the FlowGraph (which is added via parameters, not
            # using prms_channel_flow_graph_postprocess)
            for key, val in answers.items():
                if key in convert_to_vol:
                    current = val.current / (24 * 60 * 60)
                else:
                    current = val.current

                # <
                # Fill the first node, the fake disconnected node, with
                # whatever value is in the graph
                # current = np.append(current, flow_graph[key][-3:-2])
                # Fill on the downstream reaches of PRMS, where we "ingnore"
                # wh_ignore is on the extended dis/params but shouldnt make
                # a difference with original since original data in same
                # indices with new, disconnected node at the end.
                current[(wh_ignore,)] = flow_graph[key][(wh_ignore,)]
                # Fill in the last two nodes
                answers_conv_vol[key] = np.concatenate(
                    [current, flow_graph[key][-5:]]
                )

            # <<
            # there are no expected sources or sinks in this test
            answers_conv_vol["node_sink_source"] = np.concatenate(
                [val.current * zero, flow_graph["node_sink_source"][-5:]]
            )
            answers_conv_vol["node_negative_sink_source"] = np.concatenate(
                [
                    val.current * zero,
                    flow_graph["node_negative_sink_source"][-5:],
                ]
            )

            answers_conv_vol["node_storages"] = np.concatenate(
                [val.current * nan, flow_graph["node_storages"][-5:]]
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
    assert flow_graph._nodes[-4].budget is not None

    flow_graph.finalize()
    for vv in control.options["netcdf_output_var_names"]:
        da = xr.load_dataarray(tmp_path / f"{vv}.nc", concat_characters=True)
        assert da.node_maker_index[-5:].values.tolist() == check_indices
        assert da.node_maker_name[-5:].values.tolist() == check_names
        assert da.node_maker_id[-5:].values.tolist() == check_ids

    # test additional output files match and are working
    da_no = xr.load_dataarray(tmp_path / "node_outflows.nc")
    da_in = xr.load_dataarray(tmp_path / "node_upstream_inflows.nc")

    # full time check of the passthrough nodes (which is probably gratuitious
    # given all the above checks already passed.
    assert (abs(da_no[:, -1] - da_no[:, -4]) < 1e-12).all()
    assert (abs(da_in[:, -2] - da_in[:, -4]) < 1e-12).all()

    da_lr = xr.load_dataarray(tmp_path / "release.nc")
    da_ls = xr.load_dataarray(tmp_path / "spill.nc")

    da_no = da_no.where(da_no.node_maker_name == "starfit", drop=True)
    da_lr = da_lr.where(da_lr.node_maker_name == "starfit", drop=True)
    da_ls = da_ls.where(da_ls.node_maker_name == "starfit", drop=True)

    assert (da_no == da_lr + da_ls).all()


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
        "pass_through": PassThroughFlowNodeMaker(),
    }
    new_nodes_maker_names = ["starfit", "pass_through"]
    new_nodes_maker_indices = [0, 0]
    new_nodes_maker_ids = [-2, -1]
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
        new_nodes_maker_ids=new_nodes_maker_ids,
        new_nodes_flow_to_nhm_seg=new_nodes_flow_to_nhm_seg,
        graph_budget_type="error",  # move to error
        addtl_output_vars=["spill", "release"],
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
        model.initialize_netcdf(tmp_path, separate_files=False)

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

        # use the results as answers on wh_ingore and on the new nodes
        if do_compare_in_memory:
            for var in answers.values():
                var.advance()

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
                verbose=False,
            )

    # this checks that the budget was actually active for the starfit node
    assert flow_graph._nodes[-2].budget is not None

    flow_graph.finalize()

    # test single file output has extra coords and additional vars
    ds = xr.open_dataset(tmp_path / "FlowGraph.nc")
    ds_starfit = ds.where(ds.node_maker_name == "starfit", drop=True)
    assert (
        ds_starfit.node_outflows == ds_starfit.release + ds_starfit.spill
    ).all()
