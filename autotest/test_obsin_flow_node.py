import numpy as np
import pytest

from pywatershed import PRMSChannel
from pywatershed.base.adapter import Adapter, AdapterNetcdf, adapter_factory
from pywatershed.base.control import Control
from pywatershed.base.flow_graph import FlowGraph
from pywatershed.base.parameters import Parameters
from pywatershed.constants import nan, zero
from pywatershed.hydrology.obsin_node import ObsInNodeMaker
from pywatershed.hydrology.prms_channel_flow_graph import (
    HruSegmentFlowAdapter,
    PRMSChannelFlowNodeMaker,
)
from pywatershed.parameters import PrmsParameters
from pywatershed.utils.prms_streamflow_data import PRMSStreamflowData

do_compare_output_files = False
do_compare_in_memory = True
rtol = atol = 1.0e-7

# test the exchange too?
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
    if "drb_2yr:nhm_obsin" not in simulation["name"]:
        pytest.skip("Only testing obsin flow graph for drb_2yr:nhm_obsin")
    return Control.load_prms(
        simulation["control_file"], warn_unused_options=False
    )


@pytest.fixture(scope="function")
def discretization_prms(simulation):
    dis_hru_file = simulation["dir"] / "parameters_dis_hru.nc"
    dis_seg_file = simulation["dir"] / "parameters_dis_seg.nc"
    return Parameters.merge(
        Parameters.from_netcdf(dis_hru_file, encoding=False),
        Parameters.from_netcdf(dis_seg_file, encoding=False),
    )


@pytest.fixture(scope="function")
def parameters_prms(simulation, control):
    param_file = simulation["dir"] / "parameters_PRMSChannel.nc"
    return PrmsParameters.from_netcdf(param_file)


def test_prms_channel_obsin_compare_prms(
    simulation,
    control,
    discretization_prms,
    parameters_prms,
    tmp_path,
):
    control_param_file = (
        simulation["control_file"].parent / control.options["parameter_file"]
    )
    control_parameters = PrmsParameters.load(control_param_file)
    obsout_seg = control_parameters.parameters["obsout_segment"] - 1
    sf_data = PRMSStreamflowData(simulation["dir"] / "sf_data").data
    poi_inds = obsout_seg[np.where(obsout_seg >= 0)].tolist()
    npoi = len(poi_inds)
    poi_ids = discretization_prms.parameters["poi_gage_id"][(poi_inds),]
    obsin_data = sf_data[sf_data.columns.intersection(poi_ids)]
    # see test starfit flow node to see how slice off individual points
    # and re-merge using xarray
    obsin_params = Parameters.from_ds(
        discretization_prms.subset(["poi_gage_id"])
        .to_xr_ds()
        .isel(npoigages=poi_inds)
    )

    # combine PRMS lateral inflows to a single non-volumetric inflow
    output_dir = simulation["output_dir"]
    input_variables = {}
    for key in PRMSChannel.get_inputs():
        nc_path = output_dir / f"{key}.nc"
        input_variables[key] = AdapterNetcdf(nc_path, key, control)

    inflows_prms = HruSegmentFlowAdapter(parameters_prms, **input_variables)

    class GraphInflowAdapter(Adapter):
        def __init__(
            self,
            prms_inflows: Adapter,
            npoi: int,
            variable: str = "inflows",
        ):
            self._variable = variable
            self._prms_inflows = prms_inflows

            self._npoi = npoi
            self._nnodes = len(self._prms_inflows.current) + npoi
            self._current_value = np.zeros(self._nnodes) * nan
            return

        def advance(self) -> None:
            self._prms_inflows.advance()
            self._current_value[0:-npoi] = self._prms_inflows.current
            self._current_value[-npoi:] = zero  # no inflow at the obsin
            return

    inflows_graph = GraphInflowAdapter(inflows_prms, npoi)

    # FlowGraph

    node_maker_dict = {
        "prms_channel": PRMSChannelFlowNodeMaker(
            discretization_prms, parameters_prms
        ),
        "obsin": ObsInNodeMaker(obsin_params, obsin_data),
    }
    nnodes = parameters_prms.dims["nsegment"] + npoi
    node_maker_name = ["prms_channel"] * nnodes
    node_maker_name[-npoi:] = ["obsin"] * npoi
    node_maker_index = np.arange(nnodes)
    node_maker_index[-npoi:] = np.arange(npoi)
    to_graph_index = np.zeros(nnodes, dtype=np.int64)
    dis_params = discretization_prms.parameters
    to_graph_index[0:-npoi] = dis_params["tosegment"] - 1

    # this is  baroque, but we want to solve the nhm_seg identifier
    obsin_intervene_inds = (
        discretization_prms.parameters["poi_gage_segment"][(poi_inds),] - 1
    )
    # we want to intervene below this poi
    nhm_seg_intervene_below = discretization_prms.parameters["nhm_seg"][
        (obsin_intervene_inds),
    ]
    wh_intervene_below_nhm = np.where(
        np.isin(dis_params["nhm_seg"], nhm_seg_intervene_below)
    )

    wh_intervene_above_nhm = (dis_params["tosegment"] - 1)[
        wh_intervene_below_nhm
    ]

    # have to map to the graph from an index found in prms_channel
    wh_intervene_above_graph = np.where(
        (np.array(node_maker_name) == "prms_channel")
        & np.isin(node_maker_index, wh_intervene_above_nhm)
    )
    wh_intervene_below_graph = np.where(
        (np.array(node_maker_name) == "prms_channel")
        & np.isin(node_maker_index, wh_intervene_below_nhm)
    )

    to_graph_index[-npoi:] = wh_intervene_above_graph[0]
    to_graph_index[wh_intervene_below_graph] = np.arange(npoi) + (
        nnodes - npoi
    )

    params_flow_graph = Parameters(
        dims={
            "nnodes": nnodes,
        },
        coords={
            "nnodes": np.arange(nnodes),
        },
        data_vars={
            "node_maker_name": node_maker_name,
            "node_maker_index": node_maker_index,
            "to_graph_index": to_graph_index,
        },
        metadata={
            "nnodes": {"dims": ["nnodes"]},
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
        inflows=inflows_graph,
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
        # flow_graph.output()

        # check exchange
        lateral_inflow_answers.advance()
        assert (
            inflows_graph.current[0:-npoi] == lateral_inflow_answers.current
        ).all()

        if do_compare_in_memory:
            for key, val in answers.items():
                val.advance()
                if key in convert_to_vol:
                    desired = val.current / (24 * 60 * 60)
                else:
                    desired = val.current

                actual = flow_graph[key][0:-npoi]
                if key in ["outflows", "node_outflows"]:
                    actual[wh_intervene_below_nhm] = flow_graph[key][-npoi:]

                if key == "node_storage_changes":
                    actual[
                        wh_intervene_below_nhm
                    ] += -flow_graph.node_sink_source[-npoi:]

                np.testing.assert_allclose(
                    actual,
                    desired,
                    atol=atol,
                    rtol=rtol,
                )

    flow_graph.finalize()
