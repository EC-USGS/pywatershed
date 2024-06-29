from typing import Literal

import networkx as nx
import numpy as np

from pywatershed.base.accessor import Accessor
from pywatershed.base.adapter import adaptable
from pywatershed.base.conservative_process import ConservativeProcess
from pywatershed.base.control import Control
from pywatershed.constants import nan, zero
from pywatershed.parameters import Parameters


class FlowNode(Accessor):
    """The FlowNode base class.

    A FlowNode represents a spatial element of an explicit flow solution which
    does not (currently) include a head (water depth) term.

    A FlowNode is instantiated with its own (optional) data and calculates
    outflow, storage_change, and sink_source properties on subtimesteps.
    """

    def __init__(self, control: Control):
        """Initialize the FlowNode.

        Args:
          control: A Control object.
        """
        raise Exception("This must be overridden")

    def prepare_timestep(self):
        "Prepare the subtimestep for subtimestep calculations."
        raise Exception("This must be overridden")

    def calculate_subtimestep(
        self, isubstep: int, inflow_upstream: float, inflow_lateral: float
    ):
        """Calculate the subtimestep.

        Args:
          isubstep: Zero-based integer indicating the index of the current
            substep.
          inflow_upstream: The in-channel flows to this FlowNode on the current
            substep.
          inflow_lateral: The later flows to this FlowNode on the current
            substep.
        """
        raise Exception("This must be overridden")

    def advance(self):
        "Advance this FlowNode to the next timestep."
        raise Exception("This must be overridden")

    def finalize_timestep(self):
        "Finalize the current timestep at this FlowNode."
        raise Exception("This must be overridden")

    @property
    def outflow(self):
        "The average outflow of the FlowNode over the current timestep."
        raise Exception("This must be overridden")

    @property
    def outflow_substep(self):
        """The outflow of the FlowNode over the sub-timestep."""
        raise Exception("This must be overridden")

    @property
    def storage_change(self):
        "The storage change of the FlowNode at the current subtimestep."
        raise Exception("This must be overridden")

    @property
    def storage(self):
        "The storage of the FlowNode at the current subtimestep."
        raise Exception("This must be overridden")

    @property
    def sink_source(self):
        "The sink or source amount of the FlowNode at the current subtimestep."
        raise Exception("This must be overridden")


class FlowNodeMaker(Accessor):
    """FlowNodeMaker instantiates FlowNodes with their data.

    Notes:
      * initialize the the data needed to return the nodes
      * method to give a node at index
      * the node maker is only used on the init of the FlowGraph
      * does not track mass balance

    """

    def __init__(
        self,
        discretization: Parameters = None,
        parameters: Parameters = None,
    ):
        """Intitalize the FlowNodeMaker

        Args:
          discretization: Discretization data to parcel out to the FlowNodes.
          parameters: Parmeter data to parcel to the FlowNodes.
        """
        self.name = "FlowNodeMaker"
        return

    def get_node(control: Control, index: int) -> FlowNode:
        """Instantiate FlowNode at a given index.

        Args:
          control: A Control object.
          index: The index in the discretization and parameter data to use
            when instantiating a FlowNode.
        """
        raise Exception("This must be overridden")


class FlowGraph(ConservativeProcess):
    """FlowGraph manages and computes FlowNodes given by FlowNodeMakers.

    FlowGraph lets users combine FlowNodes of different kinds into a single
    mathmetical graph of flow solution. FlowNodes provide explicit solutions
    of flow (currently not involving a head term) on a single spatial unit. The
    FlowGraph allows these different flow solutions to be combined in arbitrary
    order on a mathematical graph. There are many applications, but a common
    one is to add a reservoir representation into an existing graph of flow,
    such as exists within PRMSChannel which computes a Muskingum-Mann solution
    of flow. This example is shown schematically in the following figure.

    .. |fg1| image:: /_static/flow_graph_schematic_1.png
       :align: middle
    +---------+
    |  |fg1|  |
    +---------+

    Above a node of class B is inserted into the original graph. Class B may
    have a different flow solution than class A in the original graph, but
    FlowGraph handles new nodes wherever you want to put them. FlowGraph checks
    mass balance over the graph.

    To delve a bit deeper, the relationship between FlowGraph, FlowNode, and
    FlowNodeMaker is shown in the figure below.

    .. |fg2| image:: /_static/flow_graph_schematic_2.png
       :align: middle
    +---------+
    |  |fg2|  |
    +---------+

    The figure above illustrates how FlowNodeMakers already have a certain
    kind of FlowNode class composed into them. A user instantiates each
    FlowNodeMaker by passing all the data required for all the FlowNodes.
    FlowGraph recieves instantiated FlowNodeMakers and calls them, in turn,
    instantiate the FlowNodes in the FlowGraph.

    Note that users generally do not create types of FlowNodes or
    FlowNodeMakers themselves, this is typically the work of code developers.
    But users pass parameters, an inflow Adapter, and instantiated
    FlowNodeMakers to FlowGraph. The example below shows the nuts and bolts of
    setting up a FlowGraph similar to that illustrated above, where a single
    pass-through node is inserted.

    For users specifically interested in adding new nodes into the PRMSChannel
    MuskingumMann routing solutions, there are helper functions available which
    greatly simplify the code. See the notebook
    `examples/06_flow_graph_starfit.ipynb <https://github.com/EC-USGS/pywatershed/blob/develop/examples/06_flow_graph_starfit.ipynb>`__
    which highlights both  helper functions
    :func:`prms_channel_flow_graph_to_model_dict`
    and :func:`prms_channel_flow_graph_postprocess`.

    Examples:
    ---------

    This example shows how to insert a pass-through node into a PRMSChannel
    simulation. It's a bit underwhelming because the flows on the PRMSChannel
    nodes are unaltered, but shows the full mechanism without helper functions.

    >>> import numpy as np
    >>> from tqdm.auto import tqdm
    >>> import xarray as xr
    >>> import pywatershed as pws
    >>> from pywatershed.constants import nan, zero
    >>> from pywatershed.constants import __pywatershed_root__ as pkg_root_dir
    >>> # this example requries the repository with test data previously generated
    >>> domain_dir = pkg_root_dir / "../test_data/drb_2yr"
    >>> control_file = domain_dir / "nhm.control"
    >>> control = pws.Control.load_prms(
    ...     control_file, warn_unused_options=False
    ... )
    >>> dis_hru_file = domain_dir / "parameters_dis_hru.nc"
    >>> dis_seg_file = domain_dir / "parameters_dis_seg.nc"
    >>> discretization_prms = pws.Parameters.merge(
    ...     pws.Parameters.from_netcdf(dis_hru_file, encoding=False),
    ...     pws.Parameters.from_netcdf(dis_seg_file, encoding=False),
    ... )
    >>> param_file = domain_dir / "parameters_PRMSChannel.nc"
    >>> parameters_prms = pws.parameters.PrmsParameters.from_netcdf(param_file)
    >>> # Build the parameters for the FlowGraph
    >>> nnodes = parameters_prms.dims["nsegment"] + 1
    >>> node_maker_name = ["prms_channel"] * nnodes
    >>> node_maker_name[-1] = "pass_throughs"
    >>> node_maker_index = np.arange(nnodes)
    >>> node_maker_index[-1] = 0
    >>> to_graph_index = np.zeros(nnodes, dtype=np.int64)
    >>> dis_params = discretization_prms.parameters
    >>> to_graph_index[0:-1] = dis_params["tosegment"] - 1
    >>> nhm_seg_intervene_above = 1829
    >>> wh_intervene_above_nhm = np.where(
    ...     dis_params["nhm_seg"] == nhm_seg_intervene_above
    ... )
    >>> wh_intervene_below_nhm = np.where(
    ...     (dis_params["tosegment"] - 1) == wh_intervene_above_nhm[0][0]
    ... )
    ... # have to map to the graph from an index found in prms_channel
    >>> wh_intervene_above_graph = np.where(
    ...     (np.array(node_maker_name) == "prms_channel")
    ...     & (node_maker_index == wh_intervene_above_nhm[0][0])
    ... )
    >>> wh_intervene_below_graph = np.where(
    ...     (np.array(node_maker_name) == "prms_channel")
    ...     & np.isin(node_maker_index, wh_intervene_below_nhm)
    ... )
    >>> to_graph_index[-1] = wh_intervene_above_graph[0][0]
    >>> to_graph_index[wh_intervene_below_graph] = nnodes - 1
    >>> parameters_flow_graph = pws.Parameters(
    ...     dims={
    ...         "nnodes": nnodes,
    ...     },
    ...     coords={
    ...         "node_coord": np.arange(nnodes),
    ...     },
    ...     data_vars={
    ...         "node_maker_name": node_maker_name,
    ...         "node_maker_index": node_maker_index,
    ...         "to_graph_index": to_graph_index,
    ...     },
    ...     metadata={
    ...         "node_coord": {"dims": ["nnodes"]},
    ...         "node_maker_name": {"dims": ["nnodes"]},
    ...         "node_maker_index": {"dims": ["nnodes"]},
    ...         "to_graph_index": {"dims": ["nnodes"]},
    ...     },
    ...     validate=True,
    ... )
    >>> # Get the FlowNodeMakers instantiated and named
    >>> node_maker_dict = {
    ...     "prms_channel": pws.PRMSChannelFlowNodeMaker(
    ...         discretization_prms, parameters_prms
    ...     ),
    ...     "pass_throughs": pws.PassThroughNodeMaker(),
    ... }
    >>> # Get the inputs to PRMSChannel combined, then add inputs to the
    ... # additional node using a custom Adapter.
    >>> input_variables = {}
    >>> for key in pws.PRMSChannel.get_inputs():
    ...     nc_path = domain_dir / f"output/{key}.nc"
    ...     input_variables[key] = pws.AdapterNetcdf(nc_path, key, control)
    ...
    >>> inflows_prms = pws.HruSegmentFlowAdapter(
    ...     parameters_prms, **input_variables
    ... )
    >>> class GraphInflowAdapter(pws.Adapter):
    ...     def __init__(
    ...         self,
    ...         prms_inflows: pws.Adapter,
    ...         variable: str = "inflows",
    ...     ):
    ...         self._variable = variable
    ...         self._prms_inflows = prms_inflows
    ...
    ...         self._nnodes = len(self._prms_inflows.current) + 1
    ...         self._current_value = np.zeros(self._nnodes) * nan
    ...         return
    ...
    ...     def advance(self) -> None:
    ...         self._prms_inflows.advance()
    ...         self._current_value[0:-1] = self._prms_inflows.current
    ...         self._current_value[-1] = zero  # no inflow at the pass through
    ...         return
    ...
    >>> inflows_graph = GraphInflowAdapter(inflows_prms)
    >>> # Instantiate the FlowGraph
    >>> flow_graph = pws.FlowGraph(
    ...     control,
    ...     discretization=None,
    ...     parameters=parameters_flow_graph,
    ...     inflows=inflows_graph,
    ...     node_maker_dict=node_maker_dict,
    ...     budget_type="error",
    ... )
    >>> # Save out the full timeseries of flows for all nodes
    >>> graph_seg_outflows = np.zeros([control.n_times, nnodes])
    >>> # Run the flow graph
    >>> for istep in tqdm(range(control.n_times)):
    ...     control.advance()
    ...     flow_graph.advance()
    ...     flow_graph.calculate(1.0)
    ...     graph_seg_outflows[istep, :] = flow_graph["node_outflows"]
    ...
    >>> flow_graph.finalize()
    >>> # Compare to the results of PRMSChannel run with out a pass-through
    ... # node.
    >>> prms_seg_outflows = xr.open_dataarray(
    ...     domain_dir / "output/seg_outflow.nc"
    ... )
    ... # The final node is the passthrough node, drop it from comparisons.
    >>> assert (
    ...     abs(graph_seg_outflows[:, 0:-1] - prms_seg_outflows.values) < 1e-10
    ... ).all()

    """  # noqa: E501

    def __init__(
        self,
        control: Control,
        discretization: Parameters,  # could use this, but not necsesary
        parameters: Parameters,
        inflows: adaptable,
        node_maker_dict: dict,
        budget_type: Literal["defer", None, "warn", "error"] = "defer",
        verbose: bool = None,
    ):
        """Initialize a FlowGraph.

        Args:
            control: A Control object
            discretization: Currently unused by FlowGraph but required by
              it's superclass, ConservativeProcess.
            parameters: A Parameter object with the FlowGraph parameters as
              described below.
            inflows: An adaptable of inflows to the graph, often referred to as
              "lateral" flows (not flows inside the graph).
            node_maker_dict: A dictionary of FlowNodeMaker instances with
              keys/names supplied in the parameters, e.g.
              {key1: flow_node_maker_instance, ...}.
            budget_type: one of ["defer", None, "warn", "error"] with "defer"
              being the default and defering to
              control.options["budget_type"] when
              available. When control.options["budget_type"] is not avaiable,
              budget_type is set to "warn".

        The `parameters` argument is a :class:`Parameters` object which
        contains the following data:

        * node_maker_name: A list or np.array of the FlowNodeMaker name for
          each node.
        * node_maker_index: A list or np.array of the indices to ask for from
          the associated/collated FlowNodeMaker (above) for each node
        * to_graph_index: the index of the downstream index in the FlowGraph
          with -1 indicating an outflow node. This must specify a DAG.

        The inputs inflows, node_maker_name, node_maker_index, and
        to_graph_index are collated. The order of execution of the graph is not
        the same as the supplied order, the execution order is solved from
        to_graph_index. Note that initial conditions are set by the node makers
        via their parameters.

        """
        super().__init__(
            control=control,
            discretization=discretization,
            parameters=parameters,
        )
        self.name = "FlowGraph"

        self._set_inputs(locals())
        self._set_options(locals())

        for fnm in self._node_maker_dict.values():
            assert isinstance(fnm, FlowNodeMaker)

        self._init_graph()

        # If/when FlowGraph handles nodes which dont tautologically balance
        # could allow the basis to be unit.
        self._set_budget(basis="global", unit_desc="flow rates")

        return

    @staticmethod
    def get_dimensions() -> tuple:
        return ("nnodes",)

    @staticmethod
    def get_parameters() -> tuple:
        return (
            "node_maker_name",
            "node_maker_index",
            "to_graph_index",
        )

    @staticmethod
    def get_inputs() -> tuple:
        return ("inflows",)

    @staticmethod
    def get_init_values() -> dict:
        return {
            "outflows": nan,
            "node_upstream_inflows": nan,
            "node_outflows": nan,
            "node_storage_changes": nan,
            "node_storages": nan,
            "node_sink_source": nan,
            "node_negative_sink_source": nan,
        }

    @classmethod
    def get_variables(cls) -> tuple:
        return list(cls.get_init_values().keys())

    @staticmethod
    def get_mass_budget_terms():
        """Get a dictionary of variable names for mass budget terms."""
        return {
            "inputs": ["inflows"],
            "outputs": ["outflows"],
            "storage_changes": [
                "node_storage_changes",
                "node_negative_sink_source",
            ],
        }

    def get_outflow_mask(self):
        """Get a mask indicataing on which nodes flow exits the graph."""
        return self._outflow_mask

    @property
    def outflow_mask(self):
        "A mask indicating on which nodes flow exits the graph."
        return self._outflow_mask

    def _set_initial_conditions(self) -> None:
        self._node_upstream_inflow_sub = np.zeros(self.nnodes) * nan
        self._node_upstream_inflow_acc = np.zeros(self.nnodes) * nan
        self._node_outflow_substep = np.zeros(self.nnodes) * nan
        return

    def _init_graph(self) -> None:
        params = self._params.parameters
        # where do flows exit the graph?
        self._outflow_mask = np.where(
            params["to_graph_index"] == -1, True, False
        )

        # which nodes do not have upstream nodes?
        self._headwater_nodes = set(range(self.nnodes)).difference(
            set(params["to_graph_index"])
        )

        connectivity = []
        for inode in range(self.nnodes):
            tonode = params["to_graph_index"][inode]
            if tonode < 0:
                continue
            connectivity.append(
                (
                    inode,
                    tonode,
                )
            )

        # use networkx to calculate the Directed Acyclic Graph
        self._graph = nx.DiGraph()
        self._graph.add_edges_from(connectivity)

        # Make sure the user isnt suppling disconnected nodes
        assert len(self._graph) == self.nnodes, "Disconnected nodes present."
        # should we allow disconnected nodes? seems not
        # these are handled in PRMSChannel if we change our minds

        if self.nnodes > 1:
            node_order = list(nx.topological_sort(self._graph))
        else:
            node_order = [0]

        self._node_order = np.array(node_order, dtype="int64")

        # any performance for doing a hash table up front?
        # a hash {to_seg: [from_seg_0, ..., from_seg_n]}

        # instatiate the nodes
        self._nodes = []
        for ii, (maker_name, maker_index) in enumerate(
            zip(params["node_maker_name"], params["node_maker_index"])
        ):
            self._nodes += [
                self._node_maker_dict[maker_name].get_node(
                    self.control, maker_index
                )
            ]

    def _advance_variables(self) -> None:
        for node in self._nodes:
            node.advance()

        # no prognostic variables on the graph
        return

    def calculate(self, time_length: float, n_substeps: int = 24) -> None:
        params = self._params.parameters

        for node in self._nodes:
            node.prepare_timestep()

        self._node_upstream_inflow_acc[:] = zero

        for istep in range(n_substeps):
            # This works because the first nodes calculated do
            # not have upstream reaches
            self._node_upstream_inflow_sub[:] = zero

            for inode in self._node_order:
                # The first nodes calculated dont have upstream inflows
                # Eventually pass timestep length and n_substems to nodes
                # Calculate
                self._nodes[inode].calculate_subtimestep(
                    istep,
                    self._node_upstream_inflow_sub[inode],
                    self.inflows[inode],
                )
                # Get the outflows back
                self._node_outflow_substep[inode] = self._nodes[
                    inode
                ].outflow_substep
                # Add this node's outflow its downstream node's inflow
                if params["to_graph_index"][inode] >= 0:
                    self._node_upstream_inflow_sub[
                        params["to_graph_index"][inode]
                    ] += self._node_outflow_substep[inode]

            # <
            # not sure how PRMS-specific this is
            self._node_upstream_inflow_acc += self._node_upstream_inflow_sub

        for node in self._nodes:
            node.finalize_timestep()

        self.node_upstream_inflows[:] = (
            self._node_upstream_inflow_acc / n_substeps
        )

        for ii in range(self.nnodes):
            self.node_outflows[ii] = self._nodes[ii].outflow
            self.node_storage_changes[ii] = self._nodes[ii].storage_change
            self.node_storages[ii] = self._nodes[ii].storage
            self.node_sink_source[ii] = self._nodes[ii].sink_source

        self.node_negative_sink_source[:] = -1 * self.node_sink_source

        # global mass balance term
        self.outflows[:] = np.where(
            self._outflow_mask, self.node_outflows, zero
        )

        if self.budget is not None:
            self.budget.advance()
            self.budget.calculate()

        return


def inflow_exchange_factory(
    dimension_names: tuple,
    parameter_names: tuple,
    input_names: tuple,
    init_values: dict,
    mass_budget_terms: dict,
    calculation: callable,
):
    """Create and InflowExchange class from function input.

    Args:
        dimension_names: tuple of dimension names of the InflowExchange
        parameter_names: tuple of parameter names of the InflowExchange
        input_names: tuple of input names for the InflowExchange
        init_values: dict of variable names and values for the public variables
        mass_budget_terms: dict of inputs, outputs, and storage_changes keys
            with a list of the variables in init_values from which to
            calculate mass balance.
        calculation: a function on self that performs the calculations.

    Returns:
        A Class of InflowExchange, a subclass of ConservativeProcess.

    """

    class InflowExchange(ConservativeProcess):
        def __init__(
            self,
            control: Control,
            discretization: Parameters,
            parameters: Parameters,
            budget_type: Literal[None, "warn", "error"] = None,
            verbose: bool = None,
            budget_basis="global",
            **kwargs,
        ) -> None:
            super().__init__(
                control=control,
                discretization=discretization,
                parameters=parameters,
            )
            self.name = "InflowExchange"

            self._set_inputs(locals() | kwargs)
            self._set_options(locals())

            self._set_budget(basis=budget_basis, unit_desc="flow rates")

            return

        @staticmethod
        def get_dimensions() -> tuple:
            return dimension_names

        @staticmethod
        def get_parameters() -> tuple:
            return parameter_names

        @staticmethod
        def get_inputs() -> tuple:
            return input_names

        @staticmethod
        def get_init_values() -> dict:
            return init_values

        @staticmethod
        def get_mass_budget_terms():
            return mass_budget_terms

        def _set_initial_conditions(self) -> None:
            pass

        def _advance_variables(self) -> None:
            pass

        def _calculate(self, simulation_time: float) -> None:
            calculation(self)
            return

    return InflowExchange
