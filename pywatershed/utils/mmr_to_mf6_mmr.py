import pathlib as pl

import flopy
import networkx as nx
import numpy as np

from ..base import Control, meta
from ..constants import SegmentType, fileish, zero
from ..parameters import PrmsParameters
from .mmr_to_mf6 import MmrToMf6


class MmrToMf6Mmr(MmrToMf6):
    """DEAL with the documentation later."""

    def __init__(
        self,
        control_file: fileish = None,
        param_file: fileish = None,
        control: Control = None,
        params: PrmsParameters = None,
        segment_shapefile: fileish = None,
        output_dir: fileish = pl.Path("."),
        bc_binary_files: bool = False,
        bc_flows_combine: bool = False,
        sim_name: str = "mmr_to_mf6",
        inflow_dir: fileish = None,
        inflow_from_PRMS: bool = True,
        # intial flows over ride from file?
        # length_units="meters",
        # time_units="seconds",
        save_flows: bool = False,
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        time_zone="UTC",
        write_on_init: bool = True,
        oc_options: dict = None,
        **kwargs,
    ):
        super().__init__(
            control_file=control_file,
            param_file=param_file,
            control=control,
            params=params,
            segment_shapefile=segment_shapefile,
            output_dir=output_dir,
            bc_binary_files=bc_binary_files,
            bc_flows_combine=bc_flows_combine,
            sim_name=sim_name,
            inflow_dir=inflow_dir,
            inflow_from_PRMS=inflow_from_PRMS,
            save_flows=save_flows,
            # length_units=length_units,
            # time_units=time_units,
            start_time=start_time,
            end_time=end_time,
            time_zone=time_zone,
            write_on_init=write_on_init,
        )

        if oc_options is None:
            oc_options = {}

        parameters = self.params.parameters

        # DISL

        # todo: vertices
        # Bring in segment shapefile to do this sometime
        # # only requires parameter file
        # vertices = [
        #     [0, 0.0, 0.0, 0.0],
        #     [1, 0.0, 1.0, 0.0],
        #     [2, 1.0, 0.0, 0.0],
        #     [3, 2.0, 0.0, 0.0],
        #     [4, 3.0, 0.0, 0.0],
        # ]
        # # icell1d fdc ncvert icvert
        # cell2d = [
        #     [0, 0.5, 2, 1, 2],
        #     [1, 0.5, 2, 0, 2],
        #     [2, 0.5, 2, 2, 3],
        #     [3, 0.5, 2, 3, 4],
        # ]

        # nvert turns off requirement of vertices and cell2d
        nvert = None  # len(vertices)
        vertices = None
        cell2d = None

        self._tosegment = parameters["tosegment"] - 1

        # unit-ed quantities:
        segment_units = self.units(meta.parameters["seg_length"]["units"])
        self._segment_length = parameters["seg_length"]
        self._segment_length = self._segment_length * segment_units

        _ = flopy.mf6.ModflowSwfdisl(
            self._swf,
            nodes=self._nsegment,
            nvert=nvert,
            reach_length=self._segment_length.to_base_units().magnitude,
            toreach=self._tosegment,
            idomain=1,  # ??
            vertices=vertices,
            cell2d=cell2d,
            length_units=self.length_units,
        )

        # set boundary conditions/flows
        self.set_flw()

        # EMS
        _ = flopy.mf6.ModflowEms(self._sim)

        # MMR
        # note: for specifying lake number, use fortran indexing!
        mmr_obs = {
            f"{sim_name}.mmr.obs.csv": [
                ("OUTFLOW1", "EXT-OUTFLOW", 1),
            ],
            "digits": 10,
        }

        # solve segment order - taken from the PRMSChannel initialization
        # convert prms data to zero-based

        connectivity = []
        outflow_mask = np.full((len(self._tosegment)), False)
        for iseg in range(self._nsegment):
            theseg = self._tosegment[iseg]
            if theseg < 0:
                outflow_mask[iseg] = True
                continue
            connectivity.append((iseg, theseg))

        # use networkx to calculate the Directed Acyclic Graph
        if self._nsegment > 1:
            graph = nx.DiGraph()
            graph.add_edges_from(connectivity)
            segment_order = list(nx.topological_sort(graph))
        else:
            segment_order = [0]

        # if the domain contains links with no upstream or
        # downstream reaches, we just throw these back at the
        # top of the order since networkx wont handle such nonsense
        wh_mask_set = set(np.where(outflow_mask)[0])
        seg_ord_set = set(segment_order)
        mask_not_seg_ord = list(wh_mask_set - seg_ord_set)
        if len(mask_not_seg_ord):
            segment_order = mask_not_seg_ord + segment_order
            # for pp in mask_not_seg_ord:
            #    assert (tosegment[pp] == -1) and (not pp in tosegment)

        segment_order = np.array(segment_order, dtype=int)

        # solve k_coef - taken from the PRMSChannel initialization
        # meta.get_units could be defined
        vel_units = {
            key: self.units(val)
            for key, val in meta.get_units(
                ["mann_n", "seg_slope", "seg_depth"], to_pint=True
            ).items()
        }
        vel_terms = {}
        for var in vel_units:
            vel_terms[var] = parameters[var] * vel_units[var]

        velocity = (
            (1.0 / vel_terms["mann_n"])
            * np.sqrt(vel_terms["seg_slope"])
            * vel_terms["seg_depth"] ** (2.0 / 3.0)
            * (3600.0 * self.units("seconds/hour"))
        )

        # JLM: This is a bad idea and should throw an error rather than edit
        # inputs in place during run
        # Shouldnt this ALSO BE ABOVE the velocity calculation?
        # seg_slope = np.where(
        #     parameters["seg_slope"] < 1e-7, 0.0001, parameters["seg_slope"]
        # )  # this is from prms5.2.1, it is not in prms6

        # Kcoef
        # initialize Kcoef to 24.0 for segments with zero velocities
        # this is different from PRMS, which relied on divide by zero resulting
        # in a value of infinity that when evaluated relative to a maximum
        # desired Kcoef value of 24 would be reset to 24. This approach is
        # equivalent and avoids the occurence of a divide by zero.
        Kcoef = np.full(self._nsegment, 24.0, dtype=float)
        Kcoef = Kcoef * (self._segment_length.units / velocity.units)

        # only calculate Kcoef for cells with velocities greater than zero
        idx = velocity > zero
        Kcoef[idx] = self._segment_length[idx] / velocity[idx]
        Kcoef = np.where(
            parameters["segment_type"] == SegmentType.LAKE.value,
            24.0 * Kcoef.units,
            Kcoef,
        )
        Kcoef = np.where(Kcoef.magnitude < 0.01, 0.01 * Kcoef.units, Kcoef)
        Kcoef = np.where(Kcoef.magnitude > 24.0, 24.0 * Kcoef.units, Kcoef)

        # qoutflow
        qoutflow_units = self.units(
            meta.get_units("segment_flow_init", to_pint=True)[
                "segment_flow_init"
            ]
        )
        qoutflow0 = parameters["segment_flow_init"] * qoutflow_units

        # x_coef
        x_coef_units = self.units(
            meta.get_units("x_coef", to_pint=True)["x_coef"]
        )
        x_coef = parameters["x_coef"] * x_coef_units

        _ = flopy.mf6.ModflowSwfmmr(
            self._swf,
            print_flows=True,
            save_flows=save_flows,
            observations=mmr_obs,
            iseg_order=segment_order,
            qoutflow0=qoutflow0.to_base_units().magnitude,
            k_coef=Kcoef.to_base_units().magnitude,
            x_coef=x_coef.to_base_units().magnitude,
        )

        # output control
        _ = flopy.mf6.ModflowSwfoc(
            self._swf,
            budget_filerecord=f"{self._sim_name}.bud",
            qoutflow_filerecord=f"{self._sim_name}.qoutflow",
            **oc_options,
        )

        # self.set_flw()

        # done, write if requested/default else delay
        if write_on_init:
            self.write()
