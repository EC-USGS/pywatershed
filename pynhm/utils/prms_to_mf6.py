import pathlib as pl
from warnings import warn

import flopy
import networkx as nx
import numpy as np
import xarray as xr

# try:
#     import geopandas as gpd

#     has_geopandas = True
# except ModuleNotFoundError:
#     has_geopandas = False


from ..constants import fileish, SegmentType, zero
from .parameters import PrmsParameters
from pynhm import Control

# more generally: this is PRMSMMR to SNFDISL, SNFFLW, and/or SNFMMR


class MMRToMF6:
    """Muskingum-Mann Routing (MMR) data from PRMS to MF6

    The MF6 IO descrption can be found here
    https://modflow6.readthedocs.io/en/latest/mf6io.html

    Using the PRMS parameter files and optionally the PRMS segment shapefiles,
    we can create the following MF6 input files
      Time (PRMS control file) independent:
        * sim_name.disl: the linear discretization (disl) of the stream network
        * sim_name.mmr: the data supporting MMR solution
      Time (PRMS control file) dependend:
        * sim_name.nam: a simulation for just MMR (optional output)
        * sim_name.flw: the boundary conditions to the stream network
        * sim_name.mmr.obs: the location of desired model observations

    Args:
      control_file: the filepath to the PRMS control file. This and control are

      param_file: The filepath to the domain parameters. This or params is
        required.
      params: a PrmsParameters object (already loaded from file). This or
        param_file is required.

      hru_shapefile: a shapefile for HRUS - NOT currently used
      length_units: only meters accepted for now (default)
      disu_file: the output file to write.

    TODOS:
        Do we manage unit conversions? Will mandate meters for now
        Is there a way to test that this is a conforming file in mf6?
        Non-contiguous polygons for HRUs.

    Examples
    --------

    # Ex 1. similar to autotest/test_dis_hru.py
    import pathlib as pl

    from pynhm.constants import __pynhm_root__
    from pynhm.utils import DisHru

    # not used
    # shape_file = (
    #     "/Users/jamesmcc/usgs/data/pynhm/20220209_gm_delaware_river"
    #     "/GIS_simple/HRU_subset.shp")
    param_file = (__pynhm_root__ / "../test_data/drb_2yr/myparam.param")
    disu_file = pl.Path(".") / "disu_example_file.mf6"
    dis = DisHru(param_file=param_file, disu_file=disu_file)

    """

    def __init__(
        self,
        control_file: fileish = None,
        param_file: fileish = None,
        control: Control = None,
        params: PrmsParameters = None,
        hru_shapefile: fileish = None,
        length_units: str = "meters",
        output_dir: fileish = pl.Path("."),
        sim_name: str = "mmr_to_mf6",
        inflow_dir: fileish = None,
        # intial flows over ride from file?
        **kwargs,
    ):

        # read the parameter/control files (similar but not identical logic)
        inputs_dict = {
            "param": {"obj": params, "file": param_file},
            "control": {"obj": control, "file": control_file},
        }

        for key, val_dict in inputs_dict.items():
            obj_file = val_dict["file"]
            obj = val_dict["obj"]

            if (obj_file is not None) and (obj is not None):
                if key == "param":
                    msg = "Can only specify one of param_file or params"
                else:
                    msg = "Can only specify one of control_file or control"
                raise ValueError(msg)

            elif (obj_file is None) and (obj is None):
                if key == "param":
                    msg = "Must specify (exactly) one of param_file or params"
                    raise ValueError(msg)
                else:
                    msg = (
                        "When control is not passed to MMRToMF6, it does note "
                        "create MF6 .nam, .flw, nor .obs files"
                    )
                    setattr(self, "control_file", None)
                    setattr(self, "control", None)
                    warn(msg)

            elif obj_file:
                setattr(self, f"{key}_file", obj_file)
                if key == "param":
                    setattr(self, "params", PrmsParameters.load(obj_file))
                else:
                    setattr(self, "control", Control.load(obj_file))

            else:
                setattr(self, f"{key}_file", None)

        # shorthand
        parameters = self.params.parameters

        # MF6 simulation
        sim = flopy.mf6.MFSimulation(
            sim_ws=output_dir,
            sim_name=sim_name,
            # version="mf6",
            # exe_name="mf6",
            memory_print_option="all",
        )

        # TDIS
        # time requires control file
        if hasattr(self, "control"):
            nper = self.control._n_times
            tdis_rc = [(1.0, 1, 1.0) for ispd in range(nper)]
            _ = flopy.mf6.ModflowTdis(
                sim,
                pname="tdis",
                time_units="DAYS",
                nper=nper,
                perioddata=tdis_rc,
            )

        # EMS
        _ = flopy.mf6.ModflowEms(sim)

        # SNF
        snf = flopy.mf6.ModflowSnf(sim, modelname=sim_name)

        # DISL

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

        nsegment = parameters["nsegment"]
        hru_segment = parameters["hru_segment"] - 1
        tosegment = parameters["tosegment"] - 1

        _ = flopy.mf6.ModflowSnfdisl(
            snf,
            nodes=nsegment,
            nvert=nvert,
            segment_length=parameters["seg_length"],
            tosegment=tosegment,
            idomain=1,  # ??
            vertices=vertices,
            cell2d=cell2d,
        )

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
        outflow_mask = np.full((len(tosegment)), False)
        for iseg in range(nsegment):
            theseg = tosegment[iseg]
            if theseg < 0:
                outflow_mask[iseg] = True
                continue
            connectivity.append((iseg, theseg))

        # use networkx to calculate the Directed Acyclic Graph
        if nsegment > 1:
            graph = nx.DiGraph()
            graph.add_edges_from(connectivity)
            segment_order = list(nx.topological_sort(graph))
        else:
            segment_order = [0]

        segment_order = np.array(segment_order, dtype=int)

        # solve k_coef - taken from the PRMSChannel initialization
        velocity = (
            (
                (1.0 / parameters["mann_n"])
                * np.sqrt(parameters["seg_slope"])
                * parameters["seg_depth"] ** (2.0 / 3.0)
            )
            * 60.0
            * 60.0
        )

        # JLM: This is a bad idea and should throw an error rather than edit
        # inputs in place during run
        # Shouldnt this ALSO BE ABOVE the velocity calculation?
        parameters["seg_slope"] = np.where(
            parameters["seg_slope"] < 1e-7, 0.0001, parameters["seg_slope"]
        )  # not in prms6

        # initialize Kcoef to 24.0 for segments with zero velocities
        # this is different from PRMS, which relied on divide by zero resulting
        # in a value of infinity that when evaluated relative to a maximum
        # desired Kcoef value of 24 would be reset to 24. This approach is
        # equivalent and avoids the occurence of a divide by zero.
        Kcoef = np.full(nsegment, 24.0, dtype=float)

        # only calculate Kcoef for cells with velocities greater than zero
        idx = velocity > zero
        Kcoef[idx] = parameters["seg_length"][idx] / velocity[idx]
        Kcoef = np.where(
            parameters["segment_type"] == SegmentType.LAKE.value, 24.0, Kcoef
        )
        Kcoef = np.where(Kcoef < 0.01, 0.01, Kcoef)
        Kcoef = np.where(Kcoef > 24.0, 24.0, Kcoef)

        _ = flopy.mf6.ModflowSnfmmr(
            snf,
            print_flows=True,
            observations=mmr_obs,
            iseg_order=segment_order,
            qoutflow0=parameters["segment_flow_init"],
            k_coef=Kcoef,
            x_coef=parameters["x_coef"],
        )

        # Boundary conditions
        # aggregate inflows over the contributing fluxes

        inflow_list = ["sroff_vol", "ssres_flow_vol", "gwres_flow_vol"]

        def read_inflow(vv):
            return xr.open_dataset(pl.Path(inflow_dir) / f"{vv}.nc")[vv]

        s_per_time = self.control.time_step_seconds
        inflows = sum([read_inflow(vv) for vv in inflow_list]) / (s_per_time)

        # calculate lateral flow term
        lat_inflow = np.zeros((nper, nsegment))

        for ihru in range(parameters["nhru"]):
            iseg = hru_segment[ihru]
            if iseg < 0:
                # This is bad, selective handling of fluxes is not cool,
                # mass is being discarded in a way that has to be coordinated
                # with other parts of the code.
                # This code shuold be removed evenutally.
                inflows[ihru] = zero
                continue

            lat_inflow[:, iseg] += inflows[:, ihru]

        flw_spd = {
            ispd: [[irch, lat_inflow[ispd, irch]] for irch in range(nsegment)]
            for ispd in range(nper)
        }

        _ = flopy.mf6.ModflowSnfflw(
            snf,
            print_input=True,
            print_flows=True,
            stress_period_data=flw_spd,
        )

        sim.write_simulation()

        return

    # def write(self, out_file: fileish = None):
    #     if out_file:
    #         self.disu_file = out_file

    #     if self.disu_file:
    #         mf6_file_writer(
    #             self,
    #             file_struct=disu_struct,
    #             required=required,
    #             output_file=self.disu_file,
    #         )
    #     else:
    #         print(
    #             "No output file (self.disu_file)  has been set for this "
    #             "DisHru object, no output written."
    #         )

    #     return
