import pathlib as pl
from warnings import warn

import flopy
import numpy as np
import pint

from pywatershed import Control, meta

from ..constants import fileish
from ..parameters import PrmsParameters

# try:
#     import geopandas as gpd

#     has_geopandas = True
# except ModuleNotFoundError:
#     has_geopandas = False


class MMRToMF6:
    """A base class for MMRToMMR, MMRToMCT, and MMRToDFW

    TODO: all the documentation below and for child classes.

    # this is PRMSMMR to SNFDISL, SNFFLW, and/or SNFMMR

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

      segment_shapefile: a shapefile for segments - NOT currently used
      length_units: only meters accepted for now (default)
      disu_file: the output file to write.

    TODOS:
        * Deal with time/sim being optional
        * add shapefile information, deal with non-contiguous HRUs.
        * run mf6 in the tests


    Examples
    --------


    """

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
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        time_zone="UTC",
        write_on_init: bool = True,
        **kwargs,
    ):
        self.output_dir = output_dir
        self.segment_shapefile = segment_shapefile

        self.units = pint.UnitRegistry(system="mks")
        # these are not really optional at the moment, but hope to make so soon
        # making the system definition using input units.
        self.length_units = "meters"
        self.time_units = "seconds"

        # Preliminaries
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
                    setattr(
                        self,
                        "control",
                        Control.load_prms(obj_file, warn_unused_options=False),
                    )

            else:
                setattr(self, f"{key}_file", None)

        # shorthand
        parameters = self.params.parameters

        # MF6 simulation
        self._sim = flopy.mf6.MFSimulation(
            sim_ws=str(self.output_dir),
            sim_name=sim_name,
            # version="mf6",
            # exe_name="mf6",
            memory_print_option="all",
        )

        # TDIS
        # time requires control file

        self.start_time = start_time
        self.end_time = end_time
        if start_time is None:
            self.start_time = self.control.start_time
        if end_time is None:
            self.end_time = self.control.end_time

        self.nper = (
            int((self.end_time - self.start_time) / self.control.time_step) + 1
        )

        self._timestep_s = self.control.time_step_seconds * self.units(
            "seconds"
        )
        # convert to output units in the output data structure
        perlen = (
            self._timestep_s.to_base_units().magnitude
        )  # not reused, ok2convert

        if hasattr(self, "control"):
            tdis_rc = [(perlen, 1, 1.0) for ispd in range(self.nper)]
            _ = flopy.mf6.ModflowTdis(
                self._sim,
                pname="tdis",
                time_units=self.time_units,
                start_date_time=str(self.start_time) + time_zone,
                nper=self.nper,
                perioddata=tdis_rc,
            )

        # SWF
        self._swf = flopy.mf6.ModflowSwf(self._sim, modelname=sim_name)

        # DISL

        # TODO: vertices
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

        self._nsegment = self.params.dims["nsegment"]
        self._hru_segment = parameters["hru_segment"] - 1
        self._tosegment = parameters["tosegment"] - 1

        # united quantities

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

        return

    def write(self):
        print(f"\nWriting simulation files to: {self.output_dir}")
        self._sim.write_simulation()

        return
