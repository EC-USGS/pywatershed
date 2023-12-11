import pathlib as pl
from warnings import warn

import flopy
import numpy as np
import pint
import xarray as xr

from pywatershed import Control, meta

from ..constants import fileish, zero
from ..parameters import PrmsParameters

# try:
#     import geopandas as gpd

#     has_geopandas = True
# except ModuleNotFoundError:
#     has_geopandas = False


class MmrToMf6:
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
        save_flows: bool = False,
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        time_zone="UTC",
        **kwargs,
    ):
        self._written = False
        self._sim_name = sim_name
        self._output_dir = pl.Path(output_dir)
        self.segment_shapefile = segment_shapefile
        self._save_flows = save_flows

        self._inflow_from_PRMS = inflow_from_PRMS
        self._inflow_dir = inflow_dir
        self._bc_flows_combine = bc_flows_combine
        self._bc_binary_files = bc_binary_files

        self.units = pint.UnitRegistry(system="mks")
        # these are not really optional at the moment, but hope to make so soon
        # making the system definition using input units.
        self.length_units = "meters"
        self.time_units = "seconds"

        # Preliminaries
        # read the parameter/control files (similar but not identical logic)
        inputs_dict = {
            "params": {"obj": params, "file": param_file},
            "control": {"obj": control, "file": control_file},
        }

        for key, val_dict in inputs_dict.items():
            obj_file = val_dict["file"]
            obj = val_dict["obj"]

            if (obj_file is not None) and (obj is not None):
                if key == "params":
                    msg = "Can only specify one of param_file or params"
                else:
                    msg = "Can only specify one of control_file or control"
                raise ValueError(msg)

            elif (obj_file is None) and (obj is None):
                if key == "params":
                    msg = "Must specify (exactly) one of param_file or params"
                    raise ValueError(msg)
                else:
                    msg = (
                        "When control is not passed to MMRToMF6, it does not "
                        "create MF6 .nam, .flw, nor .obs files"
                    )
                    setattr(self, "control_file", None)
                    setattr(self, "control", None)
                    warn(msg)

            elif obj_file:
                setattr(self, f"{key}_file", obj_file)
                if key == "params":
                    setattr(self, "params", PrmsParameters.load(obj_file))
                else:
                    setattr(
                        self,
                        "control",
                        Control.load_prms(obj_file, warn_unused_options=False),
                    )

            else:
                setattr(self, f"{key}_file", None)
                setattr(self, key, obj)

        # dimension
        self._nsegment = self.params.dims["nsegment"]
        self._hru_segment = self.params.parameters["hru_segment"] - 1

        # MF6 simulation
        self._sim = flopy.mf6.MFSimulation(
            sim_ws=str(self._output_dir),
            sim_name=self._sim_name,
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
        self._swf = flopy.mf6.ModflowSwf(
            self._sim, modelname=self._sim_name, save_flows=self._save_flows
        )

        return

    def set_flw(self, **kwargs):
        # Boundary conditions / FLW
        # aggregate inflows over the contributing fluxes

        # For non-binary data, this method has to be called after
        # flopy.mf6.Swfdisl is set on self._swf, that's why its a method

        parameters = self.params.parameters

        if self._inflow_from_PRMS:
            inflow_list = ["sroff", "ssres_flow", "gwres_flow"]
        else:
            inflow_list = ["sroff_vol", "ssres_flow_vol", "gwres_flow_vol"]

        # check they all have the same units before summing
        inflow_units = list(meta.get_units(inflow_list, to_pint=True).values())
        inflow_unit = inflow_units[0]
        assert [inflow_unit] * len(inflow_units) == inflow_units
        inflow_unit = self.units(inflow_unit)

        def read_inflow(vv, start_time, end_time):
            ff = xr.open_dataset(pl.Path(self._inflow_dir) / f"{vv}.nc")[vv]
            return ff.sel(time=slice(start_time, end_time))

        inflows = {
            vv: read_inflow(vv, self.start_time, self.end_time)
            for vv in inflow_list
        }

        if self._bc_flows_combine:
            inflows["combined"] = sum(inflows.values())
            for kk in list(inflows.keys()):
                if kk != "combined":
                    del inflows[kk]

        # add the units
        inflows = {kk: vv.values * inflow_unit for kk, vv in inflows.items()}

        # flopy adds one to the index, but if we write binary we have to do it
        add_one = int(self._bc_binary_files)

        for flow_name in inflows.keys():
            # if from pywatershed, inflows are already volumes in cubicfeet
            if "inch" in str(inflow_unit):  # PRMS style need hru areas
                hru_area_unit = self.units(
                    list(meta.get_units("hru_area").values())[0]
                )
                hru_area = parameters["hru_area"] * hru_area_unit
                inflows[flow_name] *= hru_area

            inflows[flow_name] /= self._timestep_s.to("seconds")

            new_inflow_unit = inflows[flow_name].units

            # calculate lateral flow term to the REACH/segment from HRUs
            lat_inflow = (
                np.zeros((self.nper, self._nsegment)) * new_inflow_unit
            )

            for ihru in range(self.params.dims["nhru"]):
                iseg = self._hru_segment[ihru]
                if iseg < 0:
                    # This is bad, selective handling of fluxes is not cool,
                    # mass is being discarded in a way that has to be
                    # coordinated
                    # with other parts of the code.
                    # This code shuold be removed evenutally.
                    inflows[flow_name][:, ihru] = zero * new_inflow_unit
                    continue

                lat_inflow[:, iseg] += inflows[flow_name][:, ihru]

            # convert to output units in the output data structure
            flw_spd = {}
            for ispd in range(self.nper):
                flw_ispd = [
                    (
                        irch + add_one,
                        lat_inflow[ispd, irch].to_base_units().magnitude,
                    )
                    for irch in range(self._nsegment)
                ]

                if self._bc_binary_files:
                    # should put the time in the file names?
                    ra = np.array(
                        flw_ispd, dtype=[("irch", "<i4"), ("q", "<f8")]
                    )
                    i_time_str = str(
                        self.control.start_time + ispd * self.control.time_step
                    )
                    bin_name = (
                        f"swf_flw_bc/"
                        f"flw_{flow_name}_{i_time_str.replace(':', '_')}.bin"
                    )
                    bin_name_pl = self._output_dir / bin_name
                    if not bin_name_pl.parent.exists():
                        bin_name_pl.parent.mkdir()
                    _ = ra.tofile(bin_name_pl)

                    flw_spd[ispd] = {
                        "filename": str(bin_name),
                        "binary": True,
                        "data": None,
                    }

                else:
                    flw_spd[ispd] = flw_ispd

            for key in ["print_input", "print_flows"]:
                if key not in kwargs.keys():
                    kwargs[key] = True

            _ = flopy.mf6.ModflowSwfflw(
                self._swf,
                save_flows=self._save_flows,
                stress_period_data=flw_spd,
                maxbound=self._nsegment + 1,
                pname=flow_name,
                **kwargs,
            )

        return

    def write(self, *args, rewrite=False, **kwargs):
        if self._written and not rewrite:
            msg = (
                "The simulation was already written, "
                "use rewrite=True to force."
            )
            warn(msg)
            return

        print(f"\nWriting simulation files to: {self._output_dir}")
        self._written = True
        return self._sim.write_simulation(*args, **kwargs)

    def run(self, *args, **kwargs):
        print(f"\nRunning simulation files in: {self._output_dir}")
        return self._sim.run_simulation(*args, **kwargs)
