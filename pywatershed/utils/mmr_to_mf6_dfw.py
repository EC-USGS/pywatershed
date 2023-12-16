import pathlib as pl

import flopy
import numpy as np

from pywatershed import Control, meta

from ..constants import fileish
from ..parameters import PrmsParameters
from .mmr_to_mf6 import MmrToMf6

# TODO
# * zhb?


class MmrToMf6Dfw(MmrToMf6):
    """Muskingum-Mann Routing to MF6 Diffusive Wave

    MMR: Muskingum-Mann Routing
    DFW: Diffusive Wave

    PRMS MMR Parameters used and how:
    * tosegment:
        Used to give the reach/segment connectivity (transformed to zero based
       for use in python and flopy). No units.
    * seg_length:
        Used for SfwDisl, the discretization of the reaches, directly. Also
        used to calculate seg_mid_elevation if not present in parameters.
        Units of meters from metadata (these are converted to units requested
        for modflow, which we currently hardcode to meters and seconds.)
    * seg_slope:
        Passed directly to SWF DFW. Also used to calculate seg_mid_elevation if
        not present in parameters. See its description below. Unitless slope.
    * hru_segment:
        Used to calculate seg_mid_elevation if not present in parameters. See
        its description below. Unitless index mapping from PRMS HRUs to
        segments.
    * hru_elev:
        Used to calculate seg_mid_elevation if not present in parameters. See
        its description below. Units are in meters (though not documented in
        the metadata as such, just "elev_units")
    * seg_width:
        This is (apparently) the NHD bank-full width present in the PRMS
        parameter files but unused in PRMS. Units are in meters.
    * mann_n:
        Mannings N coefficient for MMR. Passed directly to SWF DFW. Units in
        seconds / meter ** (1/3) according to the metadata.

    Can be user-supplied parameters but are NOT part of PRMS parameters:
    * seg_mid_elevation:
        The elevation at the mid-point of a segment. This is calculated always
        starting from one of the basin outlets. The outlet elevation is
        calculated as per chd_options["stress_period_data"], the downstream
        end of the segment is assumed to have the elevation of the lowest HRU
        which maps to this segment:
        - outlet_downstream_elev = lowest_hru_elev
        The upstream end of the outlet segment is then calculated:
        - outlet_upstream_elev =  lowest_hru_elev + seg_length * seg_slope
        The midpoint of the of the outlet is the average of upstream and
        downstream segment elevations. For all other segments above outlets,
        its downstream elevation is the same as its downstream segment's
        upstream elevation and so
        - segment_upstream_elev =
          downstream_seg_upstream_elev + seg_length * seg_slope
        And its midpoint is again the average elevation between upstream and
        downstream elevations of the segment.
    * chd_options["stress_period_data"]:
        At each outlet in the domain, this is taken as the lowest hru elevation
        over all hrus which map to this segment.

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
        sim_name: str = "mmr_to_dfw",
        inflow_dir: fileish = None,
        inflow_from_PRMS: bool = True,
        # intial flows over ride from file?
        # length_units="meters",
        # time_units="seconds",
        save_flows: bool = True,
        start_time: np.datetime64 = None,
        end_time: np.datetime64 = None,
        time_zone="UTC",
        write_on_init: bool = True,
        dfw_options: dict = None,
        dfw_griddata: dict = None,
        ims_options: dict = None,
        sto_options: dict = None,
        ic_options: dict = None,
        cxs_options: dict = None,
        oc_options: dict = None,
        flw_options: dict = None,
        chd_options: dict = None,
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
            # length_units=length_units,
            # time_units=time_units,
            save_flows=save_flows,
            start_time=start_time,
            end_time=end_time,
            time_zone=time_zone,
            write_on_init=write_on_init,
            dfw_options=dfw_options,
            dfw_griddata=dfw_griddata,
        )

        parameters = self.params.parameters
        # bottom elevations
        if "seg_mid_elevation" in parameters.keys():
            self._seg_mid_elevation = parameters["seg_mid_elevation"]
        else:
            self._seg_mid_elevation = self.calculate_seg_mid_elevations(
                check=False
            )

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

        self._nsegment = self.params.dims["nsegment"]
        self._tosegment = parameters["tosegment"] - 1

        # unit-ed quantities
        segment_units = self.units(meta.parameters["seg_length"]["units"])
        self._segment_length = parameters["seg_length"]
        self._segment_length = self._segment_length * segment_units

        _ = flopy.mf6.ModflowSwfdisl(
            self._swf,
            nodes=self._nsegment,
            nvert=nvert,
            reach_length=self._segment_length.to_base_units().magnitude,  # m
            reach_bottom=self._seg_mid_elevation,  # m
            toreach=self._tosegment,
            idomain=1,  # ??
            vertices=vertices,
            cell2d=cell2d,
            length_units=self.length_units,
        )

        if flw_options is None:
            flw_options = {}

        print("FLW")
        self.set_flw(**flw_options)

        _ = flopy.mf6.ModflowIms(self._sim, **ims_options)

        if "save_flows" not in dfw_options:
            dfw_options["save_flows"] = self._save_flows

        print("DFW")
        _ = flopy.mf6.ModflowSwfdfw(
            self._swf,
            width=parameters["seg_width"],  # this is in meters per metadata
            manningsn=parameters["mann_n"],  # seconds / meter ** (1/3)
            slope=parameters["seg_slope"],  # ratio
            **dfw_options,
        )

        if "save_flows" not in sto_options:
            sto_options["save_flows"] = self._save_flows

        print("STO")
        _ = flopy.mf6.ModflowSwfsto(
            self._swf,
            **sto_options,
        )

        print("IC")
        if ic_options is None:
            ic_options = {}

        _ = flopy.mf6.ModflowSwfic(self._swf, **ic_options)

        print("CXS")
        if cxs_options is None:
            pass
        else:
            _ = flopy.mf6.ModflowSwfcxs(self._swf, **cxs_options)

        print("OC")
        if oc_options is None:
            oc_options = {}

        _ = flopy.mf6.ModflowSwfoc(
            self._swf,
            budget_filerecord=f"{self._sim_name}.bud",
            stage_filerecord=f"{self._sim_name}.stage",
            # qoutflow_filerecord=f"{self._sim_name}.qoutflow",
            **oc_options,
        )

        print("CHD")
        if chd_options is None:
            chd_options = {}
        if "stress_period_data" not in chd_options.keys():
            if hasattr(self, "_outlet_chds"):
                chd_options["stress_period_data"] = list(
                    self._outlet_chds.items()
                )

            else:
                msg = (
                    "CHD without calculation of seg_mid_elevation is "
                    "not yet implemented."
                )
                raise NotImplementedError(msg)
                # This approach just will not be consistent with
                # randomly supplied data for seg_mid_elevation,
                # likely not possible and would have to be passed in.

        chd_options["maxbound"] = len(chd_options["stress_period_data"])
        _ = flopy.mf6.ModflowSwfchd(self._swf, **chd_options)

        return

    def calculate_seg_mid_elevations(self, check=False):
        nseg = self._nsegment
        parameters = self.params.parameters

        # the rise of the reach
        seg_dy = parameters["seg_slope"] * parameters["seg_length"]
        # elevation at upstream end of each reach
        seg_y = seg_dy * np.nan

        # all 1-based indexers for fortran
        tosegment0 = parameters["tosegment"] - 1  # move to zero-based indexing
        is_outflow = -1  # indicates outflow in zero based
        hru_seg = parameters["hru_segment"] - 1

        # not an indexer
        hru_elev = parameters["hru_elev"]

        # We want the elevation at middle of the segment
        # seg_dy is the total rise (y) of each segment.
        # so the elevation at the middle of each segment is:
        # seg_y_mid[ss] = sum(seg_dy[segs_downstream]) + seg_dy[ss] / 2

        # probably best to solve
        # seg_y[ss] = sum(seg_dy[segs_downstream]) + seg_dy[ss]
        # and then solve
        # seg_y_mid = seg_y - seg_dy/2
        # because we want to leverage already calculated segments for
        # efficiency

        # for each segment,
        # find its first downstream reach already solved or the outlet:
        # start_ind
        # this gives the starting datum as: start_seg_y[start_ind] or zero
        # (outlet).
        #  go from outlet back to segment, solving all seg_y

        self._outlet_chds = {}
        for ss in range(nseg):
            already_solved = ~np.isnan(seg_y[ss])
            if already_solved:
                continue
            segment_ind = ss
            downstream_seg_inds = []
            while not already_solved and segment_ind != is_outflow:
                # segment ind only gets added if it is NOT already solved
                downstream_seg_inds += [segment_ind]
                # advance downstream
                segment_ind = tosegment0[segment_ind]
                # check if solved
                already_solved = ~np.isnan(seg_y[segment_ind])

            _ = downstream_seg_inds.reverse()  # reverses in-place

            for ds_seg_ind in downstream_seg_inds:
                my_downstream_ind = tosegment0[ds_seg_ind]
                if my_downstream_ind == is_outflow:
                    # Use the associated hru elevation (minimum if multiple)
                    # as the height of the outlet.
                    # not a great assumption, but better than nothing.
                    # Also save these outflow elevations to specify a
                    # constant head boundary to use later
                    outlet_hrus = np.where(hru_seg == ds_seg_ind)
                    outlet_elev = hru_elev[outlet_hrus].min()
                    self._outlet_chds[ds_seg_ind] = outlet_elev
                    seg_y[ds_seg_ind] = seg_dy[ds_seg_ind] + outlet_elev
                else:
                    seg_y[ds_seg_ind] = (
                        seg_dy[ds_seg_ind] + seg_y[my_downstream_ind]
                    )

        # check?
        # Compare highest segment elevation and highest HRU elevation
        # print(
        #     f"{domain_name}:\n"
        #     "seg_y.max() / parameters['hru_elev'].max() = "
        #     f"{seg_y.max()} / {parameters['hru_elev'].max()} = "
        #     f"{seg_y.max() / parameters['hru_elev'].max()}"
        # )
        # drb_2yr:
        # seg_y.max() / parameters['hru_elev'].max() = 1141.3 / 932.0 = 1.225
        # ucb_2yr:
        # seg_y.max() / parameters['hru_elev'].max() = 3019.86 / 3804.0 = 0.794
        # could correlate hru_to_seg ... ?

        # check going downstream that all differences current - downstream are
        # seg_dy
        if check:
            for ss in range(nseg):
                my_downstream_ind = tosegment0[ss]
                if my_downstream_ind == is_outflow:
                    assert (
                        abs(seg_y[ss] - seg_dy[ss] - self._outlet_chds[ss])
                        < 1.0e-7
                    )
                    # assert abs(seg_y[ss] - seg_dy[ss]) < 1.0e-7
                else:
                    assert (
                        (seg_y[ss] - seg_y[my_downstream_ind]) - seg_dy[ss]
                    ) < 1.0e-7

        self._seg_y_mid = seg_y - (seg_dy / 2)
        return
