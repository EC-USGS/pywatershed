import numpy as np

# try:
#     import geopandas as gpd

#     has_geopandas = True
# except ModuleNotFoundError:
#     has_geopandas = False


from ..constants import fileish
from .mf6_file_writer import mf6_file_writer
from .parameters import PrmsParameters

acres_to_m2 = 4046.8564224
idt = " " * 2  # indent is 2 spaces in output file

required = [
    "nodes",
    "nja",
    "top",
    "bot",
    "area",
    "iac",
    "ja",
    "ihc",
    "cl12",
    "hwva",
]

disu_struct = {
    "options": {
        "length_units": "scalar",
        "NOGRB": "scalar",
        "xorigin": "scalar",
        "yorigin": "scalar",
        "angrot": "scalar",
        "vertical_offset_tolerance": "scalar",
    },
    "dimensions": {
        "nodes": "scalar",
        "nja": "scalar",
        "nvert": "scalar",
    },
    "griddata": {
        "top": "vector",
        "bot": "vector",
        "area": "vector",
        "idomain": "vector",
    },
    "connectiondata": {
        "iac": "vector",
        "ja": "vector",
        "ihc": "vector",
        "cl12": "vector",
        "hwva": "vector",
        "xangldegx": "vector",
    },
    "vertices": {
        "iv": "vector",
        "xv": "vector",
        "yv": "vector",
    },
    "cell2d": {
        "icell2d": "vector",
        "xc": "vector",
        "yc": "vector",
        "ncvert": "vector",
        "icvert": "vector",
    },
}


class DisHru:
    """Write the HRU spatial discretization (no cascades) to an MF6 disu file

    The MF6 IO descrption can be found here
    https://water.usgs.gov/water-resources/software/MODFLOW-6/mf6io_6.3.0.pdf
    The disu section is (currently) found on page 37.

    Args:
      param_file: The filepath to the domain parameters.
      params: a PrmsParameters object (already loaded from file)
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
        param_file: fileish = None,
        params: PrmsParameters = None,
        hru_shapefile: fileish = None,
        length_units: str = "meters",
        disu_file: fileish = None,
        **kwargs,
    ):

        # read the parameter file: currently just for HRU areas
        self._param_file = param_file
        self.params = params
        if (param_file is not None) and (params is not None):
            msg = "Can only specify one of param_file or params"
            raise ValueError(msg)
        elif (param_file is None) and (params is None):
            msg = "Must specify (exactly) one of param_file or params"
            raise ValueError(msg)
        elif param_file:
            self.params = PrmsParameters.load(param_file)

        # which of these are @properties ?

        # OPTIONS block
        if length_units.lower() != "meters":
            msg = "Only currently supporting meters"
            raise ValueError(msg)
        self.length_units = length_units.upper()
        # optional: NOGRB, xorigin, yorigin, angrot, vertical_offset_tolerance

        # DIMENSIONS block
        self.nodes = len(self.params.parameters["nhm_id"])
        self.nja = self.nodes
        # optional: nvert

        # GRIDDATA block
        self.top = 2  # 2 m? why not
        self.bot = 0
        self.area = self.params.parameters["hru_area"] * acres_to_m2
        # optional: idomain

        # CONNECTIONDATA block
        self.iac = np.ones(self.nja, dtype=np.int32)
        self.ja = np.arange(self.nja, dtype=np.int32)
        # according to the examples, ihc can just be constant
        # the shape of ihc, cli12 and hwva all correspond to the shape of ja?
        # it could be more clear in the docs, just by saying so.
        self.ihc = np.zeros(self.nja, dtype=np.int32)  # vertical only
        self.cl12 = np.zeros(self.nja, dtype=np.int32)
        self.hwva = self.area  # for vertical
        # optional: angldegx

        # VERTICES block
        # optional: iv, xv, yv
        # The shapefile information is not required
        # (Do not use area from these polygons)
        # Difficulty is that HRU polygons are not necessarily contiguous
        # 14% of the DRB HRUs are not contiguous (see notebook in examples)
        # Comment this for now.
        # if hru_shapefile:
        #     self._hru_shapefile = hru_shapefile
        #     self.hru_gdf = (
        #         gpd.read_file(self._hru_shapefile)
        #         .rename(columns={"nhru_v1_1": "nhm_id"})
        #         .set_index("nhm_id")
        #     )
        #     # need to sort the gdf to match the params?
        #     # implement if needed
        #     assert (
        #         self.hru_gdf.reset_index()["nhm_id"]
        #         == self.params.parameters["nhm_id"]
        #     ).all()

        # CELL2D block
        # optional: icell2d, xc, yc, ncvert, icvert

        # write out if file is specified
        self.disu_file = disu_file
        if self.disu_file:
            self.write()

        return

    def write(self, out_file: fileish = None):
        if out_file:
            self.disu_file = out_file

        if self.disu_file:
            mf6_file_writer(
                self,
                file_struct=disu_struct,
                required=required,
                output_file=self.disu_file,
            )
        else:
            print(
                "No output file (self.disu_file)  has been set for this "
                "DisHru object, no output written."
            )

        return
