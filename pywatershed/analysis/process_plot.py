import pathlib as pl
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon

from ..base import meta
from ..base.model import Model
from ..base.process import Process
from ..utils.optional_import import import_optional_dependency


class ProcessPlot:
    def __init__(
        self,
        gis_dir: pl.Path,
        hru_shp_file_name: str = "HRU_subset.shp",
        seg_shp_file_name: str = "Segments_subset.shp",
    ):
        gpd = import_optional_dependency("geopandas")

        self.hru_shapefile = gis_dir / hru_shp_file_name
        self.seg_shapefile = gis_dir / seg_shp_file_name

        # HRU one-time setups
        self.hru_gdf = (
            gpd.read_file(self.hru_shapefile)
            .drop("nhm_id", axis=1)
            .rename(columns={"nhru_v1_1": "nhm_id"})
            .set_index("nhm_id")
        )

        # segment one-time setup
        self.seg_gdf = gpd.read_file(self.seg_shapefile)
        # if (self.__seg_poly.crs.name
        #     == "USA_Contiguous_Albers_Equal_Area_Conic_USGS_version"):
        #     print("Overriding USGS aea crs with EPSG:5070")
        self.seg_gdf.crs = "EPSG:5070"

        self.seg_geoms_exploded = (
            self.seg_gdf.explode()
            .reset_index(level=1, drop=True)
            .drop("model_idx", axis=1)
            .rename(columns={"nsegment_v": "nhm_seg"})
            .set_index("nhm_seg")
        )

        return

    def plot(self, var_name: str, process: Process, cmap: str = None):
        var_dims = list(meta.get_vars(var_name)[var_name]["dims"])
        if "nsegment" in var_dims:
            if not cmap:
                cmap = "cool"
            return self.plot_seg_var(var_name, process, cmap)
        elif "nhru" in var_dims:
            return self.plot_hru_var(var_name, process)
        else:
            raise ValueError()

    def plot_seg_var(self, var_name: str, process: Process, cmap="cool"):
        ccrs = import_optional_dependency("cartopy.crs")

        data_df = pd.DataFrame(
            {
                "nhm_seg": process.params.coords["nhm_seg"],
                var_name: process[var_name],
            }
        ).set_index("nhm_seg")

        minx, miny, maxx, maxy = self.hru_gdf.geometry.total_bounds
        hru_geoms_exploded = self.hru_gdf.explode().reset_index(
            level=1, drop=True
        )

        aa = {}
        for yy in self.seg_gdf.crs.coordinate_operation.params:
            aa[yy.name] = yy.value
        if "9822" in self.seg_gdf.crs.coordinate_operation.method_code:
            # Albers Equal Area
            crs_proj = ccrs.AlbersEqualArea(
                central_longitude=aa["Longitude of false origin"],
                central_latitude=aa["Latitude of false origin"],
                standard_parallels=(
                    aa["Latitude of 1st standard parallel"],
                    aa["Latitude of 2nd standard parallel"],
                ),
                false_easting=aa["Easting at false origin"],
                false_northing=aa["Northing at false origin"],
            )
        elif "9802" in self.seg_gdf.crs.coordinate_operation.method_code:
            # Lambert Conformal Conic
            crs_proj = ccrs.LambertConformal(
                central_latitude=aa["Latitude of false origin"],
                central_longitude=aa["Longitude of false origin"],
                standard_parallels=(
                    aa["Latitude of 1st standard parallel"],
                    aa["Latitude of 2nd standard parallel"],
                ),
                false_easting=aa["Easting at false origin"],
                false_northing=aa["Northing at false origin"],
            )

        df_plot = self.seg_geoms_exploded.join(data_df)
        norm = Normalize(
            vmin=df_plot[var_name].min().min(),
            vmax=df_plot[var_name].max().max(),
        )

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(30, 20))
        ax = plt.axes(projection=crs_proj)
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([minx, maxx, miny, maxy], crs=crs_proj)

        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        mapper.set_array(df_plot[var_name])

        metadata = meta.get_vars(var_name)[var_name]
        plt.title("Variable: {}".format(var_name))
        plt.colorbar(mapper, shrink=0.6, label=metadata["units"], ax=ax)

        _ = plot_polygon_collection(
            ax,
            hru_geoms_exploded.geometry,
            **dict(cmap=cmap, norm=norm, linewidth=0.5, alpha=0.05),
        )

        _ = plot_line_collection(
            ax,
            df_plot.geometry,
            values=df_plot[var_name],
            **dict(cmap=cmap, norm=norm),
        )

        return

    def get_hru_var(self, var_name: str, model: Model):
        # find the process
        for proc_name, proc in model.processes.items():
            params_vars = list(set(proc.variables) | set(proc.parameters))
            if var_name in params_vars:
                process = proc
                break

        data_df = pd.DataFrame(
            {
                "nhm_id": process.params.coords["nhm_id"],
                var_name: process[var_name],
            }
        ).set_index("nhm_id")
        return data_df

    def plot_hru_var(
        self,
        var_name: str,
        process: Process,
        data: np.ndarray = None,
        data_units: str = None,
        nhm_id: np.ndarray = None,
    ):
        _ = import_optional_dependency("hvplot.pandas")

        ccrs = import_optional_dependency("cartopy.crs")

        if data is None:
            # data_df = self.get_hru_var(var_name, model)
            data_df = pd.DataFrame(
                {
                    "nhm_id": process.params.coords["nhm_id"],
                    var_name: process[var_name],
                }
            ).set_index("nhm_id")

        else:
            if nhm_id is None:
                # nhm_id = model.parameters["nhm_id"]
                raise ValueError("code needs work to handle nhm_id=None")

            data_df = pd.DataFrame(
                {
                    "nhm_id": nhm_id,
                    var_name: data,
                }
            ).set_index("nhm_id")

        plot_df = self.hru_gdf.join(data_df)

        metadata = meta.get_vars(var_name)
        if not len(metadata):
            metadata = meta.get_params(var_name)
        if len(metadata):
            metadata = metadata[var_name]
        else:
            metadata = None

        frame_height = 550
        title = f'"{var_name}"\n'
        clabel = data_units
        if metadata is not None:
            title += "\n".join(
                wrap(
                    f"{metadata['desc']}, {metadata['units']}",
                    width=frame_height / 10,
                )
            )
            clabel = f'{metadata["units"]}'

        plot = plot_df.hvplot(
            tiles=True,
            crs=ccrs.epsg(5070),
            frame_height=frame_height,
            c=var_name,
            line_width=0,
            alpha=0.75,
            # clim=stat_lims[stat],
            hover_cols=["nhm_id"],
            title=title,
            clabel=clabel,
            xlabel="Longitude (degrees East)",
            ylabel="Latitude (degrees North)",
        )
        return plot


def plot_line_collection(
    ax,
    geoms,
    values=None,
    cmap=None,
    norm=None,
    vary_width=False,
    vary_color=True,
    colors=None,
    alpha=1.0,
    linewidth=1.0,
    **kwargs,
):
    """Plot a collection of line geometries"""
    shapely = import_optional_dependency("shapely")

    lines = []
    for geom in geoms:
        a = np.asarray(geom.coords)

        if geom.has_z:
            a = shapely.geometry.LineString(zip(*geom.xy))

        lines.append(shapely.geometry.LineString(a))

    if vary_width:
        lwidths = ((values / values.max()).to_numpy() + 0.01) * linewidth
        if vary_color:
            lines = LineCollection(
                lines, linewidths=lwidths, cmap=cmap, norm=norm, alpha=alpha
            )
        else:
            lines = LineCollection(
                lines, linewidths=lwidths, colors=colors, alpha=alpha
            )
    elif vary_color:
        lines = LineCollection(
            lines, linewidth=linewidth, alpha=alpha, cmap=cmap, norm=norm
        )

    if vary_color and values is not None:
        lines.set_array(values)
        # lines.set_cmap(cmap)

    ax.add_collection(lines, autolim=True)
    ax.autoscale_view()
    return lines


def plot_polygon_collection(
    ax,
    geoms,
    values=None,
    cmap=None,
    norm=None,
    facecolor=None,
    edgecolor=None,
    alpha=1.0,
    linewidth=1.0,
    **kwargs,
):
    """Plot a collection of Polygon geometries"""
    # from https://stackoverflow.com/questions/33714050/
    #     geopandas-plotting-any-way-to-speed-things-up
    shapely = import_optional_dependency("shapely")

    patches = []

    for poly in geoms:
        a = np.asarray(poly.exterior)
        if poly.has_z:
            a = shapely.geometry.Polygon(zip(*poly.exterior.xy))

        patches.append(Polygon(a))

    patches = PatchCollection(
        patches,
        facecolor=facecolor,
        linewidth=linewidth,
        edgecolor=edgecolor,
        alpha=alpha,
        cmap=cmap,
        norm=norm,
    )

    if values is not None:
        patches.set_array(values)
        # patches.set_cmap(cmap)

    ax.add_collection(patches, autolim=True)
    ax.autoscale_view()
    return patches
