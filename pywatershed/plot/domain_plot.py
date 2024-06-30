# plot with folium
# take only GIS files on init,
# add methods for initialize/finalized model? model yaml?

import pathlib as pl
from typing import Union

import geopandas as gpd
import pandas as pd

from ..base import Parameters
from ..utils import import_optional_dependency

folium = import_optional_dependency("folium", errors="warn")
IPython = import_optional_dependency("IPython", errors="warn")


class DomainPlot:
    """Plot PRMS domain information using the folium package

    Args:
        hru_shp_file: Optional shapefile for the HRUs of the domain.
        segment_shp_file: Optional shapefile for the stream segments of the
            domain.
        hru_parameters: Optional parameter file or Parameters object for HRU
            parameters to display on tooltip/hover. Must also supply
            hru_parameter_names to select the variables names to show.
        segment_parameters: Optional parameter file or Parameters object for
            segment parmaeters to isplay on tooptip/hover. Must also supply
            segment_parameter_names to select the variables to show.
        hru_parameter_names: A list of string names to display with tooltip on
            each HRU.
        segment_parameter_names: A list of string names to display with
            tooltip on each segment.
        hru_highlight_indices: Indices (zero-based) of HRUs to highlight.
        segment_highlight_indices: Indices (zero-based) of segments to
            highlight.
        crs: int = 4326,
        start_lat: float = None,
        start_lon: float = None,
        start_zoom: int = 7,
        dislpay_plot: bool = True,
        return_plot: bool = False,


    """

    def __init__(
        self,
        hru_shp_file: Union[pl.Path, str] = None,
        segment_shp_file: Union[pl.Path, str] = None,
        hru_parameters: Union[pl.Path, Parameters] = None,
        segment_parameters: Union[pl.Path, Parameters] = None,
        hru_parameter_names: list[str] = None,
        segment_parameter_names: list[str] = None,
        hru_highlight_indices: list = None,
        segment_highlight_indices: list = None,
        crs: int = 4326,
        start_lat: float = None,
        start_lon: float = None,
        start_zoom: int = 7,
        display_plot: bool = True,
        return_plot: bool = False,
    ):
        self._hru_shp_file = hru_shp_file
        self._segment_shp_file = segment_shp_file
        # these are set by the add_parameter methods
        # self._hru_parameters = hru_parameters
        # self._segment_parameters = segment_parameters
        # self._hru_parameter_names = hru_parameter_names
        # self._segment_parameter_names = segment_parameter_names
        self._hru_highlight_indices = hru_highlight_indices
        self._segment_highlight_indices = segment_highlight_indices
        self._hru_gdf = None
        self._segment_gdf = None

        self.crs = crs
        self.start_lat = start_lat
        self.start_lon = start_lon
        self.start_zoom = start_zoom
        self.display_plot = display_plot
        self.return_plot = return_plot

        self.dom_plot = None

        if self._hru_shp_file is not None:
            self.add_hru_gdf(self._hru_shp_file)

        if self._segment_shp_file is not None:
            self.add_segment_gdf(self._segment_shp_file)

        if self._hru_gdf is not None:
            if hru_parameters is not None:
                self.add_hru_parameters(
                    hru_parameters,
                    hru_parameter_names,
                )
            else:
                self._set_default_hru_parameter_names()

        # <<<<
        if self._segment_gdf is not None:
            if segment_parameters is not None:
                self.add_segment_parameters(
                    segment_parameters,
                    segment_parameter_names,
                )
            else:
                self._set_default_segment_parameter_names()

        # could add the topo and hydro layers to arguments
        self.topo_layer = folium.TileLayer(
            tiles=(
                "https://basemap.nationalmap.gov/arcgis/rest/services/"
                "USGSTopo/MapServer/tile/{z}/{y}/{x}"
            ),
            attr="USGStopo",
            zoom_start=start_zoom,
            name="USGSTopo",
            overlay=False,
            show=True,
            control=True,
        )

        self.satellite_layer = folium.TileLayer(
            tiles=(
                "https://server.arcgisonline.com/ArcGIS/rest/services/"
                "World_Imagery/MapServer/tile/{z}/{y}/{x}"
            ),
            attr="EsriSatellite",
            name="EsriSatellite",
            overlay=True,
            show=False,
            control=True,
        )

        self.hydro_layer = folium.TileLayer(
            tiles=(
                "https://basemap.nationalmap.gov/arcgis/rest/services/"
                "USGSHydroCached/MapServer/tile/{z}/{y}/{x}"
            ),
            attr="USGSHydroCached",
            zoom_start=start_zoom,
            name="USGSHydroCached",
            overlay=True,
            control=True,
            show=False,
        )

        if self.display_plot:
            self.display()

        if self.return_plot:
            if self.dom_plot is None:
                self.make_domain_plot()
            return self.dom_plot
        else:
            return

    def display(self):
        if self.dom_plot is None:
            self.make_domain_plot()
        IPython.display.display(self.dom_plot)

    def add_hru_gdf(self, gdf):
        if self.dom_plot is not None:
            self.dom_plot = None
        if isinstance(gdf, gpd.geodataframe.GeoDataFrame):
            hru_gdf = gdf
        else:
            # assume it's a file handled by gpd.read_file
            self._hru_shp_file = gdf
            hru_gdf = gpd.read_file(gdf).fillna(0)

        # have to handle at least 2 conventions adopted over the years
        # the original naming convention and variables are absurd.
        rename_cols = {"model_hru_": "model_idx"}
        for old_col, new_col in rename_cols.items():
            if old_col in hru_gdf.columns:
                hru_gdf = hru_gdf.rename(columns={old_col: new_col})

        # this seems so unnecessary
        if "nhm_id" in hru_gdf.columns and "nhru_v1_1" in hru_gdf.columns:
            hru_gdf["nhm_id"] = hru_gdf["nhm_id"] - 1

        drop_cols = ["model_hru_", "nhru_v1_1"]
        for col in drop_cols:
            if col in hru_gdf.columns:
                hru_gdf.drop(columns=col, inplace=True)

        self._hru_gdf = hru_gdf.to_crs(self.crs)
        self._hru_parameters = None
        self._hru_parameter_names = None
        self._set_default_hru_parameter_names()

    def add_segment_gdf(self, gdf):
        if self.dom_plot is not None:
            self.dom_plot = None
        if isinstance(gdf, gpd.geodataframe.GeoDataFrame):
            segment_gdf = gdf
        else:
            # assume it's a file handled by gpd.read_file
            self._segment_shp_file = gdf
            segment_gdf = gpd.read_file(gdf).fillna(0)

        # have to handle at least 2 conventions adopted over the years
        rename_cols = {"model_seg_": "model_idx"}
        for old_col, new_col in rename_cols.items():
            if old_col in segment_gdf.columns:
                segment_gdf = segment_gdf.rename(columns={old_col: new_col})

        drop_cols = ["model_segment_"]
        for col in drop_cols:
            if col in segment_gdf.columns:
                segment_gdf.drop(columns=col, inplace=True)

        rename_cols = {"nsegment_v": "nhm_seg"}
        for old_col, new_col in rename_cols.items():
            if old_col in segment_gdf.columns:
                segment_gdf = segment_gdf.rename(columns={old_col: new_col})

        self._segment_gdf = segment_gdf.to_crs(self.crs)
        self._set_default_segment_parameter_names()

    def add_hru_parameters(self, parameters, parameter_names):
        self._hru_parameters = parameters
        self._hru_parameter_names = parameter_names
        self._hru_gdf = self._add_parameters(
            self._hru_gdf,
            self._hru_parameters,
            self._hru_parameter_names,
        )
        if self.dom_plot is not None:
            self.dom_plot = None

    def add_segment_parameters(self, parameters, parameter_names):
        self._segment_parameters = parameters
        self._segment_parameter_names = parameter_names
        self._segment_gdf = self._add_parameters(
            self._segment_gdf,
            self._segment_parameters,
            self._segment_parameter_names,
        )
        if self.dom_plot is not None:
            self.dom_plot = None

    def _set_default_hru_parameter_names(self):
        for coord_name in ["nhm_id", "nhru_v"]:
            if coord_name in self._hru_gdf.columns:
                self._hru_parameter_names = [coord_name]

    def _set_default_segment_parameter_names(self):
        for coord_name in ["nhm_seg", "nsegment_v"]:
            if coord_name in self._segment_gdf.columns:
                self._segment_parameter_names = [coord_name]

    def _get_basemap(self):
        if self.start_lat is None or self.start_lon is None:
            location = None
        else:
            location = [self.start_lat, self.start_lon]

        basemap = folium.Map(
            location=location,
            # width=self.width,  # width and height disable LayerControl
            # height=self.height,
            zoom_start=self.start_zoom,
        )
        return basemap

    def make_domain_plot(self):
        # add parameters first because that delets the self.dom_plot

        # <<<<
        self.dom_plot = self._get_basemap()
        self.topo_layer.add_to(self.dom_plot)
        self.hydro_layer.add_to(self.dom_plot)
        self.satellite_layer.add_to(self.dom_plot)

        if self._hru_gdf is not None:
            hru_map = folium.GeoJson(
                self._hru_gdf,
                style_function=self.style_function_hru,
                name="HRUs",
                control=True,
            ).add_to(self.dom_plot)
            tooltip_hru = folium.GeoJsonPopup(
                fields=self._hru_parameter_names, labels=True
            )
            hru_map.add_child(tooltip_hru)

            # TODO: highlight HRUs
        # <
        if self._segment_gdf is not None:
            segment_map = folium.GeoJson(
                self._segment_gdf,
                style_function=self.style_function_segments,
                name="Segments",
                control=True,
            ).add_to(self.dom_plot)
            tooltip_segments = folium.GeoJsonTooltip(
                fields=self._segment_parameter_names, labels=True
            )
            segment_map.add_child(tooltip_segments)

            if self._segment_highlight_indices is not None:
                self.segment_highlight_gdf = self._segment_gdf[
                    self._segment_gdf.index.isin(
                        self._segment_highlight_indices
                    )
                ]
                segment_highlight_map = folium.GeoJson(
                    self.segment_highlight_gdf,
                    style_function=self.style_function_segment_highlights,
                    name="Segment Highlights",
                    control=True,
                ).add_to(self.dom_plot)
                tooltip_segment_highlights = folium.GeoJsonTooltip(
                    fields=self._segment_parameter_names, labels=True
                )
                segment_highlight_map.add_child(tooltip_segment_highlights)

        # <<
        folium.LayerControl().add_to(self.dom_plot)
        return

    @staticmethod
    def _add_parameters(
        gdf: gpd.geodataframe.GeoDataFrame,
        parameters: Union[Parameters, pl.Path, str],
        parameter_names: list,
    ):
        if isinstance(parameters, (pl.Path, str)):
            parameters = Parameters.from_netcdf(parameters, use_xr=True)

        if parameter_names is None or len(parameter_names) == 0:
            raise ValueError("No parameter names specified for tooltip.")
        # if coord_name not in parameter_names:
        #     parameter_names += [coord_name]
        param_dict = {
            var: parameters.parameters[var] for var in parameter_names
        }
        params_df = pd.DataFrame.from_dict(param_dict)
        params_df.reset_index(inplace=True)
        params_df["model_idx"] = params_df.index + 1

        join_cols = []
        for col in params_df.columns:
            if col in gdf.columns:
                join_cols += [col]

        gdf = gpd.GeoDataFrame(
            pd.merge(params_df, gdf, on=join_cols),
            geometry="geometry",
        )

        return gdf

    @staticmethod
    def style_function_hru(xx):
        return {
            "fillColor": "gray",
            "color": "black",
            "weight": 2,
            "fill_opacity": 0.15,
        }

    @staticmethod
    def style_function_segments(xx):
        return {
            "color": "cyan",
            "weight": 2,
        }

    @staticmethod
    def style_function_segment_highlights(xx):
        return {
            "color": "orange",
            "weight": 14,
        }
