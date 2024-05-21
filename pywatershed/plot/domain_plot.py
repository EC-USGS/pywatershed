# plot with folium
# take only GIS files on init,
# add methods for initialize/finalized model? model yaml?

import pathlib as pl
from typing import Union

import folium
import geopandas as gpd
import pandas as pd
from IPython.display import display

import pywatershed as pws


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
        hru_parameters: Union[pl.Path, pws.Parameters] = None,
        segment_parameters: Union[pl.Path, pws.Parameters] = None,
        hru_parameter_names: list[str] = None,
        segment_parameter_names: list[str] = None,
        crs: int = 4326,
        start_lat: float = None,
        start_lon: float = None,
        start_zoom: int = 7,
        display_plot: bool = True,
        return_plot: bool = False,
    ):
        # TODO: make segments optional
        # TODO: make HRUs optional

        self._hru_shp_file = hru_shp_file
        self._segment_shp_file = segment_shp_file
        self._hru_parameters = hru_parameters
        self._segment_parameters = segment_parameters
        self._hru_parameter_names = hru_parameter_names
        self._segment_parameter_names = segment_parameter_names

        self.crs = crs
        self.start_lat = start_lat
        self.start_lon = start_lon
        self.start_zoom = start_zoom
        self.display_plot = display_plot
        self.return_plot = return_plot

        if self._hru_shp_file is not None:
            hru_gdf = gpd.read_file(hru_shp_file).fillna(0)
            hru_gdf.drop(columns="model_hru_", inplace=True)
            self.hru_gdf = hru_gdf.to_crs(self.crs)

        if self._segment_shp_file is not None:
            seg_gdf = gpd.read_file(segment_shp_file).fillna(0)
            seg_gdf.drop(columns="model_seg_", inplace=True)
            seg_gdf = seg_gdf.to_crs(self.hru_gdf.crs)
            self.segment_gdf = seg_gdf.to_crs(self.crs)

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

        self.make_domain_plot()

        if self.display_plot:
            display(self.dom_plot)

        if self.return_plot:
            return self.dom_plot
        else:
            return

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
        self.dom_plot = self._get_basemap()
        self.topo_layer.add_to(self.dom_plot)
        self.hydro_layer.add_to(self.dom_plot)
        self.satellite_layer.add_to(self.dom_plot)

        if self._hru_shp_file:
            if self._hru_parameters:
                self._add_parameters(self._hru_parameters)
            else:
                self._hru_parameter_names = ["nhm_id"]
            # <
            hru_map = folium.GeoJson(
                self.hru_gdf,
                style_function=self.style_function_hru,
                name="HRUs",
                control=True,
            ).add_to(self.dom_plot)
            tooltip_hru = folium.GeoJsonPopup(
                fields=self._hru_parameter_names, labels=True
            )
            hru_map.add_child(tooltip_hru)

        if self._segment_shp_file:
            if self._segment_parameters:
                self._add_parameters(self._segment_parameters)
            else:
                self._segment_parameter_names = ["nhm_seg"]

            seg_map = folium.GeoJson(
                self.segment_gdf,
                style_function=self.style_function_segments,
                name="Segments",
                control=True,
            ).add_to(self.dom_plot)
            tooltip_seg = folium.GeoJsonTooltip(
                fields=self._segment_parameter_names, labels=True
            )
            seg_map.add_child(tooltip_seg)

        # <
        folium.LayerControl().add_to(self.dom_plot)
        return

    def _add_parameters(self, parameters: Union[pws.Parameters, pl.Path, str]):
        if isinstance(parameters, (pl.Path, str)):
            parameters = pws.Parameters.from_netcdf(parameters, use_xr=True)

        if "nhm_id" in parameters.parameters.keys():
            coord_name = "nhm_id"
            gdf = self.hru_gdf
            parameter_names = self._hru_parameter_names
        elif "nhm_seg" in parameters.parameters.keys():
            coord_name = "nhm_seg"
            gdf = self.segment_gdf
            parameter_names = self._segment_parameter_names
        else:
            msg = "Coordinate names not found in parameter data"
            raise ValueError(msg)

        if parameter_names is None:
            raise ValueError("No parameter names specified for tooltip.")
        if coord_name not in parameter_names:
            parameter_names += [coord_name]
        param_dict = {
            var: parameters.parameters[var] for var in parameter_names
        }
        params_df = pd.DataFrame.from_dict(param_dict)
        params_df.reset_index(inplace=True)
        # necessary? just an index
        params_df["model_idx"] = params_df.index + 1

        gdf = gpd.GeoDataFrame(
            pd.merge(params_df, gdf, on=coord_name),
            geometry="geometry",
        )
        # gdf.set_index(coord_name, inplace=True, drop=True)

        if coord_name == "nhm_id":
            self.hru_gdf = gdf
        else:
            self.segment_gdf = gdf

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
