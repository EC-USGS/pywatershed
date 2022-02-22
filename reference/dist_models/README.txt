Please Note:

These examples are for illustrative and educational purposes only
to demonstrate various simulation options of PRMS. They are not
verified or documented, and are in no way to be used as a representation
of the water resources in the subject watersheds.

                            PRMS Sample Problems
                                 November 2021
 
Seven sample projects with input data sets are provided in this subdirectory
to verify that PRMS is correctly installed and running on the user's system.
The sample projects also may be looked at as examples of how to use various 
PRMS simulation, input, and ouput options. The Data and Parameter Files can 
be found in the 'input' subdirectory. Files ending in 'day' refer to 
pre-processed climate-data files used with the climate_by_hru module. The 
Control File(s) can be found in the 'control' subdirectory. Results for
simulations can be found in the 'output' subdirectory. Results run from
the development computer for simulations can be found in the 'output-test' 
subdirectory and are intended for comparison purposes of results in the 
'output' subdirectories. Note, the streamflow routing module is specified as
muskingum for the three Apalachicola-Chattahoochee-Flint River Basin projects
as this is a large watershed where the travel time of flow is typically
greater than a day. The strmflow routing module for the other examples is
specified as the the simple method (module strmflow) that assumes all flow 
entering the stream network leaves the basin in a single day as they are 
small watersheds.

Double-click on one nogui.bat or gui.bat files to run the model or paramtool.bat 
to bring up a spread-sheet like interface to display and edit the primary 
Parameter File. Notes, a) as some projects have multiple script files, the 
names may differ to indicate each example within a project, e.g., nogui_IDE.bat 
and nogui_XYZ.bat; b) after running the gui.bat script, users may need to enter 
"control C" to end the process. 

Sample problems:

1. sagehen: This sample model is for the Sagehen Creek Watershed and 
   is described in the original GSFLOW documentation (Markstrom and others,
   2008, USGS TM 6-D1). This sample includes script files to run three simulations:
   1) using modules temp_1sta and precip_1sta to distribute temperature and
   preciptiation, respectively, from observed station data, 2) using module
   climate_by_hru that reads predistributed climate values, and 3) using module
   map_results to write mapped output for possible use to loosely couple
   PRMS to a grid-based model requiring recharge values as input. Various output
   options are selected for each simulation. The new control parameter
   soilzone_aet_flag was set to 1, which computes soil-water ET based on
   the potential evapotranspiration (PET) rate limited by the unsatisfied PET
   rate to be consistent with other components of the actual evapotransiration
   (AET) rate can be higher on some timesteps, thus different than previous
   versions of PRMS. Note: results from the Map Results module using and input
   file (gvr.params) that specifies the intersections between the HRU map and
   output gridded map will be written to the current directory.

2. merced: This sample model is for the Merced River at Pohono Bridge near 
   Yosemite National Park, California and is described in the PRMS-IV 
   documentation (Markstrom and others, 2015, USGS TM 6-B7). This sample 
   includes script files to run one simulation using the ide_dist module
   and one using the xyz_dist module for precipation and temperature
   distribution.
   
3. acf: This sample model is for the Apalachicola-Chattahoochee-Flint 
   River Basin the watershed above the Chattahoochee River near Norcross, 
   Georgia and is described in LaFontaine and others (2013). This sample 
   illustrates running PRMS to produce output from the nhru_summary module 
   that is described in Regan and LaFontaine (2017), a statvar file that is
   described in Markstro and others (2015), and from the prms_summary and 
   basin_summary modules that are described in the Release notes. 
   
4. acfb_dyn_params: This sample model is for the Apalachicola-Chattahoochee-Flint 
   River Basin the watershed above the Chattahoochee River near Norcross, 
   Georgia and is described in LaFontaine and others (2013). This sample 
   illustrates running PRMS with a time-series of dynamic parameters input 
   to the dynamic_param_read module that is described in Regan and LaFontaine (2017).
   
5. acfb_water_use: This sample model is for the Apalachicola-Chattahoochee-Flint 
   River Basin the watershed above the Chattahoochee River near Norcross, 
   Georgia and is described in LaFontaine and others (2013). This sample 
   illustrates running PRMS with a time-series of water-use input to the 
   water_use_read module that is described in Regan and LaFontaine (2017).

6. Tazlina: A model of the Tazlina Basin in Alaske is provided for illustrative
   and educational purposes only. The glacier dynamics simulation method  is
   described in Van Beusekom and Viger (2015).

7. sagehen_restart: This is the sagehen model described in sample problem 1
   to demonstrate the use of the restart option with 16 simulations. The
   restart option is decribe in Regan and LaFontaine (2017) and Regan and 
   others (2015). The first simulation runs the entire simulation period 
   from October 1, 1980, through January 30, 1984. This simulation provides
   results for the entire simulation period. PRMS is then run 15 additional
   times to reproduce the restart simulations. The first of these simulations
   is the 'hindcast' simulation, which runs from October 1, 1980, through
   September 1, 1983. Note that for this run, the simulation 'end_time'
   has been reset to '1983,9,1,0,0,0' in the batch file; several other 
   input control parameters also are reset on the command lines in prms.bat
   for this and each of the subsequent runs.

References:

LaFontaine, J.H., Hay, L.E., Viger, R.J., Markstrom, S.L., Regan, R.S., 
Elliott, C.M., and Jones, J.W., 2013, Application of the Precipitation-
Runoff Modeling System (PRMS) in the Apalachicola–Chattahoochee–Flint 
River Basin in the southeastern United States: U.S. Geological Survey 
Scientific Investigations Report 2013–5162, 118 p., accessed October 13, 
2016, at https://pubs.usgs.gov/sir/2013/5162/.

Markstrom, S.L., Regan, R.S., Hay, L.E., Viger, R.J., Webb, R.M.T., 
Payn, R.A., and LaFontaine, J.H., 2015, PRMS-IV, the precipitation-
runoff modeling system, version 4: U.S. Geological Survey Techniques 
and Methods, book 6, chap. B7, 158 p., http://dx.doi.org/10.3133/tm6B7.

Markstrom, S.L., Niswonger, R.G., Regan, R.S., Prudic, D.E., and
Barlow, P.M., 2008, GSFLOW--Coupled Ground-water and Surface-water
FLOW model based on the integration of the Precipitation-Runoff
Modeling System (PRMS) and the Modular Ground-Water Flow Model
(MODFLOW-2005): U.S. Geological Survey Techniques and Methods
6-D1, 240 p.

Regan, R.S., and LaFontaine, J.H., 2017, Documentation of the dynamic 
parameter, water-use, stream and lake flow routing, and two summary 
output modules and updates to surface-depression storage simulation and 
initial conditions specification options with the Precipitation-Runoff 
Modeling System (PRMS): U.S. Geological Survey Techniques and Methods, 
book 6, chap. B8, 60 p., https://doi.org/10.3133/tm6B8.

Regan, R.S., Niswonger, R.G., Markstrom, S.L., and Barlow, P.M., 2015,
Documentation of a restart option for the U.S. Geological Survey coupled
groundwater and surface-water flow (GSFLOW) model: U.S. Geological Survey
Techniques and Methods, book 6, chap. D3, 19 p., http://dx.doi.org/10.3133/tm6D3.

Van Beusekom, A.E., and Viger, R.J., 2015, A glacier runoff extension to the
Precipitation Runoff Modeling System, Journal of Geophysical Research: Earth 
Science, 21 p., https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JF003789. 