:: Hail Mary attempt to keep netCDF from erroring out based on this SO entry:
::https://stackoverflow.com/questions/49317927/errno-101-netcdf-hdf-error-when-opening-netcdf-file
::
set RES=5000
set HDF5_USE_FILE_LOCKING=FALSE
set LOGFILE_DIR=../swb_logfiles
set OUTPUT_DIR=../swb_output
set DATA_DIR=../
set LOOKUP_DIR=../
set SWB_CONTROL_FILE=../swb_control_file.ctl
set OUTPUT_FILE_PREFIX=hru_1_%RES%__

set SWB2=swb2.exe
set SWBSTATS2=swbstats2.exe
