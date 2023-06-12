call set_run_variables.bat
%SWB2% --output_prefix=%OUTPUT_FILE_PREFIX%   ^
	   --output_dir=%OUTPUT_DIR%              ^
	   --logfile_dir=%LOGFILE_DIR%            ^
       --data_dir=%DATA_DIR%                  ^
       --lookup_dir=%LOOKUP_DIR%              ^
       %SWB_CONTROL_FILE%

