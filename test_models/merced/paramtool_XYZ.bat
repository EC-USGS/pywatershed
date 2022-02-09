@ECHO OFF
..\..\bin\prms -C.\control\mercdXYZ.control -print
java -cp ..\..\dist\oui4.jar oui.paramtool.ParamTool .\input\mercdXYZ.param .\control\mercdXYZ.control.par_name
ECHO.
ECHO Run complete. Please press enter to continue.
PAUSE>NUL
