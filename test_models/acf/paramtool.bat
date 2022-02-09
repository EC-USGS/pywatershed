@ECHO OFF
..\..\bin\prms -C.\control\acf.control -set print_debug -2 -print
java -cp ..\..\dist\oui4.jar oui.paramtool.ParamTool .\input\acf.params .\control\acf.control.par_name
ECHO.
ECHO Run complete. Please press enter to continue.
PAUSE>NUL
