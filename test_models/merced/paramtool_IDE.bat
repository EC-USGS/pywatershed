@ECHO OFF
..\..\bin\prms -C.\control\mercdIDE.control -print
java -cp ..\..\dist\oui4.jar oui.paramtool.ParamTool .\input\mercdIDE.param .\control\mercdIDE.control.par_name
ECHO.
ECHO Run complete. Please press enter to continue.
PAUSE>NUL
