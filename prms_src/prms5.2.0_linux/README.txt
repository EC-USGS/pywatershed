

                      PRMS - Version: 5.2.0
            USGS Precipitation-Runoff Modeling System


Any use of trade, product, or firm names is for descriptive purposes only and does
not imply endorsement by the U.S. Government. This software has been approved for
release by the U.S. Geological Survey (USGS). Although the software has been subjected
to rigorous review, the USGS reserves the right to update the software as needed pursuant
to further analysis and review. No warranty, expressed or implied, is made by the USGS
or the U.S. Government as to the functionality of the software and related material nor
shall the fact of release constitute any such warranty. Furthermore, the software is 
released on condition that neither the USGS nor the U.S. Government shall be held liable 
for any damages resulting from its authorized or unauthorized use.

PRMS version 5.2.0 is packaged for personal computers using one of the Linux or 
Microsoft Windows operating systems. An executable file compiled for 64-bit operating 
systems is provided. The Linux executable was compiled using the gfortran and gcc
compilers (8.3.0 version). The Windows executable was compiled using the Microsoft
Visual Studio Community 2019, Version 16.3.7 and the Intel Parallel Studio XE 2019
Update 5 Composer Edition and Microsoft Visual C++ 2019 compilers.

The source code and Linux Makefiles are provided to aid users in compilation
on other computers. However, no support is provided for compilation.

IMPORTANT: Users should review the file 'Release_notes_prms_5.2.0.pdf' for a 
description of, and references for, and changes that have been introduced into PRMS
with each official release; these changes may substantially affect model results.

Instructions for installation, execution, and testing of this version of
PRMS are provided below.



                            TABLE OF CONTENTS

                         A. DISTRIBUTION FILE
                         B. INSTALLING
                         C. EXECUTING THE SOFTWARE
                         D. TESTING
                         E. COMPILING


A. DISTRIBUTION FILE

The following distribution files are for use on personal computers:

         prms_5.2.0.zip for Windows-based personal computers
         prms_5.2.0_linux.zip for Linux-based personal computers

The distribution file contains:

          Executable, source code, and GUI.
          PRMS documentation.
          Five example problems.

Unzipping the distribution file creates numerous individual files contained in 
several subdirectories. The following directory structure will be created in the 
installation directory:

   |
   |--prms5.2.0                ; This file, release notes, and updated tables
   |    |--bin                 ; Compiled PRMS executable
   |    |--projects            ; Five example application models
   |        |--acf             ; scripts for ACF River using the Climate-by-HRU
                                 climate distribution module
   |            |--input       ; Data, Parameter, and Climate-by-HRU Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control File and name files for GUI written here
   |        |--acfb_dyn_params ; scripts for ACF River using the Climate-by-HRU
                                 climate distribution module and dynamic parameters
   |            |--input       ; Data, Parameter, Climate-by-HRU, and dynamic parameter
                                 Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control File and name files for GUI written here
   |        |--acfb_water_use  ; scripts for ACF River using the Climate-by-HRU
                                 climate distribution module
   |            |--input       ; Data, Parameter, Climate-by-HRU, and water use Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control File and name files for GUI written here
   |        |--merced          ; scripts for Merced River example using XYZ and
                                 IDE climate distribution modules
   |            |--input       ; Data and Parameter Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control Files and name files for GUI written here 
   |        |--sagehen         ; scripts for Sagehen River example using temp_1sta
                                 and precip_1sta climate distribution modules
   |            |--input       ; Data and Parameter Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control Files and name files for GUI written here
   |        |--Tazlina         ; scripts for Tazlina glacier example using the
                                 xyz_dist climate distribution module and simulating
 				 glacier and frozen ground dynamics
   |            |--input       ; Data and Parameter Files
   |            |--output      ; output files written here
   |            |--output-test ; output files as produced for comparison purposes
   |            |--control     ; Control Files and name files for GUI written here
   |    |--doc                 ; Documentation reports
   |    |--src                 ; Source files with Makefiles
   |        |--prms            ; Source code for PRMS Modules
   |        |--mmf             ; Source code for MMF utilities
   |    |--dist                ; Java jar files for PRMS GUI


It is recommended that no user files be kept in the prms5.2.0 directory
structure.  If you do plan to put your own files in the directory structure, 
do so only by creating additional subdirectories of the "projects" subdirectory.

Included with the release are several documents that use the Portable Document 
Format (PDF) file structure. The PDF files are readable and printable on various 
computer platforms using Acrobat Reader from Adobe. The Acrobat Reader is freely 
available from the following World Wide Web site: http://www.adobe.com/


B. INSTALLING

To make the executable version of PRMS accessible from any directory, the 
directory containing the executable (prms5.2.0/bin on Linux-based computers or 
prms5.2.0\bin on Windows-based computers) should be included in the PATH
environment variable. Also, if a prior release of PRMS is installed on your 
system, the directory containing the executable for the prior release should be 
removed from the PATH environment variable.

As an alternative, the executable file in the "bin" subdirectory can be copied 
into a directory already included in the PATH environment variable. Note, the example 
problems provided with the release (described below) have example Linux shell script or
Windows batch files that require the executable be in the "bin" subdirectory.


C. EXECUTING THE SOFTWARE

A 64-bit (prms for Linux and prms.exe for Windows) executable is provided in the "bin" 
subdirectory. After the "bin" subdirectory is included in your PATH, PRMS is 
initiated in a Command-Prompt window using the command:

      prms [Fname]

The optional Fname argument is the name of the PRMS Control File.  If no argument is used,
then PRMS will look for a Control File named "control" in the user's current directory.

The arrays in PRMS are dynamically allocated, so models are not limited
by the size of input data. However, it is best to have at least 4 MB of 
random-access memory (RAM) for model execution and more RAM for large models.
If there is less available RAM than the model requires, which depends
on the size of the application, the program will use virtual memory; however,
this can slow execution significantly. If there is insufficient memory to 
run the model, then PRMS will not initiate the beginning of the simulation; 
however, if on a Windows-based computer, the Command-Prompt window may continue 
to indicate that PRMS is executing. For this circumstance, the program must be 
terminated manually using the Windows Task Manager application.

The initial conditions files written by PRMS are unformatted files. The structure
of these files depends on the compiler and options in the code. These files likely
cannot be used on different computers that have different operating systems. For example,
an initial conditions file written on a Windows-based computer may not be usable on a
Linux-based computers. For Windows based computers, PRMS is compiled with the unformatted
file type specified as "BINARY". For Linux-based computers, PRMS is compiled with the 
unformatted file type specified as "UNFORMATTED".


D. TESTING

Five example problems with PRMS data sets are provided in the "projects" subdirectory
to verify that PRMS is correctly installed and running on the user's system.
The example problems also may be looked at as examples of how to use the program. 
See the 'Readme.txt' file in that subdirectory for a description of the five 
example problems.


E. COMPILING

The executable file provided in distribution was created using compilers as described
above. Although executable versions of the program are provided, the source code also
is provided in the "src" subdirectory so that PRMS can be recompiled if necessary.
However, the USGS cannot provide assistance to those compiling PRMS. In general, the 
requirements are a Fortran compiler, a compatible C compiler, and the knowledge
of using the compilers. Makefiles are included in the "src" subdirectories as an example
for compiling PRMS.

