# FastHenry

This project is a fork of the repo by [edilorean](https://github.com/ediloren/FastHenry2). The source code has been modified, and some Linux dependent references have been removed such as those related to time calculations. This allows the code to be compiled on Windows using MSYS     

FastHenry is the premium inductance solver originally developed at M.I.T. on Unix platform. A de-facto golden reference standard, FastHenry extracts the inductances and resistances of any arbitrary 3D conductive geometry by solving the Maxwell equations in quasi-static regime.

## Compile Instructions - Windows

1. Clone repo
1. Download [MSYS](https://www.msys2.org/) and follow istallation instructions to install the mingw-w64 GCC tools. 
1. Install [Make](https://packages.msys2.org/packages/make) package
1. Restart MSYS terminal 
1. Navigate to cloned folder 
1. Run command `make fasthenry`
1. If compiled successfully, an exe file will be created in the bin folder
1. optional: run `make clean` to clean project directory and remove object files. 
