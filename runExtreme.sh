#!/bin/bash
#ifort  -o extreme.exe -heap-arrays 1000000000 ExtremeRainfallAndTemp.F90 -I/usr/local/include/netcdf.inc -L/usr/local/lib -lnetcdff
ifort -O3 -o nexER.exe -mcmodel=large -heap-arrays 1000000000 ExtremeRainNEX.F90 -I/opt/netcdf/include/netcdf.inc -L/opt/netcdf/lib -lnetcdff

