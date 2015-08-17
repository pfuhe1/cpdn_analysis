#!/usr/local/bin/convsh
#
#set umid $argv(0)
readfile 0 atmos_restart.day 
writefile netcdf atmos_${argv}_sm.nc 4
clearall
readfile 0 region_restart.day
writefile netcdf region_${argv}_sm.nc 4
