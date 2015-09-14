#!/usr/local/bin/convsh
#
set infile [lindex $argv 0]
set outfile [lindex $argv 1]
readfile 0 ${infile}
writefile netcdf $outfile 4
