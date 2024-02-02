#!/usr/bin/bash

echo "Test decompress input netcdf files"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test18-input
dir_out=$PWD/Test18-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir $dir_out && cd $dir_out

cp $dir_in/decompress-ncc.atmos_daily.nc.copy .

#Decompress compressed netcdf file(s) into 1
decompress-ncc \
  decompress-ncc.atmos_daily.nc.copy \
  decompress-ncc_output.nc
