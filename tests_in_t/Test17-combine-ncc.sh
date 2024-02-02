#!/usr/bin/bash

echo "Test combine-ncc combines compressed netcdf files"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test17-input
dir_out=$PWD/Test17-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir $dir_out && cd $dir_out

cp $dir_in/combine-ncc.atmos_daily.nc.copy .

#Combine netcdf copy file
combine-ncc \
  combine-ncc.atmos_daily.nc.copy \
  combine-ncc_output.nc
