#!/usr/bin/bash

# test make_quick_mosaic
#First create an ocean_mosaic and ocean_topog.nc
#Make_hgrid: create ocean_hgrid"

echo "Test make_quick_mosaic"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test13-input
dir_out=$PWD/Test13-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

#Make the quick mosaic
make_quick_mosaic \
		--input_mosaic $dir_in/ocean_mosaic.nc \
		--ocean_topog $dir_in/ocean_topog.nc
