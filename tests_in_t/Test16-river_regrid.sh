#!/usr/bin/bash

# Test river_regrid: remap data from lat-lon onto C48 grid

echo "Test river_regrid runoff data from regular lat-lon grid onto cm2m grid"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test16-input
dir_out=$PWD/Test16-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

river_regrid \
	--mosaic $dir_in/grid_spec.nc \
	--river_src $dir_in/z1l_river_output_M45_tripolar_aug24.nc \
	--output river_data_C48 \
	--land_thresh 0.000001
