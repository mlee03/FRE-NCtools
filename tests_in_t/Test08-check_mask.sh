#!/usr/bin/bash

# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a sea-ice model, baltic1 experiment

echo "Test check_mask for baltic1 experiment"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test08-input
dir_out=$PWD/Test08-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir $dir_out && cd $dir_out

check_mask \
	--grid_file $dir_in/baltic1_grid_spec.nc \
	--min_pe 60 \
	--max_pe 80
