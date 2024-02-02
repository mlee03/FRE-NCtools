#!/usr/bin/bash

# Test remap_land: remap land restart files.
# Prepare takes a directory name, and generates the required netCDF input
# files for a given test.  prepare_input_data expects a single argument <dir_name>
# which is the directory name in $top_srcdir/my_tests/Test10-inpumy_tests/<dir_name>
# that contains the input data

echo "Test remap_land can remap land restart files"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test10-input
dir_out=$PWD/Test10-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

remap_land \
	--file_type land  \
	--src_mosaic $dir_in/C48_mosaic/C48_mosaic.nc \
	--dst_mosaic $dir_in/C192_mosaic/C192_mosaic.nc \
	--src_restart $dir_in/src_restart/land.res \
	--dst_restart $dir_in/land.res \
	--dst_cold_restart $dir_in/dst_cold_restart/land.res \
	--remap_file remap_file_C48_to_C192 --print_memory
