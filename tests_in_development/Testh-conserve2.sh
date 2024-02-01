#!/usr/bin/bash

source run_tools.sh
echo "Testh second order conservative remapping from C128 to lonlat 384x240"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testh-input
dir_out=$PWD/tests_in_development/Testh-conserve2-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file="atmos_4xdaily_ave"

scalar_field="wmaxup_max"

input_mosaic=$(run_generate_cubed_sphere 128)
output_mosaic=$(run_generate_lonlat 384 240)

run_fregrid_lonlat_conserve2 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field

run_fregrid_lonlat_conserve2 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field
