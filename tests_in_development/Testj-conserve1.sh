#!/usr/bin/bash

source run_tools.sh
echo "Testj first order conservative remapping from C512 to lonlat 1536x960"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testj-input
dir_out=$PWD/tests_in_development/Testj-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file="atmos_1min"

scalar_field="ps,tb,qn,qp"

input_mosaic=$(run_generate_cubed_sphere 512)
output_mosaic=$(run_generate_lonlat 1536 980)

run_fregrid_lonlat_conserve1 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field

run_fregrid_lonlat_conserve1 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field
