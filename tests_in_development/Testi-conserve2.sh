#!/usr/bin/bash

source run_tools.sh
echo "Testi second order conservative remapping from C256 to lonlat 768x480"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testi-input
dir_out=$PWD/tests_in_development/Testi-conserve2-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file="atmos_5min"

scalar_field="ps,lw,u850,v850,w850,vort850,w,temp,rainwat,liq_wat,uh25"

input_mosaic=$(run_generate_cubed_sphere 256)
output_mosaic=$(run_generate_lonlat 768 480)

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
