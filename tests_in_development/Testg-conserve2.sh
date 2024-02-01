#!/usr/bin/bash

source run_tools.sh
echo "Testg second order conservative remapping from C3072 to lonlat 11520x5760"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testf-input
dir_out=$PWD/tests_in_development/Testg-conserve2-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=""

scalar_field="CAPE_max,CIN_max,MXUPHL2_5km_max,MNUPHL2_5km_min,MAXUVV_max,MAXDVV_min"

input_mosaic=$(run_generate_cubed_sphere 3072)
output_mosaic=$(run_generate_lonlat 11520 5760)

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
