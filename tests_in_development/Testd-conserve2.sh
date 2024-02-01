#!/usr/bin/bash

source run_tools.sh
echo "Testd second order conservative remapping from C768 to lonlat 2304x1440"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testd-input
dir_out=$PWD/tests_in_development/Testd-conserve2-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=atmos_daily

scalar_field="ps_ic,ps,ua_ic,va_ic,ucomp,vcomp,vort,pv,delp,sphum"

input_mosaic=$(run_generate_cubed_sphere 768)
output_mosaic=$(run_generate_lonlat 2304 1140)

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
