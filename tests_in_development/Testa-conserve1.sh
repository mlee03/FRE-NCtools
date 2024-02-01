#!/usr/bin/bash

source run_tools.sh
echo "Testa first order conservative remapping from C96 to lonlat 288x180"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testa-input
dir_out=$PWD/Testa-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=00010101.atmos_month_aer

scalar_field=\
"zsurf,ps,temp,sphum,sulfate,sulfate_col,sm_dust,sm_dust_col,\
lg_dust,lg_dust_col,salt,salt_col,blk_crb,blk_crb_col,org_crb,\
org_crb_col,sulfate_ex_c_vs,sm_dst_ex_c_vs,lg_dst_ex_c_vs,\
blk_crb_ex_c_vs,org_crb_ex_c_vs,salt_ex_c_vs,aer_ex_c_vs,\
aer_ab_c_vs,aer_c,aer_ex_vs,aer_ab_vs"

input_mosaic=$(run_generate_cubed_sphere 96)
output_mosaic=$(run_generate_lonlat 288 180)

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
