#!/usr/bin/bash

source run_tools.sh
echo "Testb first order conservative remapping from C192 to lonlat 576x360"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testb-input
dir_out=$PWD/tests_in_development/Testb-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=19790101.atmos_daily_cmip

scalar_field=\
"rsdt,rsut,rsdscs,rsuscs,rldscs,rlutcs,rsutcs,clwvi,clivi,wap500,ta700_unmsk,\
ccb,cct,rsdsdiff,rsdscsdiff,cldnvi,ts,prw,zg_plev19,ps,huss,tasmin,tasmax,tas,\
height2m,pr,psl,sfcWind,height10m,hurs,hursmin,hursmax,clt,tslsi,prc,prsn,uas,\
vas,sfcWindmax,hfls,hfss,rlds,rlus,rsds,rsus,rlut,ta_unmsk,hur_unmsk,hus_unmsk,\
wap_unmsk,va_unmsk,ua_unmsk"

input_mosaic=$(run_generate_cubed_sphere 192)
output_mosaic=$(run_generate_lonlat 576 360)

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
