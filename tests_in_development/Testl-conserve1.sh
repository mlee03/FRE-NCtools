#!/usr/bin/bash

source run_tools.sh
echo "Testl first order conservative remapping from 1/4 tripolar to lonlat 360x180"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testl-input
dir_out=$PWD/tests_in_development/Testl-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file="00010101.ocean_month"

scalar_field=\
"Heat_PmE,LW,LwLatSens,MLD_003,MLD_003_max,MLD_003_min,ML_buoy_restrat,PRCmE,\
ssh,sss,SST,sst_max,sst_min,SW,TKE_tidal,evap,evs,ficeberg,fprec,friver,\
frunoff,heat_content_cond,heat_content_fprec,heat_content_frunoff,heat_content_lprec,\
heat_content_lrunoff,heat_content_massin,heat_content_massout,heat_content_surfwater,\
hfds,hfevapds,hfibthermds,hflso,hfrainds,hfrunoffds,hfsifrazil,hfsnthermds,hfsso,\
latent,lprec,lrunoff,net_heat_coupler,net_heat_surface,net_massin,net_massout,\
p_surf,prlq,prsn,rlntds,rsntds,salt_flux,sensible,sfdsi,sfdsi,speed,ustar,wfo"

input_mosaic=$dir_in/ocean_mosaic
output_mosaic=$(run_generate_lonlat 360 180)

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
