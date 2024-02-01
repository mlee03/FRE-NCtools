#!/usr/bin/bash

source run_tools.sh
echo "Testk first order conservative remapping from 1/2 tripolar to lonlat 360x180"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testk-input
dir_out=$PWD/Testk-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file="00010101.ocean_z_month"

scalar_field=\
"S_advection_xy,Sh_tendency_vert_remap,T_advection_xy,Th_tendency_vert_remap,
boundary_forcing_heat_tendency,boundary_forcing_salt_tendency,difvho,
frazil_heat_tendency,opottempdiff,opottemppmdiff,opottemptend,osaltdiff,
osaltpmdiff,osalttend,rhoinsitu,rhopot0,rhopot2,rsdoabsorb,so,thetao,
uh,uhGM,uhml,umo,uo,vh,vhGM,vhml,vmo,vo"

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
