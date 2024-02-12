#!/usr/bin/bash

echo "remap data from C48 to N45 regular lat-lon grid with conservative 2nd order"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test05-input
dir_out=$PWD/Test05-conserve2-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

fregrid \
  --debug \
	--input_mosaic $dir_in/C48_mosaic.nc \
	--input_file $dir_in/19800101.atmos_daily \
	--scalar_field zsurf,temp,t_surf \
	--nlon 144 \
	--nlat 90 \
	--interp_method conserve_order2 \
	--output_dir ./ \
	--output_file 19800101.atmos_daily.nc \
	--check_conserve \
	--remap_file C48_to_N45_remap.nc

fregrid \
  --debug \
	--input_mosaic $dir_in/C48_mosaic.nc \
	--input_file $dir_in/19800101.atmos_daily \
	--scalar_field zsurf,temp,t_surf \
	--nlon 144 \
	--nlat 90 \
	--interp_method conserve_order2 \
	--output_dir ./ \
	--output_file 19800101.atmos_daily.nc \
	--check_conserve \
	--remap_file C48_to_N45_remap.nc
