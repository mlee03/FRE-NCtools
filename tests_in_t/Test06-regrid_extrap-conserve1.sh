#!/usr/bin/bash

# Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation

echo "Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation conserve 1"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test06-input
dir_out=$PWD/Test06-conserve1-output

[[ -d $dir_out ]] ; rm $dir_out
mkdir -p $dir_our && cd $dir_out

make_hgrid \
		--grid_type regular_lonlat_grid \
		--nxbnd 2 \
		--nybnd 2 \
		--xbnd 0,360 \
		--ybnd -90,90 \
		--nlon 720 \
		--nlat 360 \
		--grid_name levitus_grid

make_solo_mosaic \
	--num_tiles 1 \
	--dir ./ \
	--mosaic_name levitus_mosaic \
	--tile_file levitus_grid.nc \
	--periodx 360

fregrid \
  --debug \
	--input_mosaic levitus_mosaic.nc \
	--input_file $dir_in/WOA09_ann_theta.nc \
	--scalar_field POTENTIAL_TEMP \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--extrapolate \
  --interp_method conserve_order1 \
	--dst_vgrid $dir_in/ocean_vgrid.nc \
	--check_conserve

fregrid \
  --debug \
	--input_mosaic levitus_mosaic.nc \
	--input_file $dir_in/WOA09_ann_theta.nc \
	--scalar_field POTENTIAL_TEMP \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--extrapolate \
  --interp_method conserve_order1 \
	--dst_vgrid $dir_in/ocean_vgrid.nc \
	--check_conserve
