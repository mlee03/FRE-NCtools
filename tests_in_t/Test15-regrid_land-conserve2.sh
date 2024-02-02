#!/usr/bin/bash

# Test regrid land data with cell_measures and cell_methods attribute

echo "Test regrid land data 2nd order conservative"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test15-input
dir_out=$PWD/Test15-conserve2-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

#remap static field
fregrid  \
  --debug \
	--input_mosaic $dir_in/C180_mosaic.nc  \
	--interp_method conserve_order2  \
	--nlon 144  \
	--nlat 90  \
	--input_file $dir_in/00050101.land_static  \
	--scalar_field soil_frac,lake_frac,glac_frac,area,soil_area,lake_area,glac_area \
	--output_file out_all.nc  \
	--remap_file remap_file.nc

#remap static field
fregrid  \
  --debug \
	--input_mosaic $dir_in/C180_mosaic.nc  \
	--interp_method conserve_order2  \
	--nlon 144  \
	--nlat 90  \
	--input_file $dir_in/00050101.land_static  \
	--scalar_field soil_frac,lake_frac,glac_frac,area,soil_area,lake_area,glac_area \
	--output_file out_half.nc  \
	--remap_file remap_file.nc
