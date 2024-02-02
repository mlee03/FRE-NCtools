#!/usr/bin/bash

echo "Test fregrid ocean data conserve 1"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test20-input
dir_out=$PWD/Test20-conserve1-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

#Create regular lat-lon grid (100:160, -15:15, size is 360x30)
make_hgrid  \
  --grid_type regular_lonlat_grid  \
  --nxbnd 2  \
  --nybnd 2  \
  --xbnd 0,360  \
  --ybnd -15,15  \
  --nlon 720  \
  --nlat 60  \
  --grid_name latlon_grid

#Create lat-lon mosaic
make_solo_mosaic  \
  --num_tiles 1  \
  --dir ./  \
  --mosaic_name latlon_mosaic  \
  --tile_file latlon_grid.nc

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
fregrid   \
  --debug \
	--input_mosaic $dir_in/CM2.1_mosaic.nc   \
	--input_file $dir_in/ocean_temp_salt.res.nc   \
	--scalar_field temp,salt  \
	--output_file ocean_temp_salt.res.latlon.all.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--check_conserve

fregrid   \
  --debug \
	--input_mosaic $dir_in/CM2.1_mosaic.nc   \
	--input_file $dir_in/ocean_temp_salt.res.nc   \
	--scalar_field temp,salt  \
	--output_file ocean_temp_salt.res.latlon.half.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--check_conserve
