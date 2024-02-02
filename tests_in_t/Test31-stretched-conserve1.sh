#!/usr/bin/bash

# test stretched grid

echo "Test stretched grid data stretch_factor 1.5, 2.5, 3.5 conserve 1"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test31-input
dir_out=$PWD/Test31-conserve1-output

[[ -d $dir_out [[ ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

#Make streetched grid
make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 1.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_1.5

make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 2.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_2.5

make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 3.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_3.5

#Create stretched grid mosaic
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_1.5_mosaic \
  --tile_file C256_1.5.tile1.nc,C256_1.5.tile2.nc,C256_1.5.tile3.nc,C256_1.5.tile4.nc,C256_1.5.tile5.nc,C256_1.5.tile6.nc

make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_2.5_mosaic \
  --tile_file C256_2.5.tile1.nc,C256_2.5.tile2.nc,C256_2.5.tile3.nc,C256_2.5.tile4.nc,C256_2.5.tile5.nc,C256_2.5.tile6.nc

make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_3.5_mosaic_35.4 \
  --tile_file C256_3.5.tile1.nc,C256_3.5.tile2.nc,C256_3.5.tile3.nc,C256_3.5.tile4.nc,C256_3.5.tile5.nc,C256_3.5.tile6.nc

# stretched grid
fregrid \
  --debug \
  --input_mosaic C256_1.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order1 \
  --remap_file fregrid_remap_file_1.5.nc \
  --check_conserve

# stretched grid
fregrid \
  --debug \
  --input_mosaic C256_2.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order1 \
  --remap_file fregrid_remap_file_2.5.nc \
--check_conserve

# stretched grid
fregrid \
  --debug \
  --input_mosaic C256_3.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order1 \
  --remap_file fregrid_remap_file_3.5.nc \
  --check_conserve
