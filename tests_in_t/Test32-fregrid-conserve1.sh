#!/usr/bin/bash

# test no stretched grid

echo "Test no stretched grid data"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test32-input
dir_out=$PWD/Test32-conserve1-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

#Make no stretched grid
make_hgrid \
  --grid_type gnomonic_ed \
  --nlon 512 \
  --grid_name C256

#Create no stretched grid mosaic
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_mosaic \
  --tile_file C256.tile1.nc,C256.tile2.nc,C256.tile3.nc,C256.tile4.nc,C256.tile5.nc,C256.tile6.nc

# no stretched grid
fregrid \
  --debug \
  --input_mosaic C256_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin 15.0 \
  --latEnd 65.0 \
  --lonBegin 230.0 \
  --lonEnd 310.0 \
  --remap_file fregrid_remap.nc \
  --interp_method conserve_order1 \
  --check_conserve
