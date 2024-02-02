#!/usr/bin/bash

# Test grid for multiple same level and telescoping nests

echo "Test grid for multiple telescope nests"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

out_dir=$PWD/Test27-telescope-output

[[ -d $out_dir ]] ; rm -rf $out_dir
mkdir -p $out_dir && cd $out_dir

make_hgrid \
  --grid_type gnomonic_ed \
  --nlon 96 \
  --grid_name C48_grid \
  --do_schmidt \
  --stretch_factor 1.0 \
  --target_lon -97.5 \
  --target_lat 36.5 \
  --nest_grids 3 \
  --parent_tile 2,5,7 \
  --refine_ratio 2,2,2 \
  --istart_nest 7,13,7 \
  --jstart_nest 7,7,23 \
  --iend_nest 58,68,40 \
  --jend_nest 58,68,48 \
  --halo 3 \
  --great_circle_algorithm \
  --verbose 1

make_solo_mosaic \
  --num_tiles=9 \
  --dir ./ \
  --mosaic_name C48_level1_mosaic \
  --tile_file C48.tile1.nc,C48.tile2.nc,C48.tile3.nc,C48.tile4.nc,C48.tile5.nc,C48.tile6.nc,C48.tile7.nc,C48.tile8.nc,C48.tile9.nc
