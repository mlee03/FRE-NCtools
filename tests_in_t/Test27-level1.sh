#!/usr/bin/bash

# Test grid for multiple same level and telescoping nests

echo "Test fregrid for multiple same level nests great_circle"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

out_dir=$PWD/Test27-level1-output

[[ -f $out_dir ]] ; rm -rf $out_dir
mkdir -p $out_dir ; cd $out_dir

#Make_hgrid: create three same level -level1- nests in tiles 2,5,6"
make_hgrid \
	--grid_type gnomonic_ed \
	--nlon 96 \
	--grid_name C48 \
	--nest_grids 3 \
	--parent_tile 2,5,6 \
  --refine_ratio 2,2,2 \
  --istart_nest 7,13,7 \
  --iend_nest 58,68,40 \
  --jstart_nest 7,7,23 \
  --jend_nest 58,68,48 \
  --halo 3 \
  --verbose 1\
  --great_circle_algorithm

make_solo_mosaic \
  --num_tiles=9 \
  --dir ./ \
  --mosaic_name C48_level1_mosaic \
  --tile_file C48.tile1.nc,C48.tile2.nc,C48.tile3.nc,C48.tile4.nc,C48.tile5.nc,C48.tile6.nc,C48.tile7.nc,C48.tile8.nc,C48.tile9.nc
