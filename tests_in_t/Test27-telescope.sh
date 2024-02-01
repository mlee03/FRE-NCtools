#!/usr/bin/bash

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************

# Test grid for multiple same level and telescoping nests

echo "Test grid for multiple telescope nests"

out_dir=$PWD/tests_in_t/Test27-telescope-output

rm -rf $out_dir
mkdir -p $out_dir
cd $out_dir

SECONDS=0
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
echo "**** TEST27 SECONDS TO MAKE TELESCOPE NESTS $SECONDS"


#make_solo_mosaic \
#  --num_tiles=9 \
#  --dir ./ \
#  --mosaic_name C48_level1_mosaic \
#  --tile_file C48.tile1.nc,C48.tile2.nc,C48.tile3.nc,C48.tile4.nc,C48.tile5.nc,C48.tile6.nc,C48.tile7.nc,C48.tile8.nc,C48.tile9.nc
