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

# test stretched grid

echo "Test stretched grid data lats 32.0 34.0 35.4"

dir_in=$PWD/t/Test31-input
dir_out=$PWD/t/Test31-output
mkdir -p $dir_out
cd $dir_out

#Make streetched grid
make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 2.5 \
  --target_lon 262.4 \
  --target_lat 32.0 \
  --nlon 512 \
  --grid_name C256_grid_32.0


make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 2.5 \
  --target_lon 262.4 \
  --target_lat 34.0 \
  --nlon 512 \
  --grid_name C256_grid_34.0


make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 2.5 \
  --target_lon 262.4 \
  --target_lat 35.4 \
  --nlon 512 \
  --grid_name C256_grid_35.4


#Create stretched grid mosaic
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_mosaic_32.0 \
  --tile_file C256_grid_32.0.tile1.nc,C256_grid_32.0.tile2.nc,C256_grid_32.0.tile3.nc,C256_grid_32.0.tile4.nc,C256_grid_32.0.tile5.nc,C256_grid_32.0.tile6.nc


make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_mosaic_34.0 \
  --tile_file C256_grid_34.0.tile1.nc,C256_grid_34.0.tile2.nc,C256_grid_34.0.tile3.nc,C256_grid_34.0.tile4.nc,C256_grid_34.0.tile5.nc,C256_grid_34.0.tile6.nc


make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_mosaic_35.4 \
--tile_file C256_grid_35.4.tile1.nc,C256_grid_35.4.tile2.nc,C256_grid_35.4.tile3.nc,C256_grid_35.4.tile4.nc,C256_grid_35.4.tile5.nc,C256_grid_35.4.tile6.nc


# stretched grid lats 32.0 34.0 35.4
fregrid \
  --input_mosaic C256_mosaic_32.0.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin 15.0 \
  --latEnd 65.0 \
  --lonBegin 230.0 \
  --lonEnd 310.0 \
  --remap_file fregrid_remap_file_640_by_400_32.0.nc \
  --output_file out_32.0.nc \
  --check_conserve

fregrid \
  --input_mosaic C256_mosaic_34.0.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin 15.0 \
  --latEnd 65.0 \
  --lonBegin 230.0 \
  --lonEnd 310.0 \
  --remap_file fregrid_remap_file_640_by_400_34.0.nc \
  --output_file out_34.0.nc \
  --check_conserve


fregrid \
  --input_mosaic C256_mosaic_35.4.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin 15.0 \
  --latEnd 65.0 \
  --lonBegin 230.0 \
  --lonEnd 310.0 \
  --remap_file fregrid_remap_file_640_by_400_35.4.nc \
  --output_file out_35.4.nc \
  --check_conserve
