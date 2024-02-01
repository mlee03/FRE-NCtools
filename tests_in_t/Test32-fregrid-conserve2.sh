.1;5204;0c#!/usr/bin/bash

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

# test no stretched grid

echo "Test no stretched grid data"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test32-input
dir_out=$PWD/tests_in_t/Test32-conserve2-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

#Make no stretched grid
SECONDS=0
make_hgrid \
  --grid_type gnomonic_ed \
  --nlon 512 \
  --grid_name C256
echo "**** TEST32 SECONDS TO GENERATE C256:  $SECONDS"

#Create no stretched grid mosaic
SECONDS=0
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_mosaic \
  --tile_file C256.tile1.nc,C256.tile2.nc,C256.tile3.nc,C256.tile4.nc,C256.tile5.nc,C256.tile6.nc
echo "**** TEST32 SECONDS TO GENERATE C256_MOSAIC:  $SECONDS"

# no stretched grid
SECONDS=0
fregrid \
  --debug \
  --input_mosaic C256_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin 15.0 \
  --latEnd 65.0 \
  --lonBegin 230.0 \
  --lonEnd 310.0 \
  --interp_method conserve_order2 \
  --remap_file fregrid_remap.nc \
  --check_conserve
echo "**** TEST32 SECONDS TO RUN FREGRID:  $SECONDS"
