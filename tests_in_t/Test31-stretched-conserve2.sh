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

echo "Test stretched grid data stretch_factor 1.5, 2.5, 3.5"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test31-input
dir_out=$PWD/tests_in_t/Test31-conserve2-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

#Make streetched grid
SECONDS=0
make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 1.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_1.5
echo "**** TEST31 SECONDS TO GENERATE C256_1.5:  $SECONDS"

SECONDS=0
make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 2.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_2.5
echo "**** TEST31 SECONDS TO GENERATE C256_2.5:  $SECONDS"

SECONDS=0
make_hgrid \
  --grid_type gnomonic_ed --do_schmidt \
  --stretch_factor 3.5 \
  --target_lon 285.35 \
  --target_lat 40.34 \
  --nlon 512 \
  --grid_name C256_3.5
echo "**** TEST31 SECONDS TO GENERATE C256_3.5:  $SECONDS"

#Create stretched grid mosaic
SECONDS=0
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_1.5_mosaic \
  --tile_file C256_1.5.tile1.nc,C256_1.5.tile2.nc,C256_1.5.tile3.nc,C256_1.5.tile4.nc,C256_1.5.tile5.nc,C256_1.5.tile6.nc
echo "**** TEST31 SECONDS TO GENERATE C256_1.5_mosaic:  $SECONDS"

SECONDS=0
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_2.5_mosaic \
  --tile_file C256_2.5.tile1.nc,C256_2.5.tile2.nc,C256_2.5.tile3.nc,C256_2.5.tile4.nc,C256_2.5.tile5.nc,C256_2.5.tile6.nc
echo "**** TEST31 SECONDS TO GENERATE C256_2.5_mosaic:  $SECONDS"

SECONDS=0
make_solo_mosaic \
  --num_tiles 6 \
  --dir ./ \
  --mosaic_name C256_3.5_mosaic \
  --tile_file C256_3.5.tile1.nc,C256_3.5.tile2.nc,C256_3.5.tile3.nc,C256_3.5.tile4.nc,C256_3.5.tile5.nc,C256_3.5.tile6.nc
echo "**** TEST31 SECONDS TO GENERATE C256_3.5_mosaic:  $SECONDS"

# stretched grid
SECONDS=0
fregrid \
  --debug \
  --input_mosaic C256_1.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order2 \
  --remap_file fregrid_remap_file_1.5.nc \
  --check_conserve
echo "**** TEST31 SECONDS TO RUN FREGRD 1.5 :  $SECONDS"

# stretched grid
SECONDS=0
fregrid \
  --debug \
  --input_mosaic C256_2.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order2 \
  --remap_file fregrid_remap_file_2.5.nc \
--check_conserve
echo "**** TEST31 SECONDS TO RUN FREGRD 2.5 :  $SECONDS"

# stretched grid
SECONDS=0
fregrid \
  --debug \
  --input_mosaic C256_3.5_mosaic.nc \
  --nlon 640 \
  --nlat 400 \
  --latBegin -20. \
  --latEnd 90. \
  --lonBegin -180. \
  --lonEnd 0. \
  --interp_method conserve_order2 \
  --remap_file fregrid_remap_file_3.5.nc \
  --check_conserve
echo "**** TEST31 SECONDS TO RUN FREGRD 3.5 :  $SECONDS"
