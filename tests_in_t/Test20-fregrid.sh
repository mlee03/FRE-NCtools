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

echo "Test fregrid ocean data"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test20-input
dir_out=$PWD/tests_in_t/Test20-conserve1-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

#Create regular lat-lon grid (100:160, -15:15, size is 360x30)
SECONDS=0
make_hgrid  \
  --grid_type regular_lonlat_grid  \
  --nxbnd 2  \
  --nybnd 2  \
  --xbnd 0,360  \
  --ybnd -15,15  \
  --nlon 720  \
  --nlat 60  \
  --grid_name latlon_grid
echo "**** TEST20 SECONDS TO MAKE LATLON_GRID $SECONDS"

#Create lat-lon mosaic
SECONDS=0
make_solo_mosaic  \
  --num_tiles 1  \
  --dir ./  \
  --mosaic_name latlon_mosaic  \
  --tile_file latlon_grid.nc
echo "**** TEST20 SECONDS TO MAKE LATLON_MOSAIC $SECONDS"

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
SECONDS=0
fregrid   \
  --debug \
	--input_mosaic $dir_in/CM2.1_mosaic.nc   \
	--input_file $dir_in/ocean_temp_salt.res.nc   \
	--scalar_field temp,salt  \
	--output_file ocean_temp_salt.res.latlon.all.nc   \
	--output_file ocean_temp_salt.res.latlon.all.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--check_conserve
echo "**** TEST20 SECONDS TO REMAP GENERATE TO END $SECONDS"


SECONDS=0
fregrid   \
  --debug \
	--input_mosaic $dir_in/CM2.1_mosaic.nc   \
	--input_file $dir_in/ocean_temp_salt.res.nc   \
	--scalar_field temp,salt  \
	--output_file ocean_temp_salt.res.latlon.half.nc   \
	--output_file ocean_temp_salt.res.latlon.half.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--output_mosaic latlon_mosaic.nc   \
	--check_conserve
echo "**** TEST20 SECONDS TO REMAP READ TO END $SECONDS"
