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

dir_in=$PWD/t/Test20-input
dir_out=$PWD/t/Test20-output
mkdir -p $dir_out

for ncl_file in $dir_in/*.ncl ; do
  nc_file=${ncl_file/'.ncl'/'.nc'}
  ncgen $ncl_file -O $nc_file
done


#Create regular lat-lon grid (100:160, -15:15, size is 360x180)
   make_hgrid  \
     --grid_type regular_lonlat_grid  \
     --nxbnd 2  \
     --nybnd 2  \
     --xbnd 0,360  \
     --ybnd -2,2  \
     --nlon 720  \
     --nlat 10  \
     --grid_name $dir_out/latlon_grid

   #Create lat-lon mosaic
   make_solo_mosaic  \
     --num_tiles 1  \
     --dir $dir_out  \
     --mosaic_name $dir_out/latlon_mosaic  \
     --tile_file latlon_grid.nc

#Remap data from CM2.1 ocean grid onto regular lat-lon grid.
   fregrid   \
		--input_mosaic $dir_in/CM2.1_mosaic.nc   \
		--input_file $dir_in/ocean_temp_salt.res.nc   \
		--scalar_field temp,salt  \
		--output_file $dir_out/ocean_temp_salt.res.latlon.nc   \
		--output_mosaic $dir_out/latlon_mosaic.nc   \
		--check_conserve

   for ref_nc in $PWD/t/Test20-reference/*.nc ; do
     out_nc=${ref_nc/'reference'/'output'}
     nccmp -dm $ref_nc $out_nc
   done
