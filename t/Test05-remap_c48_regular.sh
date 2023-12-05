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

echo "remap data from C48 to regular lat-lon grid"

dir_in=$PWD/t/Test05-input
dir_out=$PWD/t/Test05-output

mkdir -p $dir_out
for ncl_file in $dir_in/*.ncl ; do
  nc_file=${ncl_file/'.ncl'/'.nc'}
  ncgen $ncl_file -o $nc_file
done


fregrid \
		--input_mosaic $dir_in/C48_mosaic.nc \
		--input_file $dir_in/19800101.atmos_daily \
		--scalar_field zsurf,temp,t_surf \
		--nlon 144 \
		--nlat 90 \
		--interp_method conserve_order2 \
		--output_dir $dir_out \
		--output_file 19800101.atmos_daily.nc \
		--check_conserve \
		--remap_file $dir_out/C48_to_N45_remap.nc
