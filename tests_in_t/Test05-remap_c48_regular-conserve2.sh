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

echo "remap data from C48 to regular lat-lon grid with conservative 1rst order"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test05-input
dir_out=$PWD/tests_in_t/Test05-conserve2-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

SECONDS=0
fregrid \
  --debug \
	--input_mosaic $dir_in/C48_mosaic.nc \
	--input_file $dir_in/19800101.atmos_daily \
	--scalar_field zsurf,temp,t_surf \
	--nlon 144 \
	--nlat 90 \
	--interp_method conserve_order2 \
	--output_dir ./ \
	--output_dir ./ \
	--output_file 19800101.atmos_daily.nc \
	--output_file 19800101.atmos_daily.nc \
	--check_conserve \
	--remap_file C48_to_N45_remap.nc
echo "**** TEST05 SECONDS TO REMAP GENERATE TO END $SECONDS"

SECONDS=0
fregrid \
  --debug \
	--input_mosaic $dir_in/C48_mosaic.nc \
	--input_file $dir_in/19800101.atmos_daily \
	--scalar_field zsurf,temp,t_surf \
	--nlon 144 \
	--nlat 90 \
	--interp_method conserve_order2 \
	--output_dir ./ \
	--output_dir ./ \
	--output_file 19800101.atmos_daily.nc \
	--output_file 19800101.atmos_daily.nc \
	--check_conserve \
	--remap_file C48_to_N45_remap.nc
echo "**** TEST05 SECONDS TO REMAP READING TO END $SECONDS"
