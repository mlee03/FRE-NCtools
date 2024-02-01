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

# Test river_regrid: remap data from lat-lon onto C48 grid

echo "Test  remap runoff data from regular lat-lon grid onto cm2m grid"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test16-input
dir_out=$PWD/tests_in_t/Test16-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

SECONDS=0
river_regrid \
	--mosaic $dir_in/grid_spec.nc \
	--river_src $dir_in/z1l_river_output_M45_tripolar_aug24.nc \
	--river_src $dir_in/z1l_river_output_M45_tripolar_aug24.nc \
	--output river_data_C48 \
	--output river_data_C48 \
	--land_thresh 0.000001
echo "**** TEST16 SECONDS TO RUN RIVER_REGRID $SECONDS"
