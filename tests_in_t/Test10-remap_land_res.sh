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

# Test remap_land: remap land restart files.

# Prepare takes a directory name, and generates the required netCDF input
# files for a given test.  prepare_input_data expects a single argument <dir_name>
# which is the directory name in $top_srcdir/my_tests/Test10-inpumy_tests/<dir_name>
# that contains the input data

echo "Test remap_land can remap land restart files"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test10-input
dir_out=$PWD/tests_in_t/Test10-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

SECONDS=0
remap_land \
	--file_type land  \
	--src_mosaic $dir_in/C48_mosaic/C48_mosaic.nc \
	--dst_mosaic $dir_in/C192_mosaic/C192_mosaic.nc \
	--src_restart $dir_in/src_restart/land.res \
	--dst_restart $dir_in/land.res \
	--dst_cold_restart $dir_in/dst_cold_restart/land.res \
	--remap_file remap_file_C48_to_C192 --print_memory
echo "**** TEST010 SECONDS TO RUN REMAP_LAND $SECONDS"
