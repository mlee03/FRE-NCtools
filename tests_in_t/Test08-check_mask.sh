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

# Test check_mask: create mask_table to mask out some all-land domain
# to save processor usage for a sea-ice model, baltic1 experiment


echo "Test check_mask for baltic1 experiment"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test08-input
dir_out=$PWD/tests_in_t/Test08-output

rm -rf $dir_out
mkdir $dir_out
cd $dir_out

SECONDS=0
check_mask \
	--grid_file $dir_in/baltic1_grid_spec.nc \
	--min_pe 60 \
	--max_pe 80
echo "**** TEST08 SECONDS TO RUN CHECK_MASK"
