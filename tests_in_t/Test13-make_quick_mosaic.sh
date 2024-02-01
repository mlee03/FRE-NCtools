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

# test make_quick_mosaic

#First create an ocean_mosaic and ocean_topog.nc
#Make_hgrid: create ocean_hgrid"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test13-input
dir_out=$PWD/tests_in_t/Test13-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

#Make the quick mosaic
SECONDS=0
make_quick_mosaic \
		--input_mosaic $dir_in/ocean_mosaic.nc \
		--ocean_topog $dir_in/ocean_topog.nc
echo "**** TEST13 SECONDS TO RUN MAKE_QUICK_MOSAIC $SECONDS"
