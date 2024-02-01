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
dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test18-input
dir_out=$PWD/tests_in_t/Test18-output

rm -rf $dir_out
mkdir $dir_out
cd $dir_out

echo "decompress input netcdf files"

cp $dir_in/decompress-ncc.atmos_daily.nc.copy .

#Decompress compressed netcdf file(s) into 1
SECONDS=0
decompress-ncc \
  decompress-ncc.atmos_daily.nc.copy \
  decompress-ncc_output.nc
  decompress-ncc_output.nc
echo "**** TEST18 SECONDS TO RUN DECOMPRESS-NCC"
