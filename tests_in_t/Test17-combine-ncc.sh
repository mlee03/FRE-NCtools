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

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test17-input
dir_out=$PWD/tests_in_t/Test17-output

rm -rf $dir_out
mkdir $dir_out
cd $dir_out

echo "combine-ncc combines compressed netcdf files"

cp $dir_in/combine-ncc.atmos_daily.nc.copy .

#Combine netcdf copy file
SECONDS=0
combine-ncc \
  combine-ncc.atmos_daily.nc.copy \
  combine-ncc_output.nc
  combine-ncc_output.nc
echo "**** TEST17 SECONDS TO RUN COMBINE-NCC"
