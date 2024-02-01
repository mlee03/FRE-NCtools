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

echo "Test make_hgrid do_transform vs. FV3 generated grid spec"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test33-input
dir_out=$PWD/tests_in_t/Test33-output
mkdir -p $dir_out
cd $dir_out

gs=$dir_in/"fv3_48_grid_spec.tile"   #name prefix for the internal fv3 grid spec files.

#Ib) Call make_hgrid to generate the NCTools grids from FV3 files.
## TODO:  Note argument of --great_circle_algorithm option, which should not
## be necessary and code may be corrected in future.
SECONDS=0
make_hgrid --grid_type from_file \
	   --great_circle_algorithm \
	   --my_grid $gs"1.nc",$gs"2.nc",$gs"3.nc",$gs"4.nc",$gs"5.nc",$gs"6.nc" \
	   --grid_name  C48_ff_grid_spec
echo " **** TEST33 SECONDS TO MAKE_HGRID FROM FILE"


#II) "Analytically" create a equal distance gnomonic cubic grid using
#     the do_cube_transform option
SECONDS=0
make_hgrid \
    --grid_type gnomonic_ed \
    --nlon 96 \
    --do_cube_transform \
    --stretch_factor 3 \
    --target_lat 17.5 \
    --target_lon 172.5 \
    --grid_name C48_grid_spec
echo "**** TEST33 SECONDS TO GENERATE MAKE_HGRID DO_CUBE_TRANSFORM"
