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

dir_in=$PWD/t/Test33-input
dir_out=$PWD/t/Test33-output
mkdir -p $dir_out

gs=$dir_in/"fv3_48_grid_spec.tile"   #name prefix for the internal fv3 grid spec files.

#Ib) Call make_hgrid to generate the NCTools grids from FV3 files.
## TODO:  Note argument of --great_circle_algorithm option, which should not
## be necessary and code may be corrected in future.
make_hgrid --grid_type from_file \
	   --great_circle_algorithm \
	   --my_grid $gs"1.nc",$gs"2.nc",$gs"3.nc",$gs"4.nc",$gs"5.nc",$gs"6.nc" \
	   --grid_name $dir_out/C48_ff_grid_spec


#II) "Analytically" create a equal distance gnomonic cubic grid using
#     the do_cube_transform option
make_hgrid \
    --grid_type gnomonic_ed \
    --nlon 96 \
    --do_cube_transform \
    --stretch_factor 3 \
    --target_lat 17.5 \
    --target_lon 172.5 \
    --grid_name $dir_out/C48_grid_spec


# III)  Compare the six tile files generated "analytically" to the corresponding ones
##    generated from FV3 grids. Note that we are comparing files with origins from
## two different systems, one Fortran based and another C based, and exact matching
# will generally not be possible. The two tolerances chosen below were determined
# by running on hardware, with FV3 origin files generated on Intel Skylake (C4) and
# the NCTools analytically generated files created on various AMD and Intel hardware
# (e.g. Intel analysis nodes, AMD T5 nodes, various home computers).
for i in 1 2 3 4 5 6
do
    fv3_file=$dir_out/"C48_ff_grid_spec.tile"$i".nc"
    nct_file=$dir_out/"C48_grid_spec.tile"$i".nc"
    nccmp -d --variable=x,y,dx,dy --Tolerance=1.0e-9 $fv3_file $nct_file
    nccmp -d --variable=area --Tolerance=1.0e-6 $fv3_file $nct_file
    # TODO: angle_dx and angle_dy may be done in future.
done
