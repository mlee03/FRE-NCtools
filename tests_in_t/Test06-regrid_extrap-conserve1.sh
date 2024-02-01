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

# Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation

echo "Test remap data onto cm2m ocean grid with extrapolation and vertical interpolation"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test06-input
dir_out=$PWD/tests_in_t/Test06-conserve1-output

rm -rf $dir_out
mkdir -p $dir_out
cd $dir_out

SECONDS=0
make_hgrid \
		--grid_type regular_lonlat_grid \
		--nxbnd 2 \
		--nybnd 2 \
		--xbnd 0,360 \
		--ybnd -90,90 \
		--nlon 720 \
		--nlat 360 \
		--grid_name levitus_grid
echo "**** TEST06 SECONDS TO GENERATE LEVITUS_GRID $SECONDS"

SECONDS=0
make_solo_mosaic \
	--num_tiles 1 \
	--dir ./ \
	--mosaic_name levitus_mosaic \
	--tile_file levitus_grid.nc \
	--periodx 360
echo "**** TEST06 SECONDS TO GENERATE LEVITUS_MOSAIC $SECONDS"

SECONDS=0
fregrid \
  --debug \
	--input_mosaic levitus_mosaic.nc \
	--input_file $dir_in/WOA09_ann_theta.nc \
	--scalar_field POTENTIAL_TEMP \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--extrapolate \
  --interp_method conserve_order1 \
	--dst_vgrid $dir_in/ocean_vgrid.nc \
	--check_conserve
echo "**** TEST06 SECONDS TO REMAP GENERATE TO END $SECONDS"

SECONDS=0
fregrid \
  --debug \
	--input_mosaic levitus_mosaic.nc \
	--input_file $dir_in/WOA09_ann_theta.nc \
	--scalar_field POTENTIAL_TEMP \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_file WOA09_ann_theta_cm2g_extrap1.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--output_mosaic $dir_in/ocean_mosaic.nc \
	--extrapolate \
  --interp_method conserve_order1 \
	--dst_vgrid $dir_in/ocean_vgrid.nc \
	--check_conserve
echo "**** TEST06 SECONDS TO REMAP READ TO END $SECONDS"
