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

# Test regrid land data with cell_measures and cell_methods attribute

echo "Test regrid land data"

dir_in=$PWD/t/Test15-input
dir_out=$PWD/t/Test15-output
mkdir -p $dir_out
cd $dir_out

   fregrid \
		--input_mosaic $dir_in/C180_mosaic.nc \
		--interp_method conserve_order1 \
		--nlon 144 \
		--nlat 90 \
		--remap_file remap_file.nc

#remap static field
   fregrid  \
		--input_mosaic $dir_in/C180_mosaic.nc  \
		--interp_method conserve_order1  \
		--nlon 144  \
		--nlat 90  \
		--input_file $dir_in/00050101.land_static  \
		--scalar_field soil_frac,lake_frac,glac_frac,area,soil_area,lake_area,glac_area  \
		--output_file out.nc  \
		--remap_file remap_file.nc

# parallel call
     mpirun -n 4 fregrid_parallel \
		--input_mosaic $dir_in/C180_mosaic.nc \
		--interp_method conserve_order1 \
		--nlon 144 \
		--nlat 90 \
		--remap_file remap_file4.nc

     mpirun -n 4 fregrid_parallel  \
		--input_mosaic $dir_in/C180_mosaic.nc  \
		--interp_method conserve_order1  \
		--nlon 144  \
		--nlat 90  \
		--input_file $dir_in/00050101.land_static  \
		--scalar_field soil_frac,lake_frac,glac_frac,area,soil_area,lake_area,glac_area  \
		--output_file out_parallel.nc  \
		--remap_file remap_file4.nc

# remap other fields
# Commented this part out because the input file is too large
#   fregrid  \
#		--input_mosaic C180_mosaic.nc  \
#		--interp_method conserve_order1  \
#		--nlon 144  \
#		--nlat 90  \
#		--input_file 00050101.land_month  \
#		--scalar_field evap_land,evap_soil,evap_glac,evap_lake  \
#		--output_file 00050101.land_month.nc  \
#		--remap_file remap_file.nc
