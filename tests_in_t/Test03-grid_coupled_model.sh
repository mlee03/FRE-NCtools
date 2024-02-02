#!/usr/bin/bash

# Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)
# set setup to generate ncls from test's input directory
#Combined with Test09

echo "Test grid for coupled model (land and atmosphere are C48 and ocean is 1 degree tripolar grid)"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test03-input
dir_out=$PWD/Test03-output

mkdir -p $dir_out && cd $dir_out

#Make_hgrid: create ocean_hgrid"
make_hgrid \
	--grid_type tripolar_grid \
	--nxbnd 2 \
	--nybnd 7 \
	--xbnd -280,80 \
	--ybnd -82,-30,-10,0,10,30,90 \
	--dlon 1.0,1.0 \
	--dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 \
	--grid_name ocean_hgrid \
	--center c_cell

#Make_vgrid: create ocean_vgrid
make_vgrid \
	--nbnds 3 \
	--bnds 0.,220.,5500. \
	--dbnds 10.,10.,367.14286 \
	--center c_cell \
  --grid_name ocean_vgrid

#Make_solo_mosaic: create ocean solo mosaic
make_solo_mosaic \
	--num_tiles 1 \
	--dir ./       \
	--mosaic_name ocean_mosaic \
	--tile_file ocean_hgrid.nc \
	--periodx 360

#Make_topog: create ocean topography data
make_topog \
	--mosaic ocean_mosaic.nc \
	--topog_type realistic \
	--topog_file $dir_in/OCCAM_p5degree.nc \
	--topog_field TOPO \
	--scale_factor -1 \
	--vgrid ocean_vgrid.nc \
	--output ocean_topog.nc
	--output ocean_topog.nc

check_mask \
  --grid_file ocean_mosaic.nc \
  --ocean_topog ocean_topog.nc \
  --layout 45,72

# MPI only tests
mpirun -n 4 make_topog_parallel \
		   --mosaic ocean_mosaic.nc \
		   --topog_type realistic \
		   --topog_file $dir_in/OCCAM_p5degree.nc \
		   --topog_field TOPO \
		   --scale_factor -1 \
		   --vgrid ocean_vgrid.nc \
		   --output ocean_topog_parallel.nc
		   --output ocean_topog_parallel.nc

check_mask \
  --grid_file ocean_mosaic.nc \
  --ocean_topog ocean_topog_parallel.nc \
  --layout 45,72

#Make_hgrid: create C48 grid for atmos/land
make_hgrid \
	--grid_type gnomonic_ed \
	--nlon 96 \
	--grid_name C48_grid

#Make_solo_mosaic: create C48 solo mosaic for atmos/land
make_solo_mosaic \
	--num_tiles 6 \
	--dir ./ \
	--mosaic C48_mosaic \
	--tile_file C48_grid.tile1.nc,C48_grid.tile2.nc,C48_grid.tile3.nc,C48_grid.tile4.nc,C48_grid.tile5.nc,C48_grid.tile6.nc

#Make_coupler_mosaic: coupler_mosaic with and without parallel and compare
make_coupler_mosaic \
	--atmos_mosaic C48_mosaic.nc \
	--ocean_mosaic ocean_mosaic.nc \
	--ocean_topog  ocean_topog.nc \
	--check \
	--area_ratio_thresh 1.e-10 \
	--mosaic_name grid_spec

mpirun -n 4 make_coupler_mosaic_parallel \
		   --atmos_mosaic C48_mosaic.nc \
		   --ocean_mosaic ocean_mosaic.nc \
		   --ocean_topog  ocean_topog_parallel.nc \
       --area_ratio_thresh 1.e-10 \
		   --mosaic_name grid_spec_parallel
