#!/usr/bin/bash

echo "Test make_hgrid do_transform vs. FV3 generated grid spec"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test33-input
dir_out=$PWD/Test33-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir -p $dir_out && cd $dir_out

gs=$dir_in/"fv3_48_grid_spec.tile"   #name prefix for the internal fv3 grid spec files.

#Ib) Call make_hgrid to generate the NCTools grids from FV3 files.
## TODO:  Note argument of --great_circle_algorithm option, which should not
## be necessary and code may be corrected in future.
make_hgrid --grid_type from_file \
	   --great_circle_algorithm \
	   --my_grid $gs"1.nc",$gs"2.nc",$gs"3.nc",$gs"4.nc",$gs"5.nc",$gs"6.nc" \
	   --grid_name  C48_ff_grid_spec

#II) "Analytically" create a equal distance gnomonic cubic grid using
#     the do_cube_transform option
make_hgrid \
    --grid_type gnomonic_ed \
    --nlon 96 \
    --do_cube_transform \
    --stretch_factor 3 \
    --target_lat 17.5 \
    --target_lon 172.5 \
    --grid_name C48_grid_spec
