#!/bin/bash

export PATH=$PATH:$PWD/frenctools/bin

t/Test03-grid_coupled_model.sh
echo "3**************************************************"
t/Test05-remap_c48_regular.sh
echo "5**************************************************"
t/Test06-regrid_extrap.sh
echo "6**************************************************"
#t/Test07-runoff_regrid.sh #takes a long time
t/Test10-remap_land_res.sh
echo "10*************************************************"
#t/Test13-make_quick_mosaic.sh #figure out output location
#echo "13**************************************************"
t/Test15-regrid_land.s
echo "15*************************************************"
t/Test16-river_regrid.sh
echo "16*************************************************"
t/Test20-fregrid.sh
echo "20*************************************************"
t/Test27-multiple_nests.sh
echo "27*************************************************"
#t/Test31-fregrid_stretched.sh #takes a long time
t/Test32-fregrid_no_stretched.sh
echo "31*************************************************"
t/Test33-reference_make_hgrid.sh
echo "33*************************************************"
