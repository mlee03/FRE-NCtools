#!/bin/bash


function run_generate_cubed_sphere() {

  grid_name=C$1
  mosaic_name=C$1_mosaic

  echo "=========================================" >> LOG
  echo "`date` GENERATING $mosaic_name.nc" >> LOG
  echo "=========================================" >> LOG
  eval "make_hgrid \
        --grid_type gnomonic_ed \
        --nlon $(( 2 *$1 )) \
        --grid_name $grid_name " >> LOG 2>&1

  grid_name=$grid_name.tile

  eval "make_solo_mosaic \
         --num_tiles 6 \
         --dir ./ \
         --mosaic_name $mosaic_name \
         --tile_file ${grid_name}1.nc,${grid_name}2.nc,${grid_name}3.nc,${grid_name}4.nc,${grid_name}5.nc,${grid_name}6.nc" >> LOG 2>&1

  echo $mosaic_name.nc

}


function run_generate_lonlat() {

  grid_name=lonlat_$1x$2
  mosaic_name=lonlat_$1x$2_mosaic

  echo "=========================================" >> LOG
  echo "`date` GENERATING $mosaic_name.nc" >> LOG
  echo "=========================================" >> LOG
  eval "make_hgrid \
         --grid_type regular_lonlat_grid \
         --nxbnd 2 \
         --nybnd 2 \
         --xbnd 0,360 \
         --ybnd -90,90 \
         --nlon $(( 2 * $1 )) \
         --nlat $(( 2 * $2)) \
         --grid_name $grid_name" >> LOG 2>&1

  eval "make_solo_mosaic \
         --num_tiles 1 \
         --dir ./ \
         --mosaic_name $mosaic_name \
        --tile_file $grid_name.nc" >> LOG 2>&1

  echo $mosaic_name.nc

}

function run_fregrid_lonlat_conserve1(){

  echo "=========================================" >> LOG
  echo "`date` GENERATING $5" >> LOG
  echo "=========================================" >> LOG
  SECONDS=0
  eval "fregrid \
         --debug \
         --input_mosaic $1 \
         --output_mosaic $2 \
         --input_dir $3 \
         --output_dir $4 \
         --associated_file_dir $3 \
         --input_file $5 \
         --output_file $5 \
         --scalar_field $6 \
         --interp_method conserve_order1 \
         --remap_file remap" >> LOG 2>&1
  echo "=> SECONDS TO REMAP $1 TO $2:  $SECONDS" > LOG 2>&1

}


function run_fregrid_lonlat_conserve1(){

  echo "=========================================" >> LOG
  echo "`date` GENERATING $5" >> LOG
  echo "=========================================" >> LOG
  SECONDS=0
  eval "fregrid \
         --debug \
         --input_mosaic $1 \
         --output_mosaic $2 \
         --input_dir $3 \
         --output_dir $4 \
         --associated_file_dir $3 \
         --input_file $5 \
         --output_file $5 \
         --scalar_field $6 \
         --interp_method conserve_order2 \
         --remap_file remap" >> LOG 2>&1
  echo "=> SECONDS TO REMAP $1 TO $2:  $SECONDS" > LOG 2>&1

}
