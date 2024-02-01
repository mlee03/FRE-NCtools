#!/usr/bin/bash

source run_tools.sh
echo "Testm first order conservative remapping from 1/8 tripolar to lonlat 360x180"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testm-input
dir_out=$PWD/tests_in_development/Testm-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=""

scalar_field=""

input_mosaic=$dir_in/$input_mosaic
output_mosaic=$(run_generate_lonlat 4608 2880)

  echo "=========================================" >> LOG
  echo "`date` GENERATING $5" >> LOG
  echo "=========================================" >> LOG
  SECONDS=0
  eval "fregrid \
         --debug \
         --input_mosaic $input_mosaic \
         --output_mosaic $output_mosaic \
         --input_dir $dir_in \
         --output_dir $dir_out \
         --interp_method conserve_order1 \
         --remap_file remap" >> LOG 2>&1
  echo "=> SECONDS TO GENERATE REMAP $1 TO $2:  $SECONDS" > LOG 2>&1
