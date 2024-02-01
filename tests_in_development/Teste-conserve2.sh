#!/usr/bin/bash

source run_tools.sh
echo "Teste second order conservative remapping from C1536 to lonlat 4608x2880"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Teste-input
dir_out=$PWD/tests_in_development/Teste-conserve1-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=""

scalar_field=""

input_mosaic=$(run_generate_cubed_sphere 1536)
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
         --interp_method conserve_order2 \
         --remap_file remap" >> LOG 2>&1
  echo "=> SECONDS TO GENERATE REMAP $1 TO $2:  $SECONDS" > LOG 2>&1
