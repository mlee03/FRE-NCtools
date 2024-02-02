#!/usr/bin/bash

# Test procedure for mppncscatter is to first split the files.  Once scattered,
# run mppnccombine to re-combine the files.  The final, recombined file should
# match the unscattered file.

# The mppnccombine and mppncscatter commands should probably be tested in
# the same file, since here we assume mppnccombine is running correctly.

echo "Test mppncombine and mppnscatter"
echo "COMPILED ON `tail -n 1 $my_bin/COMPILE_HISTORY`"
echo "=============================================="

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test12-input
dir_out=$PWD/Test12-output

[[ -d $dir_out ]] ; rm -rf $dir_out
mkdir $dir_out && cd $dir_out

test_file=fv_core.res.tile1.nc

# Scatter the file
mppncscatter -i 2 -j 3 -x 2 -y 12 $dir_in/$test_file

# Combine the file:
mppnccombine -64 $dir_out/$test_file ${test_file}.????
