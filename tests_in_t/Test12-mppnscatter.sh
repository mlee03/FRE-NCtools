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

# Test procedure for mppncscatter is to first split the files.  Once scattered,
# run mppnccombine to re-combine the files.  The final, recombined file should
# match the unscattered file.

# The mppnccombine and mppncscatter commands should probably be tested in
# the same file, since here we assume mppnccombine is running correctly.

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Test12-input
dir_out=$PWD/tests_in_t/Test12-output

rm -rf $dir_out
mkdir $dir_out
cd $dir_out

test_file=fv_core.res.tile1.nc

SECONDS=0

# Scatter the file
mppncscatter -i 2 -j 3 -x 2 -y 12 $dir_in/$test_file

# Combine the file:
mppnccombine -64 $dir_out/$test_file ${test_file}.????

echo "**** TEST12 SECONDS TO RUN TEST12"
