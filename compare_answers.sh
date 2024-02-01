#!/bin/bash

ref_dir="/home/Mikyung.Lee/FRE-NCTools/test-benchmark/t"
test_dir=$PWD/t

for rdir in $ref_dir/*output ; do

  echo "***************************"
  tdir=${rdir/$ref_dir/$test_dir}

  for rfile in $rdir/*.nc ; do
    tfile=${rfile/$rdir/$tdir}
    set -x
    nccmp -Sfdm -c 10 -T 0.1 $rfile $tfile
    set +x
  done

  echo "***************************"

done
