#!/bin/bash

ref_dir="/home/Mikyung.Lee/FRE-NCTools/test-benchmark/t"
test_dir=$PWD/t
log='NCCMP_RESULTS'

if [ -f $log ] ; then ; rm $log ; fi
touch -a $log

for rdir in $ref_dir/*-output ; do

  echo "***************************"
  tdir=${rdir/$ref_dir/$test_dir}
  echo $tdir

  for rfile in $tdir ; do
    tfile=${rfile/$rdir/$tdir}
    nccmp -Sfdm -T 0.1 $rfile $tfile > $log 2>&1 &
  done

  echo "***************************"

done
