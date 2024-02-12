#!/bin/bash

install_dir="/home/Mikyung.Lee/FRE-NCTools/test-benchmark/fre-nctools/"
my_bin=$install_dir/"bin"
export PATH=$PATH:$my_bin

function run_test() {
  logfile=${1/.sh/.LOG}
  echo "`date` $1 USING $my_bin" > $logfile
  echo "`tail -n 1 $install_dir/COMPILE_HISTORY`" >> $logfile
  eval "./$1" >> $logfile 2>&1
}

if [ $# -eq 0 ] ; then
  for itest in ./Test*.sh ; do run_test $itest ; done
else
  run_test $1
fi
