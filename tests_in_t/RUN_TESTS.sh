#!/bin/bash


my_bin="/home/Mikyung.Lee/FRE-NCTools/test-benchmark/frenctools/bin"
export PATH=$PATH:$my_bin

function run_test() {
  echo "`date` $1 USING $my_bin"
  eval "./$1" > ${1/.sh/.LOG} 2>&1
}

if [ $? -eq 0 ] ; then
  for itest in ./Test*.sh ; do run_test $itest ; done
else
  run_test $1
fi
