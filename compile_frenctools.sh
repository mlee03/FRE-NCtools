#!/bin/bash

set -x

curr_dir=$PWD
build_dir=./build
if [ -z $2 ] ; then
  install_dir=$PWD/fre-nctools
else
  install_dir=$2
fi

if [ $1 == 'cpu' ] ; then
  export CC=mpicc
  export FC=mpif90
  export FCFLAGS="-g -I`nc-config --includedir` -Duse_libMPI -Wall -Minfo=all -O3"
  export CFLAGS="-g -I`nc-config --includedir` -std=c11 -Duse_libMPI -Wall -Minfo=all -O3"
  export LDFLAGS="`nc-config --libs`"
elif [ $1 == 'gpu' ] ; then
  export CC=nvc
  export FC=nvfortran
  export FCFLAGS="-g -O0 -I`nc-config --includedir` -acc -D_OPENACC "
  export CFLAGS="-g -O0 -I`nc-config --includedir` -std=c11 -Minfo=accel -acc -D_OPENACC"
  export LDFLAGS="`nc-config --libs`"
else
  export CC=nvc
  export FC=nvfortran
  export FCFLAGS="-g -O0 -I`nc-config --includedir` -acc -D_OPENACC "
  export CFLAGS="-g -O0 -I`nc-config --includedir` -std=c11 -Minfo=accel -acc -D_OPENACC"
  export LDFLAGS="`nc-config --libs`"
fi

if [ -d $build_dir ] ; then rm -rf $build_dir ; fi
if [ -d $install_dir ] ; then rm -rf $install_dir ; fi
mkdir $install_dir

autoreconf -iv
mkdir $build_dir && cd $build_dir
$curr_dir/configure --prefix=$install_dir --with-mpi
echo "COMPILING" && make > COMPILE_LOG 2>&1
make install

echo "FRESH COMPILE ON \`date\`" >> $install_dir/COMPILE_HISTORY

cat > compileme.sh <<EOF
 #!/bin/bash
 make install
 echo "RECOMPILE ON `date` >> $install_dir/COMPILE_HISTORY
EOF
chmod +x compileme.sh
