echo 'builindg FRE-NCtools conda package...'
echo "SRC_DIR / Build directory is: $PWD"
echo -e "Contents of SRC_DIR / Build directory are: \n`ls`"

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib
echo -e "PRE CONFIGURATION::\n\n"
echo -e "PATH is: \n $PATH"
#echo -e "LD_LIBRARY_PATH is: \n $LD_LIBRARY_PATH \n"


## this is sufficient
mkdir build && cd build
autoreconf -iv ../
../configure --prefix=$PREFIX --with-mpi --enable-quad-precision || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log

#export LD_LIBRARY_PATH=${PREFIX}/lib
echo -e "POST CONFIGURATION::\n\n"
echo -e "PATH is: \n $PATH"
#echo "LD_LIBRARY_PATH is:"
#echo $LD_LIBRARY_PATH
#echo ""

echo "compiling/building and installing into $PREFIX"
make install

echo "testing"
make check

echo "removing build directory"
cd ../ && rm -rf build

### to test, build-dir option, ala README
#autoreconf -i
#mkdir build && cd build
#../configure --prefix=$PREFIX || cat config.log
##../configure --prefix=$PREFIX --with-mpi || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.logecho "compiling/building"
#make
#echo "installing into $PREFIX"
#make install
