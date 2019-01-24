module unload gcc/7.2.0
module load intel/2017u4
module load intel-mpi/2017u4

cd ../

rm -r mdi/build
make -j 32 pw
