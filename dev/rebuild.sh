module unload gcc/7.2.0
module load intel/2017u4
module load intel-mpi/2017u4

cd ../

make -j 32 pw
