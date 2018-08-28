module unload gcc/7.2.0
module load intel/2017u4
module load intel-mpi/2017u4

cd ../

./configure CC=icc CXX=icpc F77=ifort F90=ifort FC=ifort CPP=cpp MPIF90=mpiifort -enable-parallel -enable-openmp
make -j 32 pw
