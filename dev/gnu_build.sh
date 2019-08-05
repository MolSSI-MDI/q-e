#module unload gcc/7.2.0
#module load intel/2017u4
#module load intel-mpi/2017u4

# This build script is intended to work in the following environment:
# conda install -c conda-forge mpi4py
# conda install -c conda-forge fortran-compiler

cd ../

rm -r mdi/build
#export MPICH_F90=gfortran
./configure CC=gcc CXX=c++ F77=gfortran F90=g95 FC=gfortran CPP=cpp MPIF90=mpif90 -enable-parallel -enable-openmp
make -j 32 pw
