cd ../

rm -r mdi/build
#./configure CC=icc CXX=icpc F77=ifort F90=ifort FC=ifort CPP=cpp MPIF90=mpiifort -enable-parallel -enable-openmp
./configure
make -j 6 pw
