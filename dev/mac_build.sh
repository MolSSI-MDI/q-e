cd ../

# copy the MDI modules out of the subtree and into the appropriate directories
cp mdi/molssi_driver_interface/mdi.c clib
cp mdi/molssi_driver_interface/mdi.h clib
cp mdi/molssi_driver_interface/mdi_f90.f90 Modules

#./configure CC=icc CXX=icpc F77=ifort F90=ifort FC=ifort CPP=cpp MPIF90=mpiifort -enable-parallel -enable-openmp
./configure
make -j 6 pw
