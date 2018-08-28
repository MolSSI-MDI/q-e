module unload gcc/7.2.0
module load intel/2017u4
module load intel-mpi/2017u4

cd ../

cp ~/qmmm/controller/qe/Modules/* Modules
cp ~/qmmm/controller/qe/PW/src/* PW/src
cp ~/qmmm/controller/qe/clib/* clib

make -j 32 pw
