# Makefile for the MolSSI Driver Interface Library

include ../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS)

#OBJS = write_upf.o read_upf.o
#QEMODS = ../Modules/libqemod.a ../FFTXlib/libqefft.a ../UtilXlib/libutil.a

TLDEPS = 

all : libmdi.a

libmdi.a : 
	- mkdir build
	- cd build; \
          cmake -Dlibtype=SHARED -Dlanguage=Fortran -Dtest_codes=ON -D CMAKE_C_COMPILER=$(CC) \
	     -D CMAKE_Fortran_COMPILER=$(MPIF90) .. -D CMAKE_C_FLAGS="-fPIC" -D CMAKE_Fortran_FLAGS="-fPIC" ; \
          make

#tldeps:
#	if test -n "$(TLDEPS)" ; then \
#	( cd ../ ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f  *.x *.o *~ *_tmp.f90 *.mod *.d *.i *.L

#include make.depend
