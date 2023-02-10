# Makefile for library coupling interface
include ../../make.inc

# location of required modules and include files
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../build/mdi_build-prefix/src/mdi_build-build/MDI_Library

TLDEPS=pw

########################################################################

all : tldeps libqemdi.so qemdi.x

qemdi.x : libqemdi.so qemdi.o
	$(LD) $(LDFLAGS) -o $@ \
	   qemdi.o libqemdi.so $(MODFLAGS) $(LIBOBJS) $(QELIBS) \
        -Wl,--no-whole-archive \
        -L../build/mdi_build-prefix/src/mdi_build-build/MDI_Library -lmdi \
        -Wl,-rpath,../build/mdi_build-prefix/src/mdi_build-build/MDI_Library \

libqemdi.so: libqemdi.o
	$(MPIF90) -shared -o libqemdi.so libqemdi.o \
        -Wl,--whole-archive \
        ../../UtilXlib/libutil.a \
        ../../FFTXlib/src/libqefft.a \
        ../../LAXlib/libqela.a \
        ../../dft-d3/libdftd3qe.a \
        ../../KS_Solvers/libks_solvers.a \
        ../../Modules/libqemod.a \
        ../../PW/src/libpw.a \
        ../../MBD/libmbd.a \
        ../../upflib/libupf.a \
        ../../XClib/xc_lib.a \
        -Wl,--no-whole-archive \
        -L../build/mdi_build-prefix/src/mdi_build-build/MDI_Library -lmdi \
        -Wl,-rpath,../build/mdi_build-prefix/src/mdi_build-build/MDI_Library \
	$(QELIBS) \
        $(F90FLAGS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L *.so

# we currently need to specify explicit dependencies here
libqemdi.o: libqemdi.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../UtilXlib/parallel_include.o

qemdi.o: qemdi.f90 libqemdi.so

# DO NOT DELETE
