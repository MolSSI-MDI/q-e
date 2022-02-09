# Makefile for library coupling interface
include ../../make.inc

# location of required modules and include files
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../CPV/src

TLDEPS=pw cp

COUPLEOBJS= \
libpwscf.o \
libcpv.o \
libqemod.o

########################################################################

all : tldeps libqecouple.a libplugin.so
#all : tldeps libplugin.so

libqecouple.a : $(COUPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

#libplugin.so: libplugin.o
#	$(MPIF90) -Wl,--whole-archive -shared -o libplugin.so libplugin.o $(F90FLAGS)

#libplugin.so: libplugin.o
#	$(MPIF90) -shared -o libplugin.so libplugin.o $(F90FLAGS) \
#        -L../../mdi/build/MDI_Library/ -Wl,-rpath,../../mdi/build/MDI_Library/ -lmdi \
#        -Wl,--whole-archive \
#        -L../../Modules -Wl,-rpath,/repo/build/q-e/Modules -lqemod \
#        -Wl,--no-whole-archive

libplugin.so: libplugin.o
	$(MPIF90) -shared -o libplugin.so libplugin.o \
        -Wl,--whole-archive \
        ../../Modules/libqemod.a \
        ../../UtilXlib/libutil.a \
        ../../FFTXlib/libqefft.a \
        ../../LAXlib/libqela.a \
        ../../dft-d3/libdftd3qe.a \
        ../../KS_Solvers/libks_solvers.a \
        ../../clib/clib.a \
        ../../PW/src/libpw.a \
        ../../FoX/lib/libFoX_common.a \
        ../../FoX/lib/libFoX_dom.a \
        ../../FoX/lib/libFoX_fsys.a \
        ../../FoX/lib/libFoX_sax.a \
        ../../FoX/lib/libFoX_utils.a \
        ../../FoX/lib/libFoX_wcml.a \
        ../../FoX/lib/libFoX_wkml.a \
        ../../FoX/lib/libFoX_wxml.a \
        -Wl,--no-whole-archive \
	$(QELIBS) \
        $(F90FLAGS)

#libplugin.so: libplugin.o
#	$(MPIF90) -Wl,--whole-archive -shared -o libplugin.so libplugin.o \
#        -fPIC -x f95-cpp-input -fopenmp -D__FFTW -D__MPI \
#        -I/repo/build/q-e//include \
#        -I/repo/build/q-e//FoX/finclude \
#        -I/repo/build/q-e//S3DE/iotk/include/ \
#        -I/repo/build/q-e//iotk/src -I/repo/build/q-e//Modules \
#        -I/repo/build/q-e//FFTXlib \
#        -I/repo/build/q-e//LAXlib \
#        -I/repo/build/q-e//UtilXlib \
#        -I/repo/build/q-e//FoX/finclude \
#        -I/repo/build/q-e//mdi/build/MDI_Library \
#        -I../../CPV/src

#libplugin.so: libplugin.o
#	$(MPIF90) -shared -o libplugin.so libplugin.o \
        -fPIC -x f95-cpp-input -fopenmp -D__FFTW -D__MPI

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

# we currently need to specify explicit dependencies here
libpwscf.o: libpwscf.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../UtilXlib/parallel_include.o \
	../../mdi/build/MDI_Library/libmdi.so

libplugin.o: libplugin.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../UtilXlib/parallel_include.o \
	../../mdi/build/MDI_Library/libmdi.so

libcpv.o: libcpv.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../Modules/mp_images.o ../../Modules/check_stop.o \
	../../Modules/io_global.o ../../UtilXlib/parallel_include.o \
	../../CPV/src/input.o \
	../../mdi/build/MDI_Library/libmdi.so

libqemod.o: libqemod.f90 ../../Modules/qmmm.o \
	../../mdi/build/MDI_Library/libmdi.so

# DO NOT DELETE
