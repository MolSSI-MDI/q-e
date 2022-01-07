# Makefile for library coupling interface
include ../../make.inc

# location of required modules and include files
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../CPV/src

TLDEPS=pw cp

COUPLEOBJS= \
libpwscf.o \
libplugin.o \
libcpv.o \
libqemod.o

########################################################################

all : tldeps libqecouple.a libqecouple.so

libqecouple.a : $(COUPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

libqecouple.so: $(COUPLEOBJS)
	$(MPIF90) -shared -o libqecouple.so $(COUPLEOBJS) $(F90FLAGS)

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
