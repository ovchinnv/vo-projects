ROOT=..
OBJS=bestfit.o
NAME=bestfit
DEPS=constants output

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

# additional objects for linking executable
LINKOBJS=constants.o output.o

include $(ROOT)/Makefile.inc

charmm: $(CHARMM_SOURCES)

#lib: force_look
#	/bin/sh -c 'for file in "bestfit.o"; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'

cbind : bestfit_c.o
bestfit_c.F90:bestfit_c.ftn $(LOCALDEPS)
	$(FPP) $(FPPFLAGS) $(DEFS) -D__BINDC $(LOCALMACROS) -P $< > $@
ctest : bestfit_c.o cbindtest.o
	$(LD) $(LDFLAGS) bestfit_c.o cbindtest.o $(OBJS) $(LOBJS) $(LIBS) $(MPI_LIBS) -o  ctest
cbindtest.o : cbindtest.c bestfit.h $(LOCALDEPS)
	$(CC) $(CFLAGS) -D__BINDC $(LOCALMACROS) -P $< -o $@
bestfit_c.o:bestfit_c.F90 bestfit.o
	$(FC) $(FFLAGS) $(DEFS) $(INCLUDES)  $< -o $@

clean::
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(ROOT)/$(NAME).stamp bestfit_c.F90 bestfit_c.o cbindtest.o fort.666
