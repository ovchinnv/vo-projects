ROOT=..
OBJS=bestfit.o
NAME=bestfit
DEPS=constants output
#
CBINDOBJS=bestfit_c.o
BINDSOURCES=$(OBJS:%.o=%_.F90) $(CBINDOBJS:%.o=%_.F90)
SBINDOBJS  =$(OBJS:%.o=%_.o)   $(CBINDOBJS:%.o=%_.o) 
DBINDOBJS  =$(OBJS:%.o=%__.o)  $(CBINDOBJS:%.o=%__.o)
CLIB=-lbestfit
#
LOCALMACROS=-include source.defs
LOCALDEPS=source.defs
# additional objects for linking test.o executable

LINKOBJS=constants.o output.o

include $(ROOT)/Makefile.inc

SWIGC=swig
SFLAGS=

charmm: $(CHARMM_SOURCES)

%_.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/source.msg
	$(FPP) $(FPPFLAGS) $(DEFS) -D__BINDC $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@

libraries : cbind tcl
cbind : bestfit.h static shared
shared : dynamic
dynamic : libbestfit.so
static : libbestfit.a
libbestfit.a: $(SBINDOBJS)
	$(AR) $(AFLAGS) $@ $(SBINDOBJS)
libbestfit.so: $(DBINDOBJS)
	$(CC) -shared $(DBINDOBJS) $(FLIBS) -o $@

tcl : libbestfit_tcl.so
libbestfit_tcl.so : bestfit.i $(DBINDOBJS)
	$(SWIGC) -D__BINDC $(SFLAGS) -tcl bestfit.i
	$(CC) -fPIC $(CFLAGS) -D__BINDC $(LOCALMACROS) bestfit_wrap.c
	$(CC) -shared bestfit_wrap.o $(DBINDOBJS) $(FLIBS) -o libbestfit_tcl.so
bestfit.i : bestfit.h
bestfit.h : bestfit.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) -D__BINDC -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@

ctest :: ctest.o static force_look
	$(LD) $(LDFLAGS) -Bstatic ctest.o $(LIBS) $(MPI_LIBS) $(CLIB) -o  ctest
	echo "==================================================="
	echo "Running C test :"
	echo "==================================================="
	time ./ctest > ctest.out 2>&1 ; cat ctest.out
	echo "==================================================="
	echo
ctest.o : ctest.c bestfit.h $(LOCALDEPS)
	$(CC) $(CFLAGS) -D__BINDC $(LOCALMACROS) -P $< -o $@

clean::
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(BINDSOURCES) $(OBJS) $(CBINDOBJ) $(SBINDOBJS) $(DBINDOBJS) $(ROOT)/$(NAME).stamp \
bestfit_c.F90 ctest.o fort.* ctest *_wrap* bestfit.h

