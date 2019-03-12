ROOT=..
OBJS=bestfit.o
NAME=bestfit
DEPS=constants multicom
#parselist : parser vectors
#multicom : parselist
#parser: output
#
CBINDOBJS=bestfit_c.o
BINDSOURCES=$(OBJS:%.o=%_.F90) $(CBINDOBJS:%.o=%_.F90)
SBINDOBJS  =$(OBJS:%.o=%_.o)   $(CBINDOBJS:%.o=%_.o) 
DBINDOBJS  =$(OBJS:%.o=%__.o)  $(CBINDOBJS:%.o=%__.o)
CLIB=-L. -lbestfit
#
LOCALMACROS=-include source.defs
LOCALDEPS=source.defs eig3s_lapack.def
#LOCALDEFS+=-fno-fast-math
# additional objects for linking test.o executable

LINKOBJS=constants.o output.o parser.o cmd.o prm.o multicom_aux.o multicom.o ivector.o parselist.o isort.o
TESTOBJS=$(OBJS)
PYTHONH=-I /usr/include/python2.7/
LPYTHON=-lpython2.7

include $(ROOT)/Makefile.inc

SWIGC=swig
SFLAGS=

charmm:: $(CHARMM_SOURCES)

%_.F90:%.ftn $(LOCALDEPS) $(GLOBALDEPS) $(ROOT)/source.msg
	$(FPP) $(FPPFLAGS) $(DEFS) -D__BINDC $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@

libraries : cbind #tcl python
cbind : bestfit.h static shared
shared : dynamic
dynamic : libbestfit.so
static : libbestfit.a
libbestfit.a: $(SBINDOBJS)
	$(AR) $(AFLAGS) $@ $(SBINDOBJS)
libbestfit.so: $(DBINDOBJS)
	$(LD) -shared $(DBINDOBJS) $(FLIBS) -o $@

tcl : libbestfit_tcl.so
libbestfit_tcl.so : bestfit.i $(DBINDOBJS)
	$(SWIGC) -D__BINDC $(SFLAGS) -tcl bestfit.i
	$(CC) $(PIC) $(CFLAGS) -D__BINDC $(LOCALMACROS) bestfit_wrap.c
	$(CC) -shared bestfit_wrap.o $(DBINDOBJS) $(FLIBS) -o $@

python : _bestfit.so
_bestfit.so : bestfit.i $(DBINDOBJS)
	$(SWIGC) -D__BINDC $(SFLAGS) -python bestfit.i
	$(CC) $(PIC) $(CFLAGS) -D__BINDC $(LOCALMACROS) $(PYTHONH) bestfit_wrap.c
	$(CC) -shared bestfit_wrap.o $(DBINDOBJS) $(FLIBS) $(LPYTHON) -o $@

bestfit.i : bestfit.h
bestfit.h : bestfit.H $(LOCALDEPS) $(GLOBALDEPS)
	$(FPP) $(FPPFLAGS) -D__BINDC -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) -P $< > $@

ctest :: ctest.o shared static force_look
	$(CC) -Bdynamic ctest.o $(LIBS) $(MPI_LIBS) $(CLIB) $(FLIBS) -o  ctest
	echo "==================================================="
	echo "Running C test :"
	echo "==================================================="
	time ./ctest > ctest.out 2>&1 ; cat ctest.out
	echo "==================================================="
	echo
ctest.o : ctest.c bestfit.h $(LOCALDEPS)
	$(CC) $(CFLAGS) -D__BINDC $(LOCALMACROS) $< -o $@

clean::
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(BINDSOURCES) $(OBJS) $(CBINDOBJ) $(SBINDOBJS) $(DBINDOBJS) $(ROOT)/$(NAME).stamp \
bestfit_c.F90 ctest.o fort.* ctest *_wrap* bestfit.h bestfit.py bestfit.pyc

