
FC=gfortran
FFLAGS=-O3 -fbounds-check -frange-check -c
FC=ifort
FFLAGS= -O3 -c -extend-source -fast -fno-math-errno -finline -ipo
#
DEFS=-D float=real*8\
     -D int=integer\
     -D bool=logical\
     -D'__DEPMOD=use const, only: ERRTOL'\
     -D'__ERRTOL=ERRTOL()*50'\
     -D'__DEBUG=.false.'\
     -D'__CHKDET=.true.'\
     -D'__DETFAILTRACE=.true.'\

MACROS=-D'__WARN(__MSG)=write(0,*) __MSG' \
       -D'__DEBUG_OUT(__FMT,...)=write(666,__FMT) __VA_ARGS__'


LD=$(FC)
LDFLAGS=

FPP=cpp
FPPFLAGS=-C -w

EXE=bfit
OBJS=bestfit.o const.o bfit.o
SOURCES=bestfit.F90 const.F90

CHARMM_SOURCES=bestfit.src


default: all
all: $(EXE)
$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)

%.o:%.F90
	$(FC) $(DEFS) $(FFLAGS) $<
%.F90:%.ftn
	$(FPP) $(DEFS) $(MACROS) $(FPPFLAGS) $< > $@
%.src:%.ftn
	$(FPP) $(FPPFLAGS) -D__CHARMM -P $< > $@

bfit.o: bestfit.o const.o

bestfit.o : bestfit.F90 const.o
bestfit.F90 : bestfit.ftn source.defs
const.o :  const.F90
const.F90 : const.ftn

charmm: $(CHARMM_SOURCES)
bestfit.src: bestfit.ftn source.defs

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
