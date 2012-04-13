ROOT=..
CHARMM_SOURCES=bestfit.src
include $(ROOT)/Makefile.inc

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

EXE=bfit
OBJS=bestfit.o const.o bfit.o
SOURCES=bestfit.F90 const.F90

default: all $(GLOBALDEPS) $(LOCALDEPS)
all: $(EXE)
$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.F90:%.ftn
	$(FPP) $(FPPFLAGS) $(DEFS) -U__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@
%.src:%.ftn
	$(FPP) $(FPPFLAGS) -D__CHARMM $(GLOBALMACROS) $(LOCALMACROS) $(MSG) -P $< > $@

bfit.o: bestfit.o const.o

bestfit.o : bestfit.F90 const.o
bestfit.F90 : bestfit.ftn source.defs
const.o :  const.F90
const.F90 : const.ftn

charmm: $(CHARMM_SOURCES) $(GLOBALDEPS) $(LOCALDEPS)
bestfit.src: bestfit.ftn 

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
