ROOT=..
SOURCES=bestfit.F90 const.F90
OBJS=bestfit.o const.o bfit.o
CHARMM_SOURCES=bestfit.src

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

EXE=bfit

default: all $(GLOBALDEPS) $(LOCALDEPS)
all: $(EXE)
$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)

bfit.o    : bestfit.o const.o
bestfit.o : bestfit.F90 const.o

charmm: $(CHARMM_SOURCES)

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
