ROOT=..
OBJS=bestfit.o
CHARMM_SOURCES=bestfit.src

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

# additional objects for linnking
LINKOBJECTS=constants.o

include $(ROOT)/Makefile.inc

EXE=bfit

default: all $(GLOBALDEPS) $(LOCALDEPS)
all: $(EXE)
$(EXE): $(OBJS) bfit.o
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) bfit.o -o $(EXE)

bfit.o    : bestfit.o constants
bestfit.o : constants

charmm: $(CHARMM_SOURCES)

lib: force_look
	/bin/sh -c 'for file in "bestfit.o"; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE) bfit.o
