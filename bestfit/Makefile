ROOT=..
OBJS=bestfit.o
NAME=bestfit
DEPS=constants output

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

# additional objects for linking executable
LINKOBJS=constants.o output.o

include $(ROOT)/Makefile.inc

EXE=test


test: $(EXE)
$(EXE): $(OBJS) bfit.o
	$(LD) $(LDFLAGS) $(OBJS) $(LOBJS) bfit.o -o $(EXE)

bfit.o    : bestfit.o

charmm: $(CHARMM_SOURCES)

#lib: force_look
#	/bin/sh -c 'for file in "bestfit.o"; do  ln -fs `pwd`/$$file ${ROOT}/lib/$$file; done'

clean::
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE) bfit.o $(ROOT)/$(NAME).stamp
