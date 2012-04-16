ROOT=..
SOURCES=multicom.F90 multicom_aux.F90
OBJS=multicom.o multicom_aux.o
CHARMM_SOURCES=multicom.src multicom_aux.src

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

default: all
all: $(OBJS)
$(OBJS):$(SOURCES)
charmm: $(CHARMM_SOURCES)
source: $(SOURCES)
clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
