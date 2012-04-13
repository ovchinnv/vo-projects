ROOT=../..
CHARMM_SOURCES=ftsm_var.src
include $(ROOT)/Makefile.inc

COMMON=..

#LOCALMACROS=-include source.defs
#LOCALDEPS=source.defs
LIBDEP=$(COMMON)/sm_common.o

OBJS=ftsm_var.o
SOURCES=ftsm_var.F90

default: all
all: $(OBJS)

charmm: $(CHARMM_SOURCES)

clean: 
	rm -f *.mod $(CHARMM_SOURCES) $(SOURCES) $(OBJS) $(EXE)
