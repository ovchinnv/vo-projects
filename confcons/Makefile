ROOT=..
DYNA=$(ROOT)/dynamol
OBJS=chirality.o confcons.o
NAME=confcons
DEPS=constants output multicom parser dynamol

LOCALMACROS=-include $(ROOT)/dmol.defs
LOCALDEPS=$(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc

charmm: $(CHARMM_SOURCES)
