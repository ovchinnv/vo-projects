ROOT=..
DYNA=$(ROOT)/dynamol
OBJS=chirality.o confcons.o
NAME=coonfcons
DEPS=constants output multicom parser dynamol

LOCALMACROS=-include source.defs -include $(ROOT)/dmol.defs
LOCALDEPS=source.defs $(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc

charmm: $(CHARMM_SOURCES)
