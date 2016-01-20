ROOT=..
DYNA=$(ROOT)/dynamol
OBJS=chirality.o confcons.o
NAME=confcons
DEPS=dynamol
#multicom: output parser
#dynamol: multicom
#dynamol: constants

# NOTE : use += vs. = to avoid including in charmm distrib per charmm_distrib script in root (a trick)
LOCALMACROS+=-include $(ROOT)/source.defs
LOCALMACROS=-include source.defs
LOCALMACROS+=-include $(ROOT)/dmol.defs

LOCALDEPS=source.defs
LOCALDEPS+=$(ROOT)/dmol.defs

include $(ROOT)/Makefile.inc

charmm:: $(CHARMM_SOURCES)
