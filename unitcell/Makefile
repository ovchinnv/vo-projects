ROOT=..
OBJS=unitcell.o
NAME=unitcell
DEPS=output constants bestfit

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

LINKOBJS=output.o

include $(ROOT)/Makefile.inc

charmm:: $(CHARMM_SOURCES)
clean ::
	rm -f _parser.tmp