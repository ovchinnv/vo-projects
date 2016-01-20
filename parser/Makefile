ROOT=..
OBJS=cmd.o prm.o parser.o
NAME=parser
DEPS=output constants mpi_stub

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

LINKOBJS=output.o

include $(ROOT)/Makefile.inc

parser.o : cmd.o prm.o
charmm:: $(CHARMM_SOURCES)
clean ::
	rm -f _parser.tmp