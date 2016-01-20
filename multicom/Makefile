ROOT=..
OBJS=multicom.o multicom_aux.o ifstack.o
NAME=multicom
DEPS=parselist mpi_stub
#parser: output
#parselist : parser vectors

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

multicom.o:multicom_aux.o
charmm:: $(CHARMM_SOURCES)
