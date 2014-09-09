ROOT=..
OBJS=multicom.o multicom_aux.o ifstack.o
NAME=multicom
DEPS=vectors output parser parselist

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

include $(ROOT)/Makefile.inc

multicom.o:multicom_aux.o
charmm:: $(CHARMM_SOURCES)
