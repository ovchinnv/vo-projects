ROOT=..
OBJS=multicom.o multicom_aux.o
NAME=multicom

LOCALMACROS=-include source.defs
LOCALDEPS=source.defs

INCLUDES=-I $(ROOT)/modules

include $(ROOT)/Makefile.inc

multicom.o:multicom_aux.o output.stamp vectors.stamp parser.stamp parselist.stamp
charmm: $(CHARMM_SOURCES)
