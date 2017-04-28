NAME=string/openmm
ROOT=../..
COMMON=..
OBJS=sm_plugin.o
HEADS=sm.h
DEPS=string dynamol mpi_stub

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include Makefile.inc

sm_plugin.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod ../atomlist.def

clean::
	rm -f $(OBJS) sm.h
