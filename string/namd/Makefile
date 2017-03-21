NAME=string/namd
ROOT=../..
COMMON=..
OBJS=sm_namd.o
HEADS=GlobalMasterSM.h
DEPS=string dynamol mpi_stub

LOCALDEFS+=-D__BINDC -D__CPLUSPLUS


include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs ../atomlist.def
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include MakeNAMD.inc

sm_namd.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod

clean::
	rm -f $(OBJS) $(HEADS)
