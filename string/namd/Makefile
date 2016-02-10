NAME=string/namd
ROOT=../..
COMMON=..
OBJS=smcv_namd.o
HEADS=GlobalMasterSMCV.h
DEPS=string dynamol mpi_compat

LOCALDEFS+=-D__BINDC -D__CPLUSPLUS

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include MakeNAMD.inc

smcv_namd.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod

clean::
	rm -f $(OBJS) smcv.h
