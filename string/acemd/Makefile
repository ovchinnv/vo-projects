NAME=string/acemd
ROOT=../..
COMMON=..
OBJS=plugin-smcv.o smcv_acemd.o
HEADS=smcv.h
DEPS=string dynamol mpi_compat

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include MakeACEMD.inc

smcv_acemd.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod

plugin-smcv.o : plugin-smcv.c aceplug.h smcv.h $(LOCALDEPS)
	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

clean::
	rm -f $(OBJS) smcv.h
