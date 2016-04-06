NAME=string/acemd
ROOT=../..
COMMON=..
OBJS=plugin-smcv.o smcv_acemd.o plugin-ftsm.o ftsm_acemd.o
HEADS=smcv.h ftsm.h
DEPS=string dynamol mpi_stub

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include MakeACEMD.inc

smcv_acemd.o ftsm_acemd.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod

plugin-smcv.o : plugin-smcv.c aceplug.h smcv.h $(LOCALDEPS)
	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

plugin-ftsm.o : plugin-ftsm.c aceplug.h ftsm.h $(LOCALDEPS)
	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

clean::
	rm -f $(OBJS) smcv.h ftsm.h
