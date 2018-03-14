NAME=string/acemd
ROOT=../..
COMMON=..
OBJS=sm_plugin_acemd.o sm_acemd.o
HEADS=sm.h
DEPS=string dynamol mpi_stub

LOCALDEFS+=-D__BINDC

include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc
include MakeACEMD.inc

sm_acemd.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod ../atomlist.def

sm_plugin_acemd.o : sm_plugin_acemd.c aceplug.h sm.h $(LOCALDEPS)
	$(CC) $(CFLAGS) $(LOCALDEFS) $(LOCALMACROS) $< -o $@

aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h

clean::
	rm -f $(OBJS) sm.h
