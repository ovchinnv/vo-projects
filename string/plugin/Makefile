NAME=string/plugin
ROOT=../..
COMMON=..
OBJS=sm_plugin.o sm_plugin_acemd.o
HEADS=struna.h
DEPS=dynamol mpi_stub

LIB=libstruna.so libstruna.a
LIBS+=$(LIB)

LOCALDEFS+=-D__BINDC

default : plugin
plugin : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS) linkobjs.inc

include linkobjs.inc
include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc

sm_plugin.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod ../atomlist.def
sm_plugin_acemd.o : aceplug.h

aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h

clean::
	rm -f $(OBJS) $(HEADS)
