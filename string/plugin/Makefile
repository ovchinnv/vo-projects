NAME=string/plugin
ROOT=../..
COMMON=..
OBJS=sm_plugin.o
HEADS=struna.h struna.hpp
DEPS=dynamol mpi_stub unitcell

LIB=libstruna.so libstruna.a
LIBS+=$(LIB)

LOCALDEFS+=-D__BINDC

include linkobjs.inc
include $(ROOT)/Makefile.inc
ifeq ($(strip($ACEMD_PLUGIN)),1)
 OBJS+=sm_plugin_acemd.o
endif


default : plugin
plugin : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS) linkobjs.inc

LOCALMACROS+=-include source.defs
LOCALMACROS+=-include $(COMMON)/source.defs
LOCALDEPS+=source.defs
LOCALDEPS+=$(COMMON)/source.defs

include $(COMMON)/Makefile.inc

$(SOURCES) : $(LOCALDEPS)

sm_plugin.o : $(COMMON)/sm_var.mod $(COMMON)/sm_config.mod ../atomlist.def
sm_plugin_acemd.o : aceplug.h

aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h

clean::
	rm -f $(OBJS) $(HEADS)
