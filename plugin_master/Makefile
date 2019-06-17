NAME=dynamo
ROOT=..
HEADS=plugin.h plugin.hpp
OBJS= plugin_master.o

DEPS=files dynamol continua output parser constants adaptive string watershell

LIB=lib$(NAME).so
LIBS+=$(LIB)

LOCALMACROS+=-include source.defs -D__BINDC 
LOCALDEPS+=source.defs

all : $(LIB)

include plugins.inc
include watershell.inc
include struna.inc
include adaptive.inc

include $(ROOT)/Makefile.inc
ifeq ($(strip($ACEMD_PLUGIN)),1)
 OBJS+=plugin_acemd.o
endif

$(LIB) : $(OBJS) $(COBJS) $(LOBJS)

$(SOURCES) : $(LOCALDEPS)

plugin_acemd.o : aceplug.h
aceplug.h:
	wget raw.githubusercontent.com/Acellera/plumed1/master/plumed-1.3/ACEMD/aceplug.h
plugin_master.F90 : parser/parser_macro2.src
clean ::
	rm -f *src
	