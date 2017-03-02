NAME=watershell
ROOT=..
HEADS=watershell.h
OBJS= watershell.o watershell_plugin.o watershell_plugin_acemd.o

DEPS=files dynamol continua output parser constants erf lu

LIB=libwatershell.so libwatershell.a
LIBS+=$(LIB)

include linkobjs.inc
include $(ROOT)/Makefile.inc

LOCALMACROS+=-include source.defs -D__BINDC -include erf_as.def #-D__erf=erfo7
LOCALDEPS+=source.defs

default : plugin
plugin : $(LIB)
$(LIB) : $(OBJS) $(COBJS) $(LOBJS) linkobjs.inc

$(SOURCES) : $(LOCALDEPS)

watershell_plugin.o : watershell.o
watershell.o :

watershell.F90 watershell_plugin.F90 : parser/parser_macro2.src
